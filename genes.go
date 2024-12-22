package genes

import (
	"database/sql"
	"fmt"
	"sort"
	"strings"

	"os"
	"path/filepath"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys"
	"github.com/rs/zerolog/log"
)

const GENE_INFO_SQL = `SELECT id, chr, start, end, strand, gene_symbol, gene_id, transcript_id, 0 as tss_dist
	FROM genes
 	WHERE level = ?1 AND (gene_symbol = ?2 OR gene_id = ?2 OR transcript_id = ?2)
	ORDER BY chr, start`

const WITHIN_GENE_SQL = `SELECT id, chr, start, end, strand, gene_symbol, gene_id, transcript_id, ?1 - tss as tss_dist
	FROM genes
 	WHERE level = ?2 AND chr= ?3 AND ((start >= ?4 AND start <= ?5) OR (end >= ?4 AND end <= ?5) OR (start <= ?4 AND end >= ?5)) 
 	ORDER BY start`

const CLOSEST_GENE_SQL = `SELECT id, chr, start, end, strand, gene_symbol, gene_id, transcript_id, ?1 - tss as tss_dist
	FROM genes
 	WHERE level = ?2 AND chr = ?3
 	ORDER BY ABS(tss - ?1)
 	LIMIT ?4`

const WITHIN_GENE_AND_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_symbol, gene_id, transcript_id, ?1 - tss as tss_dist 
	FROM genes 
 	WHERE level = ?2 AND 
	chr = ?3 AND 
	(
		(strand = '+' AND ((start - ?6 >= ?4 AND start - ?6 <= ?5) OR (end <= ?4 AND end <= ?5) OR (start - ?6 <= ?4 AND end >= ?5))) OR
		(strand = '-' AND ((start <= ?4 AND start <= $5) OR (end + ?6 >= ?4 AND end + ?6 <= ?5) OR (start <= ?4 AND end + ?6 >= ?5)))
	)
 	ORDER BY start`

const IN_EXON_SQL = `SELECT id, chr, start, end, strand, gene_symbol, gene_id, transcript_id, ?1 - tss as tss_dist
	FROM genes 
 	WHERE level = 3 AND gene_id = ?2 AND chr = ?3 AND ((start >= ?4 AND start <= ?5) OR (end >= ?4 AND end <= ?5) OR (start <= ?4 AND end >= ?5)) 
 	ORDER BY start`

const OVERLAPPING_GENES_FROM_LOCATION_SQL = `SELECT id, level, chr, start, end, strand, gene_symbol, gene_id 
	FROM genes 
	WHERE level = 1 AND chr = ?1 ((start >= ?2 AND start <= ?3) OR (end >= ?2 AND end <= ?3) OR (start <= ?2 AND end >= ?3))
	ORDER BY start`

const OVERLAPPING_TRANSCRIPTS_FROM_LOCATION_SQL = `SELECT id, level, chr, start, end, strand, transcript_id 
	FROM genes 
	WHERE level = 2 AND gene_id = ?1
	ORDER BY start`

const OVERLAPPING_EXONS_FROM_LOCATION_SQL = `SELECT id, level, chr, start, end, strand, exon_id 
	FROM genes 
	WHERE level = 3 AND transcript_id = ?1
	ORDER BY start`

// const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_symbol, gene_id, transcript_id, start - ?
// 	FROM genes
//  	WHERE level = 2 AND gene_id = ? AND chr = ? AND ? >= stranded_start - ? AND ? <= stranded_start + ?
//  	ORDER BY start ASC`

type GenomicFeature struct {
	Location     *dna.Location     `json:"loc"`
	Strand       string            `json:"strand"`
	Level        string            `json:"level"`
	GeneSymbol   string            `json:"geneSymbol,omitempty"`
	GeneId       string            `json:"geneId,omitempty"`
	TranscriptId string            `json:"transcriptId,omitempty"`
	ExonId       string            `json:"exonId,omitempty"`
	PromLabel    string            `json:"promLabel,omitempty"`
	Children     []*GenomicFeature `json:"children,omitempty"`
	//Id           uint              `json:"-"`
	TssDist int `json:"tssDist,omitempty"`
}

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func (feature *GenomicFeature) TSS() *dna.Location {
	var s uint

	if feature.Strand == "+" {
		s = feature.Location.Start
	} else {
		s = feature.Location.End
	}

	return dna.NewLocation(feature.Location.Chr, s, s)
}

type GenomicFeatures struct {
	Location *dna.Location `json:"location"`
	//Id       string        `json:"id,omitempty"`
	//Name     string        `json:"name,omitempty"`
	Level    string            `json:"level"`
	Features []*GenomicFeature `json:"features"`
}

type Level uint8

const (
	LEVEL_GENE       Level = 1
	LEVEL_TRANSCRIPT Level = 2
	LEVEL_EXON       Level = 3
)

func (level Level) String() string {
	switch level {
	case LEVEL_EXON:
		return "Exon"
	case LEVEL_TRANSCRIPT:
		return "Transcript"
	default:
		return "Gene"
	}
}

func ParseLevel(level string) Level {
	switch level {
	case "t", "tran", "transcript", "2":
		return LEVEL_TRANSCRIPT
	case "e", "ex", "exon", "3":
		return LEVEL_EXON
	default:
		return LEVEL_GENE
	}
}

type GeneDBCache struct {
	cacheMap map[string]*GeneDB
	dir      string
}

func NewGeneDBCache(dir string) *GeneDBCache {
	cacheMap := make(map[string]*GeneDB)

	files, err := os.ReadDir(dir)

	log.Debug().Msgf("---- genedb ----")

	if err != nil {
		log.Fatal().Msgf("error opening %s", dir)
	}

	log.Debug().Msgf("caching gene databases in %s...", dir)

	for _, file := range files {
		basename := file.Name()

		if strings.HasSuffix(basename, ".db") {

			name := strings.TrimSuffix(basename, filepath.Ext(basename))
			db := NewGeneDB(filepath.Join(dir, basename))
			cacheMap[name] = db

			log.Debug().Msgf("found gene database %s", name)
		}
	}

	log.Debug().Msgf("---- end ----")

	return &GeneDBCache{dir: dir, cacheMap: cacheMap}
}

func (cache *GeneDBCache) Dir() string {
	return cache.dir
}

func (cache *GeneDBCache) List() []string {

	ids := make([]string, 0, len(cache.cacheMap))

	for id := range cache.cacheMap {
		ids = append(ids, id)
	}

	sort.Strings(ids)

	return ids
}

func (cache *GeneDBCache) GeneDB(assembly string) (*GeneDB, error) {
	_, ok := cache.cacheMap[assembly]

	if !ok {
		db := NewGeneDB(filepath.Join(cache.dir, fmt.Sprintf("%s.db", assembly)))

		cache.cacheMap[assembly] = db
	}

	return cache.cacheMap[assembly], nil
}

func (cache *GeneDBCache) Close() {
	for _, db := range cache.cacheMap {
		db.Close()
	}
}

type GeneDB struct {
	db *sql.DB
	//withinGeneStmt *sql.Stmt
	//withinGeneAndPromStmt *sql.Stmt
	//inExonStmt      *sql.Stmt
	//closestGeneStmt *sql.Stmt
}

func NewGeneDB(file string) *GeneDB {
	db := sys.Must(sql.Open("sqlite3", file))

	return &GeneDB{db: db} //withinGeneStmt: sys.Must(db.Prepare(WITHIN_GENE_SQL)),
	//withinGeneAndPromStmt: sys.Must(db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL)),
	//inExonStmt:      sys.Must(db.Prepare(IN_EXON_SQL)),
	//closestGeneStmt: sys.Must(db.Prepare(CLOSEST_GENE_SQL))

}

func (genedb *GeneDB) Close() {
	genedb.db.Close()
}

func (genedb *GeneDB) OverlappingGenes(location *dna.Location) ([]*GenomicFeature, error) {

	var id uint
	var level Level
	var chr string
	var start uint
	var end uint
	var strand string
	var geneSymbol string
	var geneId string
	var transcriptId string
	var exonId string

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation

	var features = make([]*GenomicFeature, 0, 10)
	var currentGene *GenomicFeature
	var currentTranscript *GenomicFeature

	geneRows, err := genedb.db.Query(OVERLAPPING_GENES_FROM_LOCATION_SQL,
		location.Chr,
		location.Start,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer geneRows.Close()

	//var currentExon *GenomicSearchFeature

	for geneRows.Next() {
		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneSymbol, &geneId, &transcriptId, &exonId)
		err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneSymbol, &geneId)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		location := &dna.Location{Chr: chr, Start: start, End: end}

		feature := &GenomicFeature{Location: location,
			Level:      Level(1).String(),
			GeneSymbol: geneSymbol,
			GeneId:     geneId,
			Strand:     strand,
			Children:   make([]*GenomicFeature, 0, 10)}

		features = append(features, feature)
		currentGene = feature

		transcriptRows, err := genedb.db.Query(OVERLAPPING_TRANSCRIPTS_FROM_LOCATION_SQL,
			currentGene.GeneId)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database query")
		}

		defer transcriptRows.Close()

		//var currentExon *GenomicSearchFeature

		for transcriptRows.Next() {
			err := transcriptRows.Scan(&id, &level, &chr, &start, &end, &strand, &transcriptId)

			if err != nil {
				return nil, err //fmt.Errorf("there was an error with the database records")
			}

			location := &dna.Location{Chr: chr, Start: start, End: end}

			feature := &GenomicFeature{Location: location,
				Level:        Level(2).String(),
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				Strand:       strand,
				Children:     make([]*GenomicFeature, 0, 10)}
			currentGene.Children = append(currentGene.Children, feature)
			currentTranscript = feature

			exonRows, err := genedb.db.Query(OVERLAPPING_EXONS_FROM_LOCATION_SQL,
				currentTranscript.TranscriptId)

			if err != nil {
				return nil, err //fmt.Errorf("there was an error with the database query")
			}

			defer exonRows.Close()

			//var currentExon *GenomicSearchFeature

			for exonRows.Next() {
				err := exonRows.Scan(&id, &level, &chr, &start, &end, &strand, &exonId)

				if err != nil {
					return nil, err //fmt.Errorf("there was an error with the database records")
				}

				location := &dna.Location{Chr: chr, Start: start, End: end}
				feature := &GenomicFeature{Location: location,
					Level:        Level(3).String(),
					GeneSymbol:   geneSymbol,
					GeneId:       geneId,
					TranscriptId: transcriptId,
					ExonId:       exonId,
					Strand:       strand}
				currentTranscript.Children = append(currentTranscript.Children, feature)

			}
		}

	}

	return features, nil
}

func (genedb *GeneDB) GeneInfo(search string, level Level) ([]*GenomicFeature, error) {

	rows, err := genedb.db.Query(GENE_INFO_SQL,

		level,
		search)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	features, err := rowsToRecords(rows, level)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return features, nil
}

func (genedb *GeneDB) WithinGenes(location *dna.Location, level Level) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	// rows, err := genedb.withinGeneStmt.Query(
	// 	mid,
	// 	level,
	// 	location.Chr,
	// 	location.Start,
	// 	location.Start,
	// 	location.End,
	// 	location.End)

	rows, err := genedb.db.Query(WITHIN_GENE_SQL,
		mid,
		level,
		location.Chr,
		location.Start,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, level)
}

func (genedb *GeneDB) WithinGenesAndPromoter(location *dna.Location, level Level, pad uint) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	// rows, err := genedb.withinGeneAndPromStmt.Query(
	// 	mid,
	// 	level,
	// 	location.Chr,
	// 	pad,
	// 	location.Start,
	// 	pad,
	// 	location.Start,
	// 	pad,
	// 	location.End,
	// 	pad,
	// 	location.End)

	rows, err := genedb.db.Query(WITHIN_GENE_AND_PROMOTER_SQL,
		mid,
		level,
		location.Chr,
		location.Start,
		location.End,
		pad)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, level)
}

func (genedb *GeneDB) InExon(location *dna.Location, geneId string) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	// rows, err := genedb.inExonStmt.Query(
	// 	mid,
	// 	geneId,
	// 	location.Chr,
	// 	location.Start,
	// 	location.Start,
	// 	location.End,
	// 	location.End)

	rows, err := genedb.db.Query(IN_EXON_SQL,
		mid,
		geneId,
		location.Chr,
		location.Start,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, LEVEL_EXON)
}

func (genedb *GeneDB) ClosestGenes(location *dna.Location, n uint16, level Level) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	// rows, err := genedb.closestGeneStmt.Query(mid,
	// 	level,
	// 	location.Chr,
	// 	mid,
	// 	n)

	rows, err := genedb.db.Query(CLOSEST_GENE_SQL, mid,
		level,
		location.Chr,
		n)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, level)
}

func rowsToFeatures(location *dna.Location, rows *sql.Rows, level Level) (*GenomicFeatures, error) {

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	features, err := rowsToRecords(rows, level)

	if err != nil {
		return nil, err
	}

	return &GenomicFeatures{Location: location, Level: level.String(), Features: features}, nil
}

func rowsToRecords(rows *sql.Rows, level Level) ([]*GenomicFeature, error) {
	defer rows.Close()

	var id uint
	var chr string
	var start uint
	var end uint
	var strand string
	var geneId string
	var geneSymbol string
	var transcriptId string
	var d int

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	var features = make([]*GenomicFeature, 0, 10)

	for rows.Next() {
		err := rows.Scan(&id, &chr, &start, &end, &strand, &geneSymbol, &geneId, &transcriptId, &d)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		location := dna.NewLocation(chr, start, end)

		feature := GenomicFeature{
			Location:     location,
			Strand:       strand,
			GeneId:       geneId,
			GeneSymbol:   geneSymbol,
			TranscriptId: transcriptId,
			Level:        level.String(),
			TssDist:      d}

		features = append(features, &feature)
	}

	return features, nil
}
