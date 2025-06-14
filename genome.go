package genome

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

const GENE_DB_INFO_SQL = `SELECT id, genome, version FROM info`

const GENE_INFO_SQL = `SELECT g.id, g.chr, g.start, g.end, g.strand, g.gene_id, g.gene_symbol, transcript_id, gene_type, 0 as tss_dist
	FROM genes
 	WHERE (gene_symbol LIKE ?2 OR g.gene_id LIKE ?2 OR g.transcript_id LIKE ?2)
	ORDER BY g.chr, g.start
	LIMIT ?3`

// const GENE_INFO_FUZZY_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, transcript_id, gene_type, 0 as tss_dist
//  	FROM genes
//   	WHERE (gene_symbol LIKE ?2 OR g.gene_id LIKE ?2 OR g.transcript_id LIKE ?2)
//  	ORDER BY chr, start
// 	LIMIT ?3`

const WITHIN_GENE_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, transcript_id, gene_type, ?1 - tss as tss_dist
	FROM genes
 	WHERE level = ?2 AND chr= ?3 AND (start <= ?5 AND end >= ?4)
 	ORDER BY start`

const CLOSEST_GENE_SQL = `SELECT g.id, g.chr, g.start, g.end, g.strand, g.gene_symbol, g.gene_id, g.gene_type, ?1 - tss as g.tss_dist
	FROM genes AS g
 	WHERE chr = ?3
 	ORDER BY ABS(g.tss - ?1)
 	LIMIT ?4`

const CLOSEST_TRANSCRIPT_SQL = `SELECT t.id, t.chr, t.start, t.end, g.strand, g.gene_id, g.gene_symbol, t.transcript_id, t.gene_type, ?1 - g.tss as tss_dist
	FROM transcripts AS t
	INNER JOIN genes AS g ON g.id = t.gene_id
 	WHERE g.chr = ?3
 	ORDER BY ABS(g.tss - ?1)
 	LIMIT ?4`

const WITHIN_GENE_AND_PROMOTER_SQL = `SELECT t.id, t.chr, t.start, t.end, g.strand, g.gene_id, g.gene_symbol, t.transcript_id, g.gene_type, ?1 - g.tss as tss_dist 
	FROM genes AS g
	INNER JOIN transcripts AS t ON g.id = t.gene_id
 	WHERE g.chr = ?3 AND 
	(
		((t.start - ?6 <= ?5 AND t.end >= ?4) AND g.strand = '+') OR
		((t.start >= ?5 AND t.end + ?6 <= ?4) AND g.strand = '-')
	)
 	ORDER BY t.start`

const IN_EXON_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, transcript_id, gene_type, ?1 - tss as tss_dist
	FROM genes 
 	WHERE level = 3 AND gene_id = ?2 AND chr = ?3 AND (start <= ?5 AND end >= ?4)
 	ORDER BY start`

// If start is less x2 and end is greater than x1, it constrains the feature to be overlapping
// our location
const OVERLAPPING_GENES_FROM_LOCATION_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, is_canonical, gene_type
	FROM genes AS g
	WHERE g.chr = ?1 AND (g.start <= ?3 AND g.end >= ?2)`

const TRANSCRIPTS_IN_GENE_SQL = `SELECT t.id, g.chr, t.start, t.end, g.strand, t.transcript_id, g.is_canonical, g.gene_type
	FROM transcripts AS t
	INNER JOIN genes AS g ON g.id = t.gene_id
	WHERE t.gene_id = ?1
	ORDER BY t.start`

// const CANONICAL_TRANSCRIPTS_IN_GENE_SQL = `SELECT id, level, chr, start, end, strand, transcript_id, is_canonical, gene_type
// 	FROM genes
// 	WHERE level = 2 AND gene_id = ?1 AND is_canonical = 1
// 	ORDER BY start`

const EXONS_IN_TRANSCRIPT_SQL = `SELECT e.id, g.chr, e.start, e.end, g.strand, e.exon_id, g.is_canonical
	FROM exons AS e
	INNER JOIN transcripts AS t ON t.id = e.transcript_id
	INNER JOIN genes AS g ON g.id = t.gene_id
	WHERE t.transcript_id = ?1
	ORDER BY e.start`

const ID_TO_NAME_SQL = `SELECT name FROM ids WHERE id = ?1`

const MAX_GENE_INFO_RESULTS uint16 = 100

// const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, transcript_id, start - ?
// 	FROM genes
//  	WHERE level = 2 AND gene_id = ? AND chr = ? AND ? >= stranded_start - ? AND ? <= stranded_start + ?
//  	ORDER BY start ASC`

type GenomicFeature struct {
	Location     *dna.Location     `json:"loc"`
	Strand       string            `json:"strand"`
	Level        Level             `json:"level"`
	GeneSymbol   string            `json:"geneSymbol,omitempty"`
	GeneId       string            `json:"geneId,omitempty"`
	TranscriptId string            `json:"transcriptId,omitempty"`
	ExonId       string            `json:"exonId,omitempty"`
	PromLabel    string            `json:"promLabel,omitempty"`
	Children     []*GenomicFeature `json:"children,omitempty"`
	TssDist      int               `json:"tssDist,omitempty"`
	IsCanonical  bool              `json:"isCanonical"`
	GeneType     string            `json:"geneType,omitempty"`
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
	Level    Level             `json:"level"`
	Features []*GenomicFeature `json:"features"`
}

type Level string

const (
	LEVEL_GENE       Level = "gene"
	LEVEL_TRANSCRIPT Level = "transcript"
	LEVEL_EXON       Level = "exon"
)

// func (level Level) String() string {
// 	switch level {
// 	case LEVEL_EXON:
// 		return "Exon"
// 	case LEVEL_TRANSCRIPT:
// 		return "Transcript"
// 	default:
// 		return "Gene"
// 	}
// }

// func ParseLevel(level string) Level {
// 	switch level {
// 	case "t", "tran", "transcript", "2":
// 		return LEVEL_TRANSCRIPT
// 	case "e", "ex", "exon", "3":
// 		return LEVEL_EXON
// 	default:
// 		return LEVEL_GENE
// 	}
// }

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
			db := NewGeneDB(name, filepath.Join(dir, basename))
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

func (cache *GeneDBCache) List() ([]*GeneDBInfo, error) {

	ids := make([]string, 0, len(cache.cacheMap))

	for id := range cache.cacheMap {

		ids = append(ids, id)
	}

	sort.Strings(ids)

	infos := make([]*GeneDBInfo, 0, len(ids))

	for _, id := range ids {
		info, err := cache.cacheMap[id].GeneDBInfo()

		if err != nil {
			return nil, err
		}

		infos = append(infos, info)
	}

	return infos, nil
}

func (cache *GeneDBCache) GeneDB(assembly string) (*GeneDB, error) {
	_, ok := cache.cacheMap[assembly]

	if !ok {
		db := NewGeneDB(assembly, filepath.Join(cache.dir, fmt.Sprintf("%s.db", assembly)))

		cache.cacheMap[assembly] = db
	}

	return cache.cacheMap[assembly], nil
}

func (cache *GeneDBCache) Close() {
	for _, db := range cache.cacheMap {
		db.Close()
	}
}

type GeneDBInfo struct {
	Name    string `json:"name"`
	Genome  string `json:"genome"`
	Version string `json:"version"`
}

type GeneDB struct {
	db *sql.DB
	//withinGeneAndPromStmt *sql.Stmt
	name string
}

func NewGeneDB(name string, file string) *GeneDB {
	db := sys.Must(sql.Open("sqlite3", file))

	return &GeneDB{name: name, db: db} //withinGeneStmt: sys.Must(db.Prepare(WITHIN_GENE_SQL)),
	//withinGeneAndPromStmt: sys.Must(db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL))

	//inExonStmt:      sys.Must(db.Prepare(IN_EXON_SQL)),
	//closestGeneStmt: sys.Must(db.Prepare(CLOSEST_GENE_SQL))

}

func (genedb *GeneDB) Close() {
	genedb.db.Close()
}

func (genedb *GeneDB) GeneDBInfo() (*GeneDBInfo, error) {
	var id uint
	var genome string
	var version string

	err := genedb.db.QueryRow(GENE_DB_INFO_SQL).Scan(&id, &genome, &version)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return &GeneDBInfo{Name: genedb.name, Genome: genome, Version: version}, nil
}

// comprehensive gene,transcript exon search for a given location
func (genedb *GeneDB) OverlappingGenes(location *dna.Location, canonical bool) ([]*GenomicFeature, error) {

	var id uint

	var chr string
	var start uint
	var end uint
	var strand string
	var geneSymbol string
	var geneId string
	var transcriptId string
	var exonId string
	var isCanonical bool
	var geneType string

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	var features = make([]*GenomicFeature, 0, 10)
	var currentGene *GenomicFeature
	var currentTranscript *GenomicFeature

	genesSql := OVERLAPPING_GENES_FROM_LOCATION_SQL

	if canonical {
		genesSql += " AND g.is_canonical = 1"
	}

	genesSql += " ORDER BY g.start"

	geneRows, err := genedb.db.Query(genesSql,
		location.Chr,
		location.Start,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer geneRows.Close()

	//var currentExon *GenomicSearchFeature

	for geneRows.Next() {
		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &exonId)
		err := geneRows.Scan(&id, &chr, &start, &end, &strand, &geneId, &geneSymbol, &isCanonical, &geneType)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		location := &dna.Location{Chr: chr, Start: start, End: end}

		feature := &GenomicFeature{Location: location,
			Level:       LEVEL_GENE,
			GeneSymbol:  geneSymbol,
			GeneId:      geneId,
			Strand:      strand,
			IsCanonical: isCanonical,
			GeneType:    geneType,
			Children:    make([]*GenomicFeature, 0, 10)}

		features = append(features, feature)
		currentGene = feature

		transcriptRows, err := genedb.db.Query(TRANSCRIPTS_IN_GENE_SQL,
			currentGene.GeneId)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database query")
		}

		defer transcriptRows.Close()

		//var currentExon *GenomicSearchFeature

		for transcriptRows.Next() {
			err := transcriptRows.Scan(&id, &chr, &start, &end, &strand, &transcriptId, &isCanonical, &geneType)

			if err != nil {
				return nil, err //fmt.Errorf("there was an error with the database records")
			}

			location := &dna.Location{Chr: chr, Start: start, End: end}

			feature := &GenomicFeature{Location: location,
				Level:        "transcript",
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				Strand:       strand,
				IsCanonical:  isCanonical,
				GeneType:     geneType,
				Children:     make([]*GenomicFeature, 0, 10)}
			currentGene.Children = append(currentGene.Children, feature)
			currentTranscript = feature

			exonRows, err := genedb.db.Query(EXONS_IN_TRANSCRIPT_SQL,
				currentTranscript.TranscriptId)

			if err != nil {
				return nil, err //fmt.Errorf("there was an error with the database query")
			}

			defer exonRows.Close()

			//var currentExon *GenomicSearchFeature

			for exonRows.Next() {
				err := exonRows.Scan(&id, &chr, &start, &end, &strand, &exonId, &isCanonical, &geneType)

				if err != nil {
					return nil, err //fmt.Errorf("there was an error with the database records")
				}

				location := &dna.Location{Chr: chr, Start: start, End: end}
				feature := &GenomicFeature{Location: location,
					Level:        "exon",
					GeneSymbol:   geneSymbol,
					GeneId:       geneId,
					TranscriptId: transcriptId,
					ExonId:       exonId,
					Strand:       strand,
					IsCanonical:  isCanonical,
					GeneType:     geneType}
				currentTranscript.Children = append(currentTranscript.Children, feature)

			}
		}

	}

	return features, nil
}

func (genedb *GeneDB) SearchForGeneByName(search string, level Level, n uint16, fuzzy bool) ([]*GenomicFeature, error) {
	n = max(1, min(n, MAX_GENE_INFO_RESULTS))

	ret := make([]*GenomicFeature, 0, n)

	// case insensitive search
	search = strings.ToLower(search)

	if len(search) < 2 || strings.Contains(search, "chr:") {
		return ret, fmt.Errorf("%s is an invalid search term", search)
	}

	//log.Debug().Msgf("search %s", search)

	var rows *sql.Rows
	var err error

	if fuzzy && !strings.HasSuffix(search, "%") {
		search += "%"
	}

	rows, err = genedb.db.Query(GENE_INFO_SQL,
		level,
		search,
		n)

	if err != nil {
		// return empty array
		return ret, err //fmt.Errorf("there was an error with the database query")
	}

	err = rowsToRecords(rows, level, &ret)

	if err != nil {
		return ret, err //fmt.Errorf("there was an error with the database query")
	}

	//sort.Sort(SortFeatureByPos(ret))

	return ret, nil
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

	return rowsToFeatures(location, rows, LEVEL_TRANSCRIPT)
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

	var sql string
	if level == LEVEL_TRANSCRIPT || level == LEVEL_EXON {
		sql = CLOSEST_TRANSCRIPT_SQL
	} else {
		sql = CLOSEST_GENE_SQL
	}

	rows, err := genedb.db.Query(sql, mid,
		level,
		location.Chr,
		n)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, level)
}

func rowsToFeatures(location *dna.Location, rows *sql.Rows, level Level) (*GenomicFeatures, error) {
	ret := GenomicFeatures{Location: location, Level: level, Features: make([]*GenomicFeature, 0, 10)}

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	err := rowsToRecords(rows, level, &ret.Features)

	if err != nil {
		return &ret, err
	}

	return &ret, nil
}

// func (genedb *GeneDB) IdToName(id uint) string {

// 	name := "n/a"

// 	// ignore error as it is not critical
// 	// if the id is not found, we return "n/a"
// 	genedb.db.QueryRow(ID_TO_NAME_SQL, id).Scan(&name)

// 	return name
// }

func rowsToRecords(rows *sql.Rows, level Level, features *[]*GenomicFeature) error {
	defer rows.Close()

	var id uint
	var chr string
	var start uint
	var end uint
	var strand string
	var geneId string
	var geneSymbol string
	var transcriptId string
	var geneType string
	var d int
	var err error

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	//var features = make([]*GenomicFeature, 0, 10)

	for rows.Next() {
		geneId = ""
		geneSymbol = ""
		transcriptId = ""
		geneType = ""
		d = 0 // tss distance

		switch level {
		case LEVEL_TRANSCRIPT, LEVEL_EXON:
			err = rows.Scan(&id, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &geneType, &d)
		default:
			err = rows.Scan(&id, &chr, &start, &end, &strand, &geneId, &geneSymbol, &geneType)
		}

		if err != nil {
			return err //fmt.Errorf("there was an error with the database records")
		}

		location := dna.NewLocation(chr, start, end)

		feature := GenomicFeature{
			Level:        level,
			Location:     location,
			Strand:       strand,
			GeneId:       geneId,
			GeneSymbol:   geneSymbol,
			TranscriptId: transcriptId,
			GeneType:     geneType,
			TssDist:      d}

		*features = append(*features, &feature)
	}

	// enforce sorted correctly by chr and then position
	sort.Sort(SortFeatureByPos(*features))

	return nil
}

type SortFeatureByPos []*GenomicFeature

func (features SortFeatureByPos) Len() int      { return len(features) }
func (features SortFeatureByPos) Swap(i, j int) { features[i], features[j] = features[j], features[i] }
func (features SortFeatureByPos) Less(i, j int) bool {
	ci := dna.ChromToInt(features[i].Location.Chr)
	cj := dna.ChromToInt(features[j].Location.Chr)

	// on different chrs so sort by chr
	if ci != cj {
		return ci < cj
	}

	return features[i].Location.Start < features[j].Location.Start
}
