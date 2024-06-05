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

const WITHIN_GENE_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, ?1 - stranded_start
	FROM genes
 	WHERE level=?2 AND chr=?3 AND ((start <= ?4 AND end >= ?4) OR (start <= ?5 AND end >= ?5))
 	ORDER BY start ASC`

const CLOSEST_GENE_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, ?1 - stranded_start
	FROM genes
 	WHERE level = ?2 AND chr = ?3
 	ORDER BY ABS(stranded_start - ?1)
 	LIMIT ?4`

//  rows, err := genedb.closestGeneStmt.Query(mid,
// 	level,
// 	location.Chr,
// 	n)

const WITHIN_GENE_AND_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, ?1 - stranded_start 
	FROM genes 
 	WHERE level = ?2 AND chr = ?3 AND ((start - ?6 <= ?4 AND end + ?6 >= ?4) OR (start - ?6 <= ?5 AND end + ?6 >= ?5)) 
 	ORDER BY start ASC`

//  mid,
//  level,
//  location.Chr,
//  location.Start,
//  location.End,
//  pad)

const IN_EXON_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ?1 
	FROM genes 
 	WHERE level=3 AND gene_id=?2 AND chr=?3 AND ((start <= ?4 AND end >= ?4) OR (start <= ?5 AND end >= ?5)) 
 	ORDER BY start ASC`

const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? 
	FROM genes 
 	WHERE level=2 AND gene_id=? AND chr=? AND ? >= stranded_start - ? AND ? <= stranded_start + ? 
 	ORDER BY start ASC`

type GenomicFeature struct {
	Id       uint          `json:"-"`
	Location *dna.Location `json:"loc"`
	//Start      uint   `json:"start"`
	//End        uint   `json:"end"`
	Strand     string `json:"strand"`
	GeneId     string `json:"geneId"`
	GeneSymbol string `json:"geneSymbol"`
	TssDist    int    `json:"tssDist"`
	PromLabel  string `json:"promLabel"`
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
	Location string            `json:"location"`
	Level    string            `json:"level"`
	Features []*GenomicFeature `json:"features"`
}

//var ERROR_FEATURES = GenomicFeatures{Location: "", Level: "", Features: []GenomicFeature{}}

type Level uint8

const (
	LEVEL_GENE       Level = 1
	LEVEL_TRANSCRIPT Level = 2
	LEVEL_EXON       Level = 3
)

func ParseLevel(level string) Level {
	switch level {
	case "t", "transcript", "2":
		return LEVEL_TRANSCRIPT
	case "e", "exon", "3":
		return LEVEL_EXON
	default:
		return LEVEL_GENE
	}
}

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

type GeneDBCache struct {
	dir      string
	cacheMap map[string]*GeneDB
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

	return &GeneDBCache{dir, cacheMap}
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
	db                    *sql.DB
	withinGeneStmt        *sql.Stmt
	withinGeneAndPromStmt *sql.Stmt
	inExonStmt            *sql.Stmt
	closestGeneStmt       *sql.Stmt
}

func NewGeneDB(file string) *GeneDB {
	db := sys.Must(sql.Open("sqlite3", file))

	return &GeneDB{db: db,
		withinGeneStmt:        sys.Must(db.Prepare(WITHIN_GENE_SQL)),
		withinGeneAndPromStmt: sys.Must(db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL)),
		inExonStmt:            sys.Must(db.Prepare(IN_EXON_SQL)),
		closestGeneStmt:       sys.Must(db.Prepare(CLOSEST_GENE_SQL))}
}

func (genedb *GeneDB) Close() {
	genedb.db.Close()
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

	rows, err := genedb.withinGeneStmt.Query(
		mid,
		level,
		location.Chr,
		location.Start,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToRecords(location, rows, level)
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

	rows, err := genedb.withinGeneAndPromStmt.Query(
		mid,
		level,
		location.Chr,
		location.Start,
		location.End,
		pad)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToRecords(location, rows, level)
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

	rows, err := genedb.inExonStmt.Query(
		mid,
		geneId,
		location.Chr,
		location.Start,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToRecords(location, rows, LEVEL_EXON)
}

func (genedb *GeneDB) ClosestGenes(location *dna.Location, n uint16, level Level) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	// rows, err := genedb.closestGeneStmt.Query(mid,
	// 	level,
	// 	location.Chr,
	// 	mid,
	// 	n)

	rows, err := genedb.closestGeneStmt.Query(mid,
		level,
		location.Chr,
		n)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToRecords(location, rows, level)
}

func rowsToRecords(location *dna.Location, rows *sql.Rows, level Level) (*GenomicFeatures, error) {
	defer rows.Close()

	var id uint
	var chr string
	var start uint
	var end uint
	var strand string
	var geneId string
	var geneSymbol string
	var d int

	var features = make([]*GenomicFeature, 0, 10)

	for rows.Next() {
		err := rows.Scan(&id, &chr, &start, &end, &strand, &geneId, &geneSymbol, &d)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		location = dna.NewLocation(chr, start, end)

		feature := GenomicFeature{Id: id,
			Location:   location,
			Strand:     strand,
			GeneId:     geneId,
			GeneSymbol: geneSymbol,
			TssDist:    d}

		features = append(features, &feature)
	}

	return &GenomicFeatures{Location: location.String(), Level: level.String(), Features: features}, nil
}
