package genes

import (
	"database/sql"
	"fmt"
	"path/filepath"

	"github.com/antonybholmes/go-dna"
)

const WITHIN_GENE_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, ? - stranded_start
 FROM genes
 WHERE level=? AND chr=? AND ((start <= ? AND end >= ?) OR (start <= ? AND end >= ?))
 ORDER BY start ASC`

const CLOSEST_GENE_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, ? - stranded_start
 FROM genes
 WHERE level = ? AND chr = ?
 ORDER BY ABS(stranded_start - ?)
 LIMIT ?`

const WITHIN_GENE_AND_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, ? - stranded_start 
 FROM genes 
 WHERE level = ? AND chr = ? AND ((start - ? <= ? AND end + ? >= ?) OR (start - ? <= ? AND end + ? >= ?)) 
 ORDER BY start ASC`

const IN_EXON_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? 
 FROM genes 
 WHERE level=3 AND gene_id=? AND chr=? AND ((start <= ? AND end >= ?) OR (start <= ? AND end >= ?)) 
 ORDER BY start ASC`

const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? 
 FROM genes 
 WHERE level=2 AND gene_id=? AND chr=? AND ? >= stranded_start - ? AND ? <= stranded_start + ? 
 ORDER BY start ASC`

type GenomicFeature struct {
	Id         int    `json:"id"`
	Chr        string `json:"chr"`
	Start      uint   `json:"start"`
	End        uint   `json:"end"`
	Strand     string `json:"strand"`
	GeneId     string `json:"geneId"`
	GeneSymbol string `json:"geneSymbol"`
	TssDist    int    `json:"tssDist"`
}

func (feature *GenomicFeature) ToLocation() *dna.Location {
	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
}

func (feature *GenomicFeature) TSS() *dna.Location {
	var s uint

	if feature.Strand == "+" {
		s = feature.Start
	} else {
		s = feature.End
	}

	return dna.NewLocation(feature.Chr, s, s)
}

type GenomicFeatures struct {
	Location string           `json:"location"`
	Level    string           `json:"level"`
	Features []GenomicFeature `json:"features"`
}

var ERROR_FEATURES = GenomicFeatures{Location: "", Level: "", Features: []GenomicFeature{}}

type Level int

const (
	Gene       Level = 1
	Transcript Level = 2
	Exon       Level = 3
)

func ParseLevel(level string) Level {
	switch level {
	case "t", "transcript", "2":
		return Transcript
	case "e", "exon", "3":
		return Exon
	default:
		return Gene
	}
}

func (level Level) String() string {
	switch level {
	case Exon:
		return "Exon"
	case Transcript:
		return "Transcript"
	default:
		return "Gene"
	}
}

type GeneDbCache struct {
	dir   string
	cache *map[string]*GeneDb
}

func NewGeneDbCache() *GeneDbCache {
	return &GeneDbCache{}
}

func (genedbcache *GeneDbCache) Init(dir string) {
	genedbcache.dir = dir
	genedbcache.cache = new(map[string]*GeneDb)
}

func (genedbcache *GeneDbCache) Dir() string {
	return genedbcache.dir
}

func (genedbcache *GeneDbCache) Db(assembly string) (*GeneDb, error) {
	_, ok := (*genedbcache.cache)[assembly]

	if !ok {
		db, err := NewGeneDb(filepath.Join(genedbcache.dir, fmt.Sprintf("%s.db", assembly)))

		if err != nil {
			return nil, err
		}

		(*genedbcache.cache)[assembly] = db
	}

	return (*genedbcache.cache)[assembly], nil
}

func (genedbcache *GeneDbCache) Close() {
	for _, db := range *genedbcache.cache {
		db.Close()
	}
}

type GeneDb struct {
	db                    *sql.DB
	withinGeneStmt        *sql.Stmt
	withinGeneAndPromStmt *sql.Stmt
	inExonStmt            *sql.Stmt
	closestGeneStmt       *sql.Stmt
}

func NewGeneDb(file string) (*GeneDb, error) {
	db, err := sql.Open("sqlite3", file)

	if err != nil {
		return nil, err
	}

	withinGeneStmt, err := db.Prepare(WITHIN_GENE_SQL)

	if err != nil {
		return nil, err
	}

	withinGeneAndPromStmt, err := db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL)

	if err != nil {
		return nil, err
	}

	inExonStmt, err := db.Prepare(IN_EXON_SQL)

	if err != nil {
		return nil, err
	}

	closestGeneStmt, err := db.Prepare(CLOSEST_GENE_SQL)

	if err != nil {
		return nil, err
	}

	return &GeneDb{db, withinGeneStmt, withinGeneAndPromStmt, inExonStmt, closestGeneStmt}, nil
}

func (genedb *GeneDb) Close() {
	genedb.db.Close()
}

func (genedb *GeneDb) WithinGenes(location *dna.Location, level Level) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	rows, err := genedb.withinGeneStmt.Query(
		mid,
		level,
		location.Chr,
		location.Start,
		location.Start,
		location.End,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToRecords(location, rows, level)
}

func (genedb *GeneDb) WithinGenesAndPromoter(location *dna.Location, level Level, pad uint) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	rows, err := genedb.withinGeneAndPromStmt.Query(
		mid,
		level,
		location.Chr,
		pad,
		location.Start,
		pad,
		location.Start,
		pad,
		location.End,
		pad,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToRecords(location, rows, level)
}

func (genedb *GeneDb) InExon(location *dna.Location, geneId string) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	rows, err := genedb.inExonStmt.Query(
		mid,
		geneId,
		location.Chr,
		location.Start,
		location.Start,
		location.End,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToRecords(location, rows, Exon)
}

func (genedb *GeneDb) ClosestGenes(location *dna.Location, n uint16, level Level) (*GenomicFeatures, error) {
	mid := (location.Start + location.End) / 2

	rows, err := genedb.closestGeneStmt.Query(mid,
		level,
		location.Chr,
		mid,
		n)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToRecords(location, rows, level)
}

func rowsToRecords(location *dna.Location, rows *sql.Rows, level Level) (*GenomicFeatures, error) {
	defer rows.Close()

	var id int
	var chr string
	var start uint
	var end uint
	var strand string
	var geneId string
	var geneSymbol string
	var d int

	var features = []GenomicFeature{}

	for rows.Next() {
		err := rows.Scan(&id, &chr, &start, &end, &strand, &geneId, &geneSymbol, &d)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		features = append(features, GenomicFeature{Id: id, Chr: chr, Start: start, End: end, Strand: strand, GeneId: geneId, GeneSymbol: geneSymbol, TssDist: d})
	}

	return &GenomicFeatures{Location: location.String(), Level: level.String(), Features: features}, nil
}
