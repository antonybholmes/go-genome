package genome

import (
	"database/sql"
	"fmt"
	"slices"
	"sort"
	"strings"

	"os"
	"path/filepath"

	"github.com/antonybholmes/go-basemath"
	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys"
	"github.com/rs/zerolog/log"
)

const GENE_DB_INFO_SQL = `SELECT id, genome, version FROM info`

const STANDARD_GENE_FIELDS = `SELECT DISTINCT
	g.id, 
	g.feature,
	g.seqname, 
	g.start, 
	g.end, 
	g.strand, 
	g.gene_id, 
	g.gene_name, 
	g.type,
	g.is_canonical,
	g.transcript_id, 
	g.transcript_name,
	g.exon_number`

const GENE_INFO_SQL = STANDARD_GENE_FIELDS +
	`, 0 as tss_dist
	FROM gtf AS g
 	WHERE g.feature = 'gene' AND (g.gene_name LIKE ?1 OR g.gene_id LIKE ?1)`

const TRANSCRIPT_INFO_SQL = STANDARD_GENE_FIELDS +
	` FROM gtf AS g
 	WHERE g.feature = 'transcript' AND (g.gene_name LIKE ?1 OR g.gene_id LIKE ?1 OR g.transcript_id LIKE ?1)`

// const GENE_INFO_FUZZY_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_name, transcript_id, gene_type, 0 as tss_dist
//  	FROM gene
//   	WHERE (gene_name LIKE ?2 OR g.gene_id LIKE ?2 OR g.transcript_id LIKE ?2)
//  	ORDER BY chr, start
// 	LIMIT ?3`

const WITHIN_GENE_SQL = STANDARD_GENE_FIELDS +
	`, ?3 - g.tss as tss_dist
	FROM gtf AS g
 	WHERE g.feature = ?1 AND g.seqname = ?2 AND (g.start <= ?5 AND g.end >= ?4)
 	ORDER BY g.start`

const CLOSEST_GENE_SQL = STANDARD_GENE_FIELDS +
	`, ?3 - g.tss as tss_dist
	FROM gtf AS g
 	WHERE g.feature = ?1 AND g.seqname = ?2
 	ORDER BY ABS(tss_dist)
 	LIMIT ?4`

const TRANSCRIPTS_PER_GENE = STANDARD_GENE_FIELDS +
	`, ?3 - g.tss as tss_dist
	FROM gtf AS g
 	WHERE g.feature = ?1 AND g.seqname = ?2
 	ORDER BY ABS(tss_dist)
 	LIMIT ?4`

// const CLOSEST_TRANSCRIPT_SQL = `SELECT
// 	g.id,
// 	g.feature,
// 	g.seqname,
// 	t.start,
// 	t.end,
// 	g.strand,
// 	g.gene_id,
// 	g.gene_name,
// 	gt.name as gene_type,
// 	t.transcript_id,
// 	t.is_canonical,
// 	tt.name as transcript_type,
// 	?1 - g.tss as tss_dist
// 	FROM transcript AS t
// 	INNER JOIN gene AS g ON g.id = t.gene_id
// 	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
// 	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
//  	WHERE g.seqname = ?1
//  	ORDER BY ABS(g.tss - ?2)
//  	LIMIT ?3`

const WITHIN_GENE_AND_PROMOTER_SQL = STANDARD_GENE_FIELDS +
	`, ?2 - g.tss as tss_dist 
	FROM gtf AS g
 	WHERE g.feature = 'transcript' AND g.seqname = ?1 AND 
	(
		((g.start - ?5 <= ?4 AND g.end + ?6 >= ?3) AND g.strand = '+') OR
		((g.start - ?6 >= ?4 AND g.end + ?5 <= ?3) AND g.strand = '-')
	)
 	ORDER BY tss_dist`

// when annotating genes, see if position falls within an exon
const IN_EXON_SQL = STANDARD_GENE_FIELDS +
	`, ?2 - g.tss as tss_dist
	FROM gtf AS g
 	WHERE g.feature = 'exon' AND g.transcript_id = ?1 AND (g.start <= ?4 AND g.end >= ?3)
 	ORDER BY g.start, g.end DESC`

// If start is less x2 and end is greater than x1, it constrains the feature to be overlapping
// our location
// const OVERLAPPING_GENES_FROM_LOCATION_SQL = `SELECT
// 	g.id,
// 	g.seqname,
// 	g.start,
// 	g.end,
// 	g.strand,
// 	g.gene_id,
// 	g.gene_name,
// 	gt.name as gene_type
// 	FROM gene AS g
// 	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
// 	WHERE g.seqname = ?1 AND (g.start <= ?3 AND g.end >= ?2)`

// get all genes, transcripts, and exons overlapping a location
// which we can use to build a nested gene structure
const OVERLAP_LOCATION_SQL = STANDARD_GENE_FIELDS +
	` FROM gtf AS g
	WHERE g.seqname = ?1 AND (g.start <= ?3 AND g.end >= ?2)`

const OVERLAP_ORDER_BY_SQL = ` ORDER BY 
	g.gene_id,
	g.transcript_id,
	g.exon_number`

// const TRANSCRIPTS_IN_GENE_SQL = `SELECT
// 	t.id,
// 	g.seqname,
// 	t.start,
// 	t.end,
// 	g.strand,
// 	gt.name as gene_type,
// 	t.transcript_id,
// 	t.is_canonical,
// 	tt.name as transcript_type
// 	FROM transcript AS t
// 	INNER JOIN gene AS g ON g.id = t.gene_id
// 	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
// 	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
// 	WHERE t.gene_id = ?1`

// const CANONICAL_TRANSCRIPTS_IN_GENE_SQL = `SELECT id, level, chr, start, end, strand, transcript_id, is_canonical, gene_type
// 	FROM gene
// 	WHERE level = 2 AND gene_id = ?1 AND is_canonical = 1
// 	ORDER BY start`

// const EXONS_IN_TRANSCRIPT_SQL = `SELECT
// 	e.id,
// 	g.seqname,
// 	e.start,
// 	e.end,
// 	g.strand,
// 	gt.name as gene_type,
// 	e.exon_id,
// 	t.is_canonical,
// 	tt.name as transcript_type
// 	FROM exon AS e
// 	INNER JOIN transcript AS t ON t.id = e.transcript_id
// 	INNER JOIN gene AS g ON g.id = t.gene_id
// 	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
// 	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
// 	WHERE e.transcript_id = ?1`

const ID_TO_NAME_SQL = `SELECT name FROM ids WHERE id = ?1`

const MAX_GENE_INFO_RESULTS uint16 = 100

// const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_name, transcript_id, start - ?
// 	FROM gene
//  	WHERE level = 2 AND gene_id = ? AND chr = ? AND ? >= stranded_start - ? AND ? <= stranded_start + ?
//  	ORDER BY start ASC`

type GenomicFeature struct {
	Location       *dna.Location     `json:"loc"`
	Type           string            `json:"type,omitempty"`
	Feature        Feature           `json:"feature"`
	GeneId         string            `json:"geneId,omitempty"`
	GeneSymbol     string            `json:"geneSymbol,omitempty"`
	TranscriptId   string            `json:"transcriptId,omitempty"`
	TranscriptName string            `json:"transcriptName,omitempty"`
	ExonId         string            `json:"exonId,omitempty"`
	Children       []*GenomicFeature `json:"children,omitempty"`
	ExonNumber     uint              `json:"exonNumber,omitempty"`
	Id             uint              `json:"-"`
	TssDist        int               `json:"tssDist,omitempty"`
	IsCanonical    bool              `json:"isCanonical,omitempty"`
	IsLongest      bool              `json:"isLongest,omitempty"`
}

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func (feature *GenomicFeature) TSS() *dna.Location {
	var s uint

	if feature.Location.Strand == "+" {
		s = feature.Location.Start
	} else {
		s = feature.Location.End
	}

	return dna.NewStrandedLocation(feature.Location.Chr, s, s, feature.Location.Strand)
}

type GenomicFeatures struct {
	Location *dna.Location `json:"location"`
	//Id       string        `json:"id,omitempty"`
	//Name     string        `json:"name,omitempty"`
	Feature  Feature           `json:"feature"`
	Features []*GenomicFeature `json:"features"`
}

type Feature string

const (
	FEATURE_GENE       Feature = "gene"
	FEATURE_TRANSCRIPT Feature = "transcript"
	FEATURE_EXON       Feature = "exon"
)

type Level uint8

const (
	LEVEL_GENE       Level = 1
	LEVEL_TRANSCRIPT Level = 2
	LEVEL_EXON       Level = 3
)

func LevelToFeature(level Level) Feature {

	switch level {
	case LEVEL_TRANSCRIPT:
		return FEATURE_TRANSCRIPT
	case LEVEL_EXON:
		return FEATURE_EXON
	default:
		return FEATURE_GENE
	}
}

func FeatureToLevel(feature string) Level {
	switch feature {
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

		if strings.Contains(basename, "gtf") && strings.HasSuffix(basename, ".db") {

			name := strings.TrimSuffix(basename, filepath.Ext(basename))
			name = strings.TrimPrefix(name, "gtf_")
			db := NewGeneDB(name, filepath.Join(dir, basename))
			cacheMap[name] = db

			log.Debug().Msgf("found gene database %s %s", name, filepath.Join(dir, basename))
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
		//db := NewGeneDB(assembly, filepath.Join(cache.dir, fmt.Sprintf("%s.db", assembly)))
		db := NewGeneDB(assembly, filepath.Join(cache.dir, fmt.Sprintf("gtf_%s.db", assembly)))

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
	DB   *sql.DB
	File string
	//withinGeneAndPromStmt *sql.Stmt
	Name string
}

func NewGeneDB(name string, file string) *GeneDB {
	db := sys.Must(sql.Open("sqlite3", file))

	return &GeneDB{Name: name, File: file, DB: db} //withinGeneStmt: sys.Must(db.Prepare(WITHIN_GENE_SQL)),
	//withinGeneAndPromStmt: sys.Must(db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL))

	//inExonStmt:      sys.Must(db.Prepare(IN_EXON_SQL)),
	//closestGeneStmt: sys.Must(db.Prepare(CLOSEST_GENE_SQL))

}

func (genedb *GeneDB) Close() {
	genedb.DB.Close()
}

func (genedb *GeneDB) GeneDBInfo() (*GeneDBInfo, error) {
	var id uint
	var genome string
	var version string

	err := genedb.DB.QueryRow(GENE_DB_INFO_SQL).Scan(&id, &genome, &version)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return &GeneDBInfo{Name: genedb.Name, Genome: genome, Version: version}, nil
}

func (genedb *GeneDB) OverlappingGenes(location *dna.Location,
	canonicalMode bool,
	geneTypeFilter string) ([]*GenomicFeature, error) {

	var gid uint
	var feature Feature
	var chr string
	var start uint
	var end uint
	var strand string
	var geneName string
	var geneId string
	var geneType string

	var transcriptId string
	var transcriptName string
	var isCanonical bool

	var exonNumber uint

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	var features = make([]*GenomicFeature, 0, 10)
	var currentGene *GenomicFeature
	var currentTranscript *GenomicFeature
	var geneRows *sql.Rows
	var err error

	sql := OVERLAP_LOCATION_SQL

	//if canonical {
	//	sql += " AND (g.level = 1 OR g.is_canonical = 1)"
	//}

	if geneTypeFilter != "" {
		sql += " AND g.type = ?4" + OVERLAP_ORDER_BY_SQL

		geneRows, err = genedb.DB.Query(sql,
			location.Chr,
			location.Start,
			location.End,
			geneTypeFilter)
	} else {
		sql += OVERLAP_ORDER_BY_SQL

		geneRows, err = genedb.DB.Query(sql,
			location.Chr,
			location.Start,
			location.End)
	}

	//log.Debug().Msgf("overlapping genes %s", sql)

	if err != nil {
		//log.Debug().Msgf("error querying overlapping gene %s", err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer geneRows.Close()

	//var currentExon *GenomicSearchFeature

	for geneRows.Next() {
		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &exonId)
		err := geneRows.Scan(&gid,
			&feature,
			&chr,
			&start,
			&end,
			&strand,
			&geneId,
			&geneName,
			&geneType,
			&isCanonical,
			&transcriptId,
			&transcriptName,
			&exonNumber,
		)

		if err != nil {
			//log.Debug().Msgf("error reading overlapping gene rows %s", err)
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		//log.Debug().Msgf("overlap gene row %d %s %s %d-%d %s %s %s %s %t %s %s %d", gid, feature, chr, start, end, strand, geneId, geneName, geneType, isCanonical, transcriptId, transcriptName, exonNumber)

		location := dna.NewStrandedLocation(chr, start, end, strand)

		switch feature {
		case FEATURE_GENE:
			// only add a new gene if we don't already have it. We
			// assume the rows are ordered by gene id hence if the
			// id changes, we are processing a set of rows for a new gene
			if currentGene == nil || currentGene.GeneId != geneId {
				feature := &GenomicFeature{Id: gid,
					Location:   location,
					Feature:    FEATURE_GENE,
					GeneSymbol: geneName,
					GeneId:     geneId,
					//Strand:   strand,
					Type:     geneType,
					Children: make([]*GenomicFeature, 0, 10)}

				currentGene = feature

				features = append(features, currentGene)
			}
		case FEATURE_TRANSCRIPT:
			//if canonical mode only add if transcript is canonical
			// also only add if we have a current gene
			// also only add if we don't already have this transcript
			if currentGene != nil &&
				geneId == currentGene.GeneId &&
				(!canonicalMode || isCanonical) &&
				(currentTranscript == nil || currentTranscript.TranscriptId != transcriptId) {
				feature := &GenomicFeature{Id: gid,
					Location: location,
					//Strand:       strand,
					Feature:      FEATURE_TRANSCRIPT,
					GeneSymbol:   geneName,
					GeneId:       geneId,
					Type:         geneType,
					IsCanonical:  isCanonical,
					TranscriptId: transcriptId,
					Children:     make([]*GenomicFeature, 0, 10)}

				currentTranscript = feature

				currentGene.Children = append(currentGene.Children, currentTranscript)
			}
		case FEATURE_EXON:
			// only add exon if we have a current transcript and it matches
			// the transcript id
			if currentTranscript != nil &&
				currentTranscript.TranscriptId == transcriptId {
				feature := &GenomicFeature{Id: gid,
					Location:     location,
					Feature:      FEATURE_EXON,
					GeneSymbol:   geneName,
					GeneId:       geneId,
					TranscriptId: transcriptId,
					ExonNumber:   exonNumber,
					//Strand:       strand,
					Type: geneType}
				currentTranscript.Children = append(currentTranscript.Children, feature)
			}
		default:
			// do nothing
		}
	}

	//log.Debug().Msgf("overlap gene row here %v", len(features))

	return features, nil
}

func (genedb *GeneDB) SearchForGeneByName(search string,
	feature Feature,
	n uint16,
	fuzzy bool,
	canonical bool,
	geneType string) ([]*GenomicFeature, error) {
	n = max(1, min(n, MAX_GENE_INFO_RESULTS))

	// case insensitive search
	search = strings.ToLower(search)

	if len(search) < 2 || strings.Contains(search, "chr:") {
		return nil, fmt.Errorf("%s is an invalid search term", search)
	}

	//log.Debug().Msgf("search %s", search)

	var rows *sql.Rows
	var err error

	if fuzzy && !strings.HasSuffix(search, "%") {
		search += "%"
	}

	var sql string

	if feature == FEATURE_GENE {
		sql = GENE_INFO_SQL
	} else {
		sql = TRANSCRIPT_INFO_SQL
	}

	if canonical {
		sql += " AND t.is_canonical = 1"
	}

	//orderSql := " ORDER BY g.seqname, g.start, g.end, g.gene_name, g.transcript_id, g.exon_number"

	if geneType != "" {
		ids := strings.Split(strings.ToLower(geneType), ",")

		for i := range ids {
			ids[i] = strings.TrimSpace(ids[i])
		}

		placeholders := make([]string, len(ids))
		args := make([]any, len(ids)+2)

		args[0] = search
		args[1] = n

		for i, id := range ids {
			placeholders[i] = fmt.Sprintf("?%d", i+3)
			args[i+2] = id
		}

		sql += fmt.Sprintf(" AND gene_type IN (%s)", strings.Join(placeholders, ",")) + " LIMIT ?2"

		rows, err = genedb.DB.Query(sql, args...)
	} else {
		sql += " LIMIT ?2"

		rows, err = genedb.DB.Query(sql, search, n)
	}

	//log.Debug().Msgf("search %s  ", sql)

	if err != nil {
		//log.Debug().Msgf("err %s  ", err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	ret, err := rowsToRecords(rows)

	if err != nil {
		//log.Debug().Msgf("error reading rows %s", err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	// put in order of position
	SortFeaturesByPos(ret)

	//log.Debug().Msgf("search found %d features", len(ret))

	return ret, nil
}

func (genedb *GeneDB) WithinGenes(location *dna.Location, feature Feature) (*GenomicFeatures, error) {

	// rows, err := genedb.withinGeneStmt.Query(
	// 	mid,
	// 	level,
	// 	location.Chr,
	// 	location.Start,
	// 	location.Start,
	// 	location.End,
	// 	location.End)

	rows, err := genedb.DB.Query(WITHIN_GENE_SQL,
		feature,
		location.Chr,
		location.Mid(),
		location.Start,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, feature)
}

func (genedb *GeneDB) WithinGenesAndPromoter(location *dna.Location, pad5p uint, pad3p uint) (*GenomicFeatures, error) {

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

	//log.Debug().Msgf("within genes and promoter %s %d %d pad5p=%d pad3p=%d %s", location.Chr, location.Start, location.End, pad5p, pad3p, WITHIN_GENE_AND_PROMOTER_SQL)

	rows, err := genedb.DB.Query(WITHIN_GENE_AND_PROMOTER_SQL,
		location.Chr,
		location.Mid(),
		location.Start,
		location.End,
		pad5p,
		pad3p)

	if err != nil {
		//log.Debug().Msgf("error querying within gene and promoter %s", err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, FEATURE_TRANSCRIPT)
}

func (genedb *GeneDB) InExon(location *dna.Location, transcriptId string) (*GenomicFeatures, error) {

	// rows, err := genedb.inExonStmt.Query(
	// 	mid,
	// 	geneId,
	// 	location.Chr,
	// 	location.Start,
	// 	location.Start,
	// 	location.End,
	// 	location.End)

	rows, err := genedb.DB.Query(IN_EXON_SQL,
		transcriptId,
		location.Mid(),
		location.Start,
		location.End)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, FEATURE_EXON)
}

func (genedb *GeneDB) ClosestGenes(location *dna.Location, closestN uint16) ([]*GenomicFeature, error) {
	mid := (location.Start + location.End) / 2

	// rows, err := genedb.closestGeneStmt.Query(mid,
	// 	level,
	// 	location.Chr,
	// 	mid,

	// 	n)

	rows, err := genedb.DB.Query(CLOSEST_GENE_SQL, FEATURE_GENE, location.Chr, mid, closestN)

	if err != nil {
		//log.Debug().Msgf("error querying closest gene %s %s ", sql, err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	genes, err := rowsToRecords(rows)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	ids := make([]string, len(genes))
	for i, gene := range genes {
		ids[i] = gene.GeneId
	}

	placeholders := make([]string, len(ids))
	args := make([]any, len(ids)+1)

	args[0] = mid

	for i, id := range ids {
		placeholders[i] = fmt.Sprintf("?%d", i+2)
		args[i+1] = id
	}

	sql := STANDARD_GENE_FIELDS +
		`, ?1 - g.tss as tss_dist
		FROM gtf AS g
		WHERE g.feature = 'transcript'`

	sql += fmt.Sprintf(" AND g.gene_id IN (%s)", strings.Join(placeholders, ",")) + " ORDER BY g.gene_id, ABS(tss_dist)"

	//log.Debug().Msgf("closest transcripts %s", sql)

	rows, err = genedb.DB.Query(sql, args...)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	transcripts, err := rowsToRecords(rows)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	// find the closest transcript per gene

	closesMap := make(map[string]*GenomicFeature)

	for _, feature := range transcripts {
		existing, ok := closesMap[feature.GeneId]

		if !ok || basemath.AbsInt(feature.TssDist) < basemath.AbsInt(existing.TssDist) {
			closesMap[feature.GeneId] = feature
		}
	}

	closest := make([]*GenomicFeature, 0, closestN)

	for _, feature := range genes {

		transcript, ok := closesMap[feature.GeneId]

		if ok {
			closest = append(closest, transcript)
		}
	}

	return closest, nil
}

func rowsToFeatures(location *dna.Location, rows *sql.Rows, feature Feature) (*GenomicFeatures, error) {

	//log.Debug().Msgf("rowsToFeatures %s   %d %d", location.Chr, location.Start, location.End)

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	features, err := rowsToRecords(rows)

	if err != nil {
		return nil, err
	}

	ret := GenomicFeatures{Location: location, Feature: feature, Features: features}

	return &ret, nil
}

// func (genedb *GeneDB) IdToName(id uint) string {

// 	name := "n/a"

// 	// ignore error as it is not critical
// 	// if the id is not found, we return "n/a"
// 	genedb.db.QueryRow(ID_TO_NAME_SQL, id).Scan(&name)

// 	return name
// }

func rowsToRecords(rows *sql.Rows) ([]*GenomicFeature, error) {
	defer rows.Close()

	var id uint
	var feature Feature
	var chr string
	var start uint
	var end uint
	var strand string
	var geneId string
	var geneName string
	var geneType string
	var transcriptId string
	var transcriptName string
	var isCanonical bool
	var exonNumber uint

	var d int
	var err error

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	//var features = make([]*GenomicFeature, 0, 10)

	features := make([]*GenomicFeature, 0, 10)

	for rows.Next() {
		geneId = ""
		geneName = ""
		transcriptId = ""
		transcriptName = ""
		geneType = ""

		d = 0 // tss distance

		err = rows.Scan(&id,
			&feature,
			&chr,
			&start,
			&end,
			&strand,
			&geneId,
			&geneName,
			&geneType,
			&isCanonical,
			&transcriptId,
			&transcriptName,
			&exonNumber,
			&d)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		location := dna.NewLocation(chr, start, end)

		feature := GenomicFeature{
			Feature:  feature,
			Location: location,
			//Strand:       strand,
			GeneId:       geneId,
			GeneSymbol:   geneName,
			Type:         geneType,
			TranscriptId: transcriptId,
			IsCanonical:  isCanonical,
			TssDist:      d}

		features = append(features, &feature)
	}

	// enforce sorted correctly by chr and then position
	//sort.Sort(SortFeatureByPos(*features))

	return features, nil
}

// SortFeaturesByPos sorts features in place by chr, start, end
func SortFeaturesByPos(features []*GenomicFeature) {
	slices.SortFunc(features, func(a, b *GenomicFeature) int {
		ci := dna.ChromToInt(a.Location.Chr)
		cj := dna.ChromToInt(b.Location.Chr)

		// on different chrs so sort by chr
		if ci != cj {
			return int(ci) - int(cj)
		}

		// same chr so sort by position
		if a.Location.Start != b.Location.Start {
			return int(a.Location.Start) - int(b.Location.Start)
		}

		// same start so sort by end
		return int(a.Location.End) - int(b.Location.End)
	})
}

// type SortFeatureByPos []*GenomicFeature

// func (features SortFeatureByPos) Len() int      { return len(features) }
// func (features SortFeatureByPos) Swap(i, j int) { features[i], features[j] = features[j], features[i] }
// func (features SortFeatureByPos) Less(i, j int) bool {
// 	ci := dna.ChromToInt(features[i].Location.Chr)
// 	cj := dna.ChromToInt(features[j].Location.Chr)

// 	// on different chrs so sort by chr
// 	if ci != cj {
// 		return ci < cj
// 	}

// 	if features[i].Location.Start != features[j].Location.Start {
// 		return features[i].Location.Start < features[j].Location.Start
// 	}

// 	return features[i].Location.End < features[j].Location.End
// }
