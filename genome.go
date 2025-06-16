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

const GENE_INFO_SQL = `SELECT 
	g.id, 
	g.chr, 
	g.start, 
	g.end, 
	g.strand, 
	g.gene_id, 
	g.gene_symbol, 
	gt.name as gene_type, 
	0 as tss_dist
	FROM gene AS g
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
 	WHERE (g.gene_symbol LIKE ?1 OR g.gene_id LIKE ?1)`

const TRANSCRIPT_INFO_SQL = `SELECT 
	t.id, 
	g.chr, 
	t.start, 
	t.end, 
	t.strand, 
	g.gene_id, 
	g.gene_symbol, 
	g.gene_type, 
	t.transcript_id, 
	tt.name as transcript_type,
	0 as tss_dist
	FROM transcript AS t
	INNER JOIN gene AS g ON g.id = t.gene_id
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
 	WHERE (g.gene_symbol LIKE ?1 OR g.gene_id LIKE ?1 OR t.transcript_id LIKE ?1)`

// const GENE_INFO_FUZZY_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, transcript_id, gene_type, 0 as tss_dist
//  	FROM gene
//   	WHERE (gene_symbol LIKE ?2 OR g.gene_id LIKE ?2 OR g.transcript_id LIKE ?2)
//  	ORDER BY chr, start
// 	LIMIT ?3`

const WITHIN_GENE_SQL = `SELECT 
	g.id, 
	g.chr, 
	t.start, 
	t.end, 
	g.strand, 
	g.gene_id, 
	g.gene_symbol, 
	gt.name as gene_type, 
	t.transcript_id, 
	tt.name as transcript_type,
	?1 - g.tss as tss_dist
	FROM transcript AS t
	INNER JOIN gene AS g ON g.id = t.gene_id
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
 	WHERE g.chr= ?3 AND (t.start <= ?5 AND t.end >= ?4)
 	ORDER BY t.start`

const CLOSEST_GENE_SQL = `SELECT 
	g.id, 
	g.chr, 
	g.start, 
	g.end, 
	g.strand, 
	g.gene_symbol, 
	g.gene_id, 
	gt.name as gene_type, 
	?1 - g.tss as tss_dist
	FROM gene AS g
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
 	WHERE g.chr = ?1
 	ORDER BY ABS(g.tss - ?2)
 	LIMIT ?3`

const CLOSEST_TRANSCRIPT_SQL = `SELECT 
	t.id, 
	g.chr, 
	t.start, 
	t.end, 
	g.strand, 
	g.gene_id, 
	g.gene_symbol, 
	gt.name as gene_type, 
	t.transcript_id,
	t.is_canonical,
	tt.name as transcript_type,
	?1 - g.tss as tss_dist
	FROM transcript AS t
	INNER JOIN gene AS g ON g.id = t.gene_id
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
 	WHERE g.chr = ?1
 	ORDER BY ABS(g.tss - ?2)
 	LIMIT ?3`

const WITHIN_GENE_AND_PROMOTER_SQL = `SELECT 
	t.id, 
	g.chr, 
	t.start, 
	t.end, 
	g.strand, 
	g.gene_id, 
	g.gene_symbol, 
	gt.name as gene_type,
	t.transcript_id,
	t.is_canonical,
	tt.name as transcript_type,
	?1 - g.tss as tss_dist 
	FROM gene AS g
	INNER JOIN transcript AS t ON g.id = t.gene_id
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
 	WHERE g.chr = ?2 AND 
	(
		((t.start - ?5 <= ?4 AND t.end >= ?3) AND g.strand = '+') OR
		((t.start >= ?4 AND t.end + ?5 <= ?3) AND g.strand = '-')
	)
 	ORDER BY t.start, t.end DESC`

// when annotating genes, see if position falls within an exon
const IN_EXON_SQL = `SELECT
	e.id,
	g.chr,
	e.start,
	e.end,
	g.strand,
	g.gene_id,
	g.gene_symbol,
	gt.name as gene_type,
	t.transcript_id,
	t.is_canonical,
	tt.name as transcript_type,
	?1 - g.tss as tss_dist
	FROM exon AS e
	INNER JOIN transcript AS t ON t.id = e.transcript_id
	INNER JOIN gene AS g ON g.id = t.gene_id
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
 	WHERE g.gene_id = ?2 AND (e.start <= ?4 AND e.end >= ?3)
 	ORDER BY e.start, e.end DESC`

// If start is less x2 and end is greater than x1, it constrains the feature to be overlapping
// our location
// const OVERLAPPING_GENES_FROM_LOCATION_SQL = `SELECT
// 	g.id,
// 	g.chr,
// 	g.start,
// 	g.end,
// 	g.strand,
// 	g.gene_id,
// 	g.gene_symbol,
// 	gt.name as gene_type
// 	FROM gene AS g
// 	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
// 	WHERE g.chr = ?1 AND (g.start <= ?3 AND g.end >= ?2)`

// get all genes, transcripts, and exons overlapping a location
// which we can use to build a nested gene structure
const OVERLAP_LOCATION_SQL = `SELECT 
	g.id AS gid, 
	g.chr, 
	g.start AS gstart, 
	g.end AS gend, 
	g.strand AS gstrand, 
	g.gene_id, 
	g.gene_symbol, 
	gt.name AS gene_type,
	t.id AS tid,
	t.start AS tstart,
	t.end AS tend,
	t.transcript_id,
	t.is_canonical,
	tt.name AS transcript_type,
	e.id AS eid, 
	e.start AS estart, 
	e.end AS eend,  
	e.exon_id,
	e.exon_number
	FROM exon AS e
	INNER JOIN transcript AS t ON t.id = e.transcript_id
	INNER JOIN gene AS g ON g.id = t.gene_id
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
	WHERE g.chr = ?1 AND (g.start <= ?3 AND g.end >= ?2)`

const OVERLAP_ORDER_BY_SQL = ` ORDER BY 
	g.chr, 
	g.start, 
	g.end DESC,
	g.gene_id, 
	t.start, 
	t.end DESC,
	t.transcript_id, 
	e.start, 
	e.end DESC`

const TRANSCRIPTS_IN_GENE_SQL = `SELECT 
	t.id, 
	g.chr, 
	t.start, 
	t.end, 
	g.strand,  
	gt.name as gene_type,
	t.transcript_id, 
	t.is_canonical,
	tt.name as transcript_type
	FROM transcript AS t
	INNER JOIN gene AS g ON g.id = t.gene_id
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
	WHERE t.gene_id = ?1`

// const CANONICAL_TRANSCRIPTS_IN_GENE_SQL = `SELECT id, level, chr, start, end, strand, transcript_id, is_canonical, gene_type
// 	FROM gene
// 	WHERE level = 2 AND gene_id = ?1 AND is_canonical = 1
// 	ORDER BY start`

const EXONS_IN_TRANSCRIPT_SQL = `SELECT 
	e.id, 
	g.chr, 
	e.start, 
	e.end, 
	g.strand, 
	gt.name as gene_type,
	e.exon_id, 
	t.is_canonical, 
	tt.name as transcript_type
	FROM exon AS e
	INNER JOIN transcript AS t ON t.id = e.transcript_id
	INNER JOIN gene AS g ON g.id = t.gene_id
	INNER JOIN gene_type AS gt ON g.gene_type_id = gt.id
	INNER JOIN transcript_type AS tt ON t.transcript_type_id = tt.id
	WHERE e.transcript_id = ?1`

const ID_TO_NAME_SQL = `SELECT name FROM ids WHERE id = ?1`

const MAX_GENE_INFO_RESULTS uint16 = 100

// const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_symbol, transcript_id, start - ?
// 	FROM gene
//  	WHERE level = 2 AND gene_id = ? AND chr = ? AND ? >= stranded_start - ? AND ? <= stranded_start + ?
//  	ORDER BY start ASC`

type GenomicFeature struct {
	Location       *dna.Location     `json:"loc"`
	GeneType       string            `json:"geneType,omitempty"`
	Strand         string            `json:"strand"`
	Level          Level             `json:"level"`
	GeneSymbol     string            `json:"geneSymbol,omitempty"`
	GeneId         string            `json:"geneId,omitempty"`
	TranscriptId   string            `json:"transcriptId,omitempty"`
	TranscriptType string            `json:"transcriptType,omitempty"`
	ExonId         string            `json:"exonId,omitempty"`
	ExonNumber     uint              `json:"exonNumber,omitempty"`
	PromLabel      string            `json:"promLabel,omitempty"`
	Children       []*GenomicFeature `json:"children,omitempty"`
	Id             uint              `json:"-"`
	TssDist        int               `json:"tssDist,omitempty"`
	IsCanonical    bool              `json:"isCanonical"`
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
// func (genedb *GeneDB) OverlappingGenes(location *dna.Location, canonical bool, geneTypeFilter string) ([]*GenomicFeature, error) {

// 	var id uint

// 	var chr string
// 	var start uint
// 	var end uint
// 	var strand string
// 	var geneSymbol string
// 	var geneId string
// 	var transcriptId string
// 	var exonId string
// 	var isCanonical bool
// 	var geneType string

// 	// 10 seems a reasonable guess for the number of features we might see, just
// 	// to reduce slice reallocation
// 	var features = make([]*GenomicFeature, 0, 10)
// 	var currentGene *GenomicFeature
// 	var currentTranscript *GenomicFeature
// 	var geneRows *sql.Rows
// 	var err error

// 	sql := OVERLAPPING_GENES_FROM_LOCATION_SQL

// 	if geneTypeFilter != "" {
// 		sql += " AND gene_type = ?4 ORDER BY g.chr, g.start, g.end DESC"

// 		geneRows, err = genedb.db.Query(sql,
// 			location.Chr,
// 			location.Start,
// 			location.End,
// 			geneTypeFilter)
// 	} else {
// 		sql += " ORDER BY g.chr, g.start, g.end DESC"
// 		geneRows, err = genedb.db.Query(sql,
// 			location.Chr,
// 			location.Start,
// 			location.End)
// 	}

// 	if err != nil {
// 		log.Debug().Msgf("error querying overlapping gene %s", err)
// 		return nil, err //fmt.Errorf("there was an error with the database query")
// 	}

// 	defer geneRows.Close()

// 	//var currentExon *GenomicSearchFeature

// 	for geneRows.Next() {
// 		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &exonId)
// 		err := geneRows.Scan(&id, &chr, &start, &end, &strand, &geneId, &geneSymbol, &geneType)

// 		if err != nil {
// 			return nil, err //fmt.Errorf("there was an error with the database records")
// 		}

// 		location := &dna.Location{Chr: chr, Start: start, End: end}

// 		feature := &GenomicFeature{Id: id,
// 			Location:   location,
// 			Level:      LEVEL_GENE,
// 			GeneSymbol: geneSymbol,
// 			GeneId:     geneId,
// 			Strand:     strand,
// 			GeneType:   geneType,
// 			Children:   make([]*GenomicFeature, 0, 10)}

// 		currentGene = feature

// 		sql = TRANSCRIPTS_IN_GENE_SQL

// 		if canonical {
// 			sql += " AND t.is_canonical = 1"
// 		}

// 		sql += " ORDER BY t.start, t.end DESC"

// 		transcriptRows, err := genedb.db.Query(sql,
// 			currentGene.Id)

// 		if err != nil {
// 			return nil, err //fmt.Errorf("there was an error with the database query")
// 		}

// 		defer transcriptRows.Close()

// 		//var currentExon *GenomicSearchFeature

// 		for transcriptRows.Next() {
// 			err := transcriptRows.Scan(&id, &chr, &start, &end, &strand, &transcriptId, &isCanonical, &geneType)

// 			if err != nil {
// 				log.Debug().Msgf("error scanning transcript rows %s", err)
// 				return nil, err //fmt.Errorf("there was an error with the database records")
// 			}

// 			location := &dna.Location{Chr: chr, Start: start, End: end}

// 			feature := &GenomicFeature{Id: id,
// 				Location:     location,
// 				Level:        "transcript",
// 				GeneSymbol:   geneSymbol,
// 				GeneId:       geneId,
// 				TranscriptId: transcriptId,
// 				Strand:       strand,
// 				IsCanonical:  isCanonical,
// 				GeneType:     geneType,
// 				Children:     make([]*GenomicFeature, 0, 10)}

// 			currentTranscript = feature

// 			sql = EXONS_IN_TRANSCRIPT_SQL

// 			if canonical {
// 				sql += " AND t.is_canonical = 1"
// 			}

// 			sql += " ORDER BY e.start, e.end DESC"

// 			exonRows, err := genedb.db.Query(sql,
// 				currentTranscript.Id)

// 			if err != nil {
// 				return nil, err //fmt.Errorf("there was an error with the database query")
// 			}

// 			defer exonRows.Close()

// 			//var currentExon *GenomicSearchFeature

// 			for exonRows.Next() {
// 				err := exonRows.Scan(&id, &chr, &start, &end, &strand, &exonId, &isCanonical, &geneType)

// 				if err != nil {
// 					return nil, err //fmt.Errorf("there was an error with the database records")
// 				}

// 				location := &dna.Location{Chr: chr, Start: start, End: end}
// 				feature := &GenomicFeature{Id: id,
// 					Location:     location,
// 					Level:        "exon",
// 					GeneSymbol:   geneSymbol,
// 					GeneId:       geneId,
// 					TranscriptId: transcriptId,
// 					ExonId:       exonId,
// 					Strand:       strand,
// 					IsCanonical:  isCanonical,
// 					GeneType:     geneType}
// 				currentTranscript.Children = append(currentTranscript.Children, feature)
// 			}

// 			if len(currentTranscript.Children) > 0 {
// 				currentGene.Children = append(currentGene.Children, currentTranscript)
// 			}
// 		}

// 		if len(currentGene.Children) > 0 {
// 			features = append(features, currentGene)
// 		}
// 	}

// 	return features, nil
// }

func (genedb *GeneDB) OverlappingGenes(location *dna.Location, canonical bool, geneTypeFilter string) ([]*GenomicFeature, error) {

	var gid uint
	var chr string
	var gstart uint
	var gend uint
	var strand string
	var geneSymbol string
	var geneId string
	var geneType string
	var tid uint
	var tstart uint
	var tend uint
	var transcriptId string
	var isCanonical bool
	var transcriptType string
	var eid uint
	var estart uint
	var eend uint
	var exonId string
	var exonNumber uint

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	var features = make([]*GenomicFeature, 0, 10)
	var currentGene *GenomicFeature
	var currentTranscript *GenomicFeature
	var geneRows *sql.Rows
	var err error

	sql := OVERLAP_LOCATION_SQL

	if geneTypeFilter != "" {
		sql += " AND gene_type = ?4" + OVERLAP_ORDER_BY_SQL

		geneRows, err = genedb.db.Query(sql,
			location.Chr,
			location.Start,
			location.End,
			geneTypeFilter)
	} else {
		sql += OVERLAP_ORDER_BY_SQL
		geneRows, err = genedb.db.Query(sql,
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
			&chr,
			&gstart,
			&gend,
			&strand,
			&geneId,
			&geneSymbol,
			&geneType,
			&tid,
			&tstart,
			&tend,
			&transcriptId,
			&isCanonical,
			&transcriptType,
			&eid,
			&estart,
			&eend,
			&exonId,
			&exonNumber,
		)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		if currentGene == nil || currentGene.Id != gid {

			feature := &GenomicFeature{Id: gid,
				Location:   &dna.Location{Chr: chr, Start: gstart, End: gend},
				Level:      LEVEL_GENE,
				GeneSymbol: geneSymbol,
				GeneId:     geneId,
				Strand:     strand,
				GeneType:   geneType,
				Children:   make([]*GenomicFeature, 0, 10)}

			currentGene = feature

			features = append(features, currentGene)
		}

		if currentTranscript == nil || currentTranscript.TranscriptId != transcriptId {

			feature := &GenomicFeature{Id: gid,
				Location:     &dna.Location{Chr: chr, Start: tstart, End: tend},
				Level:        "transcript",
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				Strand:       strand,
				IsCanonical:  isCanonical,
				GeneType:     geneType,
				Children:     make([]*GenomicFeature, 0, 10)}

			currentTranscript = feature

			currentGene.Children = append(currentGene.Children, currentTranscript)
		}

		location := &dna.Location{Chr: chr, Start: estart, End: eend}
		feature := &GenomicFeature{Id: eid,
			Location:     location,
			Level:        "exon",
			GeneSymbol:   geneSymbol,
			GeneId:       geneId,
			TranscriptId: transcriptId,
			ExonId:       exonId,
			ExonNumber:   exonNumber,
			Strand:       strand,
			IsCanonical:  isCanonical,
			GeneType:     geneType}
		currentTranscript.Children = append(currentTranscript.Children, feature)
	}

	return features, nil
}

func (genedb *GeneDB) SearchForGeneByName(search string,
	level Level,
	n uint16,
	fuzzy bool,
	canonical bool,
	geneType string) ([]*GenomicFeature, error) {
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

	var sql string

	if level == LEVEL_GENE {
		sql = GENE_INFO_SQL
	} else {
		sql = TRANSCRIPT_INFO_SQL
	}

	if level != LEVEL_GENE && canonical {
		sql += " AND t.is_canonical = 1"
	}

	var orderSql string

	if level == LEVEL_GENE {
		orderSql = " ORDER BY g.chr, g.start, g.end DESC"
	} else {
		orderSql = " ORDER BY g.chr, t.start, t.end DESC"
	}

	if geneType != "" {
		sql += " AND gene_type = ?2" + orderSql + " LIMIT ?3"

		rows, err = genedb.db.Query(sql,
			search,
			geneType,
			n)
	} else {
		sql += orderSql + " LIMIT ?2"

		rows, err = genedb.db.Query(sql,
			search,
			n)
	}

	if err != nil {
		//log.Debug().Msgf("search %s  ", err)
		return ret, err //fmt.Errorf("there was an error with the database query")
	}

	err = rowsToRecords(rows, level, &ret)

	if err != nil {
		//log.Debug().Msgf("error reading rows %s", err)
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

func (genedb *GeneDB) WithinGenesAndPromoter(location *dna.Location, pad uint) (*GenomicFeatures, error) {
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
		location.Chr,
		location.Start,
		location.End,
		pad)

	if err != nil {
		//log.Debug().Msgf("error querying within gene and promoter %s", err)
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

	rows, err := genedb.db.Query(sql, location.Chr, mid, n)

	if err != nil {
		//log.Debug().Msgf("error querying closest gene %s %s ", sql, err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return rowsToFeatures(location, rows, level)
}

func rowsToFeatures(location *dna.Location, rows *sql.Rows, level Level) (*GenomicFeatures, error) {
	ret := GenomicFeatures{Location: location, Level: level, Features: make([]*GenomicFeature, 0, 10)}

	//log.Debug().Msgf("rowsToFeatures %s   %d %d", location.Chr, location.Start, location.End)

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
	var geneType string
	var transcriptId string
	var isCanonical bool
	var transcriptType string
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
		transcriptType = ""
		d = 0 // tss distance

		switch level {
		case LEVEL_TRANSCRIPT, LEVEL_EXON:
			err = rows.Scan(&id,
				&chr,
				&start,
				&end,
				&strand,
				&geneId,
				&geneSymbol,
				&geneType,
				&transcriptId,
				&isCanonical,
				&transcriptType,
				&d)
		default:
			err = rows.Scan(&id, &chr, &start, &end, &strand, &geneId, &geneSymbol, &geneType, &d)
		}

		if err != nil {
			return err //fmt.Errorf("there was an error with the database records")
		}

		location := dna.NewLocation(chr, start, end)

		feature := GenomicFeature{
			Level:          level,
			Location:       location,
			Strand:         strand,
			GeneId:         geneId,
			GeneSymbol:     geneSymbol,
			GeneType:       geneType,
			TranscriptId:   transcriptId,
			IsCanonical:    isCanonical,
			TranscriptType: transcriptType,
			TssDist:        d}

		*features = append(*features, &feature)
	}

	// enforce sorted correctly by chr and then position
	//sort.Sort(SortFeatureByPos(*features))

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

	if features[i].Location.Start != features[j].Location.Start {
		return features[i].Location.Start < features[j].Location.Start
	}

	return features[i].Location.End < features[j].Location.End
}
