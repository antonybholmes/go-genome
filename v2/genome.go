package v2

import (
	"database/sql"
	"fmt"
	"path/filepath"
	"strings"

	"github.com/antonybholmes/go-basemath"
	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-genome"
	"github.com/antonybholmes/go-sys"
)

type (
	V2GeneDB struct {
		db   *sql.DB
		file string
		//withinGeneAndPromStmt *sql.Stmt
		name string
	}
)

const (
	StandardGeneFieldsSql = `SELECT DISTINCT
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

	GeneInfoSql = StandardGeneFieldsSql +
		`, 0 as tss_dist
		FROM gtf AS g
 		WHERE g.feature = 'gene' AND (g.gene_name LIKE ?1 OR g.gene_id LIKE ?1)`

	TranscriptInfoSql = StandardGeneFieldsSql +
		` FROM gtf AS g
 		WHERE g.feature = 'transcript' AND (g.gene_name LIKE ?1 OR g.gene_id LIKE ?1 OR g.transcript_id LIKE ?1)`

		// const GENE_INFO_FUZZY_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_name, transcript_id, gene_type, 0 as tss_dist
		//  	FROM gene
		//   	WHERE (gene_name LIKE ?2 OR g.gene_id LIKE ?2 OR g.transcript_id LIKE ?2)
		//  	ORDER BY chr, start
		// 	LIMIT ?3`

	WithinGeneSql = StandardGeneFieldsSql +
		`, ?3 - g.tss as tss_dist
		FROM gtf AS g
		WHERE g.feature = ?1 AND g.seqname = ?2 AND (g.start <= ?5 AND g.end >= ?4)
		ORDER BY g.start`

	ClosestGeneSql = StandardGeneFieldsSql +
		`, ?3 - g.tss as tss_dist
		FROM gtf AS g
		WHERE g.feature = ?1 AND g.seqname = ?2
		ORDER BY ABS(tss_dist)
		LIMIT ?4`

	TranscriptsPerGene = StandardGeneFieldsSql +
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

	WithinGeneAndPromoterSql = StandardGeneFieldsSql +
		`, ?2 - g.tss as tss_dist 
		FROM gtf AS g
		WHERE g.feature = 'transcript' AND g.seqname = ?1 AND 
		(
			((g.start - ?5 <= ?4 AND g.end + ?6 >= ?3) AND g.strand = '+') OR
			((g.start - ?6 >= ?4 AND g.end + ?5 <= ?3) AND g.strand = '-')
		)
		ORDER BY tss_dist`

		// when annotating genes, see if position falls within an exon
	InExonSql = StandardGeneFieldsSql +
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
	OverlapLocationSql = StandardGeneFieldsSql +
		` FROM gtf AS g
		WHERE g.seqname = ?1 AND (g.start <= ?3 AND g.end >= ?2)`

	OverlapOrderBySql = ` ORDER BY 
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

	IdToNameSql = `SELECT name FROM ids WHERE id = ?1`

	// const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_name, transcript_id, start - ?
	// 	FROM gene
	//  	WHERE level = 2 AND gene_id = ? AND chr = ? AND ? >= stranded_start - ? AND ? <= stranded_start + ?
	//  	ORDER BY start ASC`

)

func NewGeneDB(assembly string, dir string) genome.GeneDB {
	file := filepath.Join(dir, fmt.Sprintf("gtf_%s.db", assembly))

	db := sys.Must(sql.Open(sys.Sqlite3DB, file))

	return &V2GeneDB{name: assembly, file: file, db: db} //withinGeneStmt: sys.Must(db.Prepare(WITHIN_GENE_SQL)),
	//withinGeneAndPromStmt: sys.Must(db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL))

	//inExonStmt:      sys.Must(db.Prepare(IN_EXON_SQL)),
	//closestGeneStmt: sys.Must(db.Prepare(CLOSEST_GENE_SQL))

}

func (genedb *V2GeneDB) Close() error {
	return genedb.db.Close()
}

func (genedb *V2GeneDB) GeneDBInfo() (*genome.GeneDBInfo, error) {
	var id int
	var g string
	var version string

	err := genedb.db.QueryRow(genome.GeneDBInfoSql).Scan(&id, &g, &version)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return &genome.GeneDBInfo{Name: genedb.name, Genome: g, Version: version}, nil
}

func (genedb *V2GeneDB) OverlappingGenes(location *dna.Location,
	level string,
	prom *dna.PromoterRegion,
	canonicalMode bool,
	geneTypeFilter string) ([]*genome.GenomicFeature, error) {

	var gid int
	var feature string
	var chr string
	var start int
	var end int
	var strand string
	var geneName string
	var geneId string
	var geneType string

	var transcriptId string
	var transcriptName string
	var isCanonical bool

	var exonNumber int

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	var features = make([]*genome.GenomicFeature, 0, 10)
	var currentGene *genome.GenomicFeature
	var currentTranscript *genome.GenomicFeature
	var geneRows *sql.Rows
	var err error

	sql := OverlapLocationSql

	//if canonical {
	//	sql += " AND (g.level = 1 OR g.is_canonical = 1)"
	//}

	if geneTypeFilter != "" {
		sql += " AND g.type = ?4" + OverlapOrderBySql

		geneRows, err = genedb.db.Query(sql,
			location.Chr(),
			location.Start(),
			location.End(),
			geneTypeFilter)
	} else {
		sql += OverlapOrderBySql

		geneRows, err = genedb.db.Query(sql,
			location.Chr(),
			location.Start(),
			location.End())
	}

	if err != nil {
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
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		location, err := dna.NewStrandedLocation(chr, start, end, strand)

		if err != nil {
			return nil, err
		}

		switch feature {
		case genome.GeneLevel:
			// only add a new gene if we don't already have it. We
			// assume the rows are ordered by gene id hence if the
			// id changes, we are processing a set of rows for a new gene
			if currentGene == nil || currentGene.GeneId != geneId {
				feature := &genome.GenomicFeature{Id: gid,
					Location:   location,
					Feature:    genome.GeneLevel,
					GeneSymbol: geneName,
					GeneId:     geneId,
					//Strand:   strand,
					Type:     geneType,
					Children: make([]*genome.GenomicFeature, 0, 10)}

				currentGene = feature

				features = append(features, currentGene)
			}
		case genome.TranscriptLevel:
			//if canonical mode only add if transcript is canonical
			// also only add if we have a current gene
			// also only add if we don't already have this transcript
			if currentGene != nil &&
				geneId == currentGene.GeneId &&
				(!canonicalMode || isCanonical) &&
				(currentTranscript == nil || currentTranscript.TranscriptId != transcriptId) {
				feature := &genome.GenomicFeature{Id: gid,
					Location: location,
					//Strand:       strand,
					Feature:    genome.TranscriptLevel,
					GeneSymbol: geneName,
					GeneId:     geneId,
					Type:       geneType,
					//IsCanonical:  &isCanonical,
					TranscriptId: transcriptId,
					Children:     make([]*genome.GenomicFeature, 0, 10)}

				currentTranscript = feature

				currentGene.Children = append(currentGene.Children, currentTranscript)
			}
		case genome.ExonLevel:
			// only add exon if we have a current transcript and it matches
			// the transcript id
			if currentTranscript != nil &&
				currentTranscript.TranscriptId == transcriptId {
				feature := &genome.GenomicFeature{Id: gid,
					Location:     location,
					Feature:      genome.ExonLevel,
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

	return features, nil
}

func (genedb *V2GeneDB) SearchForGeneByName(search string,
	feature string,
	fuzzy bool,
	canonical bool,
	geneType string,
	n int16) ([]*genome.GenomicFeature, error) {
	n = max(1, min(n, genome.MaxGeneInfoResults))

	// case insensitive search
	search = strings.ToLower(search)

	if len(search) < 2 || strings.Contains(search, "chr:") {
		return nil, fmt.Errorf("%s is an invalid search term", search)
	}

	var rows *sql.Rows
	var err error

	if fuzzy && !strings.HasSuffix(search, "%") {
		search += "%"
	}

	var sql string

	if feature == genome.GeneLevel {
		sql = GeneInfoSql
	} else {
		sql = TranscriptInfoSql
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

		rows, err = genedb.db.Query(sql, args...)
	} else {
		sql += " LIMIT ?2"

		rows, err = genedb.db.Query(sql, search, n)
	}

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	ret, err := genome.RowsToRecords(rows)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	// put in order of position
	genome.SortFeaturesByPos(ret)

	return ret, nil
}

func (genedb *V2GeneDB) WithinGenes(location *dna.Location, feature string, prom *dna.PromoterRegion) (*genome.GenomicFeatures, error) {

	// rows, err := genedb.withinGeneStmt.Query(
	// 	mid,
	// 	level,
	// 	location.Chr(),
	// 	location.Start(),
	// 	location.Start(),
	// 	location.End(),
	// 	location.End())

	rows, err := genedb.db.Query(WithinGeneSql,
		feature,
		location.Chr(),
		location.Mid(),
		location.Start(),
		location.End())

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return genome.RowsToFeatures(location, feature, rows)
}

func (genedb *V2GeneDB) WithinGenesAndPromoter(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion) ([]*genome.GenomicFeature, error) {

	//level := genome.MaxLevel(levels)

	// rows, err := genedb.withinGeneAndPromStmt.Query(
	// 	mid,
	// 	level,
	// 	location.Chr(),
	// 	pad,
	// 	location.Start(),
	// 	pad,
	// 	location.Start(),
	// 	pad,
	// 	location.End(),
	// 	pad,
	// 	location.End())

	rows, err := genedb.db.Query(WithinGeneAndPromoterSql,
		location.Chr(),
		location.Mid(),
		location.Start(),
		location.End(),
		prom.Upstream,
		prom.Downstream)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return genome.RowsToRecords(rows)
}

func (genedb *V2GeneDB) InExon(location *dna.Location,
	transcriptId string,
	prom *dna.PromoterRegion) ([]*genome.GenomicFeature, error) {

	// rows, err := genedb.inExonStmt.Query(
	// 	mid,
	// 	geneId,
	// 	location.Chr(),
	// 	location.Start(),
	// 	location.Start(),
	// 	location.End(),
	// 	location.End())

	rows, err := genedb.db.Query(InExonSql,
		transcriptId,
		location.Mid(),
		location.Start(),
		location.End())

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return genome.RowsToRecords(rows)
}

func (genedb *V2GeneDB) ClosestGenes(location *dna.Location, prom *dna.PromoterRegion, closestN int8) ([]*genome.GenomicFeature, error) {
	mid := (location.Start() + location.End()) / 2

	// rows, err := genedb.closestGeneStmt.Query(mid,
	// 	level,
	// 	location.Chr(),
	// 	mid,

	// 	n)

	rows, err := genedb.db.Query(ClosestGeneSql, genome.GeneLevel, location.Chr(), mid, closestN)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	genes, err := genome.RowsToRecords(rows)

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

	sql := StandardGeneFieldsSql +
		`, ?1 - g.tss as tss_dist
		FROM gtf AS g
		WHERE g.feature = 'transcript'`

	sql += fmt.Sprintf(" AND g.gene_id IN (%s)", strings.Join(placeholders, ",")) + " ORDER BY g.gene_id, ABS(tss_dist)"

	rows, err = genedb.db.Query(sql, args...)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	transcripts, err := genome.RowsToRecords(rows)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	// find the closest transcript per gene

	closesMap := make(map[string]*genome.GenomicFeature)

	for _, feature := range transcripts {
		existing, ok := closesMap[feature.GeneId]

		if !ok || basemath.AbsInt(feature.TssDist) < basemath.AbsInt(existing.TssDist) {
			closesMap[feature.GeneId] = feature
		}
	}

	closest := make([]*genome.GenomicFeature, 0, closestN)

	for _, feature := range genes {

		transcript, ok := closesMap[feature.GeneId]

		if ok {
			closest = append(closest, transcript)
		}
	}

	return closest, nil
}

// func (genedb *GeneDB) IdToName(id int) string {

// 	name := "n/a"

// 	// ignore error as it is not critical
// 	// if the id is not found, we return "n/a"
// 	genedb.db.QueryRow(ID_TO_NAME_SQL, id).Scan(&name)

// 	return name
// }
