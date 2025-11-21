package v1

import (
	"database/sql"
	"fmt"
	"path/filepath"
	"strings"

	"github.com/antonybholmes/go-basemath"
	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-genome"
	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-sys/log"
)

type (
	V1GeneDB struct {
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

	// WithinGeneSql = StandardGeneFieldsSql +
	// 	`, ?3 - g.tss as tss_dist
	// 	FROM gtf AS g
	// 	WHERE g.feature = ?1 AND g.seqname = ?2 AND (g.start <= ?5 AND g.end >= ?4)
	// 	ORDER BY g.start`

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

	// InGeneAndPromoterSql = `SELECT DISTINCT
	// 	g.id,
	// 	g.chr,
	// 	g.start,
	// 	g.end,
	// 	g.strand,
	// 	g.gene_id,
	// 	g.gene_symbol,
	// 	gt.name as gene_type,
	// 	CASE
	// 		WHEN strand = '-' THEN ?2 - g.end
	// 		ELSE ?2 - g.start
	// 	END AS tss_dist
	// 	FROM genes as g
	// 	JOIN gene_types AS gt ON g.gene_type_id = gt.id
	// 	WHERE g.chr = ?1 AND
	// 	(
	// 		((g.start - ?5 <= ?4 AND g.end + ?6 >= ?3) AND g.strand = '+') OR
	// 		((g.start - ?6 >= ?4 AND g.end + ?5 <= ?3) AND g.strand = '-')
	// 	)
	// 	ORDER BY tss_dist`

	// search within transcripts and their promoters

	CoreGeneSql = `SELECT DISTINCT
		g.id, 
		g.chr,
		g.start, 
		g.end, 
		g.strand, 
		g.gene_id, 
		g.gene_symbol,
		gt.name AS gene_type,
		t.transcript_id,
		t.start,
		t.end,
		t.is_canonical,
		t.is_longest,
		tt.name AS transcript_type,
		e.exon_id,
		e.start,
		e.end,
		e.exon_number,
		CASE
			WHEN g.strand = '-' THEN :mid - t.end
			ELSE :mid - t.start
    	END AS tss_dist,
		((:start <= t.start + :prom3p) AND (:end >= t.start - :prom5p) AND (g.strand = '+')) OR
			((:start <= t.end + :prom5p) AND (:end >= t.end - :prom3p) AND (g.strand = '-')) 
		AS in_promoter,
		:start <= e.end AND :end >= e.start AS in_exon,
		:start <= t.end AND :end >= t.start AS is_intragenic
		FROM genes as g
		JOIN transcripts AS t ON g.id = t.gene_id
		JOIN exons AS e ON e.transcript_id = t.id
		JOIN gene_types AS gt ON g.gene_type_id = gt.id
		JOIN transcript_types AS tt ON t.transcript_type_id = tt.id`

	InGeneAndPromoterSql = CoreGeneSql +
		` WHERE g.chr = :chr AND 
		(
			((:start <= t.end) AND (:end >= t.start - :prom5p) AND g.strand = '+') OR
			((:start <= t.end + :prom5p) AND (:end >= t.start) AND g.strand = '-')
		)
		ORDER BY g.gene_id, t.transcript_id, e.exon_number`

	InGeneSql = CoreGeneSql +
		` WHERE g.chr = :chr AND (:start <= t.end AND :end >= t.start)`

	ClosestGeneSql = CoreGeneSql +
		` WHERE 
			t.transcript_id IN (
				SELECT DISTINCT t.transcript_id
				FROM transcripts AS t
				WHERE 
					t.chr = :chr AND
					t.is_longest = 1
				ORDER BY ABS(
					CASE 
						WHEN t.strand ='-' THEN :mid - t.end
						ELSE :mid - t.start
					END
				)
				LIMIT :n			
			) AND
			t.is_longest = 1
		ORDER BY ABS(tss_dist)`

	// when annotating genes, see if position falls within an exon
	InExonSql = CoreGeneSql +
		` WHERE e.transcript_id = :transcriptId AND (e.start <= :end AND e.end >= :start)
		ORDER BY e.start, e.end DESC`

	// OverlapSql = `SELECT DISTINCT
	// 	g.id,
	// 	g.chr,
	// 	g.start,
	// 	g.end,
	// 	g.strand,
	// 	g.gene_id,
	// 	g.gene_symbol,
	// 	gt.name as gene_type,
	// 	t.transcript_id,
	// 	t.start,
	// 	t.end,
	// 	t.is_canonical,
	// 	t.is_longest,
	// 	tt.name as transcript_type,
	// 	e.exon_id,
	// 	e.start,
	// 	e.end,
	// 	e.exon_number,
	// 	CASE
	// 		WHEN strand = '-' THEN ?2 - t.end
	// 		ELSE ?2 - t.start
	// 	END AS tss_dist
	// 	FROM genes as g
	// 	JOIN transcripts AS t ON g.id = t.gene_id
	// 	JOIN exons AS e ON e.transcript_id = t.id
	// 	JOIN gene_types AS gt ON g.gene_type_id = gt.id
	// 	JOIN transcript_types AS tt ON t.transcript_type_id = tt.id
	// 	WHERE g.chr = ?1 AND (g.start <= ?3 AND g.end >= ?2)`

	OverlapSql = CoreGeneSql +
		` WHERE t.chr = :chr AND (t.start <= :end AND t.end >= :start)`

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
	// OverlapLocationSql = StandardGeneFieldsSql +
	// 	` FROM gtf AS g
	// 	WHERE g.seqname = ?1 AND (g.start <= ?3 AND g.end >= ?2)`

	OverlapOrderBySql = ` ORDER BY 
		g.gene_id,
		t.transcript_id,
		e.exon_number`

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
)

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func NewGeneDB(assembly string, dir string) genome.GeneDB {
	file := filepath.Join(dir, fmt.Sprintf("%s.db", assembly))

	db := sys.Must(sql.Open(sys.Sqlite3DB, file))

	return &V1GeneDB{name: assembly, file: file, db: db} //withinGeneStmt: sys.Must(db.Prepare(WITHIN_GENE_SQL)),
	//withinGeneAndPromStmt: sys.Must(db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL))

	//inExonStmt:      sys.Must(db.Prepare(IN_EXON_SQL)),
	//closestGeneStmt: sys.Must(db.Prepare(CLOSEST_GENE_SQL))

}

func (genedb *V1GeneDB) Close() error {
	return genedb.db.Close()
}

func (genedb *V1GeneDB) GeneDBInfo() (*genome.GeneDBInfo, error) {
	var id int
	var g string
	var version string

	err := genedb.db.QueryRow(genome.GeneDBInfoSql).Scan(&id, &g, &version)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return &genome.GeneDBInfo{Name: genedb.name, Genome: g, Version: version}, nil
}

func (genedb *V1GeneDB) OverlappingGenes(location *dna.Location,
	level string,
	prom *dna.PromoterRegion,
	canonicalMode bool,
	geneTypeFilter string) ([]*genome.GenomicFeature, error) {

	var geneRows *sql.Rows
	var err error

	query := OverlapSql

	log.Debug().Msgf("canonical mode %v gene type filter %s", canonicalMode, geneTypeFilter)

	//if canonical {
	//	sql += " AND (g.level = 1 OR g.is_canonical = 1)"
	//}

	if geneTypeFilter != "" {
		query += " AND gene_type = :geneType" + OverlapOrderBySql

		geneRows, err = genedb.db.Query(query,
			sql.Named("chr", location.Chr()),
			sql.Named("start", location.Start()),
			sql.Named("end", location.End()),
			sql.Named("mid", location.Mid()),
			sql.Named("prom5p", prom.Upstream()),
			sql.Named("prom3p", prom.Downstream()),
			sql.Named("geneType", geneTypeFilter))
	} else {
		query += OverlapOrderBySql

		geneRows, err = genedb.db.Query(query,
			sql.Named("chr", location.Chr()),
			sql.Named("start", location.Start()),
			sql.Named("end", location.End()),
			sql.Named("mid", location.Mid()),
			sql.Named("prom5p", prom.Upstream()),
			sql.Named("prom3p", prom.Downstream()))
	}

	if err != nil {
		log.Error().Msgf("error querying overlapping gene %s", err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer geneRows.Close()

	//var currentExon *GenomicSearchFeature

	// g.id,
	// 	g.chr,
	// 	e.start,
	// 	e.end,
	// 	g.strand,
	// 	g.gene_id,
	// 	g.gene_symbol,
	// 	gt.name as gene_type,
	// 	t.transcript_id,
	// 	t.is_canonical,
	// 	t.is_longest,
	// 	tt.name as transcript_type,
	// 	e.exon_id,
	// 	e.exon_number,

	return RowsToRecords(geneRows, level, canonicalMode)
}

func (genedb *V1GeneDB) SearchForGeneByName(search string,
	level string,
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

	if level == genome.GeneLevel {
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

func (genedb *V1GeneDB) WithinGenes(location *dna.Location, feature string,
	prom *dna.PromoterRegion) (*genome.GenomicFeatures, error) {

	// rows, err := genedb.withinGeneStmt.Query(
	// 	mid,
	// 	level,
	// 	location.Chr(),
	// 	location.Start(),
	// 	location.Start(),
	// 	location.End(),
	// 	location.End())

	rows, err := genedb.db.Query(InGeneSql,
		sql.Named("chr", location.Chr()),
		sql.Named("mid", location.Mid()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("prom5p", prom.Upstream()),
		sql.Named("prom3p", prom.Downstream()))

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return RowsToFeatures(location, feature, rows)
}

func (genedb *V1GeneDB) WithinGenesAndPromoter(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion) ([]*genome.GenomicFeature, error) {

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

	rows, err := genedb.db.Query(InGeneAndPromoterSql,
		sql.Named("chr", location.Chr()),
		sql.Named("mid", location.Mid()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("prom5p", prom.Upstream()),
		sql.Named("prom3p", prom.Downstream()))

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return RowsToRecords(rows, levels, false)
}

func (genedb *V1GeneDB) InExon(location *dna.Location, transcriptId string, prom *dna.PromoterRegion) ([]*genome.GenomicFeature, error) {

	// rows, err := genedb.inExonStmt.Query(
	// 	mid,
	// 	geneId,
	// 	location.Chr(),
	// 	location.Start(),
	// 	location.Start(),
	// 	location.End(),
	// 	location.End())

	rows, err := genedb.db.Query(InExonSql,
		sql.Named("transcriptId", transcriptId),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("mid", location.Mid()),
		sql.Named("prom5p", prom.Upstream()),
		sql.Named("prom3p", prom.Downstream()))

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return RowsToRecords(rows, genome.ExonLevel, false)
}

func (genedb *V1GeneDB) ClosestGenes(location *dna.Location,
	prom *dna.PromoterRegion,
	closestN int8) ([]*genome.GenomicFeature, error) {

	// rows, err := genedb.closestGeneStmt.Query(mid,
	// 	level,
	// 	location.Chr(),
	// 	mid,

	// 	n)

	rows, err := genedb.db.Query(ClosestGeneSql,
		sql.Named("chr", location.Chr()),
		sql.Named("mid", location.Mid()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("prom5p", prom.Upstream()),
		sql.Named("prom3p", prom.Downstream()),
		sql.Named("n", closestN))

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return RowsToRecords(rows, genome.GeneLevel, false)

	// genes, err := genome.RowsToRecords(rows)

	// if err != nil {
	// 	return nil, err //fmt.Errorf("there was an error with the database query")
	// }

	// ids := make([]string, len(genes))
	// for i, gene := range genes {
	// 	ids[i] = gene.GeneId
	// }

	// placeholders := make([]string, len(ids))
	// args := make([]any, len(ids)+1)

	// args[0] = mid

	// for i, id := range ids {
	// 	placeholders[i] = fmt.Sprintf("?%d", i+2)
	// 	args[i+1] = id
	// }

	// sql := StandardGeneFieldsSql +
	// 	`, ?1 - g.tss as tss_dist
	// 	FROM gtf AS g
	// 	WHERE g.feature = 'transcript'`

	// sql += fmt.Sprintf(" AND g.gene_id IN (%s)", strings.Join(placeholders, ",")) + " ORDER BY g.gene_id, ABS(tss_dist)"

	// rows, err = genedb.db.Query(sql, args...)

	// if err != nil {
	// 	return nil, err //fmt.Errorf("there was an error with the database query")
	// }

	// transcripts, err := genome.RowsToRecords(rows)

	// if err != nil {
	// 	return nil, err //fmt.Errorf("there was an error with the database query")
	// }

	// // find the closest transcript per gene

	// closesMap := make(map[string]*genome.GenomicFeature)

	// for _, feature := range transcripts {
	// 	existing, ok := closesMap[feature.GeneId]

	// 	if !ok || basemath.AbsInt(feature.TssDist) < basemath.AbsInt(existing.TssDist) {
	// 		closesMap[feature.GeneId] = feature
	// 	}
	// }

	// closest := make([]*genome.GenomicFeature, 0, closestN)

	// for _, feature := range genes {

	// 	transcript, ok := closesMap[feature.GeneId]

	// 	if ok {
	// 		closest = append(closest, transcript)
	// 	}
	// }

	// return closest, nil
}

func RowsToRecords(rows *sql.Rows, levels string, canonicalMode bool) ([]*genome.GenomicFeature, error) {
	var gid int
	var chr string
	var geneStart int
	var geneEnd int
	var strand string
	var geneSymbol string
	var geneId string
	var geneType string

	var transcriptId string
	var transcriptType string
	var transcriptStart int
	var transcriptEnd int

	var isCanonical bool
	var isLongest bool

	var exonId string
	var exonStart int
	var exonEnd int
	var exonNumber int

	var tssDist int

	var inPromoter bool
	var inExon bool
	var isIntragenic bool

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation

	var currentGene *genome.GenomicFeature
	var currentTranscript *genome.GenomicFeature
	var currentExon *genome.GenomicFeature

	var ret = make([]*genome.GenomicFeature, 0, 10)

	for rows.Next() {
		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &exonId)
		err := rows.Scan(&gid,
			&chr,
			&geneStart,
			&geneEnd,
			&strand,
			&geneId,
			&geneSymbol,
			&geneType,
			&transcriptId,
			&transcriptStart,
			&transcriptEnd,
			&isCanonical,
			&isLongest,
			&transcriptType,
			&exonId,
			&exonStart,
			&exonEnd,
			&exonNumber,
			&tssDist,
			&inPromoter,
			&inExon,
			&isIntragenic,
		)

		if err != nil {
			log.Error().Msgf("error reading overlapping gene rows %s", err)
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		// only add a new gene if we don't already have it. We
		// assume the rows are ordered by gene id hence if the
		// id changes, we are processing a set of rows for a new gene
		if strings.Contains(levels, "gene") && (currentGene == nil || currentGene.GeneId != geneId) {
			location, err := dna.NewStrandedLocation(chr, geneStart, geneEnd, strand)

			if err != nil {
				return nil, err
			}

			currentGene = &genome.GenomicFeature{Id: gid,
				Location:   location,
				Feature:    genome.GeneLevel,
				GeneSymbol: geneSymbol,
				GeneId:     geneId,
				//Strand:   strand,
				Type:     geneType,
				Children: make([]*genome.GenomicFeature, 0, 10)}

			ret = append(ret, currentGene)
		}

		if currentGene != nil {
			currentGene.InPromoter = currentGene.InPromoter || inPromoter
			currentGene.IsCanonical = currentGene.IsCanonical || isCanonical
			currentGene.IsLongest = currentGene.IsLongest || isLongest
			currentGene.InExon = currentGene.InExon || inExon
			currentGene.IsIntragenic = currentGene.IsIntragenic || isIntragenic
			currentGene.Label = genome.MakePromLabel(currentGene.InPromoter, currentGene.InExon, currentGene.IsIntragenic)

			if (basemath.AbsInt(tssDist) < basemath.AbsInt(currentGene.TssDist)) || currentGene.TssDist == 0 {
				currentGene.TssDist = tssDist
			}
		}

		//if canonical mode only add if transcript is canonical
		// also only add if we have a current gene
		// also only add if we don't already have this transcript
		if strings.Contains(levels, "transcript") &&
			(currentTranscript == nil || currentTranscript.TranscriptId != transcriptId) &&
			(!canonicalMode || isCanonical) {

			location, err := dna.NewStrandedLocation(chr, transcriptStart, transcriptEnd, strand)

			if err != nil {
				return nil, err
			}

			currentTranscript = &genome.GenomicFeature{Id: gid,
				Location: location,
				//Strand:       strand,
				Feature:      genome.TranscriptLevel,
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				IsCanonical:  isCanonical,
				IsLongest:    isLongest,
				InPromoter:   inPromoter,
				IsIntragenic: isIntragenic,
				TssDist:      tssDist,
				//InExon:       false,
				Type:     transcriptType,
				Children: make([]*genome.GenomicFeature, 0, 10)}

			if currentGene != nil {
				currentGene.Children = append(currentGene.Children, currentTranscript)
			} else {
				ret = append(ret, currentTranscript)
			}
		}

		if currentTranscript != nil {
			currentTranscript.InExon = currentTranscript.InExon || inExon
			currentTranscript.Label = genome.MakePromLabel(currentTranscript.InPromoter, currentTranscript.InExon, isIntragenic)
		}

		// only add exon if we have a current transcript and it matches
		// the transcript id

		if strings.Contains(levels, "exon") {
			location, err := dna.NewStrandedLocation(chr, exonStart, exonEnd, strand)

			if err != nil {
				return nil, err
			}

			currentExon = &genome.GenomicFeature{Id: gid,
				Location:     location,
				Feature:      genome.ExonLevel,
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				ExonId:       exonId,
				ExonNumber:   exonNumber,
				//InPromoter:   inPromoter,
				InExon: inExon,
				Label:  genome.MakePromLabel(inPromoter, inExon, true),
				//Strand:       strand,
			}

			if currentTranscript != nil {
				currentTranscript.Children = append(currentTranscript.Children, currentExon)
			} else if currentGene != nil {
				// normally add to transcript, but if no transcript, add to gene
				currentGene.Children = append(currentGene.Children, currentExon)
			} else {
				ret = append(ret, currentExon)
			}
		}
	}

	return ret, nil
}

func RowsToFeatures(location *dna.Location, levels string, rows *sql.Rows) (*genome.GenomicFeatures, error) {

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	features, err := RowsToRecords(rows, levels, false)

	if err != nil {
		log.Error().Msgf("error converting rows to features %s", err)
		return nil, err
	}

	var level = genome.MaxLevel(levels)

	ret := genome.GenomicFeatures{Location: location, Feature: level, Features: features}

	return &ret, nil
}
