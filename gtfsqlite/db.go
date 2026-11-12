package gtfsqlite

import (
	"database/sql"
	"fmt"
	"math"
	"path/filepath"
	"slices"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-genome"
	basemath "github.com/antonybholmes/go-sys/math"

	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-sys/db"
	"github.com/antonybholmes/go-sys/log"
)

type (
	GtfSqliteDB struct {
		db         *sql.DB
		annotation *genome.Annotation
		file       string
	}
)

const (
	GeneDBInfoSql = `SELECT id, 
		public_id, 
		genome, 
		assembly, 
		name, 
		file 
		FROM info LIMIT 1`

	MaxGeneInfoResults int16 = 100

	// const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_name, transcript_id, start - ?
	// 	FROM gene
	//  	WHERE level = 2 AND gene_id = ? AND chr = ? AND ? >= stranded_start - ? AND ? <= stranded_start + ?
	//  	ORDER BY start ASC`

	//GeneFeature       string = "gene"
	//TranscriptFeature string = "transcript"
	//ExonFeature       string = "exon"

	// GeneLevel       Level = 1
	// TranscriptLevel Level = 2
	// ExonLevel       Level = 3

	// search within transcripts and their promoters

	BasicLocationSql = `SELECT DISTINCT
		g.id, 
		c.name AS chr,
		g.start, 
		g.end, 
		g.strand, 
		g.gene_id, 
		g.symbol,
		gt.name AS gene_biotype,
		t.transcript_id,
		t.start,
		t.end,
		t.is_canonical,
		t.is_longest,
		ft.name AS feature_type,
		f.start,
		f.end,
		f.exon_id,
		f.exon_number
		FROM genes as g
		JOIN chromosomes AS c ON g.chr_id = c.id
		JOIN transcripts AS t ON g.id = t.gene_id
		JOIN features AS f ON f.transcript_id = t.id
		JOIN feature_types AS ft ON f.feature_type_id = ft.id
		JOIN biotypes AS gt ON g.biotype_id = gt.id`

	CoreLocationSql = `SELECT DISTINCT
		g.id, 
		c.name AS chr,
		g.start, 
		g.end, 
		g.strand, 
		g.gene_id, 
		g.symbol,
		gt.name AS gene_biotype,
		t.transcript_id,
		t.start,
		t.end,
		t.is_canonical,
		t.is_longest,
		ft.name AS feature_type,
		f.start,
		f.end,
		f.exon_id,
		f.exon_number,
		CASE
			WHEN g.strand = '+' THEN :mid - t.start
			ELSE :mid - t.end
		END AS tss_dist,
		((g.strand = '+') AND (:start <= t.start + :prom3p) AND (:end >= t.start - :prom5p)) OR
			((g.strand = '-') AND (:start <= t.end + :prom5p) AND (:end >= t.end - :prom3p)) 
		AS in_promoter,
		-- test only within exons
		f.feature_type_id = 1 AND :start <= f.end AND :end >= f.start AS in_exon,
		:start <= t.end AND :end >= t.start AS is_intragenic
	FROM genes as g
	JOIN chromosomes AS c ON g.chr_id = c.id
	JOIN transcripts AS t ON g.id = t.gene_id
	JOIN features AS f ON f.transcript_id = t.id
	JOIN feature_types AS ft ON f.feature_type_id = ft.id
	JOIN biotypes AS gt ON g.biotype_id = gt.id`

	// Within genes or their promoters
	// IntragenicSql = `SELECT DISTINCT
	// 		g.id,
	// 		c.name AS chr,
	// 		g.start,
	// 		g.end,
	// 		g.strand,
	// 		g.gene_id,
	// 		g.symbol,
	// 		gt.name AS gene_biotype,
	// 		t.transcript_id,
	// 		t.start,
	// 		t.end,
	// 		t.is_canonical,
	// 		t.is_longest,
	// 		ft.name AS feature_type,
	// 		f.start,
	// 		f.end,
	// 		e.exon_id,
	// 		e.exon_number,
	// 		CASE
	// 			WHEN g.strand = '+' THEN :mid - t.start
	// 			ELSE :mid - t.end
	// 		END AS tss_dist,
	// 		((g.strand = '+') AND (:start <= t.start + :prom3p) AND (:end >= t.start - :prom5p)) OR
	// 			((g.strand = '-') AND (:start <= t.end + :prom5p) AND (:end >= t.end - :prom3p))
	// 		AS in_promoter,
	// 		:start <= f.end AND :end >= f.start AS in_exon,
	// 		:start <= t.end AND :end >= t.start AS is_intragenic
	// 	FROM genes as g
	// 	JOIN chromosomes AS c ON g.chr_id = c.id
	// 	JOIN transcripts AS t ON g.id = t.gene_id
	// 	JOIN features AS f ON f.transcript_id = t.id
	// 	JOIN feature_types AS ft ON f.feature_type_id = ft.id
	// 	JOIN exons AS e ON f.exon_id = e.id
	// 	JOIN biotypes AS gt ON g.biotype_id = gt.id
	// 	WHERE c.name = :chr AND
	// 		(
	// 			(g.strand = '+' AND (:start <= t.end) AND (:end >= t.start - :prom5p)) OR
	// 			(g.strand = '-' AND (:start <= t.end + :prom5p) AND (:end >= t.start))
	// 		) AND
	// 		(:use_official = 0 OR g.official_gene_id IS NOT NULL)
	// 		ORDER BY g.gene_id, t.transcript_id, e.exon_number`

	IntragenicSql = CoreLocationSql +
		` WHERE c.name = :chr AND 
			(
				(g.strand = '+' AND (:start <= t.end) AND (:end >= t.start - :prom5p)) OR
				(g.strand = '-' AND (:start <= t.end + :prom5p) AND (:end >= t.start))
			) AND
			(:use_official = 0 OR g.official_gene_id IS NOT NULL)
			ORDER BY g.gene_id, t.transcript_id, e.exon_number`

	ClosestGeneSql = `WITH ranked_transcripts AS 
		(
			SELECT DISTINCT
				g.id as gene_id,
				t.id AS transcript_id,
				ABS(
					CASE 
						WHEN g.strand ='+' THEN :mid - t.start
						ELSE :mid - t.end
					END
				) as tss_dist,
				-- Rank transcripts within the SAME gene by closest distance
				ROW_NUMBER() OVER (
					PARTITION BY t.gene_id 
					ORDER BY 
						ABS(
							CASE 
								WHEN g.strand ='+' THEN :mid - t.start
								ELSE :mid - t.end
							END
						) ASC
				) AS gene_transcript_rank
			FROM transcripts AS t
			JOIN genes AS g ON t.gene_id = g.id
			JOIN chromosomes AS c ON g.chr_id = c.id
			WHERE 
				c.name = :chr AND
				-- avoid annotating to genes with an ENSG symbol as these are likely to be less well characterized 
				(:use_official = 0 OR g.official_gene_id IS NOT NULL)
		),
		closest_transcripts AS (
			-- get the ids of the closest transcripts for the closest genes
			SELECT DISTINCT
				ROW_NUMBER() OVER (ORDER BY tss_dist ASC) AS rank,
				transcript_id
			FROM 
				ranked_transcripts
			WHERE
				gene_transcript_rank = 1
			ORDER BY 
				tss_dist ASC
			LIMIT 
				:n	
		)
		SELECT DISTINCT
			g.id, 
			c.name AS chr,
			g.start, 
			g.end, 
			g.strand, 
			g.gene_id, 
			g.symbol,
			gt.name AS gene_biotype,
			t.transcript_id,
			t.start,
			t.end,
			t.is_canonical,
			t.is_longest,
			ft.name AS feature_type,
			f.start,
			f.end,
			f.exon_id,
			f.exon_number,
			CASE
				WHEN g.strand = '+' THEN :mid - t.start
				ELSE :mid - t.end
			END AS tss_dist,
			((g.strand = '+') AND (:start <= t.start + :prom3p) AND (:end >= t.start - :prom5p)) OR
				((g.strand = '-') AND (:start <= t.end + :prom5p) AND (:end >= t.end - :prom3p)) 
			AS in_promoter,
			f.feature_type_id = 1 AND :start <= f.end AND :end >= f.start AS in_exon,
			:start <= t.end AND :end >= t.start AS is_intragenic
		FROM genes as g
		JOIN chromosomes AS c ON g.chr_id = c.id
		JOIN transcripts AS t ON g.id = t.gene_id
		JOIN closest_transcripts ct ON ct.transcript_id = t.id
		JOIN features AS f ON f.transcript_id = t.id
		JOIN feature_types AS ft ON f.feature_type_id = ft.id
		JOIN biotypes AS gt ON g.biotype_id = gt.id
		ORDER BY 
			ct.rank ASC`

	// when annotating genes, see if position falls within an exon
	// InExonSql = CoreLocationSql +
	// 	` WHERE e.transcript_id = :transcriptId AND (e.start <= :end AND e.end >= :start)
	// 	ORDER BY e.start, e.end DESC`

	// do we overlap with an exon?
	InExonSql = `SELECT COUNT(DISTINCT f.id)
		FROM transcripts t
		JOIN features AS f ON f.transcript_id = t.id
		-- only consider exons
		WHERE f.feature_type_id = 1 AND f.transcript_id = :transcriptId AND (f.start <= :end AND f.end >= :start)`

	// order by gene, then transcript, then exon number, then feature type
	// so that exons come before cds and cds come before utrs, which is important for building the gene structure in memory
	BasicOverlapSql = BasicLocationSql +
		` WHERE 
			c.name = :chr AND (t.start <= :end AND t.end >= :start)
			AND (:biotype = '' OR LOWER(gt.name) = :biotype)
		ORDER BY 
			g.gene_id,
			t.transcript_id,
			f.start,
			f.end,
			f.feature_type_id`

	OverlapSql = CoreLocationSql +
		` WHERE 
			c.name = :chr AND (t.start <= :end AND t.end >= :start)
			AND (:biotype = '' OR LOWER(gt.name) = :biotype)
		ORDER BY 
			g.gene_id,
			t.transcript_id,
			f.start,
			f.end,
			f.feature_type_id`

	// InGeneSql = CoreLocationSql +
	// 	` WHERE
	// 		c.name = :chr AND (:start <= t.end AND :end >= t.start)
	// 	ORDER BY
	// 		g.gene_id,
	// 		t.transcript_id,
	// 		f.start,
	// 		f.end,
	// 		f.feature_type_id`

	OverlapOrderBySql = ` ORDER BY 
		g.gene_id,
		t.transcript_id,
		e.exon_number`

	IdToNameSql = `SELECT name FROM ids WHERE id = ?1`
)

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func NewGtfSqliteDB(dir string, annotation *genome.Annotation) *GtfSqliteDB {
	file := filepath.Join(dir, annotation.Url)

	log.Debug().Msgf("opening gene database %s", file)

	db := sys.Must(sql.Open(db.Sqlite3DB, file+db.SqliteReadOnlySuffix))

	return &GtfSqliteDB{annotation: annotation, file: file, db: db}
}

func (gdb *GtfSqliteDB) Close() error {
	return gdb.db.Close()
}

func (gdb *GtfSqliteDB) OverlappingGenes(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion,
	canonicalMode bool,
	annotationMode bool,
	biotypeFilter string) ([]*genome.GenomicFeature, error) {

	//log.Debug().Msgf("canonical mode %v gene type filter %s", canonicalMode, biotypeFilter)

	var geneRows *sql.Rows
	var err error

	//levelTable := "exons"

	// if strings.Contains(levels, "cds") {
	// 	levelTable = "cds"
	// } else if strings.Contains(levels, "utr") {
	// 	levelTable = "utrs"
	// } else {
	// 	levelTable = "exons"
	// }

	stmt := BasicOverlapSql

	if annotationMode {
		stmt = OverlapSql
	}

	//stmt = strings.Replace(stmt, "<<LEVEL>>", levelTable, 1)

	//log.Debug().Msgf("querying overlapping genes with sql %s", stmt)

	geneRows, err = gdb.db.Query(stmt,
		sql.Named("chr", location.Chr()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("mid", location.Mid()),
		sql.Named("prom5p", prom.FivePrime()),
		sql.Named("prom3p", prom.ThreePrime()),
		sql.Named("biotype", biotypeFilter))

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer geneRows.Close()

	return rowsToRecords(geneRows, levels, canonicalMode, annotationMode)
}

// func (gdb *GtfSqliteDB) WithinGenes(location *dna.Location,
// 	feature string,
// 	prom *dna.PromoterRegion) (*genome.GenomicSearchResults, error) {

// 	rows, err := gdb.db.Query(InGeneSql,
// 		sql.Named("chr", location.Chr()),
// 		sql.Named("mid", location.Mid()),
// 		sql.Named("start", location.Start()),
// 		sql.Named("end", location.End()),
// 		sql.Named("prom5p", prom.FivePrime()),
// 		sql.Named("prom3p", prom.ThreePrime()))

// 	if err != nil {
// 		return nil, err //fmt.Errorf("there was an error with the database query")
// 	}

// 	defer rows.Close()

// 	return rowsToSearchResults(location, feature, rows, true)
// }

func (gdb *GtfSqliteDB) IntragenicFeatures(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion,
	useOfficialGenes bool) ([]*genome.GenomicFeature, error) {

	rows, err := gdb.db.Query(IntragenicSql,
		sql.Named("chr", location.Chr()),
		sql.Named("mid", location.Mid()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("prom5p", prom.FivePrime()),
		sql.Named("prom3p", prom.ThreePrime()),
		sql.Named("use_official", useOfficialGenes)) // only consider genes with an official id for intragenic annotation

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	return rowsToRecords(rows, levels, false, true)
}

func (gdb *GtfSqliteDB) InExon(location *dna.Location,
	transcriptId string,
	prom *dna.PromoterRegion) (bool, error) {

	var count int

	err := gdb.db.QueryRow(InExonSql,
		sql.Named("transcriptId", transcriptId),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End())).Scan(&count)

	// rows, err := gdb.db.Query(InExonSql,
	// 	sql.Named("transcriptId", transcriptId),
	// 	sql.Named("start", location.Start()),
	// 	sql.Named("end", location.End()),
	// 	sql.Named("mid", location.Mid()),
	// 	sql.Named("prom5p", prom.FivePrime()),
	// 	sql.Named("prom3p", prom.ThreePrime()))

	if err != nil {
		return false, err //fmt.Errorf("there was an error with the database query")
	}

	return count > 0, nil
}

func (gdb *GtfSqliteDB) ClosestGenes(location *dna.Location,
	prom *dna.PromoterRegion,
	closestN int8,
	useOfficialGenes bool) ([]*genome.GenomicFeature, error) {

	///log.Debug().Msgf("querying closest genes with sql %s", ClosestGeneSql)

	rows, err := gdb.db.Query(ClosestGeneSql,
		sql.Named("chr", location.Chr()),
		sql.Named("mid", location.Mid()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("prom5p", prom.FivePrime()),
		sql.Named("prom3p", prom.ThreePrime()),
		sql.Named("n", closestN),
		sql.Named("use_official", useOfficialGenes)) // only consider genes with an official id for closest gene annotation

	if err != nil {
		return nil, err
	}

	defer rows.Close()

	return rowsToRecords(rows, genome.GeneLevel, false, true)

}

func rowsToRecords(rows *sql.Rows, levels string, canonicalMode bool, annotationMode bool) ([]*genome.GenomicFeature, error) {
	var gid int
	var chr string
	var geneStart int
	var geneEnd int
	var strand string
	var geneSymbol string
	var geneId string
	var geneBiotype string

	var transcriptId string

	var transcriptStart int
	var transcriptEnd int

	var isCanonical bool
	var isLongest bool

	var exonId string
	var featureType string
	var featureStart int
	var featureEnd int
	var exonNumber int

	var tssDist int

	var inPromoter bool
	var inExon bool
	var isIntragenic bool

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation

	var currentGene *genome.GenomicFeature
	var currentTranscript *genome.GenomicFeature

	var currentFeature *genome.GenomicFeature

	var ret = make([]*genome.GenomicFeature, 0, 10)

	//exonMap := make(map[int]*GenomicFeature)

	var err error

	for rows.Next() {
		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &exonId)
		if annotationMode {
			err = rows.Scan(&gid,
				&chr,
				&geneStart,
				&geneEnd,
				&strand,
				&geneId,
				&geneSymbol,
				&geneBiotype,
				&transcriptId,
				&transcriptStart,
				&transcriptEnd,
				&isCanonical,
				&isLongest,
				&featureType, // are an exon, cds, or utr
				&featureStart,
				&featureEnd,
				&exonId, // tie feature to an exon
				&exonNumber,
				&tssDist,
				&inPromoter,
				&inExon,
				&isIntragenic,
			)
		} else {
			// a shorter query for non-annotation mode when you just
			// want coordinates without extra checks to see if in exon or not etc.
			err = rows.Scan(&gid,
				&chr,
				&geneStart,
				&geneEnd,
				&strand,
				&geneId,
				&geneSymbol,
				&geneBiotype,
				&transcriptId,
				&transcriptStart,
				&transcriptEnd,
				&isCanonical,
				&isLongest,
				&featureType, // are an exon, cds, or utr
				&featureStart,
				&featureEnd,
				&exonId, // tie feature to an exon
				&exonNumber,
			)
		}

		if err != nil {
			log.Error().Msgf("error reading overlapping gene rows %s", err)
			return nil, err
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
				Location: location,
				Type:     genome.GeneLevel,
				Symbol:   geneSymbol,
				GeneId:   geneId,
				//Strand:   strand,
				Biotype: geneBiotype,
				TssDist: math.MaxInt,
			}

			ret = append(ret, currentGene)
		}

		// gene inherits properties from transcripts and these become true
		// if any transcript has them true
		if currentGene != nil {
			currentGene.IsCanonical = currentGene.IsCanonical || isCanonical
			//currentGene.IsLongest = currentGene.IsLongest || isLongest

			// add extra properties if in annotation mode
			if annotationMode {
				currentGene.InPromoter = currentGene.InPromoter || inPromoter
				currentGene.InExon = currentGene.InExon || inExon
				currentGene.IsIntragenic = currentGene.IsIntragenic || isIntragenic
				currentGene.Label = genome.MakePromLabel(currentGene.InPromoter, currentGene.InExon, currentGene.IsIntragenic)

			}
		}

		//if canonical mode only add if transcript is canonical
		// also only add if we have a current gene
		// also only add if we don't already have this transcript
		if strings.Contains(levels, "transcript") &&
			(currentTranscript == nil || currentTranscript.Transcript != transcriptId) &&
			(!canonicalMode || isCanonical) {

			location, err := dna.NewStrandedLocation(chr, transcriptStart, transcriptEnd, strand)

			if err != nil {
				return nil, err
			}

			// set the properties that will not change for a transcript
			currentTranscript = &genome.GenomicFeature{Id: gid,
				Location: location,
				//Strand:       strand,
				Type:       genome.TranscriptLevel,
				Symbol:     geneSymbol,
				GeneId:     geneId,
				Transcript: transcriptId,
				//Biotype:      transcriptBiotype,
				IsCanonical: isCanonical,
				//IsLongest:   isLongest,
				TssDist: tssDist,
			}

			if annotationMode {
				currentTranscript.InPromoter = inPromoter
				currentTranscript.IsIntragenic = isIntragenic
			}

			if currentGene != nil {
				if (basemath.AbsInt(tssDist) < basemath.AbsInt(currentGene.TssDist)) || currentGene.TssDist == 0 {
					currentGene.TssDist = tssDist
				}

				if currentGene.Children == nil {
					currentGene.Children = make([]*genome.GenomicFeature, 0, 10)
				}

				currentGene.Children = append(currentGene.Children, currentTranscript)
			} else {
				ret = append(ret, currentTranscript)
			}
		}

		if currentTranscript != nil {
			// these properties may be updated as we see more rows and we find
			// we are in an exon
			currentTranscript.InExon = currentTranscript.InExon || inExon
			currentTranscript.Label = genome.MakePromLabel(currentTranscript.InPromoter,
				currentTranscript.InExon,
				currentTranscript.IsIntragenic)
		}

		// only add exon if we have a current transcript and it matches
		// the transcript id

		// exons are a bit more complex because they can be tagged as exon, cds, or utr and we want
		// to capture that information in the feature level. We also want to add the exon
		// to the transcript if we have a current transcript, but if not, we add it to the gene.
		// If we don't have a gene, we just add it to the return list.
		// This is because some databases may not have transcripts and we still want to
		// capture exon information if it is available.
		if strings.Contains(levels, featureType) {
			location, err := dna.NewStrandedLocation(chr, featureStart, featureEnd, strand)

			if err != nil {
				return nil, err
			}

			currentFeature = &genome.GenomicFeature{Id: gid,
				Location:   location,
				Type:       featureType,
				Symbol:     geneSymbol,
				GeneId:     geneId,
				Transcript: transcriptId,
				Exon:       exonId,
				ExonNumber: exonNumber,
			}

			if annotationMode {
				currentFeature.InExon = inExon
				currentFeature.Label = genome.MakePromLabel(inPromoter, inExon, isIntragenic)
			}

			if currentTranscript != nil {
				// lazy initialize children slice only if we have exons to add
				if currentTranscript.Children == nil {
					currentTranscript.Children = make([]*genome.GenomicFeature, 0, 10)
				}

				currentTranscript.Children = append(currentTranscript.Children, currentFeature)
			} else if currentGene != nil {
				// normally add to transcript, but if no transcript, add to gene
				if currentGene.Children == nil {
					currentGene.Children = make([]*genome.GenomicFeature, 0, 10)
				}

				currentGene.Children = append(currentGene.Children, currentFeature)
			} else {
				// no gene or transcript, just add to return list
				ret = append(ret, currentFeature)
			}
		}
	}

	//log.Debug().Msgf("converted rows to %d features", len(ret))

	return ret, nil
}

// func rowsToSearchResults(location *dna.Location,
// 	levels string,
// 	rows *sql.Rows,
// 	annotationMode bool) (*genome.GenomicSearchResults, error) {

// 	// 10 seems a reasonable guess for the number of features we might see, just
// 	// to reduce slice reallocation
// 	features, err := rowsToRecords(rows, levels, false, annotationMode)

// 	if err != nil {
// 		log.Error().Msgf("error converting rows to features %s", err)
// 		return nil, err
// 	}

// 	return genome.FeaturesToSearchResults(location, levels, features, annotationMode)
// }

func MakeInExonsSql(query, table string, exons []int, namedArgs *[]any) string {

	inPlaceholders := make([]string, len(exons))

	for i, exon := range exons {
		ph := fmt.Sprintf("e%d", i+1)
		inPlaceholders[i] = ":" + ph
		*namedArgs = append(*namedArgs, sql.Named(ph, exon))
	}

	// table e.g. c.exon_id or u.exon_id depending on if we are building the sql for cds or utrs
	clause := table + ".exon_id IN (" + strings.Join(inPlaceholders, ",") + ")"

	return strings.Replace(query, "<<EXONS>>", clause, 1)

}

// SortFeaturesByPos sorts features in place by chr, start, end
func SortFeaturesByPos(features []*genome.GenomicFeature) {
	slices.SortFunc(features, func(a, b *genome.GenomicFeature) int {
		ci := dna.ChromToInt(a.Location.Chr())
		cj := dna.ChromToInt(b.Location.Chr())

		// on different chrs so sort by chr
		if ci != cj {
			return int(ci) - int(cj)
		}

		// same chr so sort by position
		if a.Location.Start() != b.Location.Start() {
			return int(a.Location.Start()) - int(b.Location.Start())
		}

		// same start so sort by end
		return int(a.Location.End()) - int(b.Location.End())
	})
}
