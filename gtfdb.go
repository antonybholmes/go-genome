package genome

import (
	"database/sql"
	"fmt"
	"path/filepath"
	"slices"
	"strings"

	"github.com/antonybholmes/go-basemath"
	"github.com/antonybholmes/go-dna"

	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-sys/log"
)

type (
	GtfDB struct {
		db         *sql.DB
		annotation *Annotation
		file       string
	}

	// GtfDBInfo struct {
	// 	PublicId string `json:"id"`
	// 	Genome   string `json:"genome"`
	// 	Assembly string `json:"assembly"`
	// 	Name     string `json:"name"`
	// 	File     string `json:"-"`
	// 	Id       int    `json:"-"`
	// }

	GenomicFeatures struct {
		Type     string            `json:"type"`
		Features []*GenomicFeature `json:"features"`
	}

	GenomicFeature struct {
		PublicId     string        `json:"id,omitempty"`
		Location     *dna.Location `json:"loc"`
		Type         string        `json:"type,omitempty"`
		Biotype      string        `json:"biotype,omitempty"`
		GeneId       string        `json:"geneId,omitempty"`
		GeneSymbol   string        `json:"geneSymbol,omitempty"`
		TranscriptId string        `json:"transcriptId,omitempty"`
		ExonId       string        `json:"exonId,omitempty"`
		Label        string        `json:"label,omitempty"`

		Children []*GenomicFeature `json:"children,omitempty"`
		//Utrs         []*GenomicFeature `json:"utrs,omitempty"`
		//Exons        []*GenomicFeature `json:"exons,omitempty"`
		//Transcripts  []*GenomicFeature `json:"transcripts,omitempty"`
		TssDist      int  `json:"tssDist,omitempty"`
		Id           int  `json:"-"`
		ExonNumber   int  `json:"exonNumber,omitempty"`
		InPromoter   bool `json:"inPromoter,omitempty"`
		InExon       bool `json:"inExon,omitempty"`
		IsIntragenic bool `json:"isIntragenic,omitempty"`
		IsLongest    bool `json:"isLongest,omitempty"`
		IsCanonical  bool `json:"isCanonical,omitempty"`
	}

	GenomicSearchResults struct {
		Location *dna.Location `json:"location"`
		//Id       string        `json:"id,omitempty"`
		//Name     string        `json:"name,omitempty"`
		Type     string            `json:"type"`
		Features []*GenomicFeature `json:"features"`
	}
)

const (
	GeneDBInfoSql = `SELECT id, public_id, genome, assembly, name, file FROM info LIMIT 1`

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

	AllLevels               string = "gene,transcript,exon,cds,utr"
	GeneLevel               string = "gene"
	TranscriptLevel         string = "transcript"
	ExonLevel               string = "exon"
	GeneAndTranscriptLevels string = "gene,transcript"
	GeneAndExonLevels       string = "gene,exon"
	TranscriptAndExonLevels string = "transcript,exon"

	// search within transcripts and their promoters

	BasicLocationSql = `SELECT DISTINCT
		g.id, 
		g.chr,
		g.start, 
		g.end, 
		g.strand, 
		g.gene_id, 
		g.gene_symbol,
		gt.name AS gene_biotype,
		t.transcript_id,
		t.start,
		t.end,
		t.is_canonical,
		t.is_longest,
		ft.name AS feature_type,
		f.start,
		f.end,
		e.exon_id,
		e.exon_number
		FROM genes as g
		JOIN transcripts AS t ON g.id = t.gene_id
		JOIN features AS f ON f.transcript_id = t.id
		JOIN feature_types AS ft ON f.feature_type_id = ft.id
		JOIN exons AS e ON f.exon_id = e.id
		JOIN biotypes AS gt ON g.biotype_id = gt.id`

	CoreLocationSql = `SELECT DISTINCT
		g.id, 
		g.chr,
		g.start, 
		g.end, 
		g.strand, 
		g.gene_id, 
		g.gene_symbol,
		gt.name AS gene_biotype,
		t.transcript_id,
		t.start,
		t.end,
		t.is_canonical,
		t.is_longest,
		ft.name AS feature_type,
		f.start,
		f.end,
		e.exon_id,
		e.exon_number,
		CASE
			WHEN g.strand = '-' THEN :mid - t.end
			ELSE :mid - t.start
    	END AS tss_dist,
		((:start <= t.start + :prom3p) AND (:end >= t.start - :prom5p) AND (g.strand = '+')) OR
			((:start <= t.end + :prom5p) AND (:end >= t.end - :prom3p) AND (g.strand = '-')) 
		AS in_promoter,
		:start <= f.end AND :end >= f.start AS in_exon,
		:start <= t.end AND :end >= t.start AS is_intragenic
		FROM genes as g
		JOIN transcripts AS t ON g.id = t.gene_id
		JOIN features AS f ON f.transcript_id = t.id
		JOIN feature_types AS ft ON f.feature_type_id = ft.id
		JOIN exons AS e ON f.exon_id = e.id
		JOIN biotypes AS gt ON g.biotype_id = gt.id`

	CdsSql = `SELECT DISTINCT
		e.id,
		c.id,
		c.start,
		c.end
		FROM cds AS c
		JOIN exons AS e ON e.id = c.exon_id
		WHERE
			<<EXONS>>
		ORDER BY e.id, c.start`

	UtrSql = `SELECT DISTINCT
		e.id,
		u.id,
		u.start,
		u.end
		FROM utrs AS u
		JOIN exons AS e ON e.id = u.exon_id
		WHERE
			<<EXONS>>
		ORDER BY u.id, u.start`

	InGeneAndPromoterSql = CoreLocationSql +
		` WHERE g.chr = :chr AND 
		(
			((:start <= t.end) AND (:end >= t.start - :prom5p) AND g.strand = '+') OR
			((:start <= t.end + :prom5p) AND (:end >= t.start) AND g.strand = '-')
		)
		ORDER BY g.gene_id, t.transcript_id, e.exon_number`

	InGeneSql = CoreLocationSql +
		` WHERE g.chr = :chr AND (:start <= t.end AND :end >= t.start)`

	ClosestGeneSql = CoreLocationSql +
		` WHERE 
			t.transcript_id IN (
				SELECT DISTINCT t.transcript_id
				FROM transcripts AS t
				WHERE 
					g.chr = :chr AND
					t.is_longest = 1
				ORDER BY ABS(
					CASE 
						WHEN g.strand ='-' THEN :mid - t.end
						ELSE :mid - t.start
					END
				)
				LIMIT :n			
			) AND
			t.is_longest = 1
		ORDER BY ABS(tss_dist)`

	// when annotating genes, see if position falls within an exon
	InExonSql = CoreLocationSql +
		` WHERE e.transcript_id = :transcriptId AND (e.start <= :end AND e.end >= :start)
		ORDER BY e.start, e.end DESC`

	// order by gene, then transcript, then exon number, then feature type
	// so that exons come before cds and cds come before utrs, which is important for building the gene structure in memory
	BasicOverlapSql = BasicLocationSql +
		` WHERE 
			g.chr = :chr AND (t.start <= :end AND t.end >= :start)
			AND (:biotype = '' OR LOWER(gt.name) = :biotype)
		ORDER BY 
			g.gene_id,
			t.transcript_id,
			f.start,
			f.end,
			f.feature_type_id`

	OverlapSql = CoreLocationSql +
		` WHERE 
			g.chr = :chr AND (t.start <= :end AND t.end >= :start)
			AND (:biotype = '' OR LOWER(gt.name) = :biotype)
		ORDER BY 
			g.gene_id,
			t.transcript_id,
			f.start,
			f.end,
			f.feature_type_id`

	OverlapOrderBySql = ` ORDER BY 
		g.gene_id,
		t.transcript_id,
		e.exon_number`

	IdToNameSql = `SELECT name FROM ids WHERE id = ?1`
)

var (
	ExonLevelMap = map[string]string{
		"exon": "exons",
		"cds":  "cds",
		"utr":  "utrs",
	}
)

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func NewGtfDB(dir string, annotation *Annotation) *GtfDB {
	file := filepath.Join(dir, annotation.Url)

	log.Debug().Msgf("opening gene database %s", file)

	db := sys.Must(sql.Open(sys.Sqlite3DB, file+sys.SqliteReadOnlySuffix)) //sys.SqliteReadOnlySuffix

	return &GtfDB{annotation: annotation, file: file, db: db} //withinGeneStmt: sys.Must(db.Prepare(WITHIN_GENE_SQL)),
	//withinGeneAndPromStmt: sys.Must(db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL))

	//inExonStmt:      sys.Must(db.Prepare(IN_EXON_SQL)),
	//closestGeneStmt: sys.Must(db.Prepare(CLOSEST_GENE_SQL))

}

func (gtfdb *GtfDB) Close() error {
	return gtfdb.db.Close()
}

// func (gtfdb *GtfDB) LoadGeneDBInfo() (*GtfDBInfo, error) {

// 	var info GtfDBInfo

// 	err := gtfdb.db.QueryRow(GeneDBInfoSql).Scan(&info.Id,
// 		&info.PublicId,
// 		&info.Genome,
// 		&info.Assembly,
// 		&info.Name,
// 		&info.File)

// 	if err != nil {
// 		return nil, err //fmt.Errorf("there was an error with the database query")
// 	}

// 	return &info, nil
// }

func (gtfdb *GtfDB) OverlappingGenes(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion,
	canonicalMode bool,
	annotationMode bool,
	biotypeFilter string) ([]*GenomicFeature, error) {

	//log.Debug().Msgf("canonical mode %v gene type filter %s", canonicalMode, biotypeFilter)

	var geneRows *sql.Rows
	var err error

	levelTable := "exons"

	if strings.Contains(levels, "cds") {
		levelTable = "cds"
	} else if strings.Contains(levels, "utr") {
		levelTable = "utrs"
	} else {
		levelTable = "exons"
	}

	stmt := BasicOverlapSql

	if annotationMode {
		stmt = OverlapSql
	}

	stmt = strings.Replace(stmt, "<<LEVEL>>", levelTable, 1)

	//log.Debug().Msgf("querying overlapping genes with sql %s", stmt)

	geneRows, err = gtfdb.db.Query(stmt,
		sql.Named("chr", location.Chr()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("mid", location.Mid()),
		sql.Named("prom5p", prom.Upstream()),
		sql.Named("prom3p", prom.Downstream()),
		sql.Named("biotype", biotypeFilter))

	if err != nil {
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

	return rowsToRecords(geneRows, levels, canonicalMode, annotationMode)
}

func (gtfdb *GtfDB) WithinGenes(location *dna.Location, feature string,
	prom *dna.PromoterRegion) (*GenomicSearchResults, error) {

	// rows, err := genedb.withinGeneStmt.Query(
	// 	mid,
	// 	level,
	// 	location.Chr(),
	// 	location.Start(),
	// 	location.Start(),
	// 	location.End(),
	// 	location.End())

	rows, err := gtfdb.db.Query(InGeneSql,
		sql.Named("chr", location.Chr()),
		sql.Named("mid", location.Mid()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("prom5p", prom.Upstream()),
		sql.Named("prom3p", prom.Downstream()))

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	return rowsToFeatures(location, feature, rows, true)
}

func (gtfdb *GtfDB) WithinGenesAndPromoter(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion) ([]*GenomicFeature, error) {

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

	rows, err := gtfdb.db.Query(InGeneAndPromoterSql,
		sql.Named("chr", location.Chr()),
		sql.Named("mid", location.Mid()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("prom5p", prom.Upstream()),
		sql.Named("prom3p", prom.Downstream()))

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	return rowsToRecords(rows, levels, false, true)
}

func (gtfdb *GtfDB) InExon(location *dna.Location, transcriptId string, prom *dna.PromoterRegion) ([]*GenomicFeature, error) {

	// rows, err := genedb.inExonStmt.Query(
	// 	mid,
	// 	geneId,
	// 	location.Chr(),
	// 	location.Start(),
	// 	location.Start(),
	// 	location.End(),
	// 	location.End())

	rows, err := gtfdb.db.Query(InExonSql,
		sql.Named("transcriptId", transcriptId),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("mid", location.Mid()),
		sql.Named("prom5p", prom.Upstream()),
		sql.Named("prom3p", prom.Downstream()))

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	return rowsToRecords(rows, ExonLevel, false, true)
}

func (gtfdb *GtfDB) ClosestGenes(location *dna.Location,
	prom *dna.PromoterRegion,
	closestN int8) ([]*GenomicFeature, error) {

	// rows, err := genedb.closestGeneStmt.Query(mid,
	// 	level,
	// 	location.Chr(),
	// 	mid,

	// 	n)

	rows, err := gtfdb.db.Query(ClosestGeneSql,
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

	defer rows.Close()

	return rowsToRecords(rows, GeneLevel, false, true)

	// genes, err := rowsToRecords(rows)

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

	// transcripts, err := rowsToRecords(rows)

	// if err != nil {
	// 	return nil, err //fmt.Errorf("there was an error with the database query")
	// }

	// // find the closest transcript per gene

	// closesMap := make(map[string]*GenomicFeature)

	// for _, feature := range transcripts {
	// 	existing, ok := closesMap[feature.GeneId]

	// 	if !ok || basemath.AbsInt(feature.TssDist) < basemath.AbsInt(existing.TssDist) {
	// 		closesMap[feature.GeneId] = feature
	// 	}
	// }

	// closest := make([]*GenomicFeature, 0, closestN)

	// for _, feature := range genes {

	// 	transcript, ok := closesMap[feature.GeneId]

	// 	if ok {
	// 		closest = append(closest, transcript)
	// 	}
	// }

	// return closest, nil
}

func rowsToRecords(rows *sql.Rows, levels string, canonicalMode bool, annotationMode bool) ([]*GenomicFeature, error) {
	var gid int
	var chr string
	var geneStart int
	var geneEnd int
	var strand string
	var geneSymbol string
	var geneId string
	var geneBiotype string

	var transcriptId string
	//var transcriptBiotype string
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

	var currentGene *GenomicFeature
	var currentTranscript *GenomicFeature
	//var currentExon *GenomicFeature
	var currentFeature *GenomicFeature

	var ret = make([]*GenomicFeature, 0, 10)

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

			currentGene = &GenomicFeature{Id: gid,
				Location:   location,
				Type:       GeneLevel,
				GeneSymbol: geneSymbol,
				GeneId:     geneId,
				//Strand:   strand,
				Biotype: geneBiotype,

				//Children: make([]*GenomicFeature, 0, 10)
			}

			ret = append(ret, currentGene)
		}

		// gene inherits properties from transcripts and these become true
		// if any transcript has them true
		if currentGene != nil {
			currentGene.IsCanonical = currentGene.IsCanonical || isCanonical
			currentGene.IsLongest = currentGene.IsLongest || isLongest

			// add extra properties if in annotation mode
			if annotationMode {
				currentGene.InPromoter = currentGene.InPromoter || inPromoter
				currentGene.InExon = currentGene.InExon || inExon
				currentGene.IsIntragenic = currentGene.IsIntragenic || isIntragenic
				currentGene.Label = MakePromLabel(currentGene.InPromoter, currentGene.InExon, currentGene.IsIntragenic)

				if (basemath.AbsInt(tssDist) < basemath.AbsInt(currentGene.TssDist)) || currentGene.TssDist == 0 {
					currentGene.TssDist = tssDist
				}
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

			// set the properties that will not change for a transcript
			currentTranscript = &GenomicFeature{Id: gid,
				Location: location,
				//Strand:       strand,
				Type:         TranscriptLevel,
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				//Biotype:      transcriptBiotype,
				IsCanonical: isCanonical,
				IsLongest:   isLongest,
			}

			if annotationMode {
				currentTranscript.InPromoter = inPromoter
				currentTranscript.IsIntragenic = isIntragenic
				currentTranscript.TssDist = tssDist
			}

			if currentGene != nil {
				if currentGene.Children == nil {
					currentGene.Children = make([]*GenomicFeature, 0, 10)
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
			currentTranscript.Label = MakePromLabel(currentTranscript.InPromoter,
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

			currentFeature = &GenomicFeature{Id: gid,
				Location:     location,
				Type:         featureType,
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				ExonId:       exonId,
				ExonNumber:   exonNumber,
				//InPromoter:   inPromoter,

				//Label:  MakePromLabel(inPromoter, inExon, isIntragenic),
				//Strand:       strand,
			}

			if annotationMode {
				currentFeature.InExon = inExon
				currentFeature.Label = MakePromLabel(inPromoter, inExon, isIntragenic)
			}

			if currentTranscript != nil {
				// lazy initialize children slice only if we have exons to add
				if currentTranscript.Children == nil {
					currentTranscript.Children = make([]*GenomicFeature, 0, 10)
				}

				currentTranscript.Children = append(currentTranscript.Children, currentFeature)
			} else if currentGene != nil {
				// normally add to transcript, but if no transcript, add to gene
				if currentGene.Children == nil {
					currentGene.Children = make([]*GenomicFeature, 0, 10)
				}

				currentGene.Children = append(currentGene.Children, currentFeature)
			} else {
				// no gene or transcript, just add to return list
				ret = append(ret, currentFeature)
			}
		}
	}

	log.Debug().Msgf("converted rows to %d features", len(ret))

	return ret, nil
}

func rowsToFeatures(location *dna.Location, levels string, rows *sql.Rows, annotationMode bool) (*GenomicSearchResults, error) {

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	features, err := rowsToRecords(rows, levels, false, annotationMode)

	if err != nil {
		log.Error().Msgf("error converting rows to features %s", err)
		return nil, err
	}

	var level = MaxLevel(levels)

	ret := GenomicSearchResults{Location: location, Type: level, Features: features}

	return &ret, nil
}

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

// func (cache *GeneDBCache) Close() {
// 	for _, db := range cache.cacheMap {
// 		db.Close()
// 	}
// }

// func RowsToFeatures(location *dna.Location, level string, rows *sql.Rows) (*GenomicSearchResults, error) {

// 	//log.Debug().Msgf("rowsToFeatures %s   %d %d", location.Chr, location.Start, location.End)

// 	// 10 seems a reasonable guess for the number of features we might see, just
// 	// to reduce slice reallocation
// 	features, err := RowsToRecords(rows)

// 	if err != nil {
// 		return nil, err
// 	}

// 	ret := GenomicSearchResults{Location: location, Type: level, Features: features}

// 	return &ret, nil
// }

// func (genedb *GeneDB) IdToName(id int) string {

// 	name := "n/a"

// 	// ignore error as it is not critical
// 	// if the id is not found, we return "n/a"
// 	db.QueryRow(ID_TO_NAME_SQL, id).Scan(&name)

// 	return name
// }

// func RowsToRecords(rows *sql.Rows) ([]*GenomicFeature, error) {
// 	defer rows.Close()

// 	var id int
// 	var level string
// 	var chr string
// 	var start int
// 	var end int
// 	var strand string
// 	var geneId string
// 	var geneName string
// 	var biotype string
// 	var transcriptId string
// 	var transcriptName string
// 	var isCanonical bool
// 	var exonNumber int

// 	var d int
// 	var err error

// 	// 10 seems a reasonable guess for the number of features we might see, just
// 	// to reduce slice reallocation
// 	//var features = make([]*GenomicFeature, 0, 10)

// 	features := make([]*GenomicFeature, 0, 10)

// 	for rows.Next() {
// 		geneId = ""
// 		geneName = ""
// 		transcriptId = ""
// 		transcriptName = ""
// 		biotype = ""

// 		d = 0 // tss distance

// 		err = rows.Scan(&id,
// 			&level,
// 			&chr,
// 			&start,
// 			&end,
// 			&strand,
// 			&geneId,
// 			&geneName,
// 			&biotype,
// 			&isCanonical,
// 			&transcriptId,
// 			&transcriptName,
// 			&exonNumber,
// 			&d)

// 		if err != nil {
// 			return nil, err //fmt.Errorf("there was an error with the database records")
// 		}

// 		location, err := dna.NewLocation(chr, start, end)

// 		if err != nil {
// 			return nil, err
// 		}

// 		feature := GenomicFeature{
// 			Type:     level,
// 			Location: location,
// 			//Strand:       strand,
// 			GeneId:       geneId,
// 			GeneSymbol:   geneName,
// 			Biotype:      biotype,
// 			TranscriptId: transcriptId,
// 			IsCanonical:  isCanonical,
// 			TssDist:      d}

// 		features = append(features, &feature)
// 	}

// 	// enforce sorted correctly by chr and then position
// 	//sort.Sort(SortFeatureByPos(*features))

// 	return features, nil
// }

// SortFeaturesByPos sorts features in place by chr, start, end
func SortFeaturesByPos(features []*GenomicFeature) {
	slices.SortFunc(features, func(a, b *GenomicFeature) int {
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

func MaxLevel(levels string) string {
	var level string

	switch strings.ToLower(string(levels)) {
	case "gene", "genes":
		level = GeneLevel
	case "transcript", "transcripts":
		level = TranscriptLevel
	default:
		level = ExonLevel
	}

	return level
}
