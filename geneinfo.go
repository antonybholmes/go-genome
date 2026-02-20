package genome

import (
	"database/sql"
	"fmt"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys/log"
)

//
// Represents the search by name features of the db
// i.e. looking for a gene by name to get its position
// etc.
//

const (
	GeneInfoSql = `SELECT DISTINCT
			g.id, 
			g.chr,
			g.start, 
			g.end, 
			g.strand, 
			g.gene_id, 
			g.gene_symbol,
			gt.name AS gene_type
	FROM genes as g
	JOIN gene_types AS gt ON g.gene_type_id = gt.id
	WHERE (g.gene_symbol LIKE :q OR g.gene_id LIKE :q)
	ORDER BY g.gene_symbol
	LIMIT :n`

	TranscriptInfoSql = `SELECT *
		FROM(
			SELECT DISTINCT
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
				ROW_NUMBER() OVER (PARTITION BY g.gene_id ORDER BY g.gene_id) AS rank
			FROM genes as g
			JOIN transcripts AS t ON g.id = t.gene_id
			JOIN gene_types AS gt ON g.gene_type_id = gt.id
			JOIN transcript_types AS tt ON t.transcript_type_id = tt.id 
			WHERE (g.gene_symbol LIKE :q OR g.gene_id LIKE :q OR t.transcript_id LIKE :q) {{whereClause}}
			ORDER BY g.gene_symbol
		) g
		WHERE g.rank < :n`

	ExonInfoSql = `SELECT *
		FROM(
			SELECT DISTINCT
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
				ROW_NUMBER() OVER (PARTITION BY g.gene_id ORDER BY g.gene_id) AS rank
			FROM genes as g
			JOIN transcripts AS t ON g.id = t.gene_id
			JOIN exons AS e ON e.transcript_id = t.id
			JOIN gene_types AS gt ON g.gene_type_id = gt.id
			JOIN transcript_types AS tt ON t.transcript_type_id = tt.id 
			WHERE (g.gene_symbol LIKE :q OR g.gene_id LIKE :q OR t.transcript_id LIKE :q OR e.exon_id LIKE :q) {{whereClause}}
			ORDER BY g.gene_symbol
		) g
		WHERE g.rank < :n`
)

func (gtfdb *GtfDB) SearchByName(search string,
	level string,
	fuzzyMode bool,
	canonicalMode bool,
	geneType string,
	n int16) ([]*GenomicFeature, error) {
	n = max(1, min(n, MaxGeneInfoResults))

	//log.Debug().Msgf("SearchForGeneByName for gene %s level %s fuzzy %v canonical %v geneType %s n %d",
	//	search, level, fuzzy, canonical, geneType, n)

	// case insensitive search
	search = strings.ToLower(search)

	if len(search) < 2 || strings.Contains(search, "chr:") {
		return nil, fmt.Errorf("%s is an invalid search term", search)
	}

	switch level {
	case "transcript":
		return gtfdb.searchTranscripts(search,
			canonicalMode,
			false,
			n)

	case "exon":
		return gtfdb.searchTranscripts(search,
			canonicalMode,
			true,
			n)

	default:
		log.Debug().Msgf("searching for genes with term %s", search)
		return gtfdb.searchGenes(search, n)
	}

}

// Searching for exons or transcripts uses essentially
// the same pipeline so combine into one method.
func (gtfdb *GtfDB) searchTranscripts(search string,
	canonicalMode bool,
	exonMode bool,
	n int16) ([]*GenomicFeature, error) {
	n = max(1, min(n, MaxGeneInfoResults))

	var sqlStmt string

	if exonMode {
		sqlStmt = ExonInfoSql
	} else {
		sqlStmt = TranscriptInfoSql
	}

	if canonicalMode {
		sqlStmt = strings.Replace(sqlStmt, "{{whereClause}}", "AND t.is_canonical = 1", 1)
	} else {
		sqlStmt = strings.Replace(sqlStmt, "{{whereClause}}", "", 1)
	}

	//log.Debug().Msgf("SQL: %s %s %d", sqlStmt, search, n)

	rows, err := gtfdb.db.Query(sqlStmt,
		sql.Named("q", search),
		sql.Named("n", n))

	if err != nil {
		log.Debug().Msgf("error querying gene by name %s", err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	if exonMode {
		return exonsToGeneInfoRecords(rows, canonicalMode)
	} else {
		return transcriptsToGeneInfoRecords(rows, canonicalMode)
	}

}

func (gtfdb *GtfDB) searchGenes(search string,
	n int16) ([]*GenomicFeature, error) {
	n = max(1, min(n, MaxGeneInfoResults))

	//log.Debug().Msgf("SearchForGeneByName for gene %s level %s fuzzy %v canonical %v geneType %s n %d",
	//	search, level, fuzzy, canonical, geneType, n)

	//log.Debug().Msgf("SQL: %s %s %d", sqlStmt, search, n)

	rows, err := gtfdb.db.Query(GeneInfoSql,
		sql.Named("q", search),
		sql.Named("n", n))

	if err != nil {
		log.Debug().Msgf("error querying gene by name %s", err)
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	ret, err := genesToGeneInfoRecords(rows)

	if err != nil {
		return nil, err
	}

	return ret, nil
}

// complete records for genes -> transcripts -> exons
func genesToGeneInfoRecords(rows *sql.Rows) ([]*GenomicFeature, error) {
	var gid int
	var chr string
	var geneStart int
	var geneEnd int
	var strand string
	var geneSymbol string
	var geneId string
	var geneType string

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation

	var ret = make([]*GenomicFeature, 0, 10)

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
		)

		if err != nil {
			log.Error().Msgf("error reading gene info rows %s", err)
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		// only add a new gene if we don't already have it. We
		// assume the rows are ordered by gene id hence if the
		// id changes, we are processing a set of rows for a new gene

		location, err := dna.NewStrandedLocation(chr, geneStart, geneEnd, strand)

		if err != nil {
			return nil, err
		}

		currentGene := &GenomicFeature{Id: gid,
			Location:   location,
			Feature:    GeneLevel,
			GeneSymbol: geneSymbol,
			GeneId:     geneId,
			//Strand:   strand,
			Type: geneType,
			//Children: make([]*GenomicFeature, 0, 10)
		}

		ret = append(ret, currentGene)

	}

	log.Debug().Msgf("converted rows to %d features", len(ret))

	return ret, nil
}

// complete records for genes -> transcripts -> exons
func transcriptsToGeneInfoRecords(rows *sql.Rows, canonicalMode bool) ([]*GenomicFeature, error) {
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

	var rank int

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation

	var currentGene *GenomicFeature
	var currentTranscript *GenomicFeature

	var ret = make([]*GenomicFeature, 0, 10)

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
			&rank,
		)

		if err != nil {
			log.Error().Msgf("error reading gene info rows %s", err)
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		// only add a new gene if we don't already have it. We
		// assume the rows are ordered by gene id hence if the
		// id changes, we are processing a set of rows for a new gene
		if currentGene == nil || currentGene.GeneId != geneId {
			location, err := dna.NewStrandedLocation(chr, geneStart, geneEnd, strand)

			if err != nil {
				return nil, err
			}

			currentGene = &GenomicFeature{Id: gid,
				Location:   location,
				Feature:    GeneLevel,
				GeneSymbol: geneSymbol,
				GeneId:     geneId,
				//Strand:   strand,
				Type: geneType,
				//Children: make([]*GenomicFeature, 0, 10)
			}

			ret = append(ret, currentGene)
		}

		// gene inherits properties from transcripts and these become true
		// if any transcript has them true
		if currentGene != nil {
			currentGene.IsCanonical = currentGene.IsCanonical || isCanonical
			currentGene.IsLongest = currentGene.IsLongest || isLongest

		}

		//if canonical mode only add if transcript is canonical
		// also only add if we have a current gene
		// also only add if we don't already have this transcript
		if (currentTranscript == nil || currentTranscript.TranscriptId != transcriptId) &&
			(!canonicalMode || isCanonical) {

			location, err := dna.NewStrandedLocation(chr, transcriptStart, transcriptEnd, strand)

			if err != nil {
				return nil, err
			}

			// set the properties that will not change for a transcript
			currentTranscript = &GenomicFeature{Id: gid,
				Location: location,
				//Strand:       strand,
				Feature:      TranscriptLevel,
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				IsCanonical:  isCanonical,
				IsLongest:    isLongest,
				Type:         transcriptType,
				//Children: make([]*GenomicFeature, 0, 10)
			}

			if currentGene != nil {
				if currentGene.Features == nil {
					currentGene.Features = make([]*GenomicFeature, 0, 10)
				}

				currentGene.Features = append(currentGene.Features, currentTranscript)
			} else {
				ret = append(ret, currentTranscript)
			}
		}
	}

	log.Debug().Msgf("converted rows to %d features", len(ret))

	return ret, nil
}

// complete records for genes -> transcripts -> exons
func exonsToGeneInfoRecords(rows *sql.Rows, canonicalMode bool) ([]*GenomicFeature, error) {
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

	var rank int

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation

	var currentGene *GenomicFeature
	var currentTranscript *GenomicFeature
	var currentExon *GenomicFeature

	var ret = make([]*GenomicFeature, 0, 10)

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
			&rank,
		)

		if err != nil {
			log.Error().Msgf("error reading gene info rows %s", err)
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		// only add a new gene if we don't already have it. We
		// assume the rows are ordered by gene id hence if the
		// id changes, we are processing a set of rows for a new gene
		if currentGene == nil || currentGene.GeneId != geneId {
			location, err := dna.NewStrandedLocation(chr, geneStart, geneEnd, strand)

			if err != nil {
				return nil, err
			}

			currentGene = &GenomicFeature{Id: gid,
				Location:   location,
				Feature:    GeneLevel,
				GeneSymbol: geneSymbol,
				GeneId:     geneId,
				//Strand:   strand,
				Type: geneType,
				//Children: make([]*GenomicFeature, 0, 10)
			}

			ret = append(ret, currentGene)
		}

		// gene inherits properties from transcripts and these become true
		// if any transcript has them true
		if currentGene != nil {
			currentGene.IsCanonical = currentGene.IsCanonical || isCanonical
			currentGene.IsLongest = currentGene.IsLongest || isLongest

		}

		//if canonical mode only add if transcript is canonical
		// also only add if we have a current gene
		// also only add if we don't already have this transcript
		if (currentTranscript == nil || currentTranscript.TranscriptId != transcriptId) &&
			(!canonicalMode || isCanonical) {

			location, err := dna.NewStrandedLocation(chr, transcriptStart, transcriptEnd, strand)

			if err != nil {
				return nil, err
			}

			// set the properties that will not change for a transcript
			currentTranscript = &GenomicFeature{Id: gid,
				Location: location,
				//Strand:       strand,
				Feature:      TranscriptLevel,
				GeneSymbol:   geneSymbol,
				GeneId:       geneId,
				TranscriptId: transcriptId,
				IsCanonical:  isCanonical,
				IsLongest:    isLongest,
				Type:         transcriptType,
				//Children: make([]*GenomicFeature, 0, 10)
			}

			if currentGene != nil {
				if currentGene.Features == nil {
					currentGene.Features = make([]*GenomicFeature, 0, 10)
				}

				currentGene.Features = append(currentGene.Features, currentTranscript)
			} else {
				ret = append(ret, currentTranscript)
			}
		}

		// only add exon if we have a current transcript and it matches
		// the transcript id

		location, err := dna.NewStrandedLocation(chr, exonStart, exonEnd, strand)

		if err != nil {
			return nil, err
		}

		currentExon = &GenomicFeature{Id: gid,
			Location:     location,
			Feature:      ExonLevel,
			GeneSymbol:   geneSymbol,
			GeneId:       geneId,
			TranscriptId: transcriptId,
			ExonId:       exonId,
			ExonNumber:   exonNumber,
		}

		if currentTranscript != nil {
			if currentTranscript.Features == nil {
				currentTranscript.Features = make([]*GenomicFeature, 0, 10)
			}

			currentTranscript.Features = append(currentTranscript.Features, currentExon)
		} else if currentGene != nil {
			// normally add to transcript, but if no transcript, add to gene
			currentGene.Features = append(currentGene.Features, currentExon)
		} else {
			// no gene or transcript, just add to return list
			ret = append(ret, currentExon)
		}

	}

	log.Debug().Msgf("converted rows to %d features", len(ret))

	return ret, nil
}
