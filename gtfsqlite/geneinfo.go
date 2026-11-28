package gtfsqlite

import (
	"database/sql"
	"fmt"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-genome"
	"github.com/antonybholmes/go-sys/collections"
)

const (
	// GeneInfoSql = `SELECT DISTINCT
	// 	g.id,
	// 	c.name AS chr,
	// 	g.start,
	// 	g.end,
	// 	g.strand,
	// 	g.gene_id,
	// 	g.symbol,
	// 	bt.name AS biotype
	// FROM genes as g
	// JOIN biotypes AS bt ON g.biotype_id = bt.id
	// JOIN chromosomes AS c ON g.chr_id = c.id
	// WHERE (g.symbol LIKE :symbol OR g.gene_id = :q)
	// ORDER BY g.symbol
	// LIMIT :n`

	// Search for transcripts by gene symbol, gene ID, or transcript ID
	// TranscriptInfoSql = `SELECT *
	// 	FROM(
	// 		SELECT DISTINCT
	// 			g.id,
	// 			c.name AS chr,
	// 			g.start,
	// 			g.end,
	// 			g.strand,
	// 			g.gene_id,
	// 			g.symbol,
	// 			bt.name AS biotype,
	// 			t.transcript_id,
	// 			t.start,
	// 			t.end,
	// 			t.is_canonical,
	// 			t.is_longest,
	// 			ROW_NUMBER() OVER (PARTITION BY g.gene_id ORDER BY g.gene_id) AS rank
	// 		FROM genes as g
	// 		JOIN transcripts AS t ON g.id = t.gene_id
	// 		JOIN biotypes AS bt ON g.biotype_id = bt.id
	// 		JOIN chromosomes AS c ON g.chr_id = c.id
	// 		WHERE (g.symbol LIKE :q OR g.gene_id LIKE :q OR t.transcript_id LIKE :q)
	// 			<<CANONICAL>>
	// 		ORDER BY g.symbol
	// 	) features
	// 	WHERE features.rank < :n`

	GeneInfoSql = `SELECT *
		FROM(
			SELECT DISTINCT
				g.id,
				c.name AS chr,
				g.start,
				g.end,
				g.strand,
				g.gene_id,
				g.symbol,
				bt.name AS biotype,
				ROW_NUMBER() OVER (PARTITION BY g.gene_id ORDER BY g.gene_id) AS rank
			FROM genes as g
			JOIN biotypes AS bt ON g.biotype_id = bt.id
			JOIN chromosomes AS c ON g.chr_id = c.id
			WHERE (g.symbol LIKE :q OR g.gene_id LIKE :q)
				<<CANONICAL>>
			ORDER BY g.symbol
		) features
		WHERE features.rank < :n`

	GeneTranscriptInfoSql = `SELECT *
		FROM(
			SELECT DISTINCT
				g.id,
				c.name AS chr,
				g.start,
				g.end,
				g.strand,
				g.gene_id,
				g.symbol,
				bt.name AS biotype,
				t.transcript_id,
				t.start,
				t.end,
				t.is_canonical,
				t.is_longest,
				ROW_NUMBER() OVER (PARTITION BY g.gene_id ORDER BY g.gene_id) AS rank
			FROM genes as g
			JOIN transcripts AS t ON g.id = t.gene_id
			JOIN biotypes AS bt ON g.biotype_id = bt.id
			JOIN chromosomes AS c ON g.chr_id = c.id
			WHERE (g.symbol LIKE :q OR g.gene_id LIKE :q OR t.transcript_id LIKE :q)
				<<CANONICAL>>
			ORDER BY g.symbol
		) features
		WHERE features.rank < :n`

	GeneTranscriptExonInfoSql = `SELECT *
		FROM(
			SELECT DISTINCT
				g.id,
				c.name AS chr,
				g.start,
				g.end,
				g.strand,
				g.gene_id,
				g.symbol,
				bt.name AS biotype,
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
				ROW_NUMBER() OVER (PARTITION BY g.gene_id ORDER BY g.gene_id) AS rank
			FROM genes as g
			JOIN transcripts AS t ON g.id = t.gene_id
			JOIN features AS f ON f.transcript_id = t.id
			JOIN feature_types AS ft ON f.feature_type_id = ft.id
			JOIN biotypes AS bt ON g.biotype_id = bt.id
			JOIN chromosomes AS c ON g.chr_id = c.id
			WHERE (g.symbol LIKE :q OR g.gene_id LIKE :q OR t.transcript_id LIKE :q OR f.exon_id LIKE :q)
				<<CANONICAL>>
			ORDER BY g.symbol
		) features
		WHERE features.rank < :n`

	// ExonInfoSql = `SELECT *` +
	// 	` FROM (` + BasicLocationSql +
	// 	` WHERE (g.symbol LIKE :q OR g.gene_id LIKE :q OR t.transcript_id LIKE :q OR f.exon_id LIKE :q)
	// 			<<CANONICAL>>
	// 		ORDER BY g.symbol
	// 	) features
	// 	WHERE features.rank < :n`
)

func (gdb *GtfSqliteDB) SearchByName(search string,
	levels collections.Set[genome.FeatureType],
	canonicalMode bool,
	n int16) ([]*genome.GenomicFeature, error) {
	n = max(1, min(n, MaxGeneInfoResults))

	//log.Debug().Msgf("SearchForGeneByName for gene %s level %s fuzzy %v canonical %v geneType %s n %d",
	//	search, level, fuzzy, canonical, geneType, n)

	// case insensitive search
	search = strings.ToLower(search)

	if len(search) < 2 || strings.Contains(search, "chr:") {
		return nil, fmt.Errorf("%s is an invalid search term", search)
	}

	var sqlStmt string

	// use shorter SQL statements depending on the feature level
	if levels.Has("exon") {
		sqlStmt = GeneTranscriptExonInfoSql
	} else if levels.Has("transcript") {
		sqlStmt = GeneTranscriptInfoSql
	} else {
		sqlStmt = GeneInfoSql
	}

	if canonicalMode {
		sqlStmt = strings.Replace(sqlStmt, "<<CANONICAL>>", "AND t.is_canonical = 1", 1)
	} else {
		sqlStmt = strings.Replace(sqlStmt, "<<CANONICAL>>", "", 1)
	}

	rows, err := gdb.db.Query(sqlStmt,
		sql.Named("q", search),
		sql.Named("n", n))

	if err != nil {

		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	return rowsToGeneInfoRecords(rows, levels, canonicalMode)

	// switch level {
	// case genome.FeatureTypeTranscript:
	// 	return gdb.searchTranscripts(search,
	// 		canonicalMode,
	// 		false,
	// 		n)

	// case genome.FeatureTypeExon:
	// 	return gdb.searchExons(search,
	// 		canonicalMode,
	// 		n)
	// default:
	// 	return gdb.searchGenes(search, n)
	// }

}

// Searching for exons or transcripts uses essentially
// the same pipeline so combine into one method.
// func (gdb *GtfSqliteDB) searchTranscripts(search string,
// 	canonicalMode bool,
// 	addExonInfo bool,
// 	n int16) ([]*genome.GenomicFeature, error) {
// 	n = max(1, min(n, MaxGeneInfoResults))

// 	var sqlStmt string

// 	if addExonInfo {
// 		sqlStmt = ExonInfoSql
// 	} else {
// 		sqlStmt = TranscriptInfoSql
// 	}

// 	if canonicalMode {
// 		sqlStmt = strings.Replace(sqlStmt, "<<CANONICAL>>", "AND t.is_canonical = 1", 1)
// 	} else {
// 		sqlStmt = strings.Replace(sqlStmt, "<<CANONICAL>>", "", 1)
// 	}

// 	//log.Debug().Msgf("SQL: %s %s %d", sqlStmt, search, n)

// 	rows, err := gdb.db.Query(sqlStmt,
// 		sql.Named("q", search),
// 		sql.Named("n", n))

// 	if err != nil {
// 		return nil, err //fmt.Errorf("there was an error with the database query")
// 	}

// 	defer rows.Close()

// 	if addExonInfo {
// 		return exonsToGeneInfoRecords(rows, canonicalMode)
// 	} else {
// 		return transcriptsToGeneInfoRecords(rows, canonicalMode)
// 	}

// }

// func (gdb *GtfSqliteDB) searchExons(search string,
// 	canonicalMode bool,
// 	n int16) ([]*genome.GenomicFeature, error) {
// 	n = max(1, min(n, MaxGeneInfoResults))

// 	var sqlStmt string

// 	if canonicalMode {
// 		sqlStmt = strings.Replace(ExonInfoSql, "<<CANONICAL>>", "AND t.is_canonical = 1", 1)
// 	} else {
// 		sqlStmt = strings.Replace(ExonInfoSql, "<<CANONICAL>>", "", 1)
// 	}

// 	//log.Debug().Msgf("SQL: %s %s %d", sqlStmt, search, n)

// 	rows, err := gdb.db.Query(sqlStmt,
// 		sql.Named("q", search),
// 		sql.Named("n", n))

// 	if err != nil {
// 		return nil, err //fmt.Errorf("there was an error with the database query")
// 	}

// 	defer rows.Close()

// 	return exonsToGeneInfoRecords(rows, canonicalMode)
// }

// func (gdb *GtfSqliteDB) searchGenes(search string,
// 	n int16) ([]*genome.GenomicFeature, error) {
// 	n = max(1, min(n, MaxGeneInfoResults))

// 	//log.Debug().Msgf("searching for gene: %s, n: %d, SQL: %s", search, n, GeneInfoSql)

// 	rows, err := gdb.db.Query(GeneInfoSql,
// 		sql.Named("q", search),
// 		sql.Named("symbol", search+"%"),
// 		sql.Named("n", n))

// 	if err != nil {
// 		return nil, err //fmt.Errorf("there was an error with the database query")
// 	}

// 	defer rows.Close()

// 	ret, err := genesToGeneInfoRecords(rows)

// 	if err != nil {
// 		return nil, err
// 	}

// 	return ret, nil
// }

// complete records for genes -> transcripts -> exons
// func genesToGeneInfoRecords(rows *sql.Rows) ([]*genome.GenomicFeature, error) {
// 	var gid int
// 	var chr string
// 	var geneStart int
// 	var geneEnd int
// 	var strand string
// 	var geneSymbol string
// 	var geneId string
// 	var geneType string

// 	// 10 seems a reasonable guess for the number of features we might see, just
// 	// to reduce slice reallocation

// 	var ret = make([]*genome.GenomicFeature, 0, 10)

// 	for rows.Next() {
// 		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &exonId)
// 		err := rows.Scan(&gid,
// 			&chr,
// 			&geneStart,
// 			&geneEnd,
// 			&strand,
// 			&geneId,
// 			&geneSymbol,
// 			&geneType,
// 		)

// 		if err != nil {
// 			log.Error().Msgf("error reading gene info rows %s", err)
// 			return nil, err //fmt.Errorf("there was an error with the database records")
// 		}

// 		// only add a new gene if we don't already have it. We
// 		// assume the rows are ordered by gene id hence if the
// 		// id changes, we are processing a set of rows for a new gene

// 		location, err := dna.NewStrandedLocation(chr, geneStart, geneEnd, strand)

// 		if err != nil {
// 			return nil, err
// 		}

// 		currentGene := &genome.GenomicFeature{Id: gid,
// 			//Location: location,
// 			Type:   genome.FeatureTypeGene,
// 			Symbol: geneSymbol,
// 			GeneId: geneId,
// 			//Strand:   strand,
// 			Biotype: geneType,
// 			//Children: make([]*GenomicFeature, 0, 10)
// 		}
// 		genome.AddLocationToFeature(currentGene, location)

// 		ret = append(ret, currentGene)

// 	}

// 	//log.Debug().Msgf("converted rows to %d features", len(ret))

// 	return ret, nil
// }

// complete records for genes -> transcripts -> exons
// func transcriptsToGeneInfoRecords(rows *sql.Rows, canonicalMode bool) ([]*genome.GenomicFeature, error) {
// 	var gid int
// 	var chr string
// 	var geneStart int
// 	var geneEnd int
// 	var strand string
// 	var geneSymbol string
// 	var geneId string
// 	var biotype string

// 	var transcriptId string
// 	//var transcriptType string
// 	var transcriptStart int
// 	var transcriptEnd int

// 	var isCanonical bool
// 	var isLongest bool

// 	var rank int

// 	// 10 seems a reasonable guess for the number of features we might see, just
// 	// to reduce slice reallocation

// 	var currentGene *genome.GenomicFeature
// 	var currentTranscript *genome.GenomicFeature

// 	var ret = make([]*genome.GenomicFeature, 0, 10)

// 	for rows.Next() {
// 		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &exonId)
// 		err := rows.Scan(&gid,
// 			&chr,
// 			&geneStart,
// 			&geneEnd,
// 			&strand,
// 			&geneId,
// 			&geneSymbol,
// 			&biotype,
// 			&transcriptId,
// 			&transcriptStart,
// 			&transcriptEnd,
// 			&isCanonical,
// 			&isLongest,
// 			&rank,
// 		)

// 		if err != nil {
// 			log.Error().Msgf("error reading gene info rows %s", err)
// 			return nil, err //fmt.Errorf("there was an error with the database records")
// 		}

// 		// only add a new gene if we don't already have it. We
// 		// assume the rows are ordered by gene id hence if the
// 		// id changes, we are processing a set of rows for a new gene
// 		if currentGene == nil || currentGene.GeneId != geneId {
// 			location, err := dna.NewStrandedLocation(chr, geneStart, geneEnd, strand)

// 			if err != nil {
// 				return nil, err
// 			}

// 			currentGene = &genome.GenomicFeature{Id: gid,
// 				//Location: location,
// 				Type:   genome.FeatureTypeGene,
// 				Symbol: geneSymbol,
// 				GeneId: geneId,
// 				//Strand:   strand,
// 				Biotype: biotype,
// 				//Children: make([]*GenomicFeature, 0, 10)
// 			}

// 			genome.AddLocationToFeature(currentGene, location)

// 			ret = append(ret, currentGene)
// 		}

// 		// gene inherits properties from transcripts and these become true
// 		// if any transcript has them true
// 		if currentGene != nil {
// 			currentGene.IsCanonical = currentGene.IsCanonical || isCanonical
// 			//currentGene.IsLongest = currentGene.IsLongest || isLongest

// 		}

// 		//if canonical mode only add if transcript is canonical
// 		// also only add if we have a current gene
// 		// also only add if we don't already have this transcript
// 		if (currentTranscript == nil || currentTranscript.TranscriptId != transcriptId) &&
// 			(!canonicalMode || isCanonical) {

// 			location, err := dna.NewStrandedLocation(chr, transcriptStart, transcriptEnd, strand)

// 			if err != nil {
// 				return nil, err
// 			}

// 			// set the properties that will not change for a transcript
// 			currentTranscript = &genome.GenomicFeature{Id: gid,
// 				//Location: location,
// 				//Strand:       strand,
// 				Type:         genome.FeatureTypeTranscript,
// 				Symbol:       geneSymbol,
// 				GeneId:       geneId,
// 				TranscriptId: transcriptId,
// 				IsCanonical:  isCanonical,
// 				//IsLongest:   isLongest,
// 				//Biotype:      transcriptType,
// 				//Children: make([]*GenomicFeature, 0, 10)
// 			}

// 			genome.AddLocationToFeature(currentTranscript, location)

// 			if currentGene != nil {
// 				if currentGene.Children == nil {
// 					currentGene.Children = make([]*genome.GenomicFeature, 0, 10)
// 				}

// 				currentGene.Children = append(currentGene.Children, currentTranscript)
// 			} else {
// 				ret = append(ret, currentTranscript)
// 			}
// 		}
// 	}

// 	return ret, nil
// }

// complete records for genes -> transcripts -> exons
func rowsToGeneInfoRecords(rows *sql.Rows,
	levels collections.Set[genome.FeatureType],
	canonicalMode bool) ([]*genome.GenomicFeature, error) {

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation

	var currentGene *genome.GenomicFeature
	var currentTranscript *genome.GenomicFeature
	var currentExon *genome.GenomicFeature

	var ret = make([]*genome.GenomicFeature, 0, 10)
	var err error

	record := GeneInfoRecord{}

	for rows.Next() {
		//err := geneRows.Scan(&id, &level, &chr, &start, &end, &strand, &geneId, &geneSymbol, &transcriptId, &exonId)

		// The query result will have different columns depending on the level of detail requested.
		if levels.Has(genome.FeatureTypeExon) {
			err = rows.Scan(&record.Id,
				&record.Chr,
				&record.GeneStart,
				&record.GeneEnd,
				&record.Strand,
				&record.GeneId,
				&record.GeneSymbol,
				&record.Biotype,
				&record.TranscriptId,
				&record.TranscriptStart,
				&record.TranscriptEnd,
				&record.IsCanonical,
				&record.IsLongest,
				&record.FeatureType,
				&record.FeatureStart,
				&record.FeatureEnd,
				&record.ExonId,
				&record.ExonNumber,
				&record.Rank,
			)
		} else if levels.Has(genome.FeatureTypeTranscript) {
			err = rows.Scan(&record.Id,
				&record.Chr,
				&record.GeneStart,
				&record.GeneEnd,
				&record.Strand,
				&record.GeneId,
				&record.GeneSymbol,
				&record.Biotype,
				&record.TranscriptId,
				&record.TranscriptStart,
				&record.TranscriptEnd,
				&record.IsCanonical,
				&record.IsLongest,
				&record.Rank,
			)
		} else {
			err = rows.Scan(&record.Id,
				&record.Chr,
				&record.GeneStart,
				&record.GeneEnd,
				&record.Strand,
				&record.GeneId,
				&record.GeneSymbol,
				&record.Biotype,
				&record.Rank,
			)
		}

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		// only add a new gene if we don't already have it. We
		// assume the rows are ordered by gene id hence if the
		// id changes, we are processing a set of rows for a new gene
		if levels.Has(genome.FeatureTypeGene) && (currentGene == nil || currentGene.GeneId != record.GeneId) {
			location, err := dna.NewStrandedLocation(record.Chr, record.GeneStart, record.GeneEnd, record.Strand)

			if err != nil {
				return nil, err
			}

			currentGene = &genome.GenomicFeature{Id: record.Id,
				Type:    genome.FeatureTypeGene,
				GeneId:  record.GeneId,
				Symbol:  record.GeneSymbol,
				Biotype: record.Biotype,
			}

			genome.AddLocationToFeature(currentGene, location)
			ret = append(ret, currentGene)
		}

		// gene inherits properties from transcripts and these become true
		// if any transcript has them true
		if currentGene != nil {
			currentGene.IsCanonical = currentGene.IsCanonical || record.IsCanonical
			//currentGene.IsLongest = currentGene.IsLongest || isLongest

		}

		//if canonical mode only add if transcript is canonical
		// also only add if we have a current gene
		// also only add if we don't already have this transcript
		if levels.Has(genome.FeatureTypeTranscript) &&
			(currentTranscript == nil || currentTranscript.TranscriptId != record.TranscriptId) &&
			(!canonicalMode || record.IsCanonical) {

			location, err := dna.NewStrandedLocation(record.Chr, record.TranscriptStart, record.TranscriptEnd, record.Strand)

			if err != nil {
				return nil, err
			}

			// set the properties that will not change for a transcript
			currentTranscript = &genome.GenomicFeature{Id: record.Id,
				Type: genome.FeatureTypeTranscript,
				//Symbol:       record.GeneSymbol,
				//GeneId:       record.GeneId,
				TranscriptId: record.TranscriptId,
				IsCanonical:  record.IsCanonical,
				//IsLongest:   isLongest,
				//Biotype: record.Biotype,
			}
			genome.AddLocationToFeature(currentTranscript, location)

			// if we don't have gene level features, set the gene id and symbol on the transcript
			if !levels.Has(genome.FeatureTypeGene) {
				currentTranscript.GeneId = record.GeneId
				currentTranscript.Symbol = record.GeneSymbol
				currentTranscript.Biotype = record.Biotype
			}

			if currentGene != nil {
				if currentGene.Children == nil {
					currentGene.Children = make([]*genome.GenomicFeature, 0, 10)
				}

				currentGene.Children = append(currentGene.Children, currentTranscript)
			} else {
				ret = append(ret, currentTranscript)
			}
		}

		if levels.Has(record.FeatureType) {
			// only add exon if we have a current transcript and it matches
			// the transcript id

			location, err := dna.NewStrandedLocation(record.Chr, record.FeatureStart, record.FeatureEnd, record.Strand)

			if err != nil {
				return nil, err
			}

			currentExon = &genome.GenomicFeature{Id: record.Id,
				Type: genome.FeatureTypeExon,
				//Symbol:       record.GeneSymbol,
				//GeneId:       record.GeneId,
				//TranscriptId: record.TranscriptId,
				ExonId:     record.ExonId,
				ExonNumber: record.ExonNumber,
			}
			genome.AddLocationToFeature(currentExon, location)

			if !levels.Has(genome.FeatureTypeGene) {
				currentExon.GeneId = record.GeneId
				currentExon.Symbol = record.GeneSymbol
				currentExon.Biotype = record.Biotype
			}

			if !levels.Has(genome.FeatureTypeTranscript) {
				currentExon.TranscriptId = record.TranscriptId
			}

			if currentTranscript != nil {
				if currentTranscript.Children == nil {
					currentTranscript.Children = make([]*genome.GenomicFeature, 0, 10)
				}

				currentTranscript.Children = append(currentTranscript.Children, currentExon)
			} else if currentGene != nil {
				// normally add to transcript, but if no transcript, add to gene
				currentGene.Children = append(currentGene.Children, currentExon)
			} else {
				// no gene or transcript, just add to return list
				ret = append(ret, currentExon)
			}
		}
	}

	if err := rows.Err(); err != nil {
		return nil, err
	}

	//log.Debug().Msgf("converted rows to %d features", len(ret))

	return ret, nil
}
