package gtftabix

import (
	"path/filepath"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-genome"
	"github.com/antonybholmes/go-sys/log"
	basemath "github.com/antonybholmes/go-sys/math"
)

type (
	GtfTabixDB struct {
		annotation *genome.Annotation
		file       string
	}
)

func NewGtfTabixDB(dir string, annotation *genome.Annotation) *GtfTabixDB {
	file := filepath.Join(dir, annotation.Url)

	log.Debug().Msgf("opening gene database %s", file)

	return &GtfTabixDB{annotation: annotation, file: file}
}

func (gdb *GtfTabixDB) Close() error {
	return nil
}

func (gdb *GtfTabixDB) OverlappingGenes(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion,
	canonicalMode bool,
	annotationMode bool,
	biotypeFilter string) ([]*genome.GenomicFeature, error) {

	roots, err := QueryRegion(gdb.file, location)

	if err != nil {
		return nil, err
	}

	features, err := toGenomicFeature(location,
		prom,
		roots,
		levels,
		canonicalMode,
		annotationMode)

	if err != nil {
		return nil, err
	}

	return features, nil
}

func (gdb *GtfTabixDB) ClosestGenes(location *dna.Location,
	prom *dna.PromoterRegion,
	closestN int8,
	useOfficialGenes bool) ([]*genome.GenomicFeature, error) {
	return nil, nil
}

// Returns the genes/transcripts with exons overlapping the given location
func (gdb *GtfTabixDB) InExon(location *dna.Location, transcriptId string, prom *dna.PromoterRegion) (bool, error) {
	roots, err := QueryRegion(gdb.file, location)

	if err != nil {
		return false, err
	}

	for _, root := range roots {
		root.FindFirst(func(feat *Feature) (string, bool) {
			if feat.Type == Exon && feat.TranscriptID() == transcriptId {
				exonLoc, err := dna.NewLocation(feat.Chr(), feat.Start, feat.End)
				if err != nil {
					return "", false
				}

				if exonLoc.Overlaps(location) {
					return "exon", true
				}

			}
			return "", false
		})

	}

	return false, nil
}

func (gdb *GtfTabixDB) IntragenicFeatures(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion,
	useOfficialGenes bool) ([]*genome.GenomicFeature, error) {

	maxInc := basemath.Max(prom.FivePrime(), prom.ThreePrime())

	// expand region for promoter regions if provided
	loc, err := dna.NewLocation(location.Chr(), location.Start()-maxInc, location.End()+maxInc)

	if err != nil {
		return nil, err
	}

	roots, err := QueryRegion(gdb.file, loc)

	if err != nil {
		return nil, err
	}

	features, err := toGenomicFeature(location,
		prom,
		roots,
		levels,
		false,
		true)

	if err != nil {
		return nil, err
	}

	return features, nil
}

// func (gdb *GtfTabixDB) WithinGenes(location *dna.Location,
// 	feature string,
// 	prom *dna.PromoterRegion) (*genome.GenomicSearchResults, error) {

// 	roots, err := QueryRegion(gdb.file, location)

// 	if err != nil {
// 		return nil, err
// 	}

// 	return rowsToSearchResults(location, prom, roots, feature, false, true)
// }

func (gdb *GtfTabixDB) SearchByName(search string,
	level string,
	canonicalMode bool,
	n int16) ([]*genome.GenomicFeature, error) {
	return nil, nil
}

// Convert gtf features from a file into standardized genomic features
func toGenomicFeature(location *dna.Location,
	prom *dna.PromoterRegion,
	feats []*Feature,
	levels string,
	canonicalMode bool,
	annotationMode bool) ([]*genome.GenomicFeature, error) {
	var currentGene *genome.GenomicFeature
	var currentTranscript *genome.GenomicFeature

	var currentFeature *genome.GenomicFeature

	var ret = make([]*genome.GenomicFeature, 0, 10)

	gid := 1

	for _, feat := range feats {
		chr := feat.SeqID
		geneId := feat.GeneID()
		transcriptId := feat.TranscriptID()
		exonId := feat.ExonID()
		exonNumber := feat.ExonNumber()
		start := feat.Start
		end := feat.End

		strand := feat.Strand
		geneSymbol := feat.GeneName()
		featureType := feat.Type

		biotype := feat.Biotype()

		// only add a new gene if we don't already have it. We
		// assume the rows are ordered by gene id hence if the
		// id changes, we are processing a set of rows for a new gene
		if strings.Contains(levels, "gene") && (currentGene == nil || currentGene.GeneId != geneId) {
			location, err := dna.NewStrandedLocation(chr, start, end, strand)

			if err != nil {
				return nil, err
			}

			currentGene = &genome.GenomicFeature{Id: gid,
				Location: location,
				Type:     feat.Type,
				Symbol:   geneSymbol,
				GeneId:   geneId,
				Strand:   strand,
				Biotype:  biotype,
			}

			ret = append(ret, currentGene)

			if annotationMode {
				currentGene.Label = genome.MakePromLabel(currentGene.InPromoter, currentGene.InExon, currentGene.IsIntragenic)
			}

			gid++
		}

		//if canonical mode only add if transcript is canonical
		// also only add if we have a current gene
		// also only add if we don't already have this transcript
		if strings.Contains(levels, "transcript") &&
			(currentTranscript == nil || currentTranscript.Transcript != transcriptId) {

			isCanonical := feat.IsCanonical()

			if canonicalMode && !isCanonical {
				continue
			}

			inPromoter, err := feat.InPromoter(location, prom)
			if err != nil {
				return nil, err
			}

			isIntragenic, err := feat.IsIntragenic(location)
			if err != nil {
				return nil, err
			}

			location, err := dna.NewStrandedLocation(chr, start, end, strand)

			if err != nil {
				return nil, err
			}

			tssDist, err := feat.TssDist(location)

			if err != nil {
				return nil, err
			}

			// set the properties that will not change for a transcript
			// isLongest cannot be set from tabix gtf file since this
			// information is not provided in the GTF rows
			currentTranscript = &genome.GenomicFeature{Id: gid,
				Location: location,
				//Strand:       strand,
				Type:       genome.TranscriptLevel,
				Symbol:     geneSymbol,
				GeneId:     geneId,
				Transcript: transcriptId,
				//Biotype:      transcriptBiotype,
				IsCanonical:  isCanonical,
				IsIntragenic: isIntragenic,
				InPromoter:   inPromoter,
				//IsLongest:   isLongest,
				TssDist: tssDist,
			}

			if currentGene != nil {
				currentGene.InPromoter = currentGene.InPromoter || inPromoter
				currentGene.IsIntragenic = currentGene.IsIntragenic || isIntragenic

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

			if annotationMode {
				currentTranscript.Label = genome.MakePromLabel(currentTranscript.InPromoter, currentTranscript.InExon, currentTranscript.IsIntragenic)
			}
		}

		if currentTranscript != nil {
			// these properties may be updated as we see more rows and we find
			// we are in an exon

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
			location, err := dna.NewStrandedLocation(chr, start, end, strand)

			if err != nil {
				return nil, err
			}

			inExon := feat.InExon(location)

			currentFeature = &genome.GenomicFeature{Id: gid,
				Location:   location,
				Type:       featureType,
				Symbol:     geneSymbol,
				GeneId:     geneId,
				Transcript: transcriptId,
				Exon:       exonId,
				ExonNumber: exonNumber,
				InExon:     inExon,
			}

			if currentTranscript != nil {
				currentFeature.InPromoter = currentTranscript.InPromoter
				currentFeature.IsIntragenic = currentTranscript.IsIntragenic

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

			if annotationMode {
				currentFeature.Label = genome.MakePromLabel(currentFeature.InPromoter, inExon, currentFeature.IsIntragenic)

				if currentGene != nil {
					currentGene.InExon = currentGene.InExon || inExon
					currentGene.Label = currentFeature.Label
				}

				if currentTranscript != nil {
					currentTranscript.InExon = currentTranscript.InExon || inExon
					currentTranscript.Label = currentFeature.Label
				}
			}
		}
	}

	return ret, nil
}

func rowsToSearchResults(location *dna.Location,
	prom *dna.PromoterRegion,
	feats []*Feature,
	levels string,
	canonicalMode bool,
	annotationMode bool) (*genome.GenomicSearchResults, error) {

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	features, err := toGenomicFeature(location,
		prom,
		feats,
		levels,
		canonicalMode,
		annotationMode)

	if err != nil {
		log.Error().Msgf("error converting rows to features %s", err)
		return nil, err
	}

	return genome.FeaturesToSearchResults(location, levels, features, annotationMode)
}
