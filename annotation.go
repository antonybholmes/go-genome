package genome

import (
	"fmt"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys/log"
)

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

// type ClosestGene struct {
// 	string   *GenomicFeature `json:"feature"`
// 	PromLabel string          `json:"promLabel"`
// 	TssDist   int             `json:"tssDist"`
// }

type (
	AnnotationGene struct {
		GeneId     string `json:"geneId"`
		GeneSymbol string `json:"geneSymbol"`
		PromLabel  string `json:"promLabel"`
		Strand     string `json:"strand"`
		IsPromoter bool   `json:"-"`
		IsIntronic bool   `json:"-"`
		IsExon     bool   `json:"-"`
		TssDist    int    `json:"tssDist"`
	}

	GeneAnnotation struct {
		Location *dna.Location `json:"loc"`
		// GeneIds      string            `json:"geneIds"`
		// GeneSymbols  string            `json:"geneSymbols"`
		// GeneStrands  string            `json:"geneStrands"`
		// PromLabels   string            `json:"promLabels"`
		// TSSDists     string            `json:"tssDists"`
		// Locations    string            `json:"geneLocs"`
		WithinGenes  []*GenomicFeature `json:"withinGenes"`
		ClosestGenes []*GenomicFeature `json:"closestGenes"`
	}

	// type ByAbsD []GeneProm
	// func (a ByAbsD) Len() int           { return len(a) }
	// func (a ByAbsD) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
	// func (a ByAbsD) Less(i, j int) bool { return a[i].AbsD < a[j].AbsD }

	AnnotateDb struct {
		GeneDb    GeneDB
		TSSRegion *dna.PromoterRegion
		ClosestN  int8
	}
)

const (
	Na         string = "n/a"
	Promoter   string = "promoter"
	Exonic     string = "exonic"
	Intronic   string = "intronic"
	Intergenic string = "intergenic"
	Intragenic string = "intragenic"

	GroupSeparator   string = ","
	FeatureSeparator string = "|"
)

func NewAnnotateDb(genesdb GeneDB, tssRegion *dna.PromoterRegion, closestN int8) *AnnotateDb {
	return &AnnotateDb{
		GeneDb:    genesdb,
		TSSRegion: tssRegion,
		ClosestN:  closestN,
	}
}

func (annotateDb *AnnotateDb) Annotate(location *dna.Location, levels string) (*GeneAnnotation, error) {
	//mid := location.Mid()

	// extend search area to account  for promoter
	// let search_loc: Location = Location::new(
	//     &location.chr,
	//     location.Start - annotate.tssRegion.offset_5p.abs(),
	//     location.End + annotate.tssRegion.offset_5p.abs(),
	// )?;

	//log.Debug().Msgf("Annotating location %s %s", location, annotateDb.GeneDb.File)

	level := GeneLevel

	genesWithin, err := annotateDb.GeneDb.WithinGenesAndPromoter(
		location,
		level,
		annotateDb.TSSRegion,
	)

	if err != nil {
		log.Error().Msgf("Error within genes for location %s: %v", location, err)
		return nil, err
	}

	closestGenes, err := annotateDb.GeneDb.ClosestGenes(location, annotateDb.TSSRegion, annotateDb.ClosestN)

	if err != nil {
		log.Error().Msgf("Error closest genes for location %s: %v", location, err)
		return nil, err
	}

	annotation := GeneAnnotation{
		Location: location,
		// GeneIds:     strings.Join(geneIds, OUTPUT_FEATURE_SEP),
		// GeneSymbols: strings.Join(geneNames, OUTPUT_FEATURE_SEP),
		// GeneStrands: strings.Join(geneStrands, OUTPUT_FEATURE_SEP),
		// PromLabels:  strings.Join(promLabels, OUTPUT_FEATURE_SEP),

		// TSSDists:     strings.Join(tssDists, OUTPUT_FEATURE_SEP),
		// Locations:    strings.Join(featureLocations, OUTPUT_FEATURE_SEP),
		WithinGenes:  genesWithin,
		ClosestGenes: closestGenes,
	}

	return &annotation, nil
}

func (annotateDb *AnnotateDb) ClassifyFeature(location *dna.Location, feature *GenomicFeature) (string, error) {
	mid := location.Mid()
	var start int

	if feature.Location.Strand() == "+" {
		start = feature.Location.Start() - annotateDb.TSSRegion.Upstream()
	} else {
		start = feature.Location.Start()
	}

	var end int

	if feature.Location.Strand() == "-" {
		end = feature.Location.End() + annotateDb.TSSRegion.Upstream()
	} else {
		end = feature.Location.End()
	}

	// if location.Start > end || location.End < start {
	// 	return INTERGENIC
	// }

	isPromoter := (feature.Location.Strand() == "+" && mid >= start && mid <= feature.Location.Start()+annotateDb.TSSRegion.Downstream()) ||
		(feature.Location.Strand() == "-" && mid >= feature.Location.End()-annotateDb.TSSRegion.Downstream() && mid <= end)

	exons, err := annotateDb.GeneDb.InExon(location, feature.TranscriptId, annotateDb.TSSRegion)

	if err != nil {
		return "", err
	}

	isExon := len(exons) > 0

	isIntronic := mid >= feature.Location.Start() && mid <= feature.Location.End()

	return MakePromLabel(isPromoter, isExon, isIntronic), nil
}

func GeneWithStrandLabel(id string, strand string) string {
	return fmt.Sprintf("%s:%s", id, strand)
}

func MakePromLabel(isPromoter bool, isExon bool, isIntragenic bool) string {
	labels := make([]string, 0, 3)

	if isPromoter {
		labels = append(labels, Promoter)
	}

	// favor exonic over intronic
	if isExon {
		labels = append(labels, Exonic)
	} else {
		if isIntragenic {
			labels = append(labels, Intronic)
		}
	}

	if isIntragenic {
		labels = append(labels, Intragenic)
	} else {
		labels = append(labels, Intergenic)
	}

	//slices.Sort(labels)

	return strings.Join(labels, GroupSeparator)
}
