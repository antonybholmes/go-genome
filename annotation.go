package genome

import (
	"fmt"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys/log"
)

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

// type ClosestGene struct {
// 	Feature   *GenomicFeature `json:"feature"`
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
		ClosestGenes []*AnnotationGene `json:"closestGenes"`
	}

	// type ByAbsD []GeneProm
	// func (a ByAbsD) Len() int           { return len(a) }
	// func (a ByAbsD) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
	// func (a ByAbsD) Less(i, j int) bool { return a[i].AbsD < a[j].AbsD }

	AnnotateDb struct {
		GeneDb    GeneDB
		TSSRegion *dna.PromoterRegion
		ClosestN  uint8
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

func NewAnnotateDb(genesdb GeneDB, tssRegion *dna.PromoterRegion, closestN uint8) *AnnotateDb {
	return &AnnotateDb{
		GeneDb:    genesdb,
		TSSRegion: tssRegion,
		ClosestN:  closestN,
	}
}

func (annotateDb *AnnotateDb) Annotate(location *dna.Location, levels Level) (*GeneAnnotation, error) {
	//mid := location.Mid()

	// extend search area to account  for promoter
	// let search_loc: Location = Location::new(
	//     &location.chr,
	//     location.Start - annotate.tssRegion.offset_5p.abs(),
	//     location.End + annotate.tssRegion.offset_5p.abs(),
	// )?;

	//log.Debug().Msgf("Annotating location %s %s", location, annotateDb.GeneDb.File)

	genesWithin, err := annotateDb.GeneDb.WithinGenesAndPromoter(
		location,
		levels,
		annotateDb.TSSRegion.Offset5P(),
		annotateDb.TSSRegion.Offset3P(),
	)

	log.Debug().Msgf("Found %d genes within for location %s", len(genesWithin.Features), location)

	if err != nil {
		log.Error().Msgf("Error within genes for location %s: %v", location, err)
		return nil, err
	}

	// we need the unique ids to symbols
	//idMap := make(map[string]string)
	// withinMap := make(map[string]AnnotationGene)

	// //let mut distMap: HashMap<&str, bool> = HashMap::new();

	// for _, gene := range genesWithin.Features {
	// 	id := gene.GeneId

	// 	//idMap[id] = gene.GeneName

	// 	//let labels = annotate.classify_location(location, gene);

	// 	exons, err := annotateDb.GeneDb.InExon(location, gene.TranscriptId)

	// 	if err != nil {
	// 		return nil, err
	// 	}

	// 	isExon := len(exons.Features) > 0

	// 	isPromoter := (gene.Location.Strand == "+" && mid >= gene.Location.Start-annotateDb.TSSRegion.Offset5P() &&
	// 		mid <= gene.Location.Start+annotateDb.TSSRegion.Offset3P()) ||
	// 		(gene.Location.Strand == "-" && mid >= gene.Location.End-annotateDb.TSSRegion.Offset3P() &&
	// 			mid <= gene.Location.End+annotateDb.TSSRegion.Offset5P())

	// 	isIntronic := mid >= gene.Location.Start && mid <= gene.Location.End

	// 	prom, ok := withinMap[id]

	// 	if ok {
	// 		// we test all transcripts to see if we can classify the gene
	// 		// hence we OR the tests since we want to make as many of
	// 		// them true as possible
	// 		prom.IsIntronic = prom.IsIntronic || isIntronic
	// 		prom.IsPromoter = prom.IsPromoter || isPromoter
	// 		prom.IsExon = prom.IsExon || len(exons.Features) > 0
	// 		prom.PromLabel = MakePromLabel(isPromoter, isExon, isIntronic)
	// 	} else {
	// 		withinMap[id] = AnnotationGene{
	// 			GeneId:     gene.GeneId,
	// 			GeneSymbol: gene.GeneSymbol,
	// 			PromLabel:  MakePromLabel(isPromoter, isExon, isIntronic),
	// 			Strand:     gene.Location.Strand,
	// 			//IsPromoter: isPromoter,
	// 			//IsIntronic: isIntronic,
	// 			//IsExon:     isExon,
	// 			TssDist: gene.TssDist,
	// 		}
	// 	}

	// }

	// withinGenes := make([]*AnnotationGene, 0, len(withinMap))

	// for _, gene := range withinMap {
	// 	withinGenes = append(withinGenes, &gene)
	// }

	// slices.SortFunc(withinGenes, func(a, b *AnnotationGene) int {
	// 	return int(a.TssDist) - int(b.TssDist)
	// })

	// closestAnnotations := make([]*AnnotationGene, 0, annotateDb.ClosestN)

	// if false {
	// 	closestGenes, err := annotateDb.GeneDb.ClosestGenes(location, uint16(annotateDb.ClosestN))

	// 	if err != nil {
	// 		log.Debug().Msgf("Error closest genes for location %s: %v", location, err)
	// 		return nil, err
	// 	}

	// 	for _, cg := range closestGenes {

	// 		exons, err := annotateDb.GeneDb.InExon(location, cg.TranscriptId)

	// 		if err != nil {
	// 			return nil, err
	// 		}

	// 		isExon := len(exons.Features) > 0

	// 		isPromoter := (cg.Location.Strand == "+" && mid >= cg.Location.Start-annotateDb.TSSRegion.Offset5P() &&
	// 			mid <= cg.Location.Start+annotateDb.TSSRegion.Offset3P()) ||
	// 			(cg.Location.Strand == "-" && mid >= cg.Location.End-annotateDb.TSSRegion.Offset3P() &&
	// 				mid <= cg.Location.End+annotateDb.TSSRegion.Offset5P())

	// 		isIntronic := mid >= cg.Location.Start && mid <= cg.Location.End

	// 		closestAnnotations = append(closestAnnotations, &AnnotationGene{
	// 			GeneId:     cg.GeneId,
	// 			GeneSymbol: cg.GeneSymbol,
	// 			PromLabel:  MakePromLabel(isPromoter, isExon, isIntronic),
	// 			IsPromoter: isPromoter,
	// 			IsIntronic: isIntronic,
	// 			IsExon:     isExon,
	// 			TssDist:    cg.TssDist,
	// 			Strand:     cg.Location.Strand,
	// 		})

	// 	}
	// }

	annotation := GeneAnnotation{
		Location: location,
		// GeneIds:     strings.Join(geneIds, OUTPUT_FEATURE_SEP),
		// GeneSymbols: strings.Join(geneNames, OUTPUT_FEATURE_SEP),
		// GeneStrands: strings.Join(geneStrands, OUTPUT_FEATURE_SEP),
		// PromLabels:  strings.Join(promLabels, OUTPUT_FEATURE_SEP),

		// TSSDists:     strings.Join(tssDists, OUTPUT_FEATURE_SEP),
		// Locations:    strings.Join(featureLocations, OUTPUT_FEATURE_SEP),
		WithinGenes: genesWithin.Features,
		//ClosestGenes: closestAnnotations,
	}

	return &annotation, nil
}

func (annotateDb *AnnotateDb) ClassifyFeature(location *dna.Location, feature *GenomicFeature) (string, error) {
	mid := location.Mid()
	var start uint

	if feature.Location.Strand == "+" {
		start = feature.Location.Start - annotateDb.TSSRegion.Offset5P()
	} else {
		start = feature.Location.Start
	}

	var end uint

	if feature.Location.Strand == "-" {
		end = feature.Location.End + annotateDb.TSSRegion.Offset5P()
	} else {
		end = feature.Location.End
	}

	// if location.Start > end || location.End < start {
	// 	return INTERGENIC
	// }

	isPromoter := (feature.Location.Strand == "+" && mid >= start && mid <= feature.Location.Start+annotateDb.TSSRegion.Offset3P()) || (feature.Location.Strand == "-" && mid >= feature.Location.End-annotateDb.TSSRegion.Offset3P() && mid <= end)

	exons, err := annotateDb.GeneDb.InExon(location, feature.TranscriptId)

	if err != nil {
		return "", err
	}

	isExon := len(exons.Features) > 0

	isIntronic := mid >= feature.Location.Start && mid <= feature.Location.End

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
