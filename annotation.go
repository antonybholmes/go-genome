package genome

import (
	"fmt"
	"slices"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/rs/zerolog/log"
)

const NA string = "n/a"
const PROMOTER string = "promoter"
const EXONIC string = "exonic"
const INTRONIC string = "intronic"
const INTERGENIC string = "intergenic"
const INTRAGENIC string = "intragenic"

const GROUP_SEP string = ","
const OUTPUT_FEATURE_SEP string = "|"

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

// type ClosestGene struct {
// 	Feature   *GenomicFeature `json:"feature"`
// 	PromLabel string          `json:"promLabel"`
// 	TssDist   int             `json:"tssDist"`
// }

type AnnotationGene struct {
	GeneId     string `json:"geneId"`
	GeneName   string `json:"geneName"`
	PromLabel  string `json:"promLabel"`
	Strand     string `json:"strand"`
	IsPromoter bool   `json:"-"`
	IsIntronic bool   `json:"-"`
	IsExon     bool   `json:"-"`
	TssDist    int    `json:"tssDist"`
}

type GeneAnnotation struct {
	Location *dna.Location `json:"loc"`
	// GeneIds      string            `json:"geneIds"`
	// GeneSymbols  string            `json:"geneSymbols"`
	// GeneStrands  string            `json:"geneStrands"`
	// PromLabels   string            `json:"promLabels"`
	// TSSDists     string            `json:"tssDists"`
	// Locations    string            `json:"geneLocs"`
	WithinGenes  []*AnnotationGene `json:"withinGenes"`
	ClosestGenes []*AnnotationGene `json:"closestGenes"`
}

// type ByAbsD []GeneProm
// func (a ByAbsD) Len() int           { return len(a) }
// func (a ByAbsD) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
// func (a ByAbsD) Less(i, j int) bool { return a[i].AbsD < a[j].AbsD }

type AnnotateDb struct {
	GeneDb    *GeneDB
	TSSRegion *dna.PromoterRegion
	ClosestN  uint8
}

func NewAnnotateDb(genesdb *GeneDB, tssRegion *dna.PromoterRegion, closestN uint8) *AnnotateDb {
	return &AnnotateDb{
		GeneDb:    genesdb,
		TSSRegion: tssRegion,
		ClosestN:  closestN,
	}
}

func (annotateDb *AnnotateDb) Annotate(location *dna.Location) (*GeneAnnotation, error) {
	mid := location.Mid()

	// extend search area to account  for promoter
	// let search_loc: Location = Location::new(
	//     &location.chr,
	//     location.Start - annotate.tssRegion.offset_5p.abs(),
	//     location.End + annotate.tssRegion.offset_5p.abs(),
	// )?;

	//log.Debug().Msgf("Annotating location %s %s", location, annotateDb.GeneDb.File)

	genesWithin, err := annotateDb.GeneDb.WithinGenesAndPromoter(
		location,
		annotateDb.TSSRegion.Offset5P(),
		annotateDb.TSSRegion.Offset3P(),
	)

	if err != nil {
		return nil, err
	}

	// we need the unique ids to symbols
	//idMap := make(map[string]string)
	withinMap := make(map[string]AnnotationGene)

	//let mut distMap: HashMap<&str, bool> = HashMap::new();

	for _, gene := range genesWithin.Features {
		id := gene.GeneId

		//idMap[id] = gene.GeneName

		//let labels = annotate.classify_location(location, gene);

		exons, err := annotateDb.GeneDb.InExon(location, gene.TranscriptId)

		if err != nil {
			return nil, err
		}

		isExon := len(exons.Features) > 0

		isPromoter := (gene.Location.Strand == "+" && mid >= gene.Location.Start-annotateDb.TSSRegion.Offset5P() &&
			mid <= gene.Location.Start+annotateDb.TSSRegion.Offset3P()) ||
			(gene.Location.Strand == "-" && mid >= gene.Location.End-annotateDb.TSSRegion.Offset3P() &&
				mid <= gene.Location.End+annotateDb.TSSRegion.Offset5P())

		isIntronic := mid >= gene.Location.Start && mid <= gene.Location.End

		// var d int = 0

		// if gene.Strand == "+" {
		// 	d = int(gene.Location.Start) - int(mid)
		// } else {
		// 	d = int(gene.Location.End) - int(mid)
		// }

		//fmt.Printf("%d %d %d %s\n", d, gene.Location.End, mid, gene.GeneSymbol)

		// update by inserting default case and then updating

		//absD := uint(basemath.AbsInt(gene.TssDist))

		prom, ok := withinMap[id]

		if ok {
			// we test all transcripts to see if we can classify the gene
			// hence we OR the tests since we want to make as many of
			// them true as possible
			prom.IsIntronic = prom.IsIntronic || isIntronic
			prom.IsPromoter = prom.IsPromoter || isPromoter
			prom.IsExon = prom.IsExon || len(exons.Features) > 0
			prom.PromLabel = MakePromLabel(isPromoter, isExon, isIntronic)
		} else {
			withinMap[id] = AnnotationGene{
				GeneId:    gene.GeneId,
				GeneName:  gene.GeneName,
				PromLabel: MakePromLabel(isPromoter, isExon, isIntronic),
				Strand:    gene.Location.Strand,
				//IsPromoter: isPromoter,
				//IsIntronic: isIntronic,
				//IsExon:     isExon,
				TssDist: gene.TssDist,
			}
		}

	}

	withinGenes := make([]*AnnotationGene, 0, len(withinMap))

	for _, gene := range withinMap {
		withinGenes = append(withinGenes, &gene)
	}

	slices.SortFunc(withinGenes, func(a, b *AnnotationGene) int {
		return int(a.TssDist) - int(b.TssDist)
	})

	closestAnnotations := make([]*AnnotationGene, 0, annotateDb.ClosestN)

	if false {
		closestGenes, err := annotateDb.GeneDb.ClosestGenes(location, uint16(annotateDb.ClosestN))

		if err != nil {
			log.Debug().Msgf("Error closest genes for location %s: %v", location, err)
			return nil, err
		}

		//closestGeneList := make([]*ClosestGene, 0, len(closestGenes.Features)) //[]*ClosestGene{}

		//closestGeneIndexMap := make(map[string]uint8)

		//closestMap := make(map[uint8]*AnnotationGene)

		for _, cg := range closestGenes {

			exons, err := annotateDb.GeneDb.InExon(location, cg.TranscriptId)

			if err != nil {
				return nil, err
			}

			isExon := len(exons.Features) > 0

			isPromoter := (cg.Location.Strand == "+" && mid >= cg.Location.Start-annotateDb.TSSRegion.Offset5P() &&
				mid <= cg.Location.Start+annotateDb.TSSRegion.Offset3P()) ||
				(cg.Location.Strand == "-" && mid >= cg.Location.End-annotateDb.TSSRegion.Offset3P() &&
					mid <= cg.Location.End+annotateDb.TSSRegion.Offset5P())

			isIntronic := mid >= cg.Location.Start && mid <= cg.Location.End

			closestAnnotations = append(closestAnnotations, &AnnotationGene{
				GeneId:     cg.GeneId,
				GeneName:   cg.GeneName,
				PromLabel:  MakePromLabel(isPromoter, isExon, isIntronic),
				IsPromoter: isPromoter,
				IsIntronic: isIntronic,
				IsExon:     isExon,
				TssDist:    cg.TssDist,
				Strand:     cg.Location.Strand,
			})

		}
	}

	annotation := GeneAnnotation{
		Location: location,
		// GeneIds:     strings.Join(geneIds, OUTPUT_FEATURE_SEP),
		// GeneSymbols: strings.Join(geneNames, OUTPUT_FEATURE_SEP),
		// GeneStrands: strings.Join(geneStrands, OUTPUT_FEATURE_SEP),
		// PromLabels:  strings.Join(promLabels, OUTPUT_FEATURE_SEP),

		// TSSDists:     strings.Join(tssDists, OUTPUT_FEATURE_SEP),
		// Locations:    strings.Join(featureLocations, OUTPUT_FEATURE_SEP),
		WithinGenes:  withinGenes,
		ClosestGenes: closestAnnotations,
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

func MakePromLabel(isPromoter bool, isExon bool, isIntronic bool) string {
	labels := make([]string, 0, 3)

	if isPromoter {
		labels = append(labels, PROMOTER)
	}

	// favor exonic over intronic
	if isExon {
		labels = append(labels, EXONIC)
	} else {
		if isIntronic {
			labels = append(labels, INTRONIC)
		}
	}

	if isExon || isIntronic {
		labels = append(labels, INTRAGENIC)
	} else {
		labels = append(labels, INTERGENIC)
	}

	//slices.Sort(labels)

	return strings.Join(labels, GROUP_SEP)
}
