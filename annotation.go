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
	TSSRegion *dna.TSSRegion
	ClosestN  uint8
}

func NewAnnotateDb(genesdb *GeneDB, tssRegion *dna.TSSRegion, closestN uint8) *AnnotateDb {
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

	log.Debug().Msgf("Annotating location %s %s", location, annotateDb.GeneDb.File)

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

		isPromoter := (gene.Strand == "+" && mid >= gene.Location.Start-annotateDb.TSSRegion.Offset5P() &&
			mid <= gene.Location.Start+annotateDb.TSSRegion.Offset3P()) ||
			(gene.Strand == "-" && mid >= gene.Location.End-annotateDb.TSSRegion.Offset3P() &&
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
				PromLabel: gene.PromLabel,
				Strand:    gene.Strand,
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

	// promoterMap := make(map[string]*GeneProm)

	// sort the ids by distance
	// distMap := make(map[uint][]string)

	// for id := range idMap {
	// 	d := promoterMap[id].AbsD
	// 	distMap[d] = append(distMap[d], id)
	// }

	// // now put the ids in distance order
	// distances := make([]uint, len(distMap))

	// for d := range distMap {
	// 	distances = append(distances, d)
	// }

	// slices.Sort(distances)

	// ids := make([]string, 0, len(distances))

	// for _, d := range distances {
	// 	dids := distMap[d]

	// 	// sort the gene ids
	// 	slices.Sort(dids)

	// 	ids = append(ids, dids...)
	// }

	// nids := len(transcriptsWithin.Features)
	// // arrays are always at least 1 element since if nothing is found
	// // we put NA
	// n := basemath.Max(1, nids)
	// geneIds := make([]string, n)
	// geneNames := make([]string, n)
	// geneStrands := make([]string, n)
	// tssDists := make([]string, n)
	// featureLocations := make([]string, n)
	// promLabels := make([]string, n)

	// if nids > 0 {
	// 	for i, gene := range transcriptsWithin.Features {
	// 		//p := promoterMap[id]
	// 		geneIds[i] = gene.GeneId
	// 		geneNames[i] = gene.GeneName // idMap[id] //GeneStrandLabel(idMap[id], p.Feature.Strand))
	// 		geneStrands[i] = gene.Strand
	// 		tssDists[i] = strconv.Itoa(gene.TssDist)
	// 		promLabels[i] = gene.PromLabel
	// 		featureLocations[i] = gene.Location.String()
	// 	}
	// } else {
	// 	//ids = append(ids, NA)
	// 	geneNames[0] = NA
	// 	geneStrands[0] = NA
	// 	tssDists[0] = NA
	// 	featureLocations[0] = NA
	// 	promLabels[0] = NA
	// }

	// for _, id := range ids {
	// 	p := promoterMap[id]

	// }

	// tssDists := make([]string, 0, len(ids))

	// for _, id := range ids {
	// 	p := promoterMap[id]
	// 	tssDists = append(tssDists, strconv.Itoa(p.D))
	// }

	// for _, id := range ids {
	// 	p := promoterMap[id]
	// 	featureLocations = append(featureLocations, p.Feature.ToLocation().String())
	// }

	closestGenes, err := annotateDb.GeneDb.ClosestGenes(location, 1000, FEATURE_TRANSCRIPT)

	if err != nil {
		log.Debug().Msgf("Error closest genes for location %s: %v", location, err)
		return nil, err
	}

	//closestGeneList := make([]*ClosestGene, 0, len(closestGenes.Features)) //[]*ClosestGene{}

	closestGeneIndexMap := make(map[string]uint8)

	closestMap := make(map[uint8]*AnnotationGene)

	for _, cg := range closestGenes.Features {

		// check cg.geneId in closestGeneMap
		closestIndex, ok := closestGeneIndexMap[cg.GeneId]

		if !ok {
			closestIndex = uint8(len(closestGeneIndexMap) + 1)

			if closestIndex > annotateDb.ClosestN {
				closestIndex = 0
			} else {
				closestGeneIndexMap[cg.GeneId] = closestIndex
			}
		}

		if closestIndex == 0 {
			continue
		}

		exons, err := annotateDb.GeneDb.InExon(location, cg.TranscriptId)

		if err != nil {
			return nil, err
		}

		isExon := len(exons.Features) > 0

		isPromoter := (cg.Strand == "+" && mid >= cg.Location.Start-annotateDb.TSSRegion.Offset5P() &&
			mid <= cg.Location.Start+annotateDb.TSSRegion.Offset3P()) ||
			(cg.Strand == "-" && mid >= cg.Location.End-annotateDb.TSSRegion.Offset3P() &&
				mid <= cg.Location.End+annotateDb.TSSRegion.Offset5P())

		isIntronic := mid >= cg.Location.Start && mid <= cg.Location.End

		//annotateDb.ClassifyFeature(location, cg)
		prom, ok := closestMap[closestIndex]

		if ok {
			// we test all transcripts to see if we can classify the gene
			// hence we OR the tests since we want to make as many of
			// them true as possible
			prom.IsIntronic = prom.IsIntronic || isIntronic
			prom.IsPromoter = prom.IsPromoter || isPromoter
			prom.IsExon = prom.IsExon || len(exons.Features) > 0
			prom.PromLabel = MakePromLabel(isPromoter, isExon, isIntronic)
		} else {
			closestMap[closestIndex] = &AnnotationGene{
				GeneId:     cg.GeneId,
				GeneName:   cg.GeneName,
				PromLabel:  cg.PromLabel,
				IsPromoter: isPromoter,
				IsIntronic: isIntronic,
				IsExon:     isExon,
				TssDist:    cg.TssDist,
				Strand:     cg.Strand,
			}
		}

		// closestGeneList = append(closestGeneList, &ClosestGene{Feature: &cg,
		// 	TssDist:   cg.TssDist,
		// 	PromLabel: label})
	}

	orderedClosestGenes := make([]*AnnotationGene, 0, len(closestMap))

	for i := 1; i <= int(annotateDb.ClosestN); i++ {
		g, ok := closestMap[uint8(i)]
		if ok {
			orderedClosestGenes = append(orderedClosestGenes, g)
		}
	}

	// slices.SortFunc(orderedClosestGenes, func(a, b *AnnotationGene) int {
	// 	return int(a.TssDist) - int(b.TssDist)
	// })

	annotation := GeneAnnotation{
		Location: location,
		// GeneIds:     strings.Join(geneIds, OUTPUT_FEATURE_SEP),
		// GeneSymbols: strings.Join(geneNames, OUTPUT_FEATURE_SEP),
		// GeneStrands: strings.Join(geneStrands, OUTPUT_FEATURE_SEP),
		// PromLabels:  strings.Join(promLabels, OUTPUT_FEATURE_SEP),

		// TSSDists:     strings.Join(tssDists, OUTPUT_FEATURE_SEP),
		// Locations:    strings.Join(featureLocations, OUTPUT_FEATURE_SEP),
		WithinGenes:  withinGenes,
		ClosestGenes: orderedClosestGenes,
	}

	return &annotation, nil
}

func (annotateDb *AnnotateDb) ClassifyFeature(location *dna.Location, feature *GenomicFeature) {
	mid := location.Mid()
	var start uint

	if feature.Strand == "+" {
		start = feature.Location.Start - annotateDb.TSSRegion.Offset5P()
	} else {
		start = feature.Location.Start
	}

	var end uint

	if feature.Strand == "-" {
		end = feature.Location.End + annotateDb.TSSRegion.Offset5P()
	} else {
		end = feature.Location.End
	}

	// if location.Start > end || location.End < start {
	// 	return INTERGENIC
	// }

	isPromoter := (feature.Strand == "+" && mid >= start && mid <= feature.Location.Start+annotateDb.TSSRegion.Offset3P()) || (feature.Strand == "-" && mid >= feature.Location.End-annotateDb.TSSRegion.Offset3P() && mid <= end)

	exons, err := annotateDb.GeneDb.InExon(location, feature.TranscriptId)

	if err != nil {
		return
	}

	isExon := len(exons.Features) > 0

	isIntronic := mid >= feature.Location.Start && mid <= feature.Location.End

	feature.PromLabel = MakePromLabel(isPromoter, isExon, isIntronic)
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
