package genes

import (
	"fmt"
	"slices"
	"strconv"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/rs/zerolog/log"

	"github.com/antonybholmes/go-basemath"
)

const NA string = "n/a"
const PROMOTER string = "promoter"
const EXONIC string = "exonic"
const INTRONIC string = "intronic"
const INTERGENIC string = "intergenic"

const GROUP_SEP string = ","
const OUTPUT_FEATURE_SEP string = "|"

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

// type ClosestGene struct {
// 	Feature   *GenomicFeature `json:"feature"`
// 	PromLabel string          `json:"promLabel"`
// 	TssDist   int             `json:"tssDist"`
// }

type GeneAnnotation struct {
	Location     *dna.Location     `json:"loc"`
	GeneIds      string            `json:"geneIds"`
	GeneSymbols  string            `json:"geneSymbols"`
	GeneStrands  string            `json:"geneStrands"`
	PromLabels   string            `json:"promLabels"`
	TSSDists     string            `json:"tssDists"`
	Locations    string            `json:"geneLocs"`
	ClosestGenes []*GenomicFeature `json:"closestGenes"`
}

type GeneProm struct {
	Feature    *GenomicFeature
	IsPromoter bool
	IsIntronic bool
	IsExon     bool
	AbsD       uint
	D          int
}

type AnnotateDb struct {
	GeneDb    *GeneDB
	TSSRegion *dna.TSSRegion
	N         uint16
}

func NewAnnotateDb(genesdb *GeneDB, tssRegion *dna.TSSRegion, n uint16) *AnnotateDb {
	return &AnnotateDb{
		GeneDb:    genesdb,
		TSSRegion: tssRegion,
		N:         n,
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

	genesWithin, err := annotateDb.GeneDb.WithinGenesAndPromoter(
		location,
		LEVEL_TRANSCRIPT,
		basemath.UintMax(annotateDb.TSSRegion.Offset5P(), annotateDb.TSSRegion.Offset3P()),
	)

	if err != nil {
		return nil, err
	}

	// we need the unique ids to symbols
	idMap := make(map[string]string)
	promoterMap := make(map[string]*GeneProm)

	//let mut distMap: HashMap<&str, bool> = HashMap::new();

	for _, gene := range genesWithin.Features {
		id := gene.GeneId

		idMap[id] = gene.GeneSymbol

		log.Debug().Msgf("%s %s %v", gene.GeneId, gene.GeneSymbol, gene)

		//let labels = annotate.classify_location(location, gene);

		exons, err := annotateDb.GeneDb.InExon(location, id)

		if err != nil {
			return nil, err
		}

		isExon := len(exons.Features) > 0

		isPromoter := (gene.Strand == "+" && mid >= gene.Location.Start-annotateDb.TSSRegion.Offset5P() &&
			mid <= gene.Location.Start+annotateDb.TSSRegion.Offset3P()) ||
			(gene.Strand == "-" && mid >= gene.Location.End-annotateDb.TSSRegion.Offset3P() &&
				mid <= gene.Location.End+annotateDb.TSSRegion.Offset5P())

		isIntronic := mid >= gene.Location.Start && mid <= gene.Location.End

		var d int = 0

		if gene.Strand == "+" {
			d = int(gene.Location.Start) - int(mid)
		} else {
			d = int(gene.Location.End) - int(mid)
		}

		//fmt.Printf("%d %d %d %s\n", d, gene.Location.End, mid, gene.GeneSymbol)

		// update by inserting default case and then updating

		absD := uint(basemath.AbsInt(d))

		prom, ok := promoterMap[id]

		if ok {
			// we test all transcripts to see if we can classify the gene
			// hence we OR the tests since we want to make as many of
			// them true as possible
			prom.IsIntronic = prom.IsIntronic || isIntronic
			prom.IsPromoter = prom.IsPromoter || isPromoter
			prom.IsExon = prom.IsExon || len(exons.Features) > 0

			// record the feature that is closest to our site
			if absD < prom.AbsD {
				prom.Feature = gene
				prom.D = d
				prom.AbsD = absD
			}
		} else {
			promoterMap[id] = &GeneProm{
				Feature:    gene,
				IsPromoter: isPromoter,
				IsIntronic: isIntronic,
				IsExon:     isExon,
				D:          d,
				AbsD:       absD,
			}
		}

	}

	// sort the ids by distance
	distMap := make(map[uint][]string)

	for id := range idMap {
		d := promoterMap[id].AbsD
		distMap[d] = append(distMap[d], id)
	}

	// now put the ids in distance order
	distances := make([]uint, 0, len(distMap))

	for d := range distMap {
		distances = append(distances, d)
	}

	slices.Sort(distances)

	ids := make([]string, 0, len(distances))

	for _, d := range distances {
		dids := distMap[d]

		// sort the gene ids
		slices.Sort(dids)

		ids = append(ids, dids...)
	}

	nids := len(ids)
	// arrays are always at least 1 element since if nothing is found
	// we put NA
	n := basemath.IntMax(1, nids)
	geneSymbols := make([]string, n)
	geneStrands := make([]string, n)
	tssDists := make([]string, n)
	featureLocations := make([]string, n)
	promLabels := make([]string, n)

	if nids > 0 {
		for iid, id := range ids {
			p := promoterMap[id]
			geneSymbols[iid] = idMap[id] //GeneStrandLabel(idMap[id], p.Feature.Strand))
			geneStrands[iid] = p.Feature.Strand
			tssDists[iid] = strconv.Itoa(p.D)
			promLabels[iid] = PromLabel(p.IsPromoter, p.IsExon, p.IsIntronic)
			featureLocations[iid] = p.Feature.Location.String()
		}
	} else {
		ids = append(ids, NA)
		geneSymbols[0] = NA
		geneStrands[0] = NA
		tssDists[0] = NA
		featureLocations[0] = NA
		promLabels[0] = NA
	}

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

	closestGenes, err := annotateDb.GeneDb.ClosestGenes(location, annotateDb.N, LEVEL_GENE)

	if err != nil {
		return nil, err
	}

	//closestGeneList := make([]*ClosestGene, 0, len(closestGenes.Features)) //[]*ClosestGene{}

	for _, cg := range closestGenes.Features {
		cg.PromLabel = annotateDb.ClassifyLocation(location, cg)

		// closestGeneList = append(closestGeneList, &ClosestGene{Feature: &cg,
		// 	TssDist:   cg.TssDist,
		// 	PromLabel: label})
	}

	annotation := GeneAnnotation{
		Location:    location,
		GeneIds:     strings.Join(ids, OUTPUT_FEATURE_SEP),
		GeneSymbols: strings.Join(geneSymbols, OUTPUT_FEATURE_SEP),
		GeneStrands: strings.Join(geneStrands, OUTPUT_FEATURE_SEP),
		PromLabels:  strings.Join(promLabels, OUTPUT_FEATURE_SEP),

		TSSDists:     strings.Join(tssDists, OUTPUT_FEATURE_SEP),
		Locations:    strings.Join(featureLocations, OUTPUT_FEATURE_SEP),
		ClosestGenes: closestGenes.Features,
	}

	return &annotation, nil
}

func (annotateDb *AnnotateDb) ClassifyLocation(location *dna.Location, feature *GenomicFeature) string {
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

	if location.Start > end || location.End < start {
		return INTERGENIC
	}

	isPromoter := (feature.Strand == "+" && mid >= start && mid <= feature.Location.Start+annotateDb.TSSRegion.Offset3P()) || (feature.Strand == "-" && mid >= feature.Location.End-annotateDb.TSSRegion.Offset3P() && mid <= end)

	exons, err := annotateDb.GeneDb.InExon(location, feature.GeneId)

	if err != nil {
		return ""
	}

	isExon := len(exons.Features) > 0

	isIntronic := mid >= feature.Location.Start && mid <= feature.Location.End

	return PromLabel(isPromoter, isExon, isIntronic)
}

func GeneWithStrandLabel(id string, strand string) string {
	return fmt.Sprintf("%s:%s", id, strand)
}

func PromLabel(isPromoter bool, isExon bool, isIntronic bool) string {
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

	return strings.Join(labels, GROUP_SEP)
}
