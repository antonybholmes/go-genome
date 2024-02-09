package gene

import (
	"fmt"
	"slices"
	"strconv"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-loctogene"
	"github.com/antonybholmes/go-utils"
)

const NA string = "n/a"
const PROMOTER string = "promoter"
const EXONIC string = "exonic"
const INTRONIC string = "intronic"
const INTERGENIC string = "intergenic"

const LABEL_SEP string = ","
const FEATURE_SEP string = ";"

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

type ClosestGene struct {
	Feature   *loctogene.GenomicFeature `json:"feature"`
	PromLabel string                    `json:"prom_label"`
	Dist      int                       `json:"dist"`
}

type GeneAnnotation struct {
	Location     *dna.Location  `json:"location"`
	GeneIds      string         `json:"gene_ids"`
	GeneSymbols  string         `json:"gene_symbols"`
	PromLabels   string         `json:"prom_labels"`
	Dists        string         `json:"dists"`
	Locations    string         `json:"tss"`
	ClosestGenes []*ClosestGene `json:"closest_genes"`
}

type GeneProm struct {
	Feature    *loctogene.GenomicFeature
	IsPromoter bool
	IsIntronic bool
	IsExon     bool
	AbsD       uint
	D          int
}

type AnnotateDb struct {
	GeneDb    *loctogene.LoctogeneDb
	TSSRegion *dna.TSSRegion
	N         uint16
}

func NewAnnotateDb(genesdb *loctogene.LoctogeneDb, tssRegion *dna.TSSRegion, n uint16) *AnnotateDb {
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
		loctogene.Transcript,
		utils.UintMax(annotateDb.TSSRegion.Offset5P(), annotateDb.TSSRegion.Offset3P()),
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

		// println!(
		//     "{} {} {} {} {}",
		//     gene.gene_id, gene.GeneSymbol, gene.Start, gene.End, gene.strand
		// );

		idMap[id] = gene.GeneSymbol

		//let labels = annotate.classify_location(location, gene);

		exons, err := annotateDb.GeneDb.InExon(location, id)

		if err != nil {
			return nil, err
		}

		isExon := len(exons.Features) > 0

		isPromoter := (gene.Strand == "+" && mid >= gene.Start-annotateDb.TSSRegion.Offset5P() && mid <= gene.Start+annotateDb.TSSRegion.Offset3P()) || (gene.Strand == "-" && mid >= gene.End-annotateDb.TSSRegion.Offset3P() && mid <= gene.End+annotateDb.TSSRegion.Offset5P())

		isIntronic := mid >= gene.Start && mid <= gene.End

		var d int = 0

		if gene.Strand == "+" {
			d = int(gene.Start) - int(mid)
		} else {
			d = int(gene.End) - int(mid)
		}

		//fmt.Printf("%d %d %d %s\n", d, gene.End, mid, gene.GeneSymbol)

		// update by inserting default case and then updating

		absD := uint(utils.AbsInt(d))

		prom, ok := promoterMap[id]

		if ok {
			prom.IsIntronic = prom.IsIntronic || isIntronic
			prom.IsPromoter = prom.IsPromoter || isPromoter
			prom.IsExon = prom.IsExon || len(exons.Features) > 0

			// record the feature that is closest to our site
			if absD < prom.AbsD {
				prom.Feature = &gene
				prom.D = d
				prom.AbsD = absD
			}
		} else {
			promoterMap[id] = &GeneProm{
				Feature:    &gene,
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
	var distances []uint

	for d := range distMap {
		distances = append(distances, d)
	}

	slices.Sort(distances)

	ids := []string{}

	for _, d := range distances {
		dids := distMap[d]

		slices.Sort(dids)

		ids = append(ids, dids...)
	}

	geneSymbols := []string{}

	for _, id := range ids {
		p := promoterMap[id]
		geneSymbols = append(geneSymbols, fmt.Sprintf("%s(%s)", idMap[id], p.Feature.Strand))
	}

	promLabels := []string{}

	for _, id := range ids {
		p := promoterMap[id]
		promLabels = append(promLabels, makeLabel(p.IsPromoter, p.IsExon, p.IsIntronic))
	}

	tssDists := []string{}

	for _, id := range ids {
		p := promoterMap[id]
		tssDists = append(tssDists, strconv.Itoa(p.D))
	}

	featureLocations := []string{}

	for _, id := range ids {
		p := promoterMap[id]
		featureLocations = append(featureLocations, p.Feature.ToLocation().String())
	}

	if len(ids) == 0 {
		ids = append(ids, NA)
		geneSymbols = append(geneSymbols, NA)
		tssDists = append(tssDists, NA)
		featureLocations = append(featureLocations, NA)
	}

	closestGenes, err := annotateDb.GeneDb.ClosestGenes(location, annotateDb.N, loctogene.Gene)

	if err != nil {
		return nil, err
	}

	closestGeneList := []*ClosestGene{}

	for _, cg := range closestGenes.Features {
		label := annotateDb.ClassifyLocation(location, &cg)

		if err != nil {
			return nil, err
		}

		closestGeneList = append(closestGeneList, &ClosestGene{Feature: &cg,
			Dist:      cg.Dist,
			PromLabel: label})
	}

	annotation := GeneAnnotation{
		Location:     location,
		GeneIds:      strings.Join(ids, FEATURE_SEP),
		GeneSymbols:  strings.Join(geneSymbols, FEATURE_SEP),
		PromLabels:   strings.Join(promLabels, FEATURE_SEP),
		Dists:        strings.Join(tssDists, FEATURE_SEP),
		Locations:    strings.Join(featureLocations, FEATURE_SEP),
		ClosestGenes: closestGeneList,
	}

	return &annotation, nil
}

func (annotateDb *AnnotateDb) ClassifyLocation(location *dna.Location, feature *loctogene.GenomicFeature) string {
	mid := location.Mid()
	var s uint

	if feature.Strand == "+" {
		s = feature.Start - annotateDb.TSSRegion.Offset5P()
	} else {
		s = feature.Start
	}

	var e uint

	if feature.Strand == "-" {
		e = feature.End + annotateDb.TSSRegion.Offset5P()
	} else {
		e = feature.End
	}

	if location.Start > e || location.End < s {
		return INTERGENIC
	}

	isPromoter := (feature.Strand == "+" && mid >= s && mid <= feature.Start+annotateDb.TSSRegion.Offset3P()) || (feature.Strand == "-" && mid >= feature.End-annotateDb.TSSRegion.Offset3P() && mid <= e)

	exons, err := annotateDb.GeneDb.InExon(location, feature.GeneId)

	if err != nil {
		return ""
	}

	isExon := len(exons.Features) > 0

	isIntronic := mid >= feature.Start && mid <= feature.End

	return makeLabel(isPromoter, isExon, isIntronic)
}

func makeLabel(isPromoter bool, isExon bool, isIntronic bool) string {
	labels := []string{}

	if isPromoter {
		labels = append(labels, PROMOTER)
	}

	if isExon {
		labels = append(labels, EXONIC)
	} else {
		if isIntronic {
			labels = append(labels, INTRONIC)
		}
	}

	return strings.Join(labels, LABEL_SEP)
}
