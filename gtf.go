package genome

import (
	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys"
)

type (
	FeatureType string

	GenomicFeatures struct {
		Type     FeatureType       `json:"type" msgpack:"type"`
		Features []*GenomicFeature `json:"features" msgpack:"features"`
	}

	GenomicFeature struct {
		Type         FeatureType         `json:"type,omitempty" msgpack:"type,omitempty"`
		Biotype      string              `json:"biotype,omitempty" msgpack:"biotype,omitempty"`
		Symbol       string              `json:"symbol,omitempty" msgpack:"symbol,omitempty"`
		TranscriptId string              `json:"transcriptId,omitempty" msgpack:"transcriptId,omitempty"`
		ExonId       string              `json:"exonId,omitempty" msgpack:"exonId,omitempty"`
		Strand       string              `json:"strand,omitempty" msgpack:"strand,omitempty"`
		Label        string              `json:"label,omitempty" msgpack:"label,omitempty"`
		PublicId     string              `json:"id,omitempty" msgpack:"id,omitempty"`
		GeneId       string              `json:"geneId,omitempty" msgpack:"geneId,omitempty"`
		Children     []*GenomicFeature   `json:"children,omitempty" msgpack:"children,omitempty"`
		dna.Loc      `msgpack:",inline"` // need to explicity tell it to inline props
		TssDist      int                 `json:"tssDist,omitempty" msgpack:"tssDist,omitempty"`
		Id           int                 `json:"-" msgpack:"-"`
		ExonNumber   int                 `json:"exonNumber,omitempty" msgpack:"exonNumber,omitempty"`
		Value        float64             `json:"value,omitempty" msgpack:"value,omitempty"`
		InPromoter   bool                `json:"inPromoter,omitempty" msgpack:"inPromoter,omitempty"`
		InExon       bool                `json:"inExon,omitempty" msgpack:"inExon,omitempty"`
		IsIntragenic bool                `json:"isIntragenic,omitempty" msgpack:"isIntragenic,omitempty"`
		IsCanonical  bool                `json:"isCanonical,omitempty" msgpack:"isCanonical,omitempty"`
	}

	GenomicSearchResults struct {
		Location *dna.Location     `json:"location" msgpack:"location"`
		Type     FeatureType       `json:"type" msgpack:"type"`
		Features []*GenomicFeature `json:"features" msgpack:"features"`
	}

	GtfDB interface {
		InExon(location *dna.Location,
			transcriptId string,
			prom *dna.PromoterRegion) (bool, error)

		IntragenicFeatures(location *dna.Location,
			levels sys.Set[FeatureType],
			prom *dna.PromoterRegion,
			useOfficialGenes bool) ([]*GenomicFeature, error)

		ClosestGenes(location *dna.Location,
			prom *dna.PromoterRegion,
			closestN int8,
			useOfficialGenes bool) ([]*GenomicFeature, error)

		OverlappingGenes(location *dna.Location,
			levels sys.Set[FeatureType],
			prom *dna.PromoterRegion,
			canonicalMode bool,
			annotationMode bool,
			biotypeFilter string) ([]*GenomicFeature, error)

		// SearchByName searches for genomic features by name
		SearchByName(search string,
			levels sys.Set[FeatureType],
			canonicalMode bool,
			n int16) ([]*GenomicFeature, error)

		// WithinGenes(location *dna.Location,
		// 	feature string,
		// 	prom *dna.PromoterRegion) (*GenomicSearchResults, error)

		Close() error
	}
)

const (
	// AllLevels               string = "gene,transcript,exon,cds,utr"
	// GeneLevel               string = "gene"
	// TranscriptLevel         string = "transcript"
	// ExonLevel               string = "exon"
	// GeneAndTranscriptLevels string = "gene,transcript"
	// GeneAndExonLevels       string = "gene,exon"
	// TranscriptAndExonLevels string = "transcript,exon"

	FeatureTypeGene       FeatureType = "gene"
	FeatureTypeTranscript FeatureType = "transcript"
	FeatureTypeExon       FeatureType = "exon"
	FeatureTypeCds        FeatureType = "cds"
	FeatureTypeUtr        FeatureType = "utr"
)

var (
	ExonLevelMap = map[string]string{
		"exon": "exons",
		"cds":  "cds",
		"utr":  "utrs",
	}

	FeatureLevelsAll               = sys.NewSet[FeatureType](FeatureTypeGene, FeatureTypeTranscript, FeatureTypeExon, FeatureTypeCds, FeatureTypeUtr)
	FeatureLevelsGene              = sys.NewSet[FeatureType](FeatureTypeGene)
	FeatureLevelsTranscript        = sys.NewSet[FeatureType](FeatureTypeTranscript)
	FeatureLevelsExon              = sys.NewSet[FeatureType](FeatureTypeExon)
	FeatureLevelsGeneAndTranscript = sys.NewSet[FeatureType](FeatureTypeGene, FeatureTypeTranscript)
	FeatureLevelsGeneAndExon       = sys.NewSet[FeatureType](FeatureTypeGene, FeatureTypeExon)
	FeatureLevelsTranscriptAndExon = sys.NewSet[FeatureType](FeatureTypeTranscript, FeatureTypeExon)
)

func FeaturesToSearchResults(location *dna.Location,
	levels sys.Set[FeatureType],
	features []*GenomicFeature,
	annotationMode bool) (*GenomicSearchResults, error) {

	var level = MaxLevel(levels)

	ret := GenomicSearchResults{Location: location, Type: level, Features: features}

	return &ret, nil
}

func MaxLevel(levels sys.Set[FeatureType]) FeatureType {
	var level FeatureType

	switch {
	case levels.Has(FeatureTypeGene):
		level = FeatureTypeGene
	case levels.Has(FeatureTypeTranscript):
		level = FeatureTypeTranscript
	default:
		level = FeatureTypeExon
	}

	return level
}

func AddLocationToFeature(feature *GenomicFeature, location *dna.Location) {

	if feature.Type == FeatureTypeGene {
		feature.Chr = location.Chr()
		feature.Strand = location.Strand()
	}

	feature.Start = location.Start()
	feature.End = location.End()
}
