package genome

import (
	"strings"

	"github.com/antonybholmes/go-dna"
)

type (
	GenomicFeatures struct {
		Type     string            `json:"type" msgpack:"type"`
		Features []*GenomicFeature `json:"features" msgpack:"features"`
	}

	GenomicFeature struct {
		Type         string              `json:"type,omitempty" msgpack:"type,omitempty"`
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
		Type     string            `json:"type" msgpack:"type"`
		Features []*GenomicFeature `json:"features" msgpack:"features"`
	}

	GtfDB interface {
		InExon(location *dna.Location,
			transcriptId string,
			prom *dna.PromoterRegion) (bool, error)

		IntragenicFeatures(location *dna.Location,
			levels string,
			prom *dna.PromoterRegion,
			useOfficialGenes bool) ([]*GenomicFeature, error)

		ClosestGenes(location *dna.Location,
			prom *dna.PromoterRegion,
			closestN int8,
			useOfficialGenes bool) ([]*GenomicFeature, error)

		OverlappingGenes(location *dna.Location,
			levels string,
			prom *dna.PromoterRegion,
			canonicalMode bool,
			annotationMode bool,
			biotypeFilter string) ([]*GenomicFeature, error)

		SearchByName(search string,
			level string,
			canonicalMode bool,
			n int16) ([]*GenomicFeature, error)

		// WithinGenes(location *dna.Location,
		// 	feature string,
		// 	prom *dna.PromoterRegion) (*GenomicSearchResults, error)

		Close() error
	}
)

const (
	AllLevels               string = "gene,transcript,exon,cds,utr"
	GeneLevel               string = "gene"
	TranscriptLevel         string = "transcript"
	ExonLevel               string = "exon"
	GeneAndTranscriptLevels string = "gene,transcript"
	GeneAndExonLevels       string = "gene,exon"
	TranscriptAndExonLevels string = "transcript,exon"

	FeatureTypeGene       string = "gene"
	FeatureTypeTranscript string = "transcript"
	FeatureTypeExon       string = "exon"
	FeatureTypeCds        string = "cds"
	FeatureTypeUtr        string = "utr"
)

var (
	ExonLevelMap = map[string]string{
		"exon": "exons",
		"cds":  "cds",
		"utr":  "utrs",
	}
)

func FeaturesToSearchResults(location *dna.Location,
	levels string,
	features []*GenomicFeature,
	annotationMode bool) (*GenomicSearchResults, error) {

	var level = MaxLevel(levels)

	ret := GenomicSearchResults{Location: location, Type: level, Features: features}

	return &ret, nil
}

func MaxLevel(levels string) string {
	var level string

	switch strings.ToLower(string(levels)) {
	case "gene", "genes":
		level = GeneLevel
	case "transcript", "transcripts":
		level = TranscriptLevel
	default:
		level = ExonLevel
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
