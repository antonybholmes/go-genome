package genome

import (
	"github.com/antonybholmes/go-dna"
)

type (

	// GtfDBInfo struct {
	// 	PublicId string `json:"id"`
	// 	Genome   string `json:"genome"`
	// 	Assembly string `json:"assembly"`
	// 	Name     string `json:"name"`
	// 	File     string `json:"-"`
	// 	Id       int    `json:"-"`
	// }

	GenomicFeatures struct {
		Type     string            `json:"type"`
		Features []*GenomicFeature `json:"features"`
	}

	GenomicFeature struct {
		Location   *dna.Location `json:"loc"`
		PublicId   string        `json:"id,omitempty"`
		Label      string        `json:"label,omitempty"`
		Biotype    string        `json:"biotype,omitempty"`
		GeneId     string        `json:"geneId,omitempty"`
		Symbol     string        `json:"symbol,omitempty"`
		Transcript string        `json:"transcript,omitempty"`
		Exon       string        `json:"exon,omitempty"`
		//Strand       string            `json:"strand,omitempty"`
		Type         string            `json:"type,omitempty"`
		Children     []*GenomicFeature `json:"children,omitempty"`
		TssDist      int               `json:"tssDist,omitempty"`
		Value        float64           `json:"value,omitempty"`
		Id           int               `json:"-"`
		ExonNumber   int               `json:"exonNumber,omitempty"`
		InPromoter   bool              `json:"inPromoter,omitempty"`
		InExon       bool              `json:"inExon,omitempty"`
		IsIntragenic bool              `json:"isIntragenic,omitempty"`
		IsLongest    bool              `json:"isLongest,omitempty"`
		IsCanonical  bool              `json:"isCanonical,omitempty"`
	}

	GenomicSearchResults struct {
		Location *dna.Location     `json:"location"`
		Type     string            `json:"type"`
		Features []*GenomicFeature `json:"features"`
	}

	GtfDB interface {
		InExon(location *dna.Location,
			transcriptId string,
			prom *dna.PromoterRegion) ([]*GenomicFeature, error)

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

		WithinGenes(location *dna.Location,
			feature string,
			prom *dna.PromoterRegion) (*GenomicSearchResults, error)

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
)

var (
	ExonLevelMap = map[string]string{
		"exon": "exons",
		"cds":  "cds",
		"utr":  "utrs",
	}
)
