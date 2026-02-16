package routes

import (
	"bytes"
	"encoding/csv"
	"errors"
	"fmt"
	"net/http"
	"strconv"
	"strings"

	"github.com/antonybholmes/go-dna"
	dnaroutes "github.com/antonybholmes/go-dna/routes"
	"github.com/antonybholmes/go-genome"
	"github.com/antonybholmes/go-genome/genomedb"
	basemath "github.com/antonybholmes/go-math"
	"github.com/antonybholmes/go-sys/log"
	"github.com/antonybholmes/go-web"
	"github.com/gin-gonic/gin"
)

// A GeneQuery contains info from query params.
type (
	GeneQuery struct {
		Feature  string
		Db       genome.GeneDB
		Assembly string
		GeneType string // e.g. "protein_coding", "non_coding", etc.
		// only show canonical genes
		Canonical bool
		Promoter  *dna.PromoterRegion
	}

	GenesResp struct {
		Location *dna.Location            `json:"location"`
		Features []*genome.GenomicFeature `json:"features"`
	}

	AnnotationResponse struct {
		Status int                      `json:"status"`
		Data   []*genome.GeneAnnotation `json:"data"`
	}
)

const (
	DefaultClosestN int = 5
	MaxAnnotations  int = 100
)

var (
	ErrLocationCannotBeEmpty = errors.New("location cannot be empty")
	ErrSearchTooShort        = errors.New("search too short")

	genomeNormMap = map[string]string{
		"hg19":   "gencode.v48lift37.basic.grch37",
		"grch37": "gencode.v48lift37.basic.grch37",
		"hg38":   "gencode.v48.basic.grch38",
		"grch38": "gencode.v48.basic.grch38",
		"mm10":   "gencode.vM25.basic.mm10",
	}
)

func ParseFeature(c *gin.Context) string {

	feature := web.FormatParam(c.Query("feature")) //ParseFeature(c)

	if feature == "" {
		feature = genome.AllLevels
	}

	return feature

	// feature := strings.ToLower(c.Query("feature"))

	// switch {
	// case strings.Contains(feature, "gene,transcript,exon"):
	// 	return genome.AllLevels
	// case strings.Contains(feature, "gene,transcript"):
	// 	return genome.GeneAndTranscriptLevels
	// case strings.Contains(feature, "gene,exon"):
	// 	return genome.GeneAndExonLevels
	// case strings.Contains(feature, "transcript,exon"):
	// 	return genome.TranscriptAndExonLevels
	// case strings.Contains(feature, "gene"):
	// 	return genome.GeneLevel
	// case strings.Contains(feature, "transcript"):
	// 	return genome.TranscriptLevel
	// case strings.Contains(feature, "exon"):
	// 	return genome.ExonLevel
	// default:
	// 	return genome.AllLevels
	// }

}

func parseGeneQuery(c *gin.Context) (*GeneQuery, error) {

	assembly := web.FormatParam(c.Param("assembly"))

	if assembly == "" {
		return nil, errors.New("assembly cannot be empty")
	}

	// check if assembly is valid
	if _, ok := genomeNormMap[assembly]; !ok {
		return nil, fmt.Errorf("invalid assembly: %s", assembly)
	}

	// change assembly to normalized form
	assembly = genomeNormMap[assembly]

	log.Debug().Msgf("using assembly: %s", assembly)

	//dbFile := genomeToFileMap[assembly]

	feature := ParseFeature(c)

	// switch c.Query("feature") {
	// case "exon":
	// 	feature = genome.ExonLevel
	// case "transcript":
	// 	feature = genome.TranscriptLevel
	// default:
	// 	feature = genome.GeneLevel
	// }

	canonical := strings.HasPrefix(strings.ToLower(c.Query("canonical")), "t")

	geneType := ParseGeneType(c)

	promoterRegion := ParsePromoterRegion(c)

	db, err := genomedb.GeneDB(assembly)

	if err != nil {
		return nil, fmt.Errorf("unable to open database for assembly %s %s", assembly, err)
	}

	return &GeneQuery{
			Assembly:  assembly,
			GeneType:  geneType,
			Db:        db,
			Feature:   feature,
			Canonical: canonical,
			Promoter:  promoterRegion},
		nil
}

// func GeneDBInfoRoute(c *gin.Context) {
// 	query, err := ParseGeneQuery(c, c.Param("assembly"))

// 	if err != nil {
// 		return web.ErrorReq(err)
// 	}

// 	info, _ := query.Db.GeneDBInfo()

// 	// if err != nil {
// 	// 	return web.ErrorReq(err)
// 	// }

// 	web.MakeDataResp(c, "", &info)
// }

func GenomesRoute(c *gin.Context) {
	infos, err := genomedb.GetInstance().List()

	if err != nil {
		c.Error(err)
		return
	}

	web.MakeDataResp(c, "", infos)
}

func OverlappingGenesRoute(c *gin.Context) {
	locations, err := dnaroutes.ParseLocationsFromPost(c) // dnaroutes.ParseLocationsFromPost(c)

	if err != nil {
		c.Error(err)
		return
	}

	query, err := parseGeneQuery(c)

	if err != nil {
		c.Error(err)
		return
	}

	if len(locations) == 0 {
		web.BadReqResp(c, ErrLocationCannotBeEmpty)
	}

	ret := make([]*GenesResp, 0, len(locations))

	log.Debug().Msgf("querying for %d locations ", len(locations))

	for _, location := range locations {
		features, err := query.Db.OverlappingGenes(location,
			query.Feature,
			query.Promoter,
			query.Canonical,
			query.GeneType)

		if err != nil {
			c.Error(err)
			return
		}

		ret = append(ret, &GenesResp{Location: location, Features: features})

	}

	web.MakeDataResp(c, "", &ret)
}

func SearchForGeneByNameRoute(c *gin.Context) {
	search := c.Query("q") // dnaroutes.ParseLocationsFromPost(c)

	if search == "" {
		web.BadReqResp(c, ErrSearchTooShort)
		return
	}

	fuzzyMode := c.Query("mode") == "fuzzy"

	n := web.ParseN(c, 20)

	query, err := parseGeneQuery(c)

	if err != nil {
		c.Error(err)
		return
	}

	canonical := strings.HasPrefix(strings.ToLower(c.Query("canonical")), "t")

	log.Debug().Msgf("searching for gene: %s, fuzzy: %v, canonical: %v, type: %s",
		search, fuzzyMode, canonical, query.GeneType)

	features, _ := query.Db.SearchByName(search,
		query.Feature,
		fuzzyMode,
		canonical,
		c.Query("type"),
		int16(n))

	log.Debug().Msgf("found %d features", len(features))

	// if err != nil {
	// 	return web.ErrorReq(err)
	// }

	web.MakeDataResp(c, "", &features)
}

func WithinGenesRoute(c *gin.Context) {
	locations, err := dnaroutes.ParseLocationsFromPost(c) // dnaroutes.ParseLocationsFromPost(c)

	if err != nil {
		c.Error(err)
		return
	}

	query, err := parseGeneQuery(c)

	if err != nil {
		c.Error(err)
		return
	}

	data := make([]*genome.GenomicFeatures, len(locations))

	for li, location := range locations {
		genes, err := query.Db.WithinGenes(location, query.Feature, query.Promoter)

		if err != nil {
			c.Error(err)
			return
		}

		data[li] = genes
	}

	web.MakeDataResp(c, "", &data)
}

// Find the n closest genes to a location
func ClosestGeneRoute(c *gin.Context) {
	locations, err := dnaroutes.ParseLocationsFromPost(c)

	if err != nil {
		c.Error(err)
		return
	}

	query, err := parseGeneQuery(c)

	if err != nil {
		c.Error(err)
		return
	}

	closestN := web.ParseNumParam(c, "closest", DefaultClosestN)

	data := make([]*genome.GenomicFeatures, len(locations))

	for li, location := range locations {
		genes, err := query.Db.ClosestGenes(location, query.Promoter, int8(closestN))

		if err != nil {
			c.Error(err)
			return
		}

		data[li] = &genome.GenomicFeatures{Location: location, Feature: genome.GeneLevel, Features: genes}
	}

	web.MakeDataResp(c, "", &data)
}

func ParseGeneType(c *gin.Context) string {
	geneType := c.Query("type")

	// user can specify gene type in query string, but we sanitize it
	switch {
	case strings.Contains(geneType, "protein"):
		return "protein_coding"
	default:
		return ""
	}

}

func ParsePromoterRegion(c *gin.Context) *dna.PromoterRegion {

	v := c.Query("promoter")

	if v == "" {
		return dna.DefaultPromoterRegion()
	}

	tokens := strings.Split(v, ",")

	s, err := strconv.Atoi(tokens[0])

	if err != nil {
		return dna.DefaultPromoterRegion()
	}

	e, err := strconv.Atoi(tokens[1])

	if err != nil {
		return dna.DefaultPromoterRegion()
	}

	return dna.NewPromoterRegion(s, e)
}

func AnnotateRoute(c *gin.Context) {
	locations, err := dnaroutes.ParseLocationsFromPost(c)

	if err != nil {
		c.Error(err)
		return
	}

	// limit amount of data returned per request to 1000 entries at a time
	locations = locations[0:basemath.Min(len(locations), MaxAnnotations)]

	query, err := parseGeneQuery(c)

	if err != nil {
		c.Error(err)
		return
	}

	closestN := web.ParseNumParam(c, "closest", DefaultClosestN)

	tssRegion := ParsePromoterRegion(c)

	output := web.ParseOutput(c)

	annotationDb := genome.NewAnnotateDb(query.Db, tssRegion, int8(closestN))

	data := make([]*genome.GeneAnnotation, len(locations))

	for li, location := range locations {
		annotations, err := annotationDb.Annotate(location, query.Feature)

		if err != nil {
			log.Error().Msgf("Error annotating location %s: %v", location, err)
			c.Error(err)
			return
		}

		data[li] = annotations
	}

	if output == "text" {
		tsv, err := MakeGeneTable(data, tssRegion)

		if err != nil {
			c.Error(err)
			return
		}

		c.String(http.StatusOK, tsv)
	} else {

		c.JSON(http.StatusOK, AnnotationResponse{Status: http.StatusOK, Data: data})
	}
}

func MakeGeneTable(
	data []*genome.GeneAnnotation,
	ts *dna.PromoterRegion,
) (string, error) {
	var buffer bytes.Buffer
	wtr := csv.NewWriter(&buffer)
	wtr.Comma = '\t'

	closestN := len(data[0].ClosestGenes)

	headers := make([]string, 5+4*closestN)

	headers[0] = "Location"
	headers[1] = "Gene Id"
	headers[2] = "Gene Symbol"
	headers[3] = fmt.Sprintf(
		"Relative To Gene (prom=-%d/+%dkb)",
		ts.Upstream()/1000,
		ts.Downstream()/1000)
	headers[4] = "TSS Distance"
	//headers[5] = "Gene Location"

	idx := 6
	for i := 1; i <= closestN; i++ {
		headers[idx] = fmt.Sprintf("#%d Closest Id", i)
		headers[idx] = fmt.Sprintf("#%d Closest Gene Symbols", i)
		headers[idx] = fmt.Sprintf(
			"#%d Relative To Closet Gene (prom=-%d/+%dkb)",
			i,
			ts.Upstream()/1000,
			ts.Downstream()/1000)
		headers[idx] = fmt.Sprintf("#%d TSS Closest Distance", i)
		headers[idx] = fmt.Sprintf("#%d Gene Location", i)
	}

	err := wtr.Write(headers)

	if err != nil {
		return "", err
	}

	for _, annotation := range data {
		n := len(annotation.WithinGenes)
		geneIds := make([]string, n)
		geneNames := make([]string, n)
		promLabels := make([]string, n)
		tssDists := make([]string, n)

		for i, gene := range annotation.WithinGenes {
			geneIds[i] = gene.GeneId
			geneNames[i] = gene.GeneSymbol
			promLabels[i] = gene.Label
			tssDists[i] = strconv.Itoa(gene.TssDist)

		}

		row := []string{annotation.Location.String(),
			strings.Join(geneIds, genome.FeatureSeparator),
			strings.Join(geneNames, genome.FeatureSeparator),
			strings.Join(promLabels, genome.FeatureSeparator),
			strings.Join(tssDists, genome.FeatureSeparator)}

		for _, closestGene := range annotation.ClosestGenes {
			row = append(row, closestGene.GeneId)
			row = append(row, genome.GeneWithStrandLabel(closestGene.GeneSymbol, closestGene.Location.Strand()))
			row = append(row, closestGene.Label)
			row = append(row, strconv.Itoa(closestGene.TssDist))
			//row = append(row, closestGene.Location.String())
		}

		err := wtr.Write(row)

		if err != nil {
			return "", err
		}
	}

	wtr.Flush()

	return buffer.String(), nil
}
