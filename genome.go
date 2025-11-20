package genome

import (
	"database/sql"
	"fmt"
	"slices"
	"sort"
	"strings"

	"os"
	"path/filepath"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys/log"
)

type (
	GenomicFeatures struct {
		Location *dna.Location `json:"location"`
		//Id       string        `json:"id,omitempty"`
		//Name     string        `json:"name,omitempty"`
		Feature  Feature           `json:"feature"`
		Features []*GenomicFeature `json:"features"`
	}

	GeneDBInfo struct {
		Name    string `json:"name"`
		Genome  string `json:"genome"`
		Version string `json:"version"`
	}

	Feature string
	Level   uint8

	GeneDBCache struct {
		cacheMap map[string]GeneDB
		dir      string
	}

	GenomicFeature struct {
		Location     *dna.Location `json:"loc"`
		Type         string        `json:"type,omitempty"`
		Feature      Feature       `json:"feature"`
		GeneId       string        `json:"geneId,omitempty"`
		GeneSymbol   string        `json:"geneSymbol,omitempty"`
		TranscriptId string        `json:"transcriptId,omitempty"`
		//TranscriptName string            `json:"transcriptName,omitempty"`
		ExonId      string            `json:"exonId,omitempty"`
		Children    []*GenomicFeature `json:"children,omitempty"`
		ExonNumber  uint              `json:"exonNumber,omitempty"`
		Id          uint              `json:"-"`
		TssDist     int               `json:"tssDist"`
		IsCanonical bool              `json:"isCanonical"`
		IsLongest   bool              `json:"isLongest"`
		InPromoter  bool              `json:"inPromoter"`
		InExon      bool              `json:"inExon"`
		PromLabel   string            `json:"promLabel,omitempty"`
	}

	GeneDB interface {
		GeneDBInfo() (*GeneDBInfo, error)

		Close() error

		OverlappingGenes(location *dna.Location,
			featureLevel Feature,
			canonicalMode bool,
			geneTypeFilter string) ([]*GenomicFeature, error)

		SearchForGeneByName(search string,
			featureLevel Feature,
			fuzzy bool,
			canonical bool,
			geneType string,
			n uint16) ([]*GenomicFeature, error)

		WithinGenes(location *dna.Location, featureLevel Feature) (*GenomicFeatures, error)

		WithinGenesAndPromoter(location *dna.Location, featureLevel Feature, pad5p uint, pad3p uint) (*GenomicFeatures, error)

		InExon(location *dna.Location, transcriptId string) (*GenomicFeatures, error)

		ClosestGenes(location *dna.Location, closestN uint16) ([]*GenomicFeature, error)
	}
)

const (
	GeneDBInfoSql = `SELECT id, genome, version FROM info`

	MaxGeneInfoResults uint16 = 100

	// const IN_PROMOTER_SQL = `SELECT id, chr, start, end, strand, gene_id, gene_name, transcript_id, start - ?
	// 	FROM gene
	//  	WHERE level = 2 AND gene_id = ? AND chr = ? AND ? >= stranded_start - ? AND ? <= stranded_start + ?
	//  	ORDER BY start ASC`

	GeneFeature       Feature = "gene"
	TranscriptFeature Feature = "transcript"
	ExonFeature       Feature = "exon"

	GeneLevel       Level = 1
	TranscriptLevel Level = 2
	ExonLevel       Level = 3
)

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func (feature *GenomicFeature) TSS() *dna.Location {
	var s uint

	if feature.Location.Strand == "+" {
		s = feature.Location.Start
	} else {
		s = feature.Location.End
	}

	return dna.NewStrandedLocation(feature.Location.Chr, s, s, feature.Location.Strand)
}

func LevelToFeature(level Level) Feature {

	switch level {
	case TranscriptLevel:
		return TranscriptFeature
	case ExonLevel:
		return ExonFeature
	default:
		return GeneFeature
	}
}

func FeatureToLevel(feature string) Level {
	switch feature {
	case "t", "tran", "transcript", "2":
		return TranscriptLevel
	case "e", "ex", "exon", "3":
		return ExonLevel
	default:
		return GeneLevel
	}
}

func NewGeneDBCache(dir string, creator func(assembly string, dir string) GeneDB) *GeneDBCache {
	cacheMap := make(map[string]GeneDB)

	files, err := os.ReadDir(dir)

	log.Debug().Msgf("---- genedb ----")

	if err != nil {
		log.Fatal().Msgf("error opening %s", dir)
	}

	log.Debug().Msgf("caching gene databases in %s...", dir)

	cacheMap["hg19"] = creator("hg19", dir)
	cacheMap["grch38"] = creator("grch38", dir)
	cacheMap["mm10"] = creator("mm10", dir)

	for _, file := range files {
		basename := file.Name()

		if strings.Contains(basename, "gtf") && strings.HasSuffix(basename, ".db") {

			name := strings.TrimSuffix(basename, filepath.Ext(basename))
			name = strings.TrimPrefix(name, "gtf_")

			cacheMap[name] = creator(name, dir)

			log.Debug().Msgf("found gene database %s %s", name, filepath.Join(dir, basename))
		}
	}

	log.Debug().Msgf("---- end ----")

	return &GeneDBCache{dir: dir, cacheMap: cacheMap}
}

func (cache *GeneDBCache) GeneDB(assembly string, creator func(assembly string, dir string) GeneDB) (GeneDB, error) {
	_, ok := cache.cacheMap[assembly]

	if !ok {
		//db := NewGeneDB(assembly, filepath.Join(cache.dir, fmt.Sprintf("%s.db", assembly)))
		db := creator(assembly, filepath.Join(cache.dir, fmt.Sprintf("gtf_%s.db", assembly)))

		cache.cacheMap[assembly] = db
	}

	return cache.cacheMap[assembly], nil
}

func (cache *GeneDBCache) Dir() string {
	return cache.dir
}

func (cache *GeneDBCache) List() ([]*GeneDBInfo, error) {

	ids := make([]string, 0, len(cache.cacheMap))

	for id := range cache.cacheMap {

		ids = append(ids, id)
	}

	sort.Strings(ids)

	infos := make([]*GeneDBInfo, 0, len(ids))

	for _, id := range ids {
		info, err := cache.cacheMap[id].GeneDBInfo()

		if err != nil {
			return nil, err
		}

		infos = append(infos, info)
	}

	return infos, nil
}

func (cache *GeneDBCache) Close() {
	for _, db := range cache.cacheMap {
		db.Close()
	}
}

func RowsToFeatures(location *dna.Location, feature Feature, rows *sql.Rows) (*GenomicFeatures, error) {

	//log.Debug().Msgf("rowsToFeatures %s   %d %d", location.Chr, location.Start, location.End)

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	features, err := RowsToRecords(rows)

	if err != nil {
		return nil, err
	}

	ret := GenomicFeatures{Location: location, Feature: feature, Features: features}

	return &ret, nil
}

// func (genedb *GeneDB) IdToName(id uint) string {

// 	name := "n/a"

// 	// ignore error as it is not critical
// 	// if the id is not found, we return "n/a"
// 	genedb.db.QueryRow(ID_TO_NAME_SQL, id).Scan(&name)

// 	return name
// }

func RowsToRecords(rows *sql.Rows) ([]*GenomicFeature, error) {
	defer rows.Close()

	var id uint
	var feature Feature
	var chr string
	var start uint
	var end uint
	var strand string
	var geneId string
	var geneName string
	var geneType string
	var transcriptId string
	var transcriptName string
	var isCanonical bool
	var exonNumber uint

	var d int
	var err error

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	//var features = make([]*GenomicFeature, 0, 10)

	features := make([]*GenomicFeature, 0, 10)

	for rows.Next() {
		geneId = ""
		geneName = ""
		transcriptId = ""
		transcriptName = ""
		geneType = ""

		d = 0 // tss distance

		err = rows.Scan(&id,
			&feature,
			&chr,
			&start,
			&end,
			&strand,
			&geneId,
			&geneName,
			&geneType,
			&isCanonical,
			&transcriptId,
			&transcriptName,
			&exonNumber,
			&d)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		location := dna.NewLocation(chr, start, end)

		feature := GenomicFeature{
			Feature:  feature,
			Location: location,
			//Strand:       strand,
			GeneId:       geneId,
			GeneSymbol:   geneName,
			Type:         geneType,
			TranscriptId: transcriptId,
			IsCanonical:  isCanonical,
			TssDist:      d}

		features = append(features, &feature)
	}

	// enforce sorted correctly by chr and then position
	//sort.Sort(SortFeatureByPos(*features))

	return features, nil
}

// SortFeaturesByPos sorts features in place by chr, start, end
func SortFeaturesByPos(features []*GenomicFeature) {
	slices.SortFunc(features, func(a, b *GenomicFeature) int {
		ci := dna.ChromToInt(a.Location.Chr)
		cj := dna.ChromToInt(b.Location.Chr)

		// on different chrs so sort by chr
		if ci != cj {
			return int(ci) - int(cj)
		}

		// same chr so sort by position
		if a.Location.Start != b.Location.Start {
			return int(a.Location.Start) - int(b.Location.Start)
		}

		// same start so sort by end
		return int(a.Location.End) - int(b.Location.End)
	})
}

// type SortFeatureByPos []*GenomicFeature

// func (features SortFeatureByPos) Len() int      { return len(features) }
// func (features SortFeatureByPos) Swap(i, j int) { features[i], features[j] = features[j], features[i] }
// func (features SortFeatureByPos) Less(i, j int) bool {
// 	ci := dna.ChromToInt(features[i].Location.Chr)
// 	cj := dna.ChromToInt(features[j].Location.Chr)

// 	// on different chrs so sort by chr
// 	if ci != cj {
// 		return ci < cj
// 	}

// 	if features[i].Location.Start != features[j].Location.Start {
// 		return features[i].Location.Start < features[j].Location.Start
// 	}

// 	return features[i].Location.End < features[j].Location.End
// }
