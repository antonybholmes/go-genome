package genome

import (
	"database/sql"
	"os"
	"slices"
	"sort"
	"strings"

	"path/filepath"

	"github.com/antonybholmes/go-dna"

	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-sys/log"
)

type (

	//string string
	//string string

	GeneDBCache struct {
		cacheMap map[string]*GeneDBInfo
		//dbCreator func(assembly string, dir string) GeneDB
		dir string
	}

	// GeneDB interface {
	// 	//GeneDBInfo() (*GeneDBInfo, error)

	// 	Close() error

	// 	OverlappingGenes(location *dna.Location,
	// 		level string,
	// 		prom *dna.PromoterRegion,
	// 		canonicalMode bool,
	// 		geneTypeFilter string) ([]*GenomicFeature, error)

	// 	SearchByName(search string,
	// 		level string,
	// 		fuzzy bool,
	// 		canonical bool,
	// 		geneType string,
	// 		n int16) ([]*GenomicFeature, error)

	// 	WithinGenes(location *dna.Location, level string, prom *dna.PromoterRegion) (*GenomicFeatures, error)

	// 	WithinGenesAndPromoter(location *dna.Location, level string, prom *dna.PromoterRegion) ([]*GenomicFeature, error)

	// 	InExon(location *dna.Location, transcriptId string, prom *dna.PromoterRegion) ([]*GenomicFeature, error)

	// 	ClosestGenes(location *dna.Location,
	// 		prom *dna.PromoterRegion,
	// 		closestN int8) ([]*GenomicFeature, error)
	// }
)

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func (feature *GenomicFeature) TSS() (*dna.Location, error) {
	var s int

	if feature.Location.Strand() == "+" {
		s = feature.Location.Start()
	} else {
		s = feature.Location.End()
	}

	return dna.NewStrandedLocation(feature.Location.Chr(), s, s, feature.Location.Strand())
}

// func LevelToFeature(level Level) string {

// 	switch level {
// 	case TranscriptLevel:
// 		return TranscriptFeature
// 	case ExonLevel:
// 		return ExonFeature
// 	default:
// 		return GeneFeature
// 	}
// }

// func FeatureToLevel(feature string) Level {
// 	switch feature {
// 	case "t", "tran", "transcript", "2":
// 		return TranscriptLevel
// 	case "e", "ex", "exon", "3":
// 		return ExonLevel
// 	default:
// 		return GeneLevel
// 	}
// }

func LoadGeneDBInfo(file string) (*GeneDBInfo, error) {

	db, err := sql.Open(sys.Sqlite3DB, file+sys.SqliteReadOnlySuffix)

	if err != nil {
		return nil, err //fmt.Errorf("could not open gene database file %s", file)
	}

	defer db.Close()

	var info GeneDBInfo

	err = db.QueryRow(GeneDBInfoSql).Scan(&info.Id,
		&info.Genome,
		&info.Assembly,
		&info.Version,
		&info.File)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	return &info, nil
}

func NewGeneDBCache(dir string) *GeneDBCache {
	cacheMap := make(map[string]*GeneDBInfo)

	files, err := os.ReadDir(dir)

	log.Debug().Msgf("---- genedb ----")

	if err != nil {
		log.Fatal().Msgf("error opening %s", dir)
	}

	log.Debug().Msgf("caching gene databases in %s...", dir)

	//cacheMap["hg19"] = creator("hg19", dir)
	//cacheMap["grch38"] = creator("grch38", dir)
	//cacheMap["mm10"] = creator("mm10", dir)

	for _, file := range files {
		basename := file.Name()

		if strings.Contains(basename, "gencode") && strings.HasSuffix(basename, ".db") {

			name := strings.TrimSuffix(basename, filepath.Ext(basename))

			path := filepath.Join(dir, basename)

			info, err := LoadGeneDBInfo(path)

			if err != nil {
				log.Error().Msgf("could not load gene database info from %s: %v", path, err)
				continue
			}

			cacheMap[info.Assembly] = info

			if info.Assembly == "grch37" {
				// clone and change assembly
				hg19Info := GeneDBInfo{
					Id:       info.Id,
					Genome:   info.Genome,
					Assembly: "hg19",
					Version:  info.Version,
					File:     info.File,
				}
				cacheMap["hg19"] = &hg19Info
			}

			if info.Assembly == "grch38" {
				// clone and change assembly
				hg38Info := GeneDBInfo{
					Id:       info.Id,
					Genome:   info.Genome,
					Assembly: "hg38",
					Version:  info.Version,
					File:     info.File,
				}
				cacheMap["hg38"] = &hg38Info
			}

			if info.Assembly == "grcm38" {
				// clone and change assembly
				mm10Info := GeneDBInfo{
					Id:       info.Id,
					Genome:   info.Genome,
					Assembly: "mm10",
					Version:  info.Version,
					File:     info.File,
				}
				cacheMap["mm10"] = &mm10Info
			}

			log.Debug().Msgf("found gene database %s %s", name, filepath.Join(dir, basename))
		}
	}

	log.Debug().Msgf("---- end ----")

	return &GeneDBCache{dir: dir,
		//dbCreator: creator,
		cacheMap: cacheMap,
	}
}

func (cache *GeneDBCache) GeneDB(assembly string) (*GeneDB, error) {
	//_, ok := cache.cacheMap[assembly]

	//if !ok {
	//db := NewGeneDB(assembly, filepath.Join(cache.dir, fmt.Sprintf("%s.db", assembly)))
	db := NewGeneDB(assembly, cache.dir) //   cache.dbCreator(assembly, cache.dir) //filepath.Join(cache.dir, fmt.Sprintf("gtf_%s.db", assembly)))

	//cache.cacheMap[assembly] = db
	//}

	return db, nil //cache.cacheMap[assembly], nil
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
		info := cache.cacheMap[id]

		infos = append(infos, info)
	}

	return infos, nil
}

// func (cache *GeneDBCache) Close() {
// 	for _, db := range cache.cacheMap {
// 		db.Close()
// 	}
// }

func RowsToFeatures(location *dna.Location, level string, rows *sql.Rows) (*GenomicFeatures, error) {

	//log.Debug().Msgf("rowsToFeatures %s   %d %d", location.Chr, location.Start, location.End)

	// 10 seems a reasonable guess for the number of features we might see, just
	// to reduce slice reallocation
	features, err := RowsToRecords(rows)

	if err != nil {
		return nil, err
	}

	ret := GenomicFeatures{Location: location, Feature: level, Features: features}

	return &ret, nil
}

// func (genedb *GeneDB) IdToName(id int) string {

// 	name := "n/a"

// 	// ignore error as it is not critical
// 	// if the id is not found, we return "n/a"
// 	db.QueryRow(ID_TO_NAME_SQL, id).Scan(&name)

// 	return name
// }

func RowsToRecords(rows *sql.Rows) ([]*GenomicFeature, error) {
	defer rows.Close()

	var id int
	var level string
	var chr string
	var start int
	var end int
	var strand string
	var geneId string
	var geneName string
	var geneType string
	var transcriptId string
	var transcriptName string
	var isCanonical bool
	var exonNumber int

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
			&level,
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

		location, err := dna.NewLocation(chr, start, end)

		if err != nil {
			return nil, err
		}

		feature := GenomicFeature{
			Feature:  level,
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
		ci := dna.ChromToInt(a.Location.Chr())
		cj := dna.ChromToInt(b.Location.Chr())

		// on different chrs so sort by chr
		if ci != cj {
			return int(ci) - int(cj)
		}

		// same chr so sort by position
		if a.Location.Start() != b.Location.Start() {
			return int(a.Location.Start()) - int(b.Location.Start())
		}

		// same start so sort by end
		return int(a.Location.End()) - int(b.Location.End())
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

func MaxLevel(levels string) string {
	var level string

	switch strings.ToLower(string(levels)) {
	case "gene", "genes":
		level = GeneLevel
	case "transcript", "transcripts":
		level = TranscriptLevel
	default:
		level = GeneLevel
	}

	return level
}
