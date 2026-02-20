package genome

import (
	"database/sql"
	"errors"
	"strings"

	"path/filepath"

	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-web"

	"github.com/antonybholmes/go-sys/log"
)

type (

	//string string
	//string string

	GenomeDB struct {
		//cacheMap map[string]*GtfDBInfo
		//dbCreator func(assembly string, dir string) GeneDB
		dir string
		db  *sql.DB
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

	Annotation struct {
		Id       int    `json:"-"`
		PublicId string `json:"id"`
		Genome   string `json:"genome"`
		Assembly string `json:"assembly"`
		Type     string `json:"type"`
		Name     string `json:"name"`
		Url      string `json:"-"`
	}
)

const (
	GenomesSQL = `SELECT DISTINCT
		g.id,
		g.public_id,
		g.name
		FROM genomes g
		ORDER BY g.name`

	AssembliesSql = `SELECT
		a.id,
		a.public_id,
		a.name
		FROM assemblies a
		JOIN genomes g ON asm.genome_id = g.id
		WHERE g.public_id = :genome OR LOWER(g.name) = :genome
		ORDER BY a.name`

	AnnotationsSql = `SELECT DISTINCT
		a.public_id,
		g.name AS genome,
		am.name AS assembly,
		at.name AS type,
		a.name,
		a.url,
		FROM annotations a
		JOIN assemblies am ON a.assembly_id = am.id
		JOIN genomes g ON asm.genome_id = g.id
		JOIN annotation_types at ON a.annotation_type_id = at.id
		JOIN assembly_aliases aa ON am.id = aa.assembly_id
		WHERE 
			(am.public_id = :assembly OR LOWER(aa.alias) = :assembly)
		ORDER BY a.name`

	GtfsSql = `SELECT DISTINCT
		a.id,
		a.public_id,
		a.genome,
		a.assembly,
		a.type,
		a.name,
		a.url
		FROM (
			SELECT
			a.id,
			a.public_id,
			g.name AS genome,
			asm.name AS assembly,
			at.name AS type,
			a.name,
			a.url,
			ROW_NUMBER() OVER (PARTITION BY asm.name ORDER BY a.id DESC) AS rank
			FROM annotations a
			JOIN assemblies asm ON a.assembly_id = asm.id
			JOIN genomes g ON asm.genome_id = g.id
			JOIN annotation_types at ON a.annotation_type_id = at.id
			JOIN assembly_aliases aa ON asm.id = aa.assembly_id
			WHERE
				at.name = 'GTF'
		) a
		WHERE a.rank = 1
		ORDER BY a.genome, a.name`

	AnnotationsByTypeSql = `SELECT DISTINCT
		a.public_id,
		g.name AS genome,
		asm.name AS assembly,
		at.name AS type,
		a.name,
		a.url
		FROM annotations a
		JOIN assemblies asm ON a.assembly_id = asm.id
		JOIN genomes g ON asm.genome_id = g.id
		JOIN annotation_types at ON a.annotation_type_id = at.id
		JOIN assembly_aliases aa ON am.id = aa.assembly_id
		WHERE
			LOWER(at.name) = :type
			AND (asm.public_id = :assembly OR LOWER(aa.alias) = :assembly)
		ORDER BY a.id DESC`
)

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func NewGenomeDB(dir string) *GenomeDB {

	// db, err := sql.Open("sqlite3", path)

	// if err != nil {
	// 	log.Fatal().Msgf("%s", err)
	// }

	// defer db.Close()

	return &GenomeDB{dir: dir, db: sys.Must(sql.Open(sys.Sqlite3DB, filepath.Join(dir, "genome.db"+sys.SqliteReadOnlySuffix)))}
}

func (gdb *GenomeDB) Genomes() ([]*sys.Entity, error) {

	genomes := make([]*sys.Entity, 0, 10)

	rows, err := gdb.db.Query(GenomesSQL)

	if err != nil {
		return nil, err
	}

	defer rows.Close()

	for rows.Next() {
		var genome sys.Entity

		err := rows.Scan(&genome.Id, &genome.PublicId, &genome.Name)

		if err != nil {
			return nil, err
		}

		genomes = append(genomes, &genome)
	}

	return genomes, nil
}

func (gdb *GenomeDB) Assemblies(genome string) ([]*sys.Entity, error) {

	assemblies := make([]*sys.Entity, 0, 10)

	rows, err := gdb.db.Query(AssembliesSql, sql.Named("genome", web.FormatParam(genome)))

	if err != nil {
		return nil, err
	}

	defer rows.Close()

	for rows.Next() {
		var assembly sys.Entity

		err := rows.Scan(&assembly.Id, &assembly.PublicId, &assembly.Name)

		if err != nil {
			return nil, err
		}

		assemblies = append(assemblies, &assembly)
	}

	return assemblies, nil
}

func (gdb *GenomeDB) Annotations(assembly string, annotationType string) ([]*Annotation, error) {

	namedArgs := []any{
		sql.Named("assembly", strings.ToLower(assembly)),
		sql.Named("type", strings.ToLower(annotationType)),
	}

	datasetRows, err := gdb.db.Query(AnnotationsByTypeSql, namedArgs...)

	if err != nil {
		log.Debug().Msgf("%s", err)
		return nil, err
	}

	defer datasetRows.Close()

	ret := make([]*Annotation, 0, 10)

	for datasetRows.Next() {
		var annotation Annotation

		err := datasetRows.Scan(
			&annotation.Id,
			&annotation.PublicId,
			&annotation.Genome,
			&annotation.Assembly,
			&annotation.Type,
			&annotation.Name,
			&annotation.Url)

		if err != nil {
			return nil, err
		}

		ret = append(ret, &annotation)
	}

	return ret, nil
}

// func (feature *GenomicFeature) TSS() (*dna.Location, error) {
// 	var s int

// 	if feature.Location.Strand() == "+" {
// 		s = feature.Location.Start()
// 	} else {
// 		s = feature.Location.End()
// 	}

// 	return dna.NewStrandedLocation(feature.Location.Chr(), s, s, feature.Location.Strand())
// }

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

// From the genome central db, look for the latest GTF annotation for the given assembly
func (gdb *GenomeDB) GtfDB(assembly string) (*GtfDB, error) {
	//_, ok := cache.cacheMap[assembly]

	annotations, err := gdb.Annotations(assembly, "gtf")

	if err != nil {
		return nil, err
	}

	if len(annotations) == 0 {
		return nil, errors.New("no GTF annotations found for assembly " + assembly)
	}

	// use the first annotation
	annotation := annotations[0]

	//if !ok {
	//db := NewGeneDB(assembly, filepath.Join(cache.dir, fmt.Sprintf("%s.db", assembly)))
	db := NewGtfDB(gdb.dir, annotation) //   cache.dbCreator(assembly, cache.dir) //filepath.Join(cache.dir, fmt.Sprintf("gtf_%s.db", assembly)))

	//cache.cacheMap[assembly] = db
	//}

	return db, nil //cache.cacheMap[assembly], nil
}

func (gdb *GenomeDB) Dir() string {
	return gdb.dir
}

func (gdb *GenomeDB) Gtfs() ([]*Annotation, error) {

	rows, err := gdb.db.Query(GtfsSql)

	if err != nil {
		log.Debug().Msgf("%s", err)
		return nil, err
	}

	defer rows.Close()

	ret := make([]*Annotation, 0, 10)

	for rows.Next() {
		var annotation Annotation

		err := rows.Scan(
			&annotation.Id,
			&annotation.PublicId,
			&annotation.Genome,
			&annotation.Assembly,
			&annotation.Type,
			&annotation.Name,
			&annotation.Url)

		if err != nil {
			return nil, err
		}

		ret = append(ret, &annotation)
	}

	return ret, nil
}

// func NewGeneDBCache(dir string) *GenomeDB {
// 	cacheMap := make(map[string]*GtfDBInfo)

// 	files, err := os.ReadDir(dir)

// 	log.Debug().Msgf("---- genedb ----")

// 	if err != nil {
// 		log.Fatal().Msgf("error opening %s", dir)
// 	}

// 	log.Debug().Msgf("caching gene databases in %s...", dir)

// 	//cacheMap["hg19"] = creator("hg19", dir)
// 	//cacheMap["grch38"] = creator("grch38", dir)
// 	//cacheMap["mm10"] = creator("mm10", dir)

// 	for _, file := range files {
// 		basename := file.Name()

// 		if strings.Contains(basename, "gencode") && strings.HasSuffix(basename, ".db") {

// 			name := strings.TrimSuffix(basename, filepath.Ext(basename))

// 			path := filepath.Join(dir, basename)

// 			info, err := LoadGeneDBInfo(path)

// 			if err != nil {
// 				log.Error().Msgf("could not load gene database info from %s: %v", path, err)
// 				continue
// 			}

// 			cacheMap[info.Assembly] = info

// 			if info.Assembly == "grch37" {
// 				// clone and change assembly
// 				hg19Info := GtfDBInfo{
// 					Id:       info.Id,
// 					Genome:   info.Genome,
// 					Assembly: "hg19",
// 					Name:     info.Name,
// 					File:     info.File,
// 				}
// 				cacheMap["hg19"] = &hg19Info
// 			}

// 			if info.Assembly == "grch38" {
// 				// clone and change assembly
// 				hg38Info := GtfDBInfo{
// 					Id:       info.Id,
// 					Genome:   info.Genome,
// 					Assembly: "hg38",
// 					Name:     info.Name,
// 					File:     info.File,
// 				}
// 				cacheMap["hg38"] = &hg38Info
// 			}

// 			if info.Assembly == "grcm38" {
// 				// clone and change assembly
// 				mm10Info := GtfDBInfo{
// 					Id:       info.Id,
// 					Genome:   info.Genome,
// 					Assembly: "mm10",
// 					Name:     info.Name,
// 					File:     info.File,
// 				}
// 				cacheMap["mm10"] = &mm10Info
// 			}

// 			log.Debug().Msgf("found gene database %s %s", name, filepath.Join(dir, basename))
// 		}
// 	}

// 	log.Debug().Msgf("---- end ----")

// 	return &GenomeDB{dir: dir,
// 		//dbCreator: creator,
// 		cacheMap: cacheMap,
// 	}
// }
