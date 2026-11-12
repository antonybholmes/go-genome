package genomedb

import (
	"database/sql"
	"errors"
	"path/filepath"
	"strings"
	"sync"

	"github.com/antonybholmes/go-genome"
	"github.com/antonybholmes/go-genome/gtfsqlite"
	"github.com/antonybholmes/go-genome/gtftabix"
	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-sys/db"
	"github.com/antonybholmes/go-sys/log"
	"github.com/antonybholmes/go-web"
)

type (
	GenomeDB struct {
		db   *sql.DB
		dir  string
		path string
	}
)

var (
	instance *GenomeDB
	once     sync.Once
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

	// AnnotationsSql = `SELECT DISTINCT
	// 	a.id,
	// 	a.public_id,
	// 	g.name AS genome,
	// 	asm.name AS assembly,
	// 	at.name AS type,
	// 	a.name,
	// 	a.url,
	// 	FROM annotations a
	// 	JOIN assemblies asm ON a.assembly_id = asm.id
	// 	JOIN genomes g ON asm.genome_id = g.id
	// 	JOIN annotation_types at ON a.annotation_type_id = at.id
	// 	JOIN assembly_aliases aa ON asm.id = aa.assembly_id
	// 	WHERE
	// 		(a.public_id = :id OR LOWER(aa.name) = LOWER(:id))
	// 	ORDER BY a.name`

	GtfsSql = `SELECT
		a.id,
		a.public_id,
		g.name AS genome,
		asm.name AS assembly,
		at.name AS type,
		a.name as name,
		a.url
		FROM annotations a
		JOIN assemblies asm ON a.assembly_id = asm.id
		JOIN genomes g ON asm.genome_id = g.id
		JOIN annotation_types at ON a.annotation_type_id = at.id
		WHERE LOWER(at.name) = 'gtf'
		ORDER BY asm.id`

	AnnotationsFromIdSql = `SELECT DISTINCT
		a.id,
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
		JOIN assembly_aliases aa ON asm.id = aa.assembly_id
		WHERE
			a.public_id = :id OR LOWER(aa.name) = LOWER(:id)`

	AnnotationsByTypeSql = `SELECT DISTINCT
		a.id,
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
		JOIN assembly_aliases aa ON asm.id = aa.assembly_id
		WHERE
			LOWER(at.name) = :type
			AND (asm.public_id = :assembly OR LOWER(aa.name) = LOWER(:assembly))
		ORDER BY a.id DESC`
)

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func NewGenomeDB(dbPath string) *GenomeDB {

	// db, err := sql.Open("sqlite3", path)

	// if err != nil {
	// 	log.Fatal().Msgf("%s", err)
	// }

	// defer db.Close()

	log.Debug().Msgf("initialising genome db with %s", dbPath)

	dir := filepath.Dir(dbPath)

	return &GenomeDB{dir: dir, path: dbPath, db: sys.Must(sql.Open(db.Sqlite3DB, dbPath+db.SqliteReadOnlySuffix))}

}

func (gdb *GenomeDB) Genomes() ([]*db.Entity, error) {

	genomes := make([]*db.Entity, 0, 10)

	rows, err := gdb.db.Query(GenomesSQL)

	if err != nil {
		return nil, err
	}

	defer rows.Close()

	for rows.Next() {
		var genome db.Entity

		err := rows.Scan(&genome.Id, &genome.PublicId, &genome.Name)

		if err != nil {
			return nil, err
		}

		genomes = append(genomes, &genome)
	}

	if err := rows.Err(); err != nil {
		return nil, err
	}

	return genomes, nil
}

func (gdb *GenomeDB) Assemblies(genome string) ([]*db.Entity, error) {

	assemblies := make([]*db.Entity, 0, 10)

	rows, err := gdb.db.Query(AssembliesSql, sql.Named("genome", web.FormatParam(genome)))

	if err != nil {
		return nil, err
	}

	defer rows.Close()

	for rows.Next() {
		var assembly db.Entity

		err := rows.Scan(&assembly.Id, &assembly.PublicId, &assembly.Name)

		if err != nil {
			return nil, err
		}

		assemblies = append(assemblies, &assembly)
	}

	if err := rows.Err(); err != nil {
		return nil, err
	}

	return assemblies, nil
}

func (gdb *GenomeDB) Annotation(id string) (*genome.Annotation, error) {

	namedArgs := []any{
		sql.Named("id", strings.ToLower(id)),
	}

	var annotation genome.Annotation

	err := gdb.db.QueryRow(AnnotationsFromIdSql, namedArgs...).Scan(
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

	return &annotation, nil
}

// For a given assembly and type, e.g. gtf, return all the associated
// annotation dbs
func (gdb *GenomeDB) Annotations(assembly string, annotationType string) ([]*genome.Annotation, error) {

	namedArgs := []any{
		sql.Named("assembly", strings.ToLower(assembly)),
		sql.Named("type", strings.ToLower(annotationType)),
	}

	datasetRows, err := gdb.db.Query(AnnotationsByTypeSql, namedArgs...)

	if err != nil {
		return nil, err
	}

	defer datasetRows.Close()

	ret := make([]*genome.Annotation, 0, 10)

	for datasetRows.Next() {
		var annotation genome.Annotation

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

	if err := datasetRows.Err(); err != nil {
		return nil, err
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

// Lookup a db using its public id
func (gdb *GenomeDB) GtfFromId(id string) (genome.GtfDB, error) {
	//_, ok := cache.cacheMap[assembly]

	annotation, err := gdb.Annotation(id)

	if err != nil {
		log.Debug().Msgf("Error loading annotation for id: %s, err: %s", id, err)
		return nil, err
	}

	log.Debug().Msgf("Loading GTF for annotation: %s", annotation.PublicId)

	if strings.HasSuffix(annotation.Url, ".db") {
		return gtfsqlite.NewGtfSqliteDB(gdb.dir, annotation), nil
	} else if strings.Contains(annotation.Url, "gtf") || strings.HasSuffix(annotation.Url, ".gtf") {
		return gtftabix.NewGtfTabixDB(gdb.dir, annotation), nil
	} else {
		return nil, errors.New("unsupported annotation URL: " + annotation.Url)
	}

}

// From the genome central db, look for the latest GTF annotation for the given assembly
func (gdb *GenomeDB) GtfFromAssembly(assembly string) (genome.GtfDB, error) {
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
	db := gtfsqlite.NewGtfSqliteDB(gdb.dir, annotation) //   cache.dbCreator(assembly, cache.dir) //filepath.Join(cache.dir, fmt.Sprintf("gtf_%s.db", assembly)))

	//cache.cacheMap[assembly] = db
	//}

	return db, nil //cache.cacheMap[assembly], nil
}

func (gdb *GenomeDB) Dir() string {
	return gdb.dir
}

func (gdb *GenomeDB) Gtfs() ([]*genome.Annotation, error) {

	rows, err := gdb.db.Query(GtfsSql)

	if err != nil {

		return nil, err
	}

	defer rows.Close()

	ret := make([]*genome.Annotation, 0, 10)

	for rows.Next() {
		var annotation genome.Annotation

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

	if err := rows.Err(); err != nil {
		return nil, err
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

func InitCache(path string) *GenomeDB {
	once.Do(func() {
		//instance = genome.NewGeneDBCache(dir, v2.NewGeneDB)
		instance = NewGenomeDB(path)
	})

	return instance
}

func Instance() *GenomeDB {
	return instance
}
