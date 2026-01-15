package v1

import (
	"database/sql"
	"fmt"
	"path/filepath"

	"github.com/antonybholmes/go-genome"
	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-sys/log"
)

type (
	V1GeneDB struct {
		db   *sql.DB
		file string
		//withinGeneAndPromStmt *sql.Stmt
		name string
	}
)

// func (feature *GenomicFeature) ToLocation() *dna.Location {
// 	return dna.NewLocation(feature.Chr, feature.Start, feature.End)
// }

func NewGeneDB(assembly string, dir string) genome.GeneDB {
	file := filepath.Join(dir, fmt.Sprintf("%s.db", assembly))

	log.Debug().Msgf("opening gene database %s", file)

	db := sys.Must(sql.Open(sys.Sqlite3DB, file))

	return &V1GeneDB{name: assembly, file: file, db: db} //withinGeneStmt: sys.Must(db.Prepare(WITHIN_GENE_SQL)),
	//withinGeneAndPromStmt: sys.Must(db.Prepare(WITHIN_GENE_AND_PROMOTER_SQL))

	//inExonStmt:      sys.Must(db.Prepare(IN_EXON_SQL)),
	//closestGeneStmt: sys.Must(db.Prepare(CLOSEST_GENE_SQL))

}

func (genedb *V1GeneDB) Close() error {
	return genedb.db.Close()
}
