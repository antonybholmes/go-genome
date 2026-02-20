package genomedb

import (
	"sync"

	"github.com/antonybholmes/go-genome"
)

var (
	instance *genome.GenomeDB
	once     sync.Once
)

func InitCache(dir string) *genome.GenomeDB {
	once.Do(func() {
		//instance = genome.NewGeneDBCache(dir, v2.NewGeneDB)
		instance = genome.NewGenomeDB(dir)
	})

	return instance
}

func GetInstance() *genome.GenomeDB {
	return instance
}

func Dir() string {
	return instance.Dir()
}

func GtfDB(assembly string) (*genome.GtfDB, error) {
	//return instance.GeneDB(assembly, v2.NewGeneDB)
	return instance.GtfDB(assembly)
}

func Gtfs() ([]*genome.Annotation, error) {
	return instance.Gtfs()
}
