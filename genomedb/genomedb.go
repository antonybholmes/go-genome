package genomedb

import (
	"sync"

	"github.com/antonybholmes/go-genome"
)

var (
	instance *genome.GenomeDB
	once     sync.Once
)

func InitCache(path string) *genome.GenomeDB {
	once.Do(func() {
		//instance = genome.NewGeneDBCache(dir, v2.NewGeneDB)
		instance = genome.NewGenomeDB(path)
	})

	return instance
}

func Instance() *genome.GenomeDB {
	return instance
}
