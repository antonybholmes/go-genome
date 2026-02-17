package genomedb

import (
	"sync"

	"github.com/antonybholmes/go-genome"
)

var (
	instance *genome.GeneDBCache
	once     sync.Once
)

func InitCache(dir string) *genome.GeneDBCache {
	once.Do(func() {
		//instance = genome.NewGeneDBCache(dir, v2.NewGeneDB)
		instance = genome.NewGeneDBCache(dir)
	})

	return instance
}

func GetInstance() *genome.GeneDBCache {
	return instance
}

func Dir() string {
	return instance.Dir()
}

func GeneDB(assembly string) (*genome.GeneDB, error) {
	//return instance.GeneDB(assembly, v2.NewGeneDB)
	return instance.GeneDB(assembly)
}

func List() ([]*genome.GeneDBInfo, error) {
	return instance.List()
}
