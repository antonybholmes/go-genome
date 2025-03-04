package genomedbcache

import (
	"sync"

	gene "github.com/antonybholmes/go-genome"
)

var instance *gene.GeneDBCache
var once sync.Once

func InitCache(dir string) *gene.GeneDBCache {
	once.Do(func() {
		instance = gene.NewGeneDBCache(dir)
	})

	return instance
}

func GetInstance() *gene.GeneDBCache {
	return instance
}

func Dir() string {
	return instance.Dir()
}

func GeneDB(assembly string) (*gene.GeneDB, error) {
	return instance.GeneDB(assembly)
}
