package genedbcache

import (
	"sync"

	gene "github.com/antonybholmes/go-genes"
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

func Db(assembly string) (*gene.GeneDB, error) {
	return instance.Db(assembly)
}
