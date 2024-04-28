package genedbcache

import (
	"sync"

	gene "github.com/antonybholmes/go-genes"
)

var instance *gene.GeneDbCache
var once sync.Once

func InitCache(dir string) *gene.GeneDbCache {
	once.Do(func() {
		instance = gene.NewGeneDbCache(dir)
	})

	return instance
}

func GetInstance() *gene.GeneDbCache {
	return instance
}

func Dir() string {
	return instance.Dir()
}

func Db(assembly string) (*gene.GeneDB, error) {
	return instance.Db(assembly)
}
