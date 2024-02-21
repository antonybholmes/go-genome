package genedbcache

import (
	gene "github.com/antonybholmes/go-genes"
)

var cache = gene.NewGeneDbCache()

func Init(dir string) {
	cache.Init(dir)
}

func Dir() string {
	return cache.Dir()
}

func Db(assembly string) (*gene.GeneDb, error) {
	return cache.Db(assembly)
}
