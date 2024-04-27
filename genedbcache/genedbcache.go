package genedbcache

import (
	gene "github.com/antonybholmes/go-genes"
)

var cache *gene.GeneDbCache

func InitCache(dir string) {
	cache = gene.NewGeneDbCache(dir)
}

func Dir() string {
	return cache.Dir()
}

func Db(assembly string) (*gene.GeneDB, error) {
	return cache.Db(assembly)
}
