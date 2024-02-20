package genedbcache

import (
	gene "github.com/antonybholmes/go-genes"
)

var Cache = gene.NewGeneDbCache()

func SetDir(dir string) {
	Cache.SetDir(dir)
}

func Dir() string {
	return Cache.Dir()
}

func Db(assembly string) (*gene.GeneDb, error) {
	return Cache.Db(assembly)
}
