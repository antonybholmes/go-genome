package genedbcache

import (
	"github.com/antonybholmes/go-gene"
)

var Cache = gene.NewGeneDbCache()

func Dir(dir string) {
	Cache.Dir(dir)

}
func Db(assembly string) (*gene.GeneDb, error) {
	return Cache.Db(assembly)
}
