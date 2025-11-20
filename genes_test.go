package genome_test

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-genome"
	v2 "github.com/antonybholmes/go-genome/v2"
)

func TestWithin(t *testing.T) {
	fmt.Println("Within")

	//file := fmt.Sprintf("../data/gene/%s.db", "grch38")
	db := v2.NewGeneDB("grch38", "../data/gene/")

	defer db.Close()

	location := dna.NewLocation("chr3", 187721370, 187733550)

	records, err := db.WithinGenes(location, genome.GeneFeature)

	if err != nil {
		fmt.Println(err)
		return
	}

	b, _ := json.Marshal(&records)
	fmt.Printf("%s", string(b))
}

func TestClosest(t *testing.T) {
	fmt.Println("Closest")

	//file := fmt.Sprintf("../data/gene/%s.db", "grch38")

	db := v2.NewGeneDB("grch38", "../data/gene/")

	defer db.Close()

	location := dna.NewLocation("chr3", 187721377, 187745725)

	records, err := db.ClosestGenes(location, 10)

	if err != nil {
		fmt.Println(err)
		return
	}

	b, _ := json.Marshal(&records)
	fmt.Printf("%s", string(b))
}
