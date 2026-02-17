package genome_test

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-genome"
)

func TestWithin(t *testing.T) {
	fmt.Println("Within")

	//file := fmt.Sprintf("../data/gene/%s.db", "grch38")
	db := genome.NewGeneDB("grch38", "../data/gene/")

	defer db.Close()

	location, err := dna.NewLocation("chr3", 187721370, 187733550)

	if err != nil {
		fmt.Println(err)
		return
	}

	records, err := db.WithinGenes(location, genome.GeneLevel, dna.DefaultPromoterRegion())

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

	db := genome.NewGeneDB("grch38", "../data/gene/")

	defer db.Close()

	location, err := dna.NewLocation("chr3", 187721377, 187745725)

	if err != nil {
		fmt.Println(err)
		return
	}

	records, err := db.ClosestGenes(location, dna.DefaultPromoterRegion(), 10)

	if err != nil {
		fmt.Println(err)
		return
	}

	b, _ := json.Marshal(&records)
	fmt.Printf("%s", string(b))
}
