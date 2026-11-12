package gtftabix

import (
	"bufio"
	"errors"
	"fmt"
	"math"
	"os/exec"
	"strconv"
	"strings"

	"github.com/antonybholmes/go-dna"
	basemath "github.com/antonybholmes/go-sys/math"
)

type (
	Feature struct {
		SeqID, Source, Type  string
		Start, End           int
		Score, Strand, Frame string
		Attributes           map[string][]string
		Children             []*Feature
	}
)

const (
	InitialBufferSize = 64 * 1024
	MaxBufferSize     = 1024 * 1024

	Gene       = "gene"
	Transcript = "transcript"
	Exon       = "exon"
	CDS        = "CDS"
)

func (f *Feature) Chr() string {
	return f.SeqID
}

func (f *Feature) GeneName() string {
	return f.Attr("gene_name")
}

func (f *Feature) GeneID() string {
	return f.Attr("gene_id")
}

func (f *Feature) TranscriptID() string {
	return f.Attr("transcript_id")
}

func (f *Feature) ExonID() string {
	return f.Attr("exon_id")
}

func (f *Feature) ExonNumber() int {
	v := f.Attr("exon_number")
	if v == "" {
		return 0
	}

	n, _ := strconv.Atoi(v)

	return n
}

func (f *Feature) Exons() []*Feature {
	return f.ChildrenOf("exon")
}

// InExon checks if the given location overlaps with any exon of the feature.
func (f *Feature) Biotype() string {
	v, _ := f.FindFirst(func(current *Feature) (string, bool) {
		biotype := f.Attr("gene_type")

		if biotype != "" {
			return biotype, true
		}

		return "", false
	})

	return v

	// // recursively get all exons using a stack
	// stack := []*Feature{f}

	// for len(stack) > 0 {
	// 	current := stack[len(stack)-1]
	// 	stack = stack[:len(stack)-1]

	// 	biotype := f.Attr("gene_type")

	// 	if biotype != "" {
	// 		return biotype
	// 	}

	// 	for _, child := range current.Children {
	// 		stack = append(stack, child)
	// 	}
	// }

	// return ""
}

// Checks if the given location overlaps with any exon of the feature.
// Ideally should be run on feature that is an exon to reduce search
// time, but function will work on any feature level.
func (f *Feature) InExon(location *dna.Location) bool {
	if f.Type == Exon {
		exonLoc, err := dna.NewLocation(f.SeqID, f.Start, f.End)
		if err == nil && exonLoc.Overlaps(location) {
			return true
		}

		return false
	}

	// More comprehensive search for exons in child features

	_, ok := f.FindFirst(func(current *Feature) (string, bool) {
		if current.Type != Exon {
			return "", false
		}

		exonLoc, err := dna.NewLocation(f.SeqID, current.Start, current.End)

		if err == nil && exonLoc.Overlaps(location) {
			return "", true
		}

		// continue searching in children
		return "", false
	})

	return ok

	// // recursively get all exons using a stack
	// stack := []*Feature{f}

	// for len(stack) > 0 {
	// 	current := stack[len(stack)-1]
	// 	stack = stack[:len(stack)-1]
	// 	if current.Type == Exon {
	// 		exonLoc, err := dna.NewLocation(f.Seqname, current.Start, current.End)
	// 		if err == nil && exonLoc.Overlaps(location) {
	// 			return true
	// 		}
	// 	}

	// 	for _, child := range current.Children {
	// 		stack = append(stack, child)

	// 	}
	// }

	// return false
}

// Calculates the distance from the given location to the transcription start site (TSS) of the feature.
// It searches through all transcripts and returns the minimum absolute distance to the TSS.
// Returns an error if no transcript TSS is found.
func (f *Feature) TssDist(location *dna.Location) (int, error) {
	if f.Type != Gene && f.Type != Transcript {
		return 0, errors.New("feature is not a gene or transcript")
	}

	mid := location.Mid()

	if f.Type == Transcript {
		if f.Strand == "+" {
			return mid - f.Start, nil
		} else {
			return mid - f.End, nil
		}
	}

	// search gene for transcripts

	stack := []*Feature{f}

	minTssDist := math.MaxInt
	minAbsTssDist := math.MaxInt
	var tssDist int
	var tssAbsDist int

	for len(stack) > 0 {
		current := stack[len(stack)-1]
		stack = stack[:len(stack)-1]

		if current.Type == Transcript {
			if f.Strand == "+" {
				tssDist = mid - current.Start
			} else {
				tssDist = mid - current.End
			}

			tssAbsDist = basemath.AbsInt(tssDist)

			if tssAbsDist < minAbsTssDist {
				minAbsTssDist = tssAbsDist
				minTssDist = tssDist
			}
		}

		for _, child := range current.Children {
			stack = append(stack, child)
		}
	}

	if minAbsTssDist == math.MaxInt {
		return 0, errors.New("TSS distance not found")
	}

	return minTssDist, nil
}

// InExon checks if the given location overlaps with any exon of the feature.
func (f *Feature) IsIntragenic(location *dna.Location) (bool, error) {
	if f.Type != Gene && f.Type != Transcript {
		return false, errors.New("feature is not a gene or transcript")
	}

	if f.Type == Transcript {
		exonLoc, err := dna.NewLocation(f.Chr(), f.Start, f.End)
		if err == nil && exonLoc.Overlaps(location) {
			return true, nil
		}
		return false, nil
	}

	// Search for intragenic regions in child transcripts of the gene

	_, ok := f.FindFirst(func(current *Feature) (string, bool) {
		if current.Type != Transcript {
			return "", false
		}

		exonLoc, err := dna.NewLocation(f.Chr(), current.Start, current.End)

		if err == nil && exonLoc.Overlaps(location) {
			return "", true
		}

		return "", false
	})

	return ok, nil

}

func (f *Feature) InPromoter(location *dna.Location, prom *dna.PromoterRegion) (bool, error) {
	if f.Type != Gene && f.Type != Transcript {
		return false, errors.New("feature is not a gene or transcript")
	}

	var end int

	if f.Type == Transcript {
		start := f.Start

		if f.Strand == "-" {
			start = f.End
		}

		if f.Strand == "+" {
			end = start + prom.ThreePrime()
			start -= prom.FivePrime()
		} else {
			end = start + prom.FivePrime()
			start -= prom.ThreePrime()
		}

		exonLoc, err := dna.NewLocation(f.Chr(), start, end)

		if err == nil && exonLoc.Overlaps(location) {
			return true, nil
		}

		return false, nil
	}

	// Search for promoter regions in child transcripts of the gene

	_, ok := f.FindFirst(func(current *Feature) (string, bool) {
		if current.Type != Transcript {
			return "", false
		}

		start := current.Start

		if f.Strand == "-" {
			start = current.End
		}

		if f.Strand == "+" {
			end = start + prom.ThreePrime()
			start -= prom.FivePrime()
		} else {
			end = start + prom.FivePrime()
			start -= prom.ThreePrime()
		}

		exonLoc, err := dna.NewLocation(f.Chr(), start, end)

		if err == nil && exonLoc.Overlaps(location) {
			return "", true
		}

		return "", false
	})

	return ok, nil

}

// FindFirst checks if any feature in the hierarchy satisfies the match function.
// Terminates as soon as a match is found.
// Returns a string as this is generic and for bool like functions
// we can check if the result equals "t"
func (f *Feature) FindFirst(fn func(feature *Feature) (string, bool)) (string, bool) {
	// recursively get all exons using a stack
	stack := []*Feature{f}

	for len(stack) > 0 {
		current := stack[len(stack)-1]
		stack = stack[:len(stack)-1]

		v, ok := fn(current)
		if ok {
			return v, true
		}

		for _, child := range current.Children {
			stack = append(stack, child)
		}
	}

	return "", false
}

func (f *Feature) Length() int {
	return f.End - f.Start + 1
}

func (f *Feature) CDS() []*Feature {
	return f.ChildrenOf("CDS")
}

func (f *Feature) UTRs() []*Feature {
	return f.ChildrenOf("UTR")
}

func (f *Feature) StartCodons() []*Feature {
	return f.ChildrenOf("start_codon")
}

func (f *Feature) StopCodons() []*Feature {
	return f.ChildrenOf("stop_codon")
}

func (f *Feature) Attr(key string) string {
	if v := f.Attributes[key]; len(v) > 0 {
		return v[0]
	}

	return ""
}

func (f *Feature) IsCanonical() bool {
	for _, v := range f.Attributes["tag"] {
		if v == "Ensembl_canonical" || v == "MANE_Select" || v == "appris_principal_1" {
			return true
		}
	}

	return false
}

func (f *Feature) ChildrenOf(kind string) []*Feature {
	var out []*Feature
	for _, c := range f.Children {
		if c.Type == kind {
			out = append(out, c)
		}
	}
	return out
}

// QueryRegion queries a genomic region from a GTF or GFF file indexed with tabix.
// Returns a slice of Features representing the genomic elements overlapping the specified location.
// Each Feature will be a root such as a gene or transcript, with its children attached.
func QueryRegion(gtf string, location *dna.Location) ([]*Feature, error) {

	mode := "gtf"

	if strings.Contains(gtf, "gff") {
		mode = "gff"
	}

	cmd := exec.Command("tabix", gtf, location.String())
	stdout, err := cmd.StdoutPipe()

	if err != nil {
		return nil, err
	}

	if err := cmd.Start(); err != nil {
		return nil, err
	}

	var feats []*Feature

	sc := bufio.NewScanner(stdout)
	// make a buffer of initial size  but allow it to grow
	// up to MaxBufferSize for long attribute lines
	sc.Buffer(make([]byte, 0, InitialBufferSize), MaxBufferSize)

	for sc.Scan() {
		tokens := strings.Split(sc.Text(), "\t")

		if len(tokens) < 9 {
			continue
		}

		var s, e int
		fmt.Sscan(tokens[3], &s)
		fmt.Sscan(tokens[4], &e)

		var attrs map[string][]string

		if mode == "gff" {
			attrs = ParseGFFAttrs(tokens[8])
		} else {
			attrs = ParseGTFAttrs(tokens[8])
		}

		feats = append(feats, &Feature{
			SeqID: tokens[0], Source: tokens[1], Type: tokens[2],
			Start: s, End: e,
			Score: tokens[5], Strand: tokens[6], Frame: tokens[7],
			Attributes: attrs,
		})
	}

	if err := sc.Err(); err != nil {
		return nil, err
	}

	roots := buildTree(feats)

	return roots, cmd.Wait()
}

func buildTree(feats []*Feature) []*Feature {
	genes := map[string]*Feature{}
	trans := map[string]*Feature{}

	// pass 1: parents
	for _, f := range feats {
		switch f.Type {
		case "gene":
			genes[f.GeneID()] = f
		case "transcript":
			trans[f.TranscriptID()] = f
		}
	}

	// pass 2: link children, collect roots
	var roots []*Feature

	for _, f := range feats {
		switch f.Type {
		case "gene":
			roots = append(roots, f)
		case "transcript":
			if g := genes[f.GeneID()]; g != nil {
				g.Children = append(g.Children, f)
			} else {
				roots = append(roots, f) // gene not in region
			}
		default: // exon, CDS, UTR, codons
			if t := trans[f.TranscriptID()]; t != nil {
				t.Children = append(t.Children, f)
			} else if g := genes[f.GeneID()]; g != nil {
				g.Children = append(g.Children, f)
			} else {
				roots = append(roots, f)
			}
		}
	}

	return roots
}

// GTF attributes: gene_id "ENSG..."; gene_name "TP53";
func ParseGTFAttrs(s string) map[string][]string {
	m := map[string][]string{}

	for _, kv := range strings.Split(s, ";") {
		kv = strings.TrimSpace(kv)

		if kv == "" {
			continue
		}

		k, v, ok := strings.Cut(kv, " ")

		if !ok {
			continue
		}

		m[k] = append(m[k], strings.Trim(v, `"`))
	}

	return m
}

func ParseGFFAttrs(s string) map[string][]string {
	m := map[string][]string{}

	for _, kv := range strings.Split(s, ";") {
		kv = strings.TrimSpace(kv)

		if kv == "" {
			continue
		}

		k, v, ok := strings.Cut(kv, "=")

		if !ok {
			continue
		}

		m[k] = append(m[k], v)
	}

	return m
}

// func main() {
// 	location, err := dna.ParseLocation("chr17:7661779-7687550")
// 	if err != nil {
// 		panic(err)
// 	}

// 	feats, err := QueryRegion("annotation.gtf.gz", location, "gtf")

// 	if err != nil {
// 		panic(err)
// 	}

// 	for _, f := range feats {
// 		if f.Feature == "gene" {
// 			fmt.Println(f.Attributes["gene_name"], f.Start, f.End)
// 		}
// 	}
// }
