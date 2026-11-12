package genome

import (
	"bufio"
	"fmt"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys/log"
)

type (
	Feature struct {
		Seqname, Source, Type string
		Start, End            int
		Score, Strand, Frame  string
		Attributes            map[string][]string
		Children              []*Feature
	}

	GtfTabixDB struct {
		annotation *Annotation
		file       string
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
func (f *Feature) InExon(location *dna.Location) bool {
	// recursively get all exons using a stack
	stack := []*Feature{f}

	for len(stack) > 0 {
		current := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if current.Type == Exon {
			exonLoc, err := dna.NewLocation(f.Seqname, current.Start, current.End)
			if err == nil && exonLoc.Overlaps(location) {
				return true
			}
		}

		for _, child := range current.Children {
			stack = append(stack, child)

		}
	}

	return false
}

// InExon checks if the given location overlaps with any exon of the feature.
func (f *Feature) IsIntragenic(location *dna.Location) bool {
	// recursively get all exons using a stack
	stack := []*Feature{f}

	for len(stack) > 0 {
		current := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if current.Type == Transcript {
			exonLoc, err := dna.NewLocation(f.Seqname, current.Start, current.End)
			if err == nil && exonLoc.Overlaps(location) {
				return true
			}
		}

		for _, child := range current.Children {
			stack = append(stack, child)
		}
	}

	return false
}

func (f *Feature) InPromoter(location *dna.Location, prom *dna.PromoterRegion) bool {
	// recursively get all exons using a stack
	stack := []*Feature{f}

	for len(stack) > 0 {
		current := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if current.Type == Transcript {
			start := current.Start
			 
			if f.Strand == '-' {
				start = current.End
			}

			if f.Strand == '+' {
				end = start + prom.Downstream()
				start -= prom.Upstream()
			} else {
				end = start + prom.Upstream()
				start -= prom.Downstream()
			}

			exonLoc, err := dna.NewLocation(f.Seqname, start, end)
			
			if err == nil && exonLoc.Overlaps(location) {
				return true
			}
		}

		for _, child := range current.Children {
			stack = append(stack, child)
		}
	}

	return false
}

((g.strand = '+') AND (:start <= t.start + :prom3p) AND (:end >= t.start - :prom5p)) OR
	// 			((g.strand = '-') AND (:start <= t.end + :prom5p) AND (:end >= t.end - :prom3p))

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

func NewGtfTabixDB(dir string, annotation *Annotation) *GtfTabixDB {
	file := filepath.Join(dir, annotation.Url)

	log.Debug().Msgf("opening gene database %s", file)

	return &GtfTabixDB{annotation: annotation, file: file}
}

func (gdb *GtfTabixDB) Close() error {
	return nil
}

func (gdb *GtfTabixDB) OverlappingGenes(location *dna.Location,
	levels string,
	prom *dna.PromoterRegion,
	canonicalMode bool,
	annotationMode bool,
	biotypeFilter string) ([]*GenomicFeature, error) {

	roots, err := QueryRegion(gdb.file, location)

	if err != nil {
		return nil, err
	}

	features, err := toGenomicFeature(location, roots, levels, canonicalMode, annotationMode)

	if err != nil {
		return nil, err
	}

	return features, nil
}

// Convert gtf features from a file into standardized genomic features
func toGenomicFeature(location *dna.Location,
	feats []*Feature,
	levels string,
	canonicalMode bool,
	annotationMode bool) ([]*GenomicFeature, error) {
	var currentGene *GenomicFeature
	var currentTranscript *GenomicFeature

	var currentFeature *GenomicFeature

	var ret = make([]*GenomicFeature, 0, 10)

	gid := 1

	for _, feat := range feats {
		chr := feat.Seqname
		geneId := feat.GeneID()
		transcriptId := feat.TranscriptID()
		exonId := feat.ExonID()
		exonNumber := feat.ExonNumber()
		start := feat.Start
		end := feat.End

		geneStart := feat.Start
		geneEnd := feat.End
		strand := feat.Strand
		geneSymbol := feat.GeneName()
		featureType := feat.Type
		inExon := feat.InExon(location)
		isCanonical := feat.IsCanonical()

		//geneBiotype := feat.Biotype()

		// only add a new gene if we don't already have it. We
		// assume the rows are ordered by gene id hence if the
		// id changes, we are processing a set of rows for a new gene
		if strings.Contains(levels, "gene") && (currentGene == nil || currentGene.GeneId != geneId) {
			location, err := dna.NewStrandedLocation(chr, geneStart, geneEnd, strand)

			if err != nil {
				return nil, err
			}

			currentGene = &GenomicFeature{Id: gid,
				Location: location,
				Type:     feat.Type,
				Symbol:   geneSymbol,
				GeneId:   geneId,
				//Strand:   strand,
				//Biotype: geneBiotype,
			}

			ret = append(ret, currentGene)

			gid++
		}

		// gene inherits properties from transcripts and these become true
		// if any transcript has them true
		if currentGene != nil {
			currentGene.IsCanonical = currentGene.IsCanonical || isCanonical
			//currentGene.IsLongest = currentGene.IsLongest || isLongest

			// add extra properties if in annotation mode
			if annotationMode {
				//currentGene.InPromoter = currentGene.InPromoter || inPromoter
				currentGene.InExon = currentGene.InExon || inExon
				//currentGene.IsIntragenic = currentGene.IsIntragenic || isIntragenic
				//currentGene.Label = MakePromLabel(currentGene.InPromoter, currentGene.InExon, currentGene.IsIntragenic)

				//if (basemath.AbsInt(tssDist) < basemath.AbsInt(currentGene.TssDist)) || currentGene.TssDist == 0 {
				//	currentGene.TssDist = tssDist
				//}
			}
		}

		//if canonical mode only add if transcript is canonical
		// also only add if we have a current gene
		// also only add if we don't already have this transcript
		if strings.Contains(levels, "transcript") &&
			(currentTranscript == nil || currentTranscript.Transcript != transcriptId) &&
			(!canonicalMode || isCanonical) {

			location, err := dna.NewStrandedLocation(chr, start, end, strand)

			if err != nil {
				return nil, err
			}

			// set the properties that will not change for a transcript
			currentTranscript = &GenomicFeature{Id: gid,
				Location: location,
				//Strand:       strand,
				Type:       TranscriptLevel,
				Symbol:     geneSymbol,
				GeneId:     geneId,
				Transcript: transcriptId,
				//Biotype:      transcriptBiotype,
				IsCanonical: isCanonical,
				//IsLongest:   isLongest,
			}

			// if annotationMode {
			// 	currentTranscript.InPromoter = inPromoter
			// 	currentTranscript.IsIntragenic = isIntragenic
			// 	currentTranscript.TssDist = tssDist
			// }

			if currentGene != nil {
				if currentGene.Children == nil {
					currentGene.Children = make([]*GenomicFeature, 0, 10)
				}

				currentGene.Children = append(currentGene.Children, currentTranscript)
			} else {
				ret = append(ret, currentTranscript)
			}
		}

		if currentTranscript != nil {
			// these properties may be updated as we see more rows and we find
			// we are in an exon
			currentTranscript.InExon = currentTranscript.InExon || inExon
			currentTranscript.Label = MakePromLabel(currentTranscript.InPromoter,
				currentTranscript.InExon,
				currentTranscript.IsIntragenic)
		}

		// only add exon if we have a current transcript and it matches
		// the transcript id

		// exons are a bit more complex because they can be tagged as exon, cds, or utr and we want
		// to capture that information in the feature level. We also want to add the exon
		// to the transcript if we have a current transcript, but if not, we add it to the gene.
		// If we don't have a gene, we just add it to the return list.
		// This is because some databases may not have transcripts and we still want to
		// capture exon information if it is available.
		if strings.Contains(levels, featureType) {
			location, err := dna.NewStrandedLocation(chr, start, end, strand)

			if err != nil {
				return nil, err
			}

			currentFeature = &GenomicFeature{Id: gid,
				Location:   location,
				Type:       featureType,
				Symbol:     geneSymbol,
				GeneId:     geneId,
				Transcript: transcriptId,
				Exon:       exonId,
				ExonNumber: exonNumber,
			}

			if annotationMode {
				currentFeature.InExon = inExon
				// 	currentFeature.Label = MakePromLabel(inPromoter, inExon, isIntragenic)
			}

			if currentTranscript != nil {
				// lazy initialize children slice only if we have exons to add
				if currentTranscript.Children == nil {
					currentTranscript.Children = make([]*GenomicFeature, 0, 10)
				}

				currentTranscript.Children = append(currentTranscript.Children, currentFeature)
			} else if currentGene != nil {
				// normally add to transcript, but if no transcript, add to gene
				if currentGene.Children == nil {
					currentGene.Children = make([]*GenomicFeature, 0, 10)
				}

				currentGene.Children = append(currentGene.Children, currentFeature)
			} else {
				// no gene or transcript, just add to return list
				ret = append(ret, currentFeature)
			}
		}
	}

	return ret, nil
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
			Seqname: tokens[0], Source: tokens[1], Type: tokens[2],
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
