package genome

import "github.com/antonybholmes/go-web"

type (

	//string string
	//string string

	Annotation struct {
		PublicId string `json:"id"`
		Genome   string `json:"genome"`
		Assembly string `json:"assembly"`
		Type     string `json:"type"`
		Name     string `json:"name"`
		Url      string `json:"-"`
		Id       int    `json:"-"`
	}
)

// Normalize assembly names, e.g. grch37 -> hg19 so that
// users can use either name and get the same result
func NormalizeAssembly(assembly string) string {

	assembly = web.FormatParam(assembly)

	if assembly == "grch37" {
		assembly = "hg19"
	}

	if assembly == "grcm38" {
		assembly = "mm10"
	}

	return assembly
}
