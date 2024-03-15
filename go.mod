module github.com/antonybholmes/go-genes

go 1.22.0

replace github.com/antonybholmes/go-math => ../go-math

replace github.com/antonybholmes/go-dna => ../go-dna

require github.com/antonybholmes/go-math v0.0.0-20240215163921-12bb7e52185c

require github.com/antonybholmes/go-dna v0.0.0-20240220233159-e0a18c04f799

require (
	github.com/mattn/go-colorable v0.1.13 // indirect
	github.com/mattn/go-isatty v0.0.20 // indirect
	github.com/rs/zerolog v1.32.0 // indirect
	golang.org/x/sys v0.18.0 // indirect
)
