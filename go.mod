module github.com/antonybholmes/go-genes

go 1.22.2

replace github.com/antonybholmes/go-math => ../go-math

replace github.com/antonybholmes/go-dna => ../go-dna

replace github.com/antonybholmes/go-sys => ../go-sys

require github.com/antonybholmes/go-math v0.0.0-20240215163921-12bb7e52185c

require github.com/antonybholmes/go-dna v0.0.0-20240605222456-64515d092f41

require (
	github.com/antonybholmes/go-sys v0.0.0-20240505052557-9f8864ac77aa
	github.com/rs/zerolog v1.33.0
)

require (
	github.com/mattn/go-colorable v0.1.13 // indirect
	github.com/mattn/go-isatty v0.0.20 // indirect
	golang.org/x/sys v0.22.0 // indirect
)
