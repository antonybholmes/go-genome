module github.com/antonybholmes/go-genes

go 1.22.5

replace github.com/antonybholmes/go-basemath => ../go-basemath

replace github.com/antonybholmes/go-dna => ../go-dna

replace github.com/antonybholmes/go-sys => ../go-sys

require github.com/antonybholmes/go-basemath v0.0.0-20240802221548-7773050a8f2f

require github.com/antonybholmes/go-dna v0.0.0-20240823140021-325275b6abfd

require (
	github.com/antonybholmes/go-sys v0.0.0-20240801224521-3bed2c519a83
	github.com/rs/zerolog v1.33.0
)

require (
	github.com/mattn/go-colorable v0.1.13 // indirect
	github.com/mattn/go-isatty v0.0.20 // indirect
	golang.org/x/sys v0.24.0 // indirect
)
