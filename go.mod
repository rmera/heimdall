module github.com/rmera/heimdall

go 1.22.0

require (
	github.com/rmera/chemlearn v0.0.0-20240720225553-f25e68f2bf49
	github.com/rmera/gfe v0.0.1
	github.com/rmera/gochem v0.7.1
	github.com/rmera/libsvm-go v0.0.0-20240626051442-14b6fec349c7
	github.com/rmera/scu v0.2.0
)

require (
	github.com/gonum/floats v0.0.0-20181209220543-c233463c7e82 // indirect
	github.com/gonum/internal v0.0.0-20181124074243-f884aa714029 // indirect
	github.com/skelterjohn/go.matrix v0.0.0-20130517144113-daa59528eefd // indirect
	gonum.org/v1/gonum v0.15.0 // indirect
)

replace github.com/rmera/gfe => /wrk/programs/github.com/rmera/gfe

replace github.com/rmera/chemlearn => /wrk/programs/github.com/rmera/learn
