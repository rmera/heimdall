/*
 * main.go, part of Heimdall
 *
 *
 * Copyright 2024 Raul Mera  <rmeraa{at}academicos(dot)uta(dot)cl>
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *
 */

/*To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche*/

package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/rmera/boo"
	fe "github.com/rmera/gfe"
	chem "github.com/rmera/gochem"
)

// Global variables... Sometimes, you gotta use'em
var warnings []string //this counts the total warnings top level warnings printed byt the program
var verb int

// If level is larger or equal, prints the d arguments to stderr
// otherwise, does nothing.
func LogV(level int, d ...interface{}) {
	if level <= verb {
		fmt.Fprintln(os.Stderr, d...)
	}
	if level <= 1 {
		warnings = append(warnings, fmt.Sprintln(d...))
	}

}

// If level is larger or equal, prints the d arguments to stdout
// otherwise, does nothing. I could have just made one function
// with LogV, but I think this way it's more clear when reading the call
// and we can better separate "errors" from information. The plan is that you
// redirect stout and stderr to different files when calling Bartender.
func PrintV(level int, d ...interface{}) {
	if level <= verb {
		fmt.Println(d...)
	}

}

// gets a file's extension, i.e. whatever its written after the last point/dot in the filename
func getExtension(name string) string {
	fs := strings.Split(name, ".")
	return strings.ToLower(fs[len(fs)-1])
}

var PATH = os.Getenv("HEIMROOT")
var MODEL string = PATH + "/xgbmodel1.json"
var KEYS []string = []string{"fukui+var", "c6", "c6var", "charge", "chargevar", "dipolenorm", "dipolex", "dipoley", "dipolez", "elongation", "fod", "fodvar", "fukui+", "fukui-", "fukui-var", "fukui0", "fukui0var", "hardness", "homolumogap", "planarity", "sasa"}

const DIELECTRIC float64 = 80.0
const SEPARATOR string = ","

var LMAP map[int]string = map[int]string{
	0: "Negative",
	1: "Inverted",
	2: "Positive",
}

func main() {
	charge := flag.Int("c", 0, "Charge of the molecule")
	multi := flag.Int("m", 1, "Charge of the molecule")
	verbose := flag.Int("v", 1, "Level of verbosity")

	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage:\n  %s: [flags] geomtry.pdb/.gro/.xyz  \n\nFlags:\n", os.Args[0])
		flag.PrintDefaults()
	}
	flag.Parse()

	verb = *verbose

	args := flag.Args()
	if len(args) < 1 {
		log.Fatal("Reichardt requires at least 1 argument, a file with the geometry of the file to analyze")
	}
	geoname := args[0]
	var mol *chem.Molecule
	var fmap = fe.NewFeatureMap()
	var err error
	geoname = strings.Replace(geoname, "\n", "", -1)
	extension := getExtension(geoname)
	switch extension {
	case "gro":
		mol, err = chem.GroFileRead(geoname)
	case "pdb":
		mol, err = PoorlyMadePDBFileRead(geoname) //This reads normal PDBs and LigParGen "PDBs"
	default:
		mol, err = chem.XYZFileRead(geoname)
	}
	if err != nil {
		log.Fatal("Failed to open geometry input file: " + err.Error())
	}
	if extension == "xyz" {
		mol.FillIndexes()
	}
	mol.SetCharge(*charge) //needed for the MD and the partial charges calculation
	mol.SetMulti(*multi)
	var beads []*fe.Bead
	ats := make([]int, 0, mol.Len())
	w := make([]float64, 0, mol.Len())

	//We only use per-molecule features here, so all atoms go into one "bead"
	for i := 0; i < mol.Len(); i++ {
		ats = append(ats, i)
		w = append(w, 1.0)
	}
	beads = []*fe.Bead{{Indexes: ats, Weights: w}}
	o := &fe.XTBPOptions{Vari: true, Dielectric: DIELECTRIC, CFOD: false}
	xtblines := fe.XTBProps(mol.Coords[0], mol, beads, o)
	fmap.Join(xtblines)
	hardlines := fe.XTBHardness(mol.Coords[0], mol, beads, DIELECTRIC)
	fmap.Join(hardlines)
	shapes := fe.ShapeProps(mol.Coords[0], mol, beads)
	fmap.Join(shapes)
	sa := fe.SASA(mol, beads, 1)
	fmap.Join(sa)
	//	fmt.Println(fmap.String(), fmap.SortedKeys(), KEYS) ////
	floats, _, _ := fmap.IthVector(0, KEYS) //we don't care about labels here.
	//	if err != nil {
	//		panic(err)
	//	}
	f, err := os.Open("xgbmodel1.json")
	defer f.Close()
	if err != nil {
		panic(err)
	}
	b := bufio.NewReader(f)
	m, err := boo.UnJSONMultiClass(b)
	if err != nil {
		panic(err)
	}
	probs := m.PredictSingle(floats)
	class := m.PredictSingleClass(floats)
	fmt.Printf("%s is predicted to display %s solvatochromism\n", geoname, LMAP[class])
	PrintV(2, "Probabilities:")
	PrintV(2, fmt.Sprintf("positive/inverted: %3.5f, positive/negative: %3.5f inverted/negative: %3.5f", probs[2], probs[1], probs[0]))

	//	fmt.Println(LMAP[int(math.Round(returnValue))], probvalues)
}
