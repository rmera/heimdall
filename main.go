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
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"strings"

	fe "github.com/rmera/gfe"
	chem "github.com/rmera/gochem"
	libSvm "github.com/rmera/libsvm-go"
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

//gets a file's extension, i.e. whatever its written after the last point/dot in the filename
func getExtension(name string) string {
	fs := strings.Split(name, ".")
	return strings.ToLower(fs[len(fs)-1])
}

var PATH = os.Getenv("HEIMROOT")
var MODEL string = PATH + "/RM1.svm.model"
var RANGE string = PATH + "/RM1.svm.range"

var KEPT []string = []string{"dipolenorm", "chargevar", "hardness", "planarity", "elongation", "sasa"}

const DIELECTRIC float64 = 80.0
const SEPARATOR string = ","

var LMAP map[int]string = map[int]string{
	1: "Negative",
	2: "Inverted",
	3: "Positive",
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
	beads = []*fe.Bead{&fe.Bead{Indexes: ats, Weights: w}}

	o := &fe.XTBPOptions{Vari: true, Dielectric: DIELECTRIC, CFOD: false}
	xtblines := fe.XTBProps(mol.Coords[0], mol, beads, o)
	fmap.Join(xtblines)

	hardlines := fe.XTBHardness(mol.Coords[0], mol, beads, DIELECTRIC)
	fmap.Join(hardlines)
	shapes := fe.ShapeProps(mol.Coords[0], mol, beads)
	fmap.Join(shapes)

	sa := fe.SASA(mol, beads, 1)
	fmap.Join(sa)
	PrintV(4, "Before Deleting Items", fmap.String())

	deleted := fe.ExcludeKeys(fmap.SortedKeys(), KEPT)

	PrintV(4, "Keys for deletion", deleted) /////////////
	fmap.DeleteFeatures(deleted)

	const TMPNAME string = "problem.tmp"
	tmp, err := os.Create(TMPNAME)
	PrintV(4, "Final map used:", fmap.String())
	floats, _, _ := fmap.IthVector(0, KEPT) //we don't care about labels here.
	str := fe.Floats2SVM(floats, 0)
	tmp.WriteString(str)
	tmp.Close()
	PrintV(3, "The SVM string", str)
	model := libSvm.NewModelFromFile(MODEL)
	PrintV(3, "And it keys", fe.Strings2SVM(KEPT))
	model.Probability(true)
	minmax, lu, err := libSvm.ReadRangesFromFile(RANGE)
	if err != nil {
		log.Fatalf("Unable to read from Scaling file %s. Error: %v", RANGE, err)
	}
	param := libSvm.NewParameter()
	prob, err := libSvm.NewProblem(TMPNAME, param)
	if err != nil {
		log.Fatalf("Unable to read model from tmp file %s. Error: %v", RANGE, err)
	}
	PrintV(3, "Scaling Factors", minmax, lu)
	label, line := prob.GetLine()
	line, err = fe.ScaleMap(line, minmax, lu...)
	PrintV(4, "Label read:", label)
	PrintV(4, "Data read", line)
	returnValue, vals := model.PredictValues(line)
	v32 := vals[0]
	v31 := vals[1]
	v21 := vals[1]
	fmt.Printf("%s is predicted to display %s solvatochromism\n", geoname, LMAP[int(math.Round(returnValue))])
	PrintV(2, "Decision values:")
	PrintV(2, fmt.Sprintf("positive/inverted: %3.5f, positive/negative: %3.5f inverted/negative: %3.5f", v32, v31, v21))
	PrintV(2, "Positive values favor the first of the 2 categories, in each case, while negative values favor the second one.")

	//	fmt.Println(LMAP[int(math.Round(returnValue))], probvalues)
}
