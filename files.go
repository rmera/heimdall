/*
 * files.go, part of Heimdall
 *
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

package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"regexp"
	"slices"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

//Unified function to parse the common parts of the Bartender input file.
func commonchecks(line, wanted string, linenu int, err error) string {
	//the input file format requires that the headers are given in this order
	order := []string{"BEADS", "VSITES", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER", "FEATURES"}
	if !slices.Contains(order, wanted) {
		panic("Keywork requested not in file format: " + wanted)
	}
	kwindex := slices.Index(order, wanted)
	if line == "\n" {
		return "continue"
	}
	if err != nil { //inefficient, (errs[1] can be checked once before), but clearer.
		if strings.Contains(err.Error(), "EOF") {
			return "break"
		} else {
			p := fmt.Sprintf("Failed to read line %d in input: %s", linenu, err.Error())
			panic(p)
		}
	}
	if strings.HasPrefix(line, "#") {
		return "continue" //comment line
	}

	for i, v := range order {
		if strings.HasPrefix(line, v) {
			if i > kwindex {
				return "break"
			} else if i < kwindex {
				return "continue"
			}
			return "read"
		}
	}
	return ""
}

type Bead struct {
	Indexes []int     //indexes of the atoms in the bead
	Weights []float64 //weight for each atom

}

//Parses the input file, returns an slice of functions, where the nth function will yield the coordinates of the nth bead given the current state of the atomistic molecule
//and a slice of slices of float64, where the nth slice contains, for each atom, the fraction of the
//that atom that belonging to the nth bead. Thus, a bead can contain half an atom, for instance.
func ParseInputBead(inpname string) []*Bead {
	beadslice := make([]*Bead, 0, 0)
	wslice := make([][]float64, 0, 0)
	finp, err := os.Open(inpname)
	if err != nil {
		panic("Failed to open the input file: " + err.Error())
	}
	defer finp.Close()
	inp := bufio.NewReader(finp)
	reading := false
	linenu := 0
	for {
		linenu++
		line, err := inp.ReadString('\n')
		check := commonchecks(line, "FEATURES", linenu, err)
		if check == "continue" {
			continue
		} else if check == "break" {
			break
		}
		if !reading {
			if check == "read" {
				reading = true
			}
			continue
		}
		//now the actual reading!
		pf := strings.ReplaceAll(strings.Fields(line)[1], " ", "")
		fields := strings.Split(pf, ",")
		indexes := make([]int, len(fields))
		wslice = append(wslice, make([]float64, len(fields)))
		m1 := len(wslice) - 1
		for i, v := range fields {
			var bead string
			var weight string = "1.0"
			wslice[m1][i] = 1.0
			bead = v
			if strings.Contains(v, "/") {
				info := strings.Split(v, "/")
				bead = info[0]
				weight = info[1]
			}
			wslice[m1][i], err = strconv.ParseFloat(weight, 64)
			if err != nil {
				p := fmt.Sprintf("Failed to parse the %d field in the %d line of the input file: %s", i, linenu, err.Error())
				panic(p)
			}
			wslice[m1][i] = 1 / wslice[m1][i]
			indexes[i], err = strconv.Atoi(bead)
			if err != nil {
				p := fmt.Sprintf("Failed to convert bead id: %s in line %d to int in the input file: %s", bead, linenu, err.Error())
				panic(p)
			}
			indexes[i]-- //to convert from 1-based indexes to 0-based indexes

		}
		beadslice = append(beadslice, &Bead{Indexes: indexes, Weights: wslice[m1]})

	}
	return beadslice
}

func ParseInputFeature(inpname string) []string {
	reading := false
	linenu := 0
	features := make([]string, 0, 6)
	finp, err := os.Open(inpname)
	if err != nil {
		panic("Failed to open the input file: " + err.Error())
	}
	defer finp.Close()
	inp := bufio.NewReader(finp)
	for {
		linenu++
		line, err := inp.ReadString('\n')
		check := commonchecks(line, "FEATURES", linenu, err)
		if check == "continue" {
			continue
		} else if check == "break" {
			break
		}
		if !reading && check == "read" {
			reading = true
			continue
		}

		//now the actual reading!
		pf := strings.ReplaceAll(strings.Fields(line)[1], "\n", "")
		features = append(features, pf)

	}
	return features
}

//This is a custom function that can detect and read
//the PDBs made with LigParGen, which are so wrong
//they don't actually classify as PDB at all.
//This is a pretty quick and dirty function, and I'd say quite brittle.
//It relies on LigParGen "PDBs" actually following a format where the fields
//are separated by spaces, which I am not sure is the case.
func PoorlyMadePDBFileRead(filename string) (*chem.Molecule, error) {
	fin, err := scu.NewMustReadFile(filename)
	if err != nil {
		return nil, err
	}
	symbolre := regexp.MustCompile("[a-zA-Z]+")
	first := fin.Next()
	if !strings.Contains(first, "REMARK LIGPARGEN GENERATED PDB") {
		//It it's not from LigParGen, we assume it is a reasonable PDB
		fin.Close()
		return chem.PDBFileRead(filename)
	}
	//top:=chem.NewTopology()
	//no luck, we have a LigParGen "pdb" file.
	ats := make([]*chem.Atom, 0, 10)
	tmpcoord := make([][3]float64, 0, 10)
	cont := 0
	for i := fin.Next(); i != "EOF"; i = fin.Next() {
		if !strings.HasPrefix(i, "ATOM") && !strings.HasPrefix(i, "HETATM") {
			continue
		}
		at := new(chem.Atom)
		chunks := strings.Fields(i)
		at.Symbol = symbolre.FindString(chunks[2])
		at.Name = at.Symbol
		at.MolName = chunks[3]
		at.ID = cont + 1
		at.MolID, err = strconv.Atoi(chunks[4])
		if err != nil {
			log.Printf("Couldn't obtain MolID for an atom in a LigParGen PDB: %s. Will set to 1", err.Error())
			at.MolID = 1
		}
		var coord [3]float64
		for i, c := range chunks[5:] {
			coord[i], err = strconv.ParseFloat(c, 64)
			if err != nil {
				return nil, err
			}
		}
		cont++
		tmpcoord = append(tmpcoord, coord)
		ats = append(ats, at)
	}
	coord := v3.Zeros(len(tmpcoord))
	for i, v := range tmpcoord {
		coord.Set(i, 0, v[0])
		coord.Set(i, 1, v[1])
		coord.Set(i, 2, v[2])
	}
	fin.Close()
	coords := make([]*v3.Matrix, 1)
	coords[0] = coord
	top := chem.NewTopology(0, 1, ats)
	return chem.NewMolecule(coords, top, nil)
}
