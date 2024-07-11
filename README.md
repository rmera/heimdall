# This is the Heimdall program

A Chemistry-informed machine-learning classifier for solvatochromic phenolates.

Heimdall uses chemistry-derived features obtained with [goChem](https://github.com/rmera/gochem) and [xtb](https://github.com/grimme-lab/xtb).
and a RDF-kernel Support Vector Machine, with a [locally modified](https://github.com/rmera/libsvm-go) version of the [Go implementation](https://github.com/ewalker544/libsvm-go) of [LIBSVM](https://dl.acm.org/doi/10.1145/1961189.1961199) to predict, based on the cartesian coordinates of its atoms, whether a phenolate-based solvatochromic dye will display positive, negative or
inverted solvatochromism.



Heimdall uses [xtb](https://github.com/grimme-lab/xtb). Please cite its [GFN2 methods](https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176)


LICENSE

Copyright (c) 2024 
Raul Mera  rmeraaatacademicosdotutadotcl and
Moises Dominguez moises.dominguezatusachdotcl


This program, including its documentation, 
is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation; either version 2.1 of the 
License, or (at your option) any later version.
          
This program and its documentation is distributed in the hope that 
it will be useful, but WITHOUT ANY WARRANTY; without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the GNU General Public License for more details.
                    
You should have received a copy of the GNU Lesser General 
Public License along with this program. If not, see 
<http://www.gnu.org/licenses/>. 

