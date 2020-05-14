# cubedist
Calculate the pairwise distances (rmsd) between a set of electrostatic
maps, contained in Gaussian cube files, and report these distances in
CVS and a Mega-compatible formats.

## Rationale
Given a set of electrostatic maps the pairwise distances (as rmsd) are calculated. These maps are obtained by 
a Poisson Boltzmann analysis of proteins as implemented in software such as [Delphi](http://compbio.clemson.edu/delphi) 
or [APBS](http://www.poissonboltzmann.org/) (APBS writes the maps in other format, OpenDX scalar, that need to be converted 
to a Gaussian cube, with tools such as [openbabel](http://openbabel.org/wiki/Main_Page)).

For the calculated distances to make sense the maps source proteins should be similar/related (homologs, mutants,
chimeric variants, etc) and being structurally superimposed when the electrostatic potential maps are calculated.

When the sizes of a specific pair of maps are not equal, two distances are computed and averaged for
the pair, resizing each map at a time by interpolation to match the size of the other map, as the
distance only makes sense for maps with exactly the same size.

## Usage
```
cubedist [-o OUTPUT] [--version] [-h] map...

positional arguments:
  map                  electrostatic potential maps

optional arguments:
  -o, --output OUTPUT  output files base name (default:
                       "map_distances_rmsd")
  --version            show version information and exit
  -h, --help           show this help message and exit

```

## Installation
This is not a Julia package, just a standalone script. Just download it and put it in your PATH. You
need a working Julia environment and to install the dependencies.

## Dependencies
This script depends on the following Julia packages:
* `ArgParse`
* `Distances`
* `Images`
* `Printf`

