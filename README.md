# DXseq

## Introduction

![img](//github.com/jfkimberly/DXseq/figures/dx-hp-mol-resize.gif "An atomic (left) and a simplified (right) rendition of a DX tile with hairpins.") ![img](//github.com/jfkimberly/DXseq/figures/dx-hp-resize.gif "An atomic (left) and a simplified (right) rendition of a DX tile with hairpins.")

Want to create your own *atomically* correct double-crossover tiles? Keep
reading!

**DXseq** is a program which reads in a DNA nucleotide sequence and outputs
the spatial coordinates of each atom of each nucleotide comprising that
sequence in a [Protein Data Bank (PDB) format](http://www.rcsb.org/pdb/home/home.do). Currently, the nucleotides
are stacked in B-form (other DNA conformations may be added in the future).
Once the output file is produced, molecular visualization software can be
used to visualize the tertiary structure of the DNA.

Although the atomic coordinates of any DNA structure can be produced using
this program, at present the only available command-line options are for
various types of DX tiles (hence the name, **DXseq**).

## Dependencies

**DXseq** is written in Python 2.7. In order to run the program you must have
Python installed on your computer. This software has not been tested with
Python 3.0 or higher, so I cannot guarantee it will work under those
environments. If you do not have Python installed on your computer, you may
download the latest version at [<http://www.python.org/download/>](http://www.python.org/download/).

All python package dependencies can be installed using `pip`. `pip` does not
come pre-installed with the python virtual machine, so you need to install
it manually.

1.  Install `python-pip`, if you do not have it installed (the package name
    may be `python2-pip`, since you need python 2.7).
2.  Upgrade `pip`.
    
        $ sudo pip install pip --upgrade
3.  Install all required dependencies by typing
    
        $ pip install -r requirements.txt

Presently, you just need **numpy** (something near version 1.10.1) installed.

## Usage

At the command-line prompt, type the following (you may need to replace
`python2` with `python` depending on your system setup; this program is
written in python 2.7)

       $ python2 DXseq.py -h
    usage: DXseq.py [-h] [-o {OO,OX,XO,XX}] [-c] {CR,STL,DTL,MDX}
    
    produces PDB files of different types of DX tiles with/without hairpins
    
    positional arguments:
      {CR,STL,DTL,MDX}      choose between CR, STL, DTL, and MDX tile types.
                            CR  : a single CR-type DX tile
                            STL : two DX tiles
                            DTL : three DX tiles
                            MDX : two connected regular DX tiles
    
    optional arguments:
      -h, --help            show this help message and exit
      -o {OO,OX,XO,XX}, --options {OO,OX,XO,XX}
                            OO : complementarity & geometric compatibility (default)
                            OX : complementarity & geometric incompatibility
                            XO : noncomplementarity & geometric compatibility
                            XX : noncomplementarity & geometric incompatibility
      -c, --curvature       use experimental X-ray data for each nucleotide position
                            instead of idealized values

To create a PDB file of a normal CR-type DX tile,

    $ python2 DXseq.py CR

which produces a `CR.pdb` text file looking something like the following:

    ATOM      1  P     T a   1       7.306   2.606  -2.062  1.00  0.00           P
    ATOM      2  O1P   T a   1       7.322   3.921  -2.742  1.00  0.00           O
    ATOM      3  O2P   T a   1       8.196   2.472  -0.887  1.00  0.00           O
    ATOM      4  O5*   T a   1       5.807   2.236  -1.644  1.00  0.00           O
    ATOM      5  C5*   T a   1       4.956   1.579  -2.604  1.00  0.00           C
    ATOM      6  C4*   T a   1       4.190   0.455  -1.937  1.00  0.00           C
    ATOM      7  O4*   T a   1       4.925  -0.802  -1.951  1.00  0.00           O
    ...

Here are what some of the numbers above mean (taking the first row as an
example): 

The `ATOM` in the first column identifies the particle type, the `1` in the
second column is the atom number, the `P` in the third column is the atom
type, the `T` in the fourth column is the base type, the `a` in the fifth
column is the strand number, the `1` in the sixth column is the base number,
and the following three numbers, `7.306 2.606 -2.062`, are the *xyz* spatial
coordinates.

Use your favorite molecular viewing program to create a 3D rendering of this
file. For instance a rendering of the PDB file produced by the command

    $ python2 DXseq.py STL -o OO

using the [ViewerLite](http://accelrys.com/products/collaborative-science/biovia-discovery-studio/visualization-download.php) program is shown below.

![img](//github.com/jfkimberly/DXseq/figures/stloo-mol-resize.gif) 

![img](//github.com/jfkimberly/DXseq/figures/stloo-resize.gif)

The positional arguments `STL` and `DTL` along with the optional arguments
`-o {OO, OX, XO, XX}` produce various configurations of dimer DX tiles which
we used in [this *ACS Nano* paper](http://pubs.acs.org/doi/abs/10.1021/nn201312g).

Another feature this program offers is the calculation of the intrinsic
curvature of DNA due to its flexibility. When this is turned on (via the
`-c` flag on the command-line), the structural parameters of the DNA
(*i.e.*, twist, roll, tilt, and slide) changes from idealized values to
[values calculated from 36 independent tetramers](http://www.tandfonline.com/doi/abs/10.1080/07391102.2002.10506828) based on experimental values
from the [Nucleic Acid Database](http://www.sciencedirect.com/science/article/pii/S0006349592816491). The differences between the two can be seen
below.

![img](//github.com/jfkimberly/DXseq/figures/curvature-resize.svg)

We used this model to study the relationship between DX crystal sizes and
the structural distortion arising from the flexibility of the DX tiles in
[this *Nanotechnology* paper](http://iopscience.iop.org/article/10.1088/0957-4484/22/24/245706).

## Acknowledgements

Thanks to Seungjae Kim for sharing his version of the code.
