#!/usr/bin/python
import argparse
import os
import string
from math import *

import numpy as np

from linalg import *


# globally defined constants
STACK_NUM = 42  # length of base pair/ length of cylinder
RISE_ZERO = 3.4
PI = pi/180.

RHO_ZERO = PI
TAU_ZERO = 0.
SLIDE_ZERO = 0.2
JNUM = 4
DIR = os.getcwd()
k = 1
m = 1
#junction = np.array([12, 13, 28, 29])

###############################################################################
# read-in strand sequence
parser = argparse.ArgumentParser(description="produces PDB files of different \
types of DX tiles with/without hairpins",
                                 formatter_class=argparse.RawTextHelpFormatter)


parser.add_argument("tiletype", choices=['CR', 'STL', 'DTL', 'MDX'],
                    help="""choose between CR, STL, DTL, and MDX tile types.
CR  : a single CR-type DX tile
STL : two DX tiles
DTL : three DX tiles
MDX : two connected regular DX tiles""")


parser.add_argument("-o", "--options", default='OO',
                    choices=['OO', 'OX', 'XO', 'XX'],
                    help="""OO : complementarity &\
 geometric compatibility (default)
OX : complementarity & geometric incompatibility
XO : noncomplementarity & geometric compatibility
XX : noncomplementarity & geometric incompatibility
""")


parser.add_argument("-c", "--curvature", action="store_true",
                    help="""use experimental X-ray data for each nucleotide position
instead of idealized values""")


args = parser.parse_args()

if args.tiletype == 'STL':
    tiletype = 'STL' + args.options
    tiletype2 = None

elif args.tiletype == 'DTL':
    tiletype = 'DTL' + args.options + '-1'
    tiletype2 = 'DTL' + args.options + '-2'

elif args.tiletype == 'CR':
    tiletype = 'CR'
    tiletype2 = None

elif args.tiletype == 'MDX':
    tiletype = 'MDX'
    tiletype2 = None

# remove existing pdb files if they exist
#if os.path.exists(DIR + r'/' + tiletype + '.pdb'):
#    os.remove(tiletype + '.pdb')
# 
#pdbfile = open(tiletype + ".pdb", "a")

###############################################################################
# assign strand sequences to duplexes

seq = []

with open(DIR + r"/DATA/sequence.txt") as seqfile:
    for line in seqfile:
        seq.append(line.rstrip())

dup1 = seq[seq.index(tiletype)+1]
dup2 = seq[seq.index(tiletype)+2]

if tiletype[0] == 'M':
    dup3 = seq[seq.index(tiletype)+3]
    dup4 = seq[seq.index(tiletype)+4]

# read-in hairpin sequence for ABB-type systems
if tiletype[0] == 'D':

    dup3 = seq[seq.index(tiletype)+4]
    dup4 = seq[seq.index(tiletype)+5]

    if tiletype != 'DTLOO-1':

        hp1 = seq[seq.index(tiletype)+6]
        hp2 = seq[seq.index(tiletype)+7]
        hp3 = hp1[:]
        hp4 = hp2[:]

        for i in range(STACK_NUM-len(hp1)):
            hp1 += 'T'
        for i in range(STACK_NUM-len(hp2)):
            hp2 += 'T'
        for i in range(STACK_NUM-len(hp3)):
            hp3 += 'T'
        for i in range(STACK_NUM-len(hp4)):
            hp4 += 'T'

###############################################################################
# read-in parameters(omega,rho,tau,slide) of each base pairs

omegatable = []
rhotable = []
tautable = []
slidetable = []

with open(DIR + r'/DATA/domega.txt') as domega:
    for line in domega:
        omegatable.append(float(line.rstrip()))

with open(DIR + r'/DATA/drho.txt') as drho:
    for line in drho:
        rhotable.append(float(line.rstrip()))

with open(DIR + r'/DATA/dtau.txt') as dtau:
    for line in dtau:
        tautable.append(float(line.rstrip()))

with open(DIR + r'/DATA/dslide.txt') as dslide:
    for line in dslide:
        slidetable.append(float(line.rstrip()))

# read-in atom parameters
atom_name = []
with open(DIR + r"/DATA/atoms2.txt") as atomname_f:
    for line in atomname_f:
        atom_name.append(line.rstrip())
###############################################################################


def fbp_pos():
    """returns the coordinates of the first base pair,
    same as 'coordinates' function in C except for index position
    coordinates[base pair index][atom index][x,y,z index]

    """

    coordinates = np.zeros((STACK_NUM, 164, 3))
    temp_coordinates = []

    with open(DIR + r'/DATA/AT.txt') as AT:
        for line in AT:
            temp_coordinates.append(map(float, line.rstrip().split()))

    with open(DIR + r'/DATA/GC.txt') as GC:
        for line in GC:
            temp_coordinates.append(map(float, line.rstrip().split()))

    with open(DIR + r'/DATA/CG.txt') as CG:
        for line in CG:
            temp_coordinates.append(map(float, line.rstrip().split()))

    with open(DIR + r'/DATA/TA.txt') as TA:
        for line in TA:
            temp_coordinates.append(map(float, line.rstrip().split()))

    coordinates[0] = temp_coordinates[:][:]

    return coordinates


def calculate(coordinates, omega, rho, tau, slide, axis):
    """calculates the position of all atoms for all base pairs
    with the atoms of the first base pair coordinates as the initial
    coordinates

    """

    R = np.identity(3, float)

    sum_R = np.zeros(3)
    sum_D = np.zeros(3)

    # iteration between base pairs
    for i in range(STACK_NUM-1):

        R = transform(omega, rho, tau, R, i)

        # normal vector of base pairs
        sum_R[0] = sum_R[0]+R[2, 0]
        sum_R[1] = sum_R[1]+R[2, 1]
        sum_R[2] = sum_R[2]+R[2, 2]

        # slide (D_y), slide[i] is the slide distance
        sum_D[0] = sum_D[0]+slide[i]*R[1, 0]
        sum_D[1] = sum_D[1]+slide[i]*R[1, 1]
        sum_D[2] = sum_D[2]+slide[i]*R[1, 2]

        ############################################################
        # leave the axis out for now
        for k in range(3):
            axis[i+1, k] = axis[0, 0]*R[0, k] + axis[0, 1]*R[1, k] +\
                axis[0, 2]*R[2, k]

        for k in range(3):
            axis[i+1, k] += RISE_ZERO*(sum_R[k])+sum_D[k]
        ############################################################

        # coordinates(1) = coordinates(0) * R_0
        # coordinates(2) = coordinates(1) * R_1 = coordinates(0) * R_0 * R_1
        # etc.

        for j in range(164):
            for k in range(3):
                for l in range(3):
                    coordinates[i+1, j, k] += coordinates[0, j, l]*R[l, k]

                coordinates[i+1, j, k] += RISE_ZERO*(sum_R[k])+sum_D[k]

    return coordinates, axis


def caseselect(duplex, omega, rho, tau, slide):

    duplex_c = "RABC"
    duplex_x = "RY"

    for i in range(STACK_NUM-3):
        # the direction is always from 5' to 3',
        # up one way and down the other backbone.

        duplex_c = list(duplex_c)

        if duplex[i] == 'A' or duplex[i] == 'G':
            duplex_c[0] = duplex_x[0]

        if duplex[i] == 'T' or duplex[i] == 'C':
            duplex_c[0] = duplex_x[1]

        duplex_c[1] = duplex[i+1]
        duplex_c[2] = duplex[i+2]

        if duplex[i+3] == 'A' or duplex[i+3] == 'G':
            duplex_c[3] = duplex_x[0]

        if duplex[i+3] == 'T' or duplex[i+3] == 'C':
            duplex_c[3] = duplex_x[1]

        duplex_c = ''.join(duplex_c)

        if duplex_c == "RAAR" or duplex_c == "YTTY":
            casenum = 0
        elif duplex_c == "RAAY" or duplex_c == "RTTY":
            casenum = 1
        elif duplex_c == "YAAR" or duplex_c == "YTTR":
            casenum = 2
        elif duplex_c == "YAAY" or duplex_c == "RTTR":
            casenum = 3
        elif duplex_c == "RAGR" or duplex_c == "YCTY":
            casenum = 4
        elif duplex_c == "RAGY" or duplex_c == "RCTY":
            casenum = 5
        elif duplex_c == "YAGR" or duplex_c == "YCTR":
            casenum = 6
        elif duplex_c == "YAGY" or duplex_c == "RCTR":
            casenum = 7
        elif duplex_c == "RGAR" or duplex_c == "YTCY":
            casenum = 8
        elif duplex_c == "RGAY" or duplex_c == "RTCY":
            casenum = 9
        elif duplex_c == "YGAR" or duplex_c == "YTCR":
            casenum = 10
        elif duplex_c == "YGAY" or duplex_c == "RTCR":
            casenum = 11
        elif duplex_c == "RGGR" or duplex_c == "YCCY":
            casenum = 12
        elif duplex_c == "RGGY" or duplex_c == "RCCY":
            casenum = 13
        elif duplex_c == "YGGR" or duplex_c == "YCCR":
            casenum = 14
        elif duplex_c == "YGGY" or duplex_c == "RCCR":
            casenum = 15
        elif duplex_c == "RCAR" or duplex_c == "YTGY":
            casenum = 16
        elif duplex_c == "RCAY" or duplex_c == "RTGY":
            casenum = 17
        elif duplex_c == "YCAR" or duplex_c == "YTGR":
            casenum = 18
        elif duplex_c == "YCAY" or duplex_c == "RTGR":
            casenum = 19
        elif duplex_c == "RCGR" or duplex_c == "YCGY":
            casenum = 20
        elif duplex_c == "RCGY":
            casenum = 21
        elif duplex_c == "YCGR":
            casenum = 22
        elif duplex_c == "RTAR" or duplex_c == "YTAY":
            casenum = 23
        elif duplex_c == "RTAY":
            casenum = 24
        elif duplex_c == "YTAR":
            casenum = 25
        elif duplex_c == "RATR" or duplex_c == "YATY":
            casenum = 26
        elif duplex_c == "RATY":
            casenum = 27
        elif duplex_c == "YATR":
            casenum = 28
        elif duplex_c == "RACR" or duplex_c == "YGTY":
            casenum = 29
        elif duplex_c == "RACY" or duplex_c == "RGTY":
            casenum = 30
        elif duplex_c == "YACR" or duplex_c == "YGTR":
            casenum = 31
        elif duplex_c == "YACY" or duplex_c == "RGTR":
            casenum = 32
        elif duplex_c == "RGCR" or duplex_c == "YGCY":
            casenum = 33
        elif duplex_c == "RGCY":
            casenum = 34
        elif duplex_c == "YGCR":
            casenum = 35

        if args.curvature:
            omega[i+1] = OMEGA_ZERO+omegatable[casenum]*PI
            rho[i+1] = RHO_ZERO+rhotable[casenum]*PI
            tau[i+1] = TAU_ZERO+tautable[casenum]*PI
            slide[i+1] = SLIDE_ZERO+slidetable[casenum]

        else:
            omega[i+1] = OMEGA_ZERO
            rho[i+1] = RHO_ZERO
            tau[i+1] = TAU_ZERO
            slide[i+1] = SLIDE_ZERO

    omega[0] = OMEGA_ZERO
    rho[0] = RHO_ZERO
    tau[0] = TAU_ZERO
    slide[0] = SLIDE_ZERO
    omega[STACK_NUM-2] = OMEGA_ZERO
    rho[STACK_NUM-2] = RHO_ZERO
    tau[STACK_NUM-2] = TAU_ZERO
    slide[STACK_NUM-2] = SLIDE_ZERO

    return omega, rho, tau, slide


def backbone_modify(duplex, coordinates, duplexnum):
    """modifies the coordinates of the sugar and phosphate atoms of the DNA
    so as to create the correct visual bonds in viewing programs

    """

    theta = 20.*PI
    angle = 2.5
    R = np.identity(3, float)
    STACK_NUM = 7
    distcom = np.zeros(3)
    acom = np.zeros(3)

    # iterate over all base pairs
    # from 5' -> 3' up one side of the backbone
    for i in range(STACK_NUM):

        if duplex[i] == 'A':

            theta = angle*PI

            if (
                    (duplexnum == 1 and (i == 7 or i == 8)) or
                    (duplexnum == 2 and (i == 23 or i == 24))
            ):
                theta = -angle*PI

            ############################################################
            # create unit vector 'unit' along the axis of rotation

            for j in range(3):
                distcom[j] = coordinates[i, 10, j]-coordinates[i, 11, j]

            # sum = distance(distcom, np.zeros(3))
            sum = np.linalg.norm(distcom)

            # 'unit' is a unit vector along the rotation axis
            unit = []
            for j in distcom:
                unit.append(j/sum)

            R = rot_matrix(unit, theta)

            # check 'atoms2.txt' for atoms from 0<=j<10
            for j in range(10):

                # index changed
                for k in range(3):
                    acom[k] = coordinates[i, j, k] - coordinates[i, 11, k]

                # index changed
                for k in range(3):
                    coordinates[i, j, k] =\
                        acom[0]*R[k, 0] + acom[1]*R[k, 1] +\
                        acom[2]*R[k, 2] + coordinates[i, 11, k]

        elif duplex[i] == 'G':

            theta = angle*PI

            if (
                    (duplexnum == 1 and (i == 7 or i == 8)) or
                    (duplexnum == 2 and (i == 23 or i == 24))
            ):
                theta = -angle*PI

            ############################################################
            # create unit vector 'unit' along the axis of rotation

            for j in range(3):
                distcom[j] = coordinates[i, 51, j]-coordinates[i, 52, j]

#            sum = distance(distcom, [0., 0., 0.])
            sum = np.linalg.norm(distcom)

            unit = []
            for j in distcom:
                unit.append(j/sum)

            R = rot_matrix(unit, theta)

            # check 'atoms2.txt' for atoms from 41<=j<51
            for j in range(41, 51):

                for k in range(3):
                    acom[k] = coordinates[i, j, k] - coordinates[i, 52, k]

                for k in range(3):
                    coordinates[i, j, k] =\
                        acom[0]*R[k, 0] + acom[1]*R[k, 1] +\
                        acom[2]*R[k, 2] + coordinates[i, 52, k]

        elif duplex[i] == 'C':

            theta = angle*PI

            if (
                    (duplexnum == 1 and (i == 7 or i == 8)) or
                    (duplexnum == 2 and (i == 23 or i == 24))
            ):
                theta = -angle*PI

            ############################################################
            # create unit vector 'unit' along the axis of rotation

            for j in range(3):
                distcom[j] = coordinates[i, 92, j] - coordinates[i, 93, j]

#            sum = distance(distcom, [0., 0., 0.])
            sum = np.linalg.norm(distcom)

            unit = []
            for j in distcom:
                unit.append(j/sum)

            R = rot_matrix(unit, theta)

            for j in range(82, 92):

                for k in range(3):
                    acom[k] = coordinates[i, j, k] - coordinates[i, 93, k]

                for k in range(3):
                    coordinates[i, j, k] =\
                        acom[0]*R[k, 0] + acom[1]*R[k, 1] +\
                        acom[2]*R[k, 2] + coordinates[i, 93, k]

        else:

            theta = angle*PI

            if (
                    (duplexnum == 1 and (i == 7 or i == 8)) or
                    (duplexnum == 2 and (i == 23 or i == 24))
            ):
                theta = -angle*PI

            ############################################################
            # create unit vector 'unit' along the axis of rotation

            for j in range(3):
                distcom[j] = coordinates[i, 133, j] - coordinates[i, 134, j]

#            sum = distance(distcom, [0., 0., 0.])
            sum = np.linalg.norm(distcom)

            unit = []
            for j in distcom:
                unit.append(j/sum)

            R = rot_matrix(unit, theta)

            for j in range(123, 133):

                for k in range(3):
                    acom[k] = coordinates[i, j, k] - coordinates[i, 134, k]

                for k in range(3):
                    coordinates[i, j, k] =\
                        acom[0]*R[k, 0] + acom[1]*R[k, 1] +\
                        acom[2]*R[k, 2] + coordinates[i, 134, k]

    # iterate over all base pairs
    # from 3' -> 5' up the other side of the backbone (against the grain)
    for i in range(STACK_NUM):

        if duplex[i] == 'A':

            theta = angle*PI

            if (
                    (duplexnum == 2 and (i == 7 or i == 8)) or
                    (duplexnum == 1 and (i == 23 or i == 24))
            ):
                theta = -angle*PI

            for j in range(3):
                distcom[j] = coordinates[i, 31, j]-coordinates[i, 32, j]

#            sum = distance(distcom, [0., 0., 0.])
            sum = np.linalg.norm(distcom)

            unit = []
            for j in distcom:
                unit.append(j/sum)

            R = rot_matrix(unit, theta)

            # check 'atoms2.txt' for atoms from 21<=j<31
            for j in range(21, 31):

                for k in range(3):
                    acom[k] = coordinates[i, j, k] - coordinates[i, 32, k]

                for k in range(3):
                    coordinates[i, j, k] =\
                        acom[0]*R[k, 0] + acom[1]*R[k, 1] +\
                        acom[2]*R[k, 2] + coordinates[i, 32, k]

        elif duplex[i] == 'G':

            theta = angle*PI

            if (
                    (duplexnum == 2 and (i == 7 or i == 8)) or
                    (duplexnum == 1 and (i == 23 or i == 24))
            ):
                theta = -angle*PI

            for j in range(3):
                distcom[j] = coordinates[i, 73, j]-coordinates[i, 74, j]

#            sum = distance(distcom, [0., 0., 0.])
            sum = np.linalg.norm(distcom)

            unit = []
            for j in distcom:
                unit.append(j/sum)

            R = rot_matrix(unit, theta)

            # check 'atoms2.txt' for atoms from 63<=j<73
            for j in range(63, 73):

                for k in range(3):
                    acom[k] = coordinates[i, j, k] - coordinates[i, 74, k]

                for k in range(3):
                    coordinates[i, j, k] =\
                        acom[0]*R[k, 0] + acom[1]*R[k, 1] +\
                        acom[2]*R[k, 2] + coordinates[i, 74, k]

        elif duplex[i] == 'C':

            theta = angle*PI

            if (
                    (duplexnum == 2 and (i == 7 or i == 8)) or
                    (duplexnum == 1 and (i == 23 or i == 24))
            ):
                theta = -angle*PI

            for j in range(3):
                distcom[j] = coordinates[i, 111, j]-coordinates[i, 112, j]

#            sum = distance(distcom, [0., 0., 0.])
            sum = np.linalg.norm(distcom)

            unit = []
            for j in distcom:
                unit.append(j/sum)

            R = rot_matrix(unit, theta)

            for j in range(101, 111):

                for k in range(3):
                    acom[k] = coordinates[i, j, k] - coordinates[i, 112, k]

                for k in range(3):
                    coordinates[i, j, k] =\
                        acom[0]*R[k, 0] + acom[1]*R[k, 1] +\
                        acom[2]*R[k, 2] + coordinates[i, 112, k]

        # T-base
        else:

            theta = angle*PI

            if (
                    (duplexnum == 2 and (i == 7 or i == 8)) or
                    (duplexnum == 1 and (i == 23 or i == 24))
            ):
                theta = -angle*PI

            for j in range(3):
                distcom[j] = coordinates[i, 153, j]-coordinates[i, 154, j]

#            sum = distance(distcom, [0., 0., 0.])
            sum = np.linalg.norm(distcom)

            unit = []
            for j in distcom:
                unit.append(j/sum)

            R = rot_matrix(unit, theta)

            for j in range(143, 153):

                for k in range(3):
                    acom[k] = coordinates[i, j, k] - coordinates[i, 154, k]

                for k in range(3):
                    coordinates[i, j, k] =\
                        acom[0]*R[k, 0] + acom[1]*R[k, 1] +\
                        acom[2]*R[k, 2] + coordinates[i, 154, k]

    return coordinates


def junction_builder(coord1, coord2, dup1, dup2, juncnum1, juncnum2, check):

    dist1 = 9.
    dist2 = 10.

    R = [[int(i == j) for i in range(3)] for j in range(3)]
    R2 = [[int(i == j) for i in range(3)] for j in range(3)]

    if check == 53:
        if dup1[juncnum1] == 'A':
            m = 8
        elif dup1[juncnum1] == 'G':
            m = 49
        elif dup1[juncnum1] == 'C':
            m = 90
        else:
            m = 131

        if dup2[juncnum2] == 'A':
            n = 21
        elif dup2[juncnum2] == 'G':
            n = 63
        elif dup2[juncnum2] == 'C':
            n = 101
        else:
            n = 143

    elif check == 35:
        if dup1[juncnum1] == 'A':
            m = 29
        elif dup1[juncnum1] == 'G':
            m = 71
        elif dup1[juncnum1] == 'C':
            m = 109
        else:
            m = 151

        if dup2[juncnum2] == 'A':
            n = 0
        elif dup2[juncnum2] == 'G':
            n = 41
        elif dup2[juncnum2] == 'C':
            n = 82
        else:
            n = 123

    elif check == 55:
        if dup1[juncnum1] == 'A':
            m = 8
        elif dup1[juncnum1] == 'G':
            m = 49
        elif dup1[juncnum1] == 'C':
            m = 90
        else:
            m = 131

        if dup2[juncnum2] == 'A':
            n = 0
        elif dup2[juncnum2] == 'G':
            n = 41
        elif dup2[juncnum2] == 'C':
            n = 82
        else:
            n = 123

    elif check == 33:
        if dup1[juncnum1] == 'A':
            m = 29
        elif dup1[juncnum1] == 'G':
            m = 71
        elif dup1[juncnum1] == 'C':
            m = 109
        else:
            m = 151

        if dup2[juncnum2] == 'A':
            n = 21
        elif dup2[juncnum2] == 'G':
            n = 63
        elif dup2[juncnum2] == 'C':
            n = 101
        else:
            n = 143

    x1 = []
    x2 = []

    for i in range(3):
        x1.append(coord2[juncnum2, n+4, i])
        x2.append(coord2[juncnum2, n+5, i])

    sum2 = distance(x1, x2)

    unit2 = []
    for i in range(3):
        unit2.append((x1[i]-x2[i])/sum2)

    theta2 = 1.*PI
    R2 = rot_matrix(unit2, theta2)

    while dist1 < dist2:

        dist2 = dist1

        for k in range(n, n+4):

            acom2 = []
            for l in range(3):
                acom2.append(coord2[juncnum2, k, l]-coord2[juncnum2, n+5, l])

            for l in range(3):
                coord2[juncnum2, k, l] = acom2[0]*R2[l, 0] + \
                    acom2[1]*R2[l, 1] + acom2[2]*R2[l, 2] + \
                    coord2[juncnum2, n+5, l]

        c1 = []
        c2 = []
        for i in range(3):
            c1.append(coord1[juncnum1, m, i])
            c2.append(coord2[juncnum2, n, i])
        dist1 = distance(c1, c2)

    x1 = []
    x2 = []
    for i in range(3):
        x1.append(coord2[juncnum2, n+3, i])
        x2.append(coord2[juncnum2, n+4, i])
    sum = distance(x1, x2)

    unit = []
    for i in range(3):
        unit.append((x1[i]-x2[i])/sum)

    theta = 1.*PI
    R = rot_matrix(unit, theta)

    dist1 = 9.
    dist2 = 10.

    while(dist1 < dist2):
        dist2 = dist1

        for k in range(n, n+3):

            acom2 = []
            for l in range(3):
                acom2.append(coord2[juncnum2, k, l]-coord2[juncnum2, n+4, l])

            for l in range(3):
                coord2[juncnum2, k, l] = acom2[0]*R[l, 0] +\
                    acom2[1]*R[l, 1] + acom2[2]*R[l, 2] +\
                    coord2[juncnum2, n+4, l]

        c1 = []
        c2 = []
        for i in range(3):
            c1.append(coord1[juncnum1, m, i])
            c2.append(coord2[juncnum2, n, i])
        dist1 = distance(c1, c2)


def makedxtile(coord1, coord2, dup1, dup2, junction):
    bin = 11.
    junction_bond = np.zeros(4)
    moving_distance = np.zeros(3)
    # bondrange = 1.9

    for m in range(3):

        md1 = coord1[junction[0], 8, m]
        md2 = coord2[junction[0], 29, m]
        md3 = coord1[junction[1], 8, m]
        md4 = coord2[junction[1], 29, m]
        md5 = coord1[junction[2], 29, m]
        md6 = coord2[junction[2], 8, m]
        md7 = coord1[junction[3], 29, m]
        md8 = coord2[junction[3], 8, m]

        moving_distance[m] = md1+md3+md5+md7-(md2+md4+md6+md8)

#    md_sum = distance(moving_distance, [0., 0., 0.])
    md_sum = np.linalg.norm(moving_distance)

    for n in range(3):

        moving_distance[n] *= bin/md_sum
        # print "MD\t %f\n" % moving_distance[n]

    for i in range(STACK_NUM):
        for j in range(164):
            for k in range(3):
                coord1[i, j, k] -= moving_distance[k]
                coord2[i, j, k] += moving_distance[k]

    junction_builder(
        coord1, coord2, dup1, dup2, junction[0], junction[0], 53)
    junction_builder(
        coord2, coord1, dup2, dup1, junction[1], junction[1], 35)
    junction_builder(
        coord2, coord1, dup2, dup1, junction[2], junction[2], 53)
    junction_builder(
        coord1, coord2, dup1, dup2, junction[3], junction[3], 35)
    junction_builder(
        coord1, coord1, dup1, dup1, junction[0]-1, junction[0], 55)
    junction_builder(
        coord2, coord2, dup1, dup1, junction[3]+1, junction[3], 33)
    junction_builder(
        coord2, coord2, dup2, dup2, junction[1]+1, junction[1], 55)
    junction_builder(
        coord1, coord1, dup2, dup2, junction[2]-1, junction[2], 33)

    for i in range(4):
        if i % 2 == 0:

            # print "1st\n"
            if dup1[junction[i]] == 'A':
                m = 8
            elif dup1[junction[i]] == 'G':
                m = 49
            elif dup1[junction[i]] == 'C':
                m = 90
            else:
                m = 131

            if dup2[junction[i]] == 'A':
                n = 21
            elif dup2[junction[i]] == 'G':
                n = 63
            elif dup2[junction[i]] == 'C':
                n = 101
            else:
                n = 143

        else:

            # print "2nd\n"
            if dup1[junction[i]] == 'A':
                m = 29
            elif dup1[junction[i]] == 'G':
                m = 71
            elif dup1[junction[i]] == 'C':
                m = 109
            else:
                m = 151

            if dup2[junction[i]] == 'A':
                n = 0
            elif dup2[junction[i]] == 'G':
                n = 41
            elif dup2[junction[i]] == 'C':
                n = 82
            else:
                n = 123

        x1 = []
        x2 = []
        for j in range(3):
            x1.append(coord1[junction[i], m, j])
            x2.append(coord2[junction[i], n, j])

        junction_bond[i] = distance(x1, x2)

#    while junction_bond[0]>bondrange or junction_bond[1]>bondrange or\
#            junction_bond[2]>bondrange or junction_bond[3]>bondrange:

    for i in range(STACK_NUM):
        for j in range(164):
            for k in range(3):

                coord2[i, j, k] -= 0.1*moving_distance[k]
                coord1[i, j, k] += 0.1*moving_distance[k]

    for i in range(4):
        if i == 0:

            if dup1[junction[i]] == 'A':
                m = 8
            elif dup1[junction[i]] == 'G':
                m = 49
            elif dup1[junction[i]] == 'C':
                m = 90
            else:
                m = 131

            if dup2[junction[i]] == 'A':
                n = 21
            elif dup2[junction[i]] == 'G':
                n = 63
            elif dup2[junction[i]] == 'C':
                n = 101
            else:
                n = 143

        elif i == 3:
            if dup1[junction[i]] == 'A':
                m = 29
            elif dup1[junction[i]] == 'G':
                m = 71
            elif dup1[junction[i]] == 'C':
                m = 109
            else:
                m = 151

            if dup2[junction[i]] == 'A':
                n = 0
            elif dup2[junction[i]] == 'G':
                n = 41
            elif dup2[junction[i]] == 'C':
                n = 82
            else:
                n = 123

        elif i == 2:
            if dup2[junction[i]] == 'A':
                m = 8
            elif dup2[junction[i]] == 'G':
                m = 49
            elif dup2[junction[i]] == 'C':
                m = 90
            else:
                m = 131

            if dup1[junction[i]] == 'A':
                n = 21
            elif dup1[junction[i]] == 'G':
                n = 63
            elif dup1[junction[i]] == 'C':
                n = 101
            else:
                n = 143

        else:  # i == 1
            if dup2[junction[i]] == 'A':
                m = 29
            elif dup2[junction[i]] == 'G':
                m = 71
            elif dup2[junction[i]] == 'C':
                m = 109
            else:
                m = 151

            if dup1[junction[i]] == 'A':
                n = 0
            elif dup1[junction[i]] == 'G':
                n = 41
            elif dup1[junction[i]] == 'C':
                n = 82
            else:
                n = 123

        x1 = []
        x2 = []
        for j in range(3):
            x1.append(coord1[junction[i], m, j])
            x2.append(coord2[junction[i], n, j])
        junction_bond[i] = distance(x1, x2)

    return coord1, coord2, moving_distance


def filewrite(coord, init_atomnum, fin_atomnum, strdnum, base, norm, i):

    global k, m

    if norm:
        for j in range(init_atomnum, fin_atomnum):

            pdbfile.write(('ATOM' + ' '*2 + '%5d' + ' '*2 + '%-3s' + ' '*3 +
                           '%s' + ' ' + '%1s%4d' + ' '*4 + '%8.3f'*3 + ' '*2 +
                           '1.00' + ' '*2 + '0.00' + ' '*11 + '%s\n') %
                          (k, atom_name[j], base, strdnum, m, coord[i, j, 0],
                           coord[i, j, 1], coord[i, j, 2], atom_name[j][0]))

            k += 1
        m += 1

    else:
        if base == 'A':
            base = 'T'
        elif base == 'G':
            base = 'C'
        elif base == 'C':
            base = 'G'
        elif base == 'T':
            base = 'A'

        for j in range(init_atomnum, fin_atomnum):

            pdbfile.write(('ATOM' + ' '*2 + '%5d' + ' '*2 + '%-3s' + ' '*3 +
                           '%s' + ' ' + '%1s%4d' + ' '*4 + '%8.3f'*3 + ' '*2 +
                           '1.00' + ' '*2 + '0.00' + ' '*11 + '%s\n') %
                          (k, atom_name[j], base, strdnum, m, coord[i, j, 0],
                           coord[i, j, 1], coord[i, j, 2], atom_name[j][0]))

            k += 1
        m += 1


def DXoutput(coord1, coord2, tilenum, hpcoord1, hpcoord2, strdnum):

    # map letters to numbers
    for i in range(len(string.lowercase)):
        if strdnum == string.lowercase[i]:
            strdnum2 = i+1

    if strdnum2 % 4 == 1:

        # strand 1
        # first sequence
        ibasenum = 0
        fbasenum = junction[1]
        inc = 1
        i = ibasenum

        for base in dup1[ibasenum:fbasenum:inc]:

            if base == 'A':
                filewrite(coord1, 0, 21, strdnum, base, True, i)
            elif base == 'G':
                filewrite(coord1, 41, 63, strdnum, base, True, i)
            elif base == 'C':
                filewrite(coord1, 82, 101, strdnum, base, True, i)
            elif base == 'T':
                filewrite(coord1, 123, 143, strdnum, base, True, i)
            i += inc

        # strand 1
        # second sequence
#        ibasenum = junction[0]
#        fbasenum = None
#        inc = -1
#        i = ibasenum

        ibasenum = 0
        fbasenum = junction[1]
        inc = 1
        i = ibasenum

        strdnum = '2'

        for base in dup1[ibasenum:fbasenum:inc]:

            if base == 'A':
                filewrite(coord2, 21, 41, strdnum, base, False, i)
            elif base == 'G':
                filewrite(coord2, 63, 82, strdnum, base, False, i)
            elif base == 'C':
                filewrite(coord2, 101, 123, strdnum, base, False, i)
            elif base == 'T':
                filewrite(coord2, 143, 164, strdnum, base, False, i)
            i += inc

    elif strdnum2 % 4 == 2:

        # strand 2
        # first sequence
        ibasenum = STACK_NUM-1
        fbasenum = junction[2]
        inc = -1
        i = ibasenum

        if tiletype == 'STLOX' or tiletype == 'STLXX':
            ibasenum = STACK_NUM-6
            fbasenum = junction[2]
            inc = -1
            i = ibasenum

        for base in dup1[ibasenum:fbasenum:inc]:

            if base == 'A':
                filewrite(coord1, 21, 41, strdnum, base, False, i)
            elif base == 'G':
                filewrite(coord1, 63, 82, strdnum, base, False, i)
            elif base == 'C':
                filewrite(coord1, 101, 123, strdnum, base, False, i)
            elif base == 'T':
                filewrite(coord1, 143, 164, strdnum, base, False, i)
            i += inc

        # strand 2
        # second sequence
        ibasenum = junction[3]
        fbasenum = STACK_NUM
        inc = 1
        i = ibasenum

        if (tiletype == 'STLOX' or tiletype == 'STLXX')\
           or ((tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2')
               and tilenum == 'B'):
            ibasenum = junction[3]
            fbasenum = STACK_NUM-5
            inc = 1
            i = ibasenum

        for base in dup2[ibasenum:fbasenum:inc]:

            if base == 'A':
                filewrite(coord2, 0, 21, strdnum, base, True, i)
            elif base == 'G':
                filewrite(coord2, 41, 63, strdnum, base, True, i)
            elif base == 'C':
                filewrite(coord2, 82, 101, strdnum, base, True, i)
            elif base == 'T':
                filewrite(coord2, 123, 143, strdnum, base, True, i)
            i += inc

    elif strdnum2 % 4 == 3:

        # strand 3
        # first sequence
        ibasenum = 5
        fbasenum = junction[3]
        inc = 1
        i = ibasenum

        if (tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2')\
           and tilenum == 'B':
            ibasenum = 0
            fbasenum = junction[3]-8
            inc = 1
            i = ibasenum

        elif tiletype2 == 'DTLXO-2' and tilenum == 'B':
            ibasenum = 5
            fbasenum = junction[3]-8
            inc = 1
            i = ibasenum

        for base in dup2[ibasenum:fbasenum:inc]:

            if base == 'A':
                filewrite(coord2, 0, 21, strdnum, base, True, i)
            elif base == 'G':
                filewrite(coord2, 41, 63, strdnum, base, True, i)
            elif base == 'C':
                filewrite(coord2, 82, 101, strdnum, base, True, i)
            elif base == 'T':
                filewrite(coord2, 123, 143, strdnum, base, True, i)
            i += inc


###############################################################################
# hairpin1
        if (
                tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXO-2'
                or tiletype2 == 'DTLXX-2'
        ) and tilenum == 'B':
            # 5' -> 3' up the hairpin

            ibasenum = 2
            fbasenum = 12
            inc = 1
            i = ibasenum

            for base in hp1[ibasenum:fbasenum:inc]:

                if base == 'A':
                    filewrite(hpcoord1, 0, 21, strdnum, base, True, i)
                elif base == 'G':
                    filewrite(hpcoord1, 41, 63, strdnum, base, True, i)
                elif base == 'C':
                    filewrite(hpcoord1, 82, 101, strdnum, base, True, i)
                elif base == 'T':
                    filewrite(hpcoord1, 123, 143, strdnum, base, True, i)
                i += inc

            # 5' -> 3' down the hairpin
            ibasenum = 11
            fbasenum = None
            inc = -1
            i = ibasenum

            for base in hp1[ibasenum:fbasenum:inc]:

                if base == 'A':
                    filewrite(hpcoord1, 21, 41, strdnum, base, False, i)
                elif base == 'G':
                    filewrite(hpcoord1, 63, 82, strdnum, base, False, i)
                elif base == 'C':
                    filewrite(hpcoord1, 101, 123, strdnum, base, False, i)
                elif base == 'T':
                    if i == 10 or i == 11:
                        filewrite(hpcoord1, 143, 164, strdnum, base, True, i)
                    else:
                        filewrite(hpcoord1, 143, 164, strdnum, base, False, i)
                i += inc

            # end of hairpin1
            # rest of strand
            ibasenum = junction[3]-8
            fbasenum = junction[3]
            inc = 1
            i = ibasenum

            for base in dup2[ibasenum:fbasenum:inc]:

                if base == 'A':
                    filewrite(coord2, 0, 21, strdnum, base, True, i)
                elif base == 'G':
                    filewrite(coord2, 41, 63, strdnum, base, True, i)
                elif base == 'C':
                    filewrite(coord2, 82, 101, strdnum, base, True, i)
                elif base == 'T':
                    filewrite(coord2, 123, 143, strdnum, base, True, i)
                i += inc

###############################################################################

        # strand 3
        # second sequence
        ibasenum = junction[2]
        fbasenum = 4
        inc = -1
        i = ibasenum

        for base in dup1[ibasenum:fbasenum:inc]:

            if base == 'A':
                filewrite(coord1, 21, 41, strdnum, base, False, i)
            elif base == 'G':
                filewrite(coord1, 63, 82, strdnum, base, False, i)
            elif base == 'C':
                filewrite(coord1, 101, 123, strdnum, base, False, i)
            elif base == 'T':
                filewrite(coord1, 143, 164, strdnum, base, False, i)
            i += inc

    elif strdnum2 % 4 == 0:

        # strand 4
        # first sequence
        ibasenum = STACK_NUM-6
        fbasenum = junction[0]
        inc = -1
        i = ibasenum

        if tiletype == 'STLOX' or tiletype == 'STLXX':
            ibasenum = STACK_NUM-1
            fbasenum = junction[0]
            inc = -1
            i = ibasenum

        elif (tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2')\
                and tilenum == 'B':
            ibasenum = STACK_NUM-1
            fbasenum = junction[0]+8
            inc = -1
            i = ibasenum

        elif tiletype2 == 'DTLXO-2' and tilenum == 'B':
            ibasenum = STACK_NUM-6
            fbasenum = junction[0]+8
            inc = -1
            i = ibasenum

        for base in dup2[ibasenum:fbasenum:inc]:

            if base == 'A':
                filewrite(coord2, 21, 41, strdnum, base, False, i)
            elif base == 'G':
                filewrite(coord2, 63, 82, strdnum, base, False, i)
            elif base == 'C':
                filewrite(coord2, 101, 123, strdnum, base, False, i)
            elif base == 'T':
                filewrite(coord2, 143, 164, strdnum, base, False, i)
            i += inc

###############################################################################
# hairpin2
        if (
                tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2'
                or tiletype2 == 'DTLXO-2'
        ) and tilenum == 'B':
            # 5' -> 3' up the hairpin
            ibasenum = 2
            fbasenum = 12
            inc = 1
            i = ibasenum

            for base in hp2[ibasenum:fbasenum:inc]:

                if base == 'A':
                    filewrite(hpcoord2, 0, 21, strdnum, base, True, i)
                elif base == 'G':
                    filewrite(hpcoord2, 41, 63, strdnum, base, True, i)
                elif base == 'C':
                    filewrite(hpcoord2, 82, 101, strdnum, base, True, i)
                elif base == 'T':
                    filewrite(hpcoord2, 123, 143, strdnum, base, True, i)
                i += inc

            # 5' -> 3' down the hairpin
            ibasenum = 11
            fbasenum = None
            inc = -1
            i = ibasenum

            for base in hp2[ibasenum:fbasenum:inc]:

                if base == 'A':
                    filewrite(hpcoord2, 21, 41, strdnum, base, False, i)
                elif base == 'G':
                    filewrite(hpcoord2, 63, 82, strdnum, base, False, i)
                elif base == 'C':
                    filewrite(hpcoord2, 101, 123, strdnum, base, False, i)
                elif base == 'T':
                    if i == 10 or i == 11:
                        filewrite(hpcoord2, 143, 164, strdnum, base, True, i)
                    else:
                        filewrite(hpcoord2, 143, 164, strdnum, base, False, i)
                i += inc

            # end of hairpin2
            # rest of strand
            ibasenum = junction[0]+8
            fbasenum = junction[0]
            inc = -1
            i = ibasenum

            for base in dup2[ibasenum:fbasenum:inc]:
                if base == 'A':
                    filewrite(coord2, 21, 41, strdnum, base, False, i)
                elif base == 'G':
                    filewrite(coord2, 63, 82, strdnum, base, False, i)
                elif base == 'C':
                    filewrite(coord2, 101, 123, strdnum, base, False, i)
                elif base == 'T':
                    filewrite(coord2, 143, 164, strdnum, base, False, i)
                i += inc

###############################################################################

        # strand 4
        # second sequence
        ibasenum = junction[1]
        fbasenum = STACK_NUM-5
        inc = 1
        i = ibasenum

        if tiletype == 'STLOX' or tiletype == 'STLXX':
            ibasenum = junction[1]
            fbasenum = STACK_NUM
            inc = 1
            i = ibasenum

        for base in dup1[ibasenum:fbasenum:inc]:

            if base == 'A':
                filewrite(coord1, 0, 21, strdnum, base, True, i)
            elif base == 'G':
                filewrite(coord1, 41, 63, strdnum, base, True, i)
            elif base == 'C':
                filewrite(coord1, 82, 101, strdnum, base, True, i)
            elif base == 'T':
                filewrite(coord1, 123, 143, strdnum, base, True, i)
            i += inc


if __name__ == '__main__':

    ############################################################
    # initialize parameters
    omega = np.zeros(STACK_NUM)
    rho = np.zeros(STACK_NUM)
    tau = np.zeros(STACK_NUM)
    slide = np.zeros(STACK_NUM)
    axis1 = np.zeros((STACK_NUM, 3))
    axis2 = np.zeros((STACK_NUM, 3))
    mvdist = np.zeros(3)

    junction = np.array([22, 23, 28, 29])


    ############################################################
    original_coord = fbp_pos()
    coord1 = np.copy(original_coord)
    coord2 = np.copy(original_coord)
    hpcoord1 = np.copy(original_coord)
    hpcoord2 = np.copy(original_coord)
    hpcoord3 = np.copy(original_coord)
    hpcoord4 = np.copy(original_coord)

    if tiletype[0] == 'S' or tiletype[0] == 'C':

        #######################################################################
        OMEGA_ZERO = 720.*PI/21.  # used in caseselect function
        duplexnum = 1
        omega, rho, tau, slide = caseselect(dup1, omega, rho, tau, slide)
        coordinates, axis1 = calculate(coord1, omega, rho, tau, slide, axis1)
        coord1 = backbone_modify(dup1, coordinates, duplexnum)

        OMEGA_ZERO = -720.*PI/21.
        duplexnum = 1
        omega, rho, tau, slide = caseselect(dup1, omega, rho, tau, slide)
        coordinates2, axis2 = calculate(coord2, omega, rho, tau, slide, axis2)
        coord2 = backbone_modify(dup1, coordinates2, duplexnum)

#        coordinates[base pair index][atom index][x,y,z index]

        # A tile
        tilenum = 'A'

#        for theta in np.arange(0, 45.*PI, 30.*PI):

#        for theta in np.arange(30.*PI, 45.*PI, 30.*PI):
        theta = 330.*PI
        R = np.array([[cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0,0,1]]) 
        for bp_index in range(coord2.shape[0]):
            for atom_index in range(coord2.shape[1]):
                coord2[bp_index, atom_index] = np.dot(R, coord2[bp_index, atom_index])
        theta = int(ceil(theta / PI))
        pdbfile = open('CR' + str(theta) + '.pdb', 'w')
        DXoutput(coord1, coord2, tilenum, None, None, 'a')

#        DXoutput(coord1, coord2, tilenum, None, None, 'b')
#        DXoutput(coord1, coord2, tilenum, None, None, 'c')
#        DXoutput(coord1, coord2, tilenum, None, None, 'd')

        if tiletype[0] != 'C':
            # B tile
            tilenum = 'B'
            DXoutput(coord3, coord4, tilenum, None, None, 'e')
            DXoutput(coord3, coord4, tilenum, None, None, 'f')
            DXoutput(coord3, coord4, tilenum, None, None, 'g')
            DXoutput(coord3, coord4, tilenum, None, None, 'h')
    ###########################################################################

    if tiletype[0] == 'D':
        #######################################################################
        duplexnum = 1
        omega, rho, tau, slide = caseselect(dup1, omega, rho, tau, slide)
        coordinates, axis1 = calculate(coord1, omega, rho, tau, slide, axis1)
        coord1 = backbone_modify(dup1, coordinates, duplexnum)

        duplexnum = 2
        omega, rho, tau, slide = caseselect(dup2, omega, rho, tau, slide)
        coordinates2, axis2 = calculate(coord2, omega, rho, tau, slide, axis2)
        coord2 = backbone_modify(dup2, coordinates2, duplexnum)

        coord1, coord2, mvdist =\
            makedxtile(coord1, coord2, dup1, dup2, junction)

        coord3 = np.copy(original_coord)
        coord4 = np.copy(original_coord)

        duplexnum = 1
        omega, rho, tau, slide = caseselect(dup1, omega, rho, tau, slide)
        coordinates, axis1 = calculate(coord3, omega, rho, tau, slide, axis1)
        coord3 = backbone_modify(dup1, coordinates, duplexnum)
        coord3 = rotateRz(coord3)
        coord3 = translate(coord3, tiletype2, 'up')

        duplexnum = 2
        omega, rho, tau, slide = caseselect(dup2, omega, rho, tau, slide)
        coordinates2, axis2 = calculate(coord4, omega, rho, tau, slide, axis2)
        coord4 = backbone_modify(dup2, coordinates2, duplexnum)
        coord4 = rotateRz(coord4)
        coord4 = translate(coord4, tiletype2, 'up')

        coord3, coord4, mvdist =\
            makedxtile(coord3, coord4, dup1, dup2, junction)

        coord5 = np.copy(original_coord)
        coord6 = np.copy(original_coord)

        duplexnum = 1
        omega, rho, tau, slide = caseselect(dup3, omega, rho, tau, slide)
        coordinates, axis1 = calculate(coord5, omega, rho, tau, slide, axis1)
        coord5 = backbone_modify(dup3, coordinates, duplexnum)
        coord5 = rotateRz(coord5)
        coord5 = translate(coord5, tiletype2, 'down')

        duplexnum = 2
        omega, rho, tau, slide = caseselect(dup4, omega, rho, tau, slide)
        coordinates2, axis2 = calculate(coord6, omega, rho, tau, slide, axis2)
        coord6 = backbone_modify(dup4, coordinates2, duplexnum)
        coord6 = rotateRz(coord6)
        coord6 = translate(coord6, tiletype2, 'down')

        coord5, coord6, mvdist =\
            makedxtile(coord5, coord6, dup3, dup4, junction)

        # A tile
        tilenum = 'A'
        DXoutput(coord1, coord2, tilenum, None, None, 'a')
        DXoutput(coord1, coord2, tilenum, None, None, 'b')
        DXoutput(coord1, coord2, tilenum, None, None, 'c')
        DXoutput(coord1, coord2, tilenum, None, None, 'd')

        #######################################################################
        if tiletype != 'DTLOO-1':  # add hairpins

            omega, rho, tau, slide = caseselect(hp1, omega, rho, tau, slide)
            coordinates, axis1 = calculate(
                hpcoord1, omega, rho, tau, slide, axis1)
            hpcoord1 = backbone_modify(hp1, coordinates, duplexnum)
            hpcoord1 = translate(hpcoord1, tiletype, 'hp1')

            omega, rho, tau, slide = caseselect(hp2, omega, rho, tau, slide)
            coordinates2, axis2 = calculate(
                hpcoord2, omega, rho, tau, slide, axis2)
            hpcoord2 = backbone_modify(hp2, coordinates2, duplexnum)
            hpcoord2 = translate(hpcoord2, tiletype, 'hp2')
        #######################################################################

        # B tile
        tilenum = 'B'
        DXoutput(coord3, coord4, tilenum, hpcoord1, hpcoord2, 'e')
        DXoutput(coord3, coord4, tilenum, hpcoord1, hpcoord2, 'f')
        DXoutput(coord3, coord4, tilenum, hpcoord1, hpcoord2, 'g')
        DXoutput(coord3, coord4, tilenum, hpcoord1, hpcoord2, 'h')

        if tiletype != 'DTLOO-1':  # add hairpins

            omega, rho, tau, slide = caseselect(hp3, omega, rho, tau, slide)
            coordinates, axis1 = calculate(
                hpcoord3, omega, rho, tau, slide, axis1)
            hpcoord3 = backbone_modify(hp3, coordinates, duplexnum)
            hpcoord3 = translate(hpcoord3, tiletype, 'hp3')

            omega, rho, tau, slide = caseselect(hp4, omega, rho, tau, slide)
            coordinates2, axis2 = calculate(
                hpcoord4, omega, rho, tau, slide, axis2)
            hpcoord4 = backbone_modify(hp4, coordinates2, duplexnum)
            hpcoord4 = translate(hpcoord4, tiletype, 'hp4')

        # B tile
        tilenum = 'B'
        DXoutput(coord5, coord6, tilenum, hpcoord3, hpcoord4, 'i')
        DXoutput(coord5, coord6, tilenum, hpcoord3, hpcoord4, 'j')
        DXoutput(coord5, coord6, tilenum, hpcoord3, hpcoord4, 'k')
        DXoutput(coord5, coord6, tilenum, hpcoord3, hpcoord4, 'l')

        # hairpins
        #######################################################################

    if tiletype[0] == 'M':

        duplexnum = 1
        omega, rho, tau, slide = caseselect(dup1, omega, rho, tau, slide)
        coordinates, axis1 = calculate(coord1, omega, rho, tau, slide, axis1)
        coord1 = backbone_modify(dup1, coordinates, duplexnum)

        duplexnum = 2
        omega, rho, tau, slide = caseselect(dup2, omega, rho, tau, slide)
        coordinates2, axis2 = calculate(coord2, omega, rho, tau, slide, axis2)
        coord2 = backbone_modify(dup2, coordinates2, duplexnum)

        coord1, coord2, mvdist =\
            makedxtile(coord1, coord2, dup1, dup2, junction)

        for i in range(STACK_NUM):
            for j in range(3):
                axis1[i, j] -= mvdist[j]
                axis2[i, j] += mvdist[j]

        transvec = []
        transvec.append(axis2[-5, 0] - axis1[0, 0])
        transvec.append(axis2[-5, 1] - axis1[0, 1])
        transvec.append(axis2[-5, 2] - axis1[0, 2])

        coord3 = np.copy(original_coord)
        coord4 = np.copy(original_coord)

        duplexnum = 1
        omega, rho, tau, slide = caseselect(dup3, omega, rho, tau, slide)
        coordinates, axis1 = calculate(coord3, omega, rho, tau, slide, axis1)
        coord3 = backbone_modify(dup3, coordinates, duplexnum)
        coord3 = rotateRz(coord3)
        coord3 = tiletrans(coord3, transvec)

        duplexnum = 2
        omega, rho, tau, slide = caseselect(dup4, omega, rho, tau, slide)
        coordinates2, axis2 = calculate(coord4, omega, rho, tau, slide, axis2)
        coord4 = backbone_modify(dup4, coordinates2, duplexnum)
        coord4 = rotateRz(coord4)
        coord4 = tiletrans(coord4, transvec)

        coord3, coord4, mvdist =\
            makedxtile(coord3, coord4, dup3, dup4, junction)

        # A tile
        tilenum = 'A'
        DXoutput(coord1, coord2, tilenum, None, None, 'a')
        DXoutput(coord1, coord2, tilenum, None, None, 'b')
        DXoutput(coord1, coord2, tilenum, None, None, 'c')
        DXoutput(coord1, coord2, tilenum, None, None, 'd')

        # B tile
        tilenum = 'B'
        DXoutput(coord3, coord4, tilenum, None, None, 'e')
        DXoutput(coord3, coord4, tilenum, None, None, 'f')
        DXoutput(coord3, coord4, tilenum, None, None, 'g')
        DXoutput(coord3, coord4, tilenum, None, None, 'h')
        #######################################################################

    pdbfile.write("END\n")
