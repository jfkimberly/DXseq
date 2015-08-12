from math import *


def distance(coord1, coord2):
    """calculates the norm between two points (of any dimensions)
    coord1 & coord2

    """

    dist = 0.
    for i in range(len(coord1)):
        dist += (coord1[i]-coord2[i])**2.

    return sqrt(dist)


def mm_mul(mat1, mat2):
    """returns the inner product of 2 matrices, mat1 & mat2

    """
    mat = [[0. for i in range(len(mat2[0][:]))] for j in range(len(mat1[:]))]

    for i in xrange(len(mat1[:])):
        for j in xrange(len(mat2[0][:])):
            sum = 0.
            for k in xrange(len(mat1[0][:])):
                sum += mat1[i][k]*mat2[k][j]
            mat[i][j] = sum

    return mat


def rot_matrix(unit, theta):
    """returns a rotation matrix along the 'unit' axis:
    unit[2] : 3 components of unit vector along rotation axis

    """

    R = [[0. for j in xrange(3)] for i in xrange(3)]

    R[0][0] = unit[0]*unit[0]+(1-unit[0]*unit[0])*cos(theta)
    R[0][1] = unit[0]*unit[1]*(1-cos(theta))-unit[2]*sin(theta)
    R[0][2] = unit[0]*unit[2]*(1-cos(theta))+unit[1]*sin(theta)
    R[1][0] = unit[0]*unit[1]*(1-cos(theta))+unit[2]*sin(theta)
    R[1][1] = unit[1]*unit[1]+(1-unit[1]*unit[1])*cos(theta)
    R[1][2] = unit[1]*unit[2]*(1-cos(theta))-unit[0]*sin(theta)
    R[2][0] = unit[0]*unit[2]*(1-cos(theta))-unit[1]*sin(theta)
    R[2][1] = unit[1]*unit[2]*(1-cos(theta))+unit[0]*sin(theta)
    R[2][2] = unit[2]*unit[2]+(1-unit[2]*unit[2])*cos(theta)

    return R


def transform(omega, rho, tau, R, k):
    """Transformation matrices of tilt, roll, and twist axes
    (R_tau,R_rho,R_omega), B=R_tau*R_rho*R_omega

    """
    ################################################################
    # R_tau = [[1,0,0],[0,cos(tau),sin(tau)],[0,-sin(tau),cos(tau)]]
    # R_rho = [[cos(rho),0,-sin(rho)],[0,1,0],[sin(rho),0,cos(rho)]]
    # R_ome = [[cos(ome),sin(ome),0],[-sin(ome),cos(ome),0],[0,0,1]]
    ################################################################

    A = R[:][:]
    B = [[0. for j in xrange(3)] for i in xrange(3)]

    B[0][0] = cos(omega[k])*cos(rho[k])
    B[0][1] = cos(omega[k])*sin(rho[k])*sin(tau[k])+sin(omega[k])*cos(tau[k])
    B[0][2] = -cos(omega[k])*sin(rho[k])*cos(tau[k])+sin(omega[k])*sin(tau[k])
    B[1][0] = -sin(omega[k])*cos(rho[k])
    B[1][1] = -sin(omega[k])*sin(rho[k])*sin(tau[k])+cos(omega[k])*cos(tau[k])
    B[1][2] = sin(omega[k])*sin(rho[k])*cos(tau[k])+cos(omega[k])*sin(tau[k])
    B[2][0] = sin(rho[k])
    B[2][1] = -cos(rho[k])*sin(tau[k])
    B[2][2] = cos(rho[k])*cos(tau[k])

    R = mm_mul(A, B)

    return R


def Rx(theta):
    """returns rotation matrix around x-axis"""

    Rx = [[0. for i in xrange(3)] for j in xrange(3)]

    Rx[0][0] = 1.
    Rx[1][1] = cos(theta)
    Rx[1][2] = -sin(theta)
    Rx[2][1] = sin(theta)
    Rx[2][2] = cos(theta)

    return Rx


def Ry(theta):
    """returns rotation matrix around x-axis"""

    Ry = [[0. for i in xrange(3)] for j in xrange(3)]

    Ry[0][0] = cos(theta)
    Ry[1][2] = sin(theta)
    Ry[1][1] = 1.
    Ry[2][0] = -sin(theta)
    Ry[2][2] = cos(theta)

    return Ry


def Rz(theta):
    """returns rotation matrix around z-axis"""

    Rz = [[0. for i in xrange(3)] for j in xrange(3)]

    Rz[0][0] = cos(theta)
    Rz[0][1] = -sin(theta)
    Rz[1][0] = sin(theta)
    Rz[1][1] = cos(theta)
    Rz[2][2] = 1.

    return Rz


def hp_rot(phi, theta, gamma, mat):
    """returns 'rotmat' which are rotated coordinates of 'mat'
    rotmat = Rz(phi)*Rx(theta)*Rz(gamma)*mat

    """

    rotmat = mm_mul(mm_mul(Rz(phi), mm_mul(Rx(theta), Rz(gamma))), mat)

    return rotmat


def tiletrans(coord, transvec):
    """returns the translated coordinates of 'coord' by 'transvec' """

    for i in xrange(len(coord[:])):
        for j in xrange(len(coord[0][:])):
            for k in xrange(3):
                coord[i][j][k] += transvec[k]

    return coord


def translate(coord, tiletype, direction=None):
    """returns translated coordinates of 'coord' for a given tiletype"""

    if tiletype[0] == 'S':

        for i in xrange(len(coord[:])):
            for j in xrange(len(coord[0][:])):
                coord[i][j][0] -= 30.  # 20 for exact matching
                coord[i][j][1] += 13.
                coord[i][j][2] += 125.

    if tiletype[0] == 'D':

        if direction == 'up':
            for i in xrange(len(coord[:])):
                for j in xrange(len(coord[0][:])):
                    coord[i][j][0] += 30.  # 20 for exact matching
                    coord[i][j][1] -= 3.
                    coord[i][j][2] += 125.

        if direction == 'down':
            for i in xrange(len(coord[:])):
                for j in xrange(len(coord[0][:])):
                    coord[i][j][0] -= 30.  # 20 for exact matching
                    coord[i][j][1] += 15.
                    coord[i][j][2] += 125.

        # upper B-tile bottom hairpin
        if direction == 'hp1':
            for i in xrange(len(coord[:])):
                for j in xrange(len(coord[0][:])):
                    colvec = transpose([coord[i][j]])
                    temp_coord = hp_rot(-pi/15., pi/2., pi/2., colvec)
                    coord[i][j] = [item for sublist in transpose(temp_coord)
                                   for item in sublist]

                    coord[i][j][0] += 35.  # 6.6
                    coord[i][j][1] -= 17.5
                    coord[i][j][2] += 195.

        # upper B-tile top hairpin
        if direction == 'hp2':
            for i in xrange(len(coord[:])):
                for j in xrange(len(coord[0][:])):
                    colvec = transpose([coord[i][j]])
                    temp_coord = hp_rot(-pi/10., -pi/2., pi/2., colvec)
                    coord[i][j] = [item for sublist in transpose(temp_coord)
                                   for item in sublist]

                    coord[i][j][0] += 40.  # 6.6
                    coord[i][j][1] += 3.4
                    coord[i][j][2] += 195.

        # bottom B-tile bottom hairpin
        if direction == 'hp3':
            for i in xrange(len(coord[:])):
                for j in xrange(len(coord[0][:])):
                    colvec = transpose([coord[i][j]])
                    temp_coord = hp_rot(-pi/15., pi/2., pi/2., colvec)
                    coord[i][j] = [item for sublist in transpose(temp_coord)
                                   for item in sublist]

                    coord[i][j][0] -= 25.  # 6.6
                    coord[i][j][1] -= 0.
                    coord[i][j][2] += 195.

        # bottom B-tile upper hairpin
        if direction == 'hp4':
            for i in xrange(len(coord[:])):
                for j in xrange(len(coord[0][:])):
                    colvec = transpose([coord[i][j]])
                    temp_coord = hp_rot(-pi/10., -pi/2., pi/2., colvec)
                    coord[i][j] = [item for sublist in transpose(temp_coord)
                                   for item in sublist]

                    coord[i][j][0] -= 20.
                    coord[i][j][1] += 19.4
                    coord[i][j][2] += 195.

    return coord


def transpose(lis):
    """returns transposed list of 'lis'"""

#    print zip(*lis)
    return [[row[i] for row in lis] for i in range(len(lis[0]))]


def rotateRz(coord):
    """returns rotated coordinates of 'coord' around the z-axis by PI
    i.e. coordinates for a B tile (which is Rz(pi) of an A tile)

    """

    for i in xrange(len(coord[:])):
        for j in xrange(len(coord[0][:])):
            colvec = transpose([coord[i][j]])
            colvec = mm_mul(Rx(pi/50.), colvec)
            # temp_coord = mm_mul(Ry(-pi/30.), temp_coord)
            temp_coord = mm_mul(Rz(pi), colvec)

            # flatten list; tranpose column to row vector
            coord[i][j] = [item for sublist in transpose(temp_coord)
                           for item in sublist]

    return coord


def rotRz(coord, theta):
    """returns rotated coordinates of 'coord' around the z-axis by PI
    i.e. coordinates for a B tile (which is Rz(pi) of an A tile)

    """

    for i in xrange(len(coord[:])):
        for j in xrange(len(coord[0][:])):
            colvec = transpose([coord[i][j]])
            # colvec = mm_mul(Rx(pi/50.), colvec)
            # temp_coord = mm_mul(Ry(-pi/30.), temp_coord)
            temp_coord = mm_mul(Rz(theta), colvec)

            # flatten list; tranpose column to row vector
            coord[i][j] = [item for sublist in transpose(temp_coord)
                           for item in sublist]

    return coord


def Tx(coord, dist):
    """returns 'coord' translated in the x-direction by 'dist' """

    for i in xrange(len(coord[:])):
        for j in xrange(len(coord[0][:])):
            coord[i][j][0] += dist

    return coord

if __name__ == '__main__':

    mat = [[1., 2., 3.]]
    mat = transpose(mat)
    print hp_rot(pi, pi, pi, mat)
