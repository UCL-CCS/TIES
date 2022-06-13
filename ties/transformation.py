"""

"""

import copy

import numpy
import pyqcprot as qcp


def test():
    # Setup coordinates

    frag_a = numpy.zeros((3, 7), dtype=numpy.float64)
    frag_b = numpy.zeros((3, 7), dtype=numpy.float64)
    N = 7

    frag_a[0][0] = -2.803
    frag_a[1][0] = -15.373
    frag_a[2][0] = 24.556
    frag_a[0][1] = 0.893
    frag_a[1][1] = -16.062
    frag_a[2][1] = 25.147
    frag_a[0][2] = 1.368
    frag_a[1][2] = -12.371
    frag_a[2][2] = 25.885
    frag_a[0][3] = -1.651
    frag_a[1][3] = -12.153
    frag_a[2][3] = 28.177
    frag_a[0][4] = -0.440
    frag_a[1][4] = -15.218
    frag_a[2][4] = 30.068
    frag_a[0][5] = 2.551
    frag_a[1][5] = -13.273
    frag_a[2][5] = 31.372
    frag_a[0][6] = 0.105
    frag_a[1][6] = -11.330
    frag_a[2][6] = 33.567

    frag_b[0][0] = -14.739
    frag_b[1][0] = -18.673
    frag_b[2][0] = 15.040
    frag_b[0][1] = -12.473
    frag_b[1][1] = -15.810
    frag_b[2][1] = 16.074
    frag_b[0][2] = -14.802
    frag_b[1][2] = -13.307
    frag_b[2][2] = 14.408
    frag_b[0][3] = -17.782
    frag_b[1][3] = -14.852
    frag_b[2][3] = 16.171
    frag_b[0][4] = -16.124
    frag_b[1][4] = -14.617
    frag_b[2][4] = 19.584
    frag_b[0][5] = -15.029
    frag_b[1][5] = -11.037
    frag_b[2][5] = 18.902
    frag_b[0][6] = -18.577
    frag_b[1][6] = -10.001
    frag_b[2][6] = 17.996

    # Allocate rotation array
    rot = numpy.zeros((9,), dtype=numpy.float64)

    # Calculate center of geometry
    comA = numpy.sum(frag_a, axis=1) / N
    comB = numpy.sum(frag_b, axis=1) / N

    # Center each fragment
    frag_a = frag_a - comA.reshape(3, 1)
    frag_b = frag_b - comB.reshape(3, 1)

    # Calculate rmsd and rotation matrix
    rmsd = qcp.CalcRMSDRotationalMatrix(frag_a, frag_b, N, rot, None)

    print('qcp rmsd = ', rmsd)
    print('rotation matrix:')
    print(rot.reshape((3, 3)))


    # Calculate rmsd after applying rotation
    def rmsd(a, b):
        """Returns RMSD between two coordinate sets a and b."""
        return numpy.sqrt(numpy.sum(numpy.power(a - b, 2)) / a.shape[1])


    # rotate frag_b to obtain optimal alignment
    frag_br = frag_b.T * numpy.matrix(rot.reshape((3, 3)))
    rmsd = rmsd(frag_br.T, frag_a)
    print('rmsd after applying rotation: ', rmsd)


# Calculate rmsd after applying rotation
def get_rmsd(a, b):
    """Returns RMSD between two coordinate sets a and b."""
    return numpy.sqrt(numpy.sum(numpy.power(a - b, 2)) / a.shape[1])


def rotate(coords, rotational_matrix):
    # rotate coords using the rotational matrix
    return coords.T * numpy.matrix(rotational_matrix.reshape((3, 3)))


def superimpose_coordinates(ref, mob):
    """
    Superimpose the mob structure to the ref structure and return the necessary rotational matrix etc.

    :param ref: coordinates for the reference structure
    :param mob: coordaintes for the mobile structure
    :return: (rot, rmsd, overlap_rmsd) with rot being a rotational matrix that has to be applied to the
        mobile structure in order to superimpose it with the reference structure.
        Overlap_rmsd is the rmsd value obtained after applying the rotational matrix.
    """
    if len(ref) != len(mob):
        raise Exception('Coordinates to be overlapped must have the same number of atoms')
    if ref.shape[1] != mob.shape[1] != 3:
        raise Exception('The coordinates should be have the array shape (n, 3)')

    N = len(ref)

    ref = copy.deepcopy(ref.T)
    mob_coords = copy.deepcopy(mob.T)

    # Allocate rotation array
    rotation_matrix = numpy.zeros((9,), dtype=numpy.float64)

    # Calculate center of geometry
    comA = numpy.sum(ref, axis=1) / N
    comB = numpy.sum(mob_coords, axis=1) / N

    # Center each fragment
    ref_origin = ref - comA.reshape(3, 1)
    mob_origin = mob_coords - comB.reshape(3, 1)

    # Calculate rmsd and rotation matrix
    rmsd = qcp.CalcRMSDRotationalMatrix(ref_origin, mob_origin, N, rotation_matrix, None)

    # apply the rotation to the mobile part
    rotated_mob = rotate(mob_origin, rotation_matrix)

    # move it back to its original position
    rotated_translated_mob = rotated_mob.T + comB.reshape(3, 1)

    # overlap_rmsd = get_rmsd(rotated_mob.T, ref_origin)

    return rmsd, rotated_translated_mob.T
