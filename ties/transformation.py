"""
https://github.com/synapticarbors/pyqcprot

This code code is released under the BSD 3-clause license as noted in the .pyx source code
or LICENSE file.
The original C code is copyright:
2009-2010, Pu Liu and Douglas L. Theobald
This implementation is copyright
2011, Joshua L. Adelman

  Douglas L. Theobald (2005)
  "Rapid calculation of RMSD using a quaternion-based characteristic polynomial."
  Acta Crystallographica A 61(4):478-480.

  Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2010)
  "Fast determination of the optimal rotational matrix for macromolecular superpositions."
  J. Comput. Chem. 31, 1561-1563.
"""

import copy

import numpy
from .pyqcprotext import pyqcprot


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
    rmsd = pyqcprot.CalcRMSDRotationalMatrix(frag_a, frag_b, N, rot, None)

    print('pyqcprot rmsd = ', rmsd)
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
    return coords * rotational_matrix


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
    rotation = numpy.zeros((9,), dtype=numpy.float64)

    # Calculate center of geometry
    com_ref = numpy.sum(ref, axis=1) / N
    com_mob = numpy.sum(mob_coords, axis=1) / N

    com_ref = com_ref.reshape(3, 1)
    com_mob = com_mob.reshape(3, 1)

    # Center each fragment
    ref_origin = ref - com_ref
    mob_origin = mob_coords - com_mob

    # Calculate rmsd and rotation matrix
    rmsd = pyqcprot.CalcRMSDRotationalMatrix(ref_origin, mob_origin, N, rotation, None)

    # reshape so that it can be used directly on all coordinates
    rotational_matrix = numpy.matrix(rotation.reshape((3, 3)))

    # apply the rotation to the mobile part
    rotated_mob = rotate(mob_origin.T, rotational_matrix)

    # move it back to its original position
    rotated_translated_mob = rotated_mob.T + com_mob

    # overlap_rmsd = get_rmsd(rotated_mob.T, ref_origin)

    return rmsd, (rotational_matrix, com_ref, com_mob)
