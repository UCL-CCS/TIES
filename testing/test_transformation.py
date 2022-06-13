"""
The tests focus on the transformation related to the RMSD calculation and the superimposition of molecules using the
adapted qcp package.
"""
import numpy

from ties.transformation import superimpose_coordinates


def test_transformation_default():
    """
    This test was extracted from the pyqcprot package.
    The reference value was provided for this dataset.
    """
    ref_coords = numpy.zeros((7, 3), dtype=numpy.float64)
    mob_coords = numpy.zeros((7, 3), dtype=numpy.float64)

    # create coordinates with the [atomID][xyz]
    ref_coords[0][0] = -2.803
    ref_coords[0][1] = -15.373
    ref_coords[0][2] = 24.556
    ref_coords[1][0] = 0.893
    ref_coords[1][1] = -16.062
    ref_coords[1][2] = 25.147
    ref_coords[2][0] = 1.368
    ref_coords[2][1] = -12.371
    ref_coords[2][2] = 25.885
    ref_coords[3][0] = -1.651
    ref_coords[3][1] = -12.153
    ref_coords[3][2] = 28.177
    ref_coords[4][0] = -0.440
    ref_coords[4][1] = -15.218
    ref_coords[4][2] = 30.068
    ref_coords[5][0] = 2.551
    ref_coords[5][1] = -13.273
    ref_coords[5][2] = 31.372
    ref_coords[6][0] = 0.105
    ref_coords[6][1] = -11.330
    ref_coords[6][2] = 33.567

    mob_coords[0][0] = -14.739
    mob_coords[0][1] = -18.673
    mob_coords[0][2] = 15.040
    mob_coords[1][0] = -12.473
    mob_coords[1][1] = -15.810
    mob_coords[1][2] = 16.074
    mob_coords[2][0] = -14.802
    mob_coords[2][1] = -13.307
    mob_coords[2][2] = 14.408
    mob_coords[3][0] = -17.782
    mob_coords[3][1] = -14.852
    mob_coords[3][2] = 16.171
    mob_coords[4][0] = -16.124
    mob_coords[4][1] = -14.617
    mob_coords[4][2] = 19.584
    mob_coords[5][0] = -15.029
    mob_coords[5][1] = -11.037
    mob_coords[5][2] = 18.902
    mob_coords[6][0] = -18.577
    mob_coords[6][1] = -10.001
    mob_coords[6][2] = 17.996

    # compute the rotational matrix
    rmsd, rotational_matrix, com_ref = superimpose_coordinates(ref_coords, mob_coords)

    numpy.testing.assert_allclose([rmsd], [0.719106], rtol=1e-05)
