import numpy
import scipy.sparse
import scipy.sparse.linalg

import argparse


parser = argparse.ArgumentParser(description="""Check that the coded matrix sumn multiplication is correct.

Only works on the binary file formats for matrices (see convert-to-binary.py for conversion).

Computes X = (A+B+C) (D+E+F) using scipy and then checks if ||X - actual|| / ||X|| is small.""")

parser.add_argument("matrices", metavar="M", nargs=6,
                    help="Matrices to check summed multiplication of")

parser.add_argument("actual", type=str,
                    help="Result matrix to check")

args = parser.parse_args()


def read_matrix(filename):
    with open(filename, "rb") as f:
        mnz = numpy.fromfile(f, dtype=numpy.int32, count=3)
        nz = mnz[-1]
        coords = numpy.fromfile(f, dtype=numpy.int32, count=nz*2).reshape(nz, 2)
        v = numpy.fromfile(f, dtype=numpy.float64, count=nz)

    m, n, nz = mnz

    return scipy.sparse.coo_matrix((v, tuple(coords.T)), shape=(m, n))


A, B, C, D, E, F = map(read_matrix, args.matrices)

actual = read_matrix(args.actual)

expect = numpy.dot((A + B + C), (D + E + F))

diff = expect - actual
rnorm = scipy.sparse.linalg.norm(diff) / scipy.sparse.linalg.norm(expect)

if rnorm < 1e-6:
    print("ok: ||E - A|| / ||E|| = {}".format(rnorm))
else:
    print("fail: ||E - A|| / ||E|| = {}".format(rnorm))
