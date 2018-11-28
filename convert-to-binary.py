import numpy
import argparse

parser = argparse.ArgumentParser(description="Convert the ASCII-based matrix files to a more efficient binary format")

parser.add_argument("input", type=str, help="The ASCII matrix to convert")
parser.add_argument("output", type=str, help="The name of the output file")

args = parser.parse_args()

with open(args.input, "r") as f:
    mnz = numpy.asarray(list(map(int, f.readline().split())), dtype=numpy.int32)

i, j, v = numpy.loadtxt(args.input, skiprows=1, unpack=True)

i = i.astype(numpy.int32)
j = j.astype(numpy.int32)
v = v.astype(numpy.float64)

coords = numpy.dstack([i, j])
with open(args.output, "wb") as f:
    f.write(mnz)
    f.write(coords)
    f.write(v)
