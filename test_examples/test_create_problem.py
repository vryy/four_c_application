import sys

from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.FourCApplication import *

def main():
    # print(sys.argv[1:])
    # create 4C model
    fourc_problem = FourCProblem(sys.argv[1:])
    fourc_problem.ReadInputFile("contact2D_self_saddlepoint.dat", "xxx", "")

    fourc_model = fourc_problem.GetModel()

    print(fourc_problem)

def test():
    main()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()
