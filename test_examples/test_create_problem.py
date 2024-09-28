import sys

from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.FourCApplication import *
from KratosMultiphysics.DistributedBuildersApplication import *

def main():
    # print(sys.argv[1:])
    # create 4C model
    fourc_problem = FourCProblem(sys.argv[1:])
    fourc_problem.ReadInputFile("contact2D_self_saddlepoint.dat", "xxx", "")
    # print(fourc_problem.GetDiscretizationNames())
    fourc_problem.Run()

def test():
    main()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]() # allow to run test externally by python name.py test
    else:
        main()
