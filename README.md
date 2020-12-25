# Conformalized-Mean-Curvature-Flow

These are the codes for my PhD research project on computational conformal geometry. They include my own implementation of the [conformalized mean curvature flow (cMCF)](https://arxiv.org/abs/1203.6819) using the C++ mesh data structure [OpenMesh](https://www.openmesh.org). Also, it included a new algorithm "Sphericalized cMCF" that construct homotopy of degree one maps using cMCF. <br /> For an overview, please refer to my [conference proceeding](https://arxiv.org/abs/2006.07392) <br /> For more details, please refer to the chapter two in my [dissertation](https://karrywong.github.io/files/Dissertation-compressed.pdf). 

Instructions:
  1. Install mesh library [OpenMesh](https://www.openmesh.org) and numerical solver [Eigen](http://eigen.tuxfamily.org/)
  
  2. Modify directory paths in **Makefile** and ensure correct directory paths for all *.cpp* and *.h* files from the **codes** folder
  
  3. compile to get executable **Map2Sphere.exe**
  
  4. Run **Map2Sphere.exe** should give the instructions below, e.g. the code **Map2Sphere.exe -i input.off -o output.off -f 1** applies cMCF on the mesh (input.off) and output resulting mesh (output.off).

    ===================================================================================================
    ===================================================================================================
    =                                                                                                 =
    =                                         Map2Sphere                                              =
    =                                                                                                 =
    =     This program reads in a 3D (genus-zero) surface represented by a triangular mesh and        = 
    =     maps it conformally onto S^2	                                                          =
    =                                                                                                 =
    =     Usage is:                                                                                   =
    =                 Map2Sphere.exe -i INFILE -o OUTFILE                                             =
    =     where:                                                                                      =
    =                 -i INFILE        --> Input Mesh file (usually in OFF format)                    =
    =                 -o OUTFILE       --> Ouput Mesh file in OFF format                              =
    =                 -f FLOW          --> 0: MCF   1:cMCF   2:cMCF w/ projection on sphere           =
    =                                      3: Gauss Map initializer   4: Tutte Embedding initializer  =
    =                 -w Tutte weights --> 1: cotangent weights    2: graph lapacian                  =
    =                 -max             --> Max. no. of steps (integer), default: 2^6 = 64             =
    =                 -stps            --> Step size, default: 0.01				          =
    =                 -tol             --> Tolerance for sphericity, default: 0.001                   =
    =                                                                                                 =
    ===================================================================================================
    ===================================================================================================

# Copyright
OpenMesh
Copyright (c) 2001-2019, RWTH-Aachen University
Department of Computer Graphics and Multimedia
All rights reserved.
www.openmesh.org
