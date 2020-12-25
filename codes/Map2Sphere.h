/* ===============================================================================================
   Map2Sphere: Conformal mapping of a genus zero surface represented by a triangular mesh onto
               the sphere.
	       Different methods are implemented:
		1) Conformalized Mean Curvature Flow
  
   Author:  Karry Wong
   Date:    04/04/2020
   Version: 2
   =============================================================================================== */

#ifndef _MAP2SPHERE_
#define _MAP2SPHERE_

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <ctime>

/* ===============================================================================================
   Includes associated to Eigen or CHOLMOD
   =============================================================================================== */

#include <Eigen/Core>
#include <Eigen/Sparse>

/* ===============================================================================================
   Includes associated to OpenMesh
   =============================================================================================== */

#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"


/* ===============================================================================================
   Usage
   =============================================================================================== */

static void usage(char** argv)
{
    std::cout << "\n\n" <<std::endl;
    std::cout << "     " << "==================================================================================================="<<std::endl;
    std::cout << "     " << "==================================================================================================="<<std::endl;
    std::cout << "     " << "=                                                                                                 ="<<std::endl;
    std::cout << "     " << "=                                         Map2Sphere                                              ="<<std::endl;
    std::cout << "     " << "=                                                                                                 ="<<std::endl;
    std::cout << "     " << "=     This program reads in a 3D (genus-zero) surface represented by a triangular mesh and maps it conformally onto S^2	="<<std::endl;
    std::cout << "     " << "=                                                                                                 ="<<std::endl;
    std::cout << "     " << "=     Usage is:                                                                                   ="<<std::endl;
    std::cout << "     " << "=                 Map2Sphere.exe -i INFILE -o OUTFILE                                             ="<<std::endl;
    std::cout << "     " << "=     where:                                                                                      ="<<std::endl;
    std::cout << "     " << "=                 -i INFILE        --> Input Mesh file (usually in OFF format)                    ="<<std::endl;
    std::cout << "     " << "=                 -o OUTFILE       --> Ouput Mesh file in OFF format                              ="<<std::endl;
    std::cout << "     " << "=                 -f FLOW          --> 0: MCF   1:cMCF   2:cMCF w/ projection on sphere           ="<<std::endl; 
    std::cout << "     " << "=                                      3: Gauss Map initializer   4: Tutte Embedding initializer  ="<<std::endl;
    std::cout << "     " << "=                 -w Tutte weights --> 1: cotangent weights    2: graph lapacian                  ="<<std::endl; 
    std::cout << "     " << "=                 -max             --> Max. no. of steps (integer), default: 2^6 = 64             ="<<std::endl;
    std::cout << "     " << "=                 -stps            --> Step size, default: 0.01				                           ="<<std::endl;
    std::cout << "     " << "=                 -tol             --> Tolerance for sphericity, default: 0.001                   ="<<std::endl;
    std::cout << "     " << "=                                                                                                 ="<<std::endl;
    std::cout << "     " << "==================================================================================================="<<std::endl;
    std::cout << "     " << "==================================================================================================="<<std::endl;
    std::cout << "\n\n" <<std::endl;
}

/* ===============================================================================================
   Parse Argument from command line:
   =============================================================================================== */

bool parse_args(int argc, char* argv[], std::string* INfile, std::string* OUTfile, std::string* STRflow, std::string* STRweight, std::string* STRmax, std::string* STRstps, std::string* STRtol)

{
//
// Make sure we have at least two parameters....
//
	std::string param;
	if (argc == 1) 
	{
		return false;
	} 
	else 
	{
		for (int i = 1; i < argc - 1; i = i + 2) 
		{
			param = argv[i];

			if (param == "-i") {
				*INfile = argv[i + 1];
			}
			else if (param == "-o") {
				*OUTfile = argv[i + 1];
			}
			else if (param == "-f") {
				*STRflow = argv[i + 1];
			}
      else if (param == "-w") {
        *STRweight = argv[i + 1];
      }
			else if (param == "-max") {
				*STRmax = argv[i + 1];
			}
			else if (param == "-stps") {
				*STRstps = argv[i + 1];
			}
			else if (param == "-tol") {
				*STRtol = argv[i + 1];
			}
		}
  	}
	return true;
}

#endif
