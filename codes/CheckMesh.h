/* ===============================================================================================
   CheckMesh.h: header file for CheckMesh
  
   Author:  Karry Wong
   Date:    04/04/2020
   Version: 2
   =============================================================================================== */

#ifndef _CHECKMESH_
#define _CHECKMESH_

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <fstream>

/* ===============================================================================================
   All prototypes
   =============================================================================================== */

bool mesh_info(std::string INfile, int *nvertex, int *nface);
void mesh_info_off(std::string INfile, int *nvertex, int *nface);
void mesh_info_obj(std::string INfile, int *nvertex, int *nface);
void mesh_info_ply(std::string INfile, int *nvertex, int *nface);
void mesh_info_stl(std::string INfile, int *nvertex, int *nface);

#endif
