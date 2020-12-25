/* ===============================================================================================
   CheckMesh.cpp

   Small functions that read in the number of vertices and faces from an input file.
   Those functions are designed for the following format:
        - OFF
        - OBJ
        - PLY
   Only the ASCII version of those formats is considered
  
   Author:  Karry Wong
   Date:    04/4/2020
   Version: 2
 
   =============================================================================================== */

/* ===============================================================================================
   Includes
   =============================================================================================== */

#include "CheckMesh.h"

/* ===============================================================================================
   Get number of vertices and number of faces from an input file (either OFF, OBJ, or PLY format)
   =============================================================================================== */

bool mesh_info(std::string INfile, int *nvertex, int *nface)
{
	if(INfile.substr(INfile.find_last_of(".") + 1) == "off") 
	{
		mesh_info_off(INfile, nvertex, nface);
		return true;
	}
	else if(INfile.substr(INfile.find_last_of(".") + 1) == "obj") 
	{
		mesh_info_obj(INfile, nvertex, nface);
		return true;
	}
	else if(INfile.substr(INfile.find_last_of(".") + 1) == "ply") 
	{
		mesh_info_ply(INfile, nvertex, nface);
		return true;
	}
	else
	{
		return false;
	}
		
}

/* ===============================================================================================
   Get number of vertices and number of faces from an OFF file
   =============================================================================================== */

void mesh_info_off(std::string INfile, int *nvertex, int *nface)
{
	std::ifstream infile(INfile.c_str());
	std::string line;

	std::getline(infile,line);
	std::getline(infile,line);
	std::istringstream iss(line);
	iss >> *nvertex >> *nface;

	infile.close();
}

/* ===============================================================================================
   Get number of vertices and number of faces from an OBJ file
   =============================================================================================== */

void mesh_info_obj(std::string INfile, int *nvertex, int *nface)
{
	std::ifstream infile(INfile.c_str());
	std::string line;
	std::string keyword; 

	int nv=0;
	int nf=0;
	while (infile && ! infile.eof())
	{
		std::getline(infile,line);
		std::istringstream iss(line);
		iss >> keyword;
		if(keyword == "v") nv++;
		if(keyword == "f") nf++;
	}

	*nvertex = nv;
	*nface   = nf;

	infile.close();
}
/* ===============================================================================================
   Get number of vertices and number of faces from a PLY file
   =============================================================================================== */

void mesh_info_ply(std::string INfile, int *nvertex, int *nface)
{
	std::ifstream infile(INfile.c_str());
	std::string line;
	std::string keyword=" "; 
	std::string elementName; 

	int nv,nf;
	while (keyword != "end_header")
	{
		std::getline(infile,line);
		std::istringstream iss(line);
		iss >> keyword;
		if(keyword == "element")
		{
			iss >> elementName;
			if(elementName == "vertex")
			{
				iss >>nv;
			}
			else if(elementName == "face")
			{
				iss >>nf;
			}
		}
	}

	*nvertex = nv;
	*nface   = nf;

	infile.close();
}
