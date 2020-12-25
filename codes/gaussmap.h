/* ==============================================================================================
   Gauss Map
  
   Author:  Karry
   Date:    02/29/2020

   Input: arbitrary mesh
   Output: mesh with its vertices on unit sphere
   =============================================================================================== */

#ifndef _gaussmap_
#define _gaussmap_

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include <iostream>
#include <math.h>
#include <cmath>

/* ===============================================================================================
   Includes associated to OpenMesh
   =============================================================================================== */

#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/IO/Options.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"

/* ===============================================================================================
	Driver function for Gauss Map
   =============================================================================================== */
template <typename MyMesh>

void gaussmap(MyMesh &mesh)

{
	typename MyMesh::VertexIter			v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
  typename MyMesh::FaceIter       f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());

	// gaussmap(mesh);
  	mesh.request_vertex_normals();

	// we need face normals to update the vertex normals
  	mesh.request_face_normals();

	// let the mesh update the normals
  	mesh.update_normals();

	// dispose the face normals, as we don't need them anymore
  	mesh.release_face_normals();

    // assure we have vertex normals
    if (!mesh.has_vertex_normals())
    {
    	std::cerr << "ERROR: Standard vertex property 'Normals' not available!\n";
 	  }


	for (v_it = v_begin; v_it != v_end; ++v_it)
  	{
    // std::cout << "Vertex #" << *v_it << ": " << mesh.point( *v_it )  << std::endl;
    // std::cout << "Normal: " << mesh.normal(*v_it) << std::endl;
  	mesh.set_point( *v_it, mesh.normal(*v_it) );
    // std::cout << " moved to " << mesh.point( *v_it ) << std::endl;
  	}

  	// Remove normals  
  	mesh.release_vertex_normals();
}

#endif
