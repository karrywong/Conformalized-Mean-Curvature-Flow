/* ==============================================================================================
   Scale mesh 
  
   Author:  Karry
   Date:    05/04/2020

   Input: mesh
   Output mesh with surface area equal to 4*pi
   =============================================================================================== */

#ifndef _scale_
#define _scale_

/* ===========================================================s====================================
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
	Driver function for scaling
   =============================================================================================== */

template <typename MyMesh>

void scale(MyMesh &mesh)
{

   typename MyMesh::FaceIter           f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
   typename MyMesh::VertexIter         v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
   typename MyMesh::FaceHalfedgeIter   fh_it;
   typename MyMesh::VertexHandle       v0, v1, v2;
   typename MyMesh::Point              p, p0, p1, p2;
   typename MyMesh::HalfedgeHandle     h1;

   double l0, l1, l2, area = 0.0, factor = 0.0;
   int r0, r1, r2; 

   for (f_it = f_begin; f_it != f_end; ++f_it)
   {
      fh_it = mesh.fh_iter(*f_it);   
         typename MyMesh::HalfedgeHandle h0(*fh_it);
      v0 = mesh.to_vertex_handle(h0);
      r0 = v0.idx();
      v1 = mesh.from_vertex_handle(h0);
      r1 = v1.idx();
      p0 = mesh.point(v0);
      p1 = mesh.point(v1);

      h1 = mesh.next_halfedge_handle(h0);
      v2 = mesh.to_vertex_handle(h1); 
      r2 = v2.idx();
      p2 = mesh.point(v2);

      l0 = (p1 - p2).length();
      l1 = (p2 - p0).length();
      l2 = (p0 - p1).length();

      area += 0.25 * sqrt( (l2+l1+l0) * (l0-l2+l1) * (l0+l2-l1) * (l2+l1-l0) );
   }

   factor = sqrt(4* M_PI / area);
   std::cout << "The area is " << area << " and the scaling factor is " << factor << std::endl;

   for (v_it = v_begin; v_it != v_end; ++v_it)
   {
         p = mesh.point(*v_it);
         p *= factor;
         mesh.set_point(*v_it, p);
   }

}

#endif