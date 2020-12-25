/* ==============================================================================================
   Tutte embedding
  
   Author:  Karry
   Date:    04/19/2020

   Input: arbitrary mesh
   Output: mesh flattened in a hexagon
   =============================================================================================== */

#ifndef _tutte_
#define _tutte_

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include <iostream>
#include <math.h>
#include <cmath>

/* ===============================================================================================
   Includes associated to Eigen or CHOLMOD
   =============================================================================================== */

#include <Eigen/Core>
#include <Eigen/Sparse>

/* ===============================================================================================
   Includes associated to OpenMesh
   =============================================================================================== */

#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/IO/Options.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"

//  ===============================================================================================
//    Define local acronyms
//    =============================================================================================== 

typedef Eigen::Triplet<double> T;

/* ===============================================================================================
	Function for inverse stereographic projection
   =============================================================================================== */

template <typename MyMesh>

void inv_stereo(MyMesh &mesh)
{

	typename MyMesh::VertexIter			v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
	typename MyMesh::Point 				p; 
	int n_vertices = mesh.n_vertices();
	int n_edges    = mesh.n_edges();
	int n_faces = mesh.n_faces();
	int idx;

	// double		x0[n_vertices], y0[n_vertices], z0[n_vertices];
	double 		temp = 0.0;

	for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		p = mesh.point(*v_it);
		idx = (*v_it).idx();

		temp = 1 + p[0]* p[0] + p[1]* p[1];
		p[0] = 2* p[0] / temp; 
		p[1] = 2* p[1] / temp;
		p[2] = -2 / temp + 1;

		mesh.set_point(*v_it, p);
	}

}


/* ===============================================================================================
	Driver function for Tutte embedding
   =============================================================================================== */

template <typename MyMesh>

void tuttemap(MyMesh &mesh, int &iweight)
{
	typename MyMesh::FaceIter			f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
	typename MyMesh::VertexIter			v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());

	typename MyMesh::FaceHalfedgeIter	fh_it;
	typename MyMesh::VertexHandle		v0, v1, v2;
	typename MyMesh::Point				p0, p1, p2;
	typename MyMesh::HalfedgeHandle 	h, h1;

	int n_vertices = mesh.n_vertices();
	int n_edges    = mesh.n_edges();
	int n_faces = mesh.n_faces();
	int r0, r1, r2; 

	double a0 = 0.0, a1 = 0.0, a2 = 0.0;
	double l0, l1, l2;

	// // initialization - compute discrete Gaussian curvature
	// double curvature[n_vertices];
	// for (int i = 0; i < n_vertices; i++) {curvature[i] = 2* M_PI;}

	// for (f_it = f_begin; f_it != f_end; ++f_it)
	// {
	// 	fh_it = mesh.fh_iter(*f_it);	 
	//     	typename MyMesh::HalfedgeHandle h0(*fh_it);
	// 	v0 = mesh.to_vertex_handle(h0);
	// 	r0 = v0.idx();
	// 	v1 = mesh.from_vertex_handle(h0);
	// 	r1 = v1.idx();
	// 	p0 = mesh.point(v0);
	// 	p1 = mesh.point(v1);

	// 	h1 = mesh.next_halfedge_handle(h0);
	// 	v2 = mesh.to_vertex_handle(h1); 
	// 	r2 = v2.idx();
	// 	p2 = mesh.point(v2);

	// 	l0 = (p1 - p2).length();
	// 	l1 = (p2 - p0).length();
	// 	l2 = (p0 - p1).length();

	// 	a0 = acos((l1*l1 + l2*l2 - l0*l0)/ (2.0*l1*l2));
	// 	a1 = acos((l0*l0 + l2*l2 - l1*l1)/ (2.0*l0*l2));
	// 	a2 = acos((l0*l0 + l1*l1 - l2*l2)/ (2.0*l0*l1));

	// 	curvature[r0] -= a0;
	// 	curvature[r1] -= a1;
	// 	curvature[r2] -= a2; 
	// }

	// // double* k = std::max_element(curvature, curvature + n_vertices); // maximal value in discrete curvature
	// int ind = std::max_element(curvature, curvature + n_vertices)  - curvature; // index of maximal value in discrete curvature
	int ind = n_vertices - 1;

	// the request has to be called before a vertex/face/edge can be deleted. it grants access to the status attribute
	mesh.request_face_status();
  	mesh.request_edge_status();
	mesh.request_vertex_status();

	//Delete vertex with max discrete curvature
	int valence = 0;
	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
        // if((*v_it).idx() == ind) {mesh.set_color(*v_it, typename MyMesh::Color(255,0,0));}
        if((*v_it).idx() == ind) {
        	std::cout << "Vertex to be deleted: " << *v_it << std::endl; 
        	valence = mesh.valence(*v_it);	
        	std::cout << "Its valence is " << valence << std::endl;
        	mesh.delete_vertex(*v_it, false);
        }
	}
	mesh.garbage_collection();

	//Compute boundary coordinates
	double boundary_x[valence];
	double boundary_y[valence];
	int j;

	for (j = 0; j != valence; ++ j)
	{
		boundary_x[j] = cos( (2* j* M_PI) / valence );
		boundary_y[j] = sin( (2* j* M_PI) / valence );
	}

	// Update
	if(mesh.is_trimesh()){std::cout << "Still a TriMesh!" << std::endl;}

	//Assemble sparse linear system (I - W)x  = b
	n_vertices = mesh.n_vertices();
	n_edges    = mesh.n_edges();
	n_faces    = mesh.n_faces();
	std::cout << "# Vertices: " << n_vertices << std::endl;
	std::cout << "# Edges   : " << n_edges << std::endl;
	std::cout << "# Faces   : " << n_faces << std::endl;

	std::vector<T> tripletList;
	tripletList.reserve( 2*n_edges + n_vertices );  // nonzero elements
	Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
	Eigen::VectorXd				bx(n_vertices), by(n_vertices); 
	bx.setZero();
	by.setZero();
	Eigen::VectorXd 			x0(n_vertices), y0(n_vertices), z0(n_vertices);
	z0.setZero();

	//Assemble matrix A = (I - W)
	typename MyMesh::EdgeIter			e_it, e_begin(mesh.edges_begin()), e_end(mesh.edges_end());
	f_begin = mesh.faces_begin(), f_end = mesh.faces_end();
	v_begin = mesh.vertices_begin(), v_end = mesh.vertices_end();

	typename MyMesh::VertexHandle		vf, vt;
	int idf, idt; //Matrix index
	double cos0, cos1, cos2, cot0, cot1, cot2;

	//initialization
	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
		tripletList.push_back(T( (*v_it).idx(), (*v_it).idx(), 0));  //Diagonal;
	}

	for (e_it = e_begin; e_it != e_end; ++e_it)
	{
		h = mesh.halfedge_handle(*e_it, 0);
		vf = mesh.from_vertex_handle(h);
		vt = mesh.to_vertex_handle(h);
		idf = vf.idx();
		idt = vt.idx();

		if ( !mesh.is_boundary(*e_it) ){
		tripletList.push_back(T(idf, idt, 0));  //Off-diagonal;
		tripletList.push_back(T(idt, idf, 0));
		}
	}	
	A.setFromTriplets(tripletList.begin(),tripletList.end());

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

		if (iweight == 2){
			int valence0, valence1, valence2;
			valence0 = mesh.valence(v0);	
			valence1 = mesh.valence(v1);	
			valence2 = mesh.valence(v1);

			if( !mesh.is_boundary(v0) ) {
				A.coeffRef(r0,r1) -= 1.0 / valence0;
				A.coeffRef(r0,r2) -= 1.0 / valence0;
			}

			if( !mesh.is_boundary(v1) ) {
				A.coeffRef(r1,r0) -= 1.0 / valence1;
				A.coeffRef(r1,r2) -= 1.0 / valence1;
			}		

			if( !mesh.is_boundary(v2) ) {
				A.coeffRef(r2,r0) -= 1.0 / valence2;
				A.coeffRef(r2,r1) -= 1.0 / valence2;
			}
		}

		else{
			l0 = (p1 - p2).length();
			l1 = (p2 - p0).length();
			l2 = (p0 - p1).length();

			//Cotangent angles
			cos0 = (l1*l1 + l2*l2 - l0*l0)/ (2.0*l1*l2);
			cot0 = cos0/sqrt(1 - cos0 * cos0);

			cos1 = (l0*l0 + l2*l2 - l1*l1)/ (2.0*l0*l2);
			cot1 = cos1/sqrt(1 - cos1 * cos1);

			cos2 = (l0*l0 + l1*l1 - l2*l2)/ (2.0*l0*l1);
			cot2 = (1 - cot0 * cot1)/ (cot0 + cot1); // cotangent - sum of angles

			if( !mesh.is_boundary(v0) ) {
				A.coeffRef(r0,r1) -= cot2;
				A.coeffRef(r0,r2) -= cot1;
			}

			if( !mesh.is_boundary(v1) ) {
				A.coeffRef(r1,r0) -= cot2;
				A.coeffRef(r1,r2) -= cot0;
			}		

			if( !mesh.is_boundary(v2) ) {
				A.coeffRef(r2,r0) -= cot1;
				A.coeffRef(r2,r1) -= cot0;
			}
		}		
	}

	//Normalization of each row for interior vertices
	typename MyMesh::VertexVertexIter	vv_it;
	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
		r0 = (*v_it).idx();
		if( !mesh.is_boundary(*v_it) )
		{
			for (vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
			{
				r1 = (*vv_it).idx();
				A.coeffRef(r0, r0) -= A.coeffRef(r0, r1);
			}
		}
	}

	//Boundary vertices
	typename MyMesh::HalfedgeIter  	 	h_it, h_begin(mesh.halfedges_begin()), h_end(mesh.halfedges_end());
	for (h_it = h_begin; h_it != h_end; ++h_it)
	{
		if( mesh.is_boundary(*h_it) ){
			typename MyMesh::HalfedgeHandle h0(*h_it); 
			for (int j = 0; j < valence; ++j)
			{
				v0 = mesh.from_vertex_handle(h0);
				r0 = v0.idx();
				A.coeffRef(r0, r0) += 1;
				bx[v0.idx()] = boundary_x[j];
				by[v0.idx()] = boundary_y[j];
				h0 = mesh.next_halfedge_handle(h0);
			}
			break;
		}
	}

	//std::cout << "Sparse matrix is " << A << std::endl;
	// std::cout << "Vector bx is \n" << bx << std::endl;
	// std::cout << "Vector by is \n" << by << std::endl;

	//Set up Linear system Ax = b
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//Decompose matrix A
	solver.compute(A);
	x0 = solver.solve(bx);
	y0 = solver.solve(by);

	typename MyMesh::Point       p;

	for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		p[0] = x0[(*v_it).idx()]; 
		p[1] = y0[(*v_it).idx()]; 
		p[2] = z0[(*v_it).idx()];
		mesh.set_point(*v_it, p);
	}

	// std::cout << "Vector x0 is \n" << x0 << std::endl;
	// std::cout << "Vector y0 is \n" << y0 << std::endl;

	//Stereographic projection
	inv_stereo(mesh);

	//Add back north pole vertex (0, 0, 1)
	std::vector<typename MyMesh::VertexHandle>	verh;
	typename	MyMesh::Point				np;
	typename	MyMesh::VertexHandle		v[3];

	np[0] = 0.0; np[1] = 0.0; np[2] = 1.0;

	v[2] = mesh.add_vertex(np);

	for (h_it = h_begin; h_it != h_end; ++h_it)
	{
		if ( mesh.is_boundary(*h_it) ){
			typename MyMesh::HalfedgeHandle h0(*h_it); 
				// std::cout << "Vertices are, v0: " << v[0] << ", v1: " << v[1] << ", v2: " << v[2] << std::endl;

				v[0] = mesh.from_vertex_handle(h0);
				v[1] = mesh.to_vertex_handle(h0);

				verh.clear();
				verh.push_back(v[0]); 
				verh.push_back(v[1]); 
				verh.push_back(v[2]);
				mesh.add_face(verh);

				verh.clear();
				verh.push_back(v[1]); 
				verh.push_back(v[0]); 
				verh.push_back(v[2]);
				mesh.add_face(verh);
		}
	}

/*
	if(mesh.is_trimesh()){std::cout << "Still a TriMesh!" << std::endl;}

	n_vertices = mesh.n_vertices();
	n_edges    = mesh.n_edges();
	n_faces    = mesh.n_faces();
	std::cout << "# Vertices: " << n_vertices << std::endl;
	std::cout << "# Edges   : " << n_edges << std::endl;
	std::cout << "# Faces   : " << n_faces << std::endl;
*/
}

#endif