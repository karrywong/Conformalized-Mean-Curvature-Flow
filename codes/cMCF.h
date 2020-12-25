/* ==============================================================================================
   Confomalized Mean Curvature Flow
  
   Author:  Karry
   Date:    04/04/2020

   Reference:	Conformalized Mean Curvature Flow proposed by Kazhdan, Solomon, and Ben-Chen:
   M. Kazhdan, J. Solomon, and M. Ben-Chen. "Can Mean Curvature Flow be Modified to be Non-singular?",
   Eurographics Symposium on Geometry Processing 2012. 
   =============================================================================================== */

#ifndef _cMCF_
#define _cMCF_

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>

/* ===============================================================================================
   Includes associated to Eigen or CHOLMOD
   =============================================================================================== */

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <Eigen/CholmodSupport>

/* ===============================================================================================
   Includes associated to OpenMesh
   =============================================================================================== */

#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/IO/Options.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "OpenMesh/Core/Geometry/VectorT.hh"
#include "OpenMesh/Core/Mesh/PolyConnectivity.hh"
#include "OpenMesh/Core/Mesh/TriConnectivity.hh"

/* ===============================================================================================
   Define local acronyms
   =============================================================================================== */

typedef Eigen::Triplet<double> T;

/* =============================================================================================================
   This function assembles mass matrix "mass" and stiffness matrix "stiff".
	Input: mesh		Output: mass and stiff
   ============================================================================================================= */

template <typename MyMesh>
void initial(int& iflow, MyMesh& mesh, std::vector<Eigen::Triplet<double>>& tripletList, Eigen::SparseMatrix<double>& mass, Eigen::SparseMatrix<double>& stiff, int niter, double* angdt, double ang0[], double ang1[], double ang2[])
{
	double tempang = 0.0;
	int n_faces = mesh.n_faces();
	int id, idf, idt; //Matrix index

	typename MyMesh::FaceIter			f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
	typename MyMesh::VertexIter			v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
	typename MyMesh::EdgeIter			e_it, e_begin(mesh.edges_begin()), e_end(mesh.edges_end());
	typename MyMesh::FaceHalfedgeIter   fh_it;
	typename MyMesh::HalfedgeHandle 	h, h1;
	typename MyMesh::VertexHandle		v, vf, vt;

	typename MyMesh::VertexHandle		v0, v1, v2;
	typename MyMesh::Point				p0, p1, p2;
	int r0, r1, r2; //Matrix index
	double l0, l1, l2, cos0, cos1, cos2;

	int j;

/* =============================================================================================================
   Loop over all triangles in the mesh
   ============================================================================================================= */

	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
		typename MyMesh::VertexHandle	v(*v_it);
		id = v.idx();	
		tripletList.push_back(T(id, id, 0.0));  //Diagonal;
	}

	for (e_it = e_begin; e_it != e_end; ++e_it)
	{
		h = mesh.halfedge_handle(*e_it, 0.0);
		vf = mesh.from_vertex_handle(h);
		vt = mesh.to_vertex_handle(h);
		idf = vf.idx();
		idt = vt.idx();
		tripletList.push_back(T(idf, idt, 0.0));  //Off-diagonal;
		tripletList.push_back(T(idt, idf, 0.0));
	}	
	mass.setFromTriplets(tripletList.begin(),tripletList.end());
	stiff = mass;

	makeDL(iflow, mesh, mass, stiff, niter, &tempang, ang0, ang1, ang2);	
	*angdt = tempang;

	for (f_it = f_begin, j = 0; f_it != f_end; ++f_it, ++j)
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

		cos0 = (l1*l1 + l2*l2 - l0*l0)/ (2.0*l1*l2);
		cos1 = (l0*l0 + l2*l2 - l1*l1)/ (2.0*l0*l2);
		cos2 = (l0*l0 + l1*l1 - l2*l2)/ (2.0*l0*l1);

		ang0[j] = acos(cos0);
		ang1[j] = acos(cos1);
		ang2[j] = acos(cos2);
	}
}


/* =============================================================================================================
   This function updates the mass matrix "mass" and, if needed, stiffness matrix
	Input: mesh		Output: mass and stiff
   ============================================================================================================= */

template <typename MyMesh>
void makeDL(int& iflow, MyMesh& mesh, Eigen::SparseMatrix<double>& mass, Eigen::SparseMatrix<double>& stiff, int& niter, double* angdt, double ang0[], double ang1[], double ang2[])
{
	int n_vertices = mesh.n_vertices();
	int n_faces = mesh.n_faces();
	int r0, r1, r2; //Matrix index
	double l0, l1, l2, cos0, cos1, cos2, cot0, cot1, cot2;
	double val, a_rs;	

	double a0 = 0.0, a1 = 0.0, a2 = 0.0;
	*angdt = 0.0;

	typename MyMesh::FaceIter			f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
	typename MyMesh::FaceHalfedgeIter	fh_it;
	typename MyMesh::HalfedgeHandle		h1;
	typename MyMesh::VertexHandle		v0, v1, v2;
	typename MyMesh::Point				p0, p1, p2;

	typename MyMesh::EdgeHandle			eh;
	typename MyMesh::FaceHandle			opf;

	// typename std::vector<typename MyMesh::FaceHandle>				bad;
	// typename std::vector<typename MyMesh::FaceHandle>::iterator		bad_it;
	// bad.reserve(0.0001 * n_faces);

	typename std::vector<typename MyMesh::FaceHandle>::iterator		position;

	int i = 0, j = 0;
	int count_nan = 0, count_inf = 0;

	for (f_it = f_begin, j = 0; f_it != f_end; ++f_it, ++j)
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

		val = (l2+l1+l0) * (l0-l2+l1) * (l0+l2-l1) * (l2+l1-l0);

		//Area - Heron's formula and volume
		if (val > 0){a_rs = 0.25 * sqrt(val);}
		else {a_rs = 0;}

		//Assemble mass //a_rs = a_rs/12;
		mass.coeffRef(r0,r1) += a_rs; //Off-diagonal
		mass.coeffRef(r1,r0) += a_rs;
		mass.coeffRef(r0,r2) += a_rs;
		mass.coeffRef(r2,r0) += a_rs;
		mass.coeffRef(r1,r2) += a_rs;
		mass.coeffRef(r2,r1) += a_rs;

		a_rs = 2.0* a_rs;
		mass.coeffRef(r0,r0) += a_rs; //Diagonal
		mass.coeffRef(r1,r1) += a_rs;    
		mass.coeffRef(r2,r2) += a_rs; 

		a0 = acos((l1*l1 + l2*l2 - l0*l0)/ (2.0*l1*l2));
		a1 = acos((l0*l0 + l2*l2 - l1*l1)/ (2.0*l0*l2));
		a2 = acos((l0*l0 + l1*l1 - l2*l2)/ (2.0*l0*l1));

		if(niter != 0){
			if 		( std::isnan(a0) || std::isnan(a1) || std::isnan(a2) )	{++count_nan;} 
			else if ( std::isinf(a0) || std::isinf(a1) || std::isinf(a2) )	{++count_inf;} 
			else
				{*angdt += (abs(a0 - ang0[j]) + abs(a1 - ang1[j]) + abs(a2 - ang2[j]));}
		}
		else {*angdt = 0.0;}		

		if(iflow == 0 || (iflow != 0 && niter == 0)){
		//Cotangent angles
		cos0 = (l1*l1 + l2*l2 - l0*l0)/ (2.0*l1*l2);
		cot0 = cos0/sqrt(1 - cos0 * cos0);

		cos1 = (l0*l0 + l2*l2 - l1*l1)/ (2.0*l0*l2);
		cot1 = cos1/sqrt(1 - cos1 * cos1);

		cos2 = (l0*l0 + l1*l1 - l2*l2)/ (2.0*l0*l1);
		cot2 = (1 - cot0 * cot1)/ (cot0 + cot1); // cotangent - sum of angles

		cot0 = 0.5*cot0; cot1 = 0.5*cot1; cot2 = 0.5*cot2;

		//Assemble stiff
		stiff.coeffRef(r0,r1) += cot2; //Off-diagonal
		stiff.coeffRef(r1,r0) += cot2;
		stiff.coeffRef(r0,r2) += cot1;
		stiff.coeffRef(r2,r0) += cot1;
		stiff.coeffRef(r1,r2) += cot0;
		stiff.coeffRef(r2,r1) += cot0;	

		stiff.coeffRef(r0,r0) -= (cot1 + cot2);  //Diagonal;
		stiff.coeffRef(r1,r1) -= (cot0 + cot2);
		stiff.coeffRef(r2,r2) -= (cot0 + cot1);
		}
	}
	*angdt /= (3* n_faces);

	std::cout << "Number of NaN in angdt is: " << count_nan << ", number of Inf in angdt is: " << count_inf << std::endl; 
	/*for (bad_it = bad.begin(); bad_it != bad.end(); ++bad_it){
		fh_it = mesh.fh_iter(*bad_it);
		eh = mesh.edge_handle(*fh_it);

		opf = mesh.opposite_face_handle(*fh_it);
		position = find(bad.begin(), bad.end(), opf);
		if (position != bad.end())	{bad.erase(position);}

		if (mesh.is_flip_ok(eh)){mesh.flip(eh);}
	}

	mesh.garbage_collection();
	bad.clear();
	*/
}


/* =============================================================================================================
   This function computes the center of mass
	Input: mesh		Output: center of mass
   ============================================================================================================= */
template <typename MyMesh, typename Point>
void center(MyMesh& mesh, Point& com)
{
	typename MyMesh::VertexIter  v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
	int n_vertices = mesh.n_vertices();

	for (v_it = v_begin; v_it != v_end; ++v_it){com += mesh.point(*v_it);}
	com = com / n_vertices; //center of mass
}

/* =============================================================================================================
   This function positions the mesh such that its center of mass is at origin
	Input: mesh
   ============================================================================================================= */
template <typename MyMesh>
void centermesh(MyMesh& mesh)
{
	typename MyMesh::VertexIter  v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
	typename MyMesh::Point       p, com;

	com[0] = com[1] = com [2] = 0.0; 
	center(mesh, com); //compute COM

	//Centering mesh
	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
		p = mesh.point(*v_it);
		mesh.set_point(*v_it, p - com);
	}
}

/* =============================================================================================================
   This function perform the sphere projection
	Input: mesh, com		Output: mesh
   ============================================================================================================= */
template <typename MyMesh>
void spherproj(MyMesh& mesh)
{
	typename MyMesh::VertexIter  v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
	typename MyMesh::Point	     com, pt;
	com[0] = com[1] = com [2] = 0.0; 	
	double r = 0.0, len, sss;

	int n_vertices = mesh.n_vertices();

	center(mesh, com);
	for (v_it = v_begin; v_it != v_end; ++v_it){
		r += (mesh.point(*v_it) - com).length();}
		r = r / n_vertices; //average radius

		//std::cout << "Projection onto sphere begins with radius: " << r << " and center of mass: " << com << std::endl; 
	for (v_it = v_begin; v_it != v_end; ++v_it){
		pt  = mesh.point(*v_it);
		len = (pt - com).length();
		sss = r / len;
		pt  = (pt - com) * sss + com;
		mesh.set_point(*v_it, pt);}
	std::cout << "Projection done" << std::endl;
}

/* =============================================================================================================
   This function computes the variance of radius
	Input: mesh		Output: center of mass, radius variance, updated com
   ============================================================================================================= */
template <typename MyMesh>
void rvar(MyMesh& mesh, double* var, double* spher2)
{
	typename MyMesh::VertexIter  v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
	typename MyMesh::Point p;

	*var = 0.0, *spher2 = 0.0;
	int n_vertices = mesh.n_vertices();
	double r = 0.0;

	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
		r += (mesh.point(*v_it)).length();
	}
		r = r / n_vertices; //average radius

	double temp, temp1, temp2;	
	for (v_it = v_begin; v_it != v_end; ++v_it)
	{ temp = (mesh.point(*v_it)).length();
	  temp1 = (temp - r) * (temp - r);
	  *var += temp1;

	  temp2 = abs(temp - r);
	  *spher2 += temp2;
	}
	*var = *var / n_vertices;
	*spher2 = *spher2 / n_vertices;
}

/* ===============================================================================================
	Driver function for Conformalized mean curvature flow 
   =============================================================================================== */
template <typename MyMesh>

void cMCF(int& iflow, int& iweight, MyMesh &mesh, std::vector<Eigen::Triplet<double> > tripletList, Eigen::SparseMatrix<double> D, Eigen::SparseMatrix<double> L, Eigen::SparseMatrix<double> A, double step_size, int Nitermax, double TOL, double spherTOL, std::string INfile)

{
	int n_vertices = mesh.n_vertices();
	int n_faces = mesh.n_faces();
	int n_edges = mesh.n_edges();
	typename MyMesh::VertexIter  v_it;
    typename MyMesh::FaceIter   f_it;
	typename MyMesh::Point p; //For cMCF w/ projection onto sphere
	bool proj = false;

	std::ofstream fs;
	std::string s = INfile;
	s.erase(s.find_last_of("."), std::string::npos);
	std::string filename;

	if (iflow != 4){
		filename = s + "_measure" + ".csv";
	}
	else{
		filename = s + "_measure_cmcf4" + ".csv";
	}
	double record[2*Nitermax][7];


/* ===============================================================================================
	Initialize the two sparse matrices needed by the flow: mass matrix D and stiffness matrix L
	And get initial sphericity of the mesh
   =============================================================================================== */

	int niter = 0;

	double ang0[n_faces], ang1[n_faces], ang2[n_faces];
	double area, volume, var, spher, spher2, areanew, volumenew, ss;
	double angdt = 0.0, angdtnew = 0.0;
	// double mintri0, mintri, check;

	double clrmean = 0.0;
    double clrbegin[n_edges], clrnew[n_edges];

    centermesh(mesh);
	rvar(mesh, &var, &spher2); //Compute radius variance

	typename MyMesh::Point		com;
	com[0] = com[1] = com[2] = 0.0;
	center(mesh, com);

    //Initialize D, L
	initial(iflow, mesh, tripletList, D, L, niter, &angdt, ang0, ang1, ang2);
	measure(mesh, &area, &volume, &spher);
	measure2(mesh, clrbegin);

	// std::cout << "D: " << D << std::endl;
	// std::cout << "L: " << L << std::endl;

	//mt(mesh, &mintri0);
	//measure_ang(mesh, INfile, niter);

	//std::cout << "Initial: " << niter << ", area: " << area << ", volume: " << volume <<  ", sphericity: " << spher <<  ", sphericity 2: " << spher2 << ", radius variance: " << var << ", COM: " << com << std::endl;


	if(iflow == 3 || iflow == 4)
	{
		if (iflow == 3) {
			gaussmap(mesh);
			//centermesh(mesh); 
			//CenterSphere(mesh);
		}
		if (iflow == 4) {
			std::cout << "Checking, iweight = " << iweight << std::endl;
			tuttemap(mesh, iweight); 
			CenterSphere(mesh);
		}

		scale(mesh);

		measure(mesh, &area, &volume, &spher);
		rvar(mesh, &var, &spher2); //Compute radius variance
		com[0] = com[1] = com[2] = 0.0;
		center(mesh, com); //Compute COM

		niter = 1;
		makeDL(iflow, mesh, D, L, niter, &angdtnew, ang0, ang1, ang2);
		niter = 0;	
	}


	measure2(mesh, clrnew);
	//for (t = 0; t != n_edges; t++) {std::cout << "t: " << t << ", clrbegin[t]: " << clrbegin[t] << std::endl;}

	int t;
	int count_nan = 0, count_inf = 0;
	for (t = 0; t != n_edges; t++) {
		if ( std::isnan(clrnew[t]) ) 		{++count_nan;}
		else if ( std::isinf(clrnew[t]) )	{++count_inf;}
		else								{clrmean +=  clrnew[t] / clrbegin[t];}
	}
	clrmean /= n_edges;

	std::cout << "Number of NaN in clr is: " << count_nan << ", number of Inf in clr is: " << count_inf << std::endl; 
	std::cout << "# step: " << niter << ", area: " << area << ", volume: " << volume <<  ", sphericity: " << spher <<  ", sphericity 2: " << spher2 << ", radius variance: " << var << ", clrmean: " << clrmean << ", COM: " << com << std::endl;
	//std::cout << "Area of smallest triangle: " << mintri0 << std::endl;
	record[niter][0] = spher;
	record[niter][1] = var;
	record[niter][2] = volume;
	record[niter][3] = angdt;
	record[niter][4] = area;
	record[niter][5] = clrmean;
	record[niter][6] = spher2;
	//record[niter][5] = mintri0;

/* ===============================================================================================
	Initialize left hand sides for x, y, and z
   ==============================================================================================*/

	Eigen::VectorXd x0(n_vertices), y0(n_vertices), z0(n_vertices);
	Eigen::VectorXd b(n_vertices); 
	int idx;

	for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		p = mesh.point(*v_it);
		idx = (*v_it).idx();
		x0[idx] = p[0]; y0[idx] = p[1]; z0[idx] = p[2];
	}

/* ===============================================================================================
	Define Eigen solver to be used
   =============================================================================================== */

	//Built-in Eigen solvers
//	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

/* ===============================================================================================
    Compute gradient vector
   ==============================================================================================*/
 
	//grad(mesh, D, L, x0, y0, z0, ang0, ang1, ang2, area, s);

/* ===============================================================================================
	Apply (c)MCF flow to get mesh spherical
   =============================================================================================== */
	step_size = 12* step_size; //avoid division by 12 in mass matrix D
	//volumenew = volume;
	//double v0 = 1e-3; 
	//int c0 = 0, c1 = 0;

	//Eigen::setNbThreads(4);
	while( niter < Nitermax) //std::abs(1 - spher) > TOL &&
	{
		//if(spher > 1) break;
		// Scale stiffness matrix by step size
		if(iflow == 0 || (iflow != 0 && niter == 0)) {L = L*step_size;}
		//if(abs(volumenew) < v0) {c0 = 1;}
		//if(c0 == 1 && c1 == 0)	{L = L* 50; c1 = 1;}

		//Set up Linear system Ax = b, with A = (D - sL), b = Dx, s: step size
		A = D - L;

		//Decompose matrix A
		solver.compute(A);

		//Assemble vector b and solve Ax = b for each direction (x, y, and z)
		b = D * x0; 
		x0 = solver.solve(b);

		b = D * y0; 
		y0 = solver.solve(b);

		b = D * z0;
		z0 = solver.solve(b);

		//Transfer eigen arrays into mesh coordinates
		for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			idx = (*v_it).idx();
			p[0] = x0[idx]; p[1] = y0[idx]; p[2] = z0[idx];
			mesh.set_point(*v_it, p);
			//std::cout << "x0: " << x0[idx] << ", y0: " << y0[idx] << ", z0: " << z0[idx] << std::endl;
		}


		//Project onto sphere once the mesh becomes more spherical
		if(iflow == 2 && proj)
		{
			spherproj(mesh);
		}

		niter++;
		angdtnew = 0.0;

		//Reset mass matrix (and stiffness matrix if MCF chosen)
		D *= 0;
		if(iflow == 0) {L *= 0;}

    	centermesh(mesh);
		rvar(mesh, &var, &spher2); //Compute radius variance
		com[0] = com[1] = com[2] = 0.0;
		center(mesh, com);

		// //projection after each step
 	// 	spherproj(mesh);

		makeDL(iflow, mesh, D, L, niter, &angdtnew, ang0, ang1, ang2);
		measure(mesh, &areanew, &volumenew, &spher);

		//mt(mesh, &mintri);

		/*if( c0 == 1 )*/
		ss = sqrt(area/areanew);
		//ss = sqrt(mintri0/mintri);
		//D *= (mintri0/mintri);

		D *= ss* ss;
		volumenew = ss* ss* ss* volumenew;
		areanew = ss* ss* areanew;

		// std::cout << "D: " << D << std::endl;
		// std::cout << "L: " << L/12 << std::endl;

		//if( c0 == 1 ){spher = cbrt(M_PI* 36* volumenew * volumenew)/ area;} //asuming area kept constant
		//else{spher = cbrt(M_PI* 36* volumenew * volumenew)/ areanew;} //area NOT kept constant 

		if (spher > spherTOL) {proj = true;} //Activate projection onto sphere

		for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			p = mesh.point(*v_it);
			//if( c0 == 1) {p *= ss;}	//Scaling
			p *= ss;
			mesh.set_point(*v_it, p);
		}

		measure2(mesh, clrnew);

		count_nan = 0, count_inf = 0;
		for (t = 0; t != n_edges; t++) {
		if ( std::isnan(clrnew[t]) ) 		{++count_nan;}
		else if ( std::isinf(clrnew[t]) )	{++count_inf;}
		else								{clrmean +=  clrnew[t] / clrbegin[t];}
		}
		clrmean /= n_edges;

		std::cout << "Number of NaN in clr is: " << count_nan << ", number of Inf in clr is: " << count_inf << std::endl; 
		std::cout << "# step: " << niter << ", area: " << areanew << ", volume: " << volumenew <<  ", sphericity: " << spher <<  ", sphericity 2: " << spher2 << ", radius variance: " << var << ", clrmean: " << clrmean << ", COM: " << com << ", angular dist.: " << angdtnew << std::endl;
		//std::cout << "Area of smallest triangle: " << mintri << std::endl;
		//mt(mesh, &check);
		//std::cout << "Checking: " << check << std::endl;


		//Write sphericity and radius variance into CSV		
		record[niter][0] = spher;
		record[niter][1] = var;
		record[niter][2] = volumenew;
		record[niter][3] = angdtnew;
		record[niter][4] = areanew;
		record[niter][5] = clrmean;
		record[niter][6] = spher2;
		//record[niter][5] = mintri;

		//std::cout << "D: " << D << std::endl;
 		//std::cout << "L: " << L << std::endl;

	}

	// spherproj(mesh);

	// area = 0.0; volume = 0.0; spher = 0.0; angdt = 0.0;
	// gauge(mesh, &area, &volume, &spher, &angdt, ang0, ang1, ang2);
	// rvar(mesh, &var, &spher2); //Compute radius variance

	// measure2(mesh, clrnew);
	// clrmean = 0.0;
	// for (t = 0; t != n_edges; t++) {clrmean +=  clrnew[t] / clrbegin[t];}
	// clrmean /= n_edges;	

	//mt(mesh, &mintri);
	//if (niter != 0){measure_ang(mesh, INfile, niter);}

	std::cout << "Final area: " << area << ", final volume: " << volume <<  ", final sphericity: " << spher << ", sphericity 2: " << spher2 << ", final radius variance: " << var << ", final clrmean: " << clrmean  <<  ", final COM: " << com <<", final angular dist.: " << angdt << std::endl; 

	// record[niter+1][0] = spher;
	// record[niter+1][1] = var;
	// record[niter+1][2] = volume;
	// record[niter+1][3] = angdt;
	// record[niter+1][4] = area;
	// record[niter+1][5] = clrmean;
	// record[niter][6] = spher2;

	// Write final sphericity and radius variance
	// fs.open(filename);
	// for(int j = 0; j <= niter; j++)
	// {
	// 		fs << j << "," << record[j][0] << "," << record[j][1] << "," << record[j][2] << "," << record[j][3] <<  "," << record[j][4] <<  "," << record[j][5] << "," << record[j][6] << "\n";
	// }
	// fs.close();
	
	// std::string filename2 = s + "_clr" + ".csv";
	// fs.open(filename2);
	// for (t = 0; t != n_edges; t++) {
	// 	fs << clrbegin[t] << "," << clrnew[t] << "\n";
	// }
	// fs.close();
}


#endif
