/* ===============================================================================================
   Geometry: Calculate area

   Author:  Karry Wong
   Date:    08/27/2018
   Version: 1
   =============================================================================================== */

#ifndef _GEOMETRY_
#define _GEOMETRY_

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>


/* ===============================================================================================
   Includes associated to OpenMesh
   =============================================================================================== */

#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "OpenMesh/Core/Geometry/VectorT.hh"
#include "OpenMesh/Core/Mesh/PolyConnectivity.hh"
#include "OpenMesh/Core/Mesh/TriConnectivity.hh"

/* =============================================================================================================
   This function calculates the sphericity, surface area and volume of the input mesh. Also, it calculates the minimum angle and the aspect ratio of each triangle.
	Input: mesh 		Output: area, volume, sphericity, (also possible: minimum angle, aspect ratio)
   ============================================================================================================= */

template <typename MyMesh>

void measure(MyMesh& mesh, double* area, double* volume, double* spher) //double ang[][3], double clr[]
{
	*area = 0.0; *volume = 0.0;
	int n_faces = mesh.n_faces();

	typename MyMesh::FaceIter	f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
	typename MyMesh::FaceHalfedgeIter  fh_it;
	typename MyMesh::FaceVertexIter	fv_it;
	typename MyMesh::HalfedgeHandle h0, h1;
	typename MyMesh::VertexHandle	v0, v1, v2, v11, v22, v33;
	typename MyMesh::Point		p0, p1, p2, p11, p22, p33, vec;
	typename MyMesh::Scalar		vecnorm;


	int j = 0;
	int r0, r1, r2; 
	double l0, l1, l2, p, val;

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

	    val =  (l2+l1+l0) * (l0-l2+l1) * (l0+l2-l1) * (l2+l1-l0);
		if (val <= 0){
			val = 0;
			++j;
		}
	   	*area += 0.25 * sqrt( val );
	    *volume += abs( dot(p0,cross(p1,p2)) / 6.0 ) ;

	    // ang[i][0] = acos((l1*l1 + l2*l2 - l0*l0)/ (2.0*l1*l2))*(180/M_PI);
	    // ang[i][1] = acos((l0*l0 + l2*l2 - l1*l1)/ (2.0*l0*l2))*(180/M_PI);
	    // ang[i][2] = acos((l0*l0 + l1*l1 - l2*l2)/ (2.0*l0*l1))*(180/M_PI);
/*	    p = 0.5*(l0 + l1 + l2); //aspect ratio (def. 0 =< 2r/R =< 1, r: inradius, R: circumradius)
	    AR[i] = (8*(p-l0)*(p-l1)*(p-l2))/(l0*l1*l2);
	   	++i;
*/
	}
	std::cout << "Total no. of degenerate triangles:  " << j << std::endl;

	*spher = cbrt(M_PI*36* *volume * *volume)/ *area; 	// Compute sphericity


	// typename MyMesh::EdgeIter	e_it, e_begin(mesh.edges_begin()), e_end(mesh.edges_end());
	// typename MyMesh::HalfedgeHandle hh0, hh0a, hh0b, hh1, hh1a, hh1b;
	// typename MyMesh::VertexHandle	vt0a, vf0a, vt0b, vf0b, vt1a, vf1a, vt1b, vf1b;
	// typename MyMesh::Point		pt0a, pf0a, pt0b, pf0b, pt1a, pf1a, pt1b, pf1b;
	// double len0a, len0b, len1a, len1b;
	// int id;

	// //Conformality - cross length ratio
	// for  (e_it = e_begin; e_it != e_end; ++e_it)
	// {
	//    hh0 = mesh.halfedge_handle(*e_it, 0);
	//    hh1 = mesh.halfedge_handle(*e_it, 1);

	//    hh0a = mesh.next_halfedge_handle(hh0);
	//    vt0a = mesh.to_vertex_handle(hh0a);
	//    vf0a = mesh.from_vertex_handle(hh0a);
	//    pt0a = mesh.point(vt0a);
	//    pf0a = mesh.point(vf0a);
	//    len0a = (pt0a - pf0a).length();

	//    hh0b = mesh.next_halfedge_handle(hh0a);
	//    vt0b = mesh.to_vertex_handle(hh0b);
	//    vf0b = mesh.from_vertex_handle(hh0b);
	//    pt0b = mesh.point(vt0b);
	//    pf0b = mesh.point(vf0b);
	//    len0b = (pt0b - pf0b).length();

	//    hh1a = mesh.next_halfedge_handle(hh1);
	//    vt1a = mesh.to_vertex_handle(hh1a);
	//    vf1a = mesh.from_vertex_handle(hh1a);
	//    pt1a = mesh.point(vt1a);
	//    pf1a = mesh.point(vf1a);
	//    len1a = (pt1a - pf1a).length();

	//    hh1b = mesh.next_halfedge_handle(hh1a);
	//    vt1b = mesh.to_vertex_handle(hh1b);
	//    vf1b = mesh.from_vertex_handle(hh1b);
	//    pt1b = mesh.point(vt1b);
	//    pf1b = mesh.point(vf1b);
	//    len1b = (pt1b - pf1b).length();

	//    id = (*e_it).idx();
	//    clr[id] = (len0a*len1a)/(len0b*len1b);
	// }
};

/* =============================================================================================================
   This function compute the cross-length ratio of the input mesh.
	Input: mesh 		Output: clr
   ============================================================================================================= */

template <typename MyMesh>

void measure2(MyMesh& mesh, double clr[])
{
	typename MyMesh::EdgeIter	e_it, e_begin(mesh.edges_begin()), e_end(mesh.edges_end());
	typename MyMesh::HalfedgeHandle hh0, hh0a, hh0b, hh1, hh1a, hh1b;
	typename MyMesh::VertexHandle	vt0a, vf0a, vt0b, vf0b, vt1a, vf1a, vt1b, vf1b;
	typename MyMesh::Point		pt0a, pf0a, pt0b, pf0b, pt1a, pf1a, pt1b, pf1b;
	double len0a, len0b, len1a, len1b;
	int id;

	//Conformality - cross length ratio
	for  (e_it = e_begin; e_it != e_end; ++e_it)
	{
	   hh0 = mesh.halfedge_handle(*e_it, 0);
	   hh1 = mesh.halfedge_handle(*e_it, 1);

	   hh0a = mesh.next_halfedge_handle(hh0);
	   vt0a = mesh.to_vertex_handle(hh0a);
	   vf0a = mesh.from_vertex_handle(hh0a);
	   pt0a = mesh.point(vt0a);
	   pf0a = mesh.point(vf0a);
	   len0a = (pt0a - pf0a).length();

	   hh0b = mesh.next_halfedge_handle(hh0a);
	   vt0b = mesh.to_vertex_handle(hh0b);
	   vf0b = mesh.from_vertex_handle(hh0b);
	   pt0b = mesh.point(vt0b);
	   pf0b = mesh.point(vf0b);
	   len0b = (pt0b - pf0b).length();

	   hh1a = mesh.next_halfedge_handle(hh1);
	   vt1a = mesh.to_vertex_handle(hh1a);
	   vf1a = mesh.from_vertex_handle(hh1a);
	   pt1a = mesh.point(vt1a);
	   pf1a = mesh.point(vf1a);
	   len1a = (pt1a - pf1a).length();

	   hh1b = mesh.next_halfedge_handle(hh1a);
	   vt1b = mesh.to_vertex_handle(hh1b);
	   vf1b = mesh.from_vertex_handle(hh1b);
	   pt1b = mesh.point(vt1b);
	   pf1b = mesh.point(vf1b);
	   len1b = (pt1b - pf1b).length();

	   id = (*e_it).idx();
	   clr[id] = (len0a*len1a)/(len0b*len1b);
	}
};


/* =============================================================================================================
   This function calculates write (i) three angles of each face (ii) cross length ratio into CSV
   ============================================================================================================= */

void originmeasure(double ang[][3], double clr[], int n_edges, int n_faces, std::string INfile)
{
	int i;

	std::ofstream fs1, fs2;
	std::string s = INfile;
	s.erase(s.find_last_of("."), std::string::npos);
	std::string filename1 = s + "_ang_origin" + ".csv";
	std::string filename2 = s + "_clr_origin" + ".csv";

	fs1.open(filename1);
	for(i = 0 ; i != n_faces ; i++ )
	{
	fs1 << ang[i][0] << "," << ang[i][1] << "," << ang[i][2] << "\n";
	}
	fs1.close();

	fs2.open(filename2);
	for(i = 0 ; i != n_edges ; i++ )
	{
	fs2 << clr[i] << "\n";
	}
	fs2.close();
}

/* =============================================================================================================
   This function calculates (i) angle ratio change, select maximal angular change in each triangle,
			    (ii) change in cross length ratio (clr)
			    and export them into CSV
	Input: arrays of angles, array of clr before/after (c)MCF		Output: CSV with max angle ratio for each face
   ============================================================================================================= */

void confmeasure(double ang[][3], double angn[][3], double clr[], double clrn[], int n_edges, int n_faces, std::string INfile, int iflow)
{
	int i;
	double temp1[n_faces], temp2[n_edges];
	for (i = 0; i != n_faces; i++){
		temp1[i] = std::max( std::abs(angn[i][0] - ang[i][0]) / ang[i][0], std::abs(angn[i][1] - ang[i][1]) / ang[i][1]);
		temp1[i] = std::max( temp1[i], std::abs(angn[i][2] - ang[i][2]) / ang[i][2]);}

	for (i = 0; i != n_edges; i++){
		temp2[i] = clrn[i] / clr[i];}

	std::ofstream fs1, fs2;
	std::string s = INfile;
	s.erase(s.find_last_of("."), std::string::npos);
	std::string filename1 = s + "_ang" + std::to_string(iflow) + "v.csv";
	std::string filename2 = s + "_clr" + std::to_string(iflow) + "v.csv";

	fs1.open(filename1);
	for(i = 0 ; i != n_faces ; i++ )
	{
	fs1 << temp1[i] << "\n";
	}
	fs1.close();

	fs2.open(filename2);
	for(i = 0 ; i != n_edges ; i++ )
	{
	fs2 << temp2[i] << "\n";
	}
	fs2.close();
}


/* =============================================================================================================
   This function output the identification numbers of the neighboring voxels.
	Input: mesh 		Output:  list of identification numbers
   ============================================================================================================= */

void voxelngh(int count, int xdim, int ydim, std::vector<int> *countlist)
{
	int block = (xdim)*(ydim);

	(*countlist).push_back (count);

	(*countlist).insert((*countlist).end(), {count-1, count+1});
	(*countlist).insert((*countlist).end(), {count-xdim-1, count-xdim, count-xdim+1, count+xdim-1, count+xdim, count+xdim+1});
	(*countlist).insert((*countlist).end(), {count+block-1, count+block, count+block+1});
	(*countlist).insert((*countlist).end(), {count+block-xdim-1, count+block-xdim, count+block-xdim+1});
	(*countlist).insert((*countlist).end(), {count+block+xdim-1, count+block+xdim, count+block+xdim+1});
	(*countlist).insert((*countlist).end(), {count-block-1, count-block, count-block+1});
	(*countlist).insert((*countlist).end(), {count-block-xdim-1, count-block-xdim, count-block-xdim+1});
	(*countlist).insert((*countlist).end(), {count-block+xdim-1, count-block+xdim, count-block+xdim+1});
}

/* =============================================================================================================
   This function counts the number of face intersection points
	Input: mesh 		Output: # face intersection points
   ============================================================================================================= */

template <typename MyMesh>

void checkface(MyMesh& mesh, int *j, int *f)
{
	*j =0; *f = 0;
	int n_vertices = mesh.n_vertices();
	int n_edges    = mesh.n_edges();
	int n_faces    = mesh.n_faces();

	typename MyMesh::VertexIter	v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
	typename MyMesh::EdgeIter	e_it, e_begin(mesh.edges_begin()), e_end(mesh.edges_end());
	typename MyMesh::Point		p;

	double  xmax, xmin, ymax, ymin, zmax, zmin;

	// Find the max. and min. x-, y-, z- coordinates
	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
	   p = mesh.point(*v_it);
		if (v_it == v_begin)
		{
	   		xmin = p[0]; ymin = p[1]; zmin = p[2];
	   		xmax = p[0]; ymax = p[1]; zmax = p[2];
		}
		else
		{
	  		xmin = std::min(xmin, p[0]); ymin = std::min(ymin, p[1]); zmin = std::min(zmin, p[2]);
	  		xmax = std::max(xmax, p[0]); ymax = std::max(ymax, p[1]); zmax = std::max(zmax, p[2]);
		}
	}


	double  emax = 0.0, emin = 0.0;
	typename MyMesh::HalfedgeHandle  h;
	typename MyMesh::VertexHandle	w0,w1;
	typename MyMesh::Point s0, s1;
	typename MyMesh::Scalar	 el;
	int block, dim, id, xdim, ydim, zdim, xpos, ypos, zpos;

	//Calculate the maximal edge length
	for (e_it = e_begin; e_it != e_end; ++e_it)
	{
	    typename MyMesh::EdgeHandle eh0(*e_it);
	    auto r = eh0.idx();

	    h = mesh.halfedge_handle(*e_it, 0);
	    w0 = mesh.to_vertex_handle(h);
	    w1 = mesh.from_vertex_handle(h);
	    s0 = mesh.point(w0);
	    s1 = mesh.point(w1);
	    el = (s0 - s1).length();
	    if (el > emax) emax = el;
	    if (r == 0)	 emin = el;
	    else emin = std::min(emin,el);
	}

	std::cout << "min x-coord. is : " << xmin << std::endl;
	std::cout << "max x-coord. is : " << xmax << std::endl;
	std::cout << "min y-coord. is : " << ymin << std::endl;
	std::cout << "max y-coord. is : " << ymax << std::endl;
	std::cout << "min z-coord. is : " << zmin << std::endl;
	std::cout << "max z-coord. is : " << zmax << std::endl;
	std::cout << "max. edge length: " << emax << std::endl;
	std::cout << "min. edge length: " << emin << std::endl;


	xmin = xmin - emax; xmax = xmax + emax;
	ymin = ymin - emax; ymax = ymax + emax;
	zmin = zmin - emax; zmax = zmax + emax;

	xdim = ceil((xmax - xmin)/emax);
	ydim = ceil((ymax - ymin)/emax);
	zdim = ceil((zmax - zmin)/emax);
	block = (xdim)*(ydim);
	dim = (xdim)*(ydim)*(zdim);

	//Two maps: voxel and voxlinverse
	//voxel - key: voxel, values: vertices
	//voxelinverse - key: vertex, values: voxel
	typename std::map<int, typename std::vector<typename MyMesh::VertexHandle>> voxel;
	typename std::map<typename MyMesh::VertexHandle, int> voxelinverse;

	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
	   p = mesh.point(*v_it);
	   xpos = floor((p[0] - xmin)/emax);
	   ypos = floor((p[1] - ymin)/emax);
	   zpos = floor((p[2] - zmin)/emax);
	   id = xpos + xdim*ypos + xdim*ydim*zpos;

	   if (voxel[id].empty())
	   {
		voxel[id].insert (voxel[id].begin(), *v_it);
	   }
	   else
	   {
	   	voxel[id].push_back (*v_it);
	   }

	   voxelinverse[*v_it] = id;
	}

	//Checking
	typename std::vector<typename MyMesh::VertexHandle>::iterator it;

	for (int k = 0; k < dim ; ++k)
	{
		if (!voxel[k].empty())
		{
	   	std::cout << "voxel[" << k << "] contains: " << std::endl;
			for (it = voxel[k].begin(); it != voxel[k].end(); it++)
			std::cout << ' ' << *it << std::endl;
		}
		else
		{
	   	std::cout << "voxel[" << k << "] is empty" << std::endl;
		}
	}


	typename MyMesh::FaceIter	f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
	typename MyMesh::FaceHalfedgeIter  fh_it;
	typename MyMesh::VertexEdgeIter	   ve_it;
	typename MyMesh::HalfedgeHandle  h0, h1, h00;
	typename MyMesh::VertexHandle	v0, v1, v2, vv0, vv1;
	typename MyMesh::Point		p0, p1, p2, pp, pp0, pp1, r, vec;
	typename MyMesh::Scalar		c, c1, c2, c3, t, check;
	int count1, count2, count3, count;

	std::vector<int> countlist;
	std::vector<int>::iterator itcount;
	typename std::vector<typename MyMesh::VertexHandle> list;
	typename std::vector<typename MyMesh::VertexHandle>::iterator it;

	for (f_it = f_begin; f_it != f_end; ++f_it)
	{
	    bool b = false;
	    fh_it = mesh.fh_iter(*f_it);
	    h0 = mesh.next_halfedge_handle(*fh_it);
	    v0 = mesh.to_vertex_handle(h0);
	    v1 = mesh.from_vertex_handle(h0);
	    p0 = mesh.point(v0);
	    p1 = mesh.point(v1);

	    h1 = *fh_it++;
	    v2 = mesh.from_vertex_handle(h1);
	    p2 = mesh.point(v2);
	    vec = cross(p0 - p1, p2 - p1); //surface normal
	    vec = vec/vec.norm();
	    c = dot(vec, p1);

	    //Assemble *counstlist: list of all neighboring voxels around the current face
	    count1 = voxelinverse.find(v0) -> second;
	    count2 = voxelinverse.find(v1) -> second;
	    count3 = voxelinverse.find(v2) -> second;

	    voxelngh(count1, xdim, ydim, &countlist);
	    voxelngh(count2, xdim, ydim, &countlist);
	    voxelngh(count3, xdim, ydim, &countlist);


	    std::sort(countlist.begin(), countlist.end());
	    itcount = std::unique(countlist.begin(), countlist.end());
	    countlist.resize( std::distance(countlist.begin(), itcount) );
	    //std::cout << "Countlist length: " << countlist.size() << std::endl;

	    //Assemble *list: list of all neighboring vertices around the current face
	    for (int n = 0; n < countlist.size(); n++)
	    {
		//std::cout << countlist[n] << ", " << std::endl;
		count = countlist[n];
	    	for (it = voxel[count].begin(); it != voxel[count].end(); it++)
	    	{
			if (std::find(list.begin(), list.end(), *it) == list.end())
			{
			list.push_back (*it);
			}
	    	}
	    }
	    countlist.clear();

	    for (it = list.begin(); it != list.end(); it++)
	    {
		   if (*it != v0 && *it != v1 && *it != v2 )
		   {
			//For each neigboring vertex, check all its incident edges for possible intersection with current face
			for (ve_it = mesh.ve_iter(*it); ve_it.is_valid(); ++ve_it)
			{
			// Ray-triangle intersection algorithm
		   	h00 = mesh.halfedge_handle(*ve_it, 0);
		   	vv0 = mesh.from_vertex_handle(h00);
		   	vv1 = mesh.to_vertex_handle(h00);

		   	if (vv0 != v0 && vv0 != v1 && vv0 != v2 && vv1 != v0 && vv1 != v1 && vv1 != v2)
		   	{
			pp0 = mesh.point(vv0);
		   	pp1 = mesh.point(vv1);
		   	r = pp1 - pp0;
		   	check = dot(vec,r);
		   	if ( check != 0)
		   	{
				t = (c - dot(vec,pp0))/check;

				if (t > 0 && t < 1)
				{
					pp = pp0 + t*r;
					c1 = dot(cross((p0 - p1),(pp - p1)), vec);
					c2 = dot(cross((p2 - p0),(pp - p0)), vec);
					c3 = dot(cross((p1 - p2),(pp - p2)), vec);
					if (c1 >= 0 && c2 >= 0 && c3>= 0)
					{
					b = true;
			   		++*j;
			   		//std::cout << "Intersection at face: " << *f_it << std::endl;
						if (*it == vv0)
						{if (std::find(list.begin(), list.end(), vv1) != list.end()){
							list.erase(std::find(list.begin(), list.end(), vv1));}
						}
						if (*it == vv1)
						{if (std::find(list.begin(), list.end(), vv0) != list.end()){
							list.erase(std::find(list.begin(), list.end(), vv0));}
						}*/


				}
				}
			}
		   	}
	    		}
		   }
	    }
	list.clear();
	if (b) ++*f;
	}

};

/* =============================================================================================================
   This function computes the surface area and the volume of the input mesh
	Input: mesh 		Output: surface area and volume
   ============================================================================================================= */
template <typename MyMesh>
void gauge(MyMesh& mesh, double* area, double* volume, double* spher, double* angdt, double ang0[], double ang1[], double ang2[])
{
	int n_faces = mesh.n_faces();
	double l0, l1, l2, cos0, cos1, cos2;
	double a0 = 0.0, a1 = 0.0, a2 = 0.0;

	double a_rs = 0.0, v_rs = 0.0;

	typename MyMesh::FaceIter			f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
	typename MyMesh::FaceHalfedgeIter	fh_it;
	typename MyMesh::HalfedgeHandle		h1;
	typename MyMesh::VertexHandle		v0, v1, v2;
	typename MyMesh::Point				p0, p1, p2;

	int j;

	for (f_it = f_begin, j = 0; f_it != f_end; ++f_it, ++j)
	{
		fh_it = mesh.fh_iter(*f_it);	 
	    	typename MyMesh::HalfedgeHandle h0(*fh_it);
		v0 = mesh.to_vertex_handle(h0);
		v1 = mesh.from_vertex_handle(h0);
		p0 = mesh.point(v0);
		p1 = mesh.point(v1);

		h1 = mesh.next_halfedge_handle(h0);
		v2 = mesh.to_vertex_handle(h1); 
		p2 = mesh.point(v2);

		l0 = (p1 - p2).length();
		l1 = (p2 - p0).length();
		l2 = (p0 - p1).length();

		//typename MyMesh::FaceHandle f(*f_it);
		a_rs = 0.25 * sqrt( (l2+l1+l0) * (l0-l2+l1) * (l0+l2-l1) * (l2+l1-l0) );
		v_rs = OpenMesh::dot(p1,cross(p0,p2))/6.0;
		*area += a_rs; 
		*volume += v_rs;

		a0 = acos((l1*l1 + l2*l2 - l0*l0)/ (2.0*l1*l2));
		a1 = acos((l0*l0 + l2*l2 - l1*l1)/ (2.0*l0*l2));
		a2 = acos((l0*l0 + l1*l1 - l2*l2)/ (2.0*l0*l1));
		*angdt += (abs(a0 - ang0[j]) + abs(a1 - ang1[j]) + abs(a2 - ang2[j]));
	}

	*spher = cbrt(M_PI*36* *volume * *volume)/ *area;

	*angdt /= n_faces;
}

/* =============================================================================================================
   This function return the area of the smallest triangle in the input mesh
	Input: mesh 		Output: smallest triangle
   ============================================================================================================= */
template <typename MyMesh>
void mt(MyMesh& mesh, double* st)
{
	int n_faces = mesh.n_faces();
	double l0, l1, l2;
	double a = 0.0, temp;

	typename MyMesh::FaceIter			f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
	typename MyMesh::FaceHalfedgeIter	fh_it;
	typename MyMesh::HalfedgeHandle		h1;
	typename MyMesh::VertexHandle		v0, v1, v2;
	typename MyMesh::Point				p0, p1, p2;

	*st = 0.0;
	for (f_it = f_begin; f_it != f_end; ++f_it)
	{
		fh_it = mesh.fh_iter(*f_it);	 
	    	typename MyMesh::HalfedgeHandle h0(*fh_it);
		v0 = mesh.to_vertex_handle(h0);
		v1 = mesh.from_vertex_handle(h0);
		p0 = mesh.point(v0);
		p1 = mesh.point(v1);

		h1 = mesh.next_halfedge_handle(h0);
		v2 = mesh.to_vertex_handle(h1); 
		p2 = mesh.point(v2);

		l0 = (p1 - p2).length();
		l1 = (p2 - p0).length();
		l2 = (p0 - p1).length();
		a = 0.25 * sqrt( (l2+l1+l0) * (l0-l2+l1) * (l0+l2-l1) * (l2+l1-l0) );
		if ((*f_it).idx() == 0) {temp = a;}
		else {temp = std::min(temp, a);}
	}
	*st = temp;
}

/* =============================================================================================================
   This function return the smallest and largest angles of the given mesh
	Input: mesh 		Output: list of angles 
   ============================================================================================================= */
template <typename MyMesh>
void measure_ang(MyMesh& mesh, std::string INfile, int niter)
{
	int n_faces = mesh.n_faces(), i;
	double ang[n_faces][3], area[n_faces];
	double l0, l1, l2;

	typename MyMesh::FaceIter	f_it, f_begin(mesh.faces_begin()), f_end(mesh.faces_end());
	typename MyMesh::FaceHalfedgeIter  fh_it;
	typename MyMesh::FaceVertexIter	fv_it;
	typename MyMesh::HalfedgeHandle h0, h1;
	typename MyMesh::VertexHandle	v0, v1, v2;
	typename MyMesh::Point		p0, p1, p2;

	for (f_it = f_begin; f_it != f_end; ++f_it)
	{
      	fh_it = mesh.fh_iter(*f_it);   
        	typename MyMesh::HalfedgeHandle h0(*fh_it);
      	v0 = mesh.to_vertex_handle(h0);
	    v1 = mesh.from_vertex_handle(h0);
	    p0 = mesh.point(v0);
	    p1 = mesh.point(v1);

	    h1 = mesh.next_halfedge_handle(h0);
	    v2 = mesh.to_vertex_handle(h1); 
	    p2 = mesh.point(v2);

	    l0 = (p1 - p2).length();
	    l1 = (p2 - p0).length();
	    l2 = (p0 - p1).length();

	    ang[(*f_it).idx()][0] = acos((l1*l1 + l2*l2 - l0*l0)/ (2.0*l1*l2))*(180/M_PI);
	    ang[(*f_it).idx()][1] = acos((l0*l0 + l2*l2 - l1*l1)/ (2.0*l0*l2))*(180/M_PI);
	    ang[(*f_it).idx()][2] = acos((l0*l0 + l1*l1 - l2*l2)/ (2.0*l0*l1))*(180/M_PI);

	   	area[(*f_it).idx()] = 0.25 * sqrt( (l2+l1+l0) * (l0-l2+l1) * (l0+l2-l1) * (l2+l1-l0) );
	}

	std::ofstream fs1;
	std::string s = INfile, filename1;
	s.erase(s.find_last_of("."), std::string::npos);
	if (niter == 0){
		filename1 = s + "_ang_start" + ".csv";
	}
	else{
		std::string num = std::to_string(niter);
		filename1 = s + "_ang_end_" + num + ".csv";
	}

	fs1.open(filename1);
	for(i = 0 ; i != n_faces ; i++ )
	{
	fs1 << ang[i][0] << "," << ang[i][1] << "," << ang[i][2] << ',' << area[i] << "\n";
	}
	fs1.close();
}


#endif
