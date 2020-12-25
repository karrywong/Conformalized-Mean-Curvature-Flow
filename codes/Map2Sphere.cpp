/* ===============================================================================================
   Map2Sphere: Conformal mapping of a genus zero surface represented by a triangular mesh onto
               a unified sphere using the "conformalized Mean Curvature Flow"
  
   Author:  Karry Wong
   Date:    04/04/2020
 
   =============================================================================================== */

/* ===============================================================================================
   Includes
   =============================================================================================== */

#include "Map2Sphere.h"
#include "CheckMesh.h"
#include "scale.h"
#include "Geometry.h"
#include "gaussmap.h"
#include "tuttemap.h"
#include "CenterSphere.h"
#include "cMCF.h"
#include <sys/time.h>


/* ===============================================================================================
   Define local acronyms
   =============================================================================================== */
struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;


struct timeval start, end;

/* ===============================================================================================
   Main program
   =============================================================================================== */

int main(int argc, char **argv)
{

/*	==========================================================================================
	Show usage if needed
	========================================================================================== */

	if( argc < 2 )
	{
		usage(argv);
		return -1;
	}

	std::string input = argv[1];
	if( input == "-h" || input == "-help" )
	{
		usage(argv);
		return -1;
	}

/*	==========================================================================================
	Read in all inputs (some values may have been preset)
	========================================================================================== */

	std::string INfile, OUTfile;
	std::string STRflow, STRweight, STRmax, STRstps, STRtol; //choice of flow
        
        if (!parse_args(argc, argv, &INfile, &OUTfile, &STRflow, &STRweight, &STRmax, &STRstps, &STRtol)) return 1; 

/*	==========================================================================================
	Create data structure for a new triangulated mesh    
	========================================================================================== */

	MyMesh mesh;

/*	==========================================================================================
	Read in mesh from input file
	========================================================================================== */

	OpenMesh::IO::Options ropt;
	ropt = OpenMesh::IO::Options::VertexColor;
	mesh.request_vertex_colors();

	omerr().disable();
	if ( ! OpenMesh::IO::read_mesh(mesh, INfile, ropt ) )
	{
		std::cerr << "Error: Cannot read mesh from " << INfile << std::endl;
		return 1;
	}
	omerr().enable();
	int nvertices, nfaces;
	if(! mesh_info(INfile,&nvertices, &nfaces))
	{
		std::cerr << "Error: Cannot read mesh: format not recognized (only OFF, OBJ, and PLY are accepted)." << std::endl;
		return 1; 
	}  
 
/*	==========================================================================================
	Print information about the mesh
	========================================================================================== */
	
	int n_vertices = mesh.n_vertices();
	int n_edges    = mesh.n_edges();
	int n_faces    = mesh.n_faces();
	double genus;
	genus = 1- (n_vertices - n_edges + n_faces)/2.0;

	std::cout << "# Vertices: " << n_vertices << std::endl;
	std::cout << "# Edges   : " << n_edges << std::endl;
	std::cout << "# Faces   : " << n_faces << std::endl;
	std::cout << "Genus     : " << genus << std::endl;

	// MyMesh::VertexIter	v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());

	// double ang[n_faces][3], clr[n_edges];
	// double area, volume, spher, var = 0.0;

	// std::cout <<  "Area: " << area << "\n" <<  "Volume: " << volume << "\n" << "Sphericity: " << spher << std::endl;
	// std::cout << "Min. angle: " << angmin << ",  min. aspect ratio: " << armin << std::endl;

/*	==========================================================================================
	Check face intersection
	==========================================================================================

	int j, f; 
	checkface(mesh, &j, &f);
	std::cout << "No. of intersection point: " << j << std::endl;
	std::cout << "No. of intersecting face: " << f << std::endl;
*/

/*	==========================================================================================
	Scale mesh such that its surface area is the same as a unit sphere (4*pi)
	==========================================================================================
*/

	scale(mesh);

/*  ==========================================================================================
	Conformalized Mean Curvature Flow proposed by Kazhdan, Solomon, and Ben-Chen:
	
	M. Kazhdan, J. Solomon, and M. Ben-Chen. "Can Mean Curvature Flow be Modified to be Non-singular?",
	Eurographics Symposium on Geometry Processing 2012.
    ==========================================================================================*/

	// iflow - 0: MCF,  1: cMCF,  2: cMCF w/ projection on sphere  3: Gauss Map initializer
	int			iflow;  
	if (!STRflow.empty()){iflow = std::stoi(STRflow);}
	else {iflow = 1;} //Default flow: cMCF

	if (iflow > 4){
	std::cerr << "Invalid choice of flow - 0: MCF  1:cMCF  2:cMCF w/ projection on sphere  3: Gauss Map initializer 4: Tutte embedding initializer" << std::endl;
	return 1;}

	int         iweight;
	if (!STRweight.empty()){iweight = std::stoi(STRweight);}
	else {iweight = 1;} //Default flow: cotangent weights
	if (iweight > 2){
	std::cerr << "Invalid choice of Tutte weight - 1: cotangent weights    2: graph lapacian" << std::endl;
	return 1;}

	double 			TOL;
	if(!STRtol.empty()){TOL = std::stod(STRtol);}
	else {TOL = 1e-3;} //Default TOL: 0.001

	double			spherTOL = 0.90; //For cMCF with projection onto sphere

	double			stps; 
	if(!STRstps.empty()){stps = std::stod(STRstps);}
	else {stps = 1e-2;} //Default step size: 0.01

	int 			nstpmax;	
	if(!STRmax.empty()){nstpmax = std::stoi(STRmax);}
	else {nstpmax = 64;} //Default no. of steps: 64

	//Assemble mass matrix D and stiffness matrix L
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(2* n_edges + n_vertices);  //nnz

	Eigen::SparseMatrix<double> D(n_vertices, n_vertices), L(n_vertices, n_vertices);
	Eigen::SparseMatrix<double> A(n_vertices, n_vertices);

	gettimeofday(&start, NULL);
	std::clock_t cpustart = clock();

		cMCF(iflow, iweight, mesh, tripletList, D, L, A, stps, nstpmax, TOL, spherTOL, INfile);

	std::clock_t cpuend = clock();
	gettimeofday(&end, NULL);

	double delta = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec)/ 1.e6;

	std::cout << "Elapsed time: " << delta << ", CPU Computing time: " << (float)(cpuend - cpustart) / CLOCKS_PER_SEC << std::endl;

	//MyMesh::VertexIter	v_it, v_begin(mesh.vertices_begin()), v_end(mesh.vertices_end());
	//Set vertex color for visualization
/*	for (v_it = v_begin; v_it != v_end; ++v_it)
	{
        if((*v_it).idx() >= 0 && (*v_it).idx() <= 11)
        {mesh.set_color(*v_it, MyMesh::Color(255,0,0));}
        else if ((*v_it).idx() >= 12 && (*v_it).idx() <= 41)
        {mesh.set_color(*v_it, MyMesh::Color(0,0,255));}
        else {mesh.set_color(*v_it, MyMesh::Color(192, 192, 192));}
	}
*/

	std::cout << "Conformalized MCF completed!" << std::endl;

/*	==========================================================================================
	i). Calculate angle ratio change, select maximal angular change in each triangle
	ii). Calculate the ratio change in cross length ratio (ideal = 1)
	Export in CSV
	========================================================================================== */

	//int i;
	//double test1, test2, test3, angn[n_faces][3], clrn[n_edges]; 
	//measure(mesh, &test1, &test2, &test3, angn, clrn); 
	//std::cout <<  "Area: " << test1 << ",  volume: " << test2 << ",  sphericity: " << test3 << std::endl;	

	//Write into CSV
	//originmeasure(ang, clr, n_edges, n_faces, INfile);	
	//confmeasure(ang, angn, clr, clrn, n_edges, n_faces, INfile, iflow);	

/*	==========================================================================================
	Write in csv file
	========================================================================================== 
	int count;
	double div = stps* nstpmax;
	std::ofstream fs;
	std::string s = INfile;
	s.erase(s.find_last_of("."), std::string::npos);
	std::string filename = s + "_test" + ".csv";

	fs.open(filename);
	for(count = 0; count < n_vertices ; count++)
	{
			fs << (x1[count] - x0[count])/(div) << "," << (y1[count] - y0[count])/(div) << "," << (z1[count] - z0[count])/(div) <<  "\n";
	}
*/
/*	==========================================================================================
	Write mesh in output file
	========================================================================================== */

	OpenMesh::IO::Options wopt;
	wopt = OpenMesh::IO::Options::VertexColor;
	if ( ! OpenMesh::IO::write_mesh(mesh, OUTfile, wopt) )
	{
		std::cerr << "Error: cannot write mesh to " << OUTfile << std::endl;
		return 1;
	}

	return 0;
}

static void usage(char** argv);

bool parse_args(int argc, char **argv, std::string *INfile, std::string *OUTfile, char* STRflow, char* STRmax, char* STRstps, char* STRtol);
