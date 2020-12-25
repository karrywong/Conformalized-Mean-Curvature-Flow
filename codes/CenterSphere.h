/* ==CENTERSPHERE.H ==============================================================================
 *
 *  Small functions that apply a Mobius transform on the vertices of a mesh embedded on the sphere,
 *  so that the center of gravity of those points sits at the center of the sphere
 *
 *	Apply Method recently proposed by Keenan Crane:
 *	
 *	A. Baden, K. Crane, and M. Kazhdan. "Moebius Registration", Eurographics Symp. Geom. Proc.
 *	(2018).
 *  Author:  Karry Wong
 *  Date:    04/04/2020
   =============================================================================================== */

#ifndef _CENTERSPHERE_H_
#define _CENTERSPHERE_H_

  /* ===== INCLUDES AND TYPEDEFS =================================================================
   *
   =============================================================================================== */

  /* ==CenterOfMass  ==============================================================================
  *
  *  Input:
  *	@Mesh:	mesh in OpenMesh format, with vertices on the surface of a sphere
  *  Output:
  *     COM: center of mass of vertices
  *
   =============================================================================================== */

  template <typename MyMesh>
  void CenterOfMass(MyMesh& Mesh, double *COM)
  {
        typename MyMesh::VertexIter     v_iter, v_start, v_end;
        typename MyMesh::Point point;

        v_start = Mesh.vertices_begin();
        v_end   = Mesh.vertices_end();

	double xg, yg, zg;

	xg = 0; yg = 0; zg = 0;
        for(v_iter = v_start; v_iter != v_end; ++v_iter)
        {
		point = Mesh.point(*v_iter);

		xg += point[0];
		yg += point[1];
		zg += point[2];
	}

	int n = Mesh.n_vertices();
	xg = xg/n; yg=yg/n; zg=zg/n;

	COM[0] = xg; COM[1] = yg; COM[2] = zg;

  }

  /* ==Jacobian ===================================================================================
  *
  *  Input:
  *	@Mesh:	mesh in OpenMesh format, with vertices on the surface of a sphere
  *  Output:
  *     J: Jacobian for inversion method
  *
   =============================================================================================== */

  template <typename MyMesh>
  void Jacobian(MyMesh& Mesh, double *J)
  {
        typename MyMesh::VertexIter     v_iter, v_start, v_end;
        typename MyMesh::Point point;

        v_start = Mesh.vertices_begin();
        v_end   = Mesh.vertices_end();

	for(int i = 0; i < 9; i++) J[i] = 0;

        for(v_iter = v_start; v_iter != v_end; ++v_iter)
        {
		point = Mesh.point(*v_iter);

		J[0] -= point[0]*point[0]; J[1] -= point[0]*point[1]; J[2] -= point[0]*point[2];
		J[4] -= point[1]*point[1]; J[5] -= point[1]*point[2];
		J[8] -= point[2]*point[2];
	}
	int n = Mesh.n_vertices();
	J[0] = 1.0 + J[0]/n; J[4] = 1.0 + J[4]/n; J[8] = 1.0 + J[8]/n;
	J[1]=J[1]/n; J[2]=J[2]/n; J[5]=J[5]/n;
	J[3] = J[1]; J[6] = J[2]; J[7] = J[5];
  }

  /* ==MatInv =====================================================================================
  *
  *  Input:
  *	Mat:	input 3x3 matrix (symmetric)
  *  Output:
  *	Inv:	inverse of Mat
  *
   =============================================================================================== */

  void MatInv(double *Mat, double *Inv)
  {

	Inv[0] = Mat[4]*Mat[8] - Mat[7]*Mat[5];
	Inv[1] = -(Mat[3]*Mat[8] - Mat[6]*Mat[5]);
	Inv[2] = Mat[3]*Mat[7] - Mat[6]*Mat[4];
	Inv[3] = Inv[1];
	Inv[4] = Mat[0]*Mat[8] - Mat[6]*Mat[2];
	Inv[5] = -(Mat[0]*Mat[7] - Mat[6]*Mat[1]);
	Inv[6] = Inv[2];
	Inv[7] = Inv[5];
	Inv[8] = Mat[0]*Mat[4] - Mat[3]*Mat[1];

	double det = Mat[0]*Inv[0] + Mat[1]*Inv[1] + Mat[2]*Inv[2];

	for(int i = 0; i < 9; i++)
	{
		Inv[i] = Inv[i]/det;
	}
  }

  /* ==MatVec =====================================================================================
  *
  *  Input:
  *	Mat:	input 3x3 matrix (symmetric)
  *	B:	vector (3D)
  *  Output:
  *	C:	C = Mat*B
  *
   =============================================================================================== */

  void MatVec(double *Mat, double *B, double *C)
  {
	C[0] = Mat[0]*B[0] + Mat[1]*B[1] + Mat[2]*B[2];
	C[1] = Mat[3]*B[0] + Mat[4]*B[1] + Mat[5]*B[2];
	C[2] = Mat[6]*B[0] + Mat[7]*B[1] + Mat[8]*B[2];
  }

  /* ==APPLYINVERSION ==============================================================================
  *
  *  Input:
  *	@Mesh:	mesh in OpenMesh format, with vertices on the surface of a sphere
  *     center: center for the inversion
  *  Output:
  *     @Mesh:  mesh with updated positions of the vertices, such that center of gravity
  *		of the point matches the center of the sphere
  *
   =============================================================================================== */

  template <typename MyMesh>
  void ApplyInversion(MyMesh& Mesh, double *center)
  {
        typename MyMesh::VertexIter     v_iter, v_start, v_end;
        typename MyMesh::Point point;

        v_start = Mesh.vertices_begin();
        v_end   = Mesh.vertices_end();

	double coef;
	double val, den;

	coef = 1.0 - center[0]*center[0] - center[1]*center[1] - center[2]*center[2];

        for(v_iter = v_start; v_iter != v_end; ++v_iter)
        {
		point = Mesh.point(*v_iter);

		den = (point[0]+center[0])*(point[0]+center[0]);
		den += (point[1]+center[1])*(point[1]+center[1]);
		den += (point[2]+center[2])*(point[2]+center[2]);

		val = coef/den;
		
		point[0] = val*(point[0]+center[0]) + center[0];
		point[1] = val*(point[1]+center[1]) + center[1];
		point[2] = val*(point[2]+center[2]) + center[2];

		Mesh.set_point( *v_iter, point);
	}
  }

  /* ==CENTERSPHERE ================================================================================
  *
  *  Input:
  *	@Mesh:	mesh in OpenMesh format, with vertices on the surface of a sphere
  *  Output:
  *     @Mesh:  mesh with updated positions of the vertices, such that center of gravity
  *		of the point matches the center of the sphere
  *
   =============================================================================================== */

  template <typename MyMesh>
  void CenterSphere(MyMesh& Mesh)
  {

  /* ===============================================================================================
	Apply Method recently proposed by Crane:
	
	A. Baden, K. Crane, and M. Kazhdan. "Moebius Registration", Eurographics Symp. Geom. Proc.
	(2018).
   =============================================================================================== */

	double TOL = 1.e-9; //default = 1e-8

	int istep = 0;
	double dist;
	double J[9], Jinv[9], COM[3], Center[3], dC[3];

  /* 	==========================================================================================
	Initialize center to center of the ball; corresponding Inversion is Identity
     	========================================================================================== */

	Center[0] = 0.; Center[1] = 0.; Center[2] = 0.;
	ApplyInversion(Mesh, Center);

  /* 	==========================================================================================
	Gauss - Newton iterations for finding appropriate center of inversion
     	========================================================================================== */

	int flag = 1;
	while(true)
	{

		istep++;

  /* 	==========================================================================================
		Compute center of mass; exit loop if distance to origin < TOL
     	========================================================================================== */

  		CenterOfMass(Mesh, COM);
		dist = COM[0]*COM[0] + COM[1]*COM[1] + COM[2]*COM[2];
		std::cout << "istep = " << istep << " dist = " << std::sqrt(dist) << std::endl;
		if(std::sqrt(dist) < TOL) break;

  /* 	==========================================================================================
		Compute Jacobian and its inverse
     	========================================================================================== */

		Jacobian(Mesh, J);
		MatInv(J, Jinv);

  /* 	==========================================================================================
		Compute new Inversion center
		If new center outside the unit ball, bring it back in it
     	========================================================================================== */

		MatVec(Jinv, COM, dC);
		Center[0] = -0.5*dC[0]; Center[1] = -0.5*dC[1]; Center[2] = -0.5*dC[2];
		dist = std::sqrt(Center[0]*Center[0] + Center[1]*Center[1] + Center[2]*Center[2]);
		if(dist > 1) {
			Center[0] = 0.5*Center[0]/dist; Center[1] = 0.5*Center[1]/dist; Center[2]=0.5*Center[2]/dist;
		}

  /* 	==========================================================================================
			Apply Inversion
     	========================================================================================== */

		ApplyInversion(Mesh, Center);

	}

	std::cout << "Sphere vertices centered      : Success" << std::endl;
	std::cout << "Number of steps needed        : " << istep << std::endl;
	std::cout << "Final energy ( i.e. d(CG, O)) : " << std::sqrt(dist) << std::endl;
	std::cout << " " << std::endl;

  }

#endif
