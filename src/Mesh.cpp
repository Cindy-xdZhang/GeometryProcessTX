#include "Mesh.h"
#include <tuple>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Eigen/sparse>
#include "viewer/ViewerData.h"

#include "utils/stl.h"
#include "utils/maths.h"
#include "utils/Eigen.h"
#include "utils/system.h"
#include<Eigen/IterativeLinearSolvers>	
#include<Eigen/SparseQR>
#include<Eigen/SparseLU>
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Trip;



size_t Mesh::counter = 0;

Mesh::Mesh()
{
    mID = counter++;
    mName = "";
    mRenderEdges = true;
    mRenderFlatFaces = false;
    mEdgeColor = Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
    mColorMap = colormap::ColorMapType::NUM_COLOR_MAP_TYPES;
}

bool Mesh::load(const std::string& filename)
{
    // Load a mesh
    mMesh.request_vertex_texcoords2D();
    OpenMesh::IO::Options opts(OpenMesh::IO::Options::VertexTexCoord);
    if (!OpenMesh::IO::read_mesh(mMesh, filename, opts))
    {
        std::cerr << "Error loading mesh from file " << filename << std::endl;
        return false;
    }

    // Mesh name
    mName = utils::remove_extension(utils::base_name(filename));

    // We need normals
    mMesh.request_face_normals();
    mMesh.request_vertex_normals();
    mMesh.update_normals();

    // Save original mesh
    mMeshOriginal = OpenMesh::Mesh(mMesh);

    // Success!
    return true;
}

bool Mesh::save(const std::string& filename)
{
    // Save a mesh
    if (!OpenMesh::IO::write_mesh(mMesh, filename))
    {
        std::cerr << "Error saving mesh to file " << filename << std::endl;
        return false;
    }

    // Success!
    return true;
}

std::string Mesh::name() const
{
    return mName;
}

unsigned int Mesh::id() const
{
    return mID;
}

void Mesh::updateViewerData(ViewerData& data) const
{
    // Clear viewer
    data.clear();
    data.clear_edges();

    // Convert mesh to viewer format
    Eigen::MatrixXd tmp_V;
    Eigen::MatrixXi tmp_F;
    Eigen::MatrixXi tmp_FtoF;
    Eigen::MatrixXd tmp_N;
    Eigen::MatrixXd tmp_UV;
    Eigen::MatrixXd P1, P2;
    OpenMesh::toRenderData(mMesh, tmp_V, tmp_F, tmp_FtoF, tmp_N, tmp_UV, P1, P2);

    // Plot the mesh
    data.set_mesh(tmp_V, tmp_F);
    data.FtoF = tmp_FtoF;
    data.set_normals(tmp_N);
    data.set_uv(tmp_UV);
    if (mRenderFlatFaces)
    {
        data.compute_normals();
    }
    else
    {
        data.face_based = false;
    }
    if (mRenderEdges)
    {
        Eigen::RowVector3d color;
        color << mEdgeColor[0], mEdgeColor[1], mEdgeColor[2];
        data.add_edges(P1, P2, color);
    }
    data.line_width = 2.0;
    data.show_lines = false;
    //show_texture = true;
    data.show_texture = false;

    // Colors
    data.uniform_colors(
            Eigen::Vector3d(51.0 / 255.0, 43.0 / 255.0, 33.3 / 255.0),
            Eigen::Vector3d(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0),
            Eigen::Vector3d(255.0 / 255.0, 235.0 / 255.0, 80.0 / 255.0));
}

OpenMesh::Mesh& Mesh::mesh()
{
    return mMesh;
}

OpenMesh::Mesh Mesh::mesh() const
{
    return mMesh;
}

void Mesh::clean()
{

}

void Mesh::resetMesh()
{
    mMesh = OpenMesh::Mesh(mMeshOriginal);
    clean();
}

void Mesh::setMesh(const OpenMesh::Mesh& mesh)
{
    mMesh = OpenMesh::Mesh(mesh);
    clean();
}

bool& Mesh::renderFlatFaces()
{
    return mRenderFlatFaces;
}

bool& Mesh::renderEdges()
{
    return mRenderEdges;
}

Eigen::Vector4f& Mesh::edgeColor()
{
    return mEdgeColor;
}

colormap::ColorMapType& Mesh::colormap()
{
    return mColorMap;
}

size_t Mesh::numVertices() const
{
    return mMesh.n_vertices();
}

Eigen::Vector3d Mesh::vertex(unsigned int index) const
{
    OpenMesh::Mesh::Point p = mMesh.point(mMesh.vertex_handle(index));
    return {p[0], p[1], p[2]};
}

void Mesh::setVertex(unsigned int index, const Eigen::Vector3d& v)
{
    OpenMesh::Mesh::Point p(v[0], v[1], v[2]);
    mMesh.set_point(mMesh.vertex_handle(index), p);
}

Eigen::Vector3d Mesh::faceCenter(unsigned int index) const
{
    OpenMesh::FaceHandle fh = mMesh.face_handle(index);
    OpenMesh::Mesh::Point p(0.0, 0.0, 0.0);
    for (OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.cfv_begin(fh);
         fv_it.is_valid(); ++fv_it)
    {
        p += mMesh.point(*fv_it);
    }
    p /= mMesh.valence(fh);
    return {p[0], p[1], p[2]};
}

void Mesh::setFaceCenter(unsigned int index, const Eigen::Vector3d& v)
{
    Eigen::Vector3d u = v - faceCenter(index);
    OpenMesh::Mesh::Point t(u[0], u[1], u[2]);
    OpenMesh::FaceHandle fh = mMesh.face_handle(index);
    for (OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_begin(fh);
         fv_it.is_valid(); ++fv_it)
    {
        mMesh.set_point(*fv_it, mMesh.point(*fv_it) + t);
    }
}

void Mesh::normalize()
{
    double totalArea = 0.0;
    std::vector<OpenMesh::Mesh::Point> barycenter(mMesh.n_faces());
    std::vector<double> area(mMesh.n_faces());

    // loop over faces
    for (OpenMesh::Mesh::FaceIter f_it = mMesh.faces_begin();
         f_it != mMesh.faces_end(); ++f_it)
    {
        // compute barycenter of face
        int valence = 0;
        OpenMesh::Mesh::Point center(0.0, 0.0, 0.0);
        for (OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(*f_it);
             fv_it.is_valid(); ++fv_it)
        {
            center += mMesh.point(*fv_it);
            ++valence;
        }
        barycenter[(*f_it).idx()] = center / valence;

        // compute area of face
        if (valence == 3)
        {
            OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(*f_it);
            OpenMesh::Mesh::Point v0 = mMesh.point(*fv_it);
            OpenMesh::Mesh::Point v1 = mMesh.point(*(++fv_it));
            OpenMesh::Mesh::Point v2 = mMesh.point(*(++fv_it));

            // A = 0.5 * || (v0 - v1) x (v2 - v1) ||
            double A = 0.5 * OpenMesh::length(OpenMesh::cross((v0 - v1), (v2 - v1)));
            area[(*f_it).idx()] = A;
            totalArea += area[(*f_it).idx()];
        }
        else if (valence == 4)
        {
            OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(*f_it);
            OpenMesh::Mesh::Point v0 = mMesh.point(*fv_it);
            OpenMesh::Mesh::Point v1 = mMesh.point(*(++fv_it));
            OpenMesh::Mesh::Point v2 = mMesh.point(*(++fv_it));
            OpenMesh::Mesh::Point v3 = mMesh.point(*(++fv_it));

            // A = 0.5 * || (v0 - v1) x (v2 - v1) ||
            double A012 = OpenMesh::length(OpenMesh::cross((v0 - v1), (v2 - v1)));
            double A023 = OpenMesh::length(OpenMesh::cross((v0 - v2), (v3 - v2)));
            double A013 = OpenMesh::length(OpenMesh::cross((v0 - v1), (v3 - v1)));
            double A123 = OpenMesh::length(OpenMesh::cross((v1 - v2), (v3 - v2)));
            area[(*f_it).idx()] = (A012 + A023 + A013 + A123) * 0.25;
            totalArea += area[(*f_it).idx()];
        }
        else
        {
            std::cerr << "Arbitrary polygonal faces not supported" << std::endl;
            return;
        }
    }

    // compute mesh centroid
    OpenMesh::Mesh::Point centroid(0.0, 0.0, 0.0);
    for (int i = 0; i < mMesh.n_faces(); i++)
    {
        centroid += area[i] / totalArea * barycenter[i];
    }

    // normalize mesh
    for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin();
         v_it != mMesh.vertices_end(); ++v_it)
    {
        OpenMesh::Mesh::Point p = mMesh.point(*v_it);
        p -= centroid; // subtract centroid (important for numerics)
        p /= sqrt(totalArea); // normalize to unit surface area (important for numerics)
        mMesh.set_point(*v_it, p);
    }
}

double Mesh::averageEdgeLength()
{
    double averageEdgeLength = 0.0;
    for (OpenMesh::Mesh::EdgeIter e_it = mMesh.edges_begin();
         e_it != mMesh.edges_end(); e_it++)
    {
        averageEdgeLength += mMesh.calc_edge_length(*e_it);
    }
    return averageEdgeLength / (double) mMesh.n_edges();
}

double Mesh::averageDihedralAngle()
{
    double averageDihedralAngle = 0.0;
    for (OpenMesh::Mesh::EdgeIter e_it = mMesh.edges_begin();
         e_it != mMesh.edges_end(); e_it++)
    {
        if (!mMesh.is_boundary(*e_it))
        {
            averageDihedralAngle += mMesh.calc_dihedral_angle(*e_it);
        }
    }
    return averageDihedralAngle / (double) mMesh.n_edges();
}

void Mesh::noise(double standardDeviation, NoiseDirection noiseDirection)
{
    double averageLength = averageEdgeLength();

    if (noiseDirection == NORMAL)
    {
        for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin();
             v_it != mMesh.vertices_end(); ++v_it)
        {
            OpenMesh::Mesh::Normal n = mMesh.normal(*v_it);
            double g = maths::gaussian(0, averageLength * standardDeviation);
            OpenMesh::Mesh::Point p = mMesh.point(*v_it) + n * g;
            mMesh.set_point(*v_it, p);
        }
    }

    if (noiseDirection == RANDOM)
    {
        for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin();
             v_it != mMesh.vertices_end(); ++v_it)
        {
            Eigen::Vector3d d = Eigen::randomVector();
            OpenMesh::Mesh::Normal n(d[0], d[1], d[2]);
            double g = maths::gaussian(0, averageLength * standardDeviation);
            OpenMesh::Mesh::Point p = mMesh.point(*v_it) + n * g;
            mMesh.set_point(*v_it, p);
        }
    }
}


void Mesh::LaplacianSmoothing(double lambda, int iterations, bool uniform )
{
    const int NVertices = numVertices();
   
	auto BuildLaplacianMatrix= [NVertices](const Mesh* inputMesh,bool uniform)->SpMat{
		std::vector<Trip> tripletList;
		tripletList.reserve(NVertices * 20);
        auto HalfEdgeMesh = inputMesh->mesh();

		auto computeLij_uniform = [NVertices, &HalfEdgeMesh](int i, int j)->double {

			OpenMesh::SmartVertexHandle Vi(i, &HalfEdgeMesh);
			OpenMesh::SmartVertexHandle Vj(j, &HalfEdgeMesh);
			for (auto eit : Vi.outgoing_halfedges()) {
				OpenMesh::SmartVertexHandle Vjj = eit.to();
				if (Vjj.idx() == j)
				{
					return 1.0;
				}
			}
			return 0.0;
		};

		auto computeLij_cotangent = [NVertices, &HalfEdgeMesh](int i, int j)->double {


			OpenMesh::SmartVertexHandle Vi(i, &HalfEdgeMesh);
			OpenMesh::SmartVertexHandle Vj(j, &HalfEdgeMesh);
			bool connected = false;
			OpenMesh::SmartHalfedgeHandle connectHalfEdge(-1);
			for (auto eit : Vi.outgoing_halfedges()) {
				OpenMesh::SmartVertexHandle Vjj = eit.to();
				if (Vjj.idx() == j)
				{
					connected = true;
					connectHalfEdge = eit;
					break;
				}
			}
			if (connected)
			{
				const  OpenMesh::SmartHalfedgeHandle PiPi_1 = connectHalfEdge.next();
				const  OpenMesh::SmartHalfedgeHandle PPiplus1 = connectHalfEdge.opp().next();
				const Eigen::Vector3d  VPi_1 = HalfEdgeMesh.point(PiPi_1.to());
				const  Eigen::Vector3d  VPiplus1 = HalfEdgeMesh.point(PPiplus1.to());
				const  Eigen::Vector3d  VPi = HalfEdgeMesh.point(Vj);
				const  Eigen::Vector3d  VP = HalfEdgeMesh.point(Vi);
				const Eigen::Vector3d Vec_Edge_Pi_pi_1 = VPi - VPi_1;
				const Eigen::Vector3d Vec_Edge_P_pi_1 = VP - VPi_1;
				const Eigen::Vector3d Vec_Edge_Pi_piplus1 = VPi - VPiplus1;
				const Eigen::Vector3d Vec_Edge_P_piplus1 = VP - VPiplus1;

				const double cosAlpha = Vec_Edge_Pi_pi_1.dot(Vec_Edge_P_pi_1) / (Vec_Edge_Pi_pi_1.norm() * Vec_Edge_P_pi_1.norm());
				const double sinAlpha = sqrt(1 - cosAlpha * cosAlpha);
				const double cosBeta = Vec_Edge_P_piplus1.dot(Vec_Edge_Pi_piplus1) / (Vec_Edge_Pi_piplus1.norm() * Vec_Edge_P_piplus1.norm());
				const double sinBeta = sqrt(1 - cosBeta * cosBeta);
				return (cosAlpha / sinAlpha + cosBeta / sinBeta);

			}
			else {
				return 0.0;
			}

       
		};
   
        std::function<double(int, int)> ComputeLij = computeLij_uniform;
        if(!uniform)ComputeLij = computeLij_cotangent;

        std::vector<double>SumVecs(NVertices,0.0);
		for (int i = 0; i < NVertices; i++) {
			for (int j = i + 1; j < NVertices; j++) {
				double Wij = ComputeLij(i, j);
				SumVecs[i] += Wij;
				SumVecs[j] += Wij;
				if (Wij == 0.0)
				{
					continue;
				}
				else
				{
					tripletList.emplace_back(Trip(i, j, Wij));
					tripletList.emplace_back(Trip(j, i, Wij));
				}

			}

		}

		for (int i = 0; i < NVertices; i++) {

			tripletList.emplace_back(Trip(i, i, -1.0 * SumVecs[i]));
		}

        SpMat mat(NVertices,NVertices);
		mat.setFromTriplets(tripletList.begin(), tripletList.end());
        return mat;
	};
 

	auto BuildAreaMatrix = [NVertices](const Mesh* inputMesh)->SpMat {
		std::vector<Trip> tripletList;
		tripletList.reserve(NVertices * 2);
		auto HalfEdgeMesh = inputMesh->mesh();

		auto computeInfluenceArea= [NVertices, &HalfEdgeMesh](int i)->double {

            OpenMesh::SmartVertexHandle Vi(i, &HalfEdgeMesh);
            const Eigen::Vector3d  thisPoint = HalfEdgeMesh.point(Vi);
			OpenMesh::Mesh::VertexVertexIter n_vv_it = HalfEdgeMesh.vv_iter(Vi);
			OpenMesh::Mesh::VertexVertexIter n_vv_it_last = n_vv_it;
			double sumQ = 0.0;
			Eigen::Vector3d FirstPoint = HalfEdgeMesh.point(*n_vv_it_last);
			for (n_vv_it++; n_vv_it.is_valid(); ++n_vv_it, ++n_vv_it_last) {
				Eigen::Vector3d LastPoint = HalfEdgeMesh.point(*n_vv_it_last);
				Eigen::Vector3d NextPoint = HalfEdgeMesh.point(*n_vv_it);

				Eigen::Vector3d  Edge1 = LastPoint - thisPoint;
				Eigen::Vector3d  Edge2 = NextPoint - thisPoint;
				double a = Edge1.norm();
				double b = Edge2.norm();
				double cosEdge = Edge1.dot(Edge2) / (a * b);
				//double edgeLength = sqrt(a * a + b * b - 2 * a * b * cosEdge);
                double TempQ = a * b * 0.5 * (sqrt(1 - cosEdge * cosEdge));
                sumQ += TempQ;
          
			}
			Eigen::Vector3d FinalPoint = HalfEdgeMesh.point(*n_vv_it_last);
			Eigen::Vector3d  Edge1 = FirstPoint - thisPoint;
			Eigen::Vector3d  Edge2 = FinalPoint - thisPoint;
			double a = Edge1.norm();
			double b = Edge2.norm();
			double cosEdge = Edge1.dot(Edge2) / (a * b);
			//double edgeLength = sqrt(a * a + b * b - 2 * a * b * cosEdge);
            double TempQ = a * b * 0.5 * (sqrt(1 - cosEdge * cosEdge));
            sumQ += TempQ;
        

            return sumQ;
		};

		for (int i = 0; i < NVertices; i++) {
            //todo: I've checked the problem doesn't come from the data range of computeInfluenceArea, not from re computation during iterations,
            //todo:no matter use sum length or area , sill have problem..->fixed take M not M-1
                //double Iver_VoronoiArea=1.0/( 2* computeInfluenceArea(i) ); 
          double Iver_VoronoiArea=( 2* computeInfluenceArea(i) );
				tripletList.emplace_back(Trip(i, i, Iver_VoronoiArea));
		
			}

		SpMat mat(NVertices, NVertices);
		mat.setFromTriplets(tripletList.begin(), tripletList.end());
		return mat;
	};


    auto BuildXvector = [NVertices](const Mesh* inputMesh)->Eigen::RowMatrixXd {
		Eigen::RowMatrixXd X_i(NVertices, 3);
		const auto HalfEdgeMesh = inputMesh->mesh();
		for (int i = 0; i < NVertices; i++)
		{
			Eigen::RowVector3d point = HalfEdgeMesh.point(OpenMesh::SmartVertexHandle(i, &HalfEdgeMesh));
			X_i.row(i) = point;
		}

        return X_i;
    };

	SpMat Laplacian = BuildLaplacianMatrix(this, uniform);
    SpMat M = BuildAreaMatrix(this);


	
	// Solving: u 
	SpMat A = (M - lambda * Laplacian);
	Eigen::SimplicialLDLT<SpMat> solver(A);
	Eigen::RowMatrixXd X_i = BuildXvector(this);

  
	while (iterations--)
	{
		Eigen::RowMatrixXd b = M * X_i;
		Eigen::RowMatrixXd x_iplus1 = solver.solve(b);
		X_i = x_iplus1;
	}
    

	for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin();
		v_it != mMesh.vertices_end(); ++v_it)
	{
		OpenMesh::Mesh::Point p = X_i.row(v_it->idx());
		mMesh.set_point(*v_it, p);
	}

}


void Mesh::OptimizingSmoothing(double lambda, double mu, double gama, double theta, int iterations, bool uniform) {

    const int NVertices = numVertices();



	auto BuildLaplacianMatrixTriplest = [gama,NVertices,uniform](const OpenMesh::Mesh& HalfEdgeMesh) {
		std::vector<Trip> tripletList;
		tripletList.reserve(3 * NVertices * NVertices);


		auto computeLij_uniform = [NVertices, &HalfEdgeMesh](int i, int j)->double {

			OpenMesh::SmartVertexHandle Vi(i, &HalfEdgeMesh);
			OpenMesh::SmartVertexHandle Vj(j, &HalfEdgeMesh);
			for (auto eit : Vi.outgoing_halfedges()) {
				OpenMesh::SmartVertexHandle Vjj = eit.to();
				if (Vjj.idx() == j)
				{
					return 1.0;
				}
			}
			return 0.0;
		};

		auto computeLij_cotangent = [NVertices, &HalfEdgeMesh](int i, int j)->double {


			OpenMesh::SmartVertexHandle Vi(i, &HalfEdgeMesh);
			OpenMesh::SmartVertexHandle Vj(j, &HalfEdgeMesh);
			bool connected = false;
			OpenMesh::SmartHalfedgeHandle connectHalfEdge(-1);
			for (auto eit : Vi.outgoing_halfedges()) {
				OpenMesh::SmartVertexHandle Vjj = eit.to();
				if (Vjj.idx() == j)
				{
					connected = true;
					connectHalfEdge = eit;
					break;
				}
			}
			if (connected)
			{
				const  OpenMesh::SmartHalfedgeHandle PiPi_1 = connectHalfEdge.next();
				const  OpenMesh::SmartHalfedgeHandle PPiplus1 = connectHalfEdge.opp().next();
				const Eigen::Vector3d  VPi_1 = HalfEdgeMesh.point(PiPi_1.to());
				const  Eigen::Vector3d  VPiplus1 = HalfEdgeMesh.point(PPiplus1.to());
				const  Eigen::Vector3d  VPi = HalfEdgeMesh.point(Vj);
				const  Eigen::Vector3d  VP = HalfEdgeMesh.point(Vi);
				const Eigen::Vector3d Vec_Edge_Pi_pi_1 = VPi - VPi_1;
				const Eigen::Vector3d Vec_Edge_P_pi_1 = VP - VPi_1;
				const Eigen::Vector3d Vec_Edge_Pi_piplus1 = VPi - VPiplus1;
				const Eigen::Vector3d Vec_Edge_P_piplus1 = VP - VPiplus1;

				const double cosAlpha = Vec_Edge_Pi_pi_1.dot(Vec_Edge_P_pi_1) / (Vec_Edge_Pi_pi_1.norm() * Vec_Edge_P_pi_1.norm());
				const double sinAlpha = sqrt(1 - cosAlpha * cosAlpha);
				const double cosBeta = Vec_Edge_P_piplus1.dot(Vec_Edge_Pi_piplus1) / (Vec_Edge_Pi_piplus1.norm() * Vec_Edge_P_piplus1.norm());
				const double sinBeta = sqrt(1 - cosBeta * cosBeta);
				return (cosAlpha / sinAlpha + cosBeta / sinBeta);

			}
			else {
				return 0.0;
			}


		};

		std::function<double(int, int)> ComputeLij = computeLij_uniform;
		if (!uniform)ComputeLij = computeLij_cotangent;

		std::vector<double>SumVecs(NVertices, 0.0);
		for (int i = 0; i < NVertices; i++) {
			for (int j = i + 1; j < NVertices; j++) {
				double Wij = ComputeLij(i, j) * gama;
				SumVecs[i] += Wij;
				SumVecs[j] += Wij;
				if (Wij == 0.0)
				{
					continue;
				}
				else
				{
					tripletList.emplace_back(Trip(i, j, Wij));
					tripletList.emplace_back(Trip(j, i, Wij));
				}

			}

		}

		for (int i = 0; i < NVertices; i++) {

			tripletList.emplace_back(Trip(i, i, -1.0 * SumVecs[i]));
		}
        return tripletList;
	};


    auto BuildJacobiMatrix = [ NVertices,theta,gama,lambda,uniform](const OpenMesh::Mesh& HalfEdgeMesh, const  std::vector<Trip>& LaplacianTripList) {
        std::vector<Trip> tripletList;
		tripletList.reserve(3 * NVertices * NVertices);
        //smoothingNess= laplacianX
        tripletList = LaplacianTripList;

        //promity1= (xi-x0)^2
        for (int i = NVertices; i < 2 * NVertices; i++)
        {
            tripletList.emplace_back(Trip(i, i- NVertices, lambda));
        }
        //promity2= (xi-xc)^2
        for (int i = 2 * NVertices; i < 3 * NVertices; i++)
        {
            tripletList.emplace_back(Trip(i, i-2* NVertices, theta));
        }

        SpMat mat(3 * NVertices,  NVertices);
        mat.setFromTriplets(tripletList.begin(), tripletList.end());
        return mat;
    };


	auto BuildXvector = [NVertices](const OpenMesh::Mesh& HalfEdgeMesh)->Eigen::RowMatrixXd {
		Eigen::RowMatrixXd X_i(NVertices, 3);
		for (int i = 0; i < NVertices; i++)
		{
			Eigen::RowVector3d point = HalfEdgeMesh.point(OpenMesh::SmartVertexHandle(i, &HalfEdgeMesh));
			X_i.row(i) = point;
		}

		return X_i;
	};


	auto BuildCVector = [NVertices](const  std::vector<Trip>& LaplacianTripList, Eigen::RowMatrixXd Xo, Eigen::RowMatrixXd Xc) {
		
        Eigen::RowMatrixXd d(3 * NVertices, 3); 
        d *= 0;

		SpMat L(NVertices, NVertices);
		L.setFromTriplets(LaplacianTripList.begin(), LaplacianTripList.end());
		Eigen::RowMatrixXd LX = L * Xc;

		d.block(0, 0, NVertices, 3) = LX.block(0, 0, NVertices, 3);

        d.block(NVertices, 0, NVertices, 3) = Xo.block(0, 0, NVertices, 3) -Xc.block(0, 0, NVertices, 3);

        /*for (int i = NVertices; i < 2 * NVertices; i++)
		{
            d.row(i) = Xc.row(i- NVertices) - Xo.row(i - NVertices);
		}

			for (int i = 2 * NVertices; i < 3 * NVertices; i++)
			{
			d.row(i)={0,0,0};
			}*/




		return d;
	};

	
 

    Eigen::RowMatrixXd X_o = BuildXvector(mMesh);
    Eigen::RowMatrixXd X_c = X_o;
  


	std::vector<Trip> L_trips = BuildLaplacianMatrixTriplest(mMesh);
	SpMat Jac = BuildJacobiMatrix(mMesh,L_trips);
	SpMat A = (Jac.transpose() * Jac).pruned();
	Eigen::SimplicialLDLT<SpMat> solver(A);
	Eigen::RowMatrixXd C = BuildCVector(L_trips, X_o, X_c);
	//solving:
	Eigen::RowMatrixXd b = -1.0 * Jac.transpose() * C;
	Eigen::RowMatrixXd Delta_V = solver.solve(b);
	X_c = Delta_V + X_c;
    //update
	for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin();
		v_it != mMesh.vertices_end(); ++v_it)
	{
		OpenMesh::Mesh::Point p = X_c.row(v_it->idx());
		mMesh.set_point(*v_it, p);
	}
	

	
}



void Mesh::OptimizingQuadMesh(std::vector<double> Paramters) {
    const  int iterationPerClick = Paramters[0];
         
	const int NumberOfVertices= numVertices();
	const int NumberOfFaces = mMesh.n_faces();

	//condition1: 4*nf + 1*nf ,  ,unknow 3*nv+3*nf
	//condition2: 4*nv + 1*nv  ,unknow 3*nv+3*nv
	//condition3: nf  ,unknow 3*nv
	//->Jacobi:  4*nf + 1*nf + 4*nv + 1*nv+nf   * ( 3*nv+3*nf+3*nv) (vertex coordinates( v1_x, v2_x,..., v1_y, v2_y,), face normals, per-vertex normals)
	
	//NofConstrains is unknown for now, because condition2 constrains depends on valence of vertex
	int NofConstrains = (6 * NumberOfFaces + 5 * NumberOfVertices);

	const int NofUnknows = (6 * NumberOfVertices + 3 * NumberOfFaces);


	//initialize matrix for all vertices
	Eigen::MatrixXd initialVi(NumberOfVertices, 3);
	for (int i = 0; i < NumberOfVertices; i++)
	{
		initialVi.row(i) = mMesh.point(OpenMesh::VertexHandle(i));
	}
	//initialize of face normals
	Eigen::MatrixXd initialFaceNormals(NumberOfFaces,3);
	for (int i = 0; i < NumberOfFaces; i++)
	{
		OpenMesh::FaceHandle fh = mMesh.face_handle(i);
		OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(fh);
	/*	Eigen::Vector3d n0 = mMesh.normal(*fv_it);
		Eigen::Vector3d n1 = mMesh.normal(*(++fv_it));
		Eigen::Vector3d n2 = mMesh.normal(*(++fv_it));
		Eigen::Vector3d n3 = mMesh.normal(*(++fv_it));	
		Eigen::Vector3d fnormal = (n0 + n1 + n2 + n3)/4;*/

		Eigen::Vector3d  v0 = mMesh.point(*fv_it) ;
		Eigen::Vector3d  v1 = mMesh.point( *(++fv_it));
		Eigen::Vector3d  v2 = mMesh.point (*(++fv_it));
		Eigen::Vector3d  v3 = mMesh.point (*(++fv_it));
		Eigen::Vector3d  diagonal_1 = v2 - v0;
		Eigen::Vector3d  diagonal_2 = v3 - v1;
		Eigen::Vector3d fnormal = diagonal_1.cross(diagonal_2);
	
		fnormal = fnormal.normalized();
		initialFaceNormals.row(i) = fnormal;
	}
	//initialize matrix for all vertices normals
	Eigen::MatrixXd initialVertexNormals(NumberOfVertices, 3);
	for (int i = 0; i < NumberOfVertices; i++)
	{
		initialVertexNormals.row(i) = mMesh.normal(OpenMesh::VertexHandle(i)).normalized();
	}


	auto Build_Xvector = [NofUnknows, NumberOfVertices, NumberOfFaces](const Eigen::MatrixXd& Vertexs, const Eigen::MatrixXd& currentFaceNormals,
		const Eigen::MatrixXd& VertexNormals)->Eigen::VectorXd {
			Eigen::VectorXd X_i(NofUnknows, 1);

			X_i.block(0, 0, NumberOfVertices, 1) = Vertexs.block(0, 0, NumberOfVertices, 1);
			X_i.block(NumberOfVertices, 0, NumberOfVertices, 1) = Vertexs.block(0, 1, NumberOfVertices, 1);
			X_i.block(NumberOfVertices * 2, 0, NumberOfVertices, 1) = Vertexs.block(0, 2, NumberOfVertices, 1);

			const int bas = NumberOfVertices*3;
			const int bas2 = NumberOfVertices * 3 + NumberOfFaces * 3;

			X_i.block(bas,                     0, NumberOfFaces, 1) = currentFaceNormals.block(0, 0, NumberOfFaces, 1);
			X_i.block(bas + NumberOfFaces,     0, NumberOfFaces, 1) = currentFaceNormals.block(0, 1, NumberOfFaces, 1);
			X_i.block(bas + NumberOfFaces * 2, 0, NumberOfFaces, 1) = currentFaceNormals.block(0, 2, NumberOfFaces, 1);

			X_i.block(bas2, 0, NumberOfVertices, 1) = VertexNormals.block(0, 0, NumberOfVertices, 1);
			X_i.block(bas2 + NumberOfVertices, 0, NumberOfVertices, 1) = VertexNormals.block(0, 1, NumberOfVertices, 1);
			X_i.block(bas2 + NumberOfVertices * 2, 0, NumberOfVertices, 1) = VertexNormals.block(0, 2, NumberOfVertices, 1);

			return X_i;
	};

	//groupid =0,1,2->vertex coordinates,face normals, per-vertex normals
	//xyz =0,1,2->x,y,z
	auto UnknowsAddressTrans = [NumberOfVertices, NumberOfFaces](int groupId, int inner_groupId, int xyz)->int {
		int result = 0;
		if (groupId == 0)
		{
			result = NumberOfVertices * xyz + inner_groupId;
		}
		else if (groupId == 1)
		{
			int bias = 3 * NumberOfVertices;
			result = bias + NumberOfFaces * xyz + inner_groupId;

		}
		else if (groupId == 2)
		{
			int bias = 3 * NumberOfVertices + 3 * NumberOfFaces;
			result = bias + NumberOfVertices * xyz + inner_groupId;
		}

		return result;
	};


	//todo: make Jacobi and C a function of only unknown X vector
	auto BuildJacobiMatandCVectors = [&](const Eigen::MatrixXd& Vertices,const Eigen::MatrixXd& FaceNormals,
		const Eigen::MatrixXd& VertexNormals) {
		const double w_c1 = Paramters[1];
		const double w_c2 = Paramters[2];
		const double w_c3 = Paramters[3];
		const double w_normal = Paramters[4];
		const double w_fairness = Paramters[5];
	
		std::vector<Trip> Jac_tripletList;
		Jac_tripletList.reserve(NofConstrains * NofUnknows);

		std::vector<Trip> C_tripletList;
		C_tripletList.reserve(NofConstrains * 1);
		
		//  todo:Condition1 
		for (int i = 0; i < 5*NumberOfFaces; i++)
		{

			if (i>=0&&i<NumberOfFaces)
			{
				const int Faceindex = i;
				OpenMesh::FaceHandle fh = mMesh.face_handle(Faceindex);
				Eigen::Vector3d nf = FaceNormals.row(Faceindex);
				OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(fh);

				OpenMesh::VertexHandle vh0 = *fv_it;
				OpenMesh::VertexHandle vh1 = *(++fv_it);
				OpenMesh::VertexHandle vh2 = *(++fv_it);
				OpenMesh::VertexHandle vh3 = *(++fv_it);
				Eigen::Vector3d v0 = Vertices.row(vh0.idx());
				Eigen::Vector3d  v1 = Vertices.row(vh1.idx());
				/*OpenMesh::Mesh::Point v2 = Vertices.row(vh2.idx());
				OpenMesh::Mesh::Point v3 = Vertices.row(vh3.idx());*/

				Eigen::Vector3d edge_i = v0 - v1;
				//vx
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh0.idx(),0),  w_c1* nf(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh1.idx(), 0), -w_c1* nf(0)));
				//vy
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh0.idx(), 1), w_c1 * nf(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh1.idx(), 1), -w_c1 * nf(1)));
				//vz
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh0.idx(), 2), w_c1 * nf(2)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh1.idx(), 2), -w_c1 * nf(2)));

				//nx ,ny,nz
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 0), w_c1 * edge_i(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 1), w_c1 * edge_i(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 2), w_c1 * edge_i(2)));

				//Constrain Vectors: Constrain1
				double C1 = w_c1 * edge_i.dot(nf);
				C_tripletList.emplace_back(Trip(i, 0, C1));

				continue;
			}
			if (i >= NumberOfFaces && i < 2*NumberOfFaces)
			{

				const int Faceindex = i-NumberOfFaces;
				OpenMesh::FaceHandle fh = mMesh.face_handle(Faceindex);
				Eigen::Vector3d nf = FaceNormals.row(Faceindex);
				OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(fh);

				OpenMesh::VertexHandle vh0 = *fv_it;
				OpenMesh::VertexHandle vh1 = *(++fv_it);
				OpenMesh::VertexHandle vh2 = *(++fv_it);
				OpenMesh::VertexHandle vh3 = *(++fv_it);
			
				Eigen::Vector3d  v1 = Vertices.row(vh1.idx());
				Eigen::Vector3d v2 = Vertices.row(vh2.idx());
				Eigen::Vector3d edge_i = v1 - v2;

				//vx
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh1.idx(), 0), w_c1 * nf(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh2.idx(), 0), -w_c1 * nf(0)));
				//vy
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh1.idx(), 1), w_c1 * nf(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh2.idx(), 1), -w_c1 * nf(1)));
				//vz
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh1.idx(), 2), w_c1 * nf(2)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh2.idx(), 2), -w_c1 * nf(2)));

				//nx ,ny,nz
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 0), w_c1 * edge_i(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 1), w_c1 * edge_i(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 2), w_c1 * edge_i(2)));


				//Constrain Vectors: Constrain1
				double C1 = w_c1 * edge_i.dot(nf);
				C_tripletList.emplace_back(Trip(i, 0, C1));
				continue;

			}
			if (i >= 2*NumberOfFaces && i < 3*NumberOfFaces)
			{
				const int Faceindex = i - 2*NumberOfFaces;
				OpenMesh::FaceHandle fh = mMesh.face_handle(Faceindex);
				Eigen::Vector3d nf = FaceNormals.row(Faceindex);
				OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(fh);

				OpenMesh::VertexHandle vh0 = *fv_it;
				OpenMesh::VertexHandle vh1 = *(++fv_it);
				OpenMesh::VertexHandle vh2 = *(++fv_it);
				OpenMesh::VertexHandle vh3 = *(++fv_it);

			
				Eigen::Vector3d v2 = Vertices.row(vh2.idx());
				Eigen::Vector3d  v3 = Vertices.row(vh3.idx());
				Eigen::Vector3d edge_i = v2 - v3;

				//vx
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh2.idx(), 0), w_c1* nf(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh3.idx(), 0), -w_c1 * nf(0)));
				//vy
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh2.idx(), 1), w_c1* nf(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh3.idx(), 1), -w_c1 * nf(1)));
				//vz
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh2.idx(), 2), w_c1* nf(2)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh3.idx(), 2), -w_c1 * nf(2)));

				//nx ,ny,nz
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 0), w_c1* edge_i(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 1), w_c1* edge_i(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 2), w_c1* edge_i(2)));


				//Constrain Vectors: Constrain1
				double C1 = w_c1 * edge_i.dot(nf);
				C_tripletList.emplace_back(Trip(i, 0, C1));
				continue;
			}
			if (i >= 3* NumberOfFaces && i < 4*NumberOfFaces)
			{
				const int Faceindex = i - 3 * NumberOfFaces;
				OpenMesh::FaceHandle fh = mMesh.face_handle(Faceindex);
				Eigen::Vector3d nf = FaceNormals.row(Faceindex);
				OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(fh);

				OpenMesh::VertexHandle vh0 = *fv_it;
				OpenMesh::VertexHandle vh1 = *(++fv_it);
				OpenMesh::VertexHandle vh2 = *(++fv_it);
				OpenMesh::VertexHandle vh3 = *(++fv_it);
		
				Eigen::Vector3d  v3 = Vertices.row(vh3.idx());
				Eigen::Vector3d v0 = Vertices.row(vh0.idx());
				Eigen::Vector3d edge_i = v3 - v0;

				//vx
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh3.idx(), 0), w_c1* nf(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh0.idx(), 0), -w_c1 * nf(0)));
				//vy
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh3.idx(), 1), w_c1* nf(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh0.idx(), 1), -w_c1 * nf(1)));
				//vz
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh3.idx(), 2), w_c1* nf(2)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(0, vh0.idx(), 2), -w_c1 * nf(2)));

				//nx ,ny,nz
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 0), w_c1* edge_i(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 1), w_c1* edge_i(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, Faceindex, 2), w_c1* edge_i(2)));



				//Constrain Vectors: Constrain1
				double C1 = w_c1 * edge_i.dot(nf);
				C_tripletList.emplace_back(Trip(i, 0, C1));
				continue;
			}

			//nf**2-1
			if (i >= 4 * NumberOfFaces && i < 5 * NumberOfFaces) {
				int faceindex = i - 4 * NumberOfFaces;
				Eigen::Vector3d nf = FaceNormals.row(faceindex);
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, faceindex, 0), w_normal * 2 * nf(0)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, faceindex, 1), w_normal * 2 * nf(1)));
				Jac_tripletList.emplace_back(Trip(i, UnknowsAddressTrans(1, faceindex, 2), w_normal * 2 * nf(2)));


				//Constrain Vectors: Constrain1
				double Ci = w_normal *(nf.dot(nf) - 1.0) ;
				C_tripletList.emplace_back(Trip(i, 0, Ci));
			}


		}
	

		//  todo:Condition2 
		int TotalValence = 0;
		const int ConstrainBias2_1 = 5 * NumberOfFaces;
		for (int i=0;i<NumberOfVertices;i++)
		{
			OpenMesh::SmartVertexHandle Vi(i, &mMesh);
			const Eigen::Vector3d  thisPoint = mMesh.point(Vi);
			const Eigen::Vector3d  vnf = VertexNormals.row(i);
			for (OpenMesh::Mesh::VertexVertexIter n_vv_it = mMesh.vv_iter(Vi); n_vv_it.is_valid(); ++n_vv_it) {
				const Eigen::Vector3d NeighborPoint = mMesh.point(*n_vv_it);
				const Eigen::Vector3d  edge_i = NeighborPoint - thisPoint;
				const int NeighborVid = (*n_vv_it).idx();
				const int CostrainId = ConstrainBias2_1 + TotalValence;
				//vx
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, NeighborVid, 0), w_c2* vnf(0)));
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0,	i, 0), -w_c2 * vnf(0)));
				//vy
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, NeighborVid, 1), w_c2* vnf(1)));
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, i, 1), -w_c2 * vnf(1)));
				//vz
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, NeighborVid, 2), w_c2* vnf(2)));
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, i, 2), -w_c2 * vnf(2)));

				//nx ,ny,nz
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(2, i, 0), w_c2* edge_i(0)));
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(2, i, 1), w_c2* edge_i(1)));
				Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(2, i, 2), w_c2* edge_i(2)));

				//Constrain Vectors: Constrain2
				double C2 = w_c2 * edge_i.dot(vnf);
				C_tripletList.emplace_back(Trip(CostrainId, 0, C2));

				TotalValence++;
			}
		}

		const int ConstrainBias2_2 = ConstrainBias2_1 + TotalValence;
		for (int i=0;i< NumberOfVertices;i++)
		{
			int constrain_id = ConstrainBias2_2 +i;
			const Eigen::Vector3d  vnf = VertexNormals.row(i);

			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(2, i, 0), w_normal * 2 * vnf(0)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(2, i, 1), w_normal * 2 * vnf(1)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(2, i, 2), w_normal * 2 * vnf(2)));


			//Constrain Vectors: Constrain1
			double Ci = w_normal * (vnf.dot(vnf) - 1.0);
			C_tripletList.emplace_back(Trip(constrain_id, 0, Ci));
		}

	

		//  todo:Condition3 
		const int ConstrainBias3 = ConstrainBias2_2 + NumberOfVertices;
		for (int i = 0; i < NumberOfFaces; i++)
		{
			int constrain_id = ConstrainBias3 + i;
			const int Faceindex = i;
			OpenMesh::FaceHandle fh = mMesh.face_handle(Faceindex);
			Eigen::Vector3d nf = FaceNormals.row(Faceindex);
			OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(fh);
			OpenMesh::VertexHandle vh0 = *fv_it;
			OpenMesh::VertexHandle vh1 = *(++fv_it);
			OpenMesh::VertexHandle vh2 = *(++fv_it);
			OpenMesh::VertexHandle vh3 = *(++fv_it);
			Eigen::Vector3d v0 = Vertices.row(vh0.idx());
			Eigen::Vector3d v1 = Vertices.row(vh1.idx());
			Eigen::Vector3d v2 = Vertices.row(vh2.idx());
			Eigen::Vector3d v3 = Vertices.row(vh3.idx());

			Eigen::Vector3d  diagonal_1 =  v0 -v2;
			Eigen::Vector3d  diagonal_2 =  v1 -v3;


			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh0.idx(), 0), w_c3 * 2 * diagonal_1(0)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh0.idx(), 1), w_c3 * 2 * diagonal_1(1)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh0.idx(), 2), w_c3 * 2 * diagonal_1(2)));

			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh2.idx(), 0), -w_c3 * 2 * diagonal_1(0)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh2.idx(), 1), -w_c3 * 2 * diagonal_1(1)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh2.idx(), 2), -w_c3 * 2 * diagonal_1(2)));

			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh1.idx(), 0), -w_c3 * 2 * diagonal_2(0)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh1.idx(), 1), -w_c3 * 2 * diagonal_2(1)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh1.idx(), 2), -w_c3 * 2 * diagonal_2(2)));
	
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh3.idx(), 0), w_c3 * 2 * diagonal_2(0)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh3.idx(), 1), w_c3 * 2 * diagonal_2(1)));
			Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, vh3.idx(), 2), w_c3 * 2 * diagonal_2(2)));

			//Constrain Vectors: Constrain1
			double Ci = w_c3 * (  diagonal_1.dot(diagonal_1) - diagonal_2.dot(diagonal_2));

			C_tripletList.emplace_back(Trip(constrain_id, 0, Ci));
		}
		//  todo:Fairness at every vertex, take neighbor two polylines compute difference
		int NpolyLineConstrains = 0;
		const int ConstrainBias_Fairness = 5 * NumberOfFaces + (TotalValence + NumberOfVertices) + NumberOfFaces;
		for (int i = 0; i < NumberOfVertices; i++)
		{
			OpenMesh::SmartVertexHandle Vi(i, &mMesh);
			const Eigen::Vector3d  thisPoint = mMesh.point(Vi);
			int Valience = Vi.valence();
			if (Valience ==4)
			{

				OpenMesh::Mesh::VertexVertexIter n_vv_it = mMesh.vv_iter(Vi);
				const int nV0_id = (*n_vv_it).idx();
				const Eigen::Vector3d nV0 = mMesh.point(*n_vv_it++);
				const int nV1_id = (*n_vv_it).idx();
				const Eigen::Vector3d nV1 = mMesh.point(*n_vv_it++);
				const int nV2_id = (*n_vv_it).idx();
				const Eigen::Vector3d nV2 = mMesh.point(*n_vv_it++);
				const int nV3_id = (*n_vv_it).idx();
				const Eigen::Vector3d nV3 = mMesh.point(*n_vv_it);
	

				//todo: vertical fairness
				const Eigen::Vector3d delta = nV3 + nV1 - 2 * thisPoint;
				int constrain_id = ConstrainBias_Fairness + NpolyLineConstrains;
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV3_id, 0), w_fairness * 2 * delta(0)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV3_id, 1), w_fairness * 2 * delta(1)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV3_id, 2), w_fairness * 2 * delta(2)));

				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV1_id, 0), w_fairness * 2 * delta(0)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV1_id, 1), w_fairness * 2 * delta(1)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV1_id, 2), w_fairness * 2 * delta(2)));

				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, i, 0), w_fairness * -4 * delta(0)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, i, 1), w_fairness * -4 * delta(1)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, i, 2), w_fairness * -4 * delta(2)));
				//Constrain Vectors: (v3+v1-v)
				double fairness_v = w_fairness * (delta).dot(delta);
				C_tripletList.emplace_back(Trip(constrain_id, 0, fairness_v));




				//todo: horizontal fairness Constrain Vectors: (v2+v0-v)
				constrain_id +=1;
				const Eigen::Vector3d delta2 = nV2 + nV0 - 2 * thisPoint;
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV2_id, 0), w_fairness * 2 * delta2(0)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV2_id, 1), w_fairness * 2 * delta2(1)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV2_id, 2), w_fairness * 2 * delta2(2)));

				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV0_id, 0), w_fairness * 2 * delta2(0)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV0_id, 1), w_fairness * 2 * delta2(1)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, nV0_id, 2), w_fairness * 2 * delta2(2)));

				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, i, 0), w_fairness * -4 * delta2(0)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, i, 1), w_fairness * -4 * delta2(1)));
				Jac_tripletList.emplace_back(Trip(constrain_id, UnknowsAddressTrans(0, i, 2), w_fairness * -4 * delta2(2)));
				//Constrain Vectors: (v3+v1-v)
				double fairness_h = w_fairness * (delta2).dot(delta2);
				C_tripletList.emplace_back(Trip(constrain_id, 0, fairness_h));

				NpolyLineConstrains += 2;
				}
		}
			


		//todo: get total constrains
		NofConstrains = 5 * NumberOfFaces + (TotalValence + NumberOfVertices) + NumberOfFaces + NpolyLineConstrains;


		return std::tuple(Jac_tripletList, C_tripletList);
	};







	Eigen::VectorXd X_o = Build_Xvector(initialVi, initialFaceNormals, initialVertexNormals);
	Eigen::VectorXd X_c = X_o;

	auto [JacList, CList] = BuildJacobiMatandCVectors(initialVi,initialFaceNormals, initialVertexNormals);
	SpMat Jac(NofConstrains, NofUnknows);
	Jac.setFromTriplets(JacList.begin(), JacList.end());
	SpMat C(NofConstrains, 1);
	C.setFromTriplets(CList.begin(), CList.end());



	SpMat JTJ= (Jac.transpose() * Jac).pruned();
	double maxCoeff = 0;
	for (int i=0;i<NofUnknows;i++)
	{
		maxCoeff = std::max(JTJ.coeffRef(i, i), maxCoeff);
	}
	double  lambda=1e-7;
	//JTJ+lambda*I gurantee not singular
	for (int i = 0; i < NofUnknows; i++)
	{
		JTJ.coeffRef(i, i) += lambda * maxCoeff;
	}

	Eigen::MatrixXd b = (-1.0 * Jac.transpose() * C);

	//solving:

	Eigen::SimplicialLDLT<SpMat> solver;//// decomposition failed for SimplicialLDLT
	solver.compute(JTJ);
	if (solver.info() != Eigen::Success) {
		// decomposition failed
		abort();
	}
	Eigen::VectorXd Delta_V = solver.solve(b);
	if (solver.info() != Eigen::Success) {
		// solving failed
		abort();
	}
	X_c = Delta_V + X_c;

	//update
	for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin();
		v_it != mMesh.vertices_end(); ++v_it)
	{
		Eigen::Vector3d Vx = { X_c(UnknowsAddressTrans(0,v_it->idx(),0),0),X_c(UnknowsAddressTrans(0,v_it->idx(),1),0) ,X_c(UnknowsAddressTrans(0,v_it->idx(),2),0) };
		mMesh.set_point(*v_it, Vx);
	}


}
