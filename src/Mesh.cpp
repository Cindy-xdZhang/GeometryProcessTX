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
//#include<algorithm>
//#include<execution>
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Trip;


#define FillJacobi(cid,g,i,V3d) do{\
		if (notfixed(g,i)){\
		Jac_tripletList.emplace_back(Trip((cid), (UnknowsAddressTrans(g, i, 0)), (V3d(0)))); \
		Jac_tripletList.emplace_back(Trip((cid), (UnknowsAddressTrans(g, i, 1)), (V3d(1)))); \
		Jac_tripletList.emplace_back(Trip((cid), (UnknowsAddressTrans(g, i, 2)), (V3d(2))));	}}while(0)

#define g_VertexCoordinate(fid,Xc)  Eigen::Vector3d( Xc(UnknowsAddressTrans(0, fid, 0),0), Xc(UnknowsAddressTrans(0, fid, 1),0),  Xc(UnknowsAddressTrans(0, fid, 2),0) )
#define g_FaceNormal(fid,Xc)  Eigen::Vector3d( Xc(UnknowsAddressTrans(1, fid, 0),0),Xc(UnknowsAddressTrans(1, fid, 1),0),Xc(UnknowsAddressTrans(1, fid, 2),0) )
#define g_VertexNormal(fid,Xc)  Eigen::Vector3d( Xc(UnknowsAddressTrans(2, fid, 0),0),Xc(UnknowsAddressTrans(2, fid, 1),0),Xc(UnknowsAddressTrans(2, fid, 2),0) )
#define  g_EdgeNormal(fid,Xc)  Eigen::Vector3d( Xc(UnknowsAddressTrans(3, fid, 0),0),Xc(UnknowsAddressTrans(3, fid, 1),0),Xc(UnknowsAddressTrans(3, fid, 2),0) )

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
OpenMesh::Mesh Mesh::Originmesh() 
{
	return mMeshOriginal;
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



void Mesh::OptimizingQuadMesh(std::vector<double> Paramters,std::vector<int>FixVertexIds) {
      int iterationPerClick = Paramters[0];
         
	const int NumberOfVertices= numVertices();
	const int NumberOfFaces = mMesh.n_faces();
	const int NumberOfEdges= mMesh.n_edges();
	//condition1: 4*nf + 1*nf ,  ,unknow 3*nv+3*nf
	//condition2: 4*nv + 1*nv  ,unknow 3*nv+3*nv
	//condition3: nf  ,unknow 3*nv
	// bonus: 3 ne, unknow: ne
	//->Jacobi:  4*nf + 1*nf + 4*nv + 1*nv+nf   * ( 3*nv+3*nf+3*nv) (vertex coordinates( v1_x, v2_x,..., v1_y, v2_y,), face normals, per-vertex normals)
	
	//NofConstrains is unknown for now, because condition2 constrains depends on valence of vertex
	int NofConstrains = (6 * NumberOfFaces + 5 * NumberOfVertices);

	//vx,vy,vz, fnx,fny,fnz, vnx,vny,vnz, enx,eny,enz
	const int NofUnknows = (6 * NumberOfVertices + 3 * NumberOfFaces) + 3 * NumberOfEdges;


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
	//initialize matrix for all edges normals
	Eigen::MatrixXd initialEdgeNormals(NumberOfEdges, 3);
	for (int i = 0; i < NumberOfEdges; i++)
	{
		OpenMesh::SmartEdgeHandle eh(i,&mMesh);
		Eigen::Vector3d V0 = mMesh.point(eh.v0());
		Eigen::Vector3d V1 = mMesh.point(eh.v1());
		Eigen::Vector3d N1 = mMesh.normal(eh.v1());
		Eigen::Vector3d iniEN = (V0 - V1).cross(N1);

		initialEdgeNormals.row(i) = iniEN.normalized();
	}

	auto notfixed = [&FixVertexIds](const int groupId, const int inner_groupId) ->bool {
		if (groupId == 0)
		{
			for (auto& vid : FixVertexIds) {
				if (inner_groupId == vid)
				{
					return false;
				}
			}

		}
		return true;
	};

	auto UnknowsAddressTrans = [NumberOfVertices, NumberOfEdges, NumberOfFaces](const int groupId, const int inner_groupId, const int xyz)->int {

		//groupid =0,1,2,3->vertex coordinates,face normals, per-vertex normals, per edge normal for torsion
		//xyz =0,1,2->x,y,z
		int result = -1;
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
		else if (groupId == 3)
		{
			int bias = 3 * NumberOfVertices + 3 * NumberOfFaces + 3 * NumberOfVertices;
			result = bias + NumberOfEdges * xyz + inner_groupId;
		}
		else
		{
			abort();
		}
		return result;
	};
	//build Build_Xvector
	Eigen::VectorXd X_o(NofUnknows, 1);
	{
	X_o.block(0, 0, NumberOfVertices, 1) = initialVi.block(0, 0, NumberOfVertices, 1);
	X_o.block(NumberOfVertices, 0, NumberOfVertices, 1) = initialVi.block(0, 1, NumberOfVertices, 1);
	X_o.block(NumberOfVertices * 2, 0, NumberOfVertices, 1) = initialVi.block(0, 2, NumberOfVertices, 1);
	const int bas = NumberOfVertices*3;
	const int bas2 = NumberOfVertices * 3 + NumberOfFaces * 3;
	const int bas3 = NumberOfVertices * 3 + NumberOfFaces * 3 + NumberOfVertices * 3;

	X_o.block(bas,                     0, NumberOfFaces, 1) = initialFaceNormals.block(0, 0, NumberOfFaces, 1);
	X_o.block(bas + NumberOfFaces,     0, NumberOfFaces, 1) = initialFaceNormals.block(0, 1, NumberOfFaces, 1);
	X_o.block(bas + NumberOfFaces * 2, 0, NumberOfFaces, 1) = initialFaceNormals.block(0, 2, NumberOfFaces, 1);

	X_o.block(bas2, 0, NumberOfVertices, 1) = initialVertexNormals.block(0, 0, NumberOfVertices, 1);
	X_o.block(bas2 + NumberOfVertices, 0, NumberOfVertices, 1) = initialVertexNormals.block(0, 1, NumberOfVertices, 1);
	X_o.block(bas2 + NumberOfVertices * 2, 0, NumberOfVertices, 1) = initialVertexNormals.block(0, 2, NumberOfVertices, 1);
		

	//this is edge normal for torsion free  

	X_o.block(bas3, 0, NumberOfEdges, 1) = initialEdgeNormals.block(0, 0, NumberOfEdges, 1);
	X_o.block(bas3 + NumberOfEdges, 0, NumberOfEdges, 1) = initialEdgeNormals.block(0, 1, NumberOfEdges, 1);
	X_o.block(bas3 + NumberOfEdges * 2, 0, NumberOfEdges, 1) = initialEdgeNormals.block(0, 2, NumberOfEdges, 1);
	}



	


	//todo: make Jacobi and C a function of only unknown X vector
	auto BuildJacobiMatandCVectors = [&](const Eigen::MatrixXd& X_c,
		std::vector<Trip>& Jac_tripletList, std::vector<Trip> &C_tripletList) {
		const double w_c1 = Paramters[1];
		const double w_c2 = Paramters[2];
		const double w_c3 = Paramters[3];
		const double w_normal = Paramters[4];
		const double w_fairness = Paramters[5];
		const double w_torsion = Paramters[6];
		//std::vector<Trip> Jac_tripletList;
		Jac_tripletList.reserve(NofConstrains * NofUnknows);

		//std::vector<Trip> C_tripletList;
		C_tripletList.reserve(NofConstrains * 2);

		//  todo:Condition1 
		for (int i = 0; i <  NumberOfFaces && w_c1>0; i++)
		{
			const int Faceindex = i;
			OpenMesh::FaceHandle fh = mMesh.face_handle(Faceindex);
			Eigen::Vector3d nf = g_FaceNormal(Faceindex,X_c);
			OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(fh);

			OpenMesh::VertexHandle vh0 = *fv_it;
			OpenMesh::VertexHandle vh1 = *(++fv_it);
			OpenMesh::VertexHandle vh2 = *(++fv_it);
			OpenMesh::VertexHandle vh3 = *(++fv_it);
			Eigen::Vector3d v0 =  g_VertexCoordinate(vh0.idx(), X_c);
			Eigen::Vector3d  v1 = g_VertexCoordinate(vh1.idx(), X_c);
			OpenMesh::Mesh::Point v2 = g_VertexCoordinate(vh2.idx(), X_c);
			OpenMesh::Mesh::Point v3 = g_VertexCoordinate(vh3.idx(), X_c);

			const Eigen::Vector3d edge_0 = v0 - v1;
			const Eigen::Vector3d edge_1= v1 - v2;
			const Eigen::Vector3d edge_2 = v2 - v3;
			const Eigen::Vector3d edge_3 = v3 - v0;

			const double C0 = w_c1 * edge_0.dot(nf);
			FillJacobi(i, 0, vh0.idx(), w_c1* nf);
			FillJacobi(i, 0, vh1.idx(), -w_c1 * nf);
			//nx ,ny,nz
			FillJacobi(i, 1, Faceindex, w_c1* edge_0);
			C_tripletList.emplace_back(Trip(i, 0, C0));

			const double C1 = w_c1 * edge_1.dot(nf);
			FillJacobi(i + NumberOfFaces, 0, vh1.idx(), w_c1* nf);
			FillJacobi(i + NumberOfFaces, 0, vh2.idx(), -w_c1 * nf);
			//nx ,ny,nz
			FillJacobi(i + NumberOfFaces, 1, Faceindex, w_c1* edge_1);
			C_tripletList.emplace_back(Trip(i + NumberOfFaces, 0, C1));

			const double C2 = w_c1 * edge_2.dot(nf);
			FillJacobi(i + 2 * NumberOfFaces, 0, vh2.idx(), w_c1* nf);
			FillJacobi(i + 2 * NumberOfFaces, 0, vh3.idx(), -w_c1 * nf);
			//nx ,ny,nz
			FillJacobi(i + 2 * NumberOfFaces, 1, Faceindex, w_c1* edge_2);
			C_tripletList.emplace_back(Trip(i + 2 * NumberOfFaces, 0, C2));

			const double C3 = w_c1 * edge_3.dot(nf);
			FillJacobi(i + 3 * NumberOfFaces, 0, vh3.idx(), w_c1* nf);
			FillJacobi(i + 3 * NumberOfFaces, 0, vh0.idx(), -w_c1 * nf);
			//nx ,ny,nz
			FillJacobi(i + 3 * NumberOfFaces, 1, Faceindex, w_c1* edge_3);
			C_tripletList.emplace_back(Trip(i+3*NumberOfFaces, 0, C3));


			//Constrain  nf**2-1
			double Ci = w_normal * (nf.dot(nf) - 1.0);
			C_tripletList.emplace_back(Trip(i+4* NumberOfFaces, 0, Ci));
			FillJacobi(i+4*NumberOfFaces, 1, Faceindex, w_normal * 2 * nf);

		}


		//  todo:Condition2 
		int TotalValence = 0;
		const int ConstrainBias2_1 = 5 * NumberOfFaces;
		for (int i=0;i<NumberOfVertices && w_c2>0;i++)
		{
			OpenMesh::SmartVertexHandle Vi(i, &mMesh);
			const Eigen::Vector3d  thisPoint =  g_VertexCoordinate( i,X_c);

			const Eigen::Vector3d  vnf =  g_VertexNormal(i, X_c);
			for (OpenMesh::Mesh::VertexVertexIter n_vv_it = mMesh.vv_iter(Vi); n_vv_it.is_valid(); ++n_vv_it) {
				const int NeighborVid = (*n_vv_it).idx();
				const Eigen::Vector3d NeighborPoint = g_VertexCoordinate(NeighborVid, X_c);
				const Eigen::Vector3d  edge_i = NeighborPoint - thisPoint;
				const int CostrainId = ConstrainBias2_1 + TotalValence;
				//vx
				//Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, NeighborVid, 0), w_c2* vnf(0)));
				//Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0,	i, 0), -w_c2 * vnf(0)));
				////vy
				//Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, NeighborVid, 1), w_c2* vnf(1)));
				//Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, i, 1), -w_c2 * vnf(1)));
				////vz
				//Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, NeighborVid, 2), w_c2* vnf(2)));
				//Jac_tripletList.emplace_back(Trip(CostrainId, UnknowsAddressTrans(0, i, 2), -w_c2 * vnf(2)));
				FillJacobi(CostrainId, 0, NeighborVid, w_c2* vnf);
				FillJacobi(CostrainId, 0, i,-w_c2* vnf);
				//nx ,ny,nz
				FillJacobi(CostrainId, 2, i, w_c2* edge_i);

				//Constrain Vectors: Constrain2
				double C2 = w_c2 * edge_i.dot(vnf);
				C_tripletList.emplace_back(Trip(CostrainId, 0, C2));

				TotalValence++;
			}
		}

		const int ConstrainBias2_2 = ConstrainBias2_1 + TotalValence;
		for (int i=0;i< NumberOfVertices && w_c2>0;i++)
		{
			int constrain_id = ConstrainBias2_2 +i;
			const Eigen::Vector3d  vnf = g_VertexNormal(i, X_c);
			FillJacobi(constrain_id, 2, i, w_normal * 2 * vnf);

			//Constrain Vectors: Constrain1
			double Ci = w_normal * (vnf.dot(vnf) - 1.0);
			C_tripletList.emplace_back(Trip(constrain_id, 0, Ci));
		}

	

		//  todo:Condition3 
		const int ConstrainBias3 = ConstrainBias2_2 + NumberOfVertices;
		for (int i = 0; i < NumberOfFaces && w_c3>0; i++)
		{
			int constrain_id = ConstrainBias3 + i;
			const int Faceindex = i;
			OpenMesh::FaceHandle fh = mMesh.face_handle(Faceindex);
			Eigen::Vector3d nf = g_FaceNormal(Faceindex,X_c);
			OpenMesh::Mesh::FaceVertexIter fv_it = mMesh.fv_iter(fh);
			OpenMesh::VertexHandle vh0 = *fv_it;
			OpenMesh::VertexHandle vh1 = *(++fv_it);
			OpenMesh::VertexHandle vh2 = *(++fv_it);
			OpenMesh::VertexHandle vh3 = *(++fv_it);

			const Eigen::Vector3d v0 = g_VertexCoordinate(vh0.idx(), X_c);
			const Eigen::Vector3d  v1 = g_VertexCoordinate(vh1.idx(), X_c);
			const Eigen::Vector3d  v2 = g_VertexCoordinate(vh2.idx(), X_c);
			const Eigen::Vector3d  v3 = g_VertexCoordinate(vh3.idx(), X_c);


			const Eigen::Vector3d  diagonal_1 =  v0 -v2;
			const Eigen::Vector3d  diagonal_2 =  v1 -v3;


			Eigen::Vector3d  d1 = w_c3 * 2 * diagonal_1;
			Eigen::Vector3d  d2 = w_c3 * 2 * diagonal_2;
			FillJacobi(constrain_id,0, vh0.idx(), d1);
			FillJacobi(constrain_id,0, vh2.idx(), -d1);
			FillJacobi(constrain_id, 0, vh1.idx(), -d2);
			FillJacobi(constrain_id, 0, vh3.idx(), d2);
			//Constrain Vectors: Constrain1
			double Ci = w_c3 * (  diagonal_1.dot(diagonal_1) - diagonal_2.dot(diagonal_2));

			C_tripletList.emplace_back(Trip(constrain_id, 0, Ci));
		}
		//  todo:Fairness Condition at every vertex, take neighbor two polylines compute difference
		int NpolyLineConstrains = 0;
		const int ConstrainBias_Fairness = 5 * NumberOfFaces + (TotalValence + NumberOfVertices) + NumberOfFaces;
		for (int i = 0; i < NumberOfVertices && w_fairness>0; i++)
		{
			OpenMesh::SmartVertexHandle Vi(i, &mMesh);
			const Eigen::Vector3d  thisPoint = g_VertexCoordinate(i,X_c);

			int Valience = Vi.valence();
			if (Valience ==4)
			{

				OpenMesh::Mesh::VertexVertexIter n_vv_it = mMesh.vv_iter(Vi);
				const int nV0_id = (*n_vv_it).idx(); n_vv_it++;
				const Eigen::Vector3d nV0 = g_VertexCoordinate(nV0_id, X_c);   
				const int nV1_id = (*n_vv_it).idx(); n_vv_it++;
				const Eigen::Vector3d nV1 = g_VertexCoordinate(nV1_id, X_c);
				const int nV2_id = (*n_vv_it).idx(); n_vv_it++;
				const Eigen::Vector3d nV2 = g_VertexCoordinate(nV2_id, X_c); 
				const int nV3_id = (*n_vv_it).idx(); n_vv_it++;
				const Eigen::Vector3d nV3 = g_VertexCoordinate(nV3_id, X_c);
	

				//todo: vertical fairness
				const Eigen::Vector3d delta = nV3 + nV1 - 2 * thisPoint;
				int constrain_id = ConstrainBias_Fairness + NpolyLineConstrains;

				FillJacobi(constrain_id, 0, nV3_id, w_fairness * 2 * delta);
				FillJacobi(constrain_id, 0, nV1_id, w_fairness * 2 * delta);
				FillJacobi(constrain_id, 0, i, w_fairness * -4 * delta);
				//Constrain Vectors: (v3+v1-v)
				double fairness_v = w_fairness * (delta).dot(delta);
				C_tripletList.emplace_back(Trip(constrain_id, 0, fairness_v));




				//todo: horizontal fairness Constrain Vectors: (v2+v0-v)
				constrain_id +=1;
				const Eigen::Vector3d delta2 = nV2 + nV0 - 2 * thisPoint;
			
				FillJacobi(constrain_id, 0, nV2_id, w_fairness * 2 * delta2);
				FillJacobi(constrain_id, 0, nV0_id, w_fairness * 2 * delta2);
				FillJacobi(constrain_id, 0, i, w_fairness * -4 * delta2);
				//Constrain Vectors: (v3+v1-v)
				double fairness_h = w_fairness * (delta2).dot(delta2);
				C_tripletList.emplace_back(Trip(constrain_id, 0, fairness_h));

				NpolyLineConstrains += 2;
				}
		}
			

		//todo:torsion condition
		const int ConstrainBias_torsion = ConstrainBias_Fairness + NpolyLineConstrains;
		int constrain_id = ConstrainBias_Fairness; 
		for (int i = 0; i < NumberOfEdges && w_torsion>0; i++)
		{
			OpenMesh::SmartEdgeHandle eh(i, &mMesh);
			Eigen::Vector3d V0 =  g_VertexCoordinate((eh.v0().idx()), X_c);
			Eigen::Vector3d V1 = g_VertexCoordinate((eh.v1().idx()), X_c);
			Eigen::Vector3d N0 = g_VertexNormal((eh.v0().idx()), X_c);
			Eigen::Vector3d N1 = g_VertexNormal((eh.v1().idx()), X_c);
			Eigen::Vector3d ei = V0 - V1;

			Eigen::Vector3d edgeN = g_EdgeNormal(i, X_c);

			double Ci_0 = w_torsion * (edgeN.dot(N0));
			C_tripletList.emplace_back(Trip(constrain_id, 0, Ci_0));
			//jacobi to vn
			FillJacobi(constrain_id, 2, eh.v0().idx(), w_torsion* edgeN);
			//jacobi to EN
			FillJacobi(constrain_id, 3, i, w_torsion* N0);
		

			double Ci_1 = w_torsion * (edgeN.dot(N1));
			C_tripletList.emplace_back(Trip(constrain_id + 1, 0, Ci_1));
			//jacobi to vn
			FillJacobi(constrain_id+1, 2, eh.v1().idx(), w_torsion* edgeN);
			//jacobi to EN
			FillJacobi(constrain_id+1, 3, i, w_torsion* N1);

			double Ci_2 = w_torsion * (edgeN.dot(ei));
			C_tripletList.emplace_back(Trip(constrain_id + 2, 0, Ci_2));
			//jacobi to vn
			FillJacobi(constrain_id+2, 0, eh.v0().idx(), w_torsion* edgeN);
			FillJacobi(constrain_id + 2, 0, eh.v1().idx(), -w_torsion* edgeN);
			//jacobi to EN
			FillJacobi(constrain_id+2, 3, i, w_torsion* ei);


			double Ci_n = w_normal * (edgeN.dot(edgeN) - 1.0);
			C_tripletList.emplace_back(Trip(constrain_id + 3, 0, Ci_n));
			//jacobi to EN
			FillJacobi(constrain_id + 3, 3, i, w_normal* 2*edgeN);
			constrain_id += 4;
		}


		//todo: get total constrains
		NofConstrains = 5 * NumberOfFaces + (TotalValence + NumberOfVertices) + NumberOfFaces + NpolyLineConstrains +4*NumberOfEdges;


		return; 
	};




	Eigen::VectorXd X_c = X_o;
	std::vector<Trip> JacList;
	std::vector<Trip> CList;

	std::vector<Trip>i_trips;
	i_trips.reserve(NofUnknows);
	double  lambda = 1e-6;
	for (int i = 0; i < NofUnknows; i++)
	{
		i_trips.emplace_back(Trip(i, i, lambda));
	}
	SpMat SparseIdentity(NofUnknows, NofUnknows);
	SparseIdentity.setFromTriplets(i_trips.begin(), i_trips.end());
	while (iterationPerClick--)
	{
		JacList.clear();
		CList.clear();
		BuildJacobiMatandCVectors(X_c, JacList, CList);

		SpMat Jac(NofConstrains, NofUnknows);
		Jac.setFromTriplets(JacList.begin(), JacList.end());
		SpMat C(NofConstrains, 1);
		C.setFromTriplets(CList.begin(), CList.end());

		SpMat JTJ = (Jac.transpose() * Jac).pruned();
		
		//JTJ+lambda*I gurantee not singular
	
		SpMat JTJ_I = (JTJ + SparseIdentity).pruned();
		//for (int i = 0; i < NofUnknows; i++)
		//{
		//	JTJ.coeffRef(i, i) += lambda ;
		//}

		Eigen::MatrixXd b = (-1.0 * Jac.transpose() * C);

		//solving:

		Eigen::SimplicialLDLT<SpMat> solver;//// decomposition failed for SimplicialLDLT
		solver.compute(JTJ_I);
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
	}


	
	



	//update
	for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin();
		v_it != mMesh.vertices_end(); ++v_it)
	{
		Eigen::Vector3d Vx = { X_c(UnknowsAddressTrans(0,v_it->idx(),0),0),X_c(UnknowsAddressTrans(0,v_it->idx(),1),0) ,X_c(UnknowsAddressTrans(0,v_it->idx(),2),0) };
		Eigen::Vector3d NORMAL_x = { X_c(UnknowsAddressTrans(2,v_it->idx(),0),0),X_c(UnknowsAddressTrans(2,v_it->idx(),1),0) ,X_c(UnknowsAddressTrans(2,v_it->idx(),2),0) };
		mMesh.set_point(*v_it, Vx);
		mMesh.set_normal(*v_it, NORMAL_x);
	}


}


//todo:doesn't work
void Mesh::DiagReMeshing() {
	OpenMesh::Mesh mesh;
	
	std::vector<int>Noface_Vids;
	std::map<int, int>Old_vid2newVid;
	std::vector<OpenMesh::Mesh::VertexHandle>  face_vhandles;
	auto should_add = [](const std::vector<int>& Noface_Vids, int this_vid) ->bool {
		if (std::find(Noface_Vids.begin(), Noface_Vids.end(), this_vid) != Noface_Vids.end())
		{
			return false;
		}
		return true;
	};

	for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin(); v_it != mMesh.vertices_end(); ++v_it)
	{
		int valence = mMesh.valence(*v_it);
		
		if (valence >= 3 && should_add(Noface_Vids, v_it->idx()))
		{
		
			/*const Eigen::Vector3d thisV = mMesh.point(*v_it);*/

			face_vhandles.clear();
			for (OpenMesh::Mesh::VertexVertexCCWIter vv_it = mMesh.vv_ccwbegin(*v_it); vv_it.is_valid(); ++vv_it) {
				auto it = Old_vid2newVid.find((*vv_it).idx());
				if (it != Old_vid2newVid.end())
				{
					face_vhandles.push_back(OpenMesh::Mesh::VertexHandle(it->second));
				} 
				else
				{
					const Eigen::Vector3d thisV = mMesh.point(*vv_it);
					OpenMesh::Mesh::VertexHandle vhandle = mesh.add_vertex(thisV);
					Old_vid2newVid.insert(std::make_pair((*vv_it).idx(), vhandle.idx()));
	
					face_vhandles.push_back(vhandle);
				}	

				Noface_Vids.emplace_back((*vv_it).idx());
			}
		

			mesh.add_face(face_vhandles);
		}
	 
		
		
	}

	// request vertex normals, so the mesh reader can use normal information
	// if available
	mesh.request_vertex_normals();
	mesh.request_face_normals();
	// let the mesh update the normals
	mesh.update_normals();

	
	mMesh = mesh;
}
