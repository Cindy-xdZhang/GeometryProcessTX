#include "Mesh.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Eigen/sparse>
#include "viewer/ViewerData.h"

#include "utils/stl.h"
#include "utils/maths.h"
#include "utils/Eigen.h"
#include "utils/system.h"

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



	auto BuildLaplacianMatrixTriplest = [this,gama,NVertices,uniform]( ) {
		std::vector<Trip> tripletList;
		tripletList.reserve(3 * NVertices * NVertices);

		auto HalfEdgeMesh = this->mesh();

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


    auto BuildJacobiMatrix = [ NVertices,theta,gama,lambda,uniform, this]( const  std::vector<Trip>& LaplacianTripList) {
        std::vector<Trip> tripletList;
		tripletList.reserve(3 * NVertices * NVertices);
        auto HalfEdgeMesh = this->mesh();
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


	auto BuildXvector = [NVertices, this]()->Eigen::RowMatrixXd {
		Eigen::RowMatrixXd X_i(NVertices, 3);
		const auto HalfEdgeMesh = this->mesh();
		for (int i = 0; i < NVertices; i++)
		{
			Eigen::RowVector3d point = HalfEdgeMesh.point(OpenMesh::SmartVertexHandle(i, &HalfEdgeMesh));
			X_i.row(i) = point;
		}

		return X_i;
	};


	auto BuildCVector = [this,NVertices](const  std::vector<Trip>& LaplacianTripList, Eigen::RowMatrixXd Xo, Eigen::RowMatrixXd Xc) {
		
        Eigen::RowMatrixXd d(3 * NVertices, 3); 
        d *= 0;

		SpMat L(NVertices, NVertices);
		L.setFromTriplets(LaplacianTripList.begin(), LaplacianTripList.end());
		Eigen::RowMatrixXd LX = L * Xc;

		d.block(0, 0, NVertices, 3) = LX.block(0, 0, NVertices, 3);



		for (int i = NVertices; i < 2 * NVertices; i++)
		{
            d.row(i) = Xc.row(i- NVertices) - Xo.row(i - NVertices);
		}

		/*	for (int i = 2 * NVertices; i < 3 * NVertices; i++)
			{
			d.row(i)={0,0,0};
			}*/

		return d;
	};

	
    std::vector<Trip> L_trips= BuildLaplacianMatrixTriplest();
    SpMat Jac = BuildJacobiMatrix(L_trips);
    Eigen::RowMatrixXd X_o = BuildXvector();
    Eigen::RowMatrixXd X_c = X_o;

    Eigen::RowMatrixXd C = BuildCVector(L_trips,X_o,X_c);

    //solving:
	SpMat A =  (Jac.transpose() * Jac).pruned();

	Eigen::SimplicialLDLT<SpMat> solver(A);
	Eigen::RowMatrixXd b = -1.0 * Jac.transpose() *  C ;
	Eigen::RowMatrixXd Delta_V= solver.solve(b);
    X_c = Delta_V+ X_c;

	for (OpenMesh::Mesh::VertexIter v_it = mMesh.vertices_begin();
		v_it != mMesh.vertices_end(); ++v_it)
	{
		OpenMesh::Mesh::Point p = X_c.row(v_it->idx());
		mMesh.set_point(*v_it, p);
	}
}