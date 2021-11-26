#include "viewer/widgets/MeshWidget.h"
#include <iostream>
#include "viewer/imgui/ImGuiHelpers.h"
#include "OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh"
#define mPI 3.141592653


MeshWidget::MeshWidget(bool expanded, LoaderPlugin* loader) : ViewerWidget(expanded)
{
    mAnchor = ViewerWidget::AnchorType::TopRight;
    mLoader = loader;
}

ImVec4 MeshWidget::draw(bool first, float scaling, float xSize, float ySize, float xPos, float yPos)
{
    ImGui::Begin(
            "Mesh", nullptr,
            ImGuiWindowFlags_NoSavedSettings
    );

    if (!first)
    {
        ImGui::SetWindowPos(ImVec2((xSize - xPos) * scaling, yPos * scaling), ImGuiCond_Once);
    }
    ImGui::SetWindowCollapsed(mCollapsed, ImGuiCond_Once);

    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5f);

    // Helper for setting viewport specific mesh options
    auto make_checkbox = [&](const char* label, unsigned int& option)
    {
        return ImGui::Checkbox(label,
                               [&]()
                               {
                                   return mViewer->core().is_set(option);
                               },
                               [&](bool value)
                               {
                                   return mViewer->core().set(option, value);
                               }
        );
    };

    if (mViewer->data_list.empty())
    {
        ImGui::Text("No mesh loaded!");
    }
    else if (mViewer->selected_data_index < 0)
    {
        ImGui::Text("No mesh selected!");
    }
    else
    {
        int index = mViewer->selected_data_index;

        // Reset
        if (ImGui::Button("Reset##Mesh", ImVec2(-1, 0)))
        {
            mLoader->mesh(index)->resetMesh();
            mLoader->mesh(index)->updateViewerData(mViewer->data(index));
        }

        // Save
        if (ImGui::Button("Save##Mesh", ImVec2(-1, 0)))
        {
            mViewer->open_dialog_save_mesh();
        }

        // Remove
        if (ImGui::Button("Remove##Mesh", ImVec2(-1, 0)))
        {
            mViewer->unload();
            goto end;
        }

        ImGui::Spacing();

        // Draw options
        if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Center view to mesh", ImVec2(-1, 0)))
            {
                mViewer->core().align_camera_center(mViewer->data().V, mViewer->data().F);
            }

            if (ImGui::Button("Normalize mesh", ImVec2(-1, 0)))
            {
                mLoader->mesh(index)->normalize();
                mLoader->mesh(index)->updateViewerData(mViewer->data(index));
            }

            make_checkbox("Visible", mViewer->data().is_visible);

            if (ImGui::Checkbox("Face-based", &(mViewer->data().face_based)))
            {
                mViewer->data().dirty = MeshGL::DIRTY_ALL;
            }
            if (ImGui::Checkbox("Render as flat triangles", &(mLoader->mesh(index)->renderFlatFaces())))
            {
                mLoader->mesh(index)->updateViewerData(mViewer->data(index));
            }

            make_checkbox("Show texture", mViewer->data().show_texture);

            if (ImGui::Checkbox("Invert normals", &(mViewer->data().invert_normals)))
            {
                mViewer->data().dirty |= MeshGL::DIRTY_NORMAL;
            }

            //make_checkbox("Show overlay", mViewer.data().show_overlay);
            //make_checkbox("Show overlay depth", mViewer.data().show_overlay_depth);

            if (ImGui::ColorEdit4("Line color", mLoader->mesh(index)->edgeColor().data(),
                                  ImGuiColorEditFlags_InputRGB | ImGuiColorEditFlags_PickerHueWheel))
            {
                mLoader->mesh(index)->updateViewerData(mViewer->data(index));
            }
            ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
            ImGui::DragFloat("Shininess", &(mViewer->data().shininess), 0.05f, 0.0f, 100.0f);
            ImGui::PopItemWidth();
        }

        ImGui::Spacing();

        // Noise
        if (ImGui::CollapsingHeader("Noise", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::DragScalar("sigma", ImGuiDataType_Double, &noiseStandardDeviation,
                              0.01f, &noiseStandardDeviationMin, &noiseStandardDeviationMax, "%.2f");

            ImGui::Combo("Direction", (int*) (&noiseDirection), "Normals\0Random\0\0");

            if (ImGui::Button("Add Gaussian noise", ImVec2(-1, 0)))
            {
                mLoader->mesh(index)->noise(noiseStandardDeviation, noiseDirection);
                mLoader->mesh(index)->updateViewerData(mViewer->data(index));
            }
        }

        ImGui::Spacing();





        // Overlays
        if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Checkbox("Wireframe", &(mLoader->mesh(index)->renderEdges())))
            {
                mLoader->mesh(index)->updateViewerData(mViewer->data(index));
            }

            make_checkbox("Fill", mViewer->data().show_faces);
            make_checkbox("Show vertex labels", mViewer->data().show_vertex_labels);
            make_checkbox("Show faces labels", mViewer->data().show_face_labels);
            make_checkbox("Show extra labels", mViewer->data().show_custom_labels);
        }

        ImGui::Spacing();


#ifdef PROJ1
        auto GetRotationMatrix = [](Eigen::Vector3d normalDir) {

            normalDir.normalize();

            Eigen::Vector3d tmp = Eigen::Vector3d::Ones();
            Eigen::Vector3d x_axis = normalDir.cross(tmp);
            //make sure they are linear independent
            while (x_axis.norm() == 0.0) {
                Eigen::Vector3d tmp = Eigen::Vector3d::Random();
                x_axis = normalDir.cross(tmp);
            }


            x_axis.normalize();
            auto y_axis = normalDir.cross(x_axis);
            y_axis.normalize();
            Eigen::Matrix3d mRotationI;
            mRotationI(0, 0) = x_axis(0);
            mRotationI(1, 0) = x_axis(1);
            mRotationI(2, 0) = x_axis(2);

            mRotationI(0, 1) = y_axis(0);
            mRotationI(1, 1) = y_axis(1);
            mRotationI(2, 1) = y_axis(2);

            mRotationI(0, 2) = normalDir(0);
            mRotationI(1, 2) = normalDir(1);
            mRotationI(2, 2) = normalDir(2);
            Eigen::Matrix3d RotationView = mRotationI.transpose();
            return RotationView;
        };
        //now only works for six neighbor points 
        auto FittingParaboloid = [](Eigen::MatrixXd XiYi, const int n_points) {


            Eigen::MatrixXd mA = Eigen::MatrixXd::Zero(6, 6);

            Eigen::MatrixXd mb = Eigen::MatrixXd::Zero(6, 1);

            for (int i = 0; i < n_points; i++) {
                const double  xi = XiYi(i, 0);
                const double  yi = XiYi(i, 1);
                const double  zi = XiYi(i, 2);
                const double  xi2 = XiYi(i, 0) * XiYi(i, 0);
                const double  yi2 = XiYi(i, 1) * XiYi(i, 1);
                const double  xiyi = XiYi(i, 0) * XiYi(i, 1);

                mA(0, 0) += 1.0;
                mA(0, 1) += xi;
                mA(0, 2) += yi;
                mA(0, 3) += xi2;
                mA(0, 4) += xiyi;
                mA(0, 5) += yi2;

                mA(1, 0) += xi;
                mA(1, 1) += xi2;
                mA(1, 2) += yi * xi;
                mA(1, 3) += xi2 * xi;
                mA(1, 4) += xiyi * xi;
                mA(1, 5) += yi2 * xi;

                mA(2, 0) += yi;
                mA(2, 1) += xiyi;
                mA(2, 2) += yi2;
                mA(2, 3) += xi2 * yi;
                mA(2, 4) += xiyi * yi;
                mA(2, 5) += yi2 * yi;


                mA(3, 0) += xi2;
                mA(3, 1) += xi * xi2;
                mA(3, 2) += yi * xi2;
                mA(3, 3) += xi2 * xi2;
                mA(3, 4) += xiyi * xi2;
                mA(3, 5) += yi2 * xi2;

                mA(4, 0) += xiyi;
                mA(4, 1) += xi * xiyi;
                mA(4, 2) += yi * xiyi;
                mA(4, 3) += xi2 * xiyi;
                mA(4, 4) += xiyi * xiyi;
                mA(4, 5) += yi2 * xiyi;

                mA(5, 0) += yi2;
                mA(5, 1) += xi * yi2;
                mA(5, 2) += yi * yi2;
                mA(5, 3) += xi2 * yi2;
                mA(5, 4) += xiyi * yi2;
                mA(5, 5) += yi2 * yi2;

                mb(0, 0) += zi;
                mb(1, 0) += zi * xi;
                mb(2, 0) += zi * yi;
                mb(3, 0) += zi * xi2;
                mb(4, 0) += zi * xiyi;
                mb(5, 0) += zi * yi2;

            }


            Eigen::MatrixXd mAI = mA.inverse();
            Eigen::MatrixXd x = mAI * mb;
            return x;
        };


        // Curvatures
        auto GaussianCurvatureComputation = [this, GetRotationMatrix, FittingParaboloid](const Mesh* inputMesh, const size_t SizeNeighbor,
            ViewerData& DataForDrawAsymptoticDir)-> Eigen::MatrixXd {
                //5 neighbors and 1 is the vertex of interest
                 //vertex of interest(origin) is used as the first point for fitting paraboloid
                const size_t numOfPoints = SizeNeighbor + 1;
                size_t numVertices = inputMesh->numVertices();

                /*		Eigen::MatrixXd P1 = DataForDrawAsymptoticDir.V;
                        Eigen::MatrixXd P2 = (DataForDrawAsymptoticDir.V + DataForDrawAsymptoticDir.V_normals).cast<double>();*/
                std::vector<Eigen::Vector3d>asympoticPoints;
                std::vector<Eigen::Vector3d>asympoticDir1;
                std::vector<Eigen::Vector3d>asympoticDir2;
                asympoticPoints.reserve(numVertices);
                asympoticDir1.reserve(numVertices);
                asympoticDir2.reserve(numVertices);


                auto HalfEdgeMesh = inputMesh->mesh();
                Eigen::MatrixXd   curvatures((int)numVertices, 1);


                for (auto v_it = HalfEdgeMesh.vertices_begin(); v_it != HalfEdgeMesh.vertices_end(); ++v_it) {


                    size_t vid = (*v_it).idx();
                    //OpenMesh::Mesh::Normal n = HalfEdgeMesh.normal(*v_it);
                    OpenMesh::Mesh::Normal n = HalfEdgeMesh.normal(*v_it);
                    const  Eigen::Matrix3d Rm = GetRotationMatrix(n);

                    Eigen::Vector3d thisV = HalfEdgeMesh.point(*v_it);

                    //collect neighborhoods of two round neighboors, otherwise will be two large
                    std::set<int> FirstRoundNeighborPointsId;
                    std::set<int> SecondRoundNeighborPointsId;

                    //this vertex should be ellimate from the traverse
                    FirstRoundNeighborPointsId.insert(vid);
                    SecondRoundNeighborPointsId.insert(vid);
                    for (OpenMesh::Mesh::VertexVertexIter vv_it = HalfEdgeMesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {

                        OpenMesh::ArrayKernel::VertexHandle neighborVertexHandle = *vv_it;
                        FirstRoundNeighborPointsId.insert(neighborVertexHandle.idx());
                        for (OpenMesh::Mesh::VertexVertexIter n_vv_it = HalfEdgeMesh.vv_iter(neighborVertexHandle); n_vv_it.is_valid(); ++n_vv_it) {

                            OpenMesh::ArrayKernel::VertexHandle SecondNearVHandle = (*n_vv_it);


                            SecondRoundNeighborPointsId.insert(SecondNearVHandle.idx());

                        }
                    }

                    std::set<int> mSecondRoundNeighborPointsId;
                    std::set_difference(SecondRoundNeighborPointsId.begin(), SecondRoundNeighborPointsId.end(),
                        FirstRoundNeighborPointsId.begin(), FirstRoundNeighborPointsId.end(), std::insert_iterator<std::set<int>>(mSecondRoundNeighborPointsId, mSecondRoundNeighborPointsId.begin()));

                    std::vector<Eigen::Vector3d>CollectedNeibors;
                    CollectedNeibors.reserve(numOfPoints);
                    //vertex of interest(oirgin) is used as the first point for fitting paraboloid


                    if (SecondRoundNeighborPointsId.size() < numOfPoints) {

                        std::cout << "W: Point has less neighbors than neighborSize is using random curvature.\n";
                        continue;
                    }
                    else if (FirstRoundNeighborPointsId.size() > numOfPoints) {

                        FirstRoundNeighborPointsId.erase(vid);
                        Eigen::Vector3d thisPoint = Eigen::Vector3d::Zero();
                        CollectedNeibors.emplace_back(thisPoint);

                        for (auto c : FirstRoundNeighborPointsId)
                        {
                            OpenMesh::ArrayKernel::VertexHandle vh(c);
                            Eigen::Vector3d toPoint = HalfEdgeMesh.point(vh);
                            //coordinates transform
                            toPoint = toPoint - thisV;
                            toPoint = Rm * toPoint;
                            CollectedNeibors.emplace_back(toPoint);
                        }

                    }
                    else {
                        //  first round<=neighborSize<=secondRoundNeighborPointsId
                        FirstRoundNeighborPointsId.erase(vid);
                        Eigen::Vector3d thisPoint = Eigen::Vector3d::Zero();
                        CollectedNeibors.emplace_back(thisPoint);

                        for (auto c : FirstRoundNeighborPointsId)
                        {
                            OpenMesh::ArrayKernel::VertexHandle vh(c);
                            Eigen::Vector3d toPoint = HalfEdgeMesh.point(vh);
                            //coordinates transform
                            toPoint = toPoint - thisV;
                            toPoint = Rm * toPoint;
                            CollectedNeibors.emplace_back(toPoint);
                        }
                        for (auto c : mSecondRoundNeighborPointsId)
                        {
                            OpenMesh::ArrayKernel::VertexHandle vh(c);
                            Eigen::Vector3d toPoint = HalfEdgeMesh.point(vh);
                            //coordinates transform
                            toPoint = toPoint - thisV;
                            toPoint = Rm * toPoint;
                            CollectedNeibors.emplace_back(toPoint);
                        }
                    }


                    //solve linear system for least-square regression
                    //CollectedNeibors.erase(CollectedNeibors.begin() + SizeNeighbor, CollectedNeibors.end());
                    Eigen::MatrixXd XIYI(numOfPoints, 3);
                    for (size_t i = 0; i < numOfPoints; ++i) {
                        XIYI.block(i, 0, 1, 3) << CollectedNeibors[i](0), CollectedNeibors[i](1), CollectedNeibors[i](2);
                    }


                    auto a = FittingParaboloid(XIYI, numOfPoints);
                    //compute Gaussian curvature
                    auto Gaussian_K = (4.0 * a(3) * a(5) - a(4) * a(4)) / pow((1 + a(1) * a(1) + a(2) * a(2)), 2);
                    auto H = 0.5 * ((1 + a(1) * a(1)) * 2.0 * a(5) - 2 * a(2) * a(4) * a(1) + (1 + a(2) * a(2)) \
                        * 2.0 * a(3)) / pow((1 + a(1) * a(1) + a(2) * a(2)), 1.5);
                    curvatures(vid, 0) = Gaussian_K;

                    //asympotic direction
                    if (Gaussian_K < 0.0) {
                        //double k1 = H + sqrt(H * H - Gaussian_K);
                        //double k2 = H - sqrt(H * H - Gaussian_K);
                        //double tmp = sqrt(-1.0*k1/k2);
                        //double theta = atan(tmp);
                        //double theta2 = atan(-1.0*tmp);
                        //Eigen::Vector3d thisAsymptoticDir,thisAsymptoticDir2;
                        ////local coordinate
                        //thisAsymptoticDir <<0.01* cos(theta), 0.01 * sin(theta), 0.0;
                        //thisAsymptoticDir2 << 0.01 * cos(theta2), 0.01 * sin(theta2), 0.0;
                        Eigen::Vector3d thisAsymptoticDir, thisAsymptoticDir2;
                        double delta = a(4) * a(4) - 4 * a(3) * a(5);
                        assert(delta >= 0);
                        auto ik_asym1 = (-1.0 * a(4) + sqrt(delta)) / (2.0 * a(3));
                        auto ik_asym2 = (-1.0 * a(4) - sqrt(delta)) / (2.0 * a(3));
                        thisAsymptoticDir << ik_asym1, 1.0, 0.0;
                        thisAsymptoticDir2 << ik_asym2, 1.0, 0.0;
                        thisAsymptoticDir.normalize();
                        thisAsymptoticDir2.normalize();

                        auto R = Rm.transpose();
                        thisAsymptoticDir = R * thisAsymptoticDir * 0.01;
                        thisAsymptoticDir += thisV;
                        thisAsymptoticDir2 = R * thisAsymptoticDir2 * 0.01;
                        thisAsymptoticDir2 += thisV;

                        asympoticPoints.emplace_back(thisV);//world coordinates

                        asympoticDir1.emplace_back(thisAsymptoticDir);
                        asympoticDir2.emplace_back(thisAsymptoticDir2);

                    }
                }
                int asympNums = asympoticPoints.size();
                Eigen::MatrixXd P1(asympNums, 3);
                Eigen::MatrixXd P2(asympNums, 3);
                Eigen::MatrixXd P3(asympNums, 3);
                Eigen::MatrixXd ColorRed(asympNums, 3);
                Eigen::MatrixXd ColorGreen(asympNums, 3);
                for (int i = 0; i < asympNums; i++) {
                    P1.block(i, 0, 1, 3) << asympoticPoints[i].transpose();
                    P2.block(i, 0, 1, 3) << asympoticDir1[i].transpose();
                    P3.block(i, 0, 1, 3) << asympoticDir2[i].transpose();
                    ColorRed.row(i) << 1, 0, 0;
                    ColorGreen.row(i) << 0, 1, 0;
                }



              /*  DataForDrawAsymptoticDir.add_edges(P1, P2, ColorRed);
                DataForDrawAsymptoticDir.add_edges(P1, P3, ColorGreen);
                DataForDrawAsymptoticDir.add_edges(P1, 2 * P1 - P2, ColorRed);
                DataForDrawAsymptoticDir.add_edges(P1, 2 * P1 - P3, ColorGreen);*/
                return curvatures;
        };
        // Curvatures
        auto MeanCurvatureComputation = [this, GetRotationMatrix, FittingParaboloid](const Mesh* inputMesh, const size_t SizeNeighbor)-> Eigen::MatrixXd {

            const size_t numOfPoints = SizeNeighbor + 1;
            size_t numVertices = inputMesh->numVertices();


            auto HalfEdgeMesh = inputMesh->mesh();
            Eigen::MatrixXd   curvatures((int)numVertices, 1);


            for (auto v_it = HalfEdgeMesh.vertices_begin(); v_it != HalfEdgeMesh.vertices_end(); ++v_it) {


                size_t vid = (*v_it).idx();
                //OpenMesh::Mesh::Normal n = HalfEdgeMesh.normal(*v_it);
                OpenMesh::Mesh::Normal n = HalfEdgeMesh.normal(*v_it);
                const  Eigen::Matrix3d Rm = GetRotationMatrix(n);

                Eigen::Vector3d thisV = HalfEdgeMesh.point(*v_it);

                //collect neighborhoods of two round neighboors, otherwise will be two large
                std::set<int> FirstRoundNeighborPointsId;
                std::set<int> SecondRoundNeighborPointsId;

                //this vertex should be ellimate from the traverse
                FirstRoundNeighborPointsId.insert(vid);
                SecondRoundNeighborPointsId.insert(vid);
                for (OpenMesh::Mesh::VertexVertexIter vv_it = HalfEdgeMesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {

                    OpenMesh::ArrayKernel::VertexHandle neighborVertexHandle = *vv_it;
                    FirstRoundNeighborPointsId.insert(neighborVertexHandle.idx());
                    for (OpenMesh::Mesh::VertexVertexIter n_vv_it = HalfEdgeMesh.vv_iter(neighborVertexHandle); n_vv_it.is_valid(); ++n_vv_it) {

                        OpenMesh::ArrayKernel::VertexHandle SecondNearVHandle = (*n_vv_it);


                        SecondRoundNeighborPointsId.insert(SecondNearVHandle.idx());

                    }
                }

                std::set<int> mSecondRoundNeighborPointsId;
                std::set_difference(SecondRoundNeighborPointsId.begin(), SecondRoundNeighborPointsId.end(),
                    FirstRoundNeighborPointsId.begin(), FirstRoundNeighborPointsId.end(), std::insert_iterator<std::set<int>>(mSecondRoundNeighborPointsId, mSecondRoundNeighborPointsId.begin()));

                std::vector<Eigen::Vector3d>CollectedNeibors;
                CollectedNeibors.reserve(numOfPoints);
                //vertex of interest(oirgin) is used as the first point for fitting paraboloid


                if (SecondRoundNeighborPointsId.size() < numOfPoints) {

                    std::cout << "W: Point has less neighbors than neighborSize is using random curvature.\n";
                    continue;
                }
                else if (FirstRoundNeighborPointsId.size() > numOfPoints) {

                    FirstRoundNeighborPointsId.erase(vid);
                    Eigen::Vector3d thisPoint = Eigen::Vector3d::Zero();
                    CollectedNeibors.emplace_back(thisPoint);

                    for (auto c : FirstRoundNeighborPointsId)
                    {
                        OpenMesh::ArrayKernel::VertexHandle vh(c);
                        Eigen::Vector3d toPoint = HalfEdgeMesh.point(vh);
                        //coordinates transform
                        toPoint = toPoint - thisV;
                        toPoint = Rm * toPoint;
                        CollectedNeibors.emplace_back(toPoint);
                    }

                }
                else {
                    //  first round<=neighborSize<=secondRoundNeighborPointsId
                    FirstRoundNeighborPointsId.erase(vid);
                    Eigen::Vector3d thisPoint = Eigen::Vector3d::Zero();
                    CollectedNeibors.emplace_back(thisPoint);

                    for (auto c : FirstRoundNeighborPointsId)
                    {
                        OpenMesh::ArrayKernel::VertexHandle vh(c);
                        Eigen::Vector3d toPoint = HalfEdgeMesh.point(vh);
                        //coordinates transform
                        toPoint = toPoint - thisV;
                        toPoint = Rm * toPoint;
                        CollectedNeibors.emplace_back(toPoint);
                    }
                    for (auto c : mSecondRoundNeighborPointsId)
                    {
                        OpenMesh::ArrayKernel::VertexHandle vh(c);
                        Eigen::Vector3d toPoint = HalfEdgeMesh.point(vh);
                        //coordinates transform
                        toPoint = toPoint - thisV;
                        toPoint = Rm * toPoint;
                        CollectedNeibors.emplace_back(toPoint);
                    }
                }


                //solve linear system for least-square regression
                //CollectedNeibors.erase(CollectedNeibors.begin() + SizeNeighbor, CollectedNeibors.end());
                Eigen::MatrixXd XIYI(numOfPoints, 3);
                for (size_t i = 0; i < numOfPoints; ++i) {
                    XIYI.block(i, 0, 1, 3) << CollectedNeibors[i](0), CollectedNeibors[i](1), CollectedNeibors[i](2);
                }


                auto a = FittingParaboloid(XIYI, numOfPoints);
                //compute Mean curvature
                auto Mean_K = 0.5 * ((1 + a(1) * a(1)) * 2.0 * a(5) - 2 * a(2) * a(4) * a(1) + (1 + a(2) * a(2)) * 2.0 * a(3)) / pow((1 + a(1) * a(1) + a(2) * a(2)), 1.5);



                curvatures(vid, 0) = Mean_K;
            }

            //Eigen::MatrixXd  NormalizedCurvatures = ((curvatures.array() - minCoeff) / ((maxCoeff - minCoeff) / 1.999)) - 1.0;
            return curvatures;

        };

        if (ImGui::CollapsingHeader("Project 1", ImGuiTreeNodeFlags_DefaultOpen))
        {
            static size_t sizeNeighbors = 5;
            size_t sizeNeighborsMin = 5;
            size_t sizeNeighborsMax = 20;
			ImGui::DragScalar("size of neighbors", ImGuiDataType_U32, &sizeNeighbors,
				1.0, &sizeNeighborsMin, &sizeNeighborsMax, "%d");

            if (ImGui::Button("Paint Gaussian Curvatures##Mesh", ImVec2(-1, 0)))
            {


                // TODO: compute curvatures here!
                Eigen::MatrixXd resCurvatures= GaussianCurvatureComputation(mLoader->mesh(index), sizeNeighbors, mViewer->data(index));
               
                const double caxis_min = resCurvatures.minCoeff();
				const double caxis_max = resCurvatures.maxCoeff();
				std::cout << "Max K= " << caxis_max << "Min K= " << caxis_min << std::endl;
                mViewer->data(index).set_data(resCurvatures);

            
              
            }

			if (ImGui::Button("Paint Mean Curvatures##Mesh", ImVec2(-1, 0)))
			{
				// TODO: compute curvatures here!
				Eigen::MatrixXd resCurvatures = MeanCurvatureComputation(mLoader->mesh(index), sizeNeighbors);
      
				mViewer->data(index).set_data(resCurvatures);
			}


        }

#endif
#ifdef PROJ3
        //collect K neigbors, when K <N1 return N1 neighbors; otherwise return K from N1 AND N2 neighbors 
        auto CollectNeigbors = [](const Mesh* inputMesh, const size_t SizeNeighbor, OpenMesh::Mesh::VertexIter v_it)->std::vector<int> {


            //const size_t numVertices = inputMesh->numVertices();

            std::vector<int>neighborsIDs;
            neighborsIDs.reserve(SizeNeighbor * 4);
            auto HalfEdgeMesh = inputMesh->mesh();
            const int vid = v_it->idx();

            //collect neighborhoods of two round neighboors, otherwise will be two large
            std::set<int> FirstRoundNeighborPointsId;
            std::set<int> SecondRoundNeighborPointsId;

            //this vertex should be ellimate from the traverse
            FirstRoundNeighborPointsId.insert(vid);
            SecondRoundNeighborPointsId.insert(vid);
            for (OpenMesh::Mesh::VertexVertexIter vv_it = HalfEdgeMesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {

                OpenMesh::ArrayKernel::VertexHandle neighborVertexHandle = *vv_it;
                FirstRoundNeighborPointsId.insert(neighborVertexHandle.idx());
                for (OpenMesh::Mesh::VertexVertexIter n_vv_it = HalfEdgeMesh.vv_iter(neighborVertexHandle); n_vv_it.is_valid(); ++n_vv_it) {

                    OpenMesh::ArrayKernel::VertexHandle SecondNearVHandle = (*n_vv_it);

                    SecondRoundNeighborPointsId.insert(SecondNearVHandle.idx());

                }
            }
            std::set<int> mSecondRoundNeighborPointsId;
            std::set_difference(SecondRoundNeighborPointsId.begin(), SecondRoundNeighborPointsId.end(),
                FirstRoundNeighborPointsId.begin(), FirstRoundNeighborPointsId.end(), std::insert_iterator<std::set<int>>(mSecondRoundNeighborPointsId, mSecondRoundNeighborPointsId.begin()));


            for (auto c : FirstRoundNeighborPointsId)
            {

                neighborsIDs.emplace_back(c);
            }

            if (FirstRoundNeighborPointsId.size() > SizeNeighbor) {

                /*FirstRoundNeighborPointsId.erase(vid);*/
                return neighborsIDs;

            }
            else {

                for (auto c : mSecondRoundNeighborPointsId)
                {
                    neighborsIDs.emplace_back(c);
                    if (neighborsIDs.size() >= SizeNeighbor)
                    {
                        break;
                    }
                   
                }
                return neighborsIDs;
            }

        };

        

		auto NormalCycleGaussianCurvature = [&CollectNeigbors](const Mesh* inputMesh, const size_t SizeNeighbor,
			ViewerData& DataForDrawAsymptoticDir)->Eigen::MatrixXd {
                auto HalfEdgeMesh = inputMesh->mesh();
                const int numVertices= inputMesh->numVertices();
                Eigen::MatrixXd  PerVertexK((int)numVertices, 1);

                for (auto v_it = HalfEdgeMesh.vertices_begin(); v_it != HalfEdgeMesh.vertices_end(); ++v_it) {
                    const int vid = v_it->idx();
                    const Eigen::Vector3d  thisPoint = HalfEdgeMesh.point(v_it);

                    OpenMesh::Mesh::VertexVertexIter n_vv_it = HalfEdgeMesh.vv_iter(v_it);
                    OpenMesh::Mesh::VertexVertexIter n_vv_it_last = n_vv_it;
                  
                    double sumAngles = 0.0;
                    Eigen::Vector3d FirstPoint = HalfEdgeMesh.point(*n_vv_it_last);

					for (n_vv_it++; n_vv_it.is_valid(); ++n_vv_it, ++n_vv_it_last) {
                        Eigen::Vector3d LastPoint = HalfEdgeMesh.point(*n_vv_it_last);
                        Eigen::Vector3d NextPoint = HalfEdgeMesh.point(*n_vv_it);

                        Eigen::Vector3d  Edge1 = LastPoint - thisPoint;
                        Eigen::Vector3d  Edge2 = NextPoint- thisPoint;
                        double cosEdge = Edge1.dot(Edge2) /(Edge1.norm()* Edge2.norm());
                        double Angle = acos(cosEdge);
                        sumAngles += Angle;
					}
                    Eigen::Vector3d FinalPoint = HalfEdgeMesh.point(*n_vv_it_last);
					Eigen::Vector3d  Edge1 = FirstPoint - thisPoint;
					Eigen::Vector3d  Edge2 = FinalPoint - thisPoint;
					double cosEdge = Edge1.dot(Edge2) / (Edge1.norm() * Edge2.norm());
					double Angle = acos(cosEdge);
					sumAngles += Angle;

                    PerVertexK(vid, 0) = 2 * mPI - sumAngles;

                }


                Eigen::MatrixXd  IntegratedVertexK((int)numVertices, 1);

				for (auto v_it = HalfEdgeMesh.vertices_begin(); v_it != HalfEdgeMesh.vertices_end(); ++v_it) {
					const int vid = v_it->idx();
					std::vector<int>  vns = CollectNeigbors(inputMesh, SizeNeighbor, v_it);
                    IntegratedVertexK(vid, 0) = 0;
					for (auto nvid : vns)
					{
						IntegratedVertexK(vid, 0) += PerVertexK(nvid, 0);
					}





				}






                return IntegratedVertexK;
        };

		auto NormalCycleMeanCurvature = [&CollectNeigbors](const Mesh* inputMesh, const size_t SizeNeighbor,
			ViewerData& DataForDrawAsymptoticDir)->Eigen::MatrixXd {
				auto HalfEdgeMesh = inputMesh->mesh();
				int NumEdges= HalfEdgeMesh.n_edges();
                const int numVertices = inputMesh->numVertices();
				Eigen::MatrixXd  PerEdgeM(NumEdges, 1);
                Eigen::MatrixXd  PerVertexM(numVertices, 1);
                Eigen::MatrixXd  IntegratedVertexM((int)numVertices, 1);

                for (OpenMesh::Mesh::EdgeIter e_it = HalfEdgeMesh.edges_begin(); e_it != HalfEdgeMesh.edges_end(); ++e_it) {

                    int eid = e_it->idx();

                    auto halfEdge1 = e_it->h0();
                    auto halfEdge2 = halfEdge1.opp();

                    auto point1 = halfEdge1.to();
                    auto point2 = halfEdge2.to();

                    auto next_halfedgePoint1= halfEdge1.next().to();
                    auto next_halfedgePoint2 = halfEdge2.next().to();
					
                    Eigen::Vector3d p1 = HalfEdgeMesh.point(point1);
					Eigen::Vector3d p2 = HalfEdgeMesh.point(point2);
					Eigen::Vector3d p3 = HalfEdgeMesh.point(next_halfedgePoint1);
					Eigen::Vector3d p4 = HalfEdgeMesh.point(next_halfedgePoint2);
					Eigen::Vector3d Vec_edge1 = p1 - p2;
                    Eigen::Vector3d Vec_edge2 = p3 - p1;
                    Eigen::Vector3d Vec_edge3 = p4 - p2;

					Eigen::Vector3d n1 = Vec_edge1.cross(Vec_edge2);

					Eigen::Vector3d n2 = Vec_edge3.cross(Vec_edge1);

                    double lengthEdge = Vec_edge1.norm();
					double cosEdge = n1.dot(n2) / (n1.norm() * n2.norm());
					double Angle = acos(cosEdge);
                    //TODO: Tell convex or concave
					Eigen::Vector3d vecN1 = HalfEdgeMesh.normal(point1);
					Eigen::Vector3d vecN2 = HalfEdgeMesh.normal(point2);
					Eigen::Vector3d vecP3plusP4 = 0.5 * (p3 + p4);
					Eigen::Vector3d vecM = 0.5*(p1+p2) - vecP3plusP4;
                    Eigen::Vector3d EdgeN = (vecN1 + vecN2) / 2;
                    double convex = EdgeN.dot(vecM);
                    double aflphae = 0;
                    if (convex<0)[[likely]]
                    {
                        aflphae= Angle;
                    }
                    else {
                        aflphae = -Angle;
                    }
                    PerEdgeM(eid, 0) = -0.5 * aflphae * lengthEdge;
                }


				for (auto v_it = HalfEdgeMesh.vertices_begin(); v_it != HalfEdgeMesh.vertices_end(); ++v_it) {

					const int vid = v_it->idx();
                    PerVertexM(vid, 0) = 0;
                    int valience = 0;
                    for (OpenMesh::Mesh::VertexEdgeIter ve_it= HalfEdgeMesh.ve_iter(*v_it); ve_it.is_valid(); ++ve_it)
                    {
                        int eid = ve_it->idx();
                        PerVertexM(vid, 0) += PerEdgeM(eid, 0);
                        valience++;
                    }

                    PerVertexM(vid, 0) /= (double)(valience);

				}



				for (auto v_it = HalfEdgeMesh.vertices_begin(); v_it != HalfEdgeMesh.vertices_end(); ++v_it) {
					const int vid = v_it->idx();
					std::vector<int>  vns = CollectNeigbors(inputMesh, SizeNeighbor, v_it);
                    IntegratedVertexM(vid, 0) = 0;
					for (auto nvid : vns)
					{
                        IntegratedVertexM(vid, 0) += PerVertexM(nvid, 0);
					}




				}



				return IntegratedVertexM;
		};


		if (ImGui::CollapsingHeader("Project 3", ImGuiTreeNodeFlags_DefaultOpen))
		{
			static size_t sizeNeighbors = 5;
			size_t sizeNeighborsMin = 1;
			size_t sizeNeighborsMax = 50;
            static size_t iterationPerClick = 5;
            static double slambda = 0.00001;
            static bool uniform = false;

            static double lambda = 0.1;
            static double mu = 0.01;
			static double gama = 0.001;
			static double theta = 0.000001;
         
            
            ImGui::Checkbox("Uniform", &uniform);
        
			ImGui::DragScalar("size of neighbors", ImGuiDataType_U32, &sizeNeighbors,
				1.0, &sizeNeighborsMin, &sizeNeighborsMax, "%d");
          
			ImGui::DragScalar("iterationPerClick ", ImGuiDataType_U32, &iterationPerClick,
				1.0, &sizeNeighborsMin, &sizeNeighborsMax, "%d");

           
            ImGui::InputDouble("lambda", &slambda,0.000001,0,"%.12f");


			ImGui::InputDouble("w-Ev0", &lambda, 0.1, 0, "%.8f");
			ImGui::InputDouble("w-Evc", &theta, 0.000001, 0, "%.8f");
			ImGui::InputDouble("w-Ef", &gama, 0.01, 0, "%.8f");
			ImGui::InputDouble("w-Et", &mu, 0.0001, 0, "%.8f");

			if (ImGui::Button("Discrete Gaussian Curvatures", ImVec2(-1, 0)))
			{


				// TODO: compute curvatures here!
				Eigen::MatrixXd resCurvatures = NormalCycleGaussianCurvature(mLoader->mesh(index), sizeNeighbors, mViewer->data(index));

				const double caxis_min = resCurvatures.minCoeff();
				const double caxis_max = resCurvatures.maxCoeff();
                std::cout << "Max K= "<< caxis_max <<"Min K= " << caxis_min << std::endl;
				mViewer->data(index).set_data(resCurvatures);



			}

			if (ImGui::Button("Discrete Mean Curvatures", ImVec2(-1, 0)))
			{
				// TODO: compute curvatures here!
				Eigen::MatrixXd resCurvatures = NormalCycleMeanCurvature(mLoader->mesh(index), sizeNeighbors, mViewer->data(index));

				const double caxis_min = resCurvatures.minCoeff();
				const double caxis_max = resCurvatures.maxCoeff();
				std::cout << "Max M= " << caxis_max << "Min M= " << caxis_min << std::endl;
				mViewer->data(index).set_data(resCurvatures);
			}

			if (ImGui::Button("Laplacian Smoothing", ImVec2(-1, 0)))
			{
				mLoader->mesh(index)->LaplacianSmoothing(slambda, iterationPerClick, uniform);
				mLoader->mesh(index)->updateViewerData(mViewer->data(index));

			}

			if (ImGui::Button("Energy Smoothing", ImVec2(-1, 0)))
			{
		        //double lambda, double mu, double gama, double theta, int iterations = 5, bool uniform = false
                int  iter = iterationPerClick;
                while (iter--)
                {
                    mLoader->mesh(index)->OptimizingSmoothing(lambda, mu, gama, theta, 1, uniform);
                }
				mLoader->mesh(index)->updateViewerData(mViewer->data(index));


			}



			if (ImGui::Button("render legend", ImVec2(-1, 0)))
			{
				// TODO: compute curvatures here!
                auto Npoints = mLoader->mesh(index)->mesh().n_vertices();
                Eigen::MatrixXd resCurvatures(14,1);
                const double interval = 4.0 / 6;
                for (int k=0;k<7;k++)
                {
                    int id1 = 2*k + 1;
                    int id2 = 2 * k + 2;
                    resCurvatures(id1-1, 0) = k * interval  -2.0;
                    resCurvatures(id2-1, 0) = k * interval  -2.0;
                }
               

				const double caxis_min = resCurvatures.minCoeff();
				const double caxis_max = resCurvatures.maxCoeff();
				std::cout << "Max = " << caxis_max << "Min = " << caxis_min << std::endl;
				mViewer->data(index).set_data(resCurvatures);
			}




		}

#endif
        if (ImGui::CollapsingHeader("Project 4", ImGuiTreeNodeFlags_DefaultOpen))
        { 
            size_t sizeNeighborsMin = 1;
	        size_t sizeNeighborsMax = 50;
            static size_t iterationPerClick = 5;

            static double w_c1= 0.1;
            static double w_c2 = 0.1;
            static double w_c3 = 0.1;

            static double w_normal = 0.01;
            static double w_fairness = 0.01;




            /*ImGui::Checkbox("Uniform", &uniform);

            ImGui::DragScalar("size of neighbors", ImGuiDataType_U32, &sizeNeighbors,
                1.0, &sizeNeighborsMin, &sizeNeighborsMax, "%d");
            */

			ImGui::DragScalar("iterationPerClick ", ImGuiDataType_U32, &iterationPerClick,
				1.0, &sizeNeighborsMin, &sizeNeighborsMax, "%d");


            ImGui::InputDouble("w-Con1", &w_c1, 0.0001, 0.1, "%.8f");
            ImGui::InputDouble("w-Con2", &w_c2, 0.0001, 0.1, "%.8f");
            ImGui::InputDouble("w-Con3", &w_c3, 0.0001, 0.1, "%.8f");
			ImGui::InputDouble("w-Normal", &w_normal, 0.0001, 0.1, "%.8f");
			ImGui::InputDouble("w-Fairness", &w_fairness, 0.0001, 0.1, "%.8f");
			
            if (ImGui::Button("Run Optimization", ImVec2(-1, 0)))
			{
                std::vector<double>Parameters = {(double)iterationPerClick,w_c1,w_c2,w_c3,w_normal,w_fairness};
                mLoader->mesh(index)->OptimizingQuadMesh(Parameters);
                mLoader->mesh(index)->updateViewerData(mViewer->data(index));
			}


        }



    }

	
end:

    ImGui::SetWindowSize(ImVec2(250.0f * scaling, 0.0f));

    ImVec4 dim = ImVec4(ImGui::GetWindowPos().x, ImGui::GetWindowPos().y,
                        ImGui::GetWindowSize().x, ImGui::GetWindowSize().y);

    ImGui::PopItemWidth();

    ImGui::End();

    return dim;
}
