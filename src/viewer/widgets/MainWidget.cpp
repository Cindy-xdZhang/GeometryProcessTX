#include "viewer/widgets/MainWidget.h"


MainWidget::MainWidget(bool expanded, LoaderPlugin* loader) : ViewerWidget(expanded)
{
    mAnchor = ViewerWidget::AnchorType::TopLeft;
    mLoader = loader;
}

ImVec4 MainWidget::draw(bool first, float scaling, float xSize, float ySize, float xPos, float yPos)
{
    ImGui::Begin(
            "Viewer", nullptr,
            ImGuiWindowFlags_NoSavedSettings
    );

    if (!first)
    {
        ImGui::SetWindowPos(ImVec2(xPos * scaling, yPos * scaling), ImGuiCond_Once);
    }
    ImGui::SetWindowCollapsed(mCollapsed, ImGuiCond_Once);

    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5f);

    // Viewing options
    if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Snap canonical view", ImVec2(-1, 0)))
        {
            mViewer->snap_to_canonical_quaternion();
        }

        // Zoom
        ImGui::DragFloat("Zoom", &(mViewer->core().camera_zoom), 0.05f, 0.1f, 20.0f);

        // Select rotation type
        int rotation_type = static_cast<int>(mViewer->core().rotation_type);
        static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
        static bool orthographic = true;
        if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0"))
        {
            using RT = ViewerCore::RotationType;
            auto new_type = static_cast<RT>(rotation_type);
            if (new_type != mViewer->core().rotation_type)
            {
                if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
                {
                    trackball_angle = mViewer->core().trackball_angle;
                    orthographic = mViewer->core().orthographic;
                    mViewer->core().trackball_angle = Eigen::Quaternionf::Identity();
                    mViewer->core().orthographic = true;
                }
                else if (mViewer->core().rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
                {
                    mViewer->core().trackball_angle = trackball_angle;
                    mViewer->core().orthographic = orthographic;
                }
                mViewer->core().set_rotation_type(new_type);
            }
        }

        // Orthographic view
        ImGui::Checkbox("Orthographic view", &(mViewer->core().orthographic));

        // Background
        ImGui::ColorEdit4("Background", mViewer->core().background_color.data(),
                          ImGuiColorEditFlags_InputRGB | ImGuiColorEditFlags_PickerHueWheel);
    }

    // Mesh
    if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Load mesh##Mesh", ImVec2(-1, 0)))
        {
            mViewer->open_dialog_load_mesh();
        }

        ImGui::Dummy(ImVec2(3.0f, 0.0f));
        ImGui::SameLine();

        if (mViewer->data_list.empty())
        {
            ImGui::Text("No mesh loaded!");
        }
        else
        {
            ImGui::Text("Loaded meshes:");
        }

        if (!mViewer->data_list.empty())
        {
            for (int i = 0; i < mViewer->data_list.size(); ++i)
            {
                if (mLoader->mesh(i))
                {
                    ImGui::PushID((int) i);
                    ImGui::Dummy(ImVec2(10.0f, 0.0f));
                    ImGui::SameLine();
                    bool active = mViewer->selected_data_index == i;
                    if (active)
                    {
                        ImVec4 color = ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive);
                        ImGui::PushStyleColor(ImGuiCol_Button, color);
                    }
                    if (ImGui::Button(mLoader->mesh(i)->name().c_str(), ImVec2(-1, 0)))
                    {
                        if (active)
                        {
                            mViewer->selected_data_index = -1;
                        }
                        else
                        {
                            mViewer->selected_data_index = i;
                        }
                    }
                    if (active)
                    {
                        ImGui::PopStyleColor();
                    }
                    ImGui::PopID();
                }
            }
        }
    }

    ImGui::SetWindowSize(ImVec2(250.0f * scaling, 0.0f));

    ImVec4 dim = ImVec4(ImGui::GetWindowPos().x, ImGui::GetWindowPos().y,
                        ImGui::GetWindowSize().x, ImGui::GetWindowSize().y);

    ImGui::PopItemWidth();

    ImGui::End();

    return dim;
}
