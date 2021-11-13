#ifndef VIEWER_CORE_H
#define VIEWER_CORE_H

#include "viewer/opengl/MeshGL.h"

#include <Eigen/Geometry>
#include <Eigen/Core>


class ViewerData;

class ViewerCore
{

public: // enums

    enum RotationType
    {
        ROTATION_TYPE_TRACKBALL = 0,
        ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
        ROTATION_TYPE_NO_ROTATION = 2,
        NUM_ROTATION_TYPES = 3
    };

public:

    ViewerCore();

    // Initialization
    void init();

    // Shutdown
    void shut();

    // ------------------- Camera control functions

    // Adjust the view to see the entire model
    void align_camera_center(
            const Eigen::MatrixXd& V);

    // Adjust the view to see the entire model
    void align_camera_center(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& F);

    // Determines how much to zoom and shift such that the mesh fills the unit
    // box (centered at the origin)
    static void get_scale_and_shift_to_fit_mesh(
            const Eigen::MatrixXd& V,
            float& zoom,
            Eigen::Vector3f& shift);

    // Determines how much to zoom and shift such that the mesh fills the unit
    // box (centered at the origin)
    static void get_scale_and_shift_to_fit_mesh(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& F,
            float& zoom,
            Eigen::Vector3f& shift);

    // ------------------- Drawing functions

    // Clear the frame buffers
    void clear_framebuffers();

    // Draw everything
    void draw(ViewerData& data, bool update_matrices = true);
    void draw_buffer(
            ViewerData& data,
            bool update_matrices,
            Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& R,
            Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& G,
            Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& B,
            Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& A);
    void draw_labels(
            ViewerData& data,
            const MeshGL::TextGL& labels
    );

    // Trackball angle (quaternion)
    void set_rotation_type(const RotationType& value);

    // ------------------- Option helpers

    // Set a ViewerData visualization option for this viewport
    void set(unsigned int& property_mask, bool value = true) const;

    // Unset a ViewerData visualization option for this viewport
    void unset(unsigned int& property_mask) const;

    // Toggle a ViewerData visualization option for this viewport
    void toggle(unsigned int& property_mask) const;

    // Check whether a ViewerData visualization option is set for this viewport
    [[nodiscard]] bool is_set(unsigned int property_mask) const;

public:

    bool is_visible = true;

    // Unique identifier
    unsigned int id = 1u;

    // Colors
    Eigen::Vector4f background_color;

    // Lighting
    Eigen::Vector3f light_position;
    float lighting_factor;

    RotationType rotation_type;
    Eigen::Quaternionf trackball_angle;

    // Camera parameters
    float camera_base_zoom;
    float camera_zoom;
    bool orthographic;
    Eigen::Vector3f camera_base_translation;
    Eigen::Vector3f camera_translation;
    Eigen::Vector3f camera_eye;
    Eigen::Vector3f camera_up;
    Eigen::Vector3f camera_center;
    float camera_view_angle;
    float camera_dnear;
    float camera_dfar;

    bool depth_test;

    // Animation
    bool is_animating;
    double animation_max_fps;

    // Caches the two-norm between the min/max point of the bounding box
    float object_scale;

    // Viewport size
    Eigen::Vector4f viewport;

    // Save the OpenGL transformation matrices used for the previous rendering pass
    Eigen::Matrix4f view;
    Eigen::Matrix4f proj;
    Eigen::Matrix4f norm;

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // VIEWER_CORE_H
