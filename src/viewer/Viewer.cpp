#include "viewer/Viewer.h"

#include <cmath>
#include <cstdio>
#include <chrono>
#include <thread>
#include <cassert>
#include <iostream>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "utils/system.h"
#include "utils/camera.h"


// Internal global variables used for glfw event handling
static Viewer* __viewer;
static int highdpi = 1;
static double scroll_x = 0;
static double scroll_y = 0;


static void glfw_mouse_press(GLFWwindow* /*window*/, int button, int action, int modifier)
{
    Viewer::MouseButton mb;

    if (button == GLFW_MOUSE_BUTTON_1)
    {
        mb = Viewer::MouseButton::Left;
    }
    else if (button == GLFW_MOUSE_BUTTON_2)
    {
        mb = Viewer::MouseButton::Right;
    }
    else // if (button == GLFW_MOUSE_BUTTON_3)
    {
        mb = Viewer::MouseButton::Middle;
    }

    if (action == GLFW_PRESS)
    {
        __viewer->mouse_down(mb, modifier);
    }
    else
    {
        __viewer->mouse_up(mb, modifier);
    }
}

static void glfw_error_callback(int /*error*/, const char* description)
{
    fputs(description, stderr);
}

static void glfw_char_mods_callback(GLFWwindow* /*window*/, unsigned int codepoint, int modifier)
{
    __viewer->key_pressed(codepoint, modifier);
}

static void glfw_key_callback(GLFWwindow* window, int key, int /*scancode*/, int action, int modifier)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }

    if (action == GLFW_PRESS)
    {
        __viewer->key_down(key, modifier);
    }
    else if (action == GLFW_RELEASE)
    {
        __viewer->key_up(key, modifier);
    }
    else if (action == GLFW_REPEAT)
    {
        __viewer->key_repeat(key, modifier);
    }
}

static void glfw_window_size(GLFWwindow* window, int width, int height)
{
    int w = (int) (width * highdpi);
    int h = (int) (height * highdpi);

    __viewer->post_resize(w, h);
}

static void glfw_mouse_move(GLFWwindow* /*window*/, double x, double y)
{
    __viewer->mouse_move((int) x * highdpi, (int) y * highdpi);
}

static void glfw_mouse_scroll(GLFWwindow* /*window*/, double x, double y)
{
    using namespace std;
    scroll_x += x;
    scroll_y += y;

    __viewer->mouse_scroll((float) y);
}

static void glfw_drop_callback(GLFWwindow* /*window*/, int /*count*/, const char** /*filenames*/)
{

}

[[maybe_unused]] int Viewer::launch(bool resizable /*= true*/, bool fullscreen /*= false*/, bool maximize /*= false*/,
                                    const std::string& name, int windowWidth /*= 0*/, int windowHeight /*= 0*/)
{
    if (launch_init(resizable, fullscreen, maximize, name, windowWidth, windowHeight) == EXIT_FAILURE)
    {
        launch_shut();
        return EXIT_FAILURE;
    }
    launch_rendering(true);
    launch_shut();
    return EXIT_SUCCESS;
}

int Viewer::launch_init(bool resizable, bool fullscreen, bool maximize,
                        const std::string& name, int windowWidth, int windowHeight)
{
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
    {
        fprintf(stderr, "Error: Could not initialize OpenGL context");
        return EXIT_FAILURE;
    }
    glfwWindowHint(GLFW_SAMPLES, 8);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    if (fullscreen)
    {
        GLFWmonitor* monitor = glfwGetPrimaryMonitor();
        const GLFWvidmode* mode = glfwGetVideoMode(monitor);
        window = glfwCreateWindow(mode->width, mode->height, name.c_str(), monitor, nullptr);
        windowWidth = mode->width;
        windowHeight = mode->height;
    }
    else
    {
        // Set default windows width
        if (windowWidth <= 0 && core_list.size() == 1 && core().viewport[2] > 0)
        {
            windowWidth = (int) core().viewport[2];
        }
        else if (windowWidth <= 0)
        {
            windowWidth = 1280;
        }
        // Set default windows height
        if (windowHeight <= 0 && core_list.size() == 1 && core().viewport[3] > 0)
        {
            windowHeight = (int) core().viewport[3];
        }
        else if (windowHeight <= 0)
        {
            windowHeight = 800;
        }
        window = glfwCreateWindow(windowWidth, windowHeight, name.c_str(), nullptr, nullptr);
        if (maximize)
        {
            glfwMaximizeWindow(window);
        }
    }
    if (!window)
    {
        fprintf(stderr, "Error: Could not create GLFW window");
        glfwTerminate();
        return EXIT_FAILURE;
    }
    glfwMakeContextCurrent(window);

    // Load OpenGL and its extensions
    if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
    {
        printf("Failed to load OpenGL and its extensions\n");
        return (-1);
    }
#if defined(DEBUG) || defined(_DEBUG)
    printf("OpenGL Version %d.%d loaded\n", GLVersion.major, GLVersion.minor);
    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version received: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
#endif
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    // Initialize one and only reference
    __viewer = this;

    // Register callbacks
    glfwSetKeyCallback(window, glfw_key_callback);
    glfwSetCursorPosCallback(window, glfw_mouse_move);
    glfwSetWindowSizeCallback(window, glfw_window_size);
    glfwSetMouseButtonCallback(window, glfw_mouse_press);
    glfwSetScrollCallback(window, glfw_mouse_scroll);
    glfwSetCharModsCallback(window, glfw_char_mods_callback);
    glfwSetDropCallback(window, glfw_drop_callback);

    // Handle retina displays (windows and mac)
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    int width_window, height_window;
    glfwGetWindowSize(window, &width_window, &height_window);
    highdpi = windowWidth / width_window;
    glfw_window_size(window, width_window, height_window);

    // Initialize viewer
    init();
    for (auto& core: this->core_list)
    {
        for (auto& data: this->data_list)
        {
            if (data.is_visible & core.id)
            {
                this->core(core.id).align_camera_center(data.V, data.F);
            }
        }
    }

    return EXIT_SUCCESS;
}

void Viewer::launch_rendering(bool loop)
{
    // glfwMakeContextCurrent(window);
    // Rendering loop
    bool first = true;
    const int num_extra_frames = 5;
    int frame_counter = 0;
    while (!glfwWindowShouldClose(window))
    {
        double tic = std::chrono::duration<double>(
                std::chrono::system_clock::now().time_since_epoch()).count();
        draw(first);
        if (first)
        {
            first = false;
        }
        glfwSwapBuffers(window);
        if (core().is_animating || frame_counter++ < num_extra_frames)
        {
            glfwPollEvents();
            // In microseconds
            double toc = std::chrono::duration<double>(
                    std::chrono::system_clock::now().time_since_epoch()).count();
            double duration = 1000000. * (toc - tic);
            const double min_duration = 1000000. / core().animation_max_fps;
            if (duration < min_duration)
            {
                std::this_thread::sleep_for(std::chrono::microseconds((int) (min_duration - duration)));
            }
        }
        else
        {
            glfwWaitEvents();
            frame_counter = 0;
        }
        if (!loop)
        {
            return;
        }

#ifdef __APPLE__
        static bool first_time_hack  = true;
        if(first_time_hack) {
          glfwHideWindow(window);
          glfwShowWindow(window);
          first_time_hack = false;
        }
#endif
    }
}

void Viewer::launch_shut()
{
    for (auto& data: data_list)
    {
        data.meshgl.free();
    }
    core().shut(); // Doesn't do anything
    shutdown_plugins();
    glfwDestroyWindow(window);
    glfwTerminate();
}

void Viewer::init()
{
    core().init(); // Doesn't do anything

    if (callback_init)
    {
        if (callback_init(*this))
        {
            return;
        }
    }

    init_plugins();
}

void Viewer::init_plugins()
{
    // Init all plugins
    for (auto& plugin: plugins)
    {
        plugin->init(this);
    }
}

void Viewer::shutdown_plugins()
{
    for (auto& plugin: plugins)
    {
        plugin->shutdown();
    }
}

Viewer::Viewer()
{
    window = nullptr;

    selected_data_index = -1;
    next_data_id = 1;

    core_list.emplace_back(ViewerCore());
    selected_core_index = 0;
    core_list.front().id = 1;
    next_core_id = 2;

    // Temporary variables initialization
    down = false;
    hack_never_moved = true;
    scroll_position = 0.0f;

    // C-style callbacks
    callback_init = nullptr;
    callback_pre_draw = nullptr;
    callback_post_draw = nullptr;
    callback_mouse_down = nullptr;
    callback_mouse_up = nullptr;
    callback_mouse_move = nullptr;
    callback_mouse_scroll = nullptr;
    callback_key_pressed = nullptr;
    callback_key_down = nullptr;
    callback_key_up = nullptr;
    callback_key_repeat = nullptr;

    callback_init_data = nullptr;
    callback_pre_draw_data = nullptr;
    callback_post_draw_data = nullptr;
    callback_mouse_down_data = nullptr;
    callback_mouse_up_data = nullptr;
    callback_mouse_move_data = nullptr;
    callback_mouse_scroll_data = nullptr;
    callback_key_pressed_data = nullptr;
    callback_key_down_data = nullptr;
    callback_key_up_data = nullptr;
    callback_key_repeat_data = nullptr;
}

Viewer::~Viewer() = default;

////////////////////////////////////////////////////////////////////////////////
// Mesh IO
////////////////////////////////////////////////////////////////////////////////

bool Viewer::load(const std::string& filename, bool camera)
{
    bool loaded = false;

    for (auto& plugin: plugins)
    {
        if ((loaded = plugin->load(filename)))
        {
            break;
        }
    }

    if (!loaded)
    {
        return false;
    }

    if (camera)
    {
        for (auto& core: core_list)
        {
            core.align_camera_center(data().V, data().F);
        }
    }

    for (auto& plugin: plugins)
    {
        if (plugin->post_load())
        {
            return true;
        }
    }

    return false;
}

bool Viewer::save(const std::string& filename)
{
    for (auto& plugin: plugins)
    {
        if (plugin->save(filename))
        {
            return true;
        }
    }

    return false;
}

bool Viewer::unload()
{
    for (auto& plugin: plugins)
    {
        if (plugin->unload())
        {
            return true;
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////
// Screenshots
////////////////////////////////////////////////////////////////////////////////

int Viewer::screenshot(const char* filename) const
{
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    GLsizei nrChannels = 3;
    GLsizei stride = nrChannels * width;
    stride += (stride % 4) ? (4 - stride % 4) : 0;

    GLsizei bufferSize = stride * height;
    std::vector<char> buffer(bufferSize);

    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());

    stbi_flip_vertically_on_write(true);
    return stbi_write_png(filename, width, height, nrChannels, buffer.data(), stride);
}

int Viewer::screenshot(const std::string& filename) const
{
    return screenshot(filename.c_str());
}

////////////////////////////////////////////////////////////////////////////////
// Callbacks
////////////////////////////////////////////////////////////////////////////////

bool Viewer::key_pressed(unsigned int key, int modifiers)
{
    for (auto& plugin: plugins)
    {
        if (plugin->key_pressed(key, modifiers))
        {
            return true;
        }
    }

    if (callback_key_pressed)
    {
        if (callback_key_pressed(*this, key, modifiers))
        {
            return true;
        }
    }

    return false;
}

bool Viewer::key_down(int key, int modifiers)
{
    for (auto& plugin: plugins)
    {
        if (plugin->key_down(key, modifiers))
        {
            return true;
        }
    }

    if (callback_key_down)
    {
        if (callback_key_down(*this, key, modifiers))
        {
            return true;
        }
    }

    return false;
}

bool Viewer::key_up(int key, int modifiers)
{
    for (auto& plugin: plugins)
    {
        if (plugin->key_up(key, modifiers))
        {
            return true;
        }
    }

    if (callback_key_up)
    {
        if (callback_key_up(*this, key, modifiers))
        {
            return true;
        }
    }

    return false;
}

bool Viewer::key_repeat(int key, int modifiers)
{
    for (auto& plugin: plugins)
    {
        if (plugin->key_repeat(key, modifiers))
        {
            return true;
        }
    }

    if (callback_key_repeat)
    {
        if (callback_key_repeat(*this, key, modifiers))
        {
            return true;
        }
    }

    return false;
}

void Viewer::select_hovered_core()
{
    int width_window, height_window;
    glfwGetFramebufferSize(window, &width_window, &height_window);
    for (int i = 0; i < core_list.size(); i++)
    {
        Eigen::Vector4f viewport = core_list[i].viewport;

        if (((float) current_mouse_x > viewport[0]) &&
            ((float) current_mouse_x < viewport[0] + viewport[2]) &&
            ((float) (height_window - current_mouse_y) > viewport[1]) &&
            ((float) (height_window - current_mouse_y) < viewport[1] + viewport[3]))
        {
            selected_core_index = i;
            break;
        }
    }
}

bool Viewer::mouse_down(MouseButton button, int modifier)
{
    // Remember mouse location at down even if used by callback/plugin
    down_mouse_x = current_mouse_x;
    down_mouse_y = current_mouse_y;

    for (auto& plugin: plugins)
    {
        if (plugin->mouse_down(static_cast<int>(button), modifier))
        {
            return true;
        }
    }

    if (callback_mouse_down)
    {
        if (callback_mouse_down(*this, static_cast<int>(button), modifier))
        {
            return true;
        }
    }

    down = true;

    // Select the core containing the click location.
    select_hovered_core();

    down_translation = core().camera_translation;

    // Initialization code for the trackball
    Eigen::RowVector3d center;
    if (selected_data_index < 0 || data().V.rows() == 0)
    {
        center << 0, 0, 0;
    }
    else
    {
        center = data().V.colwise().sum() / data().V.rows();
    }

    Eigen::Vector3f coord =
            utils::project(
                    Eigen::Vector3f(center(0), center(1), center(2)),
                    core().view,
                    core().proj,
                    core().viewport);
    down_mouse_z = coord[2];
    down_rotation = core().trackball_angle;

    mouse_mode = MouseMode::Rotation;

    switch (button)
    {
        case MouseButton::Left:
            if (core().rotation_type == ViewerCore::ROTATION_TYPE_NO_ROTATION)
            {
                mouse_mode = MouseMode::Translation;
            }
            else
            {
                mouse_mode = MouseMode::Rotation;
            }
            break;

        case MouseButton::Right:
            mouse_mode = MouseMode::Translation;
            break;

        default:
            mouse_mode = MouseMode::None;
            break;
    }

    return true;
}

bool Viewer::mouse_up(MouseButton button, int modifier)
{
    down = false;

    for (auto& plugin: plugins)
    {
        if (plugin->mouse_up(static_cast<int>(button), modifier))
        {
            return true;
        }
    }

    if (callback_mouse_up)
    {
        if (callback_mouse_up(*this, static_cast<int>(button), modifier))
        {
            return true;
        }
    }

    mouse_mode = MouseMode::None;

    return true;
}

bool Viewer::mouse_move(int mouse_x, int mouse_y)
{
    if (hack_never_moved)
    {
        down_mouse_x = mouse_x;
        down_mouse_y = mouse_y;
        hack_never_moved = false;
    }
    current_mouse_x = mouse_x;
    current_mouse_y = mouse_y;

    for (auto& plugin: plugins)
    {
        if (plugin->mouse_move(mouse_x, mouse_y))
        {
            return true;
        }
    }

    if (callback_mouse_move)
    {
        if (callback_mouse_move(*this, mouse_x, mouse_y))
        {
            return true;
        }
    }

    if (down)
    {
        // We need the window height to transform the mouse click coordinates
        // into viewport-mouse-click coordinates for trackball and
        // two_axis_valuator_fixed_up
        int width_window, height_window;
        glfwGetFramebufferSize(window, &width_window, &height_window);
        switch (mouse_mode)
        {
            case MouseMode::Rotation:
            {
                switch (core().rotation_type)
                {
                    default:
                        assert(false && "Unknown rotation type");
                    case ViewerCore::ROTATION_TYPE_NO_ROTATION:
                        break;
                    case ViewerCore::ROTATION_TYPE_TRACKBALL:
                        utils::trackball(
                                core().viewport(2),
                                core().viewport(3),
                                2.0,
                                down_rotation,
                                (float) down_mouse_x - core().viewport(0),
                                (float) down_mouse_y - (height_window - core().viewport(1) - core().viewport(3)),
                                (float) mouse_x - core().viewport(0),
                                (float) mouse_y - (height_window - core().viewport(1) - core().viewport(3)),
                                core().trackball_angle);
                        break;
                    case ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP:
                        utils::two_axis_evaluator_fixed_up(
                                core().viewport(2),
                                core().viewport(3),
                                2.0,
                                down_rotation,
                                (float) down_mouse_x - core().viewport(0),
                                (float) down_mouse_y - (height_window - core().viewport(1) - core().viewport(3)),
                                (float) mouse_x - core().viewport(0),
                                (float) mouse_y - (height_window - core().viewport(1) - core().viewport(3)),
                                core().trackball_angle);
                        break;
                }
                //Eigen::Vector4f snapq = core().trackball_angle;

                break;
            }

            case MouseMode::Translation:
            {
                //translation
                Eigen::Vector3f pos1 = utils::unproject(
                        Eigen::Vector3f((float) mouse_x, core().viewport[3] - mouse_y, down_mouse_z),
                        core().view, core().proj, core().viewport);
                Eigen::Vector3f pos0 = utils::unproject(
                        Eigen::Vector3f((float) down_mouse_x, core().viewport[3] - down_mouse_y, down_mouse_z),
                        core().view, core().proj, core().viewport);

                Eigen::Vector3f diff = pos1 - pos0;
                core().camera_translation = down_translation + Eigen::Vector3f(diff[0], diff[1], diff[2]);

                break;
            }
            case MouseMode::Zoom:
            {
                float delta = 0.001f * (mouse_x - down_mouse_x + mouse_y - down_mouse_y);
                core().camera_zoom *= 1 + delta;
                down_mouse_x = mouse_x;
                down_mouse_y = mouse_y;
                break;
            }

            default:
                break;
        }
    }
    return true;
}

bool Viewer::mouse_scroll(float delta_y)
{
    // Direct the scrolling operation to the appropriate viewport
    // (unless the core selection is locked by an ongoing mouse interaction).
    if (!down)
    {
        select_hovered_core();
    }
    scroll_position += delta_y;

    for (auto& plugin: plugins)
    {
        if (plugin->mouse_scroll(delta_y))
        {
            return true;
        }
    }

    if (callback_mouse_scroll)
    {
        if (callback_mouse_scroll(*this, delta_y))
        {
            return true;
        }
    }

    // Only zoom if there's actually a change
    if (delta_y != 0)
    {
        float mult = (1.0f + ((delta_y > 0) ? 1.f : -1.f) * 0.05f);
        const float min_zoom = 0.1f;
        core().camera_zoom = (core().camera_zoom * mult > min_zoom ? core().camera_zoom * mult : min_zoom);
    }
    return true;
}

void Viewer::draw(bool first)
{
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    int width_window, height_window;
    glfwGetWindowSize(window, &width_window, &height_window);

    auto highdpi_tmp = (width_window == 0 || width == 0) ? highdpi : (width / width_window);

    if (fabs(highdpi_tmp - highdpi) > 1e-8)
    {
        post_resize(width, height);
        highdpi = highdpi_tmp;
    }

    for (auto& core: core_list)
    {
        core.clear_framebuffers();
    }

    for (auto& plugin: plugins)
    {
        if (plugin->pre_draw(first))
        {
            break;
        }
    }

    if (callback_pre_draw)
    {
        if (callback_pre_draw(*this))
        {

        }
    }

    refresh();

    for (auto& plugin: plugins)
    {
        if (plugin->post_draw(first))
        {
            break;
        }
    }

    if (callback_post_draw)
    {
        if (callback_post_draw(*this))
        {

        }
    }
}

void Viewer::refresh()
{
    for (auto& core: core_list)
    {
        if (core.is_visible)
        {
            for (auto& data: data_list)
            {
                if (data.is_visible & core.id)
                {
                    core.draw(data);
                }
            }
        }
    }
}

void Viewer::resize(int w, int h)
{
    if (window)
    {
        glfwSetWindowSize(window, (int) (w / highdpi), (int) (h / highdpi));
    }
    post_resize(w, h);
}

void Viewer::post_resize(int w, int h)
{
    if (core_list.size() == 1)
    {
        core().viewport = Eigen::Vector4f(0, 0, w, h);
    }
    else
    {
        // It is up to the user to define the behavior of the post_resize() function
        // when there are multiple viewports (through the `callback_post_resize` callback)
    }
    for (auto& plugin: plugins)
    {
        plugin->post_resize(w, h);
    }
    if (callback_post_resize)
    {
        callback_post_resize(*this, w, h);
    }
}

void Viewer::snap_to_canonical_quaternion()
{
    Eigen::Quaternionf snapq = this->core().trackball_angle;
    utils::snap_to_canonical_view_quat(snapq, 1.0f, this->core().trackball_angle);
}

bool Viewer::open_dialog_load_mesh()
{
    std::string fname = utils::file_dialog_open(path);

    if (fname.length() == 0)
    {
        return false;
    }

    return this->load(fname);
}

bool Viewer::open_dialog_save_mesh()
{
    std::string fname = utils::file_dialog_save(path);

    if (fname.length() == 0)
    {
        return false;
    }

    return this->save(fname);
}

////////////////////////////////////////////////////////////////////////////////
/// Multi-mesh methods
////////////////////////////////////////////////////////////////////////////////

ViewerData& Viewer::data(int data_id /*= -1*/)
{
    int index;
    if (data_id == -1)
    {
        index = selected_data_index;
    }
    else
    {
        index = data_id;
    }

    assert((index >= 0 && index < data_list.size()) &&
           "selected_data_index or mesh_id should be in bounds");
    return data_list[index];
}

const ViewerData& Viewer::data(int data_id /*= -1*/) const
{
    int index;
    if (data_id == -1)
    {
        index = selected_data_index;
    }
    else
    {
        index = data_id;
    }

    assert((index >= 0 && index < data_list.size()) &&
           "selected_data_index or mesh_id should be in bounds");
    return data_list[index];
}

int Viewer::append_data(bool visible /*= true*/)
{
    data_list.emplace_back();
    selected_data_index = (int) data_list.size() - 1;
    data_list.back().id = next_data_id++;
    if (visible)
    {
        for (auto& core: core_list)
        {
            data_list.back().set_visible(true, core.id);
        }
    }
    else
    {
        data_list.back().is_visible = 0;
    }
    return data_list.back().id;
}

bool Viewer::erase_data(size_t index)
{
    assert((index >= 0 && index < data_list.size()) && "index should be in bounds");

    data_list[index].meshgl.free();
    data_list.erase(data_list.begin() + index);
    if (selected_data_index >= index)
    {
        selected_data_index--;
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Multi-viewport methods
////////////////////////////////////////////////////////////////////////////////

ViewerCore& Viewer::core(unsigned core_id /*= 0*/)
{
    assert(!core_list.empty() && "core_list should never be empty");
    int core_index;
    if (core_id == 0)
    {
        core_index = selected_core_index;
    }
    else
    {
        core_index = this->core_index((int) core_id);
    }
    assert((core_index >= 0 && core_index < core_list.size()) && "selected_core_index should be in bounds");
    return core_list[core_index];
}

const ViewerCore& Viewer::core(unsigned core_id /*= 0*/) const
{
    assert(!core_list.empty() && "core_list should never be empty");
    int core_index;
    if (core_id == 0)
    {
        core_index = selected_core_index;
    }
    else
    {
        core_index = this->core_index((int) core_id);
    }
    assert((core_index >= 0 && core_index < core_list.size()) && "selected_core_index should be in bounds");
    return core_list[core_index];
}

bool Viewer::erase_core(const size_t index)
{
    assert((index >= 0 && index < core_list.size()) && "index should be in bounds");
    if (core_list.size() == 1)
    {
        // Cannot remove last viewport
        return false;
    }
    core_list[index].shut(); // does nothing
    core_list.erase(core_list.begin() + index);
    if (selected_core_index >= index && selected_core_index > 0)
    {
        selected_core_index--;
    }
    return true;
}

size_t Viewer::core_index(const int id) const
{
    for (size_t i = 0; i < core_list.size(); ++i)
    {
        if (core_list[i].id == id)
        {
            return i;
        }
    }
    return 0;
}

int Viewer::append_core(Eigen::Vector4f viewport, bool append_empty /*= false*/)
{
    core_list.push_back(core()); // copies the previous active core and only changes the viewport
    core_list.back().viewport = viewport;
    core_list.back().id = next_core_id;
    next_core_id <<= 1;
    if (!append_empty)
    {
        for (auto& data: data_list)
        {
            data.set_visible(true, core_list.back().id);
            data.copy_options(core(), core_list.back());
        }
    }
    selected_core_index = core_list.size() - 1;
    return (int) core_list.back().id;
}
