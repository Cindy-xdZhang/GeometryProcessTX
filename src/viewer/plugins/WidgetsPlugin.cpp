#include "viewer/plugins/WidgetsPlugin.h"

#include <iostream>
#include <GLFW/glfw3.h>
#include <imgui/backends/imgui_impl_glfw.h>
#include <imgui/backends/imgui_impl_opengl3.h>

#include "viewer/imgui/imgui_fonts_droid_sans.h"

#include <io.h>
#if defined __unix__
#define DEVNULL "/dev/null"
#elif defined _WIN32
#define DEVNULL "nul"
#endif

WidgetsPlugin::WidgetsPlugin() : ViewerPlugin()
{
    hidpi_scaling_ = -1.0f;
    pixel_ratio_ = -1.0f;
}

WidgetsPlugin::~WidgetsPlugin()
{
    for (ViewerWidget* widget: widgets)
    {
        delete widget;
    }
    widgets.clear();
}

void WidgetsPlugin::init(Viewer* viewer)
{
    ViewerPlugin::init(viewer);

    // Setup ImGui binding
    if (viewer)
    {
        IMGUI_CHECKVERSION();
        if (!context_)
        {
            // Single global context by default, but can be overridden by the user
            static ImGuiContext* _global_context = ImGui::CreateContext();
            context_ = _global_context;
        }
        const char* glsl_version = "#version 150";
        ImGui_ImplGlfw_InitForOpenGL(viewer->window, false);
        ImGui_ImplOpenGL3_Init(glsl_version);
        ImGui::GetIO().IniFilename = nullptr;
        ImGui::StyleColorsDark();
        ImGuiStyle& style = ImGui::GetStyle();
        style.FrameRounding = 5.0f;
        reload_font();
    }

    for (ViewerWidget* widget: widgets)
    {
        widget->init(viewer);
    }
}

void WidgetsPlugin::shutdown()
{
    for (ViewerWidget* widget: widgets)
    {
        widget->shutdown();
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    // User is responsible for destroying context if a custom context is given
    // ImGui::DestroyContext(*context_);
}

bool WidgetsPlugin::post_load()
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->post_load())
        {
            //return true;
        }
    }

    return false;
}

bool WidgetsPlugin::pre_draw(bool first)
{
    glfwPollEvents();

    // Check whether window dpi has changed
    float scaling = new_hidpi_scaling();
    if (std::abs(scaling - hidpi_scaling()) > 1e-5)
    {
        reload_font();
        ImGui_ImplOpenGL3_DestroyDeviceObjects();
    }

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    for (ViewerWidget* widget: widgets)
    {
        if (widget->pre_draw(first))
        {
            return true;
        }
    }

    return false;
}

void WidgetsPlugin::draw(bool first)
{
    float x_scale, y_scale;
    glfwGetWindowContentScale(mViewer->window, &x_scale, &y_scale);
    int xsize, ysize;
    glfwGetWindowSize(mViewer->window, &xsize, &ysize);
    float x_size = float(xsize) / x_scale;
    float y_size = float(ysize) / y_scale;
    float x_spacing = 10.0f / x_scale;
    float y_spacing = 10.0f / y_scale;

    float x = x_spacing;
    float y = y_spacing;
    for (auto widget: widgets)
    {
        ImVec4 dim = widget->draw(first, menu_scaling(), x_size, y_size, x, y);
        //y += dim.w / y_scale + y_spacing;
        x += dim.z / x_scale + x_spacing;
    }
}

bool WidgetsPlugin::post_draw(bool first)
{
    draw(first);

    ImGui::Render();
#if defined _WIN32  // Dear ImGui raises error when unloading mesh
    int old = _dup(_fileno(stderr));
    FILE* file = nullptr;
    fopen_s(&file, DEVNULL, "w");
    _dup2(_fileno(file), _fileno(stderr));
    //fprintf(stderr, "ignored!\n");
#endif
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
#if defined _WIN32
    fflush(file);
    _dup2(old, _fileno(stderr));
    fclose(file);
    //fprintf(stderr, "shown\n");
#endif

    for (ViewerWidget* widget: widgets)
    {
        if (widget->post_draw(first))
        {
            return true;
        }
    }

    return false;
}

bool WidgetsPlugin::post_resize(int w, int h)
{
    if (context_)
    {
        ImGui::GetIO().DisplaySize.x = float(w);
        ImGui::GetIO().DisplaySize.y = float(h);
    }

    for (ViewerWidget* widget: widgets)
    {
        if (widget->post_resize(w, h))
        {
            return true;
        }
    }

    return false;
}

bool WidgetsPlugin::mouse_down(int button, int modifier)
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->mouse_down(static_cast<int>(button), modifier))
        {
            return true;
        }
    }

    ImGui_ImplGlfw_MouseButtonCallback(mViewer->window, button, GLFW_PRESS, modifier);
    return ImGui::GetIO().WantCaptureMouse;
}

bool WidgetsPlugin::mouse_up(int button, int modifier)
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->mouse_up(static_cast<int>(button), modifier))
        {
            return false; // should not steal mouse up
        }
    }

    //return ImGui::GetIO().WantCaptureMouse;
    // !! Should not steal mouse up
    return false;
}

bool WidgetsPlugin::mouse_move(int mouse_x, int mouse_y)
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->mouse_move(mouse_x, mouse_y))
        {
            return true;
        }
    }

    return ImGui::GetIO().WantCaptureMouse;
}

bool WidgetsPlugin::mouse_scroll(float delta_y)
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->mouse_scroll(delta_y))
        {
            return true;
        }
    }

    ImGui_ImplGlfw_ScrollCallback(mViewer->window, 0.f, delta_y);
    return ImGui::GetIO().WantCaptureMouse;
}

bool WidgetsPlugin::key_pressed(unsigned int key, int modifiers)
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->key_pressed(key, modifiers))
        {
            return true;
        }
    }

    ImGui_ImplGlfw_CharCallback(nullptr, key);
    return ImGui::GetIO().WantCaptureKeyboard;
}

bool WidgetsPlugin::key_down(int key, int modifiers)
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->key_down(key, modifiers))
        {
            return true;
        }
    }

    ImGui_ImplGlfw_KeyCallback(mViewer->window, key, 0, GLFW_PRESS, modifiers);
    return ImGui::GetIO().WantCaptureKeyboard;
}

bool WidgetsPlugin::key_up(int key, int modifiers)
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->key_up(key, modifiers))
        {
            return true;
        }
    }

    ImGui_ImplGlfw_KeyCallback(mViewer->window, key, 0, GLFW_RELEASE, modifiers);
    return ImGui::GetIO().WantCaptureKeyboard;
}

bool WidgetsPlugin::key_repeat(int key, int modifiers)
{
    for (ViewerWidget* widget: widgets)
    {
        if (widget->key_repeat(key, modifiers))
        {
            return true;
        }
    }

    ImGui_ImplGlfw_KeyCallback(mViewer->window, key, 0, GLFW_REPEAT, modifiers);
    return ImGui::GetIO().WantCaptureKeyboard;
}

void WidgetsPlugin::add(ViewerWidget* widget)
{
    widgets.push_back(widget);
}

void WidgetsPlugin::reload_font()
{
    reload_font(13);
}

void WidgetsPlugin::reload_font(int font_size)
{
    hidpi_scaling_ = new_hidpi_scaling();
    pixel_ratio_ = new_pixel_ratio();
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();
    io.Fonts->AddFontFromMemoryCompressedTTF(droid_sans_compressed_data,
                                             droid_sans_compressed_size, (float) font_size * hidpi_scaling_);
    io.FontGlobalScale = 1.0f / pixel_ratio_;
}

float WidgetsPlugin::menu_scaling() const
{
    return hidpi_scaling_ / pixel_ratio_;
}

float WidgetsPlugin::pixel_ratio() const
{
    return pixel_ratio_;
}

float WidgetsPlugin::hidpi_scaling() const
{
    return hidpi_scaling_;
}

float WidgetsPlugin::new_pixel_ratio()
{
    // Computes pixel ratio for hidpi devices
    int buf_size[2];
    int win_size[2];
    GLFWwindow* window = glfwGetCurrentContext();
    glfwGetFramebufferSize(window, &buf_size[0], &buf_size[1]);
    glfwGetWindowSize(window, &win_size[0], &win_size[1]);
    return (float) buf_size[0] / (float) win_size[0];
}

float WidgetsPlugin::new_hidpi_scaling()
{
    // Computes scaling factor for hidpi devices
    float xscale, yscale;
    GLFWwindow* window = glfwGetCurrentContext();
    glfwGetWindowContentScale(window, &xscale, &yscale);
    return 0.5f * (xscale + yscale);
}
