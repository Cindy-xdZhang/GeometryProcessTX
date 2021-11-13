#ifndef VIEWER_WIDGET_H
#define VIEWER_WIDGET_H

#include <memory>
#include <imgui/imgui.h>

#include "viewer/Viewer.h"


struct ImGuiContext;

class ViewerWidget
{

public: // enums

    enum AnchorType {TopLeft, TopRight, BottomLeft, BottomRight};

public:

    explicit ViewerWidget(bool expanded = true);

    virtual void init(Viewer* viewer);
    virtual void shutdown();

    // Callbacks
    virtual bool post_load();
    virtual bool pre_draw(bool first);
    virtual ImVec4 draw(bool first, float scaling, float xSize, float ySize, float xPos, float yPos);
    virtual bool post_draw(bool first);
    virtual bool post_resize(int width, int height);

    // Mouse IO
    virtual bool mouse_down(int button, int modifier);
    virtual bool mouse_up(int button, int modifier);
    virtual bool mouse_move(int mouse_x, int mouse_y);
    virtual bool mouse_scroll(float delta_y);

    // Keyboard IO
    virtual bool key_pressed(unsigned int key, int modifiers);
    virtual bool key_down(int key, int modifiers);
    virtual bool key_up(int key, int modifiers);
    virtual bool key_repeat(int key, int modifiers);

protected:

    // Pointer to the main Viewer class
    Viewer* mViewer;

    // Whether the widget is to start up collapsed
    bool mCollapsed;

    // How to align the widget
    ViewerWidget::AnchorType mAnchor;
};

#endif // VIEWER_WIDGET_H
