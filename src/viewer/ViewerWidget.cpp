#include "viewer/ViewerWidget.h"


ViewerWidget::ViewerWidget(bool expanded)
{
    mViewer = nullptr;
    mCollapsed = !expanded;
    mAnchor = ViewerWidget::AnchorType::TopLeft;
}

void ViewerWidget::init(Viewer* viewer)
{
    mViewer = viewer;
}

void ViewerWidget::shutdown()
{

}

bool ViewerWidget::post_load()
{
    return false;
}

bool ViewerWidget::pre_draw(bool first)
{
    return false;
}

ImVec4 ViewerWidget::draw(bool first, float scaling, float xSize, float ySize, float xPos, float yPos)
{
    return {};
}

bool ViewerWidget::post_draw(bool first)
{
    return false;
}

bool ViewerWidget::post_resize(int width, int height)
{
    return false;
}

bool ViewerWidget::mouse_down(int button, int modifier)
{
    return false;
}

bool ViewerWidget::mouse_up(int button, int modifier)
{
    return false;
}

bool ViewerWidget::mouse_move(int mouse_x, int mouse_y)
{
    return false;
}

bool ViewerWidget::mouse_scroll(float delta_y)
{
    return false;
}

bool ViewerWidget::key_pressed(unsigned int key, int modifiers)
{
    return false;
}

bool ViewerWidget::key_down(int key, int modifiers)
{
    return false;
}

bool ViewerWidget::key_up(int key, int modifiers)
{
    return false;
}

bool ViewerWidget::key_repeat(int key, int modifiers)
{
    return false;
}
