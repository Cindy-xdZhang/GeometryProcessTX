#ifndef MAIN_WIDGET_H
#define MAIN_WIDGET_H

#include "viewer/ViewerWidget.h"
#include "viewer/plugins/LoaderPlugin.h"


class MainWidget : public ViewerWidget
{

public:

    explicit MainWidget(bool expanded, LoaderPlugin* loader = nullptr);

    ImVec4 draw(bool first, float scaling, float xSize, float ySize, float xPos, float yPos) override;

private:

    LoaderPlugin* mLoader;
};

#endif // MAIN_WIDGET_H
