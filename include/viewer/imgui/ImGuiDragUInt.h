#ifndef IMGUI_DRAG_UINT_H
#define IMGUI_DRAG_UINT_H

#include <imgui/imgui.h>


namespace ImGui
{
    inline bool
    DragUInt(const char* label, unsigned int* v, float v_speed = 1.0f,
             unsigned int v_min = 0, unsigned int v_max = 0,
             const char* format = "%d", ImGuiSliderFlags flags = 0)
    {
        return DragScalar(label, ImGuiDataType_S32, v, v_speed, &v_min, &v_max, format, flags);
    }
}

#endif // IMGUI_DRAG_UINT_H
