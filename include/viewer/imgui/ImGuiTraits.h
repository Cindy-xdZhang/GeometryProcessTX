#ifndef IMGUITRAITS_H
#define IMGUITRAITS_H

#include <imgui/imgui.h>


namespace ImGui
{
    // Infer ImGuiDataType enum based on actual type
    template<typename T>
    class ImGuiDataTypeTraits
    {
        static const ImGuiDataType value; // link error
        static const char* format;
    };

    template<>
    class ImGuiDataTypeTraits<int>
    {
        static constexpr ImGuiDataType value = ImGuiDataType_S32;
        static constexpr const char* format = "%d";
    };

    template<>
    class ImGuiDataTypeTraits<unsigned int>
    {
        static constexpr ImGuiDataType value = ImGuiDataType_U32;
        static constexpr const char* format = "%u";
    };

    template<>
    class ImGuiDataTypeTraits<long long>
    {
        static constexpr ImGuiDataType value = ImGuiDataType_S64;
        static constexpr const char* format = "%lld";
    };

    template<>
    class ImGuiDataTypeTraits<unsigned long long>
    {
        static constexpr ImGuiDataType value = ImGuiDataType_U64;
        static constexpr const char* format = "%llu";
    };

    template<>
    class ImGuiDataTypeTraits<float>
    {
        static constexpr ImGuiDataType value = ImGuiDataType_Float;
        static constexpr const char* format = "%.3f";
    };

    template<>
    class ImGuiDataTypeTraits<double>
    {
        static constexpr ImGuiDataType value = ImGuiDataType_Double;
        static constexpr const char* format = "%.6f";
    };

}

#endif // IMGUITRAITS_H
