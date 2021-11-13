#ifndef IMGUIHELPERS_H
#define IMGUIHELPERS_H

#include "ImGuiTraits.h"
#include <imgui/imgui.h>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <cstddef>


// Code snippets taken from there:
// https://eliasdaler.github.io/using-imgui-with-sfml-pt2/

namespace ImGui
{
    static auto vector_getter = [](void* vec, int idx, const char** out_text)
    {
        auto& vector = *static_cast<std::vector<std::string>*>(vec);
        if (idx < 0 || idx >= static_cast<int>(vector.size()))
        {
            return false;
        }
        *out_text = vector.at(idx).c_str();
        return true;
    };

    inline bool Combo(const char* label, int* idx, std::vector<std::string>& values)
    {
        if (values.empty())
        {
            return false;
        }
        return Combo(label, idx, vector_getter,
                     static_cast<void*>(&values), static_cast<int>(values.size()));
    }

    inline bool Combo(const char* label, int* idx, std::function<const char*(int)> getter, int items_count)
    {
        auto func = [](void* data, int i, const char** out_text)
        {
            auto& getter = *reinterpret_cast<std::function<const char*(int)>*>(data);
            const char* s = getter(i);
            if (s)
            {
                *out_text = s;
                return true;
            }
            else
            {
                return false;
            }
        };
        return Combo(label, idx, func, reinterpret_cast<void*>(&getter), items_count);
    }

    inline bool ListBox(const char* label, int* idx, std::vector<std::string>& values)
    {
        if (values.empty())
        {
            return false;
        }
        return ListBox(label, idx, vector_getter,
                       static_cast<void*>(&values), static_cast<int>(values.size()));
    }

    inline bool
    InputText(const char* label, std::string& str, ImGuiInputTextFlags flags = 0, ImGuiTextEditCallback callback = NULL,
              void* user_data = NULL)
    {
        char buf[1024];
        std::fill_n(buf, 1024, 0);
        std::copy_n(str.begin(), std::min(1024, (int) str.size()), buf);
        if (ImGui::InputText(label, buf, 1024, flags, callback, user_data))
        {
            str = std::string(buf);
            return true;
        }
        return false;
    }

    template<typename T>
    inline bool SliderScalar(const char* label, T* value, T min = 0, T max = 0, const char* format = "")
    {
        const char* fmt = format;
        if (format == nullptr)
        {
            fmt = ImGuiDataTypeTraits<T>::format;
        }
        return SliderScalar(label, ImGuiDataTypeTraits<T>::value, value, &min, &max, fmt);
    }

    template<typename Getter, typename Setter>
    inline bool Checkbox(const char* label, Getter get, Setter set)
    {
        bool value = get();
        bool ret = ImGui::Checkbox(label, &value);
        set(value);
        return ret;
    }

}

#endif // IMGUIHELPERS_H
