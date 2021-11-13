#ifndef UTILS_SYSTEM_H
#define UTILS_SYSTEM_H

#include <string>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#include <Commdlg.h>
#undef max
#undef min
#define getcwd _getcwd
#define PATH_MAX 4096
#else
#include <unistd.h>
#endif


namespace utils
{
    template<class T>
    T base_name(T const& path, T const& delimiters = "/\\")
    {
        return path.substr(path.find_last_of(delimiters) + 1);
    }

    template<class T>
    T remove_extension(T const& filename)
    {
        typename T::size_type const p(filename.find_last_of('.'));
        return p > 0 && p != T::npos ? filename.substr(0, p) : filename;
    }

    inline bool folder_exists(const std::string& name)
    {
        struct stat st{};
        stat(name.c_str(), &st);
        return st.st_mode & S_IFDIR;
    }

    inline int _mkdir(const char* path)
    {
#ifdef _WIN32
        return ::_mkdir(path);
#else
#if _POSIX_C_SOURCE
        return ::mkdir(path, ACCESSPERMS);
#else
        return ::mkdir(path, 0755); // not sure if this works on mac
#endif
#endif
    }

    inline int mkdir(const char* path)
    {
        std::string current_level;
        std::string level;
        std::stringstream ss(path);

        // split path using slash as a separator
        while (std::getline(ss, level, '/'))
        {
            current_level += level; // append folder to the current level

            // create current level
            if (!level.empty())
            {
                if (!folder_exists(current_level) && _mkdir(current_level.c_str()) != 0)
                {
                    return -1;
                }
            }

            current_level += "/"; // don't forget to append a slash
        }

        return 0;
    }

    inline int mkdir(const std::string& path)
    {
        return mkdir(path.c_str());
    }

    inline std::string getWorkingDirectory()
    {
        char cwd[PATH_MAX];
        if (getcwd(cwd, sizeof(cwd)) != nullptr)
        {
            return cwd;
        }
        else
        {
            return "";
        }
    }

    inline void printProgress(float progress, int barWidth = 70)
    {
        std::cout << "[";
        int pos = (int) ((float) barWidth * progress);
        for (int i = 0; i < barWidth; ++i)
        {
            if (i < pos)
            {
                std::cout << "=";
            }
            else if (i == pos)
            {
                std::cout << ">";
            }
            else
            {
                std::cout << " ";
            }
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Dialogs
    ////////////////////////////////////////////////////////////////////////////

    // Returns a string with a path to an existing file
    // The string is returned empty if no file is selected
    // (on Linux machines, it assumes that Zenity is installed)
    //
    // Usage:
    //   std::string str = get_open_file_path();
    inline std::string file_dialog_open(const std::string& path = ".")
    {
        const int FILE_DIALOG_MAX_BUFFER = 1024;
        char buffer[FILE_DIALOG_MAX_BUFFER];
        buffer[0] = '\0';
        buffer[FILE_DIALOG_MAX_BUFFER - 1] = 'x'; // Initialize last character with a char != '\0'

#ifdef __APPLE__
        // For apple use applescript hack
    FILE * output = popen(
      "osascript -e \""
      "   tell application \\\"System Events\\\"\n"
      "           activate\n"
      "           set existing_file to choose file\n"
      "   end tell\n"
      "   set existing_file_path to (POSIX path of (existing_file))\n"
      "\" 2>/dev/null | tr -d '\n' ","r");
    if (output)
    {
      auto ret = fgets(buffer, FILE_DIALOG_MAX_BUFFER, output);
      if (ret == NULL || ferror(output))
      {
        // I/O error
        buffer[0] = '\0';
      }
      if (buffer[FILE_DIALOG_MAX_BUFFER - 1] == '\0')
      {
        // File name too long, buffer has been filled, so we return empty string instead
        buffer[0] = '\0';
      }
    }
#elif defined _WIN32
        // Use native windows file dialog box
    // (code contributed by Tino Weinkauf)

    OPENFILENAME ofn;       // common dialog box structure
    char szFile[260];       // buffer for file name

    // Initialize OPENFILENAME
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;
    ofn.lpstrFile = szFile;
    // Set lpstrFile[0] to '\0' so that GetOpenFileName does not
    // use the contents of szFile to initialize itself.
    ofn.lpstrFile[0] = '\0';
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "*.*\0";//off\0*.off\0obj\0*.obj\0mp\0*.mp\0";
    ofn.nFilterIndex = 1;
    ofn.lpstrFileTitle = NULL;
    ofn.nMaxFileTitle = 0;
    ofn.lpstrInitialDir = NULL;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    // Display the Open dialog box.
    int pos = 0;
    if (GetOpenFileName(&ofn)==TRUE)
    {
      while(ofn.lpstrFile[pos] != '\0')
      {
        buffer[pos] = (char)ofn.lpstrFile[pos];
        pos++;
      }
    }
    buffer[pos] = 0;
#else
        // For linux use zenity
        FILE* output = popen(("/usr/bin/zenity --file-selection --filename=" + path).c_str(), "r");
        if (output)
        {
            auto ret = fgets(buffer, FILE_DIALOG_MAX_BUFFER, output);
            if (ret == nullptr || ferror(output))
            {
                // I/O error
                buffer[0] = '\0';
            }
            if (buffer[FILE_DIALOG_MAX_BUFFER - 1] == '\0')
            {
                // File name too long, buffer has been filled, so we return empty string instead
                buffer[0] = '\0';
            }
        }

        // Replace last '\n' by '\0'
        if (strlen(buffer) > 0)
        {
            buffer[strlen(buffer) - 1] = '\0';
        }

#endif
        return std::string(buffer);
    }

    // Returns a string with a path to a new/existing file
    // The string is returned empty if no file is selected
    // (on Linux machines, it assumes that Zenity is installed)
    //
    // Usage:
    //   char buffer[FILE_DIALOG_MAX_BUFFER];
    //   get_save_file_path(buffer);
    inline std::string file_dialog_save(const std::string& path)
    {
        const int FILE_DIALOG_MAX_BUFFER = 1024;
        char buffer[FILE_DIALOG_MAX_BUFFER];
        buffer[0] = '\0';
        buffer[FILE_DIALOG_MAX_BUFFER - 1] = 'x'; // Initialize last character with a char != '\0'

#ifdef __APPLE__
        // For apple use applescript hack
    // There is currently a bug in Applescript that strips extensions off
    // of chosen existing files in the "choose file name" dialog
    // I'm assuming that will be fixed soon
    FILE * output = popen(
      "osascript -e \""
      "   tell application \\\"System Events\\\"\n"
      "           activate\n"
      "           set existing_file to choose file name\n"
      "   end tell\n"
      "   set existing_file_path to (POSIX path of (existing_file))\n"
      "\" 2>/dev/null | tr -d '\n' ","r");
    if (output)
    {
      auto ret = fgets(buffer, FILE_DIALOG_MAX_BUFFER, output);
      if (ret == NULL || ferror(output))
      {
        // I/O error
        buffer[0] = '\0';
      }
      if (buffer[FILE_DIALOG_MAX_BUFFER - 1] == '\0')
      {
        // File name too long, buffer has been filled, so we return empty string instead
        buffer[0] = '\0';
      }
    }
#elif defined _WIN32
        // Use native windows file dialog box
    // (code contributed by Tino Weinkauf)

    OPENFILENAME ofn;       // common dialog box structure
    char szFile[260];       // buffer for file name

    // Initialize OPENFILENAME
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;//hwnd;
    ofn.lpstrFile = szFile;
    // Set lpstrFile[0] to '\0' so that GetOpenFileName does not
    // use the contents of szFile to initialize itself.
    ofn.lpstrFile[0] = '\0';
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "";
    ofn.nFilterIndex = 1;
    ofn.lpstrFileTitle = NULL;
    ofn.nMaxFileTitle = 0;
    ofn.lpstrInitialDir = NULL;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    // Display the Open dialog box.
    int pos = 0;
    if (GetSaveFileName(&ofn)==TRUE)
    {
      while(ofn.lpstrFile[pos] != '\0')
      {
        buffer[pos] = (char)ofn.lpstrFile[pos];
        pos++;
      }
      buffer[pos] = 0;
    }
#else
        // For every other machine type use zenity
        FILE* output = popen(("/usr/bin/zenity --file-selection --save --filename=" + path).c_str(), "r");
        if (output)
        {
            auto ret = fgets(buffer, FILE_DIALOG_MAX_BUFFER, output);
            if (ret == nullptr || ferror(output))
            {
                // I/O error
                buffer[0] = '\0';
            }
            if (buffer[FILE_DIALOG_MAX_BUFFER - 1] == '\0')
            {
                // File name too long, buffer has been filled, so we return empty string instead
                buffer[0] = '\0';
            }
        }

        // Replace last '\n' by '\0'
        if (strlen(buffer) > 0)
        {
            buffer[strlen(buffer) - 1] = '\0';
        }

#endif
        return std::string(buffer);
    }
}

#endif // UTILS_SYSTEM_H
