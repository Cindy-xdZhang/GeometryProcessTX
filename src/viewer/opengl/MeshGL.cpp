#include "viewer/opengl/MeshGL.h"

#include <iostream>
#include "utils/VeraSansMono.h"


MeshGL::MeshGL()
{
    tex_filter = GL_LINEAR;
    tex_wrap = GL_REPEAT;
}

void MeshGL::init_buffers()
{
    // Mesh: Vertex Array Object & Buffer objects
    glGenVertexArrays(1, &vao_mesh);
    glBindVertexArray(vao_mesh);
    glGenBuffers(1, &vbo_V);
    glGenBuffers(1, &vbo_V_normals);
    glGenBuffers(1, &vbo_V_ambient);
    glGenBuffers(1, &vbo_V_diffuse);
    glGenBuffers(1, &vbo_V_specular);
    glGenBuffers(1, &vbo_V_uv);
    glGenBuffers(1, &vbo_F);
    glGenTextures(1, &vbo_tex);
    glGenTextures(1, &font_atlas_id);

    // Line overlay
    glGenVertexArrays(1, &vao_overlay_lines);
    glBindVertexArray(vao_overlay_lines);
    glGenBuffers(1, &vbo_lines_F);
    glGenBuffers(1, &vbo_lines_V);
    glGenBuffers(1, &vbo_lines_V_colors);

    // Point overlay
    glGenVertexArrays(1, &vao_overlay_points);
    glBindVertexArray(vao_overlay_points);
    glGenBuffers(1, &vbo_points_F);
    glGenBuffers(1, &vbo_points_V);
    glGenBuffers(1, &vbo_points_V_colors);

    // Text Labels
    vertex_labels.init_buffers();
    face_labels.init_buffers();
    custom_labels.init_buffers();

    dirty = MeshGL::DIRTY_ALL;
}

GLint MeshGL::bind_vertex_attrib_array(
        const GLuint& program_shader,
        const std::string& name,
        GLuint bufferID,
        const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& M,
        bool refresh)
{
    GLint id = glGetAttribLocation(program_shader, name.c_str());
    if (id < 0)
    {
        return id;
    }
    if (M.size() == 0)
    {
        glDisableVertexAttribArray(id);
        return id;
    }
    glBindBuffer(GL_ARRAY_BUFFER, bufferID);
    if (refresh)
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * M.size(), M.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(id, M.cols(), GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(id);
    return id;
}

void MeshGL::free_buffers()
{
    if (is_initialized)
    {
        glDeleteVertexArrays(1, &vao_mesh);
        glDeleteVertexArrays(1, &vao_overlay_lines);
        glDeleteVertexArrays(1, &vao_overlay_points);

        glDeleteBuffers(1, &vbo_V);
        glDeleteBuffers(1, &vbo_V_normals);
        glDeleteBuffers(1, &vbo_V_ambient);
        glDeleteBuffers(1, &vbo_V_diffuse);
        glDeleteBuffers(1, &vbo_V_specular);
        glDeleteBuffers(1, &vbo_V_uv);
        glDeleteBuffers(1, &vbo_F);
        glDeleteBuffers(1, &vbo_lines_F);
        glDeleteBuffers(1, &vbo_lines_V);
        glDeleteBuffers(1, &vbo_lines_V_colors);
        glDeleteBuffers(1, &vbo_points_F);
        glDeleteBuffers(1, &vbo_points_V);
        glDeleteBuffers(1, &vbo_points_V_colors);

        // Text Labels
        vertex_labels.free_buffers();
        face_labels.free_buffers();
        custom_labels.free_buffers();

        glDeleteTextures(1, &vbo_tex);
        glDeleteTextures(1, &font_atlas_id);
    }
}

void MeshGL::TextGL::init_buffers()
{
    glGenVertexArrays(1, &vao_labels);
    glBindVertexArray(vao_labels);
    glGenBuffers(1, &vbo_labels_pos);
    glGenBuffers(1, &vbo_labels_characters);
    glGenBuffers(1, &vbo_labels_offset);
    glGenBuffers(1, &vbo_labels_indices);
}

void MeshGL::TextGL::free_buffers()
{
    glDeleteBuffers(1, &vbo_labels_pos);
    glDeleteBuffers(1, &vbo_labels_characters);
    glDeleteBuffers(1, &vbo_labels_offset);
    glDeleteBuffers(1, &vbo_labels_indices);
}

void MeshGL::bind_mesh()
{
    glBindVertexArray(vao_mesh);
    glUseProgram(shader_mesh);
    bind_vertex_attrib_array(shader_mesh, "position", vbo_V, V_vbo, dirty & MeshGL::DIRTY_POSITION);
    bind_vertex_attrib_array(shader_mesh, "normal", vbo_V_normals, V_normals_vbo, dirty & MeshGL::DIRTY_NORMAL);
    bind_vertex_attrib_array(shader_mesh, "Ka", vbo_V_ambient, V_ambient_vbo, dirty & MeshGL::DIRTY_AMBIENT);
    bind_vertex_attrib_array(shader_mesh, "Kd", vbo_V_diffuse, V_diffuse_vbo, dirty & MeshGL::DIRTY_DIFFUSE);
    bind_vertex_attrib_array(shader_mesh, "Ks", vbo_V_specular, V_specular_vbo, dirty & MeshGL::DIRTY_SPECULAR);
    bind_vertex_attrib_array(shader_mesh, "texcoord", vbo_V_uv, V_uv_vbo, dirty & MeshGL::DIRTY_UV);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_F);
    if (dirty & MeshGL::DIRTY_FACE)
    {
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * F_vbo.size(), F_vbo.data(), GL_DYNAMIC_DRAW);
    }

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, vbo_tex);
    if (dirty & MeshGL::DIRTY_TEXTURE)
    {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, tex_wrap);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, tex_wrap);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, tex_filter);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, tex_filter);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex_u, tex_v, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex.data());
    }
    glUniform1i(glGetUniformLocation(shader_mesh, "tex"), 0);
    dirty &= ~MeshGL::DIRTY_MESH;
}

void MeshGL::bind_overlay_lines()
{
    bool is_dirty = dirty & MeshGL::DIRTY_OVERLAY_LINES;

    glBindVertexArray(vao_overlay_lines);
    glUseProgram(shader_overlay_lines);
    bind_vertex_attrib_array(shader_overlay_lines, "position", vbo_lines_V, lines_V_vbo, is_dirty);
    bind_vertex_attrib_array(shader_overlay_lines, "color", vbo_lines_V_colors, lines_V_colors_vbo, is_dirty);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_lines_F);
    if (is_dirty)
    {
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * lines_F_vbo.size(), lines_F_vbo.data(),
                     GL_DYNAMIC_DRAW);
    }

    dirty &= ~MeshGL::DIRTY_OVERLAY_LINES;
}

void MeshGL::bind_overlay_points()
{
    bool is_dirty = dirty & MeshGL::DIRTY_OVERLAY_POINTS;

    glBindVertexArray(vao_overlay_points);
    glUseProgram(shader_overlay_points);
    bind_vertex_attrib_array(shader_overlay_points, "position", vbo_points_V, points_V_vbo, is_dirty);
    bind_vertex_attrib_array(shader_overlay_points, "color", vbo_points_V_colors, points_V_colors_vbo, is_dirty);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_points_F);
    if (is_dirty)
    {
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * points_F_vbo.size(), points_F_vbo.data(),
                     GL_DYNAMIC_DRAW);
    }

    dirty &= ~MeshGL::DIRTY_OVERLAY_POINTS;
}

void MeshGL::init_text_rendering()
{
    // Decompress the png of the font atlas
    unsigned char font_atlas[256 * 256];
    VeraSansMono::decompress_atlas(font_atlas);

    // Bind atlas
    glBindTexture(GL_TEXTURE_2D, font_atlas_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 256, 256, 0, GL_RED, GL_UNSIGNED_BYTE, font_atlas);

    // TextGL initialization
    vertex_labels.dirty_flag = MeshGL::DIRTY_VERTEX_LABELS;
    face_labels.dirty_flag = MeshGL::DIRTY_FACE_LABELS;
    custom_labels.dirty_flag = MeshGL::DIRTY_CUSTOM_LABELS;
}

void MeshGL::bind_labels(const TextGL& labels)
{
    bool is_dirty = dirty & labels.dirty_flag;
    glBindTexture(GL_TEXTURE_2D, font_atlas_id);
    glBindVertexArray(labels.vao_labels);
    glUseProgram(shader_text);
    bind_vertex_attrib_array(shader_text, "position", labels.vbo_labels_pos, labels.label_pos_vbo, is_dirty);
    bind_vertex_attrib_array(shader_text, "character", labels.vbo_labels_characters, labels.label_char_vbo, is_dirty);
    bind_vertex_attrib_array(shader_text, "offset", labels.vbo_labels_offset, labels.label_offset_vbo, is_dirty);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, labels.vbo_labels_indices);
    if (is_dirty)
    {
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * labels.label_indices_vbo.size(),
                     labels.label_indices_vbo.data(), GL_DYNAMIC_DRAW);
    }
    dirty &= ~labels.dirty_flag;
}

void MeshGL::draw_mesh(bool solid)
{
    glPolygonMode(GL_FRONT_AND_BACK, solid ? GL_FILL : GL_LINE);

    /* Avoid Z-buffer fighting between filled triangles & wireframe lines */
    if (solid)
    {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
    }
    glDrawElements(GL_TRIANGLES, 3 * F_vbo.rows(), GL_UNSIGNED_INT, nullptr);

    glDisable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshGL::draw_overlay_lines()
{
    glDrawElements(GL_LINES, lines_F_vbo.rows(), GL_UNSIGNED_INT, nullptr);
}

void MeshGL::draw_overlay_points()
{
    glDrawElements(GL_POINTS, points_F_vbo.rows(), GL_UNSIGNED_INT, nullptr);
}

void MeshGL::draw_labels(const TextGL& labels)
{
    glDrawElements(GL_POINTS, labels.label_indices_vbo.rows(), GL_UNSIGNED_INT, nullptr);
}

void MeshGL::init()
{
    if (is_initialized)
    {
        return;
    }
    is_initialized = true;
    std::string mesh_vertex_shader_string =
            R"(#version 150
  uniform mat4 view;
  uniform mat4 proj;
  uniform mat4 normal_matrix;
  in vec3 position;
  in vec3 normal;
  out vec3 position_eye;
  out vec3 normal_eye;
  in vec4 Ka;
  in vec4 Kd;
  in vec4 Ks;
  in vec2 texcoord;
  out vec2 texcoordi;
  out vec4 Kai;
  out vec4 Kdi;
  out vec4 Ksi;

  void main()
  {
    position_eye = vec3 (view * vec4 (position, 1.0));
    normal_eye = vec3 (normal_matrix * vec4 (normal, 0.0));
    normal_eye = normalize(normal_eye);
    gl_Position = proj * vec4 (position_eye, 1.0); //proj * view * vec4(position, 1.0);"
    Kai = Ka;
    Kdi = Kd;
    Ksi = Ks;
    texcoordi = texcoord;
  }
)";

    std::string mesh_fragment_shader_string =
            R"(#version 150
  uniform mat4 view;
  uniform mat4 proj;
  uniform vec4 fixed_color;
  in vec3 position_eye;
  in vec3 normal_eye;
  uniform vec3 light_position_eye;
  vec3 Ls = vec3 (1, 1, 1);
  vec3 Ld = vec3 (1, 1, 1);
  vec3 La = vec3 (1, 1, 1);
  in vec4 Ksi;
  in vec4 Kdi;
  in vec4 Kai;
  in vec2 texcoordi;
  uniform sampler2D tex;
  uniform float specular_exponent;
  uniform float lighting_factor;
  uniform float texture_factor;
  uniform float matcap_factor;
  uniform float double_sided;
  out vec4 outColor;
  void main()
  {
    if(matcap_factor == 1.0f)
    {
      vec2 uv = normalize(normal_eye).xy * 0.5 + 0.5;
      outColor = texture(tex, uv);
    }else
    {
      vec3 Ia = La * vec3(Kai);    // ambient intensity

      vec3 vector_to_light_eye = light_position_eye - position_eye;
      vec3 direction_to_light_eye = normalize (vector_to_light_eye);
      float dot_prod = dot (direction_to_light_eye, normalize(normal_eye));
      float clamped_dot_prod = abs(max (dot_prod, -double_sided));
      vec3 Id = Ld * vec3(Kdi) * clamped_dot_prod;    // Diffuse intensity

      vec3 reflection_eye = reflect (-direction_to_light_eye, normalize(normal_eye));
      vec3 surface_to_viewer_eye = normalize (-position_eye);
      float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);
      dot_prod_specular = float(abs(dot_prod)==dot_prod) * abs(max (dot_prod_specular, -double_sided));
      float specular_factor = pow (dot_prod_specular, specular_exponent);
      vec3 Is = Ls * vec3(Ksi) * specular_factor;    // specular intensity
      vec4 color = vec4(lighting_factor * (Is + Id) + Ia + (1.0-lighting_factor) * vec3(Kdi),(Kai.a+Ksi.a+Kdi.a)/3);
      outColor = mix(vec4(1,1,1,1), texture(tex, texcoordi), texture_factor) * color;
      if (fixed_color != vec4(0.0)) outColor = fixed_color;
    }
  }
)";

    std::string overlay_vertex_shader_string =
            R"(#version 150
  uniform mat4 view;
  uniform mat4 proj;
  in vec3 position;
  in vec3 color;
  out vec3 color_frag;

  void main()
  {
    gl_Position = proj * view * vec4 (position, 1.0);
    color_frag = color;
  }
)";

    std::string overlay_fragment_shader_string =
            R"(#version 150
  in vec3 color_frag;
  out vec4 outColor;
  void main()
  {
    outColor = vec4(color_frag, 1.0);
  }
)";

    std::string overlay_point_fragment_shader_string =
            R"(#version 150
  in vec3 color_frag;
  out vec4 outColor;
  void main()
  {
    if (length(gl_PointCoord - vec2(0.5)) > 0.5)
      discard;
    outColor = vec4(color_frag, 1.0);
  }
)";

    std::string text_vert_shader =
            R"(#version 330
    in vec3 position;
    in float character;
    in float offset;
    uniform mat4 view;
    uniform mat4 proj;
    out int vPosition;
    out int vCharacter;
    out float vOffset;
    void main()
    {
      vCharacter = int(character);
      vOffset = offset;
      vPosition = gl_VertexID;
      gl_Position = proj * view * vec4(position, 1.0);
    }
)";

    std::string text_geom_shader =
            R"(#version 150 core
    layout(points) in;
    layout(triangle_strip, max_vertices = 4) out;
    out vec2 gTexCoord;
    uniform mat4 view;
    uniform mat4 proj;
    uniform vec2 CellSize;
    uniform vec2 CellOffset;
    uniform vec2 RenderSize;
    uniform vec2 RenderOrigin;
    uniform float TextShiftFactor;
    in int vPosition[1];
    in int vCharacter[1];
    in float vOffset[1];
    void main()
    {
      // Code taken from https://prideout.net/strings-inside-vertex-buffers
      // Determine the final quad's position and size:
      vec4 P = gl_in[0].gl_Position + vec4( vOffset[0]*TextShiftFactor, 0.0, 0.0, 0.0 ); // 0.04
      vec4 U = vec4(1, 0, 0, 0) * RenderSize.x; // 1.0
      vec4 V = vec4(0, 1, 0, 0) * RenderSize.y; // 1.0

      // Determine the texture coordinates:
      int letter = vCharacter[0]; // used to be the character
      letter = clamp(letter - 32, 0, 96);
      int row = letter / 16 + 1;
      int col = letter % 16;
      float S0 = CellOffset.x + CellSize.x * col;
      float T0 = CellOffset.y + 1 - CellSize.y * row;
      float S1 = S0 + CellSize.x - CellOffset.x;
      float T1 = T0 + CellSize.y;

      // Output the quad's vertices:
      gTexCoord = vec2(S0, T1); gl_Position = P - U - V; EmitVertex();
      gTexCoord = vec2(S1, T1); gl_Position = P + U - V; EmitVertex();
      gTexCoord = vec2(S0, T0); gl_Position = P - U + V; EmitVertex();
      gTexCoord = vec2(S1, T0); gl_Position = P + U + V; EmitVertex();
      EndPrimitive();
    }
)";

    std::string text_frag_shader =
            R"(#version 330
    out vec4 outColor;
    in vec2 gTexCoord;
    uniform sampler2D font_atlas;
    uniform vec3 TextColor;
    void main()
    {
      float A = texture(font_atlas, gTexCoord).r;
      outColor = vec4(TextColor, A);
    }
)";

    init_buffers();
    init_text_rendering();
    create_shader_program(
            mesh_vertex_shader_string,
            mesh_fragment_shader_string,
            {},
            shader_mesh);
    create_shader_program(
            overlay_vertex_shader_string,
            overlay_fragment_shader_string,
            {},
            shader_overlay_lines);
    create_shader_program(
            overlay_vertex_shader_string,
            overlay_point_fragment_shader_string,
            {},
            shader_overlay_points);
    create_shader_program(
            text_geom_shader,
            text_vert_shader,
            text_frag_shader,
            {},
            shader_text);
}

void MeshGL::free()
{
    const auto free = [](GLuint& id)
    {
        if (id)
        {
            destroy_shader_program(id);
            id = 0;
        }
    };

    if (is_initialized)
    {
        free(shader_mesh);
        free(shader_overlay_lines);
        free(shader_overlay_points);
        free(shader_text);
        free_buffers();
    }
}

bool MeshGL::create_shader_program(
        const std::string& geom_source,
        const std::string& vert_source,
        const std::string& frag_source,
        const std::map<std::string, GLuint>& attrib,
        GLuint& id)
{
    using namespace std;
    if (vert_source.empty() && frag_source.empty())
    {
        cerr <<
             "create_shader_program() could not create shader program,"
             " both .vert and .frag source given were empty" << endl;
        return false;
    }

    // create program
    id = glCreateProgram();
    if (id == 0)
    {
        cerr << "create_shader_program() could not create shader program." << endl;
        return false;
    }
    GLuint g = 0, f = 0, v = 0;

    if (!geom_source.empty())
    {
        // load vertex shader
        g = load_shader(geom_source, GL_GEOMETRY_SHADER);
        if (g == 0)
        {
            cerr << "geometry shader failed to compile." << endl;
            return false;
        }
        glAttachShader(id, g);
    }

    if (!vert_source.empty())
    {
        // load vertex shader
        v = load_shader(vert_source, GL_VERTEX_SHADER);
        if (v == 0)
        {
            cerr << "vertex shader failed to compile." << endl;
            return false;
        }
        glAttachShader(id, v);
    }

    if (!frag_source.empty())
    {
        // load fragment shader
        f = load_shader(frag_source, GL_FRAGMENT_SHADER);
        if (f == 0)
        {
            cerr << "fragment shader failed to compile." << endl;
            return false;
        }
        glAttachShader(id, f);
    }

    // loop over attributes
    for (const auto& ait : attrib)
    {
        glBindAttribLocation(
                id,
                ait.second,
                ait.first.c_str());
    }

    // Link program
    glLinkProgram(id);
    const auto& detach = [&id](const GLuint shader)
    {
        if (shader)
        {
            glDetachShader(id, shader);
            glDeleteShader(shader);
        }
    };
    detach(g);
    detach(f);
    detach(v);

    // print log if any
    print_program_info_log(id);

    return true;
}

bool MeshGL::create_shader_program(
        const std::string& vert_source,
        const std::string& frag_source,
        const std::map<std::string, GLuint>& attrib,
        GLuint& id)
{
    return create_shader_program("", vert_source, frag_source, attrib, id);
}

bool MeshGL::destroy_shader_program(const GLuint& id)
{
    // Don't try to destroy id == 0 (no shader program)
    if (id == 0)
    {
        fprintf(stderr, "Error: destroy_shader_program() id = %d"
                        " but must should be positive\n", id);
        return false;
    }
    // Get each attached shader one by one and detach and delete it
    GLsizei count;
    // shader id
    GLuint s;
    do
    {
        // Try to get at most *1* attached shader
        glGetAttachedShaders(id, 1, &count, &s);
        GLenum err = report_gl_error(std::string(""));
        if (GL_NO_ERROR != err)
        {
            return false;
        }
        // Check that we actually got *1*
        if (count == 1)
        {
            // Detach and delete this shader
            glDetachShader(id, s);
            glDeleteShader(s);
        }
    }
    while (count > 0);
    // Now that all of the shaders are gone we can just delete the program
    glDeleteProgram(id);
    return true;
}

GLuint MeshGL::load_shader(const std::string& src, const GLenum& type)
{
    if (src.empty())
    {
        return (GLuint) 0;
    }

    GLuint s = glCreateShader(type);
    if (s == 0)
    {
        fprintf(stderr, "Error: load_shader() failed to create shader.\n");
        return 0;
    }
    // Pass shader source string
    const char* c = src.c_str();
    glShaderSource(s, 1, &c, nullptr);
    glCompileShader(s);
    // Print info log (if any)
    print_shader_info_log(s);
    return s;
}

GLenum MeshGL::report_gl_error(const std::string& id)
{
    // http://stackoverflow.com/q/28485180/148668

    // gluErrorString was deprecated
    const auto gluErrorString = [](GLenum errorCode) -> const char*
    {
        switch (errorCode)
        {
            default:
                return "unknown error code";
            case GL_NO_ERROR:
                return "no error";
            case GL_INVALID_ENUM:
                return "invalid enum";
            case GL_INVALID_VALUE:
                return "invalid value";
            case GL_INVALID_OPERATION:
                return "invalid operation";
#ifndef GL_VERSION_3_0
                case GL_STACK_OVERFLOW:
                    return "stack overflow";
                case GL_STACK_UNDERFLOW:
                    return "stack underflow";
                case GL_TABLE_TOO_LARGE:
                    return "table too large";
#endif
            case GL_OUT_OF_MEMORY:
                return "out of memory";
#ifdef GL_EXT_framebuffer_object
                case GL_INVALID_FRAMEBUFFER_OPERATION_EXT:
                    return "invalid framebuffer operation";
#endif
        }
    };

    GLenum err = glGetError();
    if (GL_NO_ERROR != err)
    {
        fprintf(stderr, "GL_ERROR: ");
        fprintf(stderr, "%s%s\n", id.c_str(), gluErrorString(err));
    }
    return err;
}

void MeshGL::print_shader_info_log(const GLuint& obj)
{
    GLint infoLogLength = 0;
    GLint charsWritten = 0;
    char* infoLog;

    // Get shader info log from opengl
    glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &infoLogLength);

    // Only print if there is something in the log
    if (infoLogLength > 0)
    {
        infoLog = (char*) malloc(infoLogLength);
        glGetShaderInfoLog(obj, infoLogLength, &charsWritten, infoLog);
        printf("%s\n", infoLog);
        ::free(infoLog);
    }
}

void MeshGL::print_program_info_log(const GLuint& obj)
{
    GLint infoLogLength = 0;
    GLint charsWritten = 0;
    char* infoLog;

    // Get program info log from opengl
    glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &infoLogLength);

    // Only print if there is something in the log
    if (infoLogLength > 0)
    {
        infoLog = (char*) malloc(infoLogLength);
        glGetProgramInfoLog(obj, infoLogLength, &charsWritten, infoLog);
        printf("%s\n", infoLog);
        ::free(infoLog);
    }
}
