
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

class Shader
{
public:
  struct ShaderProgramSource
  {
    std::string vertex_shader_source;
    std::string fragment_shader_source;
  };

public:
  Shader(const std::string & file_path)
  {
    ShaderProgramSource source = ParseShader(file_path);
    ID = CreateShader(source.vertex_shader_source, source.fragment_shader_source);
  }

  void use() 
  { 
      glUseProgram(ID); 
  }

  void setBool(const std::string &name, bool value) const
  {         
      glUniform1i(glGetUniformLocation(ID, name.c_str()), (int)value); 
  }

  void setInt(const std::string &name, int value) const
  { 
      glUniform1i(glGetUniformLocation(ID, name.c_str()), value); 
  }

  void setFloat(const std::string &name, float value) const
  { 
      glUniform1f(glGetUniformLocation(ID, name.c_str()), value); 
  }

  ShaderProgramSource ParseShader(const std::string & file_path)
  {
    enum class ShaderType
    {
      NONE = -1, VERTEX = 0, FRAGMENT=1
    };

    std::ifstream stream(file_path);
    std::stringstream ss[2];
    std::string line;
    ShaderType type;
    while(getline(stream, line))
    {
      if(line.find("#shader") != std::string::npos)
      {
        if(line.find("vertex")!=std::string::npos)
          type = ShaderType::VERTEX;
        if(line.find("fragment")!=std::string::npos)
          type = ShaderType::FRAGMENT;
      }
      else
      {
        ss[(uint8_t)type] << line << "\n";
      }
    }
    return {ss[0].str(), ss[1].str()};
  }

  static uint32_t CompileShader(uint32_t type, const std::string & source)
  {
    uint32_t shader = glCreateShader(type);
    const char * src = source.c_str();
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);

    /** 检查是否有编译错误 */
    int success;
    char infoLog[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
      glGetShaderInfoLog(shader, 512, NULL, infoLog);
      std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    return shader;
  }

  static uint32_t CreateShader(const std::string & vshader, const std::string & fshader)
  {
    uint32_t program = glCreateProgram();
    uint32_t vs = CompileShader(GL_VERTEX_SHADER, vshader);
    uint32_t fs = CompileShader(GL_FRAGMENT_SHADER, fshader);

    glAttachShader(program, vs);
    glAttachShader(program, fs);

    glLinkProgram(program);
    glValidateProgram(program);

    glDeleteShader(vs);
    glDeleteShader(fs);
    return program;
  }
public:
  uint32_t ID;
};


class View
{
public:
  View() {}

  GLFWwindow * init(uint32_t width = 800, uint32_t height = 600)
  {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); /**< 告诉 GLFW GL 的主版本号 */
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3); /**< 告诉 GLFW GL 的次版本号 */
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    /** 创建一个窗口，800 宽 600 高，名字叫 LeanOpenGL */
    GLFWwindow* window = glfwCreateWindow(width, height, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
      std::cout << "Failed to create GLFW window" << std::endl;
      glfwTerminate();
      return nullptr;
    }
    glfwMakeContextCurrent(window);

    /** 告诉 GLAD OpenGL 的函数指针位置 */
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
      std::cout << "Failed to initialize GLAD" << std::endl;
      return nullptr;
    }
    return window;
  }
};

