#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <fstream>
#include <sstream>

struct ShaderProgramSource
{
  std::string vertex_shader_source;
  std::string fragment_shader_source;
};

static ShaderProgramSource ParseShader(const std::string & file_path)
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

void framebuffer_size_callback(GLFWwindow * , int width, int height)
{
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void draw_a_triangle()
{
}

int main(int , char ** argv)
{
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); /**< 告诉 GLFW GL 的主版本号 */
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3); /**< 告诉 GLFW GL 的次版本号 */
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  /** 创建一个窗口，800 宽 600 高，名字叫 LeanOpenGL */
  GLFWwindow* window = glfwCreateWindow(800, 600, "LearnOpenGL", NULL, NULL);
  if (window == NULL)
  {
    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window);

  /** 告诉 GLAD OpenGL 的函数指针位置 */
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  ShaderProgramSource source = ParseShader(argv[1]);
  uint32_t shader = CreateShader(source.vertex_shader_source, source.fragment_shader_source);
  glUseProgram(shader);

  /** set up vertex data (and buffer(s)) and configure vertex attributes */
  float vertices[] = {
      -0.5f, -0.5f, 0.0f, // left  
       0.5f, -0.5f, 0.0f, // right 
       0.5f,  0.5f, 0.0f,  // top   
      -0.5f,  0.5f, 0.0f  // top   
  }; 

  unsigned int VBO, VAO;
  // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
  glGenVertexArrays(1, &VAO);
  glBindVertexArray(GL_ELEMENT_ARRAY_BUFFER, VAO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, VAO)

  glGenBuffers(1, &VBO);

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);


  /** 定义视图窗口 */
  glViewport(0, 0, 800, 600);

  /** 定义回调函数 */
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

  /** 循环渲染 */
  while(!glfwWindowShouldClose(window))
  {
    processInput(window);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    /** 绘制窗口 */
    glfwSwapBuffers(window);

    /** 检测有没有发生鼠标移动，键盘输入等事件 */
    glfwPollEvents();
  }
  glDeleteProgram(shader);
  glfwTerminate();
  return 0;
}


