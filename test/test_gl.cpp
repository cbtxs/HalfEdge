#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "GL/view.h"

void framebuffer_size_callback(GLFWwindow * , int width, int height)
{
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void draw_a_triangle(Shader & shader)
{
  float timeValue = glfwGetTime();
  float greenValue = (std::sin(timeValue) / 2.0f) + 0.5f;
  int vertexColorLocation = glGetUniformLocation(shader.ID, "ourColor");
  glUniform4f(vertexColorLocation, 0.0f, greenValue, 0.0f, 1.0f);
  glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, 0);
}

int main(int , char ** argv)
{
  View view;
  GLFWwindow * window = view.init();
  if(!window)
  {
    return -1;
  }

  /** set up vertex data (and buffer(s)) and configure vertex attributes */
  float vertices[] = {
    // 位置              // 颜色
     0.5f, -0.5f, 0.0f,  1.0f, 0.0f, 0.0f,   // 右下
    -0.5f, -0.5f, 0.0f,  0.0f, 1.0f, 0.0f,   // 左下
     0.0f,  0.5f, 0.0f,  0.0f, 0.0f, 1.0f    // 顶部
  };

  uint32_t indices[] = 
  {
    0, 1, 2,
    0, 2, 3
  };

  /** 定义顶点缓冲数组 */
  unsigned int VBO, VAO;
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);

  glBindVertexArray(VAO);

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

  // 位置属性
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);

  // 颜色属性
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3* sizeof(float)));
  glEnableVertexAttribArray(1);

  /** 定义单元缓冲数组 */
  uint32_t EBO;
  glGenBuffers(1, &EBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

  Shader shader(argv[1]);
  shader.use();

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

    shader.use();

    glBindVertexArray(VAO);
    //glDrawArrays(GL_TRIANGLES, 0, 3);
    
    draw_a_triangle(shader);

    /** 绘制窗口 */
    glfwSwapBuffers(window);

    /** 检测有没有发生鼠标移动，键盘输入等事件 */
    glfwPollEvents();
  }
  glDeleteProgram(shader.ID);
  glfwTerminate();
  return 0;
}


