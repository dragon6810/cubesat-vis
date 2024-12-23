#include <window/window.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

GLFWwindow* makewindow(void)
{
    GLFWwindow *win;
    GLenum err;

    assert(win = glfwCreateWindow(1280, 720, "TRIAD Vis", NULL, NULL));

    glfwMakeContextCurrent(win);

    err = glewInit();
    if(err != GLEW_OK)
    {
        printf("failed to init glew: \"%s\".\n", glewGetErrorString(err));
        abort();
    }

    glViewport(0, 0, 1280, 720);

    return win;
}
