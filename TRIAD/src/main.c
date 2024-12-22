#include <stdio.h>

#include <Geomag/Geomag.h>
#include <window/window.h>

int main(int argc, char** argv)
{
    float in, out;
    GLFWwindow *win;

    in = 1.502;
    printf("sine of %f:\n", in);
    out = arm_sin_f32(in);
    printf("%f.\n", out);

    windowinginit();
    win = makewindow();

    while(!glfwWindowShouldClose(win))
    {
        glfwSwapBuffers(win);
        glfwPollEvents();
    }

    glfwDestroyWindow(win);
    glfwTerminate();

    return 0;
}
