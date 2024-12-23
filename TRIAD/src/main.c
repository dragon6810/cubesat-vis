#include <stdio.h>

#include <Geomag/Geomag.h>
#include <window/window.h>

typedef float vec3_t[3];

vec3_t vec3_origin = {0, 0, 0};

void VectorCopy(vec3_t dest, vec3_t v)
{
    int i;

    for(i=0; i<3; i++)
        dest[i] = v[i];
}

void VectorAdd(vec3_t dest, vec3_t a, vec3_t b)
{
    int i;

    for(i=0; i<3; i++)
        dest[i] = a[i] + b[i];
}

void VectorSubtract(vec3_t dest, vec3_t a, vec3_t b)
{
    int i;

    for(i=0; i<3; i++)
        dest[i] = a[i] - b[i];
}

void VectorScale(vec3_t dest, vec3_t v, float s)
{
    int i;

    for(i=0; i<3; i++)
        dest[i] = v[i] * s;
}

float VectorLength(vec3_t v)
{
    int i;

    float len;

    for(i=0, len=0; i<3; i++)
        len += v[i] * v[i];

    return sqrtf(len);
}

void VectorNormalize(vec3_t dest, vec3_t v)
{
    int i;

    float len;

    len = 1.0 / VectorLength(v);

    for(i=0; i<3; i++)
        dest[i] = v[i] * len;
}

float VectorDot(vec3_t a, vec3_t b)
{
    int i;

    float dot;

    for(i=0, dot=0; i<3; i++)
        dot += a[i] * b[i];

    return dot;
}

void VectorCross(vec3_t dest, vec3_t a, vec3_t b)
{
    vec3_t tempa, tempb;

    VectorCopy(tempa, a);
    VectorCopy(tempb, b);

    dest[0] = tempa[1] * tempb[2] - tempa[2] * tempb[1];
    dest[1] = tempa[2] * tempb[0] - tempa[0] * tempb[2];
    dest[2] = tempa[0] * tempb[1] - tempa[1] * tempb[0];
}

int VectorCmp(vec3_t a, vec3_t b)
{
    int i;

    for(i=0; i<3; i++)
    {
        if(a[i] == b[i])
            continue;

        return i + 1;
    }

    return 0;
}

void drawcyl(vec3_t start, vec3_t end, vec3_t col)
{
    const int segs = 8;
    const float r = 16;

    int i, j;

    vec3_t dir, delta;
    vec3_t cur, p, vx, vy;
    vec3_t up, right;
    float theta;

    VectorSubtract(delta, end, start);
    VectorNormalize(dir, delta);
    
    VectorCopy(up, vec3_origin);
    if(dir[2] > 0.99)
        up[0] = 1;
    else
        up[2] = 1;

    glColor3f(col[0], col[1], col[2]);

    VectorCross(right, up, dir);
    VectorCross(up, right, dir);

    glBegin(GL_QUAD_STRIP);
    for(i=0; i<=segs; i++)
    {
        theta = (float) i / (float) segs * M_PI * 2;

        VectorCopy(cur, start);
        for(j=2; j--; VectorAdd(cur, cur, delta))
        {
            VectorScale(vx, right, cosf(theta) * r);
            VectorScale(vy, up, sinf(theta) * r);
            VectorCopy(p, cur);
            VectorAdd(p, p, vx);
            VectorAdd(p, p, vy);

            glVertex3f(p[0], p[1], p[2]);
        }
    }
    glEnd();
}

void drawcone(vec3_t pos, vec3_t dir, vec3_t col)
{
    const int segs = 8;
    const float len = 96;
    const float r = 48;

    int i;

    vec3_t up, right;
    vec3_t vx, vy;
    vec3_t p, n;
    float theta;

    VectorCopy(up, vec3_origin);
    if(dir[2] > 0.99)
        up[0] = 1;
    else
        up[2] = 1;

    glColor3f(col[0], col[1], col[2]);

    VectorCross(right, up, dir);
    VectorCross(up, right, dir);

    VectorScale(dir, dir, len);

    glBegin(GL_QUAD_STRIP);
    for(i=0; i<=segs; i++)
    {
        theta = (float) i / (float) segs * M_PI * 2;

        VectorScale(vx, right, cosf(theta) * r);
        VectorScale(vy, up, sinf(theta) * r);

        VectorAdd(p, vx, vy);
        VectorAdd(p, p, pos);
        glVertex3f(p[0], p[1], p[2]);

        VectorAdd(p, pos, dir);
        glVertex3f(p[0], p[1], p[2]);
    }
    glEnd();
}

void drawarrow(vec3_t start, vec3_t end, vec3_t col)
{
    vec3_t dir;

    drawcyl(start, end, col);
    
    VectorSubtract(dir, end, start);
    VectorNormalize(dir, dir);
    drawcone(end, dir, col);
}

void drawball(vec3_t pos, float r, vec3_t col)
{
    const int usegs = 32;
    const int vsegs = 32;
    const float ambient = 0.5;

    int i, j;
    float v0, v1;
    float u;
    vec3_t n;
    float mul;

    glPushMatrix();
    glTranslatef(pos[0], pos[1], pos[2]);
    glColor3f(col[0], col[1], col[2]);
    
    for(i=0; i<vsegs; i++)
    {
        v0 = (float) i / (float) vsegs * M_PI;
        v1 = (float) (i+1) / (float) vsegs * M_PI;
    
        glBegin(GL_QUAD_STRIP);
        for(j=0; j<=usegs; j++)
        {
            u = (float)j / (float) usegs * M_PI * 2;

            mul = -cos(v0);
            if(mul < 0)
                mul = 0;
            mul = mul + ambient;
            glColor3f(col[0] * mul, col[1] * mul, col[2] * mul);
            glVertex3f(r * sin(v0) * cos(u), r * sin(v0) * sin(u), r * -cos(v0));
            
            mul = -cos(v1);
            if(mul < 0)
                mul = 0;
            mul = mul + ambient;
            glColor3f(col[0] * mul, col[1] * mul, col[2] * mul);
            glVertex3f(r * sin(v1) * cos(u), r * sin(v1) * sin(u), r * -cos(v1));
        }
        glEnd();
    }

    glPopMatrix();
}

void render(void)
{
    vec3_t col;
    vec3_t end;

    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 16.0 / 9.0, 1, 10000);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt
    (
        4096, 4096, 4096,
        0, 0, 0,
        0, 0, 1
    );

    glColor3f(1, 1, 1);
    glBegin(GL_QUADS);
    glVertex3f(1024, 1024, 0);
    glVertex3f(-1024, 1024, 0);
    glVertex3f(-1024, -1024, 0);
    glVertex3f(1024, -1024, 0);
    glEnd();

    col[0] = 1;
    col[1] = 0;
    col[2] = 0;

    drawball(vec3_origin, 1024, col);

    end[0] = 0;
    end[1] = 0;
    end[2] = 2048;
    col[1] = 1;
    drawarrow(vec3_origin, end, col);
}

int main(int argc, char** argv)
{
    GLFWwindow *win;

    windowinginit();
    win = makewindow();

    while(!glfwWindowShouldClose(win))
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        render();

        glfwSwapBuffers(win);
        glfwPollEvents();
    }

    glfwDestroyWindow(win);
    glfwTerminate();

    return 0;
}
