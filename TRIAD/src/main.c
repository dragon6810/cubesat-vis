#include <stdio.h>
#include <assert.h>
#include <stdbool.h>

#include <VecUtils.h>

#include <Geomag/Geomag.h>
#include <window/window.h>
#include <CoordinateConversions/CoordinateConversions.h>

typedef float vec3_t[3];

vec3_t vec3_origin = {0, 0, 0};

GLFWwindow *win;

float sampleelevation = 0;

float camtheta = 45, camphi = 45;
const float camdistance = 32768;

float lastmousex, lastmousey;
float mousex, mousey;
const float sensitivity = 128.0;
bool mouseoverui;
float windowaspect;

int dontdrawmag = 1;

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

void drawline(vec3_t start, vec3_t end, vec3_t col)
{
    const float stroke = 16;

    glLineWidth(stroke);
    glColor3f(col[0], col[1], col[2]);
    
    glBegin(GL_LINES);
    glVertex3f(start[0], start[1], start[2]);
    glVertex3f(end[0], end[1], end[2]);
    glEnd();
    
    glColor3f(1, 1, 1);
    glLineWidth(1);
}

void drawring(vec3_t center, vec3_t a, vec3_t b, vec3_t col, float start, float end)
{
    const int segs = 64;
    const float stroke = 16;

    int i, j;

    vec3_t cura, curb, p;
    float theta;

    glColor3f(col[0], col[1], col[2]);

    glLineWidth(stroke);
    glBegin(GL_LINES);

    for(i=0; i<segs; i++)
    {
        for(j=0; j<2; j++)
        {
            theta = M_PI * 2.0 / (float) segs * (float) (i + j);
            theta = theta * (end - start) / (M_PI * 2.0) + start;

            VectorScale(cura, a, cosf(theta));
            VectorScale(curb, b, sinf(theta));

            VectorCopy(p, center);
            VectorAdd(p, p, cura);
            VectorAdd(p, p, curb);

            glVertex3f(p[0], p[1], p[2]);
        }
    }

    glEnd();
    glLineWidth(1);

    glColor3f(1, 1, 1);
}

void drawcyl(vec3_t start, vec3_t end, vec3_t col)
{
    const int segs = 8;
    const float r = 128;
    const float ambient = 0.5;

    int i, j;

    vec3_t dir, delta;
    vec3_t cur, p, vx, vy;
    vec3_t up, right;
    vec3_t n;
    float theta;
    float mul;

    VectorSubtract(delta, end, start);
    VectorNormalize(dir, delta);
    
    VectorCopy(up, vec3_origin);
    if(dir[2] > 0.99)
        up[0] = 1;
    else
        up[2] = 1;

    glColor3f(col[0], col[1], col[2]);

    VectorCross(right, up, dir);
    VectorNormalize(right, right);
    VectorCross(up, right, dir);

    glBegin(GL_QUAD_STRIP);
    for(i=0; i<=segs; i++)
    {
        theta = (float) i / (float) segs * M_PI * 2;

        VectorCopy(cur, start);
        for(j=2; j--; VectorAdd(cur, cur, delta))
        {
            VectorScale(vx, right, cosf(theta));
            VectorScale(vy, up, sinf(theta));
            VectorAdd(n, vx, vy);
            mul = n[2];
            mul += ambient;
            glColor3f(col[0] * mul, col[1] * mul, col[2] * mul);

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
    const float len = 512;
    const float r = 256;
    const float ambient = 0.5;

    int i;

    vec3_t up, right;
    vec3_t vx, vy;
    vec3_t p, n;
    float theta;
    float mul;

    VectorCopy(up, vec3_origin);
    if(dir[2] > 0.99)
        up[0] = 1;
    else
        up[2] = 1;

    glColor3f(col[0], col[1], col[2]);

    VectorCross(right, up, dir);
    VectorNormalize(right, right);
    VectorCross(up, right, dir);

    VectorScale(dir, dir, len);

    glBegin(GL_QUAD_STRIP);
    for(i=0; i<=segs; i++)
    {
        theta = (float) i / (float) segs * M_PI * 2;

        VectorScale(vx, right, cosf(theta));
        VectorScale(vy, up, sinf(theta));
        VectorAdd(n, vx, vy);
        VectorScale(dir, dir, 1/len);
        VectorAdd(n, n, dir);
        VectorScale(dir, dir, len);
        VectorNormalize(n, n);

        mul = n[2];
        if(mul < 0)
            mul = 0;
        mul += ambient;
        glColor3f(col[0] * mul, col[1] * mul, col[2] * mul);

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

void drawmagentry(vec3_t entry, vec3_t pos)
{
    const float vecscale = 0.02;

    vec3_t start, end;
    vec3_t scaled;
    vec3_t col;

    col[0] = 1;
    col[1] = 0.5;
    col[2] = 0;

    VectorCopy(start, pos);
    VectorScale(scaled, entry, vecscale);
    VectorAdd(end, scaled, pos);

    drawarrow(start, end, col); 
}

// Generates and draws samples in fibonacci sphere around origin with distance d
void drawmagfield(void)
{
    const int nsamples = 512;
    const float phi = M_PI * (float) (3.0 - sqrtf(5));
    const float d = 6371 + sampleelevation;

    int i;

    vec3_t in, out;
    Vec3D_t vin;
    Vec3D_t vout;
    time_t t;
    float theta;
    float y, r;

    time(&t);
    for(i=0; i<nsamples; i++)
    {
        y = 1.0 - ((float) i / (float) (nsamples - 1)) * 2;
        r = sqrtf(1.0 - y * y);

        theta = (float) i * phi;

        in[0] = d * r * cosf(theta);
        in[1] = d * r * sinf(theta);
        in[2] = d * y;

        Coord_ECEFToECI(t, in, in);

        vin.X = in[0];
        vin.Y = in[1];
        vin.Z = in[2];

        assert(Geomag_GetMagEquatorial(&t, &vin, &vout) == FR_OK);
        out[0] = vout.X;
        out[1] = vout.Y;
        out[2] = vout.Z;

        drawmagentry(out, in);
    }
}

void updatecam(void)
{
    float deltax, deltay;

    if(mouseoverui)
        return;

    if(glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_1) == GLFW_RELEASE)
        return;

    deltax = (mousex - lastmousex) * windowaspect;
    deltay = mousey - lastmousey;

    camtheta -= deltax * sensitivity;
    camphi += deltay * sensitivity;

    if(camphi > 89)
        camphi = 89;
    if(camphi < -89)
        camphi = -89;
}

struct nk_context *nukcontext;

void render(void)
{
    vec3_t col;
    vec3_t end;
    vec3_t ringa, ringb;
    vec3_t lnstart, lnend;
    vec3_t campos;

    nk_glfw3_new_frame();

    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 16.0 / 9.0, 10, 100000);
   
    updatecam();

    lastmousex = mousex;
    lastmousey = mousey;

    campos[0] = cosf(DEG2RAD(camtheta)) * cosf(DEG2RAD(camphi)) * camdistance;
    campos[1] = sinf(DEG2RAD(camtheta)) * cosf(DEG2RAD(camphi)) * camdistance;
    campos[2] = sinf(DEG2RAD(camphi)) * camdistance;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt
    (
        campos[0], campos[1], campos[2],
        0, 0, 0,
        0, 0, 1
    );

    col[0] = 0.5;
    col[1] = 0.7;
    col[2] = 1.0;

    drawball(vec3_origin, 6371, col);

    // Equator

    col[0] = 1.0;
    col[1] = 0.0;
    col[2] = 0.7;
    VectorCopy(ringa, vec3_origin);
    VectorCopy(ringb, vec3_origin);
    ringa[0] = 6500;
    ringb[1] = 6500;
    drawring(vec3_origin, ringa, ringb, col, 0, M_PI * 2.0);

    // Prime Meridian

    col[0] = 0.1;
    col[1] = 0.6;
    col[2] = 0.2;
    VectorCopy(ringa, vec3_origin);
    VectorCopy(ringb, vec3_origin);
    ringa[0] = 6500;
    ringb[2] = 6500;
    Coord_ECEFToECI(time(NULL), ringa, ringa);
    drawring(vec3_origin, ringa, ringb, col, -M_PI_2, M_PI_2);

    // Earth's pole
    col[0] = 0.2;
    col[1] = 0.2;
    col[2] = 1.0;
    VectorCopy(lnstart, vec3_origin);
    VectorCopy(lnend, vec3_origin);
    lnstart[2] = -12000;
    lnend[2] = 12000;
    drawline(lnstart, lnend, col);

    // X Axis (ECI)
    col[0] = 1.0;
    col[1] = 0.0;
    col[2] = 0.0;
    VectorCopy(lnstart, vec3_origin);
    VectorCopy(lnend, vec3_origin);
    lnend[0] = 12000;
    drawline(lnstart, lnend, col);

    // Y Axis (ECI)
    col[0] = 0.0;
    col[1] = 1.0;
    col[2] = 0.0;
    VectorCopy(lnstart, vec3_origin);
    VectorCopy(lnend, vec3_origin);
    lnend[1] = 12000;
    drawline(lnstart, lnend, col);

    if(!dontdrawmag)
        drawmagfield();

    if (nk_begin(nukcontext, "mag vis", nk_rect(32, 32, 280, 300),
                     NK_WINDOW_BORDER)) 
    {
        nk_layout_row_static(nukcontext, 30, 256, 1);
        
        nk_label(nukcontext, "magnetic field altitude in km:", NK_TEXT_LEFT);
        nk_slider_float(nukcontext, 0, &sampleelevation, 8192, 16);

        nk_checkbox_label(nukcontext, "draw magnetic field", &dontdrawmag);

        mouseoverui = nk_window_is_hovered(nukcontext) || nk_item_is_any_active(nukcontext);
    }
    nk_end(nukcontext);

    nk_glfw3_render(NK_ANTI_ALIASING_OFF);
}

static void cursorposcallback(GLFWwindow* win, double x, double y)
{
    int w, h;

    glfwGetWindowSize(win, &w, &h);

    mousex = x / (float) w;
    mousey = y / (float) h;
}

int main(int argc, char** argv)
{
    int w, h;
    Vec3D_t in, out, expect;
    struct nk_font_atlas *atlas;
    struct nk_font *font;

    printf("\n================================ TRIAD Visualization/Testing ================================\n\n");

    printf("running tests...\n");

    Coord_TestConversions();
    Geomag_RunTests("WMM2025_TestValues.txt");

    printf("all tests passed.\n\n");

    windowinginit();
    win = makewindow();

    glfwGetWindowSize(win, &w, &h);
    windowaspect = (float) w / (float) h;

    glfwSetCursorPosCallback(win, cursorposcallback);

    nukcontext = nk_glfw3_init(win, NK_GLFW3_INSTALL_CALLBACKS);
    nk_glfw3_font_stash_begin(&atlas);
    font = nk_font_atlas_add_default(atlas, 14, 0);
    nk_glfw3_font_stash_end();
    nukcontext->style.font = &font->handle;

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
