#include <stdio.h>

#include <arm_math.h>

int main(int argc, char** argv)
{
    float in, out;

    in = 1.502;
    printf("sine of %f:\n", in);
    out = arm_sin_f32(in);
    printf("%f.\n", out);

    return 0;
}
