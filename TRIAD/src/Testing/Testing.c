#include <Testing/Testing.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define SBLD  "\033[1m"
#define SNRM  "\033[0m"

static int Testing_VectorCmp(Vec3D_t a, Vec3D_t b, float epsilon)
{
    int i;

    for(i=0; i<3; i++)
    {
        if(fabsf(a.Vec[i] - b.Vec[i]) > epsilon)
            return false;
    }

    return true;
}

void Testing_TestVectorImpl(Vec3D_t expect, Vec3D_t in, Vec3D_t out, float epsilon, const char* test, const char* file, int line)
{
    if(Testing_VectorCmp(out, expect, epsilon))
        return;

    printf("\n%s================ TEST FAILURE ================\n\n", KYEL);
    
    printf("%s%s %sFailed at %s%s:%d:\n\n", KGRN, test, KRED, KGRN, file, line);
    
    printf("%sInput value:     { %f, %f, %f }.\n", KGRN, in.X, in.Y, in.Z);
    printf("%sExpected output: { %f, %f, %f }.\n", KGRN, expect.X, expect.Y, expect.Z);
    printf("%sOutput value:    { %f, %f, %f }.\n", KYEL, out.X, out.Y, out.Z);
    
    puts(KNRM);
    abort();
}