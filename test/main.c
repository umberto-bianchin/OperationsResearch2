#include <stdio.h>
#include <cplex.h>

int main(int argc, char** argv)
{
    printf("CPX_INFBOUND: %f\n", CPX_INFBOUND);
    return 0;
}