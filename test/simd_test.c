#include <stdio.h>

int main() {
    printf( "sse: %s\n", __builtin_cpu_supports("sse") ? "Yes" : "No" );
    printf( "sse2: %s\n", __builtin_cpu_supports("sse2")  ? "Yes" : "No" );
    printf( "avx: %s\n", __builtin_cpu_supports("avx") ? "Yes" : "No" );
    printf( "avx2: %s\n", __builtin_cpu_supports("avx2") ? "Yes" : "No" );
    printf( "avx512f: %s\n", __builtin_cpu_supports("avx512f") ? "Yes" : "No" );
}
