// Online C++ compiler to run C++ program online
#include <iostream>
#include <string.h>
int rand();
int rseed = 0;

#define RAND_MAX ((1U << 31) - 1)

inline int rand() {
    return rseed = (rseed * 1103515245 + 12345) & RAND_MAX;
}


int main() {
    // Write C++ code here
    char values[100];
    std::cin >> values;

    return 0;
}