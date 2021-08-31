#include <iostream>

int main() {

    // allocate on the heap -- valgrind does bounds checking here
    double *a = new double[10];

    // allocate on the stack -- valgrind does not do stack bounds checking
    double b[10];

    // valgrind catches
    a[10] = 1.0;

    // valgrind does not catch
    b[10] = 1.0;

    // valgrind catches
    //std::cout << b[0] << std::endl;

    return 1;

}
