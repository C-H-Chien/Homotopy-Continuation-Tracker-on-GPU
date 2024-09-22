#ifndef ANY_CLASS_HPP
#define ANY_CLASS_HPP

#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <functional>
#include <iostream>
#include <iomanip>


class TEST_CLASS {
    

    //> Constructor
    TEST_CLASS() {};
    
    //> Destructor
    ~TEST_CLASS();

    int add(int x, int y) { return x + y; }
    int multiply(int x, int y) { return x * y; }
    int invoke(int x, int y, std::function<int(int, int)> func);

};

#endif
