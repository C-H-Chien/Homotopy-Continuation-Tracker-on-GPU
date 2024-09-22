#ifndef ANY_CLASS_CPP
#define ANY_CLASS_CPP

#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <functional>
#include <iostream>
#include <iomanip>
#include "any_class.hpp"

int TEST_CLASS::invoke(int x, int y, std::function<int(int, int)> func) {
    return func(x, y);
}

#endif
