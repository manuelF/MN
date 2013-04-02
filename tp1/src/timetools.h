#ifndef _UTILS
#define _UTILS

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <sys/time.h>
#define TIMEUNIT timespec
const timespec diff(const timespec& start, const timespec& end); //saca la dif de tiempo
#include "timetools.cpp"

#endif