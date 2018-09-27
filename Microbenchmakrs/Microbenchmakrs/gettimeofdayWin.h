#ifndef gettimeofday_h
#define gettimeofday_h

#include <windows.h>

int gettimeofday(struct timeval * tp, struct timezone * tzp);
void timersub(const timeval* tvp, const timeval* uvp, timeval* vvp);

#endif