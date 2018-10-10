#ifndef gettimeofday_h
#define gettimeofday_h

// MSVC defines this in winsock2.h!?
typedef struct timeval {
    long tv_sec;
    long tv_usec;
} timeval;

int gettimeofday(struct timeval *tp, struct timezone *tzp);

#endif