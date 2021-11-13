#ifndef TIMER_H
#define TIMER_H

#include <cstdlib>
#include <cstdint>

#if defined(__WIN32__) || defined(_WIN32) || defined(WIN32) || defined(__WINDOWS__) || defined(__TOS_WIN__)
#include <winsock2.h>
#include <windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif

#if defined(__WIN32__) || defined(_WIN32) || defined(WIN32) || defined(__WINDOWS__) || defined(__TOS_WIN__)
inline int gettimeofday(struct timeval* tp, struct timezone* tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    // This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
    // until 00:00:00 January 1, 1970 
    static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    time = ((uint64_t)file_time.dwLowDateTime);
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec = (long)((time - EPOCH) / 10000000L);
    tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
    return 0;
}
#endif // gettimeofday

#ifndef timersub
#define timersub(a, b, result) \
        do { \
                (result)->tv_sec = (a)->tv_sec - (b)->tv_sec; \
                (result)->tv_usec = (a)->tv_usec - (b)->tv_usec; \
                if ((result)->tv_usec < 0) { \
                        --(result)->tv_sec; \
                        (result)->tv_usec += 1000000; \
                } \
        } while (0)
#endif // timersub

namespace timer
{
    inline timeval StartTimer()
    {
        struct timeval timerStart{};
        gettimeofday(&timerStart, nullptr);
        return timerStart;
    }

    // time elapsed in ms
    inline double GetTimer(timeval timerStart)
    {
        struct timeval timerStop{}, timerElapsed{};
        gettimeofday(&timerStop, nullptr);
        timersub(&timerStop, &timerStart, &timerElapsed);
        return timerElapsed.tv_sec * 1000.0 + timerElapsed.tv_usec / 1000.0;
    }

#if defined(__WIN32__) || defined(_WIN32) || defined(WIN32) || defined(__WINDOWS__) || defined(__TOS_WIN__)

    inline void delay( unsigned long ms )
    {
        Sleep( ms );
    }

#else  /* presume POSIX */

    inline void delay(unsigned long ms)
    {
        usleep(ms * 1000);
    }

#endif
}

#endif // TIMER_H
