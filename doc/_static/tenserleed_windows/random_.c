#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* The C library: srandom/random are POSIX, not in the Windows CRT */
#ifdef _WIN32
  #define STD_SRANDOM srand
  #define STD_RANDOM  rand
#else
  #define STD_SRANDOM srandom
  #define STD_RANDOM random
#endif

/* The Fortran compiler: Intel Fortran on Windows expects upper-case
   names and __stdcall, gfortran (also on Windows) lower-case names
   with a trailing underscore. Check the C compiler, not the platform:
   gcc goes with gfortran, icx/cl with Intel Fortran. */
#if defined(_WIN32) && !defined(__GNUC__)
  #define RANDOMINIT_NAME RANDOMINIT
  #define RANDOM_NAME RANDOM
  #define FCALL __stdcall
#else
  #define RANDOMINIT_NAME randominit_
  #define RANDOM_NAME random_
  #define FCALL
#endif

void FCALL RANDOMINIT_NAME(unsigned int *Init)
{
  time_t t;
  if (*Init==0)
    STD_SRANDOM((unsigned) time(&t));
  else
    STD_SRANDOM(*Init);
}

int FCALL RANDOM_NAME(void)
{return STD_RANDOM()%1000;}
