#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32
  #define STD_SRANDOM srand
  #define STD_RANDOM  rand
  #define RANDOMINIT_NAME RANDOMINIT
  #define RANDOM_NAME RANDOM
  #define FCALL __stdcall
#else
  #define STD_SRANDOM srandom
  #define STD_RANDOM random
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
