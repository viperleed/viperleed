#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* for use with search.v103 and later */

void randominit_(unsigned int *Init)
{
  time_t t;
  if (*Init==0)
    srandom((unsigned) time(&t));
  else
    srandom(*Init);
}

int random_(void)
{return random()%1000;}
