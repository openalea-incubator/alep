// General symbol definition for adel SEPTO 


#ifndef _SYMBOLS_H
#define _SYMBOLS_H

#include <stdio.h> //printf

#ifdef LPFGPRINTF
#define show Printf
#else
#define show printf
#endif

#define Centi 0.01
#define Meter 1.
#define hour (1. / 24.)
#define epsillon 0.000005

#define SQUARE(x) ((x)*(x))
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define ABS(A) ((A) >= 0 ? (A) : -(A))
#define isEqual(A,B) ( ABS((A) - (B)) < epsillon )
#define zapSmall(A) ( isEqual(A,0) ? 0 : (A)  )
#define PI 3.14116
#define rad(deg) ((deg) * PI / 180.)
#define deg(rad) ((rad) / PI * 180.)


#define TABLENGTH(tab) (int) (sizeof tab / sizeof *tab)

#endif
