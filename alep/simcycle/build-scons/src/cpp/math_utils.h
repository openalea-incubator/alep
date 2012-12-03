//                       Fonctions math complementaires




#ifndef _MATH_UTILS_H
#define _MATH_UTILS_H

#include <math.h>
#include <stdlib.h> //rand

// arrondi un nombre de facon probabiliste : la partie decimale donnant la proba d'arrodir a l'entier superieur

double roundP(double x) {
  double N = floor(x);
  double p = x - N;
  double r = rand() / (double) RAND_MAX;
  if (r <= p)
    N += 1.;
  return(N);
}

//interpolateur lineaire. 
//Si extrapolate = 1, alors extrapole au dela de la derniere valeur avec la derniere pente

double approx_ex(double *x,double *y,int Nbp,double xout,int extrapolate) {
  int i=0;
  double yout;
  double extrate = extrapolate ? (y[Nbp -1] - y[Nbp - 2]) / (x[Nbp - 1] - x[Nbp - 2]) : 0;

  if (xout <= x[0])
    yout = y[0];
  else if (xout >=x [Nbp - 1])
    yout = y[Nbp - 1] + extrate * (xout - x[Nbp -1]);
  else {
    while(xout >= x[i])
      i++;
    //show("\npos in Tc=%d\n",i);
    yout = y[i - 1] + (y[i] - y[i-1]) / (x[i] - x[i - 1]) * (xout - x[i - 1]);
  }

  return(yout);
}

double approx(double x[],double y[],int Nbp,double xout) {
  return approx_ex(x,y,Nbp,xout,0);
}

//interpolateur simple avec 2 points

float ylin(float xout,float x1,float x2,float y1,float y2) {
	float yout;
	
  if (xout<=x1)
    yout=y1;
  else if (xout>=x2)
    yout=y2;
  else
    yout=y1+(y2-y1)/(x2-x1)*(xout-x1);

  return(yout);
}


#endif
