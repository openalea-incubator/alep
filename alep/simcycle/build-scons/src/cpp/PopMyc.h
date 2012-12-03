//Structures et Fonctions pour le stockage et la manipulation de population de champignon
//Une pop est decrite par des classes  d'age (de 0 a nbClass-1)
//Pour chaque classe d'age on a une quantité Q (ici des cm2) et un effectif N

//utilisable avec deniere classe plus courte (pour simuler sortie au bon moment). dans ce cas, on peut avoir une 'retention' artefactuelle dans la derniere classe pour la classe venant d'y arriver si le time buffer > largeur deniere classe

#ifndef _POPMYC_H
#define _POPMYC_H


#include <stdio.h>
#include <stdlib.h> //rand,calloc
#include <string.h> //strcat
#include <math.h>

#include "symbols.h"

// PopMyc types (for printing)
enum popTypes {pTLes,pTUdin,pTSpo,pTHlat};

//structure generique d'une pop

typedef struct POPMYC {
  //variables d'etat
  int type;//Type de pop (defini N et Q)
  int nbClass;//nombre de classes. Toujours >=1, meme si agemin = agemax.
  float agemin;//borne inferieure des classes d'age de la pop. Par convention, une classe est du type [ageinf,agesup[
  float agemax;//age reel auquel sortent les classes de la pop, peut etre inferieur à agemin+nbclass*largclass quand la duree de residence dans la pop ne tombe pas sur un multiple de largclass
  float *Q;//pointeur sur les surfaces des classes (attention!! : doit pointer sur une zone de nbClass float)
  float *N;//pointeur sur les effectifs des classes (attention!! : doit pointer sur une zone de nbClass float)
  float *decT;//decalage temporel de la classe par rapport au TimeBuffer (lors de l'arrivee de la classe). Est utilisé pour garantir temps de residence (decT est dans l'intervalle [0,largClass[).
  float largClass;//"largeur" d'une classe d'age
  //variable pour la gestion de la structure
  int pC0;//position dans le tableau Q de la classe 0 (autre classes sont à sa droite)
  float TimeBuffer;//Tampon de temps gerant les moment ou se produisent les changements de classe (Time Buffer est dans l'intervalle [0,largClass[)
} PopMyc;

//structure de description d'une classe (Q,age,N)

typedef struct CLASSMYC {
  float Q;
  float age;
  float N;
} ClassMyc;

//structure de description d'une baree dans la pop

typedef struct BARPOP {
  float Q;
  float decT;
  int pos;
  float N;
} BarPop;

//Piles de priorites/transformation pour classMyc/BarPop. ilim donne l'indice de la dernier classe tranformee

typedef struct STACKCMYC {
  ClassMyc *stack;
  int nclass;
  int ilim;
} StackCmyc;

typedef struct STACKBPOP {
  BarPop *stack;
  int nclass;
  int ilim;
} StackBpop;

// fonction pour l'initialisation

ClassMyc newCmyc_p(float age,float Q,float N) {
  ClassMyc Cmyc;
  Cmyc.age = age;
  Cmyc.Q = Q;
  Cmyc.N = N;
  return(Cmyc);
}

ClassMyc newCmyc(void) {
  ClassMyc Cmyc;
  Cmyc.age = 0;
  Cmyc.Q = 0;
  Cmyc.N = 0;
  return(Cmyc);
}

void ResetCmyc(ClassMyc *Cmyc) {
  Cmyc->age = 0;
  Cmyc->Q = 0;
  Cmyc->N = 0;
}

PopMyc newPopMyc(float agemin,float agemax,float largClass,int popType) {
  PopMyc pop;
  pop.type = popType;
  pop.nbClass = (int) ceil((agemax - agemin) / largClass);
  if (pop.nbClass < 1)
    pop.nbClass = 1;
  pop.agemin = agemin;
  pop.agemax = agemax;
  pop.Q = (float*) calloc(pop.nbClass, sizeof(float));
  pop.N = (float*) calloc(pop.nbClass, sizeof(float));
  pop.decT = (float*) calloc(pop.nbClass, sizeof(float));
  pop.largClass = largClass;
  pop.pC0 = 0;//position dans Q de la classe 0
  pop.TimeBuffer = 0;
  return(pop);
};

void copyPopMyc(PopMyc *source, PopMyc* dest) {
  if (source->nbClass != dest->nbClass) {
    show("<PopMyc.h> copyPopMyc : Mismatch in populations dimensions : nothing copied !");
    return;
  }

  dest->type = source->type;
  dest->nbClass = source->nbClass;
  dest->agemin = source->agemin;
  dest->agemax = source->agemax;
  dest->largClass = source->largClass;
  dest->pC0 = source->pC0;
  dest->TimeBuffer = source->TimeBuffer;
  int i;
  for (i = 0; i < source->nbClass; i++) {
    dest->Q[i] = source->Q[i];
    dest->N[i] = source->N[i];
    dest->decT[i] = source->decT[i];
  }

  return;
}

//vide une pop de ses effectifs
void ResetPop(PopMyc *pop) {
  int i;
  for(i = 0; i < pop->nbClass; i++) {
    pop->Q[i] = 0;
    pop->N[i] = 0;
    pop->decT[i] = 0;
  }
  pop->pC0 = 0;
  pop->TimeBuffer = 0;
}

int posC(int numC,PopMyc *pop);
//Remet la classe numc a zero
void ResetClass(int numC,PopMyc *pop) {
  int pos = posC(numC, pop);
  pop->Q[pos] = 0;
  pop->N[pos] = 0;
  pop->decT[pos] = 0;
}

//Remet la classe en position pos a zero
void ResetBpop(int pos,PopMyc *pop) {
  pop->Q[pos] = 0;
  pop->N[pos] = 0;
  pop->decT[pos] = 0;
}


StackCmyc newStackCmyc(int size) {
  StackCmyc s;
  s.stack = (ClassMyc*) malloc(size * sizeof(ClassMyc));
  s.nclass = 0;
  s.ilim = -1;
  return(s);
}

StackBpop newStackBpop(int size) {
  StackBpop s;
  s.stack = (BarPop*) malloc(size * sizeof(BarPop));
  s.nclass = 0;
  s.ilim = -1;
  return(s);
}


//fonction pour liberer les pointeurs

void freePop(PopMyc *pop) {
  free(pop->Q);
  free(pop->N);
  free(pop->decT);
}

void freeStackCmyc(StackCmyc *s) {
  free(s->stack);
}

void freeStackBpop(StackBpop *s) {
  free(s->stack);
}


//fonction qui retourne la position de la classe numC dans les tableaux d'une population pop, -1 si pas dans la pop 

int posC(int numC,PopMyc *pop) {
  int p;
  if (numC >= pop->nbClass || numC < 0) {
    show("\n<PopMyc.h> : posC :erreur recherche d'une classe qui n'est pas dans la pop\n");
    p = -1;
  }
  else {
    p = pop->pC0 + numC;
    p = (p > (pop->nbClass - 1)) ? (p - pop->nbClass) : p;
  }
  return(p);
}

//fonction qui retourne le numero de la classe incluant l'age dans population pop, -1 si pas dans la pop 

int numCAge(float age,PopMyc *pop) {
  int numC;
  if ( (pop->agemax - pop->agemin) < pop->largClass && age < pop->largClass )
    numC = 0;
  else if (age >= pop->agemax || (age < pop->agemin)) {
    show("<PopMyc.h> : numCAge : Warning : No class containing age found in the pop");
    numC = -1;
  }
  else 
    numC = (int) ((age - pop->agemin) / pop->largClass);
  return(numC);
}


//retourne la borne inf de la classe numC
float borneInf(int numC,PopMyc *pop) {
  return (pop->agemin + numC * pop->largClass);
}


//fonction retournant le decT d'un age dans une pop

float decTAge(float age,PopMyc *pop) {
  float dec = 0.;
  int numC = numCAge(age,pop);
  if (numC > 0)
    dec = age - borneInf(numC, pop) - pop->TimeBuffer;
  return(dec);
}


void getdecT(PopMyc *pop,float *decT) {
  int i;
  for (i = 0; i < pop->nbClass; i++) 
    decT[i] = pop->decT[posC(i,pop)];
}


// Fonction retournant l'age de la classe numC de la population pop
float ageC(int numC,PopMyc *pop) {
  //show("binf = %f, tb =%f, numC =%d, posC=%d, decT=%f\n",borneInf(numC, pop),pop.TimeBuffer,numC,posC(numC,pop),pop.decT[posC(numC, pop)]);
  return(borneInf(numC, pop) + pop->TimeBuffer + pop->decT[posC(numC, pop)]);
}



// Fonction retournant la borne max d'une pop
float bornemax(PopMyc *pop) {
  return(pop->agemin + pop->nbClass * pop->largClass);
}

//Fonction retournant la classe numC de la population pop
ClassMyc getCmyc(int numC, PopMyc *pop) {
  ClassMyc Cmyc;
  int pos = posC(numC, pop);
  Cmyc.Q = pop->Q[pos];
  Cmyc.N = pop->N[pos];
  Cmyc.age = ageC(numC, pop);
  return(Cmyc);
}

//Fonction retournant la barre numC de la population pop
BarPop getBarPop(int numC, PopMyc *pop) {
  BarPop bp;
  bp.pos = posC(numC, pop);
  bp.Q = pop->Q[bp.pos];
  bp.N = pop->N[bp.pos];
  bp.decT = pop->decT[bp.pos];
  return(bp);
}

//Fonction modifiant la barre bp de la population pop
void PutBarPop(BarPop bp,PopMyc *pop) {
  if (bp.pos < 0 || bp.pos > (pop->nbClass - 1)) 
    show("\n<PopMyc.h> : putBarPop :erreur : essai de modif d'une barre qui n'est pas dans la pop\n");
  else if (pop->decT[bp.pos] != bp.decT)
    show("\n<PopMyc.h> : putBarPop :Warning : essai de modif d'une barre qui n'a pas le meme decT\n");
  pop->Q[bp.pos] = bp.Q;
  pop->N[bp.pos] = bp.N;
  pop->decT[bp.pos] = bp.decT;
}


//Fonction retournant la Quantité totale contenu dans une pop

float Qtot(PopMyc *pop) {
  int i;
  float s = 0;
  for (i = 0; i < pop->nbClass; i++) 
    s += pop->Q[i];
  return(s);
}

//Fonction retournant l'effectif total contenu dans une pop

float Ntot(PopMyc *pop) {
  int i;
  float s = 0;
  for (i = 0; i < pop->nbClass; i++) 
    s += pop->N[i];
  return(s);
}


// idem dans une stack

float Ntot_sc(StackCmyc *sc) {
  int i;
  float s = 0;
  for (i = 0; i < sc->nclass; i++) 
    s += sc->stack[i].N;
  return(s);
}

// surface "disponible" (non protegee par ilim) dans une stack

float Qdisp_sc(StackCmyc *sc) {
  int i;
  float s = 0;
  for (i = sc->ilim + 1; i < sc->nclass; i++) 
    s += sc->stack[i].Q;
  return(s);
}

float Qdisp_sbp(StackBpop *sc) {
  int i;
  float s = 0;
  for (i = sc->ilim + 1; i < sc->nclass; i++) 
    s += sc->stack[i].Q;
  return(s);
}

float NQtot(PopMyc *pop) {
  int i;
  float s = 0;
  for (i = 0; i < pop->nbClass; i++) 
    s += pop->N[i] * pop->Q[i];
  return(s);
}

//Nombre classes non vides (Q > 0) de la pop

int NFilledC(PopMyc *pop) {
  int n = 0;
  int i;
  for (i = 0; i < pop->nbClass; i++) 
    if (pop->Q[i] > 0)
      n++;
  return(n);
}

//rempli le tableau Q avec pop.Q dans l'ordre
void getQ(PopMyc *pop,float *Q) {
  int i;
  for (i = 0; i < pop->nbClass; i++) 
    Q[i] = pop->Q[posC(i,pop)];
}

void getQd(PopMyc *pop,double *Q) {
  int i;
  for (i = 0; i < pop->nbClass; i++) 
    Q[i] = (double) pop->Q[posC(i,pop)];
}

//rempli le tableau N anvec pop.N dans l'ordre
void getN(PopMyc *pop,float *N) {
  int i;
  for (i = 0; i < pop->nbClass; i++) 
    N[i] = pop->N[posC(i,pop)];
}

//retourne l'effectif d'une pop entre deux age inclus)
float NtotEntre(PopMyc *pop,float agemin,float agemax) {
  int i;
  ClassMyc Cmyc;
  float s = 0;
  for (i = 0; i < pop->nbClass; i++) {
    Cmyc = getCmyc(i,pop);
    if (Cmyc.age >= agemin && Cmyc.age <= agemax)
      s += Cmyc.N;
  }
  return(s);
}

//retourne la quantité d'une pop entre deux age inclus)
float QtotEntre(PopMyc *pop,float agemin,float agemax) {
  int i;
  ClassMyc Cmyc;
  float s = 0;
  for (i = 0; i < pop->nbClass; i++) {
    Cmyc = getCmyc(i,pop);
    if (Cmyc.age >= agemin && Cmyc.age <= agemax)
      s += Cmyc.Q;
  }
  return(s);
}

//retourne la quantité d'une pop fois lenombre entre deux age inclus)
float NQtotEntre(PopMyc *pop,float agemin,float agemax) {
  int i;
  ClassMyc Cmyc;
  float s=0;
  for (i = 0; i < pop->nbClass; i++) {
    Cmyc = getCmyc(i,pop);
    if (Cmyc.age >= agemin && Cmyc.age <= agemax)
      s += Cmyc.Q * Cmyc.N;
  }
  return(s);
}


// retourne une pile de priorite des classes de la pop = tableau contenant les classes non vides, de la plus vielle a la plus jeune

StackCmyc PopAsStack(PopMyc *pop) {
  StackCmyc s = newStackCmyc(pop->nbClass);
  int i;
  int n = 0;
  ClassMyc Cmyc;
  for (i = pop->nbClass - 1; i >= 0 ; i--) {
    Cmyc = getCmyc(i,pop);
    if (Cmyc.N > 0) {
      s.stack[n] = Cmyc;
      n++;
    }
  }
  s.nclass = n;
  return(s);
}

StackBpop PopAsStack_bp(PopMyc *pop) {
  StackBpop s = newStackBpop(pop->nbClass);
  int i;
  int n = 0;
  BarPop Bpop;
  for (i = pop->nbClass - 1; i >= 0 ; i--) {
    Bpop = getBarPop(i,pop);
    if (Bpop.N > 0) {
      s.stack[n] = Bpop;
      n++;
    }
  }
  s.nclass = n;
  return(s);
}


//fonction pour l'affichage


void getType(char * type,PopMyc *pop) {
  sprintf(type, pop->type == pTLes ? "Lesions" : (pop->type == pTUdin ? "Udin" : (pop->type == pTSpo ? "Sporulating area" : "unknown" )));
}

char *strClass(ClassMyc *Cmyc) {
  char str[100];
  sprintf(str,"ClassMyc: age=%g, Q=%g, N=%.0f\n",Cmyc->age,Cmyc->Q,Cmyc->N);
  return(strdup(str));
}

void PrintClass(ClassMyc Cmyc) {
  char* str = strClass(&Cmyc);
  show(str);
}

//ajouter nb classes non vides et faire description longue dans Printpop. pour str ne fare que le resume


// Representation string de la pop

char* strPop(PopMyc *pop) {
  char str[200];
  char type[20];
  getType(type,pop);
  sprintf(str,"\nPopMyc object: type: %s, agemin: %.2f, agemax: %.2f, nbClass: %d, largClass: %.2f\nQtot= %.2f, Ntot= %.2f (in %d non empty classes)\nInternals: pC0=%d, TimeBuffer= %f\n",type, pop->agemin,pop->agemax,pop->nbClass,pop->largClass,Qtot(pop),Ntot(pop),NFilledC(pop),pop->pC0,pop->TimeBuffer);
  return(strdup(str));
}
    

void PrintPop(PopMyc pop) {
  char* str = strPop(&pop);
  show(str);
  if (Qtot(&pop) > 0) {
    int i;
    float bmin, bmax;
    show("PopMyc contents:\n");
    ClassMyc Cmyc;
    for (i = 0; i < pop.nbClass; i++) {
      Cmyc = getCmyc(i,&pop);
      bmin = pop.agemin + i * pop.largClass;
      if (i == pop.nbClass - 1)
	bmax = bornemax(&pop);
      else
	bmax = pop.agemin + ( i + 1) * pop.largClass;
      if (Cmyc.Q > 0) {
	if (pop.type == pTUdin)
	  show("Class %d ([%.2f;%.2f[) : %g Udin (age=%.2f)\n",i ,bmin, bmax, Cmyc.N, Cmyc.age);
	else if (pop.type == pTSpo)
	  show("Class %d ([%.2f;%.2f[) : %g cm2 (agedd  = %.2f, %.2f sporulation events)\n",i,bmin, bmax, Cmyc.Q, Cmyc.N, Cmyc.age);
	else if (pop.type == pTHlat)
	  show("Class %d ([%.2f;%.2f[) : %.0f hour since incubation (hour favorable to latency  = %.0f, agedd = %.2f)\n",i,bmin, bmax, Cmyc.Q, Cmyc.N, Cmyc.age);
	else
	  show("Class %d ([%.2f;%.2f[) : %g cm2 (%.2f lesions,age=%.2f)\n",i,bmin, bmax, Cmyc.Q, Cmyc.N, Cmyc.age);
      } else
	show("Class %d ([%.2f;%.2f[) : empty\n", i, bmin, bmax);
    }
  } else
    show("PopMyc currently empty\n\n");
}


//ajout d'une classe C dans la population pop. Renvoie -1 si ajout impossible, et le numero de la classe dans la pop sinon. 

int AddCtoPop(ClassMyc C, PopMyc *pop) {
  int numC = numCAge(C.age, pop);
  if (numC < 0) 
    show("<PopMyc.h> :AddCtoPop : Error Class number not found in the pop. Class not added\n");
  else if (C.Q > 0) {
    int pos = posC(numC, pop);
    (pop->decT)[pos] = (C.Q * decTAge(C.age, pop) + (pop->Q)[pos] * (pop->decT)[pos]) / (C.Q + (pop->Q)[pos]);
    (pop->Q)[pos] += C.Q;
    (pop->N)[pos] += C.N;
  } 
  return(numC);
}

//ajout d'une classe C dans la population pop et maj en parrallele de variables annexe dans la popannexe (moyenne pondere par Q). Renvoie -1 si ajout impossible, et le numero de la classe dans la pop sinon. 

int AddCtoPops(ClassMyc C, float Nbis, float Qbis, PopMyc *pop, PopMyc *popbis) {
  int numC = numCAge(C.age, pop);
  if (numC < 0) 
    show("<PopMyc.h> :AddCtoPop : Error Class number not found in the pop. Class not added\n");
  else if (C.Q > 0) {
    int pos = posC(numC, pop);
    float wClass = C.Q / (C.Q + (pop->Q)[pos]);
    float wPop = (pop->Q)[pos]  / (C.Q + (pop->Q)[pos]);
    (pop->decT)[pos] = wClass * decTAge(C.age, pop) + wPop * (pop->decT)[pos];
    (pop->Q)[pos] += C.Q;
    (pop->N)[pos] += C.N;
    (popbis->decT)[pos] = (pop->decT)[pos];
    (popbis->N)[pos] = wClass * Nbis + wPop * (popbis->N)[pos];
    (popbis->Q)[pos] = wClass * Qbis + wPop * (popbis->Q)[pos];
  }
  return(numC);
}


float propEntiere(float N, float prop) {
  float Nent = floorf(N * prop);
  float frac = N * prop - Nent;
  double r = rand() / (double) RAND_MAX;
  if (r <= frac)
    Nent += 1.;
  //show("propEntiere : N=%f p=%f frac=%f r=%g Nent=%f\n",N,prop,frac,r,Nent);
  return(Nent);
}


//elimine des Q pour reduire la quantité Q totale de la pop,renvoie la quantité enlevé : on enleve des surfaces a chaque effectif en gardant le meme effectif
//ne modifie pas le timebuffer (a faire un jour)
//adapté pour popChlo (avec effectifs nul et surface cumulé dans Q) 

float RemoveQtoPop(float Q, PopMyc *pop) {
  float Qt = Qtot(pop);
  float Qr = MIN(Qt, Q);
  if (Qt != 0) {
    float txred = 1 - Qr / Qt;
    int i;
    for (i = 0; i < pop->nbClass; i++) 
      if ((pop->Q)[i] > 0) {
	(pop->Q)[i] *= txred;
      }
  }
  return(Qr);
}

float RemoveNtoPop(float S, PopMyc *pop) {//on enleve la meme surface a chaque classe en supprimant des effectifs (MAIS changements effectif depend classe puisque leurs tailles sont differentes)
  float Qt = Qtot(pop);
  float Qr = MIN(Qt, S);
  if (Qt != 0) {
    float Sn = Qr / Ntot(pop);//surface a enlever par classe
    int i;
    float n;
    for (i = 0; i < pop->nbClass; i++)
      if ((pop->Q)[i] > 0) {
	n = MIN((pop->N)[i], (Sn * (pop->N)[i] / (pop->Q)[i])); //nombre d'inc a enlever par classe
	(pop->Q)[i] -= n * (pop->Q)[i] / (pop->N)[i];
	(pop->N)[i] -= n;
      }
  }
  return(Qr);
}



//Tire une classe au hasard dans la pop en fonctions des effectifs.renvoie -1 si la pop est vide

int TireCbyNd(PopMyc *pop, int debug) {
  int numC = -1;
  int i; 
  float Nt = Ntot(pop);
  if (Nt < 1)
    show("<PopMyc> : TireCbyN : Not enougth N in the pop (Nt=%f)!\n",Nt);
  else {
    float *Nrcum = (float*) calloc(pop->nbClass, sizeof(float));
    getN(pop, Nrcum);
    float cum = 0;
    for (i = 0; i < pop->nbClass; i++) {
      cum += Nrcum[i] / Nt;
      Nrcum[i] = cum;
    }
    double r = rand() / (double) RAND_MAX;
    while(r == 0)
      r = rand() / (double) RAND_MAX;
    if (debug)
      show("TireCbyN: r=%g\n",r);
    for (i = 0; i < pop->nbClass; i++)
      if (r <= Nrcum[i])
	break;
    numC = i;
    free(Nrcum);
  }
  return numC;
}

int TireCbyN(PopMyc *pop) {
  return TireCbyNd(pop,0);
}
      
//enleve une quantité Q d'une pop en enlevant des effectifs au hasard, au prorata des effectifs de classes. 
//Donne un resultat a une lesion pres (arrondi inferieur ou superieur). Le parametre p est la proba de killer la lesion faisant l'arrondi. 
//Renvoie la suface effectivement killée 

float RemoveQbyNd(float Q, PopMyc *pop,float pKillLast,int debug) {

  if (debug) {
    show("RemoveQbyN : pop\n");
    PrintPop(*pop);
  }

  int numc;
  BarPop bp;
  float sc;
  float Qdone = 0;
  float Qt = Qtot(pop);

  if (Q >= Qt - epsillon) {//on considere Qt = Q et on vide tout
    Qdone = Qt;
    ResetPop(pop);
  } else {
    while (Q > epsillon) {//boucle on enleve une lesion en exces
      if (debug)
	show("RemoveQbyN : Ntot=%f Qtot=%f Q=%f\n",Ntot(pop),Qtot(pop),Q);
      numc=TireCbyNd(pop,debug);
      if (debug)
	show("RemoveQbyN : numc a enlever=%d\n",numc);
      if (numc < 0) {
	show(" <Erreur!> PopMyc.h: RemoveQbyN : impossible d'enever la totalité de la demande !\n");
	exit(0);
      }
      bp = getBarPop(numc,pop);
      if (debug)
	show("RemoveQbyN : bp.Q=%f bp.N=%f\n",bp.Q,bp.N);
      sc = bp.Q / bp.N;
      bp.Q -= sc;
      bp.N -= 1;
      Q -= sc;
      Qdone += sc;
      PutBarPop(bp,pop);
    }
  }
  
  if (Q < 0) {//il y a un arrondi
    double r = rand() / (double) RAND_MAX;
    if (debug)
      show("r=%f\n",r);
    if (r >= pKillLast) {
      bp = getBarPop(numc,pop);
      bp.Q += sc;
      bp.N += 1;
      Qdone -= sc;
      PutBarPop(bp,pop);
    }
    else
      if (debug)
	show("RemoveQbyN : elimination de la derniere les au hasard\n");
  }    
  return(Qdone);
}

float RemoveQbyN(float Q, PopMyc *pop,float pKillLast) {
  return RemoveQbyNd(Q, pop, pKillLast, 0);
}

// elimines des lesions dans une stack, de la moins prioritaire a la dernier dispo (ie avec i > ilim).

float killInStack(StackCmyc *st, float S2kill) {
  float Skilled = 0;
  float Nkill;
  float Sles;//surface d'une lesion
  ClassMyc Cmyc;
  int i;
  int nc = st->nclass; 

  for (i = st->nclass - 1 ; i > st->ilim ; i--) 
    if (S2kill > 0) {
      Cmyc = st->stack[i];
      if (Cmyc.N * Cmyc.Q > 0) {
	if (S2kill >= Cmyc.Q) {
	  nc--;// equivaut a supression
	  Skilled += Cmyc.Q;
	  S2kill = MAX(0,S2kill - Skilled);
	} else {
	  Nkill = propEntiere(Cmyc.N, S2kill / Cmyc.Q);
	  if (Nkill >= Cmyc.N) {
	    nc--;// supression
	    Skilled += Cmyc.Q;
	    S2kill = MAX(0,S2kill - Skilled);
	  } else {
	    st->stack[i].N -= Nkill;
	    st->stack[i].Q -= Nkill * Sles;
	    Skilled += Nkill * Sles;
	    S2kill = MAX(0,S2kill - Skilled);
	  }
	}
      }
    }
  st->nclass = nc;
  return(Skilled);
}

//version stackbp 

float killInStackBp(StackBpop *st, float S2kill) {
  int i;
  StackCmyc stc = newStackCmyc(st->nclass);
  stc.ilim = st->ilim;
  stc.nclass = st->nclass;
  for (i = 0; i < st->nclass ; i++) {
    stc.stack[i].Q = st->stack[i].Q;
    stc.stack[i].N = st->stack[i].N;
  }
  float Skilled = killInStack(&stc, S2kill);
  st->ilim = stc.ilim;
  st->nclass = stc.nclass;
  for (i = 0; i < stc.nclass ; i++) {
    st->stack[i].Q = stc.stack[i].Q;
    st->stack[i].N = stc.stack[i].N;
  }
  freeStackCmyc(&stc);
  return(Skilled);
}
 
ClassMyc ClassEqQ(PopMyc *pop) {//renvoie une classeMyc ayant l'age moyen des classes pondérée par les surfaces, l'effectif total et la surface totale.
  ClassMyc Ceq, Cmyc;
  Ceq.Q = 0;
  Ceq.age = 0;
  Ceq.N = 0;
  int i;
  for (i = 0; i < pop->nbClass; i++) {
    Cmyc = getCmyc(i,pop);
    Ceq.Q += Cmyc.Q;
    Ceq.N += Cmyc.N;
    Ceq.age += Cmyc.age * Cmyc.Q;
  }
  if (Ceq.Q != 0) 
    Ceq.age /= Ceq.Q;
  return(Ceq);
}
    

// fonction faisant vieillir une pop d'un temp dt, met a jour la pop et renvoie une pop sortante. age pop sortante = age compté depuis la sortie
PopMyc AgePop(PopMyc *pop,float dt) {

  ClassMyc Cmyc;
  //pop sortante
  PopMyc popout;
  if ((pop->agemax - pop->agemin) < pop->largClass)
    popout = newPopMyc(0, pop->largClass, pop->largClass, pop->type);
  else
    popout = newPopMyc(0, pop->agemax - pop->agemin, pop->largClass, pop->type);
  
  int shift,i;//decalage d'indice pour mettre a jour la pop

  //___etape 1 shift et maj TimeBuffer
  shift = (int) ((dt + pop->TimeBuffer) / pop->largClass);
  pop->TimeBuffer = dt + pop->TimeBuffer - shift * pop->largClass;
  popout.TimeBuffer = pop->TimeBuffer;
  //show("dt=%f, shift=%d, TimeBuffer=%f\n", dt,shift,popout.TimeBuffer);


  //___etape 2 transferts pop ->  popout

  //Cas 2.1 : tout sort, on remet tous les effectifs a 0, et on vide dans popout qui a un agemin modifie et pC0 transmis
  float minage;
  if (pop->nbClass > 1)
    minage = MIN(ageC(0, pop), ageC(1, pop));//si dect <0
  else
    minage = ageC(0, pop);
  if (minage + shift * pop->largClass >= pop->agemax) {
    for (i = 0; i < pop->nbClass; i++) {
      popout.Q[i] = pop->Q[i];      
      pop->Q[i] = 0;
      popout.N[i] = pop->N[i];      
      pop->N[i] = 0;
      popout.decT[i] = pop->decT[i];      
      pop->N[i] = 0;
    }
    popout.pC0 = pop->pC0;
    popout.agemin = (shift - pop->nbClass) * pop->largClass + bornemax(pop) - pop->agemax;//on corrige de la duree plus courte de passage de derniere classe
    popout.agemax += popout.agemin;
  } 
  else {
    //Cas 2.2 une partie seulement peut sortir
    if (shift > 0) {
      //2.2.1 deplacemnt des classes. Classes potentiellement sortantes = les shift dernieres classes
      for (i = 0; i < shift; i++) {
	Cmyc = getCmyc(pop->nbClass - 1 - i, pop);
	if (Cmyc.N > 0 || Cmyc.Q > 0) {//sinon, tout est deja a jour
	  Cmyc.age += shift * pop->largClass;//pop n'a pas encore vieilli
	  if (Cmyc.age >= pop->agemax) {
	    Cmyc.age -= pop->agemax;
	    AddCtoPop(Cmyc, &popout);
	    ResetClass(pop->nbClass - 1 - i, pop);
	  }
	}
      }
    }

    //_____etape 2.2.2: deplacement des classes de pop : changement d'indice pc0
    pop->pC0 -= shift;
    pop->pC0 = (pop->pC0 < 0) ? (pop->pC0 + (shift / pop->nbClass + 1) * pop->nbClass) : pop->pC0;
    //___etape 2.2.3 : controle des sorties des deux nouvelles dernieres classes de pop
    for (i = 0; i < MIN(2,pop->nbClass); i++) {
	Cmyc = getCmyc(pop->nbClass - 1 - i, pop);
	if (Cmyc.N > 0 || Cmyc.Q > 0) {//sinon, tout est deja a jour
	  if (Cmyc.age >= pop->agemax) {
	    Cmyc.age -= pop->agemax;
	    AddCtoPop(Cmyc, &popout);
	    ResetClass(pop->nbClass - 1 - i, pop);
	  }
	}
    }
  }//fin etape2

  //si il n'y a plus rien dans  pop, le TimeBuffer est reinitialisé
  if (Ntot(pop) <= 0 && Qtot(pop) <= 0) pop->TimeBuffer = 0;
  return(popout);
}

//fonction simulant AgePop et renvoyant popout, mais en laissant intact la pop recu en argument

/* PopMyc SimAgePop(PopMyc *pop,float dt) { */
/*   PopMyc popout; */
/*   float *pQin,*pNin; */
/*   int i;  */
/*   //changement d'adresse du pointeur Q et recopie du tableau */
/*   pQin = pop->Q; */
/*   pNin = pop->N; */
/*   pop.Q = (float*) calloc(sizeof(float), pop.nbClass); */
/*   pop.N = (float*) calloc(sizeof(float), pop.nbClass); */
/*   for (i = 0; i < pop.nbClass; i++) { */
/*     pop.Q[i] = pQin[i]; */
/*     pop.N[i] = pNin[i]; */
/*   } */
/*   popout = AgePop(&pop,dt); */
/*   free(pop.Q); */
/*   free(pop.N); */
/*   return(popout); */
/* } */

//croisssance d'une classe de Q a Q*multmax, limité par une taille max de lésion et un taux de satisfaction txsat de la surface a creer
//renvoie le delta
float expandClass(ClassMyc *C,float multmax,float Slmax,float txsat) {
  float Qini = C->Q;
  float Qmax = MIN(C->Q * multmax, C->N * Slmax);
  C->Q = Qini + (Qmax - Qini) * txsat;
  return(C->Q - Qini);
}

//fonction calculant la demande de croisssance d'une pop de Sini à Sini*mult.Limitation par taille max de lesion

float calcdemC(PopMyc *pop,float mult,float Slmax) {
  float dem = 0;
  int i;
  ClassMyc Cmyc;
  for (i = 0; i < pop->nbClass; i++) {
    Cmyc = getCmyc(i,pop);
    dem += expandClass(&Cmyc, mult, Slmax, 1);
  }
  return(dem);
}

//fonction faisant croitre en surface une pop de S a S*multmax avec taille limite de lésion et tx de stisfaction des demandes de croissance
//si la croissance est limité alors certaines (ou toutes) les classes ont une croisssance inferieure a multmax
//renvoie le total de surface nouvellement crée

float expandPop(PopMyc *pop,float multmax,float Slmax,float txsat) {
  float newS = 0;
  ClassMyc Cmyc;
  int i;
  for (i = 0; i < pop->nbClass; i++) {
    Cmyc = getCmyc(i, pop);
    newS += expandClass(&Cmyc, multmax, Slmax, txsat);
    pop->Q[posC(i,pop)] = Cmyc.Q;
  }
  return(newS);
}

//Renvoie le Svert a enlever d'une demande a enlever au hasard dans le vert et dans une pop. On decoupe la demande en paquet de taille la surface moyeen des surface dans la pop, et on affecte les paquets au vert et a la pop (pour eviter d'appliquer RemoveQbyN si Svert >>>> Spop

/* float randKillVert(float Dem,float Svert,PopMyc pop) { */
/*   float Sk=0;//surface a killer */
/*   int i;  */
/*   float Nt=Ntot(pop); */
/*   if (Nt<=0) { */
/*     show("<PopMyc> : randKillVert : plus d'effectif !\n"); */
/*     Sk=Dem; */
/*   } */
/*   else { */
/*     float Spop=Qtot(pop); */
/*     float Selem=Spop/Nt; */
/*     int nbpak=(int)floor(Dem/Selem); */
/*     float plim=Svert/(Svert+Spop); */
/*     double r; */
/*     for (i=0;i<nbpak;i++) { */
/*       r=rand()/(double)RAND_MAX; */
/*       if (r<plim) */
/* 	Sk+=Selem; */
/*     } */
/*     //dernier paquet */
/*     if (Dem-nbpak*Selem>0) { *//*       r=rand()/(double)RAND_MAX; */
/*       if (r<plim) */
/* 	Sk+=Dem-nbpak*Selem; */
/*     } */
/*     return(Sk); */
/*   } */
/* } */
   
#endif
