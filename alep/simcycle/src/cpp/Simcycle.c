// ajout code pour wrapping avec ctypes


//code brut septo3D

#include "PopMyc.h"
#include "Lesions.h"
#include "Meteo.h"
#include "Secteur.h"

// Passage a restypes simple (void) des fonctions retournant des structures.


// ______________________________PopMyc.h

void newCmyc_void(ClassMyc *Cmyc, float age, float Q, float N) {
  *Cmyc = newCmyc_p(age,Q,N);
}

void newPopMyc_void(PopMyc *pop,float agemin,float agemax,float largClass,int popType) {
  *pop = newPopMyc(agemin,agemax,largClass,popType);
}

void newStackCmyc_void(StackCmyc *s,int size) {
  *s = newStackCmyc(size);
}

void newStackBpop_void(StackBpop *s,int size) {
  *s = newStackBpop(size);
}

void getCmyc_void(ClassMyc *Cmyc,int numC, PopMyc *pop) {
  *Cmyc = getCmyc(numC,pop); 
}

void getBarPop_void(BarPop *bp,int numC, PopMyc *pop) {
  *bp = getBarPop(numC,pop);
}

void PopAsStack_void(StackCmyc *s, PopMyc *pop) {
  *s = PopAsStack(pop);
}

void PopAsStack_bp_void(StackBpop *s, PopMyc *pop) {
  *s = PopAsStack_bp(pop);
}

void ClassEqQ_void(ClassMyc *Cmyc,PopMyc *pop) {
  *Cmyc = ClassEqQ(pop);
}

void AgePop_void(PopMyc *popout, PopMyc *pop, float dt) {
  PopMyc p = AgePop(pop,dt);
  popout = &p;
}

// ______________________________Lesions.h

void getParCycle_void(ParCycle *p) {
  *p = getParCycle();
}

void newLes_void(Lesions *les, ParCycle *par) {
  *les = newLes(par);
}

void getPopUdin_void(PopMyc *p, Lesions *les, float h) {
  *p = getPopUdin(les,h);
}

void calcGerm_void(ClassMyc *Cmyc, Lesions *les,int dth, float *PPFD, float *Rh, int debug) {
  *Cmyc = calcGerm(les,dth,PPFD,Rh,debug);
}

void getNewChlo_void(PopMyc *p, Lesions *les, int hdeb,int hfin, float *T, float *PPFD, float *Rh, float eradicant_efficacy) {
  *p = getNewChlo(les,hdeb,hfin,T,PPFD,Rh,eradicant_efficacy);
}

int getnLesRec() {
  return(NLESREC);
}
