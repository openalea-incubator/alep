// bibliotheque de fonction pour gerer des secteurs de feuille malade
//"modèle maladie"
//C. Fournier 2006

#ifndef _SECTEUR_H
#define _SECTEUR_H

#include <stdlib.h>

#include "Lesions.h"
#include "math_utils.h"

//etat des surface et des effectifs de pop 
//nombre d'etat
#define nbEtatS 6
const char* EtatS[nbEtatS]={"SVert","SInc","SChlo","SSpo","SEmpty","SSen"};
enum iEtatS {iVert,iInc,iChlo,iSpo,iEmpty,iSen};
int ColEtatS[nbEtatS]={2,7,3,4,5,6};//couleurs

#define nbEtatP 3
const char* namesEtatP[nbEtatP]={"NInc","NChloNec","NNec"};
enum iEtatP {iNInc,iNChloNec,iNNec};

//variables de retour de la dispersion


//enregistrements
// SvertUdin: Surface verte en place au moment de l'interception des Udin du pas de temps(ie Svert avant modification eventuelle par la maladie si cette phase a lieu apres dans le pas de temps)
// UdFeu,UdinFeu,UdinVert,UdinTotVert: nombres d'UD (avant lavage),d'UDIN tombées sur les secteurs (tout type de surface confondu),d'Udin tombées sur le vert (avant germinations) 

enum SummaryType {somme, moypond, mini, maxi};
#define nSectRec 9
enum SectRec {shmin, shmax, sNbEvSpl, sEcla, sEclin, sSvertUdin, sUdFeu, sUdinFeu, sUdinVert};
const char *SectRecNames[nSectRec]={"hmin", "hmax", "NbEvSpl", "Ecla", "Eclin", "SvertUdin", "UdFeu", "UdinFeu", "UdinVert"};
int SectRecST[nSectRec]={mini, maxi, moypond, somme, somme, somme, somme, somme, somme};

    
typedef struct SECTEUR {
  
  //Etat des surfaces visibles
  float Svert;//surface verte
  float Sdsen;//Surface morte par senescence naturelle
  Lesions Les;//lesions
  float SectRec[nSectRec];//enregistrements
} Secteur;

//raz des enregistrements
void resetSectRec(Secteur *s) {
  int i;
  for (i = 0; i < nSectRec; i++)
    s->SectRec[i] = 0;
}

// enregistrement d'une variable
void setSectRec(Secteur *s,int key,float value) {
  s->SectRec[key] = value;
}

void addSectRec(Secteur *s, int key, float value) {
  s->SectRec[key] += value;
}

void showSectRec(Secteur s) {
  show("Enregistrement du secteur :");
  int i;
  for (i = 0; i < nSectRec; i++)
    show(" %s=%g", SectRecNames[i], s.SectRec[i]);
  show("\n");
}
//fonctions pour l'initialisation

Secteur newSect(ParCycle *par) {
  Secteur sect;
  sect.Les = newLes(par);
  sect.Svert = 0;
  sect.Sdsen = 0;
  resetSectRec(&sect);
  return(sect);
}

void resetSect(Secteur *s) {
  s->Svert = 0;
  s->Sdsen = 0;
  resetSectRec(s);
  ResetLes(&(s->Les));
}


void freeSect(Secteur *s) {
  freeLes(&(s->Les));
}

//fonction pour la gestion d'un secteur

//maj l'exposition (visibilité) d'un secteur lors de l'extension
void exposeSect(Secteur *s,float deltaS) {
    s->Svert += deltaS;
}


 //surface totale visible du secteur
float Stot(Secteur s) {
  return(s.Svert + s.Sdsen + SLesTot(&(s.Les)));
}

// developpement des lesions

void devLesSect_d(Secteur *s, int dth,float *T,float *PPFD, float *Rh, int debug) {
  devLes_d(&(s->Les),&(s->Svert),dth,T,PPFD,Rh,debug);
}

void devLesSect(Secteur *s, int dth,float *T,float *PPFD, float *Rh) {
  devLes_d(&(s->Les),&(s->Svert),dth,T,PPFD,Rh,0);
}

//senescence d'un secteur. Reduit la demande de la surace efectivement senescee
//demSen : demande de surface a faire senescer (senescence nouvelle)

void senSect(Secteur *s,float *demSen) { 

  if (*demSen > 0) {
    float Svert = s->Svert;
    float SLesVert = SLesOnGreen(&(s->Les));
    float Srectot = MIN(*demSen, Svert + SLesVert);//surface totale concernee par le passage de la senescence non influencée par la presence de lesions
    if (Srectot > 0) {
      // recouvrement lesions
      float Storec = Srectot * SLesVert / (Svert + SLesVert); 
      senLes(Storec,&(s->Les),&(s->Svert),&(s->Sdsen));
      *demSen -= Storec;
      // recouvrement du vert
      Storec = Srectot * Svert / (Svert + SLesVert);
      s->Svert -= Storec;
      s->Sdsen += Storec;
      *demSen -= Storec;
    }
  }
}

//production eclin lors des pluies

void RainOnSect(Secteur *s, double Eclas, double *Eclins, double *QSpores) {
  *Eclins = 0;
  *QSpores = 0;
  float Sspo = SLesSpo(&(s->Les));
  if (Sspo > 0) 
    RainOnLes(&(s->Les),Eclas * Sspo / Stot(*s), Eclins, QSpores);
  // enregistrements secteur
  addSectRec(s, sNbEvSpl, 1);
  addSectRec(s, sEcla, Eclas);
  addSectRec(s, sEclin, *Eclins);
}

//depot udin

void DepositOnSect(Secteur *s, double Ud, double Udin, float tpluie, double *Udvert, double *Udinvert) {
  float fvert = 0;
  if (Stot(*s) > 0) 
    fvert = s->Svert / Stot(*s);
  *Udvert = (double) roundP(Ud * fvert);
  *Udinvert = (double) roundP(Udin * fvert);
  addSectRec(s,sUdFeu,*Udvert); 
  addSectRec(s,sUdinFeu,*Udinvert); 
  addUdin(&(s->Les),*Udinvert,*Udinvert,tpluie);
}

//infection d'un secteur par une quantité de surface sporulante Q
void infectSect(Secteur *s, float Q) {

  float S2kill = Q;
  float Skilled = 0;
  float Sk;
  
  //diminution des surfaces affectées
  if (s->Svert >= S2kill) {
    s->Svert -= S2kill;
    Skilled = S2kill;
  } else {
    show("<Secteur.h> infectSect : <! Warning !> : surface d'infection forcée (Q=%g)> S vert disponible (Svert=%g), trying to adjust with Inc...\n", Q, s->Svert);
    S2kill -= s->Svert;
    s->Svert = 0;
     if (SLesInc(&(s->Les)) >= S2kill)
       Skilled += RemoveQtoPop(S2kill, &(s->Les.PopInc)); // revient a une 'decroissance des surface des lesions inc
     else {
       show("<Secteur.h> infectSect : <! Warning !> : surface d'infection forcée > S vert + Sinc disponible, trying to adjust with Chlo...\n");
       S2kill -= SLesInc(&(s->Les));
       ResetPop(&(s->Les.PopInc));
       if (SLesChlo(&(s->Les)) >= S2kill)
	 Skilled += RemoveQtoPop(S2kill, &(s->Les.PopChlo));
       else {
	 show("<Secteur.h> infectSect : <! Warning !> : surface d'infection forcée (Q=%g)> S vert + Sinc + SChlo disponible (Svert+Sinc+SChlo=%g), scenario d'infection modifie...\n", Q, Q - S2kill);
	 Skilled += SLesChlo(&(s->Les));
	 float S;
	 videChlo(&(s->Les), &S);
       }
     }
  }
  
  //creation des lesions: 
  addQLes(&(s->Les), Skilled);
  if (VERBOSE>0) 
    PrintPop(s->Les.PopSpo);
}

void getEtatS(Secteur s,float eS[]) {
  //rempli eS avec les surface du secteur
  eS[iVert] = s.Svert;
  eS[iInc] = SLesInc(&(s.Les));
  eS[iChlo] = SLesChlo(&(s.Les));
  eS[iSpo] = SLesSpo(&(s.Les));
  eS[iEmpty] = s.Les.SEmpLes;
  eS[iSen] = s.Sdsen;
}

void getEtatP(Secteur s,float eP[]) {
  //rempli eS avec les surface du secteur
  eP[iNInc] = NLesInc(&(s.Les));
  eP[iNChloNec] = NLesChloNec(&(s.Les));
  eP[iNNec] = NLesNec(&(s.Les));
}


void calcPropSurf(Secteur s,float *p) {
  //rempli un vecteur des proportions relatives occupées par les surfaces du secteur
  float St = Stot(s);
  getEtatS(s, p);
  int i;
  if (St > 0) 
    for (i = 0; i < nbEtatS; i++)
      p[i] /= St;
  else
    for (i = 0; i < nbEtatS; i++)
      p[i] = 0;
}

//  reset enregistrement lors de la  dispersion

void newSectRecDisp(Secteur *s, float hbase, float htop) {

  // reset sect rec et rec caracteristiques structurales des secteurs
    resetSectRec(s);
    setSectRec(s,shmin,hbase);
    setSectRec(s,shmax,htop);
    setSectRec(s,sSvertUdin, s->Svert);
}

      
//fonction pour l'affichage

void PrintSect(Secteur s) {
  int i;
  float p[nbEtatS];
  float stot = Stot(s);
  calcPropSurf(s, p);
  show("EtatS (St=%.2f) -> ",stot);
  for (i = 0; i < nbEtatS; i++)
    show("%s : %.3g ", EtatS[i], p[i] * stot);
  show(" (dont %.2g cm2 de lesions recouvertes)\n", s.Les.SLesInSen);
}

//sortie EtatS d'un secteur

void writeEtatS(FILE *pf,Secteur *s) {

  int i;

  float eS[nbEtatS];
  float eP[nbEtatP];
  float *rS;
  float *lS;

  getEtatS(*s,eS);
  getEtatP(*s,eP);
  rS = s->SectRec;
  lS = s->Les.LesRec;

  fprintf(pf,"%g\t%g",Stot(*s),NLesTot(&(s->Les)));
  for (i=0;i<nbEtatS;i++) 
    fprintf(pf,"\t%g",eS[i]);
  for (i=0;i<nbEtatP;i++) 
    fprintf(pf,"\t%g",eP[i]);
  for (i=0;i<nSectRec;i++) 
    fprintf(pf,"\t%g",rS[i]);
  for (i=0;i<s->Les.nLesRec;i++) 
    fprintf(pf,"\t%g",lS[i]);
  fprintf(pf,"\n");
}

void showEtatSect(Secteur *s) {

int i;

  float eS[nbEtatS];
  float eP[nbEtatP];
  float *rS;
  float *lS;

  getEtatS(*s,eS);
  getEtatP(*s,eP);
  rS = s->SectRec;

  show("Stot: %g Ntot: %g\n",Stot(*s),NLesTot(&(s->Les)));
  for (i=0;i<nbEtatS;i++) 
    show("%s: %g ",EtatS[i],eS[i]);
 show("\n");
  for (i=0;i<nbEtatP;i++) 
    show("%s: %g ",namesEtatP[i],eP[i]);
show("\n");
  for (i=0;i<nSectRec;i++) 
    show("%s: %g ",SectRecNames[i],rS[i]);
  show("\n");

}

#endif 
