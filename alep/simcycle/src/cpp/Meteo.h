// Fonctions pour la lecture /manipulation des données meteo


// Pour les pluies, on se ramène a des evenements ponctuels : la pluie sera traite a un pas de temps (celui qui contient son heure moyenne ponderee par l'intensite). Cependant, les udin deposees durant un pas de temps formeront une seule cohorte, d'age egal au temps moyen des pluies(pondere par les intensite) du pas de temps

// Pour les periodes humides, on retient dans un pas de temps toutes les periodes humides permettant potentiellement au moins une germination dans ce pas de temps (ie n'etant pas terminees au debut du pas de temps et totalisant le nombre d'heures favorables requis AVANT la fin du pas de temps) 

#ifndef _METEO_H
#define _METEO_H

#include <stdio.h>
#include <stdlib.h> //exit,calloc
#include <math.h> //floor
#include <string.h> //strcpy,strcat

#include "symbols.h"
#include "RefParameter.h"

#define MAX_SIZE_METEODB (366 * 24)
#define MAXPLUIES 365
#define IMINPLUVIO 0.2


// Acces variables Parameter.h

// path du fichier meteo

char* getMeteoPath() {
  char path[100];
  strcpy(path,METEODIR);
  strcat(path,"/");
  strcat(path,METEOFILE);
  return(strdup(path));
}

// idem si appel depuis repertoire test

char* getTestMeteoPath() {
  char path[100];
  strcpy(path,".");
  strcat(path,METEODIR);
  strcat(path,"/");
  strcat(path,METEOFILE);
  return(strdup(path));
}

// structure gestion variables meteo

typedef struct METEODATA {

  // nombre d'entree dans les  tableaux
  int size;
  // annee
  int *annee;
  //numero de jour
  int *numJ;
  // numero de jour since 1er octobre
  int *Jsim;
  // heure hhmm
  int *hhmm;
  //temperature air
  float *T;
  // Relative humidity
  float *Rh;
  // Precipitations (mm)
  float *P;
  // PAR photon flux density (micromol/m2/s)
  float *PPFD;
  //windspeed (m.s-1)
  float *windspeed;
  // Rain event mean intensity (mm/h) : a filtered copy of Precipitations, conserving one non null value per rain event and positioned at the median time of the rain (mean hour of rain  weighted by the precipitation ammount).
  float *IRain;
  // Rain event duration (hours)
  int *DRain;

} MeteoData;


// constructeur

MeteoData newMeteoData(int size) {
  MeteoData m;
  m.size = size;
  m.annee = (int*) calloc(size,sizeof(int));
  m.numJ = (int*) calloc(size,sizeof(int));
  m.Jsim = (int*) calloc(size,sizeof(int));
  m.hhmm = (int*) calloc(size,sizeof(int));
  m.T = (float*) calloc(size,sizeof(float));
  m.Rh = (float*) calloc(size,sizeof(float));
  m.P = (float*) calloc(size,sizeof(float));
  m.PPFD = (float*) calloc(size,sizeof(float));
  m.windspeed = (float*) calloc(size,sizeof(float));
  m.IRain = (float*) calloc(size,sizeof(float));
  m.DRain = (int*) calloc(size,sizeof(int));
  if (m.annee != NULL && m.numJ != NULL && m.Jsim != NULL && m.hhmm != NULL && m.T != NULL && m.Rh != NULL && m.P != NULL && m.PPFD != NULL && m.windspeed != NULL && m.IRain !=NULL && m.DRain !=NULL)
    return(m);
  else {
    show("<Meteo.h> newMeteoData : can't allocate memory for database, exiting..\n");
    exit(0);
  }
}

// destructeur

void delMeteoData(MeteoData *meteo) {
  free(meteo->annee);
  free(meteo->numJ);
  free(meteo->Jsim);
  free(meteo->hhmm);
  free(meteo->T);
  free(meteo->Rh);
  free(meteo->P);
  free(meteo->PPFD);
  free(meteo->windspeed);
  free(meteo->IRain);
  free(meteo->DRain);
}


// construction meteodata partir d'un ou plusieur fichiers


// constructeur colone IRain, DRain discretisation des pluies sur une valeur horaire. mmMin : valeur mini (1 basule) du pluvio. 

void discretiseRain(MeteoData *m, float mmMin) {
  
  float mm;
  int D,ipluie;
  float t;

  int i = 0;
  while (i < m->size) {
    if (m->P[i] <= mmMin)
      i++;
    else {//une nouvelle pluie demarre
      mm = 0;
      D = 0;
      t = 0;
      while (m->P[i] > mmMin && i < m->size) {
	mm += m->P[i];
	D += 1;
	t += i * m->P[i];//ponderation indice par intensite
	i++;
      }
      ipluie = (int) round(t/mm);
      m->IRain[ipluie] = mm / D;
      m->DRain[ipluie] = D;
    }
  }
}

// constructeur a partir d'une base meteo 

MeteoData createMeteo(char *path) {
  MeteoData m = newMeteoData(MAX_SIZE_METEODB);
  FILE *pfmeteo=fopen(path,"r");

  if (pfmeteo==NULL) {
    show("<meteo.h> : createMeteo(): erreur, fichier meteo (%s) introuvable !!",path);
    exit(0);
  }
  
  char header[100];
  fgets(header,100,pfmeteo);
  int i=0;
  float PPFDmolh;
  while (fscanf(pfmeteo,"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n",&(m.annee[i]),&(m.numJ[i]),&(m.Jsim[i]),&(m.hhmm[i]),&PPFDmolh,&(m.T[i]),&(m.Rh[i]),&(m.windspeed[i]),&(m.P[i])) != EOF) {
    m.PPFD[i] = MAX(0,PPFDmolh * 3600. / 1e6); 
    i++;
  }
  fclose(pfmeteo);
  m.size = i;
  // filling IRain
  discretiseRain(&m, IMINPLUVIO);
  return(m);
}


// extraction meteo depuis meteodb pour une periode donnee.Meteo doit pointer vers un meteoData correctement initialise
void getMeteo(MeteoData *meteodb, MeteoData *meteo,int jdeb, int hdeb) {

  while (hdeb > 24) {
    jdeb ++;
    hdeb -= 24;
  }
      
  if (hdeb == 0) {
    jdeb--;
    hdeb = 24;
  }

  int ideb = (jdeb - meteodb->Jsim[0]) * 24 + (hdeb * 100 - meteodb->hhmm[0]) / 100;
  int D = meteo->size;
  
  if (ideb < 0 || ideb + D > meteodb->size) {
    show("<Meteo.h>  GetMeteo : can't find period in meteodb (Jdeb = %d, hdeb = %d, Duration = %d\n",jdeb,hdeb,D);
    exit(0);
  }

  int i;
  for (i = 0; i < D; i++) {
    meteo->annee[i] = meteodb->annee[ideb + i];
    meteo->numJ[i] = meteodb->numJ[ideb + i];
    meteo->Jsim[i] = meteodb->Jsim[ideb + i];
    meteo->hhmm[i] = meteodb->hhmm[ideb + i];
    meteo->T[i] = meteodb->T[ideb + i];
    meteo->Rh[i] = meteodb->Rh[ideb + i];
    meteo->P[i] = meteodb->P[ideb + i];
    meteo->PPFD[i] = meteodb->PPFD[ideb + i];
    meteo->windspeed[i] = meteodb->windspeed[ideb + i];
    meteo->IRain[i] = meteodb->IRain[ideb + i];
    meteo->DRain[i] = meteodb->DRain[ideb + i];
  }
}

// Renvoie le nombre et les indices du tableau meteo pour lesquels la pluie est non nulle 
void findRainIndex(MeteoData *meteo, int *npluies, int *ipluie) {

  int np = 0;
  int i;

  for (i = 0; i < meteo->size; i++) 
    if (meteo->IRain[i] > 0) {
      ipluie[np] = i;
      np++;
    }
  *npluies = np;
}
  
//force la meteo

#define DEFINE_SET(METEOVAR)					\
  void set##METEOVAR(MeteoData *m, int ideb, int ifin, float v) {\
    int i;\
    ifin = MIN(ifin,m->size);\
    for (i = ideb; i < ifin; i++)\
      m->METEOVAR[i] = v;\
  }

DEFINE_SET(T);
DEFINE_SET(Rh);
DEFINE_SET(PPFD);
DEFINE_SET(P);

void setPluie(MeteoData *m, int ipluie, float I, int D) {
  if (ipluie > m->size) {
    show("<Meteo.h> : setPluie: Can't set pluie at requested index : index outside meteo data size !!!\n");
    return;
  }
  int i;
  
  for (i = 0; i < m->size ; i++) {
    if (i == ipluie) {
      m->IRain[i] = 0;
      m->DRain[i] = 0;
    } else {
      m->IRain[i] = I;
      m->DRain[i] = D;
    }
  }
      
}


// Date median d'un tableau meteo

float Dmed(MeteoData *m) {
  return(m->Jsim[0] + (m->size - 1) * 1. / 24 / 2);
}
// 



// string representation

char* strMeteo(MeteoData *m, int maxlg) {

  int np;
  int ip[MAXPLUIES];
  findRainIndex(m,&np,ip);

  maxlg = MIN(maxlg,m->size);
  char str[100 + maxlg * 100];
  int i;
  
  sprintf(str,"MeteoData : %d entries (first %d shown), median Date: %.1f\n", m->size, maxlg,Dmed(m));
  char msg[100];
  sprintf(msg,"Rain events : %d",np);
  strcat(str,msg);
  if (np > 0) {
    strcat(str," (fist five indices :");
    for (i = 0; i < MIN(np,5) ; i++) {
      sprintf(msg," %d",ip[i]);
      strcat(str,msg);
    }
    strcat(str,"...)");
  }
  strcat(str,"\n");
  strcat(str,"annee\tnumJ\tJsim\thhmm\tPPFD\tT\tRh\twind\tRain\tIRain\tDRain\n");
  char lg[100];
  for (i = 0; i < maxlg; i++) {
    sprintf(lg,"\n%d\t%d\t%d\t%d\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%d",m->annee[i],m->numJ[i],m->Jsim[i],m->hhmm[i],m->PPFD[i],m->T[i],m->Rh[i],m->windspeed[i],m->P[i],m->IRain[i],m->DRain[i]);
    strcat(str,lg);
  }
  strcat(str,"\n");
  return(strdup(str));
}

void showMeteo(MeteoData m) {
  char *msg = strMeteo(&m,10);
  show(msg);
}

#endif
