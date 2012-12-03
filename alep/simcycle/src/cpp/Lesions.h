//          Modèle de dévelopment des lesions sur un secteur



// Idee compet entre goute pour occuper l'espace : introduire taille de goute et du coup 'fusion de goutes"

//to do : avoir une pop (et non une classe) de newudin (age = heure du depot dans le pas de temps), pour pouvoir gerer mieux les germinations apres les pluies

//curseur biotrophe/necrotrophe Calcul du tx de satisfaction (actuellement concerne le vert, mais pourrait aussi concerner l'espace dispo sur le mort quand les lesions pouseront en necrotrophe)

//Hypothèses pour le développement des lesions: 

// -1 udin produit 1 lesion la ou elle est tombée
// -croissance des lesions de epsilon à Slmin lors de l'incubation
// -croisssance des lesions de Slmin a Slmax lineaire des chloroses et des necroses par translation
// En cas  d'insuffisance de Surface verte, on empeche les sorties d'incubation et on reduit la croisssance des lesions
// Pour les lesions en incubation, la croissance en surface est theorique : elle permet de leur affecter une surface et de gerer la competition. Dans la realite, le champignon s'etendrait tres rapidement dans la feuille (surface inconnue), puis la densite de mycellium se densifierait localement et finirait par produire une chlorose de taille mesurable. 


// le modele considere des pop de cohortes  de lesions classees par stades (1 cohorte = 1 ensemble de lesions du meme age et du meme secteur) 


#ifndef _LESIONS_H
#define _LESIONS_H

#include <stdio.h>
#include <stdlib.h> //rand,calloc
#include <math.h>

#include "symbols.h"
#include "PopMyc.h"

#include "RefParameter.h" //Param cycle + MAXSUBSTEP

//definition des enregistrements (ou records = variables de sortie intermediaires)
//NtotUdin: nombre d'udin fertiles present au debut du pas de temps + nouveaux (=potentiel de germination du pas de temps) 
//NIncRecSen :nombre d'incubation recouvertes par la senescence durant le pas de temps
//SIncRecSen,SChloRecSen: surface en incubation (ou S Chlo) recouverte par la senecence durant le pas de temps
//NIncRecChlo: nombre d'incubations recouvertes par les chloroses durant le pas de temps
//SIncRecChlo: surface en incubation  recouverte par les chlorose durant le pas de temps
//NIncRecInc,SIncRecInc: nombre et surface des inc eliminées par la croisssance des inc
//newNInc,newSInc : nouvelles incubations du pas de temps (apres maladie)
#define NLESREC 11
enum LesRec {NtotUdin,newNInc,newSInc,NIncRecSen,SIncRecSen,SChloRecSen,NIncRecChlo,SIncRecChlo,NIncRecInc,SIncRecInc,iSLesInSen};
const char *LesRecNames[]={"NtotUdin","newNInc","newSInc","NIncRecSen","SIncRecSen","SChloRecSen","NIncRecChlo","SIncRecChlo","NIncRecInc","SIncRecInc","SLesInSen"};


// ****************** Definition des parametres du cycle

typedef struct PARCYCLE {

  // --- Germination et penetration
  //taux de perte des Udin quand elles ne germent pas (fraction de l'effectif initial perdu par jour)
  float txPerteUdin;
  // % Humidité air mini requis pour avancer le processus de germination (= permettant  98 % humidite sur les feuilles)
  float HrGerm;
  // Rayonnement incident (micromolPAR) maxi pour avancer le processus de germination ( ie permettant 98 % humidite sur les feuilles)
  float PARgerm;
  //Dure minimale continu à HRGERM pour obtenir une germination (heures)
  float DminGerm;
  //proba de reussite (penetration + incubation) pour les UDIN 
  float pinc;

  // --- Parametrisation des reponses a la temperature post-germmination
  // T(°C) de base
  float TbasePath;
  // T(°C) max
  float TmaxPath;
  // intervalle degree-jour definissant une cohorte de meme age (largeur classe d'age)
  float Lcdd;

  // --- Incubation
  // Flag(0/1) pour la prise en compte de l'effet humidite (0 = avec effet humidite, 1 = sans effet humidite (latence fixe))
  int FixedLat;
  // Duree (dd) d'incubation (si fixe) ou duree max d'incubation (si fonction humidité)
  float DInc;
  // Parametre de reponse de la duree d'incubation a l'humidite (Rapilly 74) : del2(duree mini d'incub en condition humide) et fdhr
  float del2;
  float fdhr;
  // Surface de la lesion au debut de l'incubation
  float SIncMin;
  // taux de croissance des lesions en incubation
  float txInc;

  // --- Phase Chlorotique
  //Duree (dd) de la phase chlorotique
  float DChlo;
  //Surface minimale d'une lesion chlorotique (cm2)
  float Slmin;
  //Surface maximale d'une lesion chlorotique (cm2)
  float Slmax;
  //Taux croissance des lesions chlorotique en cm2/dd
  float txCroi;
  
  // --- Sporulation / Emissions goutelletes infectieuses (Eclins)

  //nombre d'evenements splashants pour arriver au vidage d'une lesion (3 pour rapilly)
  int nbEvSplaMax;
  // proportion de gouttelettes non infectantes en raison de contenu en spores trop faible. Une valeur par evenement splashant
  float *pNoSpo;
    //proportion de gouttelettes non infectante en raison de l'evaporation (Rapilly, 1976)
  double pEv;


  //--- Production de Spores (utilisé pour info/test) : n'affecte pas le calcul des Eclins
  // densite de stomates sur les surfaces sporulantes
  float dstom;
  //Proportion de stomates infectés par des pycnides dans une lesion
  float pStoInf;
  //Nbre total de spores produits par une pycnide (splashées+residuelles, cf Eyal 1971)
  int NbTotSpore;
  //proportion de spores résiduelles par pycnide (Eyal)
  float pSporeRes;
  //proportion de spore splashées (parmis les splashable) par les evenement splashant (il en faut autant que NBEVSPLA)
  float *pSpoEv;

  //condition initiale epidemie
  float LaiSpoSol;
  float TTinfSol; 

} ParCycle;

// Creation parametre septo a partir des valeurs de Parameter.h

// 'lecture' des parametre donnees sous forme de Tableau
float pNoSpoInPar[NBEVSPLAMAX] = {PNOSPO};
float pSpoEvInPar[NBEVSPLAMAX] = {PSPOEV};

ParCycle getParCycle() {
  ParCycle p;
  p.txPerteUdin = TXPERTEUDIN;
  p.HrGerm = HRGERM;
  p.PARgerm = PARGERM;
  p.DminGerm = DMINGERM;
  p.pinc = PINC;
  p.TbasePath = TBASEMYC;
  p.TmaxPath = TMAXMYC;
  p.Lcdd = LARGCLASS;
  p.FixedLat = FIXED_LAT;
  p.DInc = (FIXED_LAT == 1 ? DINC : DINCMAX);
  p.del2 = DINCMIN;
  p.fdhr = FDHR;
  p.SIncMin = MININC;
  p.txInc = (FIXED_LAT == 1 ? (MINLES - MININC) / DINC : (MINLES - MININC) / DCROIINC);
  p.DChlo = (FIXED_LAT == 1 ? DLAT - DINC : DCHLO);
  p.Slmin = MINLES;
  p.Slmax = MAXLES;
  p.txCroi = TXCROI;
  p.nbEvSplaMax = NBEVSPLAMAX;
  p.pNoSpo = (float*) calloc(NBEVSPLAMAX, sizeof(float));
  p.pEv = PEV;
  p.dstom = DSTOM;
  p.pStoInf = PSTOINF;
  p.NbTotSpore = NBTOTSPORE;
  p.pSporeRes = PSPORERES;
  p.pSpoEv = (float*) calloc(NBEVSPLAMAX, sizeof(float));
  p.LaiSpoSol = LAISPOSOL;
  p.TTinfSol = TTINFSOL;
  // recopie Tableaux
  int i;
  for (i = 0; i < NBEVSPLAMAX; i++) {
    p.pNoSpo[i] = pNoSpoInPar[i];
    p.pSpoEv[i] = pSpoEvInPar[i];
  }  
  return(p);
}

ParCycle newParCycle(void) {
	ParCycle p;
	p.pNoSpo = (float*) calloc(NBEVSPLAMAX, sizeof(float));
	p.pSpoEv = (float*) calloc(NBEVSPLAMAX, sizeof(float));
	return(p);
}	


void freeParCycle(ParCycle *p) {
  free(p->pNoSpo);
  free(p->pSpoEv);
}	

void copyParCycle(ParCycle *source, ParCycle *dest) {
  dest->txPerteUdin = source->txPerteUdin;
  dest->HrGerm = source->HrGerm;
  dest->PARgerm = source->PARgerm;
  dest->DminGerm = source->DminGerm;
  dest->pinc = source->pinc;
  dest->TbasePath = source->TbasePath;
  dest->TmaxPath = source->TmaxPath;
  dest->Lcdd = source->Lcdd;
  dest->FixedLat = source->FixedLat;
  dest->DInc = source->DInc;
  dest->del2 = source->del2;
  dest->fdhr = source->fdhr;
  dest->SIncMin = source->SIncMin;
  dest->txInc = source->txInc;
  dest->DChlo = source->DChlo;
  dest->Slmin = source->Slmin;
  dest->Slmax = source->Slmax;
  dest->txCroi = source->txCroi;
  dest->nbEvSplaMax = source->nbEvSplaMax;
  dest->pEv = source->pEv;
  dest->dstom = source->dstom;
  dest->pStoInf = source->pStoInf;
  dest->NbTotSpore = source->NbTotSpore;
  dest->pSporeRes = source->pSporeRes;
  dest->pSpoEv = source->pSpoEv;
  dest->LaiSpoSol = source->LaiSpoSol;
  dest->TTinfSol = source->TTinfSol;
  // recopie Tableaux
  int i;
  for (i = 0; i < dest->nbEvSplaMax; i++) {
    dest->pNoSpo[i] = source->pNoSpo[i];
    dest->pSpoEv[i] = source->pSpoEv[i];
	}
   return;
}




// ****************** Definition des Lesions

// NewUdin : Cohorte des  Unites de dissemination a deposer durant le prochain pas de temps (prochain appel a calcgerm)
//  (N) = Effectif de la cohorte (sans les pertes liees a txPerteUdin)
//  (Q) = Estimation du nombre de spores total contenues dans ces udin (pas utilisé)
//  (age) = duree (j) entre le debut du pas de temps et le depot 

// PopUdin : Unites de dissemination (goutes infectieuses capables de produire 1 lesion) deja deposees
//  (N) = Effectif initial de la cohorte (sans les pertes liees a txPerteUdin)
//  (Q) = Estimation du nombre de spores total contenues dans ces udin (pas utilisé)
//  (age) = Nombre de jour depuis le depot 
//  (borne sup) = Duree (j) de vie max des udin (1 / TxPerteUdin)


// PopInc : lesions en incubation 
//  (N) = effectif total de la cohorte
//  (Q) = surface totale de la cohorte (N * surface individuelle)
//  (age) = temps thermique depuis la germination
//  (borne sup) = duree (dd) d'incubation ou duree max d'incubation (si incubation variable avec humidite)

// PopHlat : structure utilisee pour le calcul des temps de latence variable uniquement. Similaire a PopInc et evoluant de la meme facon, mais stockant les informations supplementaires necessaires au calcul: 
// (N) = nombre d'heures ecoulees depuis le debut d'incubation
// (Q) = nombre d'heures favorables a la latence ecoulees depuis le debut de l'incubation

// PopCroi : lesions en croissance (avec production de chlorose)
//  (N) = effectif total de la cohorte
//  (Q) = surface d'une lesion de la cohorte (chloroses + necroses)
//  (age) = temps thermique depuis la fin d'incubation
//  (borne sup) = durée de croissance max des chloroses ((Slmax - Slmin) / txcroi)

// PopMat : lesions en maturation (sans production de nouvelles chloroses, mais avec continuation de la transformation chloroses -> surfaces sporulante)
//  (N) = effectif total de la cohorte
//  (Q) = surface d'une lesion de la cohorte (Chloroses + necroses). Vaut Slmax pour les lesions n'ayant pas subi de competition, moins sinon.
//  (age) = temps thermique depuis la derniere chloroses produite
//  (borne sup) = duree de la phase chlorotique (= duree necessaire a la maturation de la derniere chlorose produite)

// A la fin de la maturation, l'effectif des lesions est ajoute a NbLesAd

// Les proportions respectives de surfaces chlorotiques, sporulantes et vides a l'interieur des lesions ne sont pas calculees lesion par lesion, mais globalement pour l'ensemble des lesions du secteur. On respecte cependant la structuration en cohorte d'age pour ces tissus. Deux 'populations' stokent ces infos: 

// PopChlo : cohortes de surfaces chlorotiques du secteur classees par age
//   (N) Nombre de lesions en croissance de meme age que la cohorte (= 0 si lesions en croisssance plus agees)
//   (Q) Surface totale de la cohorte pour le secteur
//   (age) temps thermique depuis l'apparition de la chlorose
//   (borne sup) duree de la phase chlorotique

// PopSpo : cohortes de surfaces sporulantes
//   (N) : <!> age (dd depuis naissance) de la plus vielle surface de la cohorte (pour calcul surfaces lesions par lesions)
//   (Q) : surface totale de cette cohorte pour le secteur
//   (age) Nombre cumulé d'episode d'emission de mucilage ( = nbre de pluies recues)
//   (borne sup) Nbre max d'episode sporulation avant vidage

// Les surfaces ne pouvant plus emettre de spores sont ajoutees a SEmpLes.


typedef struct LESIONS {
  //Etat de la maladie (cf ci dessus)
  ClassMyc NewUdin;
  PopMyc PopUdin,PopInc,PopHlat,PopCroi,PopMat,PopChlo,PopSpo;
  int nbLesAd;
  float SEmpLes;
  
  //variable pour calcul
  //surfaces des lesions inclues dans la surface occupee par senescence naturelle dans temoin = surface necrotique (spo et empty) car surface chloro et inc recouvertes deviennet de la sene nat 
  //Si curseur : on doit rajouter Schlo et Sinc, pourait etre un tableau distinguant Schlo(dans le cas 'chlo necrotrophe'),Snec et Sempty
  float SLesInSen;
  // nombre d'heure favorable a la germination deja accumulees par popUdin 
  int nhGerm;  
  // indicateur pour gestion iteration.
  ///dTTmax is the duration of the longuest admissible timestep, that ensure the cohorts to make at max one phase change. It is equal to the duration of the shortest phase. 
  //skipXXX flags allows to skip one or several phase if phase duration is zero, independently of dTTmax (0 = don't skip, 1 = skip)
  
  float dTTmax; 
  int skipInc;
  int skipChlo;
  int skipCroi;

  //enregistrements de variables intermédiaires pour analyse en sortie
  int nLesRec;
  float LesRec[NLESREC];
  //parametres du cycle
  ParCycle par;
} Lesions;



// ****************** fonctions pour l'initialisation/mise a jour lors d'ajout


//raz des enregistrements

void resetLesRec(Lesions *les) {
  int i;
  for (i=0;i<les->nLesRec;i++)
    les->LesRec[i]=0;
}


// initialisation lesion

Lesions newLes(ParCycle *par) {
 
  Lesions les;
  
  les.NewUdin = newCmyc();

  //PopUdin: description de la phase entre depot des UDIN et germination (temps = j)
  float supUdin = (par->txPerteUdin > 0 ? (float) ceil(1. / par->txPerteUdin) : 1.);
  les.PopUdin=newPopMyc(0, supUdin, 1., pTUdin);
  //PopInc et PopHlat: description de la phase d'incubation (temps = thermal time)
  les.PopInc = newPopMyc(0, par->DInc, par->Lcdd, pTLes);
  if (par->FixedLat == 0)
    les.PopHlat = newPopMyc(0, par->DInc, par->Lcdd, pTHlat);
  else
	les.PopHlat = newPopMyc(0,0,par->Lcdd,pTHlat);
  //PopChlo :description de la phase de chlorose (temps = thermal time)
  les.PopChlo=newPopMyc(0, par->DChlo, par->Lcdd, pTLes);
 //PopCroi : front de croisance des lesions chlorotique (temps = age depuis Slmin)
  float supCroi = (par->txCroi > 0 ? (par->Slmax - par->Slmin) * 1. / par->txCroi : 0);
  les.PopCroi=newPopMyc(0, supCroi, par->Lcdd, pTLes);
  //PopMat : maturation apres fin de croisssance
  les.PopMat=newPopMyc(0, par->DChlo, par->Lcdd, pTLes);
  //PopSpo : phase sporulante (temps = nbr evenement splash)
  les.PopSpo=newPopMyc(0, par->nbEvSplaMax, 1., pTSpo);

  les.nbLesAd=0;
  les.SEmpLes=0;

  les.SLesInSen=0;
  les.nhGerm=0; 

  //skip Phase flags: 0 = don't skip, 1 = skip
  les.skipInc = 0;
  les.skipChlo = 0;
  les.skipCroi = 0;   
  float Dmin = 9999;
  if (par->DInc <= 0)
    les.skipInc = 1;
  else 
    Dmin = MIN(Dmin,par->DInc);
  if (par->DChlo <= 0)
    les.skipChlo = 1;
  else 
    Dmin = MIN(Dmin,par->DChlo);
  if (supCroi <= 0) 
    les.skipCroi = 1;
  else
    Dmin = MIN(Dmin,par->DInc);
  // longuest admissible timestep
  les.dTTmax = Dmin;
  
  les.nLesRec = NLESREC;
  resetLesRec(&les);
  les.par = newParCycle();
  copyParCycle(par,&(les.par));
  return(les);
}

//reset

void ResetLes(Lesions *les) {
  ResetCmyc(&(les->NewUdin));
  ResetPop(&(les->PopUdin));
  ResetPop(&(les->PopInc));
  if (les->par.FixedLat == 0)
    ResetPop(&(les->PopHlat));
  ResetPop(&(les->PopCroi));
  ResetPop(&(les->PopMat));
  ResetPop(&(les->PopChlo));
  ResetPop(&(les->PopSpo));
  les->nbLesAd = 0;
  les->SEmpLes = 0;
  les-> SLesInSen = 0;
  les->nhGerm = 0;
  resetLesRec(les);  
  
}

// Add a cohorts of incubationg lesions to Chlo. age = agedd since transition to Chlo, N = number of lesiions, Q = total surface of the cohort

void addCtoChlo(Lesions *les, ClassMyc Cmyc) {

  if (Cmyc.N > 0) {
     
    // Surfaces
    if (les->skipChlo != 1)
      AddCtoPop(Cmyc,&(les->PopChlo));
    else {//DChlo == 0
      BarPop Bpop = getBarPop(0,&(les->PopSpo));
      //maj age = age de la plus vieile transition
      Bpop.N = MAX(Bpop.N, Cmyc.age); 
      Bpop.Q += Cmyc.Q;
      PutBarPop(Bpop, &(les->PopSpo));
    }
    
    // Lesions
    Cmyc.Q /= Cmyc.N;//passage surf tot -> surf par lesion
    if (les->skipCroi != 1) 
      AddCtoPop(Cmyc,&(les->PopCroi)); 
    else {
      if (les->skipChlo != 1) 
	AddCtoPop(Cmyc, &(les->PopMat));
      else //skip Croi & skip Chlo
	les->nbLesAd += (int) Cmyc.N;
    }
  }
}

// Add a cohort of germinating lesions to Inc. age = agedd since start of incubation, N = number of lesions, Q = total surface of the cohort
//nh, nhlat : number of hour(hour favorale to latency) elapsed since germination

void addCtoInc(Lesions *les,ClassMyc Cmyc, int nh, int nhlat) {
  if (Cmyc.N > 0) {
    if (les->skipInc != 1) {
      int n;
      if (les->par.FixedLat == 1)
	n = AddCtoPop(Cmyc, &(les->PopInc));
      else 
	n = AddCtoPops(Cmyc, nh, nhlat, &(les->PopInc), &(les->PopHlat));
      if (n < 0) 
	show("\a\n<Lesions.h> :Erreur : impossible d'ajouter Cmyc à PopInc : age Classe=%g, agemax pop=%g\n",Cmyc.age,les->PopInc.agemax);
    } else 
      addCtoChlo(les, Cmyc);
  }
}

  //ajout de nouvelles Udins dans NewUdin, lors de l'interception.age = delai(jour) entre le debut du pas de temps et le depot, N = Nombres d'udin, Q = nombres de spores

void addCtoUdin(Lesions *les,ClassMyc Cmyc) {
  if (Cmyc.N > 0) {
    float oldageUdin = les->NewUdin.age;
    float oldNudin = les->NewUdin.N;
    les->NewUdin.Q += Cmyc.Q;
    les->NewUdin.N += Cmyc.N;
    les->NewUdin.age = (oldageUdin * oldNudin + Cmyc.age * Cmyc.N) / (Cmyc.N + oldNudin);
  }
}

void addUdin(Lesions *les, float N, float Q, float age) {
  ClassMyc Cmyc = newCmyc_p(age,Q,N);
  addCtoUdin(les, Cmyc);
}

//set manuel de sporulant (alias setSpoLes)
void addQLes(Lesions *les,float Q) {
  ClassMyc Cmyc;
  Cmyc.Q = Q;
  Cmyc.N = MAX(1, floorf(Cmyc.Q / les->par.Slmax));
  Cmyc.age = 0;
  AddCtoPop(Cmyc,&(les->PopSpo));
  les->nbLesAd += (int) Cmyc.N;
}


// ******************  Fonctions pour destruction

// fonction liberation memoire allouee aux pop

void freeLes(Lesions *les) {
  freePop(&(les->PopUdin));
  freePop(&(les->PopInc));
  freePop(&(les->PopHlat));	
  freePop(&(les->PopCroi));
  freePop(&(les->PopMat));
  freePop(&(les->PopChlo));
  freePop(&(les->PopSpo));
}

// ****************** Affichage


PopMyc getPopUdin(Lesions *les,float h);
float SLesNec(Lesions *les);

void AfficheLesion(Lesions *les) {
  show("\n************Contenu de l'objet Lesion : \n\nnbLesAd=%d SLesInSen=%f SEmpLes=%f, SLesNec=%f nhGerm=%d\n",les->nbLesAd,les->SLesInSen,les->SEmpLes,SLesNec(les),les->nhGerm);
  show("\n************NewUdin:\n");
  PrintClass(les->NewUdin);
  show("\n************PopUdin:\n");
  PopMyc p = getPopUdin(les,0);
  PrintPop(p);
  freePop(&p);
  show("\n************PopInc:\n");
  PrintPop(les->PopInc);
  if (les->par.FixedLat == 0) {
    show("\n************PopHlat:\n");
    PrintPop(les->PopHlat);
  }
  show("\n************PopChlo:\n");
  PrintPop(les->PopChlo);
  show("\n************PopCroi:\n");
  PrintPop(les->PopCroi);
  show("\n************PopMat:\n");
  PrintPop(les->PopMat);
  show("\n************PopSpo:\n");
  PrintPop(les->PopSpo);
}


//affichage des variable d'enregistrement

void AfficheLesRec(Lesions *les) {
  show("\n***********Enregistrements de la lesion:");
  int i;
  for (i=0;i<les->nLesRec;i++)
    show("%s : %g, ",LesRecNames[i],les->LesRec[i]);
  show("\n");
}

// ****************** Get/Set

// get/set variables d'enregistrement

void setLesRec(Lesions *les,int key,float value) {
  les->LesRec[key]=value;
}

float getLesRec(Lesions *les,int key) {
  return(les->LesRec[key]);
}

void getLesRecKey(int key, char *str) {
  sprintf(str,LesRecNames[key]);
}
// get statistiques sur etat des lesions

float SLesInc(Lesions *les) {  
  return(Qtot(&(les->PopInc)));  
}

float NLesInc(Lesions *les) {
  return(Ntot(&(les->PopInc)));
}

float SLesSpo(Lesions *les) {
  return(Qtot(&(les->PopSpo)));
}

float SLesNec(Lesions *les) {
  return(Qtot(&(les->PopSpo)) + les->SEmpLes);
}

float SLesChlo(Lesions *les) {
  return(Qtot(&(les->PopChlo)));
}

float NLesChloNec(Lesions *les) {
  return(Ntot(&(les->PopCroi)) + Ntot(&(les->PopMat)));
}

float NLesNec(Lesions *les) {
 return(les->nbLesAd);
}

float SLesTot(Lesions *les) {
  return(SLesNec(les) + SLesChlo(les) + SLesInc(les));
}

float NLesTot(Lesions *les) {
  return(NLesInc(les) + NLesChloNec(les) + NLesNec(les));
}


//surface des lesions sur le vert
float SLesOnGreen(Lesions *les) {
  return SLesInc(les) + SLesChlo(les) + (SLesNec(les) - les->SLesInSen);//a changer si chloroses deviennet necrotrophe
}

//************** Reponses environnement



// Reponse a la temperature (temps thermique)

float Lesions_somT(int ideb, int ifin, float *T, ParCycle *par) {
  float sT=0;
  int i;
  for (i = ideb; i <= ifin; i++)
    if (T[i] <= par->TmaxPath)
      sT += MAX(0, T[i] - par->TbasePath);
  return(sT * 1. / 24.);
}


// nombre d'heures contigue favorable a la germination entre [debut ideb,  ifin[ dans des tableaux meteo

int nhGerm(int ideb, int ifin, float *PPFD, float *Rh, ParCycle *par) {
  int i = ideb;
  while (Rh[i] >= par->HrGerm && PPFD[i] <= par->PARgerm && i < ifin) 
      i++;
  return(i - ideb);
}

// Total nombre d'heures favorable a la latence (contigue ou non)  entre [ideb et ifin]

int nhLat(int ideb, int ifin, float *PPFD, float *Rh, ParCycle *par) {
  int res = 0;
  int i;
  for (i = ideb; i <= ifin; i++)
    if (Rh[i] >= par->HrGerm && PPFD[i] <= par->PARgerm)
      res++;
  return res;
}

// Duree  temps de latence dd fonction du ratio nbheure / nbheure favorable a la latence
// formule Rapilly 1974)

double Dlat(float nbh,float nbhlat,ParCycle *par) {
  double ddr = 0;
  if (nbh > 0 & nbhlat > 0)
    ddr = log10(nbhlat / nbh);
  return par->del2 + ddr / par->fdhr * par->TmaxPath;
}

//Taux de satisfaction (<1 si demande > offre, 1 sinon)

float TxSat(float offre,float demande) {
  float txsat = 1;
  if (demande != 0) 
    txsat = MAX(0, MIN(1, offre/demande));
  return txsat;
}


//Renvoi l'effectif fertile d'une cohorte d'udin, h heures apres/avant le debut du pas de temps. age  = temps (j) ecoule depuis le depot au debut du pas de temps

float fertUdin(ClassMyc *Cudin,float h,ParCycle *p) {
  float res = Cudin->N;
  if (Cudin->age + h / 24. > 0)// h et/ou age peuvent etre negatif
    res = propEntiere(Cudin->N, 1 - p->txPerteUdin * (Cudin->age + h / 24.)); 
  return(res);
}


//************ Germination / traitement des udins


//renvoie une copie de popudin,en ayant effectue les corrections de perte d'udin

PopMyc getPopUdin(Lesions *les,float h) {
  PopMyc udin = duplicatePopMyc(&(les->PopUdin));
  BarPop bp;
  ClassMyc Cmyc;
  float newN;
  int i;
  for (i = 0; i < udin.nbClass; i ++) {
    bp = getBarPop(i, &udin);
    if (bp.N > 0) {
      Cmyc = getCmyc(i,&udin);
      newN = fertUdin(&Cmyc,h,&(les->par));
      bp.Q = bp.Q / bp.N * newN;
      bp.N = newN;
      PutBarPop(bp,&udin);
    }
  }
  return(udin);
}
  
     
// calcul des germination pour les prochaines dth heures, mise a jour des Udins et nhGerm.
// variables meteo pour les prochaines dht heures: PPFD et Rh 
//Renvoie une cohorte avec N =  effectif d'udin ayant aquis Dmin heures germantes, Q = unused (spores) et age = delai (j) debut timestep -> debut germination (fin Dmin). 

ClassMyc calcGerm(Lesions *les,int dth, float *PPFD, float *Rh, int debug) {

  int i,n;
  float hdeb,hfin,Ngerm;
  ClassMyc Cmyc;
  ClassMyc Cgerm = newCmyc();
  
  //calcul et enregistrement du potentiel de germination
  float maxGerm = les->NewUdin.N;
  // on prends les effectifs au debut de la deniere periode humide, ie il y a nhGerm heures ou depuis depot si inferieur
  for (i = 1; i < les->PopUdin.nbClass; i++) {
    Cmyc = getCmyc(i,&(les->PopUdin));
    if (Cmyc.N > 0) {
      hdeb = - MIN(les->nhGerm, Cmyc.age * 24);
      maxGerm += fertUdin(&Cmyc,hdeb,&(les->par));
    }
  }
  setLesRec(les,NtotUdin,maxGerm);

  if (debug)
    show("\nCalcGerm : nhGerm = %d (DminGerm = %.0f), maxGerm (#viable udins at beginning of time step) = %.0f\n", les->nhGerm,les->par.DminGerm,  maxGerm);

  // calcul germinations
  float DminGerm = les->par.DminGerm;
  int ideb = 0;
  int Dgerm = 0;
  while (ideb + MIN(Dgerm,DminGerm) < dth) {// sinon germination debut pas de temps suivant
    Dgerm = nhGerm(ideb, dth, PPFD, Rh, &(les->par)) ;
    if (Dgerm + les->nhGerm >= DminGerm && maxGerm > 0) {//conditions a minima
      if (debug)
	show("CalcGerm: found potential germinating period starting at %d, lasting %d hours (given nhGerm = %d)\n",ideb,Dgerm, les->nhGerm);
      for (i = -1; i < les->PopUdin.nbClass; i++) {
	if (i == -1) {
	  Cmyc = les->NewUdin;
	  hdeb = (float) MIN(dth,MAX(ideb,les->NewUdin.age * 24));//age newUdin = delai appel calcGerm -> depot
	  hfin = MIN(dth, hdeb + MIN(Dgerm, DminGerm) - 1);// -1 car hdeb inclus
	  Cmyc.age = - les->NewUdin.age;// passage delai debut pas temps -> depot a temps depuis le depot au debut du pas de temps
	  if (debug)
	    if(Cmyc.N > 0)
	      show("CalcGerm: found newUdin (N = %.0f, age =%g, hdeb = %.0f, hgerm = %.0f)...\n",Cmyc.N, Cmyc.age, hdeb,hfin);
	} else {
	  Cmyc = getCmyc(i,&(les->PopUdin));
	  if (ideb == 0)
	    hdeb = - MIN(les->nhGerm, Cmyc.age * 24);// hdeb negative car avant pas de temps en cours
	  else
	    hdeb = (float) ideb;
	  hfin = hdeb + MIN(Dgerm, DminGerm) - 1;
	  if (debug)
	    if (Cmyc.N > epsillon)
	      show("CalcGerm: found Udin (N = %g, age =%g, hdeb = %.0f, hgerm = %.0f)...\n",Cmyc.N, Cmyc.age, hdeb,hfin);
	}
	if (Cmyc.N > epsillon && (hfin - hdeb + 1) >= DminGerm) {
	  Ngerm = fertUdin(&Cmyc, hdeb,&(les->par));// Tx perte s'applique au debut de la periode favorable
	  if (Ngerm > 0) {
	    if (debug)
	      show("CalcGerm: Cohort will germinate at %.0f hour (Ngerm=%.0f), ie age Cgerm = %.2f, updating Cgerm...\n",hfin,Ngerm,hfin / 24.);
	    Cgerm.age = (Cgerm.age * Cgerm.N + hfin / 24. * Ngerm) / (Cgerm.N + Ngerm);
	    Cgerm.N += Ngerm;
	    Cgerm.Q += Cmyc.Q / Cmyc.N * Ngerm;
	    if (i >= 0)// classe issue de popUdin
	      ResetClass(i, &(les->PopUdin));
	    else
	      les->NewUdin = newCmyc();
	  }
	}
      }
    }
    ideb += MAX(1,Dgerm);
  }

  
  //Vieillissement udin
  PopMyc p = AgePop(&(les->PopUdin), dth / 24.);
  freePop(&p);
  

  // ajout / vieillissement NewUdin
  if (les->NewUdin.N > 0) {//la germination n'a pas eu lieu
    les->NewUdin.age = dth / 24. - les->NewUdin.age;//passage age = delai debut pas de temps - depot a age = temps depuis le depot
    if (les->NewUdin.age >= 0) {//si temps depuis le depot positif a la fin du pas de temps
      AddCtoPop(les->NewUdin,&(les->PopUdin));
      les->NewUdin = newCmyc();
    } else {
      les->NewUdin.age = - les->NewUdin.age; // on repasse en temps debut prochain appel calcGerm -> depot
    }
  }

  //maj nhGerm
  // on se replace au moment de la derniere iter
  ideb -= MAX(1,Dgerm);
  if (ideb == 0)
    les->nhGerm += Dgerm;//periode avant pas de temps toujours en cours
  else
    les->nhGerm = Dgerm;

  return(Cgerm);
}
 

// fonction decomposant un pas de temps pour accomoder deux periodes (avant/apres germination) si germination a et redecomposant eventuellement en sous iterations pour etre compatibles avec dTTmax.igerm donne le numero de sous-iteration APRES laquelle il faut transformer les germinations en incub (-1 si pas de germ)

void splitTimeStep(Lesions *les,int dth,ClassMyc *Cgerm, float *T, int *nsubstep, int *igerm,int ideb[],int ifin[],float dTT[]) {
  
  int hgerm = -1;

  if (Cgerm->N > 0) {
    hgerm = (int) round(Cgerm->age * 24);
    if (hgerm >= dth)
      show("<Lesions.h> splitTimeStep : Warning !! hgerm > timestep duration : CalcGerm is to be checked\n");
    hgerm = MIN(dth - 1,hgerm);
  }

  //decomposition
  int hdeb,hfin;
  float dTTit;
  int nstep = 0;
  int h = 0;
  
  while (h < dth) {
    hdeb = h;
    if (hdeb < hgerm)
      hfin = hgerm;
    else
      hfin = dth - 1;
    dTTit = Lesions_somT(hdeb,hfin,T,&(les->par));
    while (dTTit > les->dTTmax && (hfin - hdeb) > 1) {
	hfin--;
	dTTit = Lesions_somT(hdeb,hfin,T,&(les->par));
    }
    ideb[nstep] = hdeb;
    ifin[nstep] = hfin;
    dTT[nstep] = dTTit;
    if (dTTit > les->dTTmax)
      show("<Lesions.h> : Warning !! can't acomodate dttmax with hourly substeps\n");
    nstep++;
    h = hfin + 1;
  }
  
  *nsubstep = nstep;
  *igerm = hgerm;
}

void showIter(int nstep, int igerm, int *ideb, int *ifin, float* dTT,float dTTmax) {

  show("\n Decomposition de l'iteration (substeps = %d, igerm = %d, dTTmax = %.1f)\n",nstep, igerm,dTTmax);
  show("\niter\tideb\tifin\tdTT");
  int i;
  for (i = 0; i < nstep; i++)
    show("\n%d\t%d\t%d\t%.1f",i,ideb[i],ifin[i],dTT[i]);
  show("\n");

}
 
// *********** Gestion incub


// Create new incubations induced by a Germinating cohort. No growth is computed as this function should be called at the end of a sub-iteration (see splitTimeStep)

void GermToInc(Lesions *les, ClassMyc *Cgerm, float *Svert, float protectant_efficacy, int debug) {
  if (debug)
    show("GermInc : try to make germination of Cgerm (N=%.0f, Q= %.0f), pinc = %g\n", Cgerm->N, Cgerm->Q, les->par.pinc); 
  // potentiel de germination : PINC * Cgerm
  ClassMyc newInc = *Cgerm;
  newInc.age = 0;
  newInc.N = propEntiere(Cgerm->N, les->par.pinc * (1 - MAX(0, MIN(1, protectant_efficacy))));
  //surface : SIncMin ou Slmin
  if (les->skipInc != 1)
    newInc.Q = newInc.N * les->par.SIncMin;
  else
    newInc.Q = newInc.N * les->par.Slmin;
  // calcul tx de satisfaction
  float txsat = TxSat(*Svert, newInc.Q);
  if (txsat < 1) {
    newInc.N = propEntiere(newInc.N, txsat);
    newInc.Q = newInc.N * les->par.SIncMin;
  }
  //maj des surfaces
  if (newInc.N > 0) {
    if (debug)
      show("Germination succeed ! NewInc : N = %.0f, Q = %g\n",newInc.N, newInc.Q); 
    *Svert -= newInc.Q;
    addCtoInc(les,newInc,0,0);
  }
}

//delta SInc d'une lesions durant dtTT.Q : taille courrante. 
//Si duree incub repond a l'humidite, lorsque duree incub < DCroiInc, on peut avoir des sufaces<Slmin en fin d'incub (mais pas superieures)
//dans tous les cas, les nouvelles lesions auront pour taille Slmin
float dSInc(float dtTT,float Q,ParCycle *p) {
  return(MIN(p->txInc * dtTT, MAX(0, p->Slmin - Q)));
}

//liberation d'espace pour la croissance des incubs

float freeSpaceIncC(Lesions *les, StackBpop *IncC, float *Svert,float SDem,float *cumSrec,int debug) {

  float SvertDisp = *Svert;
  float SincCDisp = Qdisp_sbp(IncC);// lesions ne s'etant pas encore developpees
  float SDisp = SvertDisp + SincCDisp;
  float txsat = TxSat(SDisp, SDem);

  float Sko = 0;
  float Skv = 0;

  if (txsat > 0) {
    // recouvrement dans incC
    Sko = MIN(SincCDisp, SDem * SincCDisp / SDisp * txsat);
    if (Sko > 0)
      Sko = killInStackBp(IncC, Sko);
    // recouvrement du vert
    if (SDem > Sko) 
      Skv = MIN(SvertDisp,SDem - Sko);
    // gestion excedent eventuel => retransformation en vert
    Skv = MIN(SDem,Sko + Skv) - Sko;
    *Svert -= Skv;
  }
  *cumSrec += Sko; 
  return(Sko + Skv);
}

// Croissance des incubs 

void CroiInc(Lesions *les, float *Svert, float dTT, float *NIncRec, float *SIncRec, int debug) {

  float Srec =0;

  // classes candidates a la croissance
  StackBpop incC = PopAsStack_bp(&(les->PopInc));
  int nincC0 = incC.nclass;
  
  // Croissance, de la plus ancienne a la plus recente
  if (incC.nclass > 0) {
    int i;
    BarPop Bpop;
    float SDem,Sdone;
    for (i = 0; i < incC.nclass; i++) {
      Bpop = incC.stack[i];
      SDem = 0;
      incC.ilim = i;// une lesion ne peut s'etendre sur elle meme
      if (Bpop.N > 0)
	SDem = dSInc(dTT,Bpop.Q / Bpop.N,&(les->par)) * Bpop.N;//SDem = 0 si plus en croisssance
      if (SDem > 0) {
	Sdone = freeSpaceIncC(les, &incC, Svert, SDem, &Srec, debug);
	incC.stack[i].Q += Sdone;
      }
    }

    // maj PopInc/PopHlat
    for (i = 0 ; i < incC.nclass ; i++)
      PutBarPop(incC.stack[i], &(les->PopInc));
    // effacement classes disparues
    if (incC.nclass < nincC0)
      for (i = incC.nclass + 1; i < nincC0; i++) {
	ResetBpop(incC.stack[i].pos,&(les->PopInc));
	ResetBpop(incC.stack[i].pos,&(les->PopHlat));
      }
	
  }

  *NIncRec += nincC0 - incC.nclass;
  *SIncRec += Srec;
  freeStackBpop(&incC);
}

// *********** Gestion Chloroses

//delta SChlo d'une lesions durant dtTT.Q : taille courrante

float dSChlo(float dtTT,float Q,ParCycle *p) {
  return(MIN(p->txCroi *dtTT, MAX(0, p->Slmax - Q)));
}

//liberation d'espace (Vert/Inc) pour la croissance. retourne l'espace liberee

float freeSpace(Lesions *les, float *Svert,float SDem,int debug) {

  float SvertDisp = *Svert;
  float SIncDisp = Qtot(&(les->PopInc));
  float SDisp = SvertDisp + SIncDisp;
  float txsat = TxSat(SDisp, SDem);

  float Sk = 0;
  float Skv = 0;

  if (txsat > 0) {
    
    // tentative recouvrement dans PopInc
    Sk = MIN(SIncDisp, SDem * SIncDisp / SDisp * txsat);
    if (Sk > 0)
      Sk = RemoveQbyNd(Sk, &(les->PopInc), SIncDisp / SDisp, debug);// Sk : surface effectivement gagnee sur les inc
    // recouvrement du vert
    if (SDem > Sk) {
      // ajustement txsat (gestion arrrondi kill inc favorable)
      SDisp = SvertDisp + Sk;
      txsat = TxSat(SDisp, SDem);
      Skv = MIN(SvertDisp, SDem * SvertDisp / SDisp * txsat);
      // gestion excedent eventuel
      Skv = MIN(SDem,Sk + Skv) - Sk;
    }
    *Svert -= Skv;
  }
  
  return(Sk + Skv);
}


void CroiChlo(Lesions *les,float *Svert, float dTT,int debug) {
  
  if (Qtot(&(les->PopCroi)) > 0){

    ClassMyc Cmyc;
    BarPop Bpop;
    int i;

    // calcul demande de surface
    float SDem = 0;
    for (i = 0; i < les->PopCroi.nbClass ; i++) {
      Cmyc = getCmyc(i, &(les->PopCroi));
      SDem += dSChlo(dTT, Cmyc.Q, &(les->par)) * Cmyc.N;
    }
    // liberation de la surface sur Inc + vert
    float Sdone = freeSpace(les, Svert, SDem, debug);

    // Realisation croissance
    if (Sdone == 0) {// plus de croissance possible, on vide tout
      Cmyc = ClassEqQ(&(les->PopCroi));
      ResetPop(&(les->PopCroi));
      if (les->skipChlo != 1) {
	Cmyc.age=0;
	AddCtoPop(Cmyc,&(les->PopMat));
      } else // DChlo = 0
	les->nbLesAd += (int) Cmyc.N;
      
    } else {// Sdone != 0, croissance
 
      //Creation nouvellle surface Chlo
      float newS;
      for (i = 0; i < les->PopCroi.nbClass ; i++) {
	Bpop = getBarPop(i, &(les->PopCroi));
	if (Bpop.N > 0) {
	  //maj Q dans popCroi
	  newS = dSChlo(dTT, Bpop.Q, &(les->par)) * Bpop.N * Sdone / SDem;
	  Bpop.Q += newS / Bpop.N;
	  PutBarPop(Bpop, &(les->PopCroi));
	  // maj popChlo
	  Cmyc = newCmyc_p(dTT, newS, 0);
	  if (les->skipChlo != 1)
	    AddCtoPop(Cmyc, &(les->PopChlo));
	  else {
	    Bpop = getBarPop(0,&(les->PopSpo));
	    //maj age = age de la plus vieile transition
	    Bpop.N = MAX(Bpop.N, Cmyc.age); 
	    Bpop.Q += Cmyc.Q;
	    PutBarPop(Bpop, &(les->PopSpo));
	  }
	}
      }

      //Vieillissement PopCroi
      PopMyc out = AgePop(&(les->PopCroi),dTT);
      for (i = 0 ; i < out.nbClass ; i++) {
	Cmyc = getCmyc(i,&out);
	if (Cmyc.N > 0) {
	  if (les->skipChlo != 1)
	    AddCtoPop(Cmyc,&(les->PopMat));
	  else
	    les->nbLesAd += (int) Cmyc.N;
	}
      }
      freePop(&out);
    }
  } else //Qtot(PopCroi = 0)
    if (debug)
      show("CroiChlo : PopCroi Empty !\n");
}

// renvoi une pop de Inc candidates a entree en Chlo. Met a jour age des Inc restante (sans evolution de leur croissance)

PopMyc getNewChlo(Lesions *les, int hdeb,int hfin, float *T, float *PPFD, float *Rh, float eradicant_efficacy) {
  
  int i;
  ClassMyc Cmyc;
  BarPop Bpop;
  

  float dTT = Lesions_somT(hdeb,hfin,T,&(les->par));
  if (eradicant_efficacy > 0) {	
		dTT *= (1 - MAX(0, MIN(1, eradicant_efficacy)));	// a tester que les fonctions acceptent bien  toute un dTTit = 0
	}
  PopMyc outInc = AgePop(&(les->PopInc), dTT);
    
  //----gestion cas dependance a l'humidite
  if (les->par.FixedLat < 1) {

    int dth = hfin - hdeb;
    int dhlat = nhLat(hdeb,hfin,PPFD,Rh,&(les->par));
      
    les->PopHlat.TimeBuffer = les->PopInc.TimeBuffer;
    AgePop(&(les->PopHlat), dTT);
      
    if (Qtot(&outInc) > 0) {
      show("\n<Lesions.h> Error: Parameter DINCMAX too small to keep all kesions in latency!!!! : to be adjusted !!\n");
      exit(0);
    }
  
    //filtrage des fin d'incubation reele
    for (i = 0; i < (les->PopInc).nbClass ; i++) {
      Cmyc = getCmyc(i, &(les->PopInc));
      if (Cmyc.Q > 0) {
	Bpop = getBarPop(i, &(les->PopHlat));
	Bpop.N += dth;
	Bpop.Q += dhlat;
	PutBarPop(Bpop, &(les->PopHlat));
	if (Cmyc.age >= Dlat(Bpop.N,Bpop.Q,&(les->par))) {
	  //show("fin latence age= %.2f\n",Cmyc.age);
	  Cmyc.age -= Dlat(Bpop.N,Bpop.Q,&(les->par));
	  AddCtoPop(Cmyc, &outInc);
	  ResetClass(i, &(les->PopInc));
	  ResetClass(i, &(les->PopHlat));
	}
      }
    }
  }

  return(outInc);
}



// liberation d'espace pour realisation des Chloroses, combinant Svert, Stot(PopInc) et outInc
//renvoie la surface liberee (<= Sdem)

float freeSpaceNewChlo(Lesions *les, StackCmyc *outInc, float *Svert,float SDem,float *cumSrec,int debug) {


  float SvertDisp = *Svert;
  float SpopIncDisp = Qtot(&(les->PopInc));
  float SoutIncDisp = Qdisp_sc(outInc);// surface encore a l'etat d'incub dans outInc
  
  float SDisp = SvertDisp + SpopIncDisp + SoutIncDisp;
  float txsat = TxSat(SDisp, SDem);

  float Sk = 0;
  float Sko = 0;
  float Skv = 0;

  if (txsat > 0) {
    
    // tentative recouvrement dans PopInc
    Sk = MIN(SpopIncDisp, SDem * SpopIncDisp / SDisp * txsat);
    if (Sk > 0)
      Sk = RemoveQbyNd(Sk, &(les->PopInc), SpopIncDisp / SDisp, debug);// Sk : surface effectivement gagnee sur les inc
    // recouvrement dans outInc
    if (SDem > Sk) {
      // ajustement txsat (gestion arrrondi RemoveQbyN favorable)
      SDisp = SvertDisp + Sk + SoutIncDisp;
      txsat = TxSat(SDisp, SDem);
      Sko = MIN(SoutIncDisp, SDem * SoutIncDisp / SDisp * txsat);
      if (Sko > 0)
	Sko = killInStack(outInc, Sko);
      // recouvrement du vert
      if (SDem > (Sk + Sko)) {
	// ajustement txsat (gestion arrrondi killInStack favorable)
	SDisp = SvertDisp + Sk + Sko;
	txsat = TxSat(SDisp, SDem);
	Skv = MIN(SvertDisp, SDem * SvertDisp / SDisp * txsat);
	// gestion excedent eventuel Skilled (Skv devient negatif): on retransforme en vert, si besoin
	Skv = MIN(SDem,Sk + Sko + Skv) - (Sk + Sko);
      }
    }
    *Svert -= Skv;
  }
  *cumSrec += Sk + Sko; 
  
  return(Sk + Sko + Skv);
}
 



// Creation nouvelles Chlo pour candidates issues de Inc.

void IncToChlo(Lesions *les, float *Svert, int hdeb, int hfin, float *T, float *PPFD, float *Rh, float eradicant_efficacy, float *NIncRec, float *SIncRec, int debug) {
  
  float Srec = 0;
  float Nrec = 0;

  // Classes candidates sortie d'incubation, maj age des inc restant en incubation
  PopMyc popout = getNewChlo(les, hdeb, hfin, T, PPFD, Rh, eradicant_efficacy);
  StackCmyc outInc = PopAsStack(&popout);
  freePop(&popout);

  if (outInc.nclass > 0) {

    ClassMyc Cmyc;
    BarPop Bpop;
    int i,j; 
    float dTT,SDem,Sdone;

    float N0 = Ntot(&(les->PopInc)) + Ntot_sc(&outInc);
    
    // ************ 1 . On transforme progressivement outInc de la plus vieille cohorte a la plus recente
    for (i = 0; i < outInc.nclass ; i++) { 
      
      // pas de temps = delta temps avant transition next classe. age de la classe sinon
      Cmyc = outInc.stack[i];
      if (i == outInc.nclass - 1)// plus d'autre classe
		dTT = Cmyc.age;
      else
		dTT = Cmyc.age - outInc.stack[i + 1].age;
      
	  if (eradicant_efficacy > 0) {	
		dTT *= (1 - MAX(0, MIN(1, eradicant_efficacy)));
	   } 
	  
      // ___________maj Croissance des lesions ayant deja fait la transition (if any)
      if (i > 0 && les->skipCroi != 1) //sinon aucune cohorte chlorotique en croissance
	for (j = 0; j < i ; j++) {// on parcours uniquement des cohortes deja chlorotiques
	  SDem = 0;
	  if (outInc.stack[j].N > 0) 
	    SDem = dSChlo(dTT, outInc.stack[j].Q / outInc.stack[j].N, &(les->par)) * outInc.stack[j].N;
	  Sdone = freeSpaceNewChlo(les, &outInc, Svert, SDem, &Srec, debug);// cohorte i pas encore chlorotique
	  outInc.stack[j].Q += Sdone;
	}
      
      
      // ___________Transition Inc->Chlo
      Cmyc = outInc.stack[i];//rafraichissement cas ou classe a changee en raison des croissances
      if (Cmyc.N > 0) {
	//Transformation SInc de la lesion
	outInc.stack[i].N = MIN(Cmyc.N,floorf(Cmyc.Q / les->par.Slmin));
	outInc.stack[i].Q = outInc.stack[i].N * les->par.Slmin;
	outInc.ilim = i;
	//Cmyc garde le nombre de lesion restant a faire, et la surface Inc non utilise (car < Slmin)
	Cmyc.N -= outInc.stack[i].N;
	Cmyc.Q -= outInc.stack[i].Q;
	// essai de transition en colonisant espace libre
	if (Cmyc.N > 0) {
	  SDem = les->par.Slmin * Cmyc.N - Cmyc.Q;
	  Sdone = freeSpaceNewChlo(les, &outInc, Svert, SDem, &Srec, debug);
	  if (Sdone >= SDem - epsillon) {
	    outInc.stack[i].N += Cmyc.N;
	    outInc.stack[i].Q = outInc.stack[i].N * les->par.Slmin;
	  } else {// Tout ne peux pas pas sortir
	    float newN = MIN(Cmyc.N,floorf(Sdone / les->par.Slmin));
	    outInc.stack[i].N += Cmyc.N;
	    outInc.stack[i].Q = outInc.stack[i].N * les->par.Slmin;
	    // on restaure la surfInc restante en S vert si pas utilisee
	    if (newN == 0)
	      *Svert += Cmyc.Q;
	  }
	}
      }
      
      
      //____________Croissance depuis la transition
      if (les->skipCroi != 1) {
	SDem = 0;
	if (outInc.stack[i].N > 0) 
	  SDem = dSChlo(dTT, outInc.stack[i].Q / outInc.stack[i].N, &(les->par)) * outInc.stack[i].N;
	Sdone = freeSpaceNewChlo(les, &outInc, Svert, SDem, &Srec, debug);
	outInc.stack[i].Q += Sdone;
      }
    }

    // *********  2. Creation des Surfaces
    for (i = 0; i < outInc.nclass; i++) 
      addCtoChlo(les,outInc.stack[i]);
    float Nfin = Ntot(&(les->PopInc)) + Ntot_sc(&outInc);
    Nrec = N0 - Nfin;
  }

  *NIncRec += Nrec;
  *SIncRec += Srec;
  freeStackCmyc(&outInc);
}


// Vieillissement Chlo, passage chlorose -> necrose (spo)

void ChloToSpo(Lesions *les, float dTT) {
  ClassMyc Cmyc;
  BarPop Bpop;
  // maj PopChlo
  PopMyc outChlo = AgePop(&(les->PopChlo), dTT);
  // Ajout Spo
  Cmyc = ClassEqQ(&outChlo);
  if (Cmyc.Q > 0) {
    Bpop = getBarPop(0,&(les->PopSpo));
    //maj age = age de la plus vieile transition
    Bpop.N = MAX(Bpop.N, Cmyc.age); 
    Bpop.Q += Cmyc.Q;
    PutBarPop(Bpop, &(les->PopSpo));
  }
  freePop(&outChlo);
  // mise a jour des lesions en maturation
  PopMyc outMat = AgePop(&(les->PopMat), dTT);
  Cmyc = ClassEqQ(&outMat);
  les->nbLesAd += (int) Cmyc.N;
  freePop(&outMat);
}


//mise a jour age  popSpo (ou N = age)

void AgeSpo(PopMyc *PopSpo, float dTT) {
  int i;
  BarPop Bpop;
  for (i = 0; i < PopSpo->nbClass; i++) {
	Bpop = getBarPop(i,PopSpo);
	if (Bpop.Q > 0)
	  Bpop.N += dTT;
	PutBarPop(Bpop, PopSpo);
  }
}


// Production de spores / vidage des lesions lors des pluies
// EclaLes: nombre de goutelettes produites par la pluie sur la lesion


void RainOnLes(Lesions *les, double EclaLes, double *Eclins, double *QSpores) {
  
  *Eclins = 0;
  *QSpores = 0;
  
  int i;
  ClassMyc Cmyc;
  double ecla;
  double Sspo = SLesSpo(les);

  if (Sspo > 0) {
    for (i = 0; i < les->par.nbEvSplaMax ; i++) {
      Cmyc = getCmyc(i,&(les->PopSpo));
      ecla = Cmyc.Q / Sspo * EclaLes; 
      *Eclins += ecla * MAX(0, 1 - les->par.pEv - les->par.pNoSpo[i]);
      *QSpores += Cmyc.Q / Sspo * les->par.dstom * les->par.pStoInf * les->par.NbTotSpore * (1 - les->par.pSporeRes) * les->par.pSpoEv[i];
    }
    PopMyc popout = AgePop(&(les->PopSpo),1);
    les->SEmpLes += Qtot(&popout);
    freePop(&popout);
  }
}


// cond ini via LaiSpoSol

double getEclinSoil(float TTem,double EclaSol, ParCycle *p) {
  double res = 0;
  double LAISpo =  (TTem >= p->TTinfSol ? 0 : p->LaiSpoSol);
  res = MIN(1,LAISpo) * EclaSol * MAX(0, 1  - p->pEv - p->pNoSpo[0]);
  return(res);
}


//Macro simulation developpement d'une lesion avec un espace (vert) libre donnee durant dth heures. 
// Var meteo horaires durant dth: Temperature (T, °C), PAR Photon Flux Density (PPFD, micromol.m2.s-1), Relative humidity (Rh, %)
// the cycle respond to pesticide via efficacy of products concept (see milne et all); Efficay = 1 means max efficacy)
// On neglige ici le fait que, sur les bordures du secteurs, les lesions ont la possibilité de s'etendre sur le secteur voisin. Ceci est valable si la surface des bordures (largeur secteur*sqrt(Slmax)/2) est petite devant la surface du secteur
//on met a jour Necroses puis Chloroses puis Incubations puis pour gerer les priorites d'acces au vert

void devLes_dp(Lesions *les,float *Svert,int dth,float *T,float *PPFD, float *Rh, float protectant_efficacy, float eradicant_efficacy, int debug) {

  //variables pour enregistrement
  float sumNIncRecChlo = 0;
  float sumSIncRecChlo = 0;
  float sumNIncRecInc = 0;
  float sumSIncRecInc = 0;
  float sumnewNInc = 0;
  float sumnewSInc = 0;

  float NInc0 = Ntot(&(les->PopInc));
  float SInc0 = Qtot(&(les->PopInc));

  if (debug)
    show("controle devLes deb : Svert=%f, Slestot=%f(SInc=%f)\n",*Svert,SLesTot(les),SLesInc(les));


  //_________________Calcul Germination + heure germination , maj Udins
 
  ClassMyc Cgerm = calcGerm(les,dth,PPFD,Rh, debug);
  if (debug)
    show("\n Calc germ new germ N=%f (age=%f)\n",Cgerm.N,Cgerm.age);
  
  //_________________Calcul sous-iterations avec stop pour gemination et respect dTTmax

  int niter;
  int igerm;
  int ideb[MAXSUBSTEP];
  int ifin[MAXSUBSTEP];
  float dTT[MAXSUBSTEP];
  splitTimeStep(les,dth,&Cgerm,T,&niter,&igerm,ideb,ifin,dTT);
  if (debug)
    showIter(niter, igerm, ideb, ifin, dTT,les->dTTmax);

  //________________Boucle sous-iterations 
  
  int it,hdeb_it,hfin_it;
  float dTTit;

  for (it = 0; it < niter ; it++) {
    
    dTTit =  dTT[it];
    hdeb_it = ideb[it];
    hfin_it = ifin[it];
    
	if (eradicant_efficacy > 0) {	
		dTTit *= (1 - MAX(0, MIN(1, eradicant_efficacy)));	// a tester que les fonctions acceptent bien  toute un dTTit = 0
	}	
	
	
    if (debug)
      show("\n\n _____ sous iteration %d (hdeb: %d, hfin : %d, dTT: %.1f)\n",it,hdeb_it, hfin_it, dTTit);

    // ******* Vieillissement Spo (ne consomme pas de vert)
    AgeSpo(&(les->PopSpo),dTTit);
    

    // ******* 1.Vieillissement Chlo, Transition Chloroses -> Necroses (ne consomme pas de vert)

    if (les->skipChlo != 1) // Dchlo non null
      ChloToSpo(les,dTTit);

    if (debug)
      show("controle devLes fin 1 : Svert=%f, Slestot=%f(SInc=%f,Ninc=%f)\n",*Svert,SLesTot(les),SLesInc(les),Ntot(&(les->PopInc)));
    
    // ******* 2 Croissance des Chloroses    

    if (les->skipCroi !=1) // DCroi non null, PopCroi non vide
      CroiChlo(les, Svert,dTTit,debug);
    //enregistrement des recouvrements
    sumNIncRecChlo += NInc0 - Ntot(&(les->PopInc));
    sumSIncRecChlo += SInc0 - Qtot(&(les->PopInc));
    
    if (debug)
      show("controle devLes fin 2 : Svert=%f, Slestot=%f(SInc=%f,Ninc=%f)\n",*Svert,SLesTot(les),SLesInc(les),Ntot(&(les->PopInc)));

    // ******* 3 Vieillissement Incub, Transition Incub -> Chlorose/Croissance

    if (les->skipInc != 1) // DInc non null
      IncToChlo(les, Svert, hdeb_it, hfin_it, T, PPFD, Rh, eradicant_efficacy, &sumNIncRecChlo,&sumSIncRecChlo,debug);

    if (debug)
      show("controle devLes fin 3 : Svert=%f, Slestot=%f(SInc=%f,Ninc=%f)\n",*Svert,SLesTot(les),SLesInc(les),Ntot(&(les->PopInc)));

    // ******* 4 Croisssance des incubs  
 
    if (les->skipInc != 1) 
      CroiInc(les, Svert, dTTit, &sumNIncRecInc, &sumSIncRecInc, debug);
    
   
    // ******* 5 Nouvelles incubs  (TRansfo germinations -> incubs) 
  
    if (hfin_it == igerm) {
      NInc0 = Ntot(&(les->PopInc));
      SInc0 = Qtot(&(les->PopInc));
      GermToInc(les, &Cgerm, Svert, protectant_efficacy, debug);
      sumnewNInc += Ntot(&(les->PopInc)) - NInc0;
      sumnewSInc += Qtot(&(les->PopInc)) - SInc0;
    }

    // ********* gestion Svert <= epsillon
    if (*Svert < epsillon)
      *Svert = 0;
    if (debug)
      show("controle devLes fin sous iteration: Svert=%f, Slestot=%f(SInc=%f,NInc=%f)\n",*Svert,SLesTot(les),SLesInc(les),Ntot(&(les->PopInc)));
    
  }//fin sous iteration
  
  setLesRec(les, NIncRecChlo, sumNIncRecChlo);
  setLesRec(les, SIncRecChlo, sumSIncRecChlo);
  setLesRec(les, NIncRecInc, sumNIncRecInc);
  setLesRec(les, SIncRecInc, sumSIncRecInc);
  setLesRec(les, newNInc, sumnewNInc);
  setLesRec(les, newSInc, sumnewSInc);
  setLesRec(les, iSLesInSen, les->SLesInSen);
}

void devLes_d(Lesions *les,float *Svert,int dth,float *T,float *PPFD, float *Rh, int debug){
  devLes_dp(les,Svert,dth,T,PPFD,Rh,0, 0, debug);
}

void devLes(Lesions *les,float *Svert,int dth,float *T,float *PPFD, float *Rh){
  devLes_d(les,Svert,dth,T,PPFD,Rh,0);
}


//vidage des chlo en fin de senescence
void videChlo(Lesions *les,float *Svert) {
  ResetPop(&(les->PopCroi));
  ResetPop(&(les->PopMat));
  ClassMyc Cmyc = ClassEqQ(&(les->PopChlo));
  *Svert+=Cmyc.Q;
  ResetPop(&(les->PopChlo));
}
 
//fonction gerant l'evolution des surfaces des lesions (et modifications induites sur vert et surface senecente naturelle) lors du passage de la senescence naturelle. 

void senLes(float StoRec,Lesions *les,float *Svert,float *Sdsen) {
  if (StoRec > 0) {
    //calcul recouvrerements par type de lesion
    float rInc,rChlo,rNec,STot;
    STot = SLesOnGreen(les);
    rInc = SLesInc(les) / STot * StoRec;
    rChlo = SLesChlo(les) / STot * StoRec;
    rNec = (SLesNec(les) - les->SLesInSen) / STot * StoRec;
    // maj necrose : passage dans SLesInSen. Pas d'augmentation de Sdsen=> sene se propage plus vite que dans tem si il y a des lesions.
    les->SLesInSen += rNec;
    //chloroses  : on convertit en sene nat
    RemoveQtoPop(rChlo,&(les->PopChlo));
    *Sdsen += rChlo;
    if (SLesChlo(les) <= epsillon) {
      videChlo(les,Svert);
      rChlo += epsillon;
    }
    //enregistrement
    setLesRec(les,SChloRecSen,rChlo);
    //incubations
    float ninc0=Ntot(&(les->PopInc));
    float sinc0 = SLesInc(les);
    float Sdone = RemoveQbyN(rInc,&(les->PopInc),0);// on ne tue pas la derniere lesion
    *Sdsen += Sdone;
    
    //enregistrements
    setLesRec(les,NIncRecSen,ninc0-Ntot(&(les->PopInc)));
    setLesRec(les,SIncRecSen,sinc0 - SLesInc(les));
  }
}





//histo des dim de lesions 
void AfficheDimLes(Lesions *les,int nbClass) {
  float Slmin = les->par.Slmin;
  float Slmax = les->par.Slmax;
  float DChlo = les->par.DChlo;
  //affiche un histogramme des dimensions de lesion(vides et spo), 
  PopMyc DimLes = newPopMyc(Slmin, Slmax, (Slmax - Slmin) / nbClass,pTLes);
  ClassMyc Cmyc;
  float SChlo, NChlo;
  int i;
  //Calcul surface chlorotique des lesions en croissance =Chloroses tot-lesions 100%chlorotique (age<DChlo)
  SChlo = Qtot(&(les->PopChlo)) - NQtotEntre(&(les->PopCroi), 0., DChlo);
  //calcul nombre pondéré de lesion participant a surface chlorotique (poid 1 pour popCroi, poids degressif pour popMat)
  NChlo = NtotEntre(&(les->PopCroi), DChlo, bornemax(&(les->PopCroi)));
  for (i = 0; i < les->PopMat.nbClass; i++) {
    Cmyc = getCmyc(i,&(les->PopMat));
    if (DChlo > 0)
      NChlo += Cmyc.N * MAX(0,(1. - Cmyc.age / DChlo));
  }
  //remplissage de DimLes : Q=surface moyenne de lesions dans la barre,N= effectif
  if (NChlo>0 && SLesNec(les)>0) {
    //popcroi
    for (i=0;i<les->PopCroi.nbClass;i++) {
      Cmyc=getCmyc(i,&(les->PopCroi));
      if (Cmyc.age>=DChlo) {
	Cmyc.Q=Cmyc.Q-SChlo/NChlo;
	Cmyc.age=Cmyc.Q;
	AddCtoPop(Cmyc,&DimLes);
      }
    }
    //popmat
    for (i=0;i<les->PopMat.nbClass;i++) {
      Cmyc = getCmyc(i,&(les->PopMat));
      //show("Contenu de popmat:\n");PrintClass(Cmyc);
      if (DChlo > 0 && NChlo > 0)
	Cmyc.Q = MAX(0, Cmyc.Q - (1. - Cmyc.age / DChlo) * SChlo / NChlo);
      else
	Cmyc.Q = 0;
      Cmyc.age=Cmyc.Q;
      //show("ajout depuis popmat:\n");PrintClass(Cmyc);
      AddCtoPop(Cmyc,&DimLes);
    }
  }
  //ajout lesions mature
  AddCtoPop(newCmyc_p(Slmax-epsillon,Slmax,les->nbLesAd),&DimLes);
  
  show("\n************Stats lesions:\n");
  show("\nSChlo=%f,NChlo=%f,SLesNec=%f, Sles non affectée=%f\n",SChlo,NChlo,SLesNec(les),SLesNec(les)-NQtot(&DimLes));
  PrintPop(DimLes);
  freePop(&DimLes);
}	 


#endif 
