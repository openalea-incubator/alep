#ifndef _PARAMETER_H
#define _PARAMETER_H
//fichier parametres 14 auout 2007 (10 sects+couche 1cm)

//*********************************************parametre generaux et options

  //parametre de controle du 'bavardage' (0 = aucune info)
#define VERBOSE 0
//controlle tallage (commenter ou non la ligne)
//#define NOTILLER
//arret a last ini ? (0 non, 1 oui)
#define STOPLASTINI 1
//prefixe des fichiers de sortie
#define PREFIXOUT "Output"
// Active (1)/desactive (0) l'ecriture des fichiers disp
#define WRITE_OUTDISP_FILE 1
// Active (1) /desactive (0) la lecture des parametres d'analyse de sensibilite
#define SENSITIVITY 0

 //ecchelle pour la largeur du vumetre (nombre de spore/m2 pou avoir la largeur 1)
#define SCALESPO 1000
//option output toalea
//#define ALEASTR

//ordre de deroulement d'une iteration : maladie doit imperativement etre preced'interception pour avoir une gestion coherente des udin (age udin apres interception = age au debut du pas de temps (eventuellement negatif pour les udin du pas de temps), viellissement udin durant maladie)
#define PHYSIO 1
#define MALADIE 4
#define EMISSION 2
#define INTERCEPTION 3
#define TOALEA 5
#define INITIALISATION 9
#define MSGPHASE "Croissance et Senescence des Feuilles","Emission et redistribution des Spores","Interception des Spores","Developpement de la maladie","Ecriture Chaine Alea"
#define FIRSTPHASE 1

#ifdef ALEASTR
#define NBPHASES 5
#define LASTPHASE 5
#else
#define NBPHASES 4
#define LASTPHASE 4
#endif


//******************************************** Pas de temps nb iter et meteo

//pas de temps (en j) d'une iteration
#define TIMESTEP 1
//pas de temps mini 1h pour sous iterations maladies
#define MAXSUBSTEP (TIMESTEP * 24)
//duree de la simul (j)
#define STEPS 250
//Temperature (suposée constante)
#define TEMP 25
//max d'evenements de pluie par pas de temps
#define MAXEV 6


//Jour, heure debut simulation en Jsim
//1=1 octobre
#define JDEB 5
//heure en hhmm
#define HDEB 100

//repertoire meteo
//*********************************************************************
//ATTENTION si utilisation depuis L-Studio, il faut un chemin absolu !
//*******************************************************************
#define METEODIR "./Meteo"
  //fichiers meteo : format Année Jour hhmm Temp(°C) Humidité(%) T Vaisala(°C) PAR (micromol/m2/s) vent(m/s) Pluie(mm)
#define METEOFILE "meteo99-00.txt"
// Option climat different de emergence 1 a phtEndmeteodeb
#define PHTENDMETEODEB 7
#define METEOFILEDEB METEOFILE



  //*****************************************conditions initiales

  //on infecte uniquement les NBFEUINF premieres feuilles du brin maitre, a hauteur de PINFINI % de la surface. 
  //L'infection se fait en FRACINF fois, en infectant a chaque fois 1/FRACINF eme de la feuille. Ces infections se font succesivement de la pointe vers la base, DELINF dd apres la sortie de la zone de feuille concernee. 
  //Le dernier parametre permet de choisir le moment du declenchement de l'attente  pour les zones (0 signifie des la sortie du debut de la zone, 0.5 signifie au mment de la sortie du milieude la zone, et 1 apres la sortie de la fin de la zone)

  //nombre de feuilles (bm) subissant l'infection forcée
#define NBFEUINF 3
  //proportion de la surface des feuille à infecter
#define PINFINI 0.2
#define DLATINF DLAT
  //Fractionnement : Nbre d'infections pour atteindre pinf (= nbre de zones de la feuilles).
#define FRACINF NBSECT
#define NBZINF FRACINF
  //position relative dans la zone ([0-1]) du repere declenchant le delai d'attente lors de sa visibilite
#define PDELINF 0.5

// alternative
// duree (TT since emergence) durant laquelle le sol infecte
#define TTINFSOL 0
// LAI sporulant du sol (m2/m2)
#define LAISPOSOL 30e-4

//********************************************options du modèle de dispersion

  //Active (1)/desactive(0) la prise en compte des tiges dans la dispersion
#define SEND_STEM_TO_DISP 1
  //Active (1) / desactive (0) l'homogeneisation forcee des couches
#define HOMOGENEOUS_CANOPY 0

//nombre de direction pour le calcul de l'incident et du splash
#define NBH 15
//nombre de classe d'angle du lai (doit etre egal à nbh,cf dispersion.h)
#define NBCANG NBH
  //hauteur d'une couche (cm)
#define HC 1.


  // ---------------------parametres du splashing

//parametre d'agregation de la vegetation (Baret et al. 1993)
#define LAMBDA0 1.56

  // -- modelisation de la pluie

//hauteurs (angle/horizontale en degres) min (hmax=90,ie on echantillone de hmin a la verticale) pour les directions a echantillonner pou les calculs de coef de projection.On distingue pluie et splash
#define HMINPLUIE 85.
  //intensité minimale de pluie (mm/h) pour qu'il y est splash (Rapilly, 1976)
#define IMINPLUIE 0.5

// -- production goutelettes infectieuses
 
  //Proportionalite entre le Nombre total d'eclabousssure produite par une pluie et son intensite (mm/h) (Rapilly, 1976)
#define PFA 6.19e7
  //proportion de goutte non infectante en raison de l'evaporation (Rapilly, 1976)
#define PEV 0.4

  //-- dispersion

#define HMINSPLASH 10.
//nombre de points decrivant la proportion cumulée de gouttes splashées en fonction de la distance a la source 
#define NBPDROP 6
//Hauteur (cm) des cumuls
#define HDROP 0,5.,10.,15.,20.,22.5
//Proportion Cumulé de gouttes a hdrop 
#define CUMDROP 1.,0.55,0.24,0.09,0.02,0
// facteur multiplicatif des hauteurs
#define SCALEDROP 1.

  // -- depot et evolution udin
//durée de persistence de FracUdin (en j apres la contamination ayant généré un reste)
#define DPERFRAC 1

//******************************************** Parametres du modèle maladie

  // ------ Germination et penetration

  // taux de perte des Udin quand elles ne germent pas (fraction de l'effectif initial perdu par jour)
#define TXPERTEUDIN 0.2
  // % Humidité air mini requis pour avancer le processus de germination (= permettant  98 % humidite sur les feuilles)
#define HRGERM 85.
  // Rayonnement incident (micromolPAR/m2/s) maxi pour avancer le processus de germination ( ie permettant 98 % humidite sur les feuilles)
  //140 wm2 definit ensoleimment * 0.2174 : conversion vers micromol/m2/s-1
#define PARGERM (140./0.2174)
  // Dure minimale continu à HRGERM pour obtenir une germination (heures)
#define DMINGERM 10.
  // Proba de reussite (penetration + incubation) pour les UDIN
#define PINC 1.

  // ------- Latence (incub + chlorose)

  // intervalle degree-jour definissant une cohorte de meme age durant la latence (largeur classe d'age)
#define LARGCLASS 10.

  // Parametrisation des reponses a la temperature des durees de latence et des processus de croissance des lesions
  //T(°C) de base du mycelium 
#define TBASEMYC 0
  //T(°C) maximale (pas d'accumulation de temps thermique au dela)
#define TMAXMYC 25

  // Parametrisation de la latence 

// Active (0) / Desactive (1) la parametrisation avec prise en compte de la reponse de la duree d'incubation aux conditions d'humdite (Rapilly 74)
#define FIXED_LAT 1

// Parametrisation si FIXED_LAT = 1
// Duree (dd) d'incubation et de latence
#define DINC 220. 
#define DLAT 330.

// Parametrisation si FIXED_LAT = 0
  // Duree (dd) de l'incubation en condition humide (parametre del2 de Rapilly), pris en compte si FIXED_LAT = 0
#define DINCMIN 155.
  // duree max de la phase d'incubation (dd)
#define DINCMAX 500.
// duree de la croissance des incub (de MININC a MINLES)
#define DCROIINC 220.
  // facteur de correction de la duree d'incub par les condition d'humidite
#define FDHR (-0.121)
  // Duree de la phase chlorotique (dd)
#define DCHLO 110.


// ------ Croissance des lésions

//surface min affectée au inc apres penetration
#define MININC 0.001
  //Surface mini d'une chlorose en sortie d'incubation (cm2)
#define MINLES .03
  //surface max d'une lesion issue d'une spore (cm2)
#define MAXLES .3
  //Taux croissance des lésions chlorotiques en cm2/dd
#define TXCROI 0.0006



  // ------ Vidage des lesions par les pluies

  //nombre d'événements splashants pour arriver au vidage d'un pycnide
#define NBEVSPLAMAX 3
  //proportion de gouttes non infectantes produite par unite de surface sporulante (Rapilly, 1976)
// (en raison d'un contenu de moins de 5 spores)  
//capacite de production d'udin = 1 - PnoSpo = proportion d'udin parmis les goutes emises 
#define PNOSPO 0.24,0.24,0.24


// ----- Estimation nombre de spores (deprecated)

  // Fractionnement du vidage entre evenements splashants 
//PSPOEV 0.5,0.2,0.1,0.1,0.1
//#define PSPOEV 0.6,0.275,0.12
#define PSPOEV PNOSPO
 //densité stomatique des feuilles (/cm2)
#define DSTOM 5000.
  //Proportion de stomates infectés par des pycnides dans une lesion
#define PSTOINF 0.5
  //Nbre total de spores produits par une pycnide (total de ce qui sera splashé par l'ensemble des evenement de pluie jusqu'au vidage de la lésion + spores residuelles dans la lesion "vidée", cf Eyal)
#define NBTOTSPORE 4000
  //prportion de spores résiduelles par pycnide (Eyal)
#define PSPORERES .2
 //Tx de production de spores (spores/jEqpluie/cm2 lesion)
#define TXSPO 2000


 
  //*********************************************Parametre du modèle plante

  // ---- Nb secteurs et representation

//discretisation des feuilles : nbre de secteurs
#define NBSECT 5
//parametre de forme de feuille
#define ALPHA -2.3
//nb pol min par feuille
#define NBPOLMIN 5

  // ---- Shoot geometry

// Active (1), Desactive(0) l'usage de phi0 dans dimT pour le controle de la forme des limbes
#define USE_PHI0 0

// nombre de courbe deffinissant l'evolution de la forme
#define NBLEAFSTAGES 10
// nombre de points par courbe
#define NBPLST 50
//Tableaux s,x,y definissant les shapes
#define LEAFSHAPES 0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0,0.01,0.01,0.02,0.02,0.03,0.03,0.03,0.04,0.04,0.05,0.05,0.06,0.06,0.07,0.07,0.08,0.08,0.09,0.09,0.1,0.1,0.11,0.11,0.12,0.12,0.13,0.13,0.13,0.14,0.14,0.15,0.15,0.16,0.16,0.17,0.17,0.18,0.18,0.19,0.19,0.2,0.2,0.21,0.21,0.22,0.22,0.23,0.23,0.23,-0.01,0.01,0.03,0.05,0.07,0.09,0.11,0.13,0.15,0.17,0.19,0.21,0.23,0.25,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0,0.01,0.01,0.02,0.02,0.03,0.03,0.04,0.04,0.05,0.05,0.06,0.06,0.07,0.07,0.08,0.09,0.09,0.1,0.1,0.11,0.12,0.12,0.13,0.14,0.14,0.15,0.16,0.16,0.17,0.18,0.18,0.19,0.2,0.21,0.21,0.22,0.23,0.24,0.25,0.26,0.26,0.27,0.28,0.29,0.3,0.31,0.31,0.32,0.33,-0.01,0.01,0.03,0.05,0.07,0.09,0.11,0.13,0.15,0.17,0.19,0.21,0.23,0.25,0.27,0.29,0.31,0.33,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.75,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.89,0.91,0.93,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0,0.01,0.01,0.02,0.03,0.03,0.04,0.04,0.05,0.05,0.06,0.07,0.07,0.08,0.09,0.09,0.1,0.1,0.11,0.12,0.12,0.13,0.14,0.14,0.15,0.16,0.17,0.17,0.18,0.19,0.2,0.21,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,-0.01,0.01,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.64,0.66,0.68,0.7,0.72,0.73,0.75,0.77,0.79,0.8,0.82,0.84,0.86,0.87,0.89,0.91,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0,0.01,0.01,0.02,0.03,0.04,0.05,0.06,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0,0.01,0.03,0.05,0.07,0.09,0.11,0.12,0.14,0.16,0.18,0.2,0.22,0.23,0.25,0.27,0.29,0.31,0.33,0.34,0.36,0.38,0.4,0.42,0.43,0.45,0.47,0.49,0.51,0.52,0.54,0.56,0.58,0.6,0.61,0.63,0.65,0.67,0.68,0.7,0.72,0.74,0.76,0.77,0.79,0.81,0.83,0.84,0.86,0.88,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,-0.01,0.01,0.02,0.03,0.05,0.06,0.07,0.09,0.1,0.11,0.13,0.14,0.15,0.17,0.18,0.19,0.21,0.22,0.23,0.25,0.26,0.27,0.29,0.3,0.31,0.33,0.34,0.35,0.37,0.38,0.39,0.41,0.42,0.43,0.45,0.46,0.47,0.49,0.5,0.51,0.53,0.54,0.55,0.57,0.58,0.59,0.61,0.62,0.63,0.65,-0.01,0.01,0.03,0.04,0.06,0.08,0.1,0.11,0.13,0.15,0.17,0.19,0.2,0.22,0.24,0.26,0.27,0.29,0.31,0.32,0.34,0.36,0.37,0.39,0.41,0.42,0.44,0.45,0.47,0.48,0.5,0.51,0.53,0.54,0.56,0.57,0.58,0.6,0.61,0.62,0.63,0.65,0.66,0.67,0.68,0.7,0.71,0.72,0.73,0.75,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,-0.02,0,0.01,0.03,0.05,0.06,0.08,0.1,0.11,0.13,0.15,0.16,0.18,0.2,0.21,0.23,0.25,0.26,0.28,0.3,0.31,0.33,0.35,0.36,0.38,0.4,0.41,0.43,0.45,0.46,0.48,0.5,0.52,0.53,0.55,0.57,0.58,0.6,0.62,0.64,0.65,0.67,0.69,0.71,0.72,0.74,0.76,0.78,0.79,0.81,0.04,0.05,0.06,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.57,0.58,0.59,0.6,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0,0.01,0.03,0.04,0.06,0.08,0.09,0.11,0.12,0.14,0.16,0.17,0.19,0.21,0.22,0.24,0.25,0.27,0.29,0.3,0.32,0.34,0.36,0.37,0.39,0.41,0.42,0.44,0.46,0.48,0.49,0.51,0.53,0.55,0.57,0.58,0.6,0.62,0.64,0.66,0.68,0.69,0.71,0.73,0.75,0.77,0.79,0.81,0.82,0.84,0.01,0.03,0.04,0.06,0.08,0.1,0.12,0.13,0.15,0.17,0.18,0.2,0.22,0.23,0.24,0.26,0.27,0.28,0.29,0.3,0.3,0.31,0.31,0.32,0.32,0.32,0.32,0.32,0.31,0.31,0.3,0.3,0.29,0.29,0.28,0.27,0.27,0.26,0.25,0.24,0.23,0.22,0.21,0.21,0.2,0.19,0.18,0.17,0.16,0.15,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0.02,0.03,0.05,0.07,0.08,0.1,0.12,0.13,0.15,0.17,0.19,0.2,0.22,0.24,0.25,0.27,0.29,0.3,0.32,0.34,0.35,0.37,0.39,0.41,0.42,0.44,0.46,0.47,0.49,0.51,0.52,0.54,0.56,0.58,0.59,0.61,0.63,0.64,0.66,0.68,0.69,0.71,0.73,0.74,0.76,0.78,0.8,0.81,0.83,0.85,0.04,0.05,0.06,0.08,0.09,0.1,0.12,0.13,0.14,0.15,0.17,0.18,0.19,0.2,0.2,0.21,0.22,0.22,0.23,0.23,0.24,0.24,0.24,0.24,0.23,0.23,0.23,0.22,0.22,0.21,0.2,0.19,0.18,0.17,0.16,0.14,0.13,0.11,0.1,0.08,0.07,0.05,0.04,0.02,0,-0.01,-0.03,-0.04,-0.06,-0.08,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,0,0,0.01,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.19,0.21,0.22,0.24,0.25,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.38,0.39,0.4,0.41,0.41,0.41,0.42,0.42,0.42,0.42,0.43,0.43,0.43,0.44,0.44,0.43,0.43,0.43,0.42,0.41,0.4,0.01,0.02,0.04,0.05,0.06,0.07,0.07,0.08,0.07,0.07,0.06,0.05,0.04,0.03,0.01,0,-0.02,-0.04,-0.05,-0.07,-0.09,-0.1,-0.12,-0.14,-0.16,-0.18,-0.2,-0.21,-0.23,-0.25,-0.27,-0.29,-0.31,-0.33,-0.35,-0.37,-0.39,-0.41,-0.43,-0.45,-0.47,-0.49,-0.51,-0.53,-0.55,-0.57,-0.59,-0.61,-0.63,-0.65


  // ---- Shoot dynamic

//  1 is for old style parameterisation of axeT and phenT (see at the end of file)
#define OLDSTYLE 1
// Active (1) / Desactive (0) le blocage de la senescence
#define NOSEN 0

  //Densité de semis (pl/m2)
#define DENSITE 250

  //nombre d'axes simulés (bm ou talles avec caracteristiques differentes)
#define NBAXES (NBTILW + 1)
  // axe Table (empty if OLDSTYLE=1)
// NbF : Nombre de phytomere sur les axes 
// Abond : Abondances des axes dans la pop au moment de l'emission  (axe density = freq * densite)
// hw0 : longueur du cornet au moment de l'emergence de l'axe
// Ploss : proportion of axes that will die
// Axeid : identifiant axe (pour la sortie);
// type : type de l'axe (0 = bm, ou > 0 tiller), pour le calcul de stats bm
//StartReg/End Reg : Thermal time at beginning and end of senescence of axes

#define AXE_T

  // ----- Growth of axes (ADEL)

  // setting the option to 0 uses linear thermal time. 1 means use compensated thermal time
#define USE_COMPENSATED_THERMAL_TIME 0
  //T(°C) de base/max de la plante (si linear TT)
#define TBASEPLANTE 0
#define TMAXPLANTE 50
  //parametre temps thermique compense
#define TCREF_TC 12
#define K_TC 3.8088e10
#define EAR_TC 8899.6
#define DSR_TC 68
#define DHR_TC 20736

  // relative phyllochronic time (Phyll-numf) at start of rapid leaf extension
#define RPHLEAFSTART -.4
//relative phyllochronic time at end of leaf extension
#define RPHLEAFEND 1.6

// Nombre de points pour le tableau des phyllochrone (nbphy + "0")
#define NBPHEN (NBPHY + 1)
// Phen Table  (empty if OLDSTYLE =1)
// Tip = date emergence Tip
// Sen = date full senescence
// Disp = date disparition
#define PHEN_T

// dimensions

  //hauteur max des tiges (col de la feuille la plus haute) en cm
#define HMAX 100.
  //longueur du plus grand limvbe
#define LBMAX 30.
  //nombre de Phyto le plus grand
#define NBDIM NBPHY
// tableau dimT (empty si old sytle)
// L : Longueur de limbe
// G : Longueur Gaine
// w : largeur limbe
// E : longeur entrenoeud
// diam : diametre tige
// Phi0 : angle a la base
// StemRate : coef modificateur vitesse tige
// Senpattern : indice du pattern de senescence (tableaux xrssi,yrssi)
#define DIM_T

  // syncho sene individuelle des feuilles avec ssi (nouvelles valeurs par defaut depuis juillet 2010 = compil Miniere/picardie/Rim)

#define NNRSSI (NBFDELSEN + 2)
#define NRSSI 3
#define XRSSI {-3,-1.1,0},{-2,-1,0.9},{-2,-1,0.7},{-1.7,-.9,0.4},{-1.5,-1,0.2},{-1,-1,0}
#define YRSSI {0,.25,1},{0,.25,1},{0,.15,1},{0,.12,1},{0,.12,1},{0,0,1}

// Old defaut si YRSSI pas rempli :
//Nombre de feuille avec rssi special (defaut = lineaire de -1 a 0)
#define NBFDELSEN 4
#define T1DELSEN -1.,-.8,-.5,-1.
#define SENRATE1 0.07
#define T2DELSEN 0.2,0.5,0.7,0.


// Parametres specifiques pour SimCycle

#define SC_TIMESTEP 1
#define SC_STEPS 150
#define SC_TEMP 10.
#define SC_SVERT 10.
#define SC_UDINSTART 1.
#define SC_AUTOINF 0.
#define SC_VERBOSE 0
#define SC_ENDSEN 1500
//duree de la sene rapide pour la flag leaf ~ 200 dd (avec parametrisation srtandard du modele)
#define SC_DURSEN (DTTENDSEN / NFVFLO)
//flag pour activer le modele de pluie (0 = pas de pluie, 1 = avec pluie)
#define SC_MAKEPLUIE 0
//intensite de la pluie (mm/h, 1.8 = moyenne grignon 98-99)
#define SC_IPLUIE 1.8
//dd premiere pluie
#define SC_FIRSTPLUIE 300.
//frequence des pluies : intervalle dd entre pluies
#define SC_FREQPLUIE 150.
//nombre de pluies
#define SC_NPLUIE 5
//heure de la pluie dans le pas de temps
#define SC_HPLUIE 12
//options pour randomiser les simulations (1 = avec randomisation, 0 = sans randomisation (=> avec repetabilite)
#define SC_RANDOMIZE 0


// ---- Old-style parametrisation (deprecated although compatible)

// Shoot dynamic

  //nombre de vagues de tallage (potentiel=6)
#define NBTILW 3
  //nombre de feuille sur les axes
#define NFVALUES 11,8,7,7
  //temps phyllochronique entre emergence f1 bm et emergence F1 de la ieme vague (valeurs types poue les vagues 1 a 6 : 3.05,3.6,4.6,6,7,8
#define DECPHYL 3.05,3.6,4.6
//nombre moyen d'axe (bm+talles) par vague et par plante(potentiel généralement observé= 1,1,1,2,4,3,3) lors de l'emission. Cumul = nbmax talles presentes sur une plante
#define NBAXEDEB 1.,1.,1.,1.2
//Nombre moyen d'axe par plante et par vague a la recolte
#define NBAXEFIN .95,.8,.8,.16
//delai (phyllochronic time) entre derniere vague d'emission et debut disparition des talles
#define DELREGT 2.
//durée (phyllochronic time) entre debut et fin de la regression des talles
#define DURREGT 3.

// development

#define PHYL 110.0
#define PHYLINI PHYL
#define TTENDINI ((LASTINI - 1) * PHYLINI)
#define TTFLAG TTENDINI + (NBPHY - LASTINI) * PHYL

// SSI
#define NFVVEG 3.5
#define NFVFLO 2.7
#define NFVEND 0.
#define DPHYLFLO 5.2
#define DTTENDSEN 550.

#define NBPHT 4
#define XPHT -PHYLINI, 0, TTENDINI, TTFLAG
#define YPHT 0, 1, LASTINI, NBPHY


/* Shoot senescence index */

#define NBSSI 4
#define XSSI (NFVVEG + 1.6 - 1) * PHYL,(NBPHY + 1.6 - 1) * PHYL,(NBPHY + 1.6 - 1 + DPHYLFLO) * PHYL,(NBPHY + 1.6 - 1 + DPHYLFLO) * PHYL + DTTENDSEN
#define YSSI 0,NBPHY-NFVVEG,NBPHY-NFVFLO,NBPHY


 // ----- Plant Morphology


//nb phyto
#define NBPHY 11
  //variables par num de phy (L,G,lw,E,Phi_0,stemrate factor) et axe 
#define PHYTOVALUES {8.125,6.35,7.825,8.825,9.25,9.5,11.65,12,9.35,13.1,15.175,16.8125,10,17.125,19.85,21.45,11.4,21.45,24.5,26.05,13.7,25,29.05,29,16.55,29.2,26.6,24.15,19.8,24.75,-999,-999,25.175,-999,-999,-999,28.8,-999,-999,-999,24.1,-999,-999,-999},{0.3,0.2625,0.3,0.35,0.325,0.4,0.465,0.5,0.4,0.7,0.715,0.675,0.45,0.825,1,0.89,0.55,1.1,1.15,0.935,0.75,1.2,1.25,1.1875,1,1.35,1.65,1.5,1.2,1.65,-999,-999,1.28,-999,-999,-999,1.425,-999,-999,-999,1.8,-999,-999,-999},{3,3,3.5,3.9,3.05,4.025,5.25,6.1,3.05,6.5,9,9.5,3.4,10.2,11.2,11.6,4.2,11.5,13,13.025,6.225,13.85,15.95,14.45,9.125,16.6,17.3,16.6,12,18.0875,-999,-999,14.2,-999,-999,-999,17.2,-999,-999,-999,18.675,-999,-999,-999},{0,0,0,0,0,0,0.1,0.1,0,0.1,1.55,1,0,2.3,5.6,5.55,0,6.45,8.95,9.8,0.1,10.3,12,10.475,1.9,12.75,14.9,15.2,6.1,16.05,30.2,26.3,9.675,31.55,-999,-999,14.45,-999,-999,-999,16.95,-999,-999,-999},{25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25},{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}

// Facteurs de modification (pour analyse sensibilité)

//rang de la derniere feuille ne subissant pas les modifs scale_xx et phylini
#define LASTINI 1
//facteur d'echelle pour la longueur des feuilles
#define SCALE_BL 1.0
//facteur multiplicatif de l'angle des feuilles
#define SCALE_PHI0 1.0
  //facteur d'echelle pour les dimensions d'entrenoeud
#define SCALE_EN 1.0
//facteur multiplicatif de la vitesse d'allongement des entrenoeuds
#define SCALE_STEMRATE 1.0
//facteur multiplicatif de la dimension d'entrenoeud, et divisant la vitesse (ie plus grand a vitess egale)
#define SCALE_ENNOTRATE 1.0
//facteur d'echelle longueur gaines
#define SCALE_SH 1.0
//facteur d'echelle pour les largeurs de feuilles
#define SCALE_LW 1.0
//facteur d'echelle pour les surfaces de feuilles (ettenue de facon homogene longueur et largeur)
#define SCALE_SURF 1.0
//facteur d'echelle tige (en et gaine)
#define SCALE_ENSH 1.0
//facteur d'echelle pour manipulation surface derniere et avant derniere feuille a LAI constant (avnt derniere feuille compense)
#define SCALE_LAST 1.
//lambda
#define SCALE_LAMBDA 1.



#endif


