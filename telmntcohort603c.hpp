/* *************************************************************
TELMNTCOHORT603C.HPP - Container class to hold land cover cohort 
                      characteristics

20051119 - Created by DWK
20070514 - DWK added int standage 
20070830 - DWK added double wfpsoff
20070917 - DWK added double nonsolc and nonsoln
 
************************************************************* */

#ifndef ELMNTCOHORT603C_H
#define ELMNTCOHORT603C_H

#include "temconsts602.hpp"
#include "bioms423.hpp"

struct ElmntCohort60
{      
  // Index for agricultural crop type used for
  //   parameterization
  int agcmnt;

  // Number of growing degree days in croplands
  double aggrowdd;

  // Current value of decomposition parameter in croplands 
  double agkd;

  // Index to indicate whether cohort was in
  //   agriculture during the previous year
//  int agprevstate;

  int agprvstate;

  // Index to indicate whether cohort is in
  //   agriculture during the current year
  int agstate;

  double c2n;

  // Area covered by cohort (square kilometers)
  //   during current year 
  double chrtarea;

  // Index for vegetation community type used for
  //   parameterization
  int cmnt;

  double cneven;

  Biomass convrtflx;

  double cropprveetmx;
  double cropprvleafmx;
  double cropprvpetmx;

  Biomass cropResidue;

  double croptopt;

  // Index for current vegetation type
  int currentveg;

  // Counter of months since disturbance
  int distmnthcnt;

  // Flag to indicate that a disturbance has occurred
  //   in a particular year.  Different values represent
  //   different types of disturbances:
  //   0 = no disturbance
  //   1 = sustained disturbed state 
  //         (e.g. conversion to agriculture or urban state)
  //   2 = timber harvest
  //   3 = fire
  //   4 = insect & disease
  int disturbflag;
     
  // Month in which disturbance occurred
  // (e.g. 1 = January, ...., 12 = December)
  int disturbmonth;

  // Monthly soil temperature at 10 cm
  double dst10[CYCLE];

  // Maximum monthly estimated evapotranspiration
  //   during the current year
  double eetmx;

  // Index to indicate if crops are fertilized (= 1)
  //   or not (= 0)
  int fertflag;

  // Counter of months since fire disturbance
  int firemnthcnt;

  // Monthly reintroduction of N volatilized in fires
  //   back to the soil (as N deposition)
  double firendep;

  Biomass formPROD10;
  Biomass formPROD100;

  double fprevozone;

  // Fire return frequency
  double FRI;

  Biomass initPROD1[CYCLE];
  Biomass initPROD10[10];
  Biomass initPROD100[100];

 // Index to indicate if crops are irrigated (= 1)
  //   or not (= 0)
  int irrgflag;

  // Current value of decomposition parameter in 
  //   natural ecosystems
  double kd;

  double natprveetmx;
  double natprvleafmx;
  double natprvpetmx;
  double natseedC;
  double natseedSTRN;
  double natseedSTON;
  double natsoil;
  double nattopt;
  double natyreet;
  double natyrpet;

  double newleafmx;

  double newtopt;

  // "Non-reactive soil organic carbon
  double nonsolc;

  // "Non-reactive soil organic nitrogen
  double nonsoln;

  double nretent;
  double nsretent;
  double nvretent;

  // Maximum monthly potential evapotranspiration 
  //   during current year
  double petmx;

  int potveg;

  // Air temperature that occurred 2 months ago
  double prev2tair;
  
  // Atmospheric CO2 concentration during the previous month
  double prevco2;

  Biomass prevCropResidue;

  // Soil temperature at 10cm during the previous month
  double prevdst10;

  Biomass prevPROD1;
  Biomass prevPROD10;
  Biomass prevPROD100;

  // Snowpack during the previous month
  double prevspack;

  // Air temperature of the previous month
  double prevtair;
  // rain of the previous month
  double prevrain;
  double prevmdc;

  // Unnormalized relative leaf area of the previous month
  double prevunrmleaf;

  // Value of y[] during previous time step
  double prevy[MAXSTATE];

  //double prod10par;

  //double prod100par;

  int productYear;

  // Area covered by cohort (square kilometers)
  //   during previous year 
  long prvchrtarea;

  // Crop net primary production during the previous month
  double prvcropnpp;

  // Maximum monthly estimated evapotranspiration
  //   during the previous year
  double prveetmx;

  // Maximum relative leaf area of cohort
  //   during the previous year 
  double prvleafmx;

  // Maximum monthly potential evapotranspiration
  //   during the previous year
  double prvpetmx;

  int qc;

  //double sconvert;

  Biomass sconvrtflx;

  Biomass slash;
  Biomass deadwood;

  double slashpar;

  // Source cohort for current cohort
  int srcCohort;

  // Age of cohort
  int standage;

  //forest actual age after disturbance
  int forestage;

  double STMdx9[MAXNODES];

  int STMis9;

  double STMsmass9;

  double STMt9[MAXNODES];

  double STMwater9[MAXNODES];

  double STMweight9[MAXSNODES];

  double STMx9[MAXNODES];

  double STMxfa9[MAXNODES];

  double STMxfb9[MAXNODES];

  int subtype;

  int tillflag;

  double topt;

  int tqc;

  double vconvert;

  Biomass vconvrtflx;

  double vrespar;
  
  double wfpsoff;

  double y[MAXSTATE];

  double yrltrc;

  double yrltrn;
  double leafmortpar;
  double stemmortpar;
  double rootmortpar;
  double leafslash;
  double leafconv;
  double rootslash;
  double rootconv;
  double stemslash;
  double stemconv;
  double standdead;
  double deadconv;
  double deadslash;
  double sconvert;
  double prod10par;
  double prod100par;
};

#endif
