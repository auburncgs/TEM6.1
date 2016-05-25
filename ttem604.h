/* *************************************************************
****************************************************************
TTEM604.H - Terrestrial Ecosystem Model Version 6.0
****************************************************************

Modifications:

20030712 - DWK created by modifying ttem50b1.h to incorporate
           open-N cycle dynamics into TEM
20030712 - DWK changed include from atms43.h to atms51.h
20030712 - DWK changed include from tveg50b1.h to tveg51.h
20030712 - DWK changed include from tsoil50b1.h to tsoil51.h
20030712 - DWK changed include from tmcrb50b1.h to tmcrb51.h
20030712 - DWK changed include from humnact50b1.h to humnact51.h
20030712 - DWK changed TTEM50 to TTEM51
20030712 - DWK changed I_MNUP to I_NIMM in temkey
20030712 - DWK added I_TNDEP, I_SNFIX, I_ANFIX, I_GMIN, I_AMMN,
           I_NOP, I_N2OP, I_NO3P, I_NTRF, I_N2P, I_DNTRF,
           I_NH3FLX, I_NOFLX, I_N2OFLX, I_N2FLX, I_LCH,
           I_AGSNFX and I_TOTSNFX to temkey
20030712 - DWK changed Atmosphere43 atms to Atmosphere51 atms
20030712 - DWK changed Tveg50 veg to Tveg51 veg
20030712 - DWK changed Tsoil50 soil to Tsoil51 soil
20030712 - DWK changed Tmicrobe50 microbe to Tmicrobe51 microbe
20030712 - DWK changed Humanact50 ag to Humanact51 ag
20030712 - DWK changed include ttem50b1.cpp to ttem51.cpp at bottom
           of file
20030718 - DWK added I_NH4, I_NO3, I_NH4DEP, I_NO3DEP, I_INNH4UP,
           I_INNO3UP, I_VNH4UP, and I_VNO3UP to temkey
20030718 - DWK changed const int MAXESTAT = 20 to
           const int MAXESTAT = 22
20030718 - DWK changed const int NUMTEM = NUMEQ + 26 to
           NUMTEM = NUMEQ + 28
20030719 - DWK added I_LCHDOC, I_RSOIL and I_LCHDON to temkey
20030719 - DWK changed I_LCH to I_LCHNO3 in temkey
20030719 - DWK changed NUMEEQ from "MAXESTAT+76" to "MAXESTAT+101"
20030719 - DWK added public double rsoil[CYCLE] and double
           yrrsoil
20030719 - DWK added I_NATVEGN, I_NATSTRN, I_NATSTON, I_NATULF,
           I_NATLEAF, I_NATLAI, I_NATFPC,  I_NATINGPP, I_NATGPP,
           I_NATINNPP, I_NATGPR, I_NATRVMNT, I_NATRVGRW,
           I_NATINNUP, I_AINNH4UP, I_NINNH4UP, I_AINNO3UP,
           I_NINNO3UP, I_NATVNUP, I_AGVNH4UP, I_NVNH4UP,
           I_AGVNO3UP, I_NVNO3UP, I_NATVSUP, I_NATVLUP,
           I_NATVNMBL, I_NVNRSB and I_NATLTRN to enum temkey
20030719 - DWK changed I_TOTNPP to I_NATNPP in temkey
20030719 - DWK changed TOTVEGC to NATVEGC in temkey
20030719 - DWK changed const int NUMTEM = NUMEQ + 28 to
           NUMTEM = NUMEQ + 56
20030731 - DWK added I_ABIMMB to enum temkey
20030731 - DWK changed  NUMEEQ = MAXESTAT + 101 to
           NUMEEQ = MAXESTAT + 102
20030819 - DWK added public int climitflg
20030826 - DWK added I_DOC, I_DON, I_DOCP and I_DONP to
           enum temkey
20030826 - DWK changed MAXESTAT = 22 to MAXESTAT = 24
20030826 - DWK changed NUMEEQ = MAXESTAT + 102 to
           NUMEEQ = MAXESTAT + 104
20030902 - DWK added I_NH4IMM, and I_NO3IMM to
           enum temkey
20030902 - DWK changed NUMEEQ = MAXESTAT + 104 to
           NUMEEQ = MAXESTAT + 106
20030909 - DWK changed I_ABIMMB to I_AIMMNH4 and added I_AIMMNO3
           to enum temkey
20030909 - DWK changed NUMEEQ = MAXESTAT + 106 to
           NUMEEQ = MAXESTAT + 107
20031015 - DWK added I_CO2G, I_CO2W, I_HCO3, I_CO3, and I_ALK to
           enum tem
20031015 - DWK changed MAXESTATE = 24 to MAXESTATE = 29
20031015 - DWK added I_ERDPOC, I_DSVCO2, I_LCHCO2, I_HCO3P,
           I_LCHHCO3, I_CO3P, I_LCHCO3, I_LCHALK and
           I_NTCB, I_ISOPREN, I_TERPEN, I_ORVOC, I_OVOC,
           I_VOC, I_ERDPON, I_ABVGPR and I_ROOTGPR to enum tem
20031015 - DWK changed NUMEEQ = MAXESTAT + 107 to
           NUMEEQ = MAXESTAT + 124
20031015 - DWK added double ntcs[CYCLE] and double yrntcs
20031019 - DWK changed char version[20], char sitename[80],
           char developer[80] and char description[80] to
           string version, string sitename, string developer
           and string description, respectively
20040201 - DWK changed include from ttem51.cpp to ttem511.cpp
           at bottom of file
20040203 - DWK changed include from tmcrb51.h to tmcrb512.h
20040203 - DWK changed include from ttem511.cpp to ttem512.cpp
           at bottom of file
20040220 - DWK changed include from ttem512.cpp to ttem513.cpp
           at bottom of file
20040228 - DWK changed class TTEM51 to TTEM60
20040228 - DWK changed include from tmcrb512.h to tmcrb60.h
20040228 - DWK changed Tmicrobe51 microbe to Tmicrobe60 microbe
20040228 - DWK changed include from ttem513.cpp to ttem60.cpp
           at bottom of file
20040707 - DWK changed include from tsoil51.h to tsoil60.h
20040707 - DWK changed Tsoil51 soil to Tsoil60 soil
20040707 - DWK changed include from tveg51.h to tveg60.h
20040707 - DWK changed Tveg51 veg to Tveg60 veg
20040707 - DWK changed include from humnact51.h to humnact60.h
20040707 - DWK changed Humanact51 ag to Humanact60 ag
20040707 - DWK changed char predstr[NUMTEM][9] to 
           string predstr[NUMTEM]
20040716 - DWK changed include from atms51.h to atms60.h
20040716 - DWK changed public Atmosphere51 atms to
           Atmosphere60 atms
20040922 - DWK added public parameters double no3cut[MAXCMNT],
           double no31a[MAXCMNT], double no31b[MAXCMNT],
           double no32a[MAXCMNT], and double no32b[MAXCMNT]
20040922 - DWK deleted public parameters double no3a[MAXCMNT]
           and double no3b[MAXCMNT]
20041003 - DWK added I_NTNS to enum temkey
20041003 - DWK changed global const int NUMEEQ = MAXESTAT + 124 
           to const int NUMEEQ = MAXESTAT + 125
20041003 - DWK added public double ntns[CYCLE] and 
           double yrntns to represent net terrestrial nitrogen
           storage
20050408 - DWK changed include from ttem60.cpp to ttem601.cpp
20050409 - DWK changed global const NUMEEQ from 
           NUMEEQ = MAXESTAT + 125 to NUMEEQ = MAXESTAT + 127                                                          
20050409 - DWK added I_CDCMP and INDCMP to enum temkey
20050409 - DWK changed include from tmcrb60.h to tmcrb601.h
20051117 - DWK added include temconsts602.hpp
20051117 - DWK changed include from atms60.h to atms602.h
20051117 - DWK changed include from tveg60.h to tveg602.h
20051117 - DWK changed include from tsoil60.h to tsoil602.h
20051117 - DWK changed include from tmcrb601.h to tmcrb602.h
20051117 - DWK changed include from humnact60.h to humnact602.h
20051117 - DWK changed include from odeint425.h to odeint602.h
20051117 - DWK changed Odeint4 to Odeint60 in class declaration
20051117 - DWK deleted ttem601.cpp from bottom of file
20051118 - DWK added functions and variables from class Odeint60
20051118 - DWK deleted include odeint60.h
20051118 - DWK deleted inheritance of class Odeint60
20051118 - DWK changed public virtual int boundcon() to
           public int bouncon()
20051118 - DWK changed public virtual int delta() to
           public void delta()
20051123 - DWK deleted public equilibrium() and stepyr()
20051123 - DWK changed protected int initFlag to 
           public int initFlag
20051123 - DWK deleted public function setMonth()
20051123 - DWK added public function updateYearSummary()
20051123 - DWK deleted public function resetODEflux()
20051123 - DWK added public functions getY() and setY()
20051124 - DWK changed double nirr[CYCLE] and double tair[CYCLE] 
           to const double& nirr and const double& tair in
           function call to setEquilEvap() 
20051124 - DWK added const double& co2 in function call to
           setEquilC2N()
20051124 - DWK deleted const int& pdm from function call of 
           deltaxclm(), massbal() and monthxclm()
20051126 - DWK changed public double bnfix[CYCLE] to 
           double bnfix
20051126 - DWK changed public double nep[CYCLE] to double nep
20051126 - DWK changed public double ntcs[CYCLE] to double ntcs 
20051126 - DWK changed public double ntns[CYCLE] to double ntns
20051126 - DWK changed public double rsoil[CYCLE] to 
           double rsoil
20051126 - DWK changed public double totalc[CYCLE] to 
           double totalc
20051130 - DWK deleted public function yearlyTransient()
20051201 - DWK deleted public function deltaxclm()
20051201 - DWK deleted const double& tgppopt from function call
           to monthxclm()
20051201 - DWK deleted private double gv
20051201 - DWK deleted private double ksoil
20051201 - DWK deleted private double temp
20051202 - DWK deleted private double respq10
20051202 - DWK deleted private double dq10
20051202 - DWK deleted public function monthxclm()
20051202 - DWK renamed public function ECDsetELMNTstate() to
           ECDsetODEstate()
20051202 - DWK renamed public function setELMNTflux() to
           resetMonthlyELEMNTFluxes()
20051202 - DWK added const int& pdm to 
           resetMonthlyELEMNTFluxes()                                                                                      
20051202 - DWK renamed public function resetYrFlux() to
           resetYrFluxes()
20051207 - DWK deleted public function resetVegEcd()
20051208 - DWK added private function pcdisplayMonth()
20051208 - DWK added private function pcdisplayODEerr()
20051208 - DWK added private function pcdisplayDT()
20051208 - DWK added public int topwind and int calwind
20051208 - DWK added enum scykey, enum snykey, enum swykey,
           enum sstykey, and enum sgykey and associated 
           public functions
20051208 - DWK added public functions getOptionalCflx(),
           getOptionalNflx(), getOptionalSoilTemp(), 
           getOptionalTraceGas() and getOptionalWflx()
20051208 - DWK added global const int WSY, const int CSY,
           const int NSY, const int STSY, and const int GSY
20051208 - DWK added public scykey scy[CSY], snykey sny[NSY],
           swykey swy[WSY], sstykey ssty[STSY], and 
           sgykey sgy[GSY]
20051208 - DWK added public functions displayOptionalCflx(),
           displayOptionalNflx(), displayOptionalSoilTemp(),
           displayOptionalTraceGas() and displayOptionalWflx()
20051208 - DWK added public int adapttol, int intbomb and
           int tolbomb 
20060609 - DWK added GET_FIRENDEP to enum snykey
20060909 - DWK changed include from tsoil602.h to tsoil603.h
20061108 - DWK added private function urbanDynamics()
20061112 - DWK added private function pastureDynamics()
20070210 - DWK changed include from tveg602.h to tveg603.h
20070210 - DWK changed include from humnact602.h to humnact603.h
20070829 - DWK changed include from tveg603.h to tveg603b.h
20070830 - DWK changed include from tsoil603b.h to tsoil603c.h
20090127 - DWK changed include from tsoil603c.h to tsoil604.h
20090128 - DWK changed include from tveg603b.h to tveg604.h
20090128 - DWK changed include from humnact603b.h to humnact604.h                                   

****************************************************************

References:

VERSION 4.1

Tian, H., J.M. Melillo, D.W. Kicklighter, A.D. McGuire and J.
  Helfrich.  1999. The sensitvity of terrestrial carbon storage to
  historical climate variability and atmospheric CO2 in the United
  States.  Tellus 51B: 414-452.

VERSION 4.2

McGuire, A.D., S. Sitch, J.S. Clein, R. Dargaville, G. Esser, J. Foley,
  M. Heimann, F. Joos, J. Kaplan, D.W. Kicklighter, R.A. Meier, J.M.
  Melillo, B. Moore III, I.C. Prentice, N. Ramankutty, T. Reichenau,
  A. Schloss, H. Tian, L.J. Williams and U. Wittenberg (2001) Carbon
  balance of the terrestrial biosphere in the twentieth century:
  analyses of CO2, climate and land use effects with four process-
  based ecosystem models.  Global Biogeochemical Cycles 15: 183-206.

Tian, H., J.M. Melillo, D.W. Kicklighter, S. Pan, J. Liu, A.D. McGuire
  and B. Moore III (2003) Regional carbon dynamics in monsoon Asia
  and its implications for the global carbon cycle. Global and
  Planetary Change 37: 201-217.

VERSION 4.3

Felzer, B., D. Kicklighter, J. Melillo, C. Wang, Q. Zhuang, and 
  R. Prinn (2004) Effects of ozone on net primary production and 
  carbon sequestration in the conterminous United States using a 
  biogeochemistry model. Tellus 56B: 230-248.

VERSION 5.0

Zhuang, Q., A.D. McGuire, J.M. Melillo, J.S. Clein, R.J.
  Dargaville, D.W. Kicklighter, R.B. Myneni, J. Dong, V.E.
  Romanovsky, J. Harden and J.E. Hobbie (2003) Carbon cycling in
  extratropical terrestrial ecosystems of the Northern Hemisphere
  during the 20th century: a modeling analysis of the influences
  of soil thermal dynamics. Tellus 55B: 751-776.

Runge - Kutta - Fehlberg (RKF) adaptive integrator:

Cheney, W., and D. Kincaid.  1985.  Numerical mathematics and
  computing.  Second edition.  Brooks/ Col Publishing Co.  Monterey,
  CA.  pp. 325-328.

****************************************************************
************************************************************** */

#ifndef TTEM604_H
#define TTEM604_H
 
#include "temconsts602.hpp"

#ifdef CALIBRATE_TEM
  // Additional global constants

  const int WSY = 5;
  const int CSY = 6;
  const int NSY = 3;
  const int STSY = 2;
  const int GSY = 4;
#endif


// Objects describing basic components of the ecosystem

#include "bioms423.hpp"   // TTEM uses Biomass class


// Objects describing the structure of the ecosystem

#include "atms602.h"    // TTEM uses Atmosphere60 class
#include "tveg604.h"   // TTEM uses Tveg60 class
#include "tsoil604.h"   // TTEM uses Tsoil60 class
#include "tmcrb604.h"   // TTEM uses Tmicrobe60 class

// Objects describing the effects of human activities on the ecosystem

#include "humnact604.h" // TTEM uses Humnact60 class


class TTEM60
{

  public:

     TTEM60();

     enum temkey
     {
       I_VEGC,     I_SOLC,     I_DOC,
       I_STRN,     I_STON,     I_SOLN,     I_DON,   I_NH4,
       I_NO3,      I_SM,       I_RGRW,     I_SGRW, I_DEADWOODC,   I_DEADWOODN,
       I_FOLC,     I_STEMC,    I_CROOTC,   I_FROOTC, I_CWD, I_AGR,
       I_AGL,      I_BGR,      I_BGL,    I_MDC,   I_AGE, I_SOC, //added by cgs2014

       I_UNRMLF,   I_LEAF,     I_LAI,      I_FPC,

       I_INGPP,    I_GPP,      I_FOZONE,   I_FINDOZONE,I_INNPP,    
       I_NPP,      I_GPR,      I_RVMNT,    I_RVGRW,    I_ABVGPR,   
       I_ROOTGPR,  I_LTRC,     I_CDCMP,    I_RH,       I_DOCP,     
       I_LCHDOC,   I_ERDPOC,
             

       I_AGFRTN,   I_BNFIX,    I_SNFIX,    I_ANFIX,    I_INNUP,
       I_INNH4UP,  I_INNO3UP,  I_VNUP,     I_VNH4UP,   I_VNO3UP,   
       I_VSUP,     I_VLUP,     I_VNMBL,    I_VNRSRB,   I_LTRN,     
       I_NDCMP,    I_DONP,     I_GMIN,     I_NH4IMM,   I_NIMM,     
       I_NMIN,     I_AIMMNH4,  I_AIMMNO3,  I_AMMN,     I_NTRF,     
       I_NO3P,     I_NOP,      I_N2OP,     I_N2P,      I_DNTRF,
       I_NH3FLX,   I_NOFLX,    I_N2OFLX,   I_N2FLX,    I_LCHNO3,   
       I_LCHDON,   I_ERDPON,

       I_AGIRRIG,  I_INEET,    I_EET,      I_RPERC,    I_SPERC,    
       I_RRUN,     I_SRUN,     

       I_NSOLC,    I_TSOLC,    I_TOTEC,    I_TOTC, 
       
       I_VEGN,     I_NSOLN,    I_TSOLN,    I_AVLN,
       
       I_SNWPCK,   I_AVLW,     I_VSM,      I_PCTP,     
       
       I_RSOIL,    I_NEP,      I_NCE,      I_NTCB, 
       
       I_NINP,     I_NLST,     I_NTNB,
       
       I_PET,      I_SNWINF,   I_WYLD,
       
       I_AGPRDC,   I_PROD10C,  I_PROD100C, I_TOTPRDC,  
       
       I_RESIDC,   I_AGSTUBC,     
       
       I_AGPRDN,   I_PROD10N,  I_PROD100N, I_TOTPRDN,
       
       I_RESIDN,   I_AGSTUBN,
       
       I_CNVRTC,   I_VCNVRTC,  I_SCNVRTC,  I_SLASHC,   I_CFLX,
       I_FFLC, I_FFDC, I_SWFC, I_DWFC,
       
       I_CNVRTN,   I_VCNVRTN,  I_SCNVRTN,  I_SLASHN,   I_NRETNT,
       I_NVRTNT,   I_NSRTNT,

       I_AGFPRDC,  I_AGFPRDN,  I_FRESIDC,  I_FRESIDN,  I_AGPRDFC,
       I_AGPRDFN,  I_RESIDFC,  I_RESIDFN,  
       
       I_PRDF10C,  I_PRDF10N,  I_PRD10FC,  I_PRD10FN,  I_PRDF100C,  
       I_PRDF100N, I_PRD100FC, I_PRD100FN, I_TOTFPRDC, I_TOTFPRDN, 
       I_TOTPRDFC, I_TOTPRDFN,

       I_CROPC,    I_NATVEGC,  I_CROPN,    I_NATVEGN,  I_CSTRN,
       I_NATSTRN,  I_CSTON,    I_NATSTON,

       I_CROPULF,  I_NATULF,   I_CROPLEAF, I_NATLEAF, I_CROPLAI,
       I_NATLAI,   I_CROPFPC,  I_NATFPC,

       I_AGINGPP,  I_NATINGPP, I_AGGPP,    I_NATGPP,   I_AGINNPP,
       I_NATINNPP, I_AGNPP,    I_NATNPP,   I_AGGPR,    I_NATGPR,
       I_AGRVMNT,  I_NATRVMNT, I_AGRVGRW,  I_NATRVGRW, I_AGLTRC,
       I_NATLTRC,

       I_AGSNFX,   I_NATSNFX,  I_AGINNUP,  I_NATINNUP, I_AINNH4UP,
       I_NINNH4UP, I_AINNO3UP, I_NINNO3UP, I_AGVNUP,   I_NATVNUP,
       I_AGVNH4UP, I_NVNH4UP,  I_AGVNO3UP, I_NVNO3UP,  I_AGVSUP,
       I_NATVSUP,  I_AGVLUP,   I_NATVLUP,  I_AGVNMBL,  I_NATVNMBL,
       I_AGVNRSRB, I_NVNRSRB,  I_AGLTRN,   I_NATLTRN,

       I_TSOIL,    I_DST0,     I_DST5,     I_DST10,    I_DST20,
       I_DST50,    I_DST100,   I_DST200,  I_DST300, I_FRONTD,   I_THAWBE,
       I_THAWEND,  I_THAWPCT,  I_ACTLAYER,

       I_CO2G,     I_CO2W,     I_HCO3,     I_RHCO3,    I_ALK, 
       I_CO2DISS,  I_LCHCO2,   I_HCO3P,    I_RHCO3P,   I_LCHHCO3,  
       I_LCHALK, 
       
       I_ISOPREN,  I_TERPEN,   I_ORVOC,    I_OVOC,     I_VOC
     };

     #ifdef CALIBRATE_TEM    
       enum scykey { NOCKEY,       GET_LEAF,    GET_LAI,     
                     GET_FPC,      GET_NSOLC,   GET_TSOLC,

                     GET_INGPP,    GET_GPP,     GET_INNPP,
                     GET_NPP,      GET_GPR,     GET_RVMNT,
                     GET_RVGRW,    GET_LTRC,    GET_AGSTUBC,
                     GET_CDCMP,    GET_RH,      GET_RSOIL,   
                     GET_LCHDOC,   GET_NEP,     GET_NTCB,

                     GET_D40,      GET_FOZONE,


                     GET_CNVRTC,   GET_VCNVRTC, GET_SCNVRTC,
                     GET_SLASHC,   GET_CFLX,    GET_NCE,

                     GET_AGPRDC,   GET_PROD10C, GET_PROD100C,
                     GET_RESIDC,

                     GET_AGFPRDC,  GET_PRDF10C, GET_PRDF100C,
                     GET_FRESIDC,

                     GET_AGPRDFC,  GET_PRD10FC, GET_PRD100FC,
                     GET_TOTPRDFC, GET_RESIDFC };

       enum snykey { NONKEY,       GET_NH4,     GET_NO3,
                     GET_DON,      GET_NSOLN,   GET_TSOLN,

                     GET_NINP,     GET_TNDEP,   GET_NH4DEP,
                     GET_NO3DEP,   GET_AGFRTN,  GET_BNFIX,
                     GET_SNFIX,    GET_ANFIX,   GET_INNUP,   
                     GET_INNH4UP,  GET_INNO3UP, GET_VNUP,    
                     GET_VNH4UP,   GET_VNO3UP,  GET_VSUP,    
                     GET_VLUP,     GET_VNMBL,   GET_VNRSRB,  
                     GET_LTRN,     GET_AGSTUBN, GET_NDCMP,   
                     GET_DONP,     GET_GMIN,    GET_NH4IMM,  
                     GET_NO3IMM,   GET_NIMM,    GET_NMIN, 
                     GET_LCHNO3,   GET_LCHDON,  GET_NLST,    
                     GET_NTNB,

                     GET_CNVRTN,   GET_VCNVRTN, GET_SCNVRTN,
                     GET_SLASHN,   GET_NRETNT,  GET_NVRTNT,
                     GET_NSRTNT,

                     GET_AGPRDN,   GET_PROD10N, GET_PROD100N,
                     GET_RESIDN,

                     GET_AGFPRDN,  GET_PRDF10N, GET_PRDF100N,
                     GET_FRESIDN,
                     GET_AGPRDFN,  GET_PRD10FN, GET_PRD100FN,
                     GET_TOTPRDFN, GET_RESIDFN, GET_FIRENDEP,

                     GET_L2SN };

       enum swykey { NOWKEY,    GET_SH2O,   GET_PCTP,   GET_VSM,

                     GET_RAIN,  GET_SNWFAL, GET_SNWINF, GET_AGIRRIG,
                     GET_PET,   GET_INEET,  GET_EET,    GET_RPERC,
                     GET_SPERC, GET_RRUN,  GET_SRUN,    GET_WYLD };

       enum sstykey { NOSTKEY,     GET_FRONTD,  GET_ACTLAY,
                      GET_THAWPCT, GET_THWBEG1, GET_THWEND1,  
                      GET_THWBEG2, GET_THWEND2 };

       enum sgykey { NOGKEY,      GET_AIMMNH4,  GET_AMMN,   
                     GET_AIMMNO3, GET_NO3P,     GET_GNLST,
                     GET_NOP,     GET_N2OP,    GET_N2P,      
                     GET_NH3FLX,  GET_NOFLX,   GET_N2OFLX,   
                     GET_N2FLX,   GET_GNMIN };
     #endif

/* **************************************************************
			Public Functions
************************************************************** */

     void askODE( ofstream& rflog1 );

     #ifdef CALIBRATE_TEM
       void displayOptionalCflx( const scykey& s );
       void displayOptionalNflx( const snykey& s );
       void displayOptionalSoilTemp( const sstykey& s );
       void displayOptionalTraceGas( const sgykey& s );
       void displayOptionalWflx( const swykey& s );
     #endif
     
     int ecdqc( const int& dcmnt );

     void ECDsetODEstate( const int& pdcmnt,
                          const double& psiplusc );
 
     void getco2( void );

     void getcropecd( const int& dv, const string& agecd );


     #ifdef CALIBRATE_TEM
       double getOptionalCflx( const int& optflx );

       double getOptionalNflx( const int& optflx );

       double getOptionalSoilTemp( const int& optflx );

       double getOptionalTraceGas( const int& optflx );

       double getOptionalWflx( const int& optflx );
     #endif


     void getsitecd( const int& numcmnt, ofstream& rflog1 );

     void getsitecd( const int& dv, const string& ecd );
     
     void initrun( ofstream& rflog1 );

     int monthlyTransient( const int& outyr,
                           const int& pdm,
                           const double& outtol,
                           ofstream& rflog1 );

     #ifdef CALIBRATE_TEM
       inline scykey& next( scykey& s )
       {
         return s = (GET_RESIDFC == s) ? GET_LEAF : scykey( s+1 );
       }

       inline sgykey& next( sgykey& s )
       {
         return s = (GET_GNMIN == s) ? GET_AIMMNH4 : sgykey( s+1 );
       }

       inline snykey& next( snykey& s )
       {
         return s = (GET_L2SN == s) ? GET_NH4 : snykey( s+1 );
       }

       inline sstykey& next( sstykey& s )
       {
         return s = (GET_THWEND2 == s) ? GET_FRONTD : sstykey( s+1 );
       }

       inline swykey& next( swykey& s )
       {
         return s = (GET_WYLD == s) ? GET_SH2O : swykey( s+1 );
       }

       inline scykey& prev( scykey& s )
       {
         return s = (GET_LEAF == s) ? GET_RESIDFC : scykey( s-1 );
       }

       inline sgykey& prev( sgykey& s )
       {
         return s = (GET_AIMMNH4 == s) ? GET_GNMIN : sgykey( s-1 );
       }

       inline snykey& prev( snykey& s )
       {
         return s = (GET_NH4 == s) ? GET_L2SN : snykey( s-1 );
       }

       inline sstykey& prev( sstykey& s )
       {
         return s = (GET_FRONTD == s) ? GET_THWEND2 : sstykey( s-1 );
       }

       inline swykey& prev( swykey& s )
       {
         return s = (GET_SH2O == s) ? GET_WYLD : swykey( s-1 );
       }
     #endif

     void resetMonthlyELMNTFluxes( void );

     void resetYrFluxes( void );

     void setELMNTecd( const int& pdcmnt,
                       const double& psiplusc );

     void setEquilC2N( const int& pdcmnt,
                       const double& co2 );

     void setEquilEvap( const double& nirr,
                        const double& tair,
                        const int& pdm );

     void setPrevState( void );
     
     int stepmonth( const int& dyr,
                    const int& pdm,
                    int& intflag,
                    const double& tol,
                    ofstream& rflog1 );

     int testEquilibrium( void );

     void updateYearSummary( void );
     
     void writesitecd( ofstream& fout, const int& dcmnt );


     // "Get" and "Set" private variables and parameters   
     
     // bnfix **************************************************
     
     inline double getBNFIX( void ) { return bnfix; }


     // doccut ***************************************************
     
     inline double getDOCCUT( const int& pcmnt )
     {
       return doccut[pcmnt];
     }

     inline void setDOCCUT( const double& pdoccut, 
                            const int& pcmnt )
     {
       doccut[pcmnt] = pdoccut;
     }

     // doc1a ***************************************************
     
     inline double getDOC1A( const int& pcmnt )
     {
       return doc1a[pcmnt];
     }

     inline void setDOC1A( const double& pdoc1a, 
                           const int& pcmnt )
     {
       doc1a[pcmnt] = pdoc1a;
     }

     // doc1b ***************************************************
     
     inline double getDOC1B( const int& pcmnt )
     {
       return doc1b[pcmnt];
     }

     inline void setDOC1B( const double& pdoc1b, 
                           const int& pcmnt )
     {
       doc1b[pcmnt] = pdoc1b;
     }


     // doc2a ***************************************************
     
     inline double getDOC2A( const int& pcmnt )
     {
       return doc2a[pcmnt];
     }

     inline void setDOC2A( const double& pdoc2a, 
                           const int& pcmnt )
     {
       doc2a[pcmnt] = pdoc2a;
     }

     // doc2b ***************************************************
     
     inline double getDOC2B( const int& pcmnt )
     {
       return doc2b[pcmnt];
     }

     inline void setDOC2B( const double& pdoc2b, 
                           const int& pcmnt )
     {
       doc2b[pcmnt] = pdoc2b;
     }


      // doncut ***************************************************
     
     inline double getDONCUT( const int& pcmnt )
     {
       return doncut[pcmnt];
     }

     inline void setDONCUT( const double& pdoncut, 
                            const int& pcmnt )
     {
       doncut[pcmnt] = pdoncut;
     }


     // don1a ***************************************************
     
     inline double getDON1A( const int& pcmnt )
     {
       return don1a[pcmnt];
     }

     inline void setDON1A( const double& pdon1a, 
                           const int& pcmnt )
     {
       don1a[pcmnt] = pdon1a;
     }

     // don1b ***************************************************
     
     inline double getDON1B( const int& pcmnt )
     {
       return don1b[pcmnt];
     }

     inline void setDON1B( const double& pdon1b, 
                           const int& pcmnt )
     {
       don1b[pcmnt] = pdon1b;
     }


     // don2a ***************************************************
     
     inline double getDON2A( const int& pcmnt )
     {
       return don2a[pcmnt];
     }

     inline void setDON2A( const double& pdon2a, 
                           const int& pcmnt )
     {
       don2a[pcmnt] = pdon2a;
     }

     // don2b ***************************************************
     
     inline double getDON2B( const int& pcmnt )
     {
       return don2b[pcmnt];
     }

     inline void setDON2B( const double& pdon2b, 
                           const int& pcmnt )
     {
       don2b[pcmnt] = pdon2b;
     }


     // maxit ***************************************************
     
     inline int getMAXIT( void )
     {
       return maxit;
     }

     inline void setMAXIT( const int& pmaxit )
     {
       maxit = pmaxit;
     }


     // maxitmon*************************************************
     
     inline long getMAXITMON( void )
     {
       return maxitmon;
     }

     inline void setMAXITMON( const long& pmaxitmon )
     {
       maxitmon = pmaxitmon;
     }


    // nce ****************************************************
     
     inline double getNCE( void ) { return nce; }


     // nep ****************************************************
     
     inline double getNEP( void ) { return nep; }


     // nh4a ***************************************************
     
     inline double getNH4A( const int& pcmnt )
     {
       return nh4a[pcmnt];
     }

     inline void setNH4A( const double& pnh4a, 
                          const int& pcmnt )
     {
       nh4a[pcmnt] = pnh4a;
     }

     // nh4b ***************************************************
     
     inline double getNH4B( const int& pcmnt )
     {
       return nh4b[pcmnt];
     }

     inline void setNH4B( const double& pnh4b, 
                          const int& pcmnt )
     {
       nh4b[pcmnt] = pnh4b;
     }

      // no3cut ***************************************************
     
     inline double getNO3CUT( const int& pcmnt )
     {
       return no3cut[pcmnt];
     }

     inline void setNO3CUT( const double& pno3cut, 
                            const int& pcmnt )
     {
       no3cut[pcmnt] = pno3cut;
     }


     // no31a **************************************************
     
     inline double getNO31A( const int& pcmnt )
     {
       return no31a[pcmnt];
     }

     inline void setNO31A( const double& pno31a, 
                           const int& pcmnt )
     {
       no31a[pcmnt] = pno31a;
     }

     // no31b **************************************************
     
     inline double getNO31B( const int& pcmnt )
     {
       return no31b[pcmnt];
     }

     inline void setNO31B( const double& pno31b, 
                           const int& pcmnt )
     {
       no31b[pcmnt] = pno31b;
     }

     // no32a **************************************************
     
     inline double getNO32A( const int& pcmnt )
     {
       return no32a[pcmnt];
     }

     inline void setNO32A( const double& pno32a, 
                           const int& pcmnt )
     {
       no32a[pcmnt] = pno32a;
     }

     // no32b **************************************************
     
     inline double getNO32B( const int& pcmnt )
     {
       return no32b[pcmnt];
     }

     inline void setNO32B( const double& pno32b, 
                           const int& pcmnt )
     {
       no32b[pcmnt] = pno32b;
     }

     // ntcb ***************************************************
     
     inline double getNTCB( void ) { return ntcb; }


     // ntnb ***************************************************
     
     inline double getNTNB( void ) { return ntnb; }


     // prevy **************************************************
     
     inline double getPREVY( const int& i ) { return prevy[i]; }

     inline void setPREVY( const double& pprevy, const int& i )
     {
       prevy[i] = pprevy;
     }


     // rsoil **************************************************
     
     inline double getRSOIL( void ) { return rsoil; }


     // solca **************************************************
     
     inline double getSOLCA( const int& pcmnt )
     {
       return solca[pcmnt];
     }

     inline void setSOLCA( const double& psolca, 
                           const int& pcmnt )
     {
       solca[pcmnt] = psolca;
     }


     // solcb **************************************************
     
     inline double getSOLCB( const int& pcmnt )
     {
       return solcb[pcmnt];
     }

     inline void setSOLCB( const double& psolcb, 
                           const int& pcmnt )
     {
       solca[pcmnt] = psolcb;
     }


     // solna **************************************************
     
     inline double getSOLNA( const int& pcmnt )
     {
       return solna[pcmnt];
     }

     inline void setSOLNA( const double& psolna, 
                           const int& pcmnt )
     {
       solna[pcmnt] = psolna;
     }


     // solnb **************************************************
     
     inline double getSOLNB( const int& pcmnt )
     {
       return solnb[pcmnt];
     }

     inline void setSOLNB( const double& psolnb, 
                           const int& pcmnt )
     {
       solnb[pcmnt] = psolnb;
     }


     // stona **************************************************
     
     inline double getSTONA( const int& pcmnt )
     {
       return stona[pcmnt];
     }

     inline void setSTONA( const double& pstona, 
                           const int& pcmnt )
     {
       stona[pcmnt] = pstona;
     }


     // stonb **************************************************
     
     inline double getSTONB( const int& pcmnt )
     {
       return stonb[pcmnt];
     }

     inline void setSTONB( const double& pstonb, 
                           const int& pcmnt )
     {
       stonb[pcmnt] = pstonb;
     }


     // strna **************************************************
     
     inline double getSTRNA( const int& pcmnt )
     {
       return strna[pcmnt];
     }

     inline void setSTRNA( const double& pstrna, 
                           const int& pcmnt )
     {
       strna[pcmnt] = pstrna;
     }


     // strnb **************************************************
     
     inline double getSTRNB( const int& pcmnt )
     {
       return strnb[pcmnt];
     }

     inline void setSTRNB( const double& pstrnb, 
                           const int& pcmnt )
     {
       strnb[pcmnt] = pstrnb;
     }


     // totalc *************************************************
     
     inline double getTOTALC( void ) { return totalc; }


     // vegca **************************************************
     
     inline double getVEGCA( const int& pcmnt )
     {
       return vegca[pcmnt];
     }

     inline void setVEGCA( const double& pvegca, 
                           const int& pcmnt )
     {
       vegca[pcmnt] = pvegca;
     }


     // vegcb **************************************************
     
     inline double getVEGCB( const int& pcmnt )
     {
       return vegcb[pcmnt];
     }

     inline void setVEGCB( const double& pvegcb, 
                           const int& pcmnt )
     {
       vegcb[pcmnt] = pvegcb;
     }


     // y ******************************************************
          
     inline double getY( const int& i ) { return y[i]; }

     inline void setY( const double& py, const int& i )
     {
       y[i] = py;
     }

     inline void setLAT( const double& tlatitude )
     {
        latitude = tlatitude;
     }
     inline void setLON( const double& tlongitude )
     {
        longitude = tlongitude;
     }


/* *************************************************************
			 Public Variables
************************************************************** */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif

     // Input Ecd files with calibration data
     
     #ifdef CALIBRATE_TEM
       string soilfile;
       string rootfile;
       string vegfile;
       string leaffile;
       string mcrvfile;
       string agfile;
       string snowfile;
       string slayerfile;
 //      string tsoilfile;
     #endif
     
     #ifdef CALIBRATE_TEM
       int adapttol;
     #endif
     
     Humanact60 ag;

     Atmosphere60 atms;

     static int avlnflag;
     
     static int baseline;

     #ifdef CALIBRATE_TEM
       int calwind;
     #endif
     
     int climitflg;

     static double ctol;

     static int diffyr;

     // Counter of months since disturbance
     int distmnthcnt;
     
     // Flag to indicate that a disturbance has occurred
     //   in a particular year.  Different values represent
     //   different types of disturbances:
     //   0 = no disturbance
     //   1 = sustained disturbed state 
     //         (e.g., conversion to agriculture or urban areas)
     //   2 = timber harvest
     //   3 = fire
     //   4 = insect & disease
     //   5 = land conversion among natural vegetation types
     int disturbflag;
     
     // Month in which disturbance occurred
     // (e.g. 1 = January, ...., 12 = December)
     int disturbmonth;   

     double elev;              // elevation (m)

     static int endeq;

     static int modstartyr; //model simulation start year
     static int outputstartyr;

     static int outputendyr;

     static int equil;

     // index for vegetation type (ez = veg.temveg-1)
     int ez;

     // Counter of months since fire disturbance
     double firemnthcnt;

     static int initbase;

     int initFlag;

     static double inittol;

     #ifdef CALIBRATE_TEM
       int intbomb;
     #endif

     static int intflag;

     static int maxnrun;

     static int maxyears;

     Tmicrobe60 microbe;

     static int moistlim;

     int nattempt;

     static int nfeed;

     double nfert;

     static double ntol;

     static int o3flag;

     int predflag;

     vector<string> predstr;

     int qualcon[MAXRTIME];

     int retry;

     static int rheqflag;

     static int runsize;

     #ifdef CALIBRATE_TEM
       scykey scy[CSY];
       sgykey sgy[GSY];
       snykey sny[NSY];
       sstykey ssty[STSY];
       swykey swy[WSY];
     #endif
     
     Tsoil60 soil;


     static int strteq;

     double tol;

     #ifdef CALIBRATE_TEM
       int tolbomb;
    
       int topwind;
     #endif

     int totyr;

     Tveg60 veg;

     static int wrtyr;

     static double wtol;

     double yrbnfix;  // (g N / (sq. meter * year))

     double yrnce;

     double yrnep;    // (g C / (sq. meter * year))

     double yrntcb;   // (g C / (sq. meter * year))

     double yrntnb;   // (g C / (sq. meter * year))
     
     double yrrsoil;  // (g C / (sq. meter * year))

     double yrtotalc;

     // Site ECD variables

     string version;
     string sitename;
     string developer;
     double sitecol;
     double siterow;
     long updated;
     string description;
     double randomnum;
     int newstandage;


/* **************************************************************
		     Protected Variables
************************************************************** */

//  protected:

     int dbugflg;


  private:

/* **************************************************************
		 Private Functions
************************************************************** */

     int adapt( const int& numeq, 
                double pstate[], 
                const double& ptol, 
                const int& pdm );

     int boundcon( double ptstate[],
                   double err[],
                   const double& ptol );

     void delta( const int& dm,
                 double pstate[],
                 double pdstate[] );

     void cropDynamics( const int& dm, double pstate[] );

     void getenviron( const int& pdyr,
                      const int& pdm,
                      ofstream& rflog1 );

     void massbal( double y[NUMEQ], 
                   double prevy[MAXSTATE] );

     void natvegDynamics( const int& dm, double pstate[] );

     void pastureDynamics( const int& dm, double pstate[] );

     #ifdef CALIBRATE_TEM
       void pcdisplayDT( const double& tottime,
                         const double& deltat );

       void pcdisplayMonth( const int& dyr,
                            const int& dm );

       void pcdisplayODEerr( const int& test, 
                             double pstate[] );

	#endif
     void resetODEflux( void );

     void rkf( const int& numeq, 
               double pstate[], 
               double& pdt, 
               const int& pdm );
                 
     void step( const int& numeq, 
                double pstate[], 
                double pdstate[], 
                double ptstate[], 
                double& pdt );

     void urbanDynamics( const int& dm, double pstate[] );
     double omitzero(const double& data1);

/* **************************************************************
		     Private Variables
************************************************************** */

     // Biological Nitrogen Fixation      
     double bnfix;    // (g N / (sq. meter * month))

     ifstream fecd[MAXCMNT];

     // Net Carbon Exchange (NCE)
     double nce;

     // Net Ecosystem Production
     double nep;      // (g C / (sq. meter * month))

     // Net Terrestrial Carbon Balance (NTCB)
     double ntcb;     // (g C / (sq. meter * month))

     // Net Terrestrial Nitrogen Balance (NTNB)
     double ntnb;     // (g C / (sq. meter * month))

     // Values of ODE state variables for previous month 
     double prevy[MAXSTATE];

     // Soil respiration (heterotropic + root respiration)
     double rsoil;    // (g C / (sq. meter * month))

     // Total carbon storage (veg.plant + soil.org)
     double totalc;   // (g C/sq. meter)

     // Values of ODE state variables for current month
     double y[NUMEQ];

    
     // Adaptive integrator variables

     int blackhol;
     static int maxit;
     static long maxitmon;

     int syint;
     int test;

     double dum4[NUMEQ];
     double error[NUMEQ];
     
     double dum5[NUMEQ];
     double dumy[NUMEQ];

     double ydum[NUMEQ];
     double dy[NUMEQ];
     double yprime[NUMEQ];
     double rk45[NUMEQ];

     double f11[NUMEQ];
     double f2[NUMEQ];
     double f13[NUMEQ];
     double f3[NUMEQ];
     double f14[NUMEQ];
     double f4[NUMEQ];
     double f15[NUMEQ];
     double f5[NUMEQ];
     double f16[NUMEQ];
     double f6[NUMEQ];
     double f111[NUMEQ]; //added by cgs2014


     static double  a1;

     static double  a3;
     static double a31;
     static double a32;

     static double  a4;
     static double a41;
     static double a42;
     static double a43;

     static double  a5;
     static double a51;
     static double a52;
     static double a53;
     static double a54;

     static double  b1;
     static double  b3;
     static double  b4;
     static double  b5;

     static double  b6;
     static double b61;
     static double b62;
     static double b63;
     static double b64;
     static double b65;

     double dummy;
     static double tempvara1;


/* **************************************************************
			Private Parameters
************************************************************** */

     double vegca[MAXCMNT];
     double vegcb[MAXCMNT];

     double strna[MAXCMNT];
     double strnb[MAXCMNT];

     double stona[MAXCMNT];
     double stonb[MAXCMNT];

     double solca[MAXCMNT];
     double solcb[MAXCMNT];

     double solna[MAXCMNT];
     double solnb[MAXCMNT];

     double doccut[MAXCMNT];
     double doc1a[MAXCMNT];
     double doc1b[MAXCMNT];
     double doc2a[MAXCMNT];
     double doc2b[MAXCMNT];

     double doncut[MAXCMNT];
     double don1a[MAXCMNT];
     double don1b[MAXCMNT];
     double don2a[MAXCMNT];
     double don2b[MAXCMNT];

     double nh4a[MAXCMNT];
     double nh4b[MAXCMNT];

     double no3cut[MAXCMNT];
     double no31a[MAXCMNT];
     double no31b[MAXCMNT];
     double no32a[MAXCMNT];
     double no32b[MAXCMNT];
     double latitude;
     double longitude;


 };

#endif
