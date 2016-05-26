/* *************************************************************
****************************************************************
TTEM604.CPP - Terrestrial Ecosystem Model Version 6.0
****************************************************************

Modifications:

20030712 - DWK created by modifying ttem50b1.cpp to incorporate
           open-N cycle dynamics into TEM
20030712 - DWK changed TTEM50:: to TTEM51::
20030712 - DWK commented out "o3flag = 0" in cropDynamics() and
           natvegDynamics()
20030712 - DWK added veg.o3para[], veg.o3parb[] and veg.o3parc[]
           to ecdqc(), getcropecd() and getsitecd()
20030712 - DWK added TOTNDEP, VNFIX, MNFIX, GROSMIN, NIMMOB, NOPROD,
           N2OPROD, NO3PROD, NITRIF, AMMNVOL, DNITRIF, LEACH,
           CRPVNFIX and TOTVNFIX to constructor
20030719 - DWK added RSOIL, LEACHDOC, LEACHDON to TTEM51()
20030719 - DWK renamed LEACH as LEACHNO3 in TTEM51()
20030719 - DWK added I_NATVEGN, I_NATSTRN, I_NATSTON, I_NATULF,
           I_NATLEAF, I_NATLAI, I_NATFPC,  I_NATINGPP, I_NATGPP,
           I_NATINNPP, I_NATGPR, I_NATRVMNT, I_NATRVGRW,
           I_NATINNUP, I_AINNH4UP, I_NINNH4UP, I_AINNO3UP,
           I_NINNO3UP, I_NATVNUP, I_AGVNH4UP, I_NVNH4UP,
           I_AGVNO3UP, I_NVNO3UP, I_NATVSUP, I_NATVLUP,
           I_NATVNMBL, I_NVNRSB and I_NATLTRN to TTEM51()
20030719 - DWK changed I_TOTNPP to I_NATNPP in TTEM51()
20030719 - DWK changed TOTVEGC to NATVEGC in TTEM51()
20030719 - DWK changed I_AVLN to I_NH4 in boundcon()
20030719 - DWK added I_NO3 to boundcon()
20030719 - DWK changed I_VNUP to I_VNH4UP in boundcon()
20030719 - DWK added I_VNO3UP to boundcon()
20030720 - DWK added veg.nfix[] to cropDynamics() and
           natvegDynamics()
20030720 - DWK added microbe.ammnvol[], microbe.nitrif[],
           microbe.no3prod[], microbe.noprod[], microbe.n2oprod[],
           microbe.n2prod[], microbe.denitrif[], soil.nh3flux[],
           soil.noflux[], soil.n2oflux[], soil.n2flux[],
           soil.leachNO3[] and soil.leachDOM[] to cropDynamics(),
           delta() and natvegDynamics(), setELMNTflux() and setMonth()
20030720 - DWK replaced kda with kd1a and kd2a in ecdqc(),
           getcropecd(), getsitecd() and writesitecd()
20030720 - DWK replaced kdb with kd1b and kd2b in ecdqc(),
           getcropecd(), getsitecd() and writesitecd()
20030720 - DWK added kdcut to ecdqc(), getcropecd() and getsitecd()
20030720 - DWK changed avlna[] to nh4a[], changed avlnb[] to
           nh4b[], and added no3a[] and no3b[] to
           ecdqc(), getcropecd(), getsitecd() and writesitecd()
20030720 - DWK changed nmaxcut[] to nupnh4cut[], changed nmax1a[]
           to nupnh41a, changed nmax1b[] to nupnh41b[], changed
           nmax2a[] to nupnh42a, changed nmax2b[] to nupnh42b[],
           and added nupno3cut[], nupno31a[], nupno31b[],
           nupno32a[] and nupno32b[] to ecdqc(), getcropecd(),
           getsitecd() and writesitecd()
20030720 - DWK added veg.rroot[], soil.DOCpar[],
           soil.lchDOMpar[], veg.nfixpara[], veg.nfixparb[],
           microbe.nfixpara[], microbe.nfixparb[],
           microbe.ammnpar[], microbe.maxntrfCN[],
           microbe.minntrfCN[], soil.lchNO3parcut[],
           soil.lchNO3par1a[], soil.lchNO3par1b[],
           soil.lchNO3par2a[], soil.lchNO3par2b[] and
           soil.DONpar[] to ecdqc(), getcropecd(), getsitecd()
           and writesitecd()
20030720 - DWK changed microbe.nupa[] to microbe.nimma[] and
           microbe.nupb[] to microbe.nimmb[] in ecdqc(),
           getcropecd(), getsitecd() and writesitecd()
20030720 - DWK deleted microbe.nfixpar[] and soil.nloss[] from
           ecdqc(), getcropecd(), getsitecd() and writesitecd()
20030720 - DWK changed yravln to yravln.total and added
           yravln.nh4, and yravln.no3 to resetYrFlux() and
           setMonth()
20030720 - DWK changed yrinnup to yrinnup.total and added
           yrinnup.nh4, and yrinnup.no3 to resetYrFlux() and
           setMonth()
20030720 - DWK changed yrnup to yrnup.total and added
           yrnup.nh4, and yrnup.no3 to resetYrFlux() and
           setMonth()
20030720 - DWK added yrrsoil, soil.yrlchDOM.carbon,
           atms.yrndep.total, atms.yrndep.nh4, atms.yrndep.no3,
           veg.yrnfix, microbe.yrnfix, microbe.yrgnmin,
           microbe.yrimmb, microbe.yrammnvol, microbe.yrnitrif,
           microbe.yrno3prd, microbe.yrnoprd, microbe.yrn2oprd,
           microbe.yrn2prd, microbe.yrdenitrif, soil.yrnh3flx,
           soil.yrnoflx, soil.yrn2oflx, soil.yrn2flx,
           soil.yrlchNO3 and soil.yrlchDOM.nitrogen to
           resetYrFlux() and setMonth()
20030720 - DWK deleted calculation of veg.nmax and added
           calculation of veg.nupnh4, veg.nupno3 and
           soil.lchNO3par to setELMNTecd()
20030720 - DWK modified "1 == baseline" statements in
           stepmonth() and stepyr()
200g.setTOTEC( (y[I_VEGC]
             + y[I_SOLC]
             + soil.getNSOLC()
             + y[I_DOC]) );
30722 - DWK added microbe.setNitrification() to cropDynamics()
           and natvegDynamics()
20030722 - DWK added conditional statements for y{I_LCHDOC],
           y[I_SNFIX], y[I_ANFIX], y[I_LCHDON], y[I_TNDEP],
           y[I_NH4DEP], y[I_NO3DEP], y[I_AMMN], y[I_NTRF],
           y[I_NO3P], y[I_NOP], y[I_N2OP], y[I_N2P], y[I_DNTRF],
           y[I_NH3FLX], y[I_NOFLX], y[I_N2OFLX], y[I_N2FLX],
           y[I_LCHNO3] and y[I_LCHDON] to massbal()
20030731 - DWK added I_ABIMMB to TTEM51()
20030731 - DWK added soil.abioticNimmob[] to delta(),
           setELMNTflux() and setMonth()
20030731 - DWK added soil.yrabNimmob[] to resetYrflux() and
           setMonth()
20030731 - DWK added y[I_ABIMMB] to massbal()
20030731 - DWK added soil.setDOMleaching() to cropDynamics() and
           natvegDynamics()
20030819 - DWK added veg.yrnfix to testEquilibrium()
20030826 - DWK added I_DOC, I_DON, I_DOCP and I_DONP to
           TTEM51()
20030826 - DWK added I_DON to boundcon()
20030826 - DWK added microbe.relntrf and
           microbe.DOMprod[].nitrogen to cropDynamics() and
           natvegDynamics()
20030826 - DWK changed calculation of soil.leachDOM[].nitrogen
           in cropDynamics() and natvegDynamics()
20030826 - DWK added I_DOC, I_DON, I_DOCP and I_DONP to delta()
           and massbal()
20030826 - DWK added I_DOC and I_DON to ECDsetElmntstate()
20030826 - DWK added soil.yrDOM and microbe.yrDOMprd to
           resetYrFlux() and setMonth()
20030826 - DWK added microbe.DOMprod to setELMNTflux() and
           setMonth()
20030902 - DWK added microbe.grossnmin[], microbe.immnh4[] and
           microbe.immno3[] to delta() and natvegDynamics()
20030902 - DWK added y[I_GMIN], y[I_NH4IMM] and y[I_NO3IMM]
           to TTEM51() and massbal()
20030902 - DWK added microbe.nh4imma[], microbe.nh4immb[],
           microbe.no3imma[] and microbe.no3immb[] to
           getcropecd() and getsitecd()
20031014 - BSF added eet/pet to veg.gppxio3() function call in
           cropDynamics() and natvegDynamics()
20031014 - BSF added half compounding to fprevozone in stepyr() and
           stepmonth()
20031015 - DWK added I_CO2G, I_CO2W, I_HCO3, I_RHCO3, I_ALK,
           I_ERDPOC, I_DSVCO2, I_LCHCO2, I_HCO3P,
           I_LCHHCO3, I_RHCO3P, I_LCHALK and
           I_NTCS, I_ISOPREN, I_TERPEN, I_ORVOC, I_OVOC,
           I_VOC and I_ERDPON to TTEM51()
20031016 - DWK changed missing values in ecdqc() from -999.9
           to -9999.9
20031020 - DWK modified conditions associated with "availn == 0"
           in cropDynamics() and natvegDynamics()
20031204 - DWK added veg.cmnt to function call 
           of microbe.setTraceGasProduction() in cropDynamics()
           and natvegDynamics()
20031204 - DWK added  microbe.tgmpar[dcmnt] to ecdqc(), 
           getcropecd() and getsitecd()
20031231 - DWK deleted initialization of veg.prvleafmx from 
           setELMNTecd()
20031231 - DWK added call to setELMNTecd() in monthlyTransient()                       
20040202 - DWK changed calculation of microbe.netnmin to be the
           difference of microbe.grossnmin[] - microbe.immnh4[]
           in cropDynamics() and natvegDynamics()
20040202 - DWK changed calculation of y[I_NMIN] to be the
           difference of y[I_GMIN] - y[I_NH4IMM] in massbal()
20040202 - DWK added microbe.immno3[] to calculation of
           pdstate[I_SOLN] in delta()
20040203 - DWK modified the calculation of microbe.relntrf in
           natvegDynamics() and cropDynamics()
20040203 - DWK modified the calculation of microbe.nitrif in
           natvegDynamics() and cropDynamics()
20040203 - DWK changed microbe.maxntrfCN[] to be 
           microbe.initntrf[] in getcropecd() and getsitecd()
20040203 - DWK changed microbe.minntrfCN[] to be 
           microbe.allntrf[] in getcropecd() and getsitecd()
20040206 - DWK changed calculation of microbe.relntrf in 
           natvegDynamics() and cropDynamics()
20040228 - DWK changed TTEM51:: to TTEM60::
20040228 - DWK changed double nh4imma[MAXCMNT] to 
           nh4imm1a[MAXCMNT] in functions
20040228 - DWK changed double nh4immb[MAXCMNT] to 
           nh4imm1b[MAXCMNT] in functions
20040228 - DWK added double nh4immcut[MAXCMNT], 
           double nh4imm2a[MAXCMNT] and double nh4imm2b[MAXCMNT]
           to functions
20040228 - DWK deleted double immno3[CYCLE], double yrimmno3, 
           double no3imma[MAXCMNT] and double no3immb[MAXCMNT]
           from functions()
20040707 - DWK deleted xtair and xtairflg from function call to
           soil.setMonthlySoilConditions() in getenviron()           
20040718 - DWK changed soil.leachNO3[] and soil.leachDOM[] to 
           be a function of (soil.rperc[] + soil.sperc[] in
           cropDynamics() and natvegDynamics()
20040919 - DWK added ntrfparcut[], ntrfpar1a[], ntrfpar1b[],
           ntrfpar2a[], ntrfpar2b[], DONpara[] and DONparb[] to
           getcropecd() and getsitecd()
20040919 - DWK added calculation of microbe.ntrfpar and 
           soil.DONpar to setELMNTecd()
20040922 - DWK added no3cut[], no31a[], no31b[], no32a[] and 
           no32b[] to ecdqc(), getcropecd() and getsitecd()
20040922 - DWK modified calculation of soil.avail[dm].no3 in 
           ECDsetELMNTstate()
20041003 - DWK added I_NTNS to TTEM60()
20041003 - DWK added ntns[] and yrntns to functions
20050409 - DWK added I_CDCMP and I_NDCMP to TTEM60()
20050409 - DWK added microbe.DOMprod[dm].carbon and 
           soil.erodePOM[dm].carbon to pdstate[I_SOLC] in 
           delta()
20050409 - DWK added pdstate[I_CDCMP], pdstate[I_DOC], 
           pdstate[I_DOCP], pdstate[I_LCHDOC], 
           pdstate[I_ERDPOC], pdstate[I_NTCS], pdstate[I_NDCMP],
           pdstate[I_ERDPON], and pdstate[I_NTNS] to delta()
20050409 - DWK added microbe.DOMprod[dm].nitrogen and 
           soil.erodePOM[dm].nitrogen to pdstate[I_SOLN] in 
           delta()
20050409 - DWK deleted microbe.DOMprod[dm].nitrogen from 
           pdstate[I_NH4] in delta()
20051117 - DWK added include ttem602.h and standard includes
20051117 - DWK changed Odeint4 to Odeint60 in TTEM60::TTEM60()
20051118 - DWK added functions adapt(), ask() and rfk() formerly
           in the RKF adaptive integrator
20051118 - DWK deleted inheritance of Odeint60
20051123 - DWK deleted stepyr()
20051123 - DWK changed stepyr() to stepmonth() in equilibrium()
           and yearlyTransient()
20051123 - DWK deleted setMonth()
20051123 - DWK added updateYearSummary()
20051123 - DWK deleted resetODEflux()
20051201 - DWK deleted deltaxclm()
20051202 - DWK deleted monthxclm()
20051207 - DWK deleted resetVegEcd()
20051208 - DWK added pcdisplayDT(), pcdisplayMonth() and 
           pcdisplayODEerr()
20051208 - DWK added getOptionalCflx(), getOptionalNflx(), 
           getOptionalSoilTemp(), getOptionalTraceGas() and 
           getOptionalWflx()
20051208 - DWK added public functions displayOptionalCflx(),
           displayOptionalNflx(), displayOptionalSoilTemp(),
           displayOptionalTraceGas() and displayOptionalWflx()
20060909 - DWK changed include from ttem602.h to ttem603.h
20060909 - DWK added soil.getThawedReactiveSOMProp() to 
           cropDynamics() and natvegDynamics()
20060909 - DWK added soil.getPROPREACTA() and 
           soil.getPROPREACTA() to ecdqc()
20060909 - DWK added soil.getPROPREACTA() and 
           soil.getPROPREACTA() to getcropecd() and getsitecd()
20061108 - DWK added private function urbanDynamics()
20061112 - DWK added private function pastureDynamics()
20070829 - DWK added veg.updateFoliage() to natvegDynamics(),
           cropDynamics(), pastureDynamics(), and urbanDynamics()        
20070830 - DWK changed include from ttem603b.h to ttem603c.h
20090127 - DWK changed include from ttem603c.h to ttem604.h

****************************************************************
************************************************************** */

#include<cstdio>

  using std::printf;

#include<iostream>

  using std::cout;
  using std::cin;
  using std::ios;
  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#ifdef BORLAND_CPP
  #include<stdlib>
#else 
  #include<cstdlib>
#endif

  using std::exit;
  using std::atof;
  using std::atoi;

#include<cmath>

  using std::exp;
  using std::fabs;
  using std::modf;
  using std::pow;
#include<math.h>
#include <float.h>
#include <limits>

#include<vector>

  using std::vector;
    
#include<string>
  
  using std::string;

// #define CALIBRATE_TEM
using namespace std;
#ifdef CALIBRATE_TEM
  #include<conio.h>
  #include<cctype>
    using std::toupper;
#endif

#include "ttem604.h"

/* **************************************************************
************************************************************** */

// Initialization of static members

int TTEM60::avlnflag = 0;
int TTEM60::nfeed = 0;
int TTEM60::rheqflag = 0;
int TTEM60::moistlim = 0;
int TTEM60::o3flag = 0;
int TTEM60::initbase = 0;
int TTEM60::baseline = 0;
int TTEM60::intflag = 0;

int TTEM60::maxnrun = 0;
int TTEM60::equil = 0;
int TTEM60::runsize = 0;
int TTEM60::maxyears = 0;
int TTEM60::strteq = 0;
int TTEM60::endeq = 0;
int TTEM60::modstartyr = 0;
int TTEM60::outputstartyr = 0;
int TTEM60::outputendyr = 0;
int TTEM60::diffyr = 0;
int TTEM60::wrtyr = 0;
const double mpe = 0.000001;

double TTEM60::ctol = 1.0;
double TTEM60::ntol = 0.02;
double TTEM60::wtol = 0.01;

// Initialization of adaptive integrator variables

double TTEM60::inittol = 0.01;
int TTEM60::maxit = 20;
long TTEM60::maxitmon = 2000;

double TTEM60::a1 =   0.115740741;

double   TTEM60::a3 =   0.548927875;
double  TTEM60::a31 =   0.09375;
double  TTEM60::a32 =   0.28125;

double   TTEM60::a4 =   0.535331384;
double  TTEM60::a41 =   0.879380974;
double  TTEM60::a42 =  -3.277196177;
double  TTEM60::a43 =   3.320892126;

double   TTEM60::a5 =  -0.20;
double  TTEM60::a51 =   2.032407407;
double  TTEM60::a52 =  -8.0;
double  TTEM60::a53 =   7.173489279;
double  TTEM60::a54 =  -0.2058966866;

double   TTEM60::b1 =   0.118518519;
double   TTEM60::b3 =   0.518986355;
double   TTEM60::b4 =   0.50613149;
double   TTEM60::b5 =  -0.18;
double   TTEM60::b6 =   0.036363636;
double  TTEM60::b61 =  -0.296296296;
double  TTEM60::b62 =   2.0;
double  TTEM60::b63 =  -1.381676413;
double  TTEM60::b64 =   0.45297271;
double  TTEM60::b65 =  -0.275;
double  TTEM60::tempvara1 =  1.00;

/* **************************************************************
************************************************************** */
/*
template<typename T>
bool is_infinite (const T &value)
{
	T max_value = std::numeric_limits<T>::max();
	T min_value = -max_value;
	return !(min_value <= value && value <= max_value);
}

template<typename T>
bool is_nan(const T &value)
{
	return value != value;
}

template<typename T>
bool is_valid(const T &value)
{
	return ! is_infinite(value) && ! is_nan(value);
}
*/

TTEM60::TTEM60() : predstr( NUMTEM )
{
  tol = inittol;
  syint = 1;  
  totyr = -99;
  rheqflag = 0;

#ifdef CALIBRATE_TEM
  soilfile = "tsoil43d.ecd";
  rootfile = "trotz43d.ecd";
  vegfile = "tveg604.ecd";
  leaffile = "tleaf43d.ecd";
  mcrvfile = "tmcrv604.ecd";
  agfile = "ag603.ecd";
  snowfile = "tspack50b2.ecd";
  slayerfile = "tslayer602.ecd";

//  soilfile = "sasssoil603.ecd";
//  rootfile = "sassrotz603.ecd";
//  vegfile = "sassveg604.ecd";
//  leaffile = "sassleaf603.ecd";
//  mcrvfile = "sassmcrv604.ecd";
//  agfile = "sassag603.ecd";
//  snowfile = "sassspack603.ecd";
//  slayerfile = "sassslayer603.ecd";

  rheqflag = 0;

  scy[0] = GET_LAI;
  scy[1] = GET_NPP;
  scy[2] = GET_CDCMP;
  scy[3] = GET_RH;
  scy[4] = GET_LCHDOC;
  scy[5] = GET_NTCB;

  sny[0] = GET_VNUP;
  sny[1] = GET_NDCMP;
  sny[2] = GET_LCHDON;

  swy[0] = GET_RAIN;
  swy[1] = GET_SNWINF;
  swy[2] = GET_VSM;
  swy[3] = GET_PET;
  swy[4] = GET_EET;

  ssty[0] = GET_ACTLAY;
  ssty[1] = GET_THAWPCT;

  sgy[0] = GET_NOFLX;
  sgy[1] = GET_N2OFLX;
  sgy[2] = GET_N2FLX;
#endif
  
// Identify potential output variables from TEM

// Ecosystem carbon pools determined by the integrator**********

  // vegetation carbon
  predstr.at( I_VEGC ) = "VEGC";  
  
  // reactive soil organic carbon      
  predstr.at( I_SOLC ) = "SOILORGC";   
  
  // dissolved organic carbon
  predstr.at( I_DOC ) = "DOC";

// Ecosystem nitrogen pools determined by the integrator********

  // vegetation structural nitrogen
  predstr.at( I_STRN ) = "VSTRUCTN";

  // vegetation labile nitrogen
  predstr.at( I_STON ) = "VSTOREN";

  // reactive soil organic nitrogen
  predstr.at( I_SOLN ) = "SOILORGN";  

  // dissolved organic nitrogen
  predstr.at( I_DON ) = "DON";

  // soil ammonium    
  predstr.at( I_NH4 ) = "SOILNH4";   
  
  // soil nitrate
  predstr.at( I_NO3 ) = "SOILNO3";   


// Ecosystem water pools determined by the integrator***********

  // soil moisture
  predstr.at( I_SM ) = "SMOIST";       

  // groundwater pool resulting from rainfall
  predstr.at( I_RGRW ) = "RGRNDH2O";

  // groundwater pool resulting from snow melt
  predstr.at( I_SGRW ) = "SGRNDH2O";

  //standing deadwood C and N
  predstr.at( I_DEADWOODC ) = "DEADWOODC";
  predstr.at( I_DEADWOODN ) = "DEADWOODN";
  predstr.at( I_FOLC ) = "FOLC";
  predstr.at( I_STEMC ) = "STEMC";
  predstr.at( I_CROOTC ) = "CROOTC";
  predstr.at( I_FROOTC ) = "FROOTC";
  predstr.at( I_CWD ) = "CWD";
  predstr.at( I_AGR ) = "AGR";
  predstr.at( I_AGL ) = "AGL";
  predstr.at( I_BGR ) = "BGR";
  predstr.at( I_BGL ) = "BGL";
  predstr.at( I_MDC ) = "MDC";
  predstr.at( I_AGE ) = "AGE";
  predstr.at( I_SOC ) = "SOC";

/// Phenology variables determined by the integrator************

  // un-normalized relative phenology variable
  predstr.at( I_UNRMLF ) = "UNRMLEAF";

  // normalized relative phenology variable (0 - 1.0)
  predstr.at( I_LEAF ) = "LEAF";

  // leaf area index
  predstr.at( I_LAI ) = "LAI";    
  
  // foliar projected cover
  predstr.at( I_FPC ) = "FPC";    


// Carbon fluxes for ecosystems ********************************

  // GPP not limited by nutrient availability
  predstr.at( I_INGPP ) = "VEGINGPP";

  // gross primary production
  predstr.at( I_GPP ) = "GPP";    

  // Direct ozone effects 
  predstr.at( I_FOZONE ) = "FOZONE";
  
  // Indirect ozone effects 
  predstr.at( I_FINDOZONE ) = "FINDOZON";

 // NPP not limited by nutrient availability
  predstr.at( I_INNPP ) = "VEGINNPP";

  // net primary production
  predstr.at( I_NPP ) = "NPP";      

  // gross plant respiration
  predstr.at( I_GPR ) = "GPR";      

  // vegetation maintenance respiration
  predstr.at( I_RVMNT ) = "RVMAINT";

  // vegetation growth respiration
  predstr.at( I_RVGRW ) = "RVGRWTH";

  // aboveground vegetation respiration
  predstr.at( I_ABVGPR ) = "ABVGPR";

  // root respiration
  predstr.at( I_ROOTGPR ) = "ROOTRESP";

  // litterfall carbon
  predstr.at( I_LTRC ) = "LTRC";   
  
  // carbon in decomposition of organic matter
  predstr.at( I_CDCMP ) = "CDECOMP";   

  // heterotrophic respiration
  predstr.at( I_RH ) = "RH";       

  // dissolved organic carbon production
  predstr.at( I_DOCP ) = "DOCPROD";

  // DOC leaching losses
  predstr.at( I_LCHDOC ) = "LEACHDOC";   

  // erosion of particulate organic carbon
  predstr.at( I_ERDPOC ) = "ERODEPOC";


// Nitrogen fluxes for ecosystems determined by the integrator

  // nitrogen fertilization
  predstr.at( I_AGFRTN ) = "AGFERTN";

  // Biological nitrogen fixation
  predstr.at( I_BNFIX ) = "BIOLNFIX";

  // Symbiotic nitrogen fixation
  predstr.at( I_SNFIX ) = "SYMNFIX";

  // Asymbiotic nitrogen fixation
  predstr.at( I_ANFIX ) = "ASYMNFIX";

  // VEGNUP not limited by carbon availability
  predstr.at( I_INNUP ) = "VEGINNUP";

  // VNH4UP not limited by carbon availability
  predstr.at( I_INNH4UP ) = "VINNH4UP";

  // VNO3UP not limited by carbon availability
  predstr.at( I_INNO3UP ) = "VINNO3UP";

  // nitrogen uptake by vegetation
  predstr.at( I_VNUP ) = "VEGNUP";

  // ammonium uptake by vegetation
  predstr.at( I_VNH4UP ) = "VEGNH4UP";

  // nitrate uptake by vegetation
  predstr.at( I_VNO3UP ) = "VEGNO3UP";

  // vegetation nitrogen uptake for structural components
  predstr.at( I_VSUP ) = "VEGSUP";

  // vegetation nitrogen uptake for labile components
  predstr.at( I_VLUP ) = "VEGLUP";

  // nitrogen mobilization by vegetation
  predstr.at( I_VNMBL ) = "VNMOBIL";

  // nitrogen resorption by vegetation
  predstr.at( I_VNRSRB ) = "VNRESORB";

  // litterfall nitrogen from vegetation
  predstr.at( I_LTRN ) = "LTRN";

  // nitrogen in decomposition of organic matter
  predstr.at( I_NDCMP ) = "CDECOMP";   

  // dissolved organic nitrogen production
  predstr.at( I_DONP ) = "DONPROD";

  // gross nitrogen mineralization
  predstr.at( I_GMIN ) = "GROSNMIN";

  // ammonium immobilization
  predstr.at( I_NH4IMM ) = "NH4IMMOB";

  // nitrate immobilization
//  predstr.at( I_NO3IMM ) = "NO3IMMOB";

  // total nitrogen immobilization
  predstr.at( I_NIMM ) = "NIMMOBIL";

  // net nitrogen mineralization
  predstr.at( I_NMIN ) = "NETNMIN";

  // abiotic NH4 immobilization
  predstr.at( I_AIMMNH4 ) = "AIMMNH4";

  // abiotic NO3 immobilization
  predstr.at( I_AIMMNO3 ) = "AIMMNO3";

  // ammonia volatilization
  predstr.at( I_AMMN ) = "AMMNVOL";

  // nitrification
  predstr.at( I_NTRF ) = "NITRIF";     
  
  // Nitrate production
  predstr.at( I_NO3P ) = "NO3PROD";    
  
  // Microbial nitric oxide production
  predstr.at( I_NOP ) = "NOPROD";      
  
  // Microbial nitrous oxide production
  predstr.at( I_N2OP ) = "N2OPROD";    
  
  // microbial dinitrogen production
  predstr.at( I_N2P ) = "N2PROD";      
  
  // Denitrification
  predstr.at( I_DNTRF ) = "DENITRIF";  

  // Soil ammonia flux
  predstr.at( I_NH3FLX ) = "NH3FLUX";
  
  // Soil nitric oxide flux
  predstr.at( I_NOFLX ) = "NOFLUX";
  
  // Soil nitrous oxide flux
  predstr.at( I_N2OFLX ) = "N2OFLUX";
  
  // Soil dinitrogen flux
  predstr.at( I_N2FLX ) = "N2FLUX";
  
  // Nitrate leaching losses
  predstr.at( I_LCHNO3 ) = "LEACHNO3";  
  
  // Dissolved organic nitrogen leaching losses
  predstr.at( I_LCHDON ) = "LEACHDON";  

  // Erosion of particulate organic nitrogen
  predstr.at( I_ERDPON ) = "ERODEPON";


// Water fluxes determined by the integrator********************

  // Irrigation
  predstr.at( I_AGIRRIG ) = "IRRIGATE"; 

  // Initial estimated evapotranspiration
  predstr.at( I_INEET ) = "INEET";

  // estimated evapotranspiration
  predstr.at( I_EET ) = "EET";

  // percolation of rainwater through soil profile
  predstr.at( I_RPERC ) = "RPERC";

  // percolation of snowmelt through soil profile
  predstr.at( I_SPERC ) = "SPERC";

  // runoff of rainwater
  predstr.at( I_RRUN ) = "RRUN";
  
  // runoff of snowmelt
  predstr.at( I_SRUN ) = "SRUN";        


// Other ecosystem carbon pools ********************************

  // nonreactive soil organic carbon      
  predstr.at(I_NSOLC ) = "NONSOLC";   

  // total soil organic carbon      
  predstr.at( I_TSOLC ) = "TOTSOLC";

  // total carbon pool found in ecosystem excluding products
  predstr.at( I_TOTEC ) = "TOTEC";

  // total carbon
  predstr.at( I_TOTC ) = "TOTALC";     

// Other ecosystem nitrogen pools ******************************

  // total nitrogen stored in vegetation
  predstr.at( I_VEGN ) = "VEGN";

  // nonreactive soil organic nitrogen
  predstr.at( I_NSOLN ) = "NONSOLN";

  // total soil organic nitrogen
  predstr.at( I_TSOLN ) = "TOTSOLN";

  // soil available nitrogen
  predstr.at( I_AVLN ) = "AVAILN";


// Other ecosystem water pools ******************************

  predstr.at( I_SNWPCK ) = "SNOWPACK";  // snowpack

  // available soil moisture
  predstr.at( I_AVLW ) = "AVAILW";

  // volumetric soil moisture
  predstr.at( I_VSM ) = "VSM";      

  // soil moisture expressed as percent total porosity
  predstr.at( I_PCTP ) = "PCTP";



// Other carbon fluxes for ecosystems **************************

  // soil respiration
  predstr.at( I_RSOIL ) = "RSOIL";  

  // net ecosystem production
  predstr.at( I_NEP ) = "NEP";     

  // net carbon exchange of ecosystem with atmosphere
  predstr.at( I_NCE ) = "NCE";

  // net terrestrial carbon balance
  predstr.at( I_NTCB ) = "NTCB";


// Other nitrogen fluxes for ecosystems ************************

  // total nitrogen inputs into ecosystem
  predstr.at( I_NINP ) = "NINPUT";

  // Total nitrogen losses from ecosystems
  predstr.at( I_NLST ) = "NLOST";

  // net terrestrial carbon balance
  predstr.at( I_NTNB ) = "NTNB";


// Other water fluxes ******************************************

  // potential evapotranspiration
  predstr.at( I_PET ) = "PET";

  // infiltration into the soil of water from snowmelt
  predstr.at( I_SNWINF ) = "SNOWINF";

  // water yield
  predstr.at( I_WYLD ) = "H2OYIELD";    


// Carbon stocks in products ***********************************

  // carbon in agricultural products
  predstr.at( I_AGPRDC ) = "AGPRODC";

  // carbon pool of products that decompose in 10 years
  predstr.at( I_PROD10C ) = "PROD10C";

  // carbon pool of products that decompose in 100 years
  predstr.at( I_PROD100C ) = "PROD100C";

  // carbon in all product pools
  predstr.at( I_TOTPRDC ) = "TOTPRODC";


// Carbon stocks in crop residue and stubble********************

  // carbon in crop residue
  predstr.at( I_RESIDC ) = "RESIDC";

  // stubble carbon
  predstr.at( I_AGSTUBC ) = "CRPSTUBC";   


// Nitrogen stocks in products *********************************

  // nitrogen in agricultural products
  predstr.at( I_AGPRDN ) = "AGPRODN";

  // nitrogen pool of products that decompose in 10 years
  predstr.at( I_PROD10N ) = "PROD10N";

  // nitrogen pool of products that decompose in 100 years
  predstr.at( I_PROD100N ) = "PROD100N";

  // nitrogen in all product pools
  predstr.at( I_TOTPRDN ) = "TOTPRODN";


// Nitrogen stocks in crop residue and stubble******************

  // nitrogen in crop residue
  predstr.at( I_RESIDN ) = "RESIDN";

  // stubble nitrogen
  predstr.at( I_AGSTUBN ) = "CRPSTUBN";


// Carbon fluxes associated with agricultural conversion *******

  // carbon loss from the ecosystem during conversion
  predstr.at( I_CNVRTC ) = "CONVERTC";

  // carbon loss from vegetation during conversion
  predstr.at( I_VCNVRTC ) = "VCONVRTC";

  // carbon loss from soils during conversion
  predstr.at( I_SCNVRTC ) = "SCONVRTC";

  // carbon associated with slash left after conversion
  predstr.at( I_SLASHC ) = "SLASHC";

  // carbon flux from ecosystem (NEP+CONVERTC)
  predstr.at( I_CFLX ) = "CFLUX";
  //carbon flux from burnt litter
  predstr.at( I_FFLC ) = "FFLC";
  //carbon flux from burnt soil duff (soil organic layer)
  predstr.at( I_FFDC ) = "FFDC";
  //carbon flux from burnt standing deadwood
  predstr.at( I_SWFC ) = "SWFC";
  //carbon fluxes from burnt coarse woody debris
  predstr.at( I_DWFC ) = "DWFC";
  predstr.at( I_AGE ) = "AGE";
  predstr.at( I_SOC ) = "SOC";

// Nitrogen fluxes associated with agricultural conversion *****

  // nitrogen loss from the ecosystem during conversion
  predstr.at( I_CNVRTN ) = "CONVERTN";

  // nitrogen loss from vegetation during conversion
  predstr.at( I_VCNVRTN ) = "VCONVRTN";

  // nitrogen loss from soils during conversion
  predstr.at( I_SCNVRTN ) = "SCONVRTN";

  // nitrogen associated with slash left after conversion
  predstr.at( I_SLASHN ) = "SLASHN";

  // Total organic N mineralized and retained in ecosystem
  //   after disturbance
  predstr.at( I_NRETNT ) = "NRETENT";

  // Vegetation N mineralized and retained in ecosystem
  //   after disturbance
  predstr.at( I_NVRTNT ) = "NVRETENT";

  // Soil organic N mineralized and retained in ecosystem
  //   after disturbance
  predstr.at( I_NSRTNT ) = "NSRETENT";


// Carbon and nitrogen fluxes to/from products *****************

  // carbon loss to formation of agricultural products
  predstr.at( I_AGFPRDC ) = "AGFPRODC";

  // nitrogen loss to formation of agricultural products
  predstr.at( I_AGPRDN ) = "AGFPRODN";

  // carbon loss to crop residue
  predstr.at( I_FRESIDC ) = "FRESIDC";

  // nitrogen loss to crop residue
  predstr.at( I_FRESIDN ) = "FRESIDN";

  // carbon loss to resulting from decomposition of agricultural
  //   products
  predstr.at( I_AGPRDFC ) = "AGPRODFC";

  // nitrogen loss resulting from decomposition of agricultural
  //   products
  predstr.at( I_AGPRDFN ) = "AGPRODFN";

  // carbon loss from crop residue
  predstr.at( I_RESIDFC ) = "RESIDFC";

  // nitrogen loss from crop residue
  predstr.at( I_RESIDFN ) = "RESIDFN";

  // carbon loss to formation of products that decompose in
  //  10 years
  predstr.at( I_PRDF10C ) = "PRDF10C";

  // nitrogen loss to formation of products that decompose in
  //   10 years
  predstr.at( I_PRDF10N ) = "PRDF10N";

  // carbon loss resulting from decomposition of PROD10C
  predstr.at( I_PRD10FC ) = "PRD10FC";

  // nitrogen loss resulting from decomposition of PROD10N
  predstr.at( I_PRD10FN ) = "PRD10FN";

  // carbon loss to formation of products that decompose in
  //  100 years
  predstr.at( I_PRDF100C ) = "PRDF100C";

  // nitrogen loss to formation of products that decompose in
  //   100 years
  predstr.at( I_PRDF100N ) = "PRDF100N";

  // carbon loss resulting from decomposition of PROD100C
  predstr.at( I_PRD100FC ) = "PRD100FC";

  // nitrogen loss resulting from decomposition of PROD100N
  predstr.at( I_PRD100FN ) = "PRD100FN";

  // carbon loss to the formation of all products
  predstr.at( I_TOTFPRDC ) = "TOTFPRDC";

  // nitrogen loss to the formation of all products
  predstr.at( I_TOTFPRDN ) = "TOTFPRDN";

  // carbon loss resulting from decomposition of all products
  predstr.at( I_TOTPRDFC ) = "TOTPRDFC";

  // nitrogen loss resulting from decomposition of all products
  predstr.at( I_TOTPRDFN ) = "TOTPRDFN";


// Agro-Ecosystem carbon and nitrogen pools *********************

  // crop carbon
  predstr.at( I_CROPC ) = "CROPC";       

  // carbon in natural vegetation
  predstr.at( I_NATVEGC ) = "NATVEGC";

  // crop nitrogen
  predstr.at( I_CROPN ) = "CROPN";       

  // nitrogen in natural vegetation
  predstr.at( I_NATVEGN ) = "NATVEGN";

  // crop structural N
  predstr.at( I_CSTRN ) = "CROPSTRN";    

  // structural N in natural vegetation
  predstr.at( I_NATSTRN ) = "NATSTRN";

  // crop labile N
  predstr.at( I_CSTON ) = "CROPSTON";    

  // labile N stored in natural vegetation
  predstr.at( I_NATSTON ) = "NATSTON";


// Crop phenology **********************************************

  // unnormalized leaf in crops
  predstr.at( I_CROPULF ) = "CRPUNMLF";

  // unnormalized leaf in natural vegetation
  predstr.at( I_NATULF ) = "NATUNMLF";

  // leaf of crops
  predstr.at ( I_CROPLEAF ) = "CROPLEAF";

  // leaf of natural vegetation
  predstr.at( I_NATLEAF ) = "NATLEAF";

  // leaf area index (LAI) of crops
  predstr.at( I_CROPLAI ) = "CROPLAI";

  // leaf area index (LAI) of natural vegetation
  predstr.at( I_NATLAI ) = "NATLAI";

  // foliar projected cover (FPC) of crops
  predstr.at( I_CROPFPC ) = "CROPFPC";

  // foliar projected cover (FPC) of natural vegetation
  predstr.at( I_NATFPC ) = "NATFPC";


// Additional carbon fluxes for agro-ecosystems *****************

  // GPP of crops not limited by nutrient availability
  predstr.at( I_AGINGPP ) = "CRPINGPP";

  // GPP of natural vegetation not limited by 
  //   nutrient availability
  predstr.at( I_NATINGPP ) = "NATINGPP";

  // gross primary production (GPP) of crops
  predstr.at( I_AGGPP ) = "CROPGPP";

  // gross primary production of natural vegetation
  predstr.at( I_NATGPP ) = "NATGPP";

  // NPP of crops not limited by nutrient availability
  predstr.at( I_AGINNPP ) = "CRPINNPP";

  // NPP of natural vegetation not limited by 
  //   nutrient availability
  predstr.at( I_NATINNPP ) = "NATINNPP";

  // net primary production (NPP) of crops
  predstr.at( I_AGNPP ) = "CROPNPP";

  // net primary production (NPP) of natural vegetation
  predstr.at( I_NATNPP ) = "NATNPP";

  // gross plant respiration of crops
  predstr.at( I_AGGPR ) = "CROPGPR";

  // gross plant respiration of natural vegetation
  predstr.at( I_NATGPR ) = "NATGPR";

  // maintenance respiration of crop plants
  predstr.at( I_AGRVMNT ) = "CRPRMNT";

  // maintenance respiration of natural vegetation
  predstr.at( I_NATRVMNT ) = "NATRVMNT";

  // growth respiration of crop plants
  predstr.at( I_AGRVGRW ) = "CRPRGRW";

  // growth respiration of natural vegetation
  predstr.at( I_NATRVGRW ) = "NATRGRW";

  // litterfall carbon from crops
  predstr.at( I_AGLTRC ) = "CROPLTRC";

  // litterfall carbon from natural vegetation
  predstr.at( I_NATLTRC ) = "NATLTRC";

  // Additional nitrogen fluxes for agro-ecosystems ************

  // symbiotic nitrogen fixation with crops
  predstr.at( I_AGSNFX ) = "CRPSNFIX";

  // symbiotic nitrogen fixation with natural vegetation
  predstr.at( I_NATSNFX ) = "NATSNFIX";

  // nitrogen uptake by crops not limited by carbon availability
  predstr.at( I_AGINNUP ) = "CRPINNUP";

  // nitrogen uptake by natural vegetation not limited by carbon
  //   availability
  predstr.at( I_NATINNUP ) = "NATINNUP";

  // ammonium uptake by crops not limited by carbon availability
  predstr.at( I_AINNH4UP ) = "CINNH4UP";

  // ammonium uptake vy natural vegetation not limited by carbon
  //  availability
  predstr.at( I_NINNH4UP ) = "NINNH4UP";

  // nitrate uptake by crops not limited by carbon availability
  predstr.at( I_AINNO3UP ) = "CINNO3UP";

  // nitrate uptake by natural vegetation not limited by carbon
  //   availability
  predstr.at( I_NINNO3UP ) = "NINNO3UP";

  // nitrogen uptake by crops
  predstr.at( I_AGVNUP ) = "CROPNUP";

  // nitrogen uptake by natural vegetation
  predstr.at( I_NATVNUP ) = "NATVNUP";

  // ammonium uptake by crops
  predstr.at( I_AGVNH4UP ) = "CRPNH4UP";

  // ammonium uptake by natural vegetation
  predstr.at( I_NVNH4UP ) = "NVNH4UP";

  // nitrate uptake by crops
  predstr.at( I_AGVNO3UP ) = "CRPNO3UP";

  // nitrate uptake by natural vegetation
  predstr.at( I_NVNO3UP ) = "NVNO3UP";

  // nitrogen uptake for structural components of crops
  predstr.at( I_AGVSUP ) = "CROPSUP";

  // nitrogen uptake for structural components of natural
  //  vegetation
  predstr.at( I_NATVSUP ) = "NATVSUP";

  // nitrogen uptake for labile components of crops
  predstr.at( I_AGVLUP ) = "CROPLUP";

  // nitrogen uptake for labile components of natural vegetation
  predstr.at( I_NATVLUP ) = "NATVLUP";

  // nitrogen mobilization by crops
  predstr.at( I_AGVNMBL ) = "CRPNMOBL";

  // nitrogen mobilization by natural vegetation
  predstr.at( I_NATVNMBL ) = "NATVNMBL";

  // nitrogen resorption by crops
  predstr.at( I_AGVNRSRB ) = "CRPNRSRB";

  // nitrogen resorption by natural vegetation
  predstr.at( I_NVNRSRB ) = "NVNRSRB";

  // litterfall nitrogen from crops
  predstr.at( I_AGLTRN ) = "CROPLTRN";

  //litterfall nitrogen from natural vegetation
  predstr.at( I_NATLTRN ) = "NATLTRN";

// Soil thermal variables **************************************

  // Soil temperature for 0 - 20 cm
  predstr.at( I_TSOIL ) = "TSOIL";
  
  // Soil temperature at the ground surface
  predstr.at( I_DST0 ) = "DST0";
  
  // Soil temperature at 5 cm depth
  predstr.at( I_DST5 ) = "DST5";
  
  // Soil temperature at 10 cm depth
  predstr.at( I_DST10 ) = "DST10";
  
  // Soil temperature at 20 cm depth
  predstr.at( I_DST20 ) = "DST20";
  
  // Soil temperature at 50 cm depth
  predstr.at( I_DST50 ) = "DST50";
  
  // Soil temperature at 100 cm depth
  predstr.at( I_DST100 ) = "DST100";
  
  // Soil temperature at 200 cm depth
  predstr.at( I_DST200 ) = "DST200";
  // Soil temperature at 300 cm depth
  predstr.at( I_DST300 ) = "DST300";
  
  predstr.at( I_FRONTD ) = "FRONTD";
  
  predstr.at( I_THAWBE ) = "THAWBE";
  
  predstr.at( I_THAWEND ) = "THAWEND";

// Percentage of month with thawed ground **********************

  predstr.at( I_THAWPCT ) = "THAWPCT";

// Percentage of month with thawed ground **********************

  predstr.at( I_ACTLAYER ) = "ACTLAYER";

  // gaseous carbon dioxide in soils
  predstr.at( I_CO2G ) = "SOILCO2G";

  // dissolved carbon dioxide in soils
  predstr.at( I_CO2W ) = "SOILCO2W";

  // dissolved bicarbonate in soils
  predstr.at( I_HCO3 ) = "SOILHCO3";

  // dissolved carbonate in soils
  predstr.at( I_RHCO3 ) = "SOLRHCO3";

  // alkalinity in soils
  predstr.at( I_ALK ) = "SOILALK";

  // dissolution of soil carbon dioxide
  predstr.at( I_CO2DISS ) = "CO2DISS";

  // leaching of dissolved carbon dioxide from soils
  predstr.at( I_LCHCO2 ) = "LCHCO2";

  // production of dissolved bicarbonate in soils
  predstr.at( I_HCO3P ) = "HCO3PRD";

  // production of dissolved carbonate in soils
  predstr.at( I_RHCO3P ) = "RHCO3PRD";

  // leaching of dissolved bicarbonate from soils
  predstr.at( I_LCHHCO3 ) = "LCHHCO3";

  // leaching of alkalinity from soils
  predstr.at( I_LCHALK ) = "LCHALK";

  // isoprene emmissions from vegetation
  predstr.at( I_ISOPREN ) = "ISOPRENE";

  // monoterpene emmissions from vegetation
  predstr.at( I_TERPEN ) = "TERPENE";

  // other reactive volatile organic carbon emmissions from
  //   vegetation
  predstr.at( I_ORVOC ) = "ORVOC";

  // other volatile organic carbon emmissions from
  //   vegetation
  predstr.at( I_OVOC ) = "OVOC";

  // volatile organic carbon emmissions from vegetation
  predstr.at( I_VOC ) = "VOC";


  dbugflg = 0;

};

/* **************************************************************
************************* Functions *****************************
************************************************************** */


/* *************************************************************
************************************************************* */

int TTEM60::adapt( const int& numeq, 
                   double pstate[], 
                   const double& ptol, 
                   const int& pdm )
{

  int i;
  double ipart;
  double fpart;
  double time = ZERO;
  double dt = 1.0;
  int mflag = 0;
  long nintmon = 0;
  double oldstate[NUMEQ];
  blackhol = 0;

  while( time != 1.0 )
  {
    test = REJECT;
    if( 1 == syint )
    {
    	while ( test != ACCEPT )
    	{
    		//if (veg.cmnt==4&& longitude ==-89.75 && latitude == 51.5) cout <<"soc0: "<<pstate[I_SOC]<<endl;
    		rkf( numeq, pstate, dt, pdm );
    		//if (veg.cmnt==4&& longitude ==-89.75 && latitude == 51.5) cout <<"soc3: "<<pstate[I_SOC]<<endl;
    		test = boundcon( dum4, error, ptol );
    		//cout <<"MDC0: "<<pstate[I_MDC]<<" ptol: "<<ptol<<" test: "<<test<<endl;
          	#ifdef CALIBRATE_TEM
    		if( test != ACCEPT )
    		{
    			// Display ODE errors to DOS screen

    			pcdisplayODEerr( test, pstate );
    		}
          	 #endif

    		if( dt <= pow( 0.5, maxit ) )
    		{
    			test = ACCEPT;
    			mflag = 1;
        
    			if( 0 == nintmon )
    			{
    				for( i = 0; i < numeq; ++i )
    				{
    					oldstate[i] = pstate[i];
    				}
    			}
	
    			++nintmon;
    		}//end of if dt
    		//cout <<"MDC1: "<<pstate[I_MDC]<<" test: "<<test<<endl;
    		if ( test == ACCEPT)
    		{
    			for( i = 0; i < numeq; ++i ) { pstate[i] = dum4[i]; }
        
    			time += dt;

          		#ifdef CALIBRATE_TEM
    			// Display time updates to the DOS screen

    			pcdisplayDT( time, dt );
          		#endif

    			fpart = modf( (0.01 + (time/(2.0*dt))), &ipart );
        
    			if( fpart < 0.1 && dt < 1.0 ) { dt *= 2.0; }
    		}
    		else { dt *= 0.5; }
    		//cout <<"MDC2: "<<pstate[I_MDC]<<" dt: "<<dt<<endl;
    		if( nintmon == maxitmon )
    		{
    			time = 1.0;
    			blackhol = 1;
        
    			for( i = 0; i < numeq; ++i )
    			{
    				pstate[i] = oldstate[i];
    			}
    		}
    		//cout <<"MDC3: "<<pstate[I_MDC]<<" nintmon: "<<nintmon<<" maxitmon: "<<maxitmon<<endl;

    	} //end of while test
    }    // end of if syint

  }      //end of while time
  return mflag;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM60::askODE( ofstream& rflog1 )
{


/* **************************************************************
	      Parameters for Adaptive Integrator
************************************************************** */

#ifdef PMODE
  *fgo >> inittol;
  *fgo >> maxit;
  *fgo >> maxitmon;

#else
  cout << endl << "Enter the proportional tolerance for the integrator: ";
  cin >> inittol;

  cout << "Enter the maximum number of iterations in the integrator: ";
  cin >> maxit;

  cout << "Enter the maximum number of times in a month that the" << endl;
  cout << "integrator can reach the maximum number of iterations: ";
  cin >> maxitmon;
#endif

  rflog1 << endl;
  rflog1 << "Enter the maximum number of iterations in the integrator: ";
  rflog1 << maxit << endl;

  rflog1 << endl;
  rflog1 << "Enter the proportional tolerance for the integrator: ";
  rflog1 << inittol << endl;

  rflog1 << endl;
  rflog1 << "Enter the maximum number of times in a month that the" << endl;
  rflog1 << "integrator can reach the maximum number of iterations: ";
  rflog1 << maxitmon << endl;

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

int TTEM60::boundcon( double ptstate[],
                      double err[],
                      const double& ptol )
{

  int test = ACCEPT;

// Check carbon and nitrogen state variables

  if( err[I_VEGC] > fabs( ptol * ptstate[I_VEGC] ) )
  {
    return test = temkey( I_VEGC )+1;
  }

  if( err[I_SOLC] > fabs( ptol * ptstate[I_SOLC] ) )
  {
    return test = temkey( I_SOLC )+1;
  }

  if( err[I_DOC] > fabs( ptol * ptstate[I_DOC] ) )
  {
    return test = temkey( I_DOC )+1;
  }

  /*
  if( 1 == nfeed
      && err[I_DEADWOODN] > fabs( ptol * ptstate[I_DEADWOODN ) )
  {
    return test = temkey( I_DEADWOODN )+1;
  }
  */

//  if( err[I_CO2G] > fabs( ptol * ptstate[I_CO2G] ) )
//  {
//    return test = temkey( I_CO2G )+1;
//  }

//  if( err[I_CO2W] > fabs( ptol * ptstate[I_CO2W] ) )
//  {
//    return test = temkey( I_CO2W )+1;
//  }

//  if( err[I_HCO3] > fabs( ptol * ptstate[I_HCO3] ) )
//  {
//    return test = temkey( I_HCO3 )+1;
//  }

//  if( err[I_RHCO3] > fabs( ptol * ptstate[I_RHCO3] ) )
//  {
//    return test = temkey( I_RHCO3 )+1;
//  }

  if( 1 == nfeed
      && err[I_STRN] > fabs( ptol * ptstate[I_STRN] ) )
  {
    return test = temkey( I_STRN )+1;
  }

  if( 1 == nfeed
      && err[I_STON] > fabs( ptol * ptstate[I_STON] ) )
  {
    return test = temkey( I_STON )+1;
  }

  if( 1 == nfeed
      && err[I_SOLN] > fabs( ptol * ptstate[I_SOLN] ) )
  {
    return test = temkey( I_SOLN )+1;
  }

  if( 1 == nfeed
      && err[I_DON] > fabs( ptol * ptstate[I_DON] ) )
  {
    return test = temkey( I_DON )+1;
  }

  if( 1 == nfeed
      && err[I_NH4] > fabs( ptol * ptstate[I_NH4] ) )
  {
    return test = temkey( I_NH4 )+1;
  }

  if( 1 == nfeed
      && err[I_NO3] > fabs( ptol * ptstate[I_NO3] ) )
  {
    return test = temkey( I_NO3 )+1;
  }
  if( err[I_DEADWOODC] > fabs( ptol * ptstate[I_DEADWOODC] ) )
  {
    return test = temkey( I_DEADWOODC )+1;
  }

  if( err[I_DEADWOODN] > fabs( ptol * ptstate[I_DEADWOODN] ) )
  {
    return test = temkey( I_DEADWOODN )+1;
  }
  if( err[I_GPP] > fabs( ptol * ptstate[I_GPP] ) )
  {
    return test = temkey( I_GPP )+1;
  }

  if( err[I_NPP] > fabs( ptol * ptstate[I_NPP] ) )
  {
    return test = temkey( I_NPP )+1;
  }

  if( err[I_RVMNT] > fabs( ptol * ptstate[I_RVMNT] ) )
  {
    //return test = temkey( I_RVMNT )+1;
  }

  if( 1 == nfeed
      && err[I_VNH4UP] > fabs( ptol * ptstate[I_VNH4UP] ) )
  {
    //return test = temkey( I_VNH4UP )+1;
  }

  if( 1 == nfeed
      && err[I_VNO3UP] > fabs( ptol * ptstate[I_VNO3UP] ) )
  {
    //return test = temkey( I_VNO3UP )+1;
  }

  if( 1 == nfeed
      && err[I_VSUP] > fabs( ptol * ptstate[I_VSUP] ) )
  {
    //return test = temkey( I_VSUP )+1;
  }

  if( 1 == nfeed
      && err[I_VNMBL] > fabs( ptol * ptstate[I_VNMBL] ) )
  {
    //return test = temkey( I_VNMBL )+1;
  }


  // Check water state variables

  if( err[I_SM]  > fabs( ptol * ptstate[I_SM] ) )
  {
    return test = temkey( I_SM )+1;
  }

  if( err[I_RGRW]  > fabs( ptol * ptstate[I_RGRW] ) )
  {
    //return test = temkey( I_RGRW )+1;
  }

  if( err[I_SGRW]  > fabs( ptol * ptstate[I_SGRW] ) )
  {
    //return test = temkey( I_SGRW )+1;
  }

  if( err[I_RPERC]  > fabs( ptol * ptstate[I_RPERC] ) )
  {
    //return test = temkey( I_RPERC )+1;
  }

  if( err[I_EET]  > fabs( ptol * ptstate[I_EET] ) )
  {
    //return test = temkey( I_EET )+1;
  }


  return test;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM60::cropDynamics( const int& pdm, double pstate[] )
{
  double avlNH4;
  double avlNO3;

//  int irrgflag = 1;

  double newpctp;
  double newsh2o;
  double newvsm;

  int perennial = 0;

  double propReact;
  
  soil.updateHydrology( elev,
                        atms.getTAIR(),
                        atms.getPREVTAIR(),
                        atms.getPREV2TAIR(),
                        atms.getRAIN(),
                        atms.getPET(),
                        pstate[I_SM],
                        pstate[I_RGRW],
                        pstate[I_SGRW],
                        ag.irrgflag,
                        ag.irrigate,
                        pdm );
                          
  
  // Assume no agricultural N fertilization (ag.fertn) and no
  //  nutrients resulting from agricultural conversion

  ag.fertn = ZERO;
  microbe.setAMMNVOL( ZERO );


  // Determine biological N fixation and partition it between
  //   symbiotic N fixation (veg.nfix) and nonsymbiotic N
  //   fixation (microbe.nfix)

  if( soil.getEET() > 0.1 )
  {
    bnfix = veg.getNFIXPARA( ag.cmnt ) * soil.getEET()
            + veg.getNFIXPARB( ag.cmnt );
  }
  else { bnfix = ZERO; }

  microbe.setNFIX( (bnfix * microbe.getNFIXPAR( ag.cmnt )) ) ;
  
  if( microbe.getNFIX() < ZERO ) 
  { 
    microbe.setNFIX( ZERO ); 
  }

  veg.setNFIX( (bnfix - microbe.getNFIX()) );

  if( veg.getNFIX() < ZERO ) { veg.setNFIX( ZERO ); }
  
  // Reduce EET if vegetation is not mature
  
  soil.setEVAPORATION( soil.getEET() 
                       * (1.0 - veg.getPROPTRANS( ag.cmnt)) );
					   
  //veg.updateFoliage( ag.cmnt, pstate[I_VEGC], soil.getEET() );
  veg.updateFoliage( ag.cmnt, pstate[I_FOLC], pstate[I_VEGC], soil.getEET() );
  soil.setEET( veg.getTRANSPIRATION() + soil.getEVAPORATION() );


  // Assume wetlands are wetter by the wfpsoff for determining
  //   moisture effects on vegetation and microbes

  newpctp = soil.getPCTP() + soil.getWFPSOFF();
  
  newsh2o = (newpctp * soil.getTOTPOR() * soil.getROOTZ()) 
		    / (100.0 * soil.getACTLAYER()); 

  newvsm = newsh2o / (soil.getROOTZ() * 1000.0);


  soil.setKH2O( newvsm, moistlim );

  // Get proportion of unfrozen organic matter in rooting zone
  
  propReact = soil.getThawedReactiveSOMProp( veg.cmnt );
  //if (veg.cmnt==7 && veg.getCURRENTVEG() == 50) cout <<" pdm0: "<<pdm<<" soiln: "<<" psnh4: "<<pstate[I_NH4]<<" psno3: "<<pstate[I_NO3]<<pstate[I_SOLN]<<endl;

  avlNH4 = pstate[I_NH4] * soil.getACTLAYER() / soil.getROOTZ();
  avlNO3 = pstate[I_NO3] * soil.getACTLAYER() / soil.getROOTZ();
  

  // Note: Microbes are assumed to be acting on "old" carbon
  //   (i.e. natural vegetation - veg.cmnt) rather than 
  //   "new" carbon associated with crops (i.e. ag.cmnt)

  microbe.updateDynamics( veg.cmnt,
                          soil.getPCTFLDCAP(),
                          soil.getACTLAYER(),
                          (pstate[I_SOLC] * propReact),
                          pstate[I_AGR],
                          pstate[I_AGL],
                          pstate[I_BGR],
                          pstate[I_BGL],
                          pstate[I_CWD],
                          (pstate[I_SOLN] * propReact),
                          newsh2o,
                          newvsm,
                          avlNH4,
                          (atms.getNH4DEP()+ ag.getNRETENT()),
                          moistlim,
                          ag.tillflag,
                          ag.getTILLFACTOR( veg.cmnt ),
                          soil.getKH2O() );
  //if (veg.cmnt==7 && veg.getCURRENTVEG() == 50) cout <<" avnh4_0.2: "<<avlNH4<< "avno3: "<<avlNO3<<" psnh4: "<<pstate[I_NH4]<<" psno3: "<<pstate[I_NO3]<<" don: "<<pstate[I_DON]<<" soiln: "<<pstate[I_SOLN]<<" vsm: "<<newvsm<<" h2o: "<<newsh2o<<endl;

   //if (pstate[I_NO3] != pstate[I_NO3]) exit(-1);
  if( ag.getGROWDD() >= GDDSEED )
  {
    if( 0 == moistlim )
    {
      if( ag.getCROPPRVPETMX() < atms.getPET() )
      {
        atms.setPRVPETMX( atms.getPET() );
      }
      else
      {
        atms.setPRVPETMX( ag.getCROPPRVPETMX() );
      }
    }
    else
    {
      if( ag.getCROPPRVEETMX() < soil.getEET() )
      {
        soil.setPRVEETMX( soil.getEET() );
      }
      else
      {
        soil.setPRVEETMX( ag.getCROPPRVEETMX() );
      }
    }    

    veg.updateDynamics( latitude,
    					longitude,
    					pdm,
    					ag.cmnt,
                        atms.getCO2(),
                        atms.getAOT40(),
						atms.getTAIR(),
                        atms.getNDEP(),
                        ag.getNRETENT(),
                        atms.getPAR(),
                        atms.getPET(),
                        atms.getPRVPETMX(),
                        soil.getEET(),
                        soil.getPRVEETMX(),
                        pstate[I_VEGC],
                        pstate[I_FOLC],
                        pstate[I_STEMC],
                        pstate[I_CROOTC],
                        pstate[I_FROOTC],
                        pstate[I_DEADWOODC],
                        pstate[I_STRN],
                        pstate[I_STON],
                        pstate[I_DEADWOODN],
                        newsh2o,
                        avlNH4,
                        avlNO3,
                        moistlim,
                        nfeed,
                        o3flag,
                        ag.state,
                        perennial,
                        ag.fertflag,
                        soil.getKH2O(),
                        microbe.getNETNMIN(),
                        microbe.getAMMNVOL(),
                        microbe.getNITRIF(),
                        microbe.getNO3PROD(),
                        ag.fertn,
                        pstate[I_AGE]);
  }                    
  else
  {
    // No crop plants exist - set all monthly fluxes to zero
    
    veg.resetMonthlyFluxes();
  }

  //if (veg.cmnt==7 && veg.getCURRENTVEG() == 50) cout <<" avnh4_0.5: "<<avlNH4<< "avno3: "<<avlNO3<<" psnh4: "<<pstate[I_NH4]<<" psno3: "<<pstate[I_NO3]<<" don: "<<pstate[I_DON]<<" soiln: "<<pstate[I_SOLN]<<endl;

  // Determine carbon and nitrogen leaching losses
  
  soil.updateLeachingLosses( veg.cmnt,
                             (pstate[I_DOC] * propReact), 
                             (pstate[I_DON] * propReact), 
                             avlNO3, 
                             (pstate[I_SM] * soil.getACTLAYER() / soil.getROOTZ()) );
  

  if ( soil.getLEACHDOC() > (pstate[I_DOC]  
       + microbe.getDOCPROD()) )
  {
    soil.setLEACHDOC( (pstate[I_DOC] + microbe.getDOCPROD()) );
  }

  if ( soil.getLEACHDON() > (pstate[I_DON] 
         + microbe.getDONPROD()) )
  {
    soil.setLEACHDON( (pstate[I_DON] + microbe.getDONPROD()) );
  }

  
  // Determine loss of POC through erosion
  
  soil.setERODEPOC( ZERO );


  // Determine loss of PON through erosion
  
  soil.setERODEPON( ZERO );
  
  
  // Determine trace gas production based on microbe.no3prod,
  //   ag.fertn and water-filled pore space

  microbe.setTraceGasProduction( veg.cmnt, 
                                 newpctp, 
                                 ag.fertn );
  //if (veg.cmnt==7 && veg.getCURRENTVEG() == 50) cout <<" avnh4_1: "<<avlNH4<< "avno3: "<<avlNO3<<" leachdon: "<<soil.getLEACHDON()<<" leachno3: "<<soil.getLEACHNO3()<<endl;


  // Limit ecosystem N losses to total NO3 available

  if( soil.getLCHNO3PAR() > ZERO )
  {    
    if ( soil.getLEACHNO3() > (avlNO3 + atms.getNO3DEP()
         + microbe.getNO3PROD() + ag.fertn - microbe.getDENITRIF() 
         - veg.getNO3UPTAKE()) )
    {
      if ( (avlNO3 + atms.getNO3DEP() + microbe.getNO3PROD() 
           + ag.fertn - microbe.getDENITRIF() 
           - veg.getNO3UPTAKE()) > ZERO )
      {    
        soil.setLEACHNO3( (avlNO3
                           + atms.getNO3DEP()
                           + microbe.getNO3PROD()
                           + ag.fertn
                           - microbe.getDENITRIF()
                           - veg.getNO3UPTAKE()) );
      }
      else
      {
        soil.setLEACHNO3( ZERO );

        if( (avlNO3 + atms.getNO3DEP() + microbe.getNO3PROD() 
            + ag.fertn - veg.getNO3UPTAKE()) > ZERO )
        {
          microbe.setDENITRIF( (avlNO3
                               + atms.getNO3DEP()
                               + microbe.getNO3PROD()
                               + ag.fertn
                               - veg.getNO3UPTAKE()) );
        }
        else
        {
          microbe.setDENITRIF( ZERO );
          microbe.setN2PROD( ZERO );
          microbe.setN2OPROD( ZERO );

          veg.setNO3UPTAKE( (avlNO3
                             + atms.getNO3DEP()
                             + microbe.getNO3PROD()
                             + ag.fertn) );
        }
      }
    }
  }
  else     
  {
    // If no leaching losses, limit ecosystem N losses to 
    //   total NO3 available for denitrification
    
    if( microbe.getDENITRIF() > (avlNO3 + atms.getNO3DEP() 
         + ag.fertn + microbe.getNO3PROD() - veg.getNO3UPTAKE()) )
    {

      if( (avlNO3 + atms.getNO3DEP() + ag.fertn 
          + microbe.getNO3PROD() - veg.getNO3UPTAKE()) > ZERO )
      {
        microbe.setDENITRIF( (avlNO3
                             + atms.getNO3DEP()
                             + ag.fertn
                             + microbe.getNO3PROD()
                             - veg.getNO3UPTAKE()) );
      }
      else
      {
        microbe.setDENITRIF( ZERO );
        microbe.setN2PROD( ZERO );
        microbe.setN2OPROD( ZERO );

        veg.setNO3UPTAKE( (avlNO3
                           + atms.getNO3DEP()
                           + ag.fertn
                           + microbe.getNO3PROD()) );
      }

                                 
      if( microbe.getN2OPROD() > ZERO )
      {
        microbe.setN2PROD( microbe.getDENITRIF() - microbe.getN2OPROD() );
      
        if( microbe.getN2PROD() < ZERO )
        {
          microbe.setN2PROD( ZERO );
          microbe.setN2OPROD( microbe.getDENITRIF() );
        }
      }
      else { microbe.setN2PROD( microbe.getDENITRIF() ); }           
    }    
  }
 //if (veg.cmnt==7 && veg.getCURRENTVEG() == 50) cout <<" avnh4_2: "<<avlNH4<< "avno3: "<<avlNO3<<" leachdon: "<<soil.getLEACHDON()<<" leachno3: "<<soil.getLEACHNO3()<<endl;

 if (avlNO3!=avlNO3) exit(-1);
  // Determine trace gas fluxes based on ammonia volatilization
  //   (microbe.ammnvol), NO production (microbe.noprod),
  //   N2O production (microbe.n2oprod) and N2 production
  //   (microbe.n2prod)

  soil.setTraceGasFluxes( microbe.getAMMNVOL(),
                          microbe.getNOPROD(),
                          microbe.getN2OPROD(),
                          microbe.getN2PROD() );


  if( 0 == avlnflag )
  {
    // Do not allow changes in available nitrogen in soil
    //   (used during calibration process)

    // Equilibrate Ammonium (NH4) pool

    microbe.setAMMNVOL( (atms.getNH4DEP()
                        + microbe.getGROSSNMIN()
                        - microbe.getIMMNH4()
                        - veg.getNH4UPTAKE()
                        - microbe.getNITRIF()) );

    soil.setNH3FLUX( microbe.getAMMNVOL() );
    
    
    // Equilibrate nitrate (NO3) pool

    soil.setLEACHNO3( (atms.getNO3DEP()
                      + microbe.getNO3PROD()
                      - microbe.getDENITRIF()
                      - veg.getNO3UPTAKE()) );

    // Equilibrate DON pool

    soil.setLEACHDON( microbe.getDONPROD() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM60::delta( const int& pdm,
                    double pstate[],
                    double pdstate[] )
{
//  double tmprootz;

  //if (veg.cmnt == 7 && longitude ==-93.5 && latitude == 41.25 && veg.getCURRENTVEG()==50) cout <<"stn1: "<<pstate[I_STRN]<<endl;

  switch( ag.state )
  {
    case 1:  cropDynamics( pdm, pstate ); break;
    case 2:  pastureDynamics( pdm, pstate ); break;
    case 3:  urbanDynamics( pdm, pstate ); break;
    default: natvegDynamics( pdm, pstate );
  }
  //if (veg.cmnt == 7 && longitude ==-93.5 && latitude == 41.25 && veg.getCURRENTVEG()==50) cout <<"stn2: "<<pstate[I_STRN]<<endl;
  //if (longitude ==-127.75 && latitude == 55.25 && ag.forestage==0 && pdm == 8) cout <<" mon: "<<pdm<<" forestage: "<<ag.forestage<<" gpp: "<<veg.getGPP()<<" disturbflag: "<<disturbflag<<endl;

  // Describe monthly changes to carbon pools and fluxes for ODE 
  //   state variables (i.e., pdstate)

  double leaffallrate =veg.getLCFALL(veg.cmnt);
  double stemfallrate=veg.getSCFALL(veg.cmnt);
  double crootfallrate=veg.getCCFALL(veg.cmnt);
  double frootfallrate=veg.getFCFALL(veg.cmnt);

  if (veg.getifwoody(veg.cmnt) == 1)
  {
	  leaffallrate = ( 1+veg.getMORTCOEFA(veg.cmnt) * exp(veg.getMORTCOEFB(veg.cmnt)*pstate[I_AGE]) ) * veg.getLCFALL(veg.cmnt);
	  stemfallrate = ( 1+veg.getMORTCOEFA(veg.cmnt) * exp(veg.getMORTCOEFB(veg.cmnt)*pstate[I_AGE]) ) * veg.getSCFALL(veg.cmnt);
	  crootfallrate = ( 1+veg.getMORTCOEFA(veg.cmnt) * exp(veg.getMORTCOEFB(veg.cmnt)*pstate[I_AGE]) ) * veg.getCCFALL(veg.cmnt);
	  frootfallrate = ( 1+veg.getMORTCOEFA(veg.cmnt) * exp(veg.getMORTCOEFB(veg.cmnt)*pstate[I_AGE]) ) * veg.getFCFALL(veg.cmnt);
  }
  double leafltc = leaffallrate * pstate[I_FOLC];
  double stemltc = stemfallrate * pstate[I_STEMC];
  double crootltc = crootfallrate * pstate[I_CROOTC];
  double frootltc = frootfallrate * pstate[I_FROOTC];

  double totltfallc=leafltc+stemltc+crootltc+frootltc;
  double totltfalln=0.0;
  totltfalln=totltfallc*pstate[I_STRN]/omitzero(pstate[I_VEGC])*veg.getNFALL(veg.cmnt)/omitzero(veg.getCFALL(veg.cmnt));


  // Carbon pools in ecosystems

  pdstate[I_VEGC] = veg.getGPP() 
                    - veg.getGPR() 
                    - totltfallc;
                    //- veg.getLTRFALC();

  //added to deal with when totNPP <0.0 and leafc, rootc & stemc are less than respective NPP

  pdstate[I_FOLC] = veg.getNPP_leaf()
                    - leafltc;
  pdstate[I_STEMC] = veg.getNPP_stem()
                    - stemltc;
  pdstate[I_CROOTC] = veg.getNPP_croot()
                    - crootltc;
  pdstate[I_FROOTC] = veg.getNPP_froot()
                    - frootltc;

  //if (latitude==39.25 && longitude==-121.75) cout <<" pdm: "<<pdm<< "leafc: "<<pstate[I_FOLC]<< " pd_folc: "<<pdstate[I_FOLC]<<" leafnpp: "<<veg.getNPP_leaf()<<" leafltc: "<<leafltc<<" gpp: "<<veg.getGPP()<<" gpr: "<<veg.getGPR()<<" npp: "<<veg.getNPP()<<endl;
  double pdalloc=pdstate[I_FOLC]+pdstate[I_STEMC]+pdstate[I_CROOTC]+pdstate[I_FROOTC];
  double allocltrc = leafltc+stemltc+crootltc+frootltc;


  double fallwoodc=veg.getwoodfall(veg.cmnt) * pstate[I_DEADWOODC];
  double fallwoodn=veg.getwoodfall(veg.cmnt) * pstate[I_DEADWOODN];

  if (fallwoodc<0.0001) fallwoodc=0.0;
  if (fallwoodn<0.0001) fallwoodn=0.0;

  pdstate[I_SOLC] = //veg.getLTRFALC()
		  	  	  	totltfallc
                    //+ ag.getSLASHC()
                    + fallwoodc
					//+ ag.getdeadslashc()
                    //- ag.getSCONVRTFLXC()
                    - microbe.getRH()
                    - microbe.getDOCPROD()
                    - soil.getERODEPOC();
  //pdstate[I_SOC]  = microbe.getSOIL_LITTER()
		  	  	  	  //-microbe.getRH_SOIL();
  //if (veg.cmnt==4) cout <<"soilc: "<<pstate[I_SOC]<<" pdsoilc: "<<pdstate[I_SOC]<<" soillitter: "<<microbe.getSOIL_LITTER()<<" rh_soil: "<<microbe.getRH_SOIL()<<endl;
  //if (longitude == -127.75 && latitude ==55.25 && disturbflag==3 && pdm==8) cout <<" ysoilc: "<<pstate[I_SOLC]<<" pdsoilc: "<<pdstate[I_SOLC]<<" getltrc: "<<veg.getLTRFALC()<<" slashc: "<<ag.getSLASHC()<<" fallwoodc: "<<fallwoodc<<" woodc: "<<pstate[I_DEADWOODC]
          //<<" sconvertC: "<<ag.getSCONVRTFLXC()<<" rh: "<<microbe.getRH()<< endl;

  if (veg.getifwoody(veg.cmnt)==1) //forest
  {
  pdstate[I_CWD] = stemltc * 0.7 //coarse woody debris
		  	  	   //+ag.getstem2slash() * 0.7 //disturbed stem cwd
		  	  	   //+ag.getdeadslashc() * 0.8
		  	  	   +fallwoodc*0.8
                   -microbe.getRH_CWD();
                   //-ag.getdwfc();
  pdstate[I_AGR] = stemltc * 0.2 //aboveground resistant litter
                   //+veg.getCFALL(veg.cmnt) * pstate[I_FOLC] * 0.1
		  	  	   //+ ag.getstem2slash()*0.3//disturbed stem litter
		  	  	   //+ ag.getdeadslashc()*0.2
		  	  	   + fallwoodc*0.2
		  	  	   -microbe.getRH_AGR();
		  	  	   //-pstate[I_AGR]* ag.getfflc() / omitzero(pstate[I_AGR] + pstate[I_AGL]) ;
  pdstate[I_AGL] = leafltc //aboveground labile litter
                   + stemltc * 0.1
                   //+ag.getleaf2slash()//disturbed leaf litter
                   -microbe.getRH_AGL();
                   //-ag.getfflc() * pstate[I_AGL] / omitzero(pstate[I_AGR] + pstate[I_AGL]) ;
  pdstate[I_BGR] = crootltc * 0.8 //belowground resistant litter
                   //+ag.getroot2slash()*0.8 //disturbed root
		  	  	   -microbe.getRH_BGR();
  pdstate[I_BGL] = frootltc //belowground labile litter
                   //+ag.getroot2slash()*0.2 //disturbed root
                   +crootltc* 0.2
                   -microbe.getRH_BGL();
  }

  if (veg.getifwoody(veg.cmnt)==2) //shrub
  {
	  pdstate[I_CWD] = stemltc* 0.2 //coarse woody debris for woodland and shrubland
                   //+ ag.getstem2slash() * 0.5//disturbed stem litter
                   +fallwoodc*0.8
			  	   -microbe.getRH_CWD();
                   //-ag.getdwfc();
	  pdstate[I_AGR] = stemltc * 0.6 //aboveground resistant litter
			       //+ ag.getstem2slash()*0.5//disturbed stem litter
			       +fallwoodc*0.2
			  	   -microbe.getRH_AGR();
			       //-ag.getfflc() * pstate[I_AGR] / omitzero(pstate[I_AGR] + pstate[I_AGL]);
	  pdstate[I_AGL] = leafltc //aboveground labile litter
	                   + stemltc * 0.2
	                   //+ag.getleaf2slash()//disturbed leaf litter
	                   -microbe.getRH_AGL();
	                   //-ag.getfflc() * pstate[I_AGL] / omitzero(pstate[I_AGR] + pstate[I_AGL]) ;
	  pdstate[I_BGR] = crootltc * 0.8 //belowground resistant litter
	                   //+ag.getroot2slash()*0.8 //disturbed root
			  	  	   -microbe.getRH_BGR();
	  pdstate[I_BGL] = frootltc //belowground labile litter
	                   //+ag.getroot2slash()*0.2 //disturbed root
	                   +crootltc* 0.2
	                   -microbe.getRH_BGL();
  }
  if (veg.getifwoody(veg.cmnt)==0) //herbaceous
  {
	  pdstate[I_CWD] = stemltc* 0.0 //no coarse woody debris for herbaceous plants
                   //+ ag.getstem2slash() * 0.0//disturbed stem litter
                   +fallwoodc*0.8
			  	   -microbe.getRH_CWD();
                   //-ag.getdwfc();
	  pdstate[I_AGR] = stemltc * 0.5 //aboveground resistant litter
			       //+ ag.getstem2slash()*1.0//disturbed stem litter
			       +fallwoodc*0.2
			  	   -microbe.getRH_AGR();
			       //-ag.getfflc() * pstate[I_AGR] / omitzero(pstate[I_AGR] + pstate[I_AGL]);
	  pdstate[I_AGL] = leafltc //aboveground labile litter
	                   + stemltc * 0.5
	                   //+ag.getleaf2slash()//disturbed leaf litter
	                   -microbe.getRH_AGL();
	                   //-ag.getfflc() * pstate[I_AGL] / omitzero(pstate[I_AGR] + pstate[I_AGL]) ;
	  pdstate[I_BGR] = crootltc * 0.6 //belowground resistant litter
		  	  	   //+ag.getroot2slash()*0.3
                   -microbe.getRH_BGR();
	  pdstate[I_BGL] = frootltc //belowground labile litter
                   +crootltc* 0.4
                   //+ag.getroot2slash()*0.7
                   -microbe.getRH_BGL();
  }

  //if (latitude==39.25&& longitude==-121.75) cout <<"pdm: "<<pdm<<" agl: "<<pstate[I_AGL]<<" agr: "<<pstate[I_AGR]<<" fflc: "<<ag.getfflc()<<" DWFC: "<<ag.getdwfc() <<" swfc: "<<ag.getswfc()<<" vconvert: "<<ag.getVCONVRTFLXC()<<" slashc: "<<ag.getSLASHC()<<" rootmort: "<<ag.rootmortpar<<" leafmort: "<<ag.leafmortpar<<" leafc: "<<pstate[I_FOLC]<<" stemc: "<<pstate[I_STEMC]<<" rootc: "<<pstate[I_CROOTC] + pstate[I_FROOTC]<<endl;
  pdstate[I_DOC] = microbe.getDOCPROD()
                   - soil.getLEACHDOC();

  pdstate[I_DEADWOODC] = //ag.getdeadwoodc()
						  - fallwoodc;

  //cout <<"mon: "<<pdm<<" cmnt: "<<veg.cmnt<<" getfall: "<<veg.getwoodfall(veg.cmnt)<<" fallwoodc: "<<fallwoodc<<" pd_deadwoodc: " <<pdstate[I_DEADWOODC]<<" deadwoodc: "<<pstate[I_DEADWOODC]<<endl;
  //if (veg.cmnt==7 && y[I_DEADWOODC]>0.0) cout <<"pdm: "<<pdm<<" deadwoodc: "<<y[I_DEADWOODC]<<" getdeadc: "<<ag.getdeadwoodc()<<" fallc: "<<fallwoodc<<" deadslash: "<<ag.getdeadslashc()<<" swfc: "<<ag.getswfc()<<endl;
  pdstate[I_DEADWOODN] = //ag.getdeadwoodn()
						 - fallwoodn;
		 	 	 	 	 //- ag.getswfc() * pstate[I_DEADWOODN]/omitzero(pstate[I_DEADWOODC])
		 	 	 	 	 //- ag.getdeadslashn();

  if (pdstate[I_DEADWOODC]<0.0001 && pdstate[I_DEADWOODC]>-0.0001) pdstate[I_DEADWOODC]=0.0;
  if (pdstate[I_DEADWOODN]<0.0001 && pdstate[I_DEADWOODN]>-0.0001) pdstate[I_DEADWOODN]=0.0;


  //if (disturbflag >1) cout <<"cmnt: " <<veg.cmnt<<" pdm: " <<pdm<<" disturbflag: " <<disturbflag<<" disturbmon: "<<disturbmonth <<" totaldeadc: " <<y[I_DEADWOODC] <<" deadwoodc: "<<ag.getdeadwoodc()<<" deadltc: "<<veg.getdeadwoodltc()<<" slashc: " <<ag.getSLASHC()<< endl;

//  pdstate[I_CO2G] = microbe.getRH()
//                    + veg.getROOTRESP()
//                    - rsoil
//                    - soil.getDISSCO2();

//  pdstate[I_CO2W] = soil.getDISSCO2()
//                    - soil.getLEACHCO2()
//                    - soil.getFORNHCO3();

//  pdstate[I_HCO3] = soil.getFORMHCO3()
//                    + soil.getFORMRHCO3()
//                    - soil.getLEACHHCO3;

//  pdstate[I_RHCO3] = -1.0 * soil.getFORMRHCO3();


  // Nitrogen pools in ecosystems

  pdstate[I_STRN] = veg.getSUPTAKE()
		  	  	  	- totltfalln
                    //- veg.getLTRFALN()
                    - veg.getNRESORB()
                    + veg.getNMOBIL()
                    + veg.getNFIX();

  pdstate[I_STON] = veg.getLUPTAKE()
                    + veg.getNRESORB() 
                    - veg.getNMOBIL();
  //if (veg.cmnt == 7) cout <<"ltrfn1: "<< veg.getLTRFALN()<<" slashn: "<<ag.getSLASHN()<<" falln: "<<fallwoodn<<" deadslashn: "<< " don: "<<microbe.getDONPROD()<<endl;

  pdstate[I_SOLN] =
		  	  	  	totltfalln
                    //+ ag.getSLASHN()
                    + fallwoodn
                    //+ ag.getdeadslashn()
                    //- ag.getSCONVRTFLXN()
                    //- ag.getNSRETENT()
                    - microbe.getNETNMIN()
                    - microbe.getDONPROD()
                    - soil.getERODEPON()
                    + microbe.getNFIX();
  //if (veg.cmnt == 7 && longitude ==-93.5 && latitude == 41.25 && veg.getCURRENTVEG()==50) cout <<"stn3: "<<pstate[I_STRN]<<endl;
  pdstate[I_DON] = microbe.getDONPROD()
                   - soil.getLEACHDON();
  //if (veg.cmnt == 7 && veg.getCURRENTVEG()==50) cout <<"mon1: "<<pdm<<" pdNo3: "<<pdstate[I_NO3] <<" psNO3: "<<pstate[I_NO3]<<" pdNH4: "<<pdstate[I_NH4] <<" psnh4: "<<pstate[I_NH4]<<" donprod: "<< microbe.getDONPROD()<<" leachdon: "<<soil.getLEACHDON()<<endl;
  pdstate[I_NH4] = atms.getNH4DEP()
				   + ag.getNSRETENT()
                   + ag.getFIRENDEP()
                   + microbe.getGROSSNMIN()
                   - microbe.getIMMNH4()
                   - veg.getNH4UPTAKE()
                   - microbe.getAMMNVOL()
                   - microbe.getNITRIF();
                   //+ microbe.getNFIX() //add microbial Nfix to available N rather than SOLN. modified by cgs2014


  pdstate[I_NO3] = atms.getNO3DEP()
                   + ag.fertn
                   + microbe.getNO3PROD()
                   - microbe.getDENITRIF()
                   - veg.getNO3UPTAKE()
                   - soil.getLEACHNO3();
  //if (veg.cmnt == 7 && veg.getCURRENTVEG()==50) cout <<"mon2: "<<pdm<<" pdNo3: "<<pdstate[I_NO3] <<" psNO3: "<<pstate[I_NO3]<<" pdNH4: "<<pdstate[I_NH4] <<" psnh4: "<<pstate[I_NH4]<<" ltrfn2: "<< veg.getLTRFALN()<<" nmin: "<<microbe.getGROSSNMIN()<<" immnh4: "<<microbe.getIMMNH4()<<" nh4up: "<< veg.getNH4UPTAKE()<<" ammn: "<<microbe.getAMMNVOL()<< " nitrif: "<<microbe.getNITRIF()<<endl;
  //if (pstate[I_NO3]!=pstate[I_NO3]) exit(-1);

  // Water pools

  pdstate[I_SM] = soil.getSNOWINF()
                    + atms.getRAIN()
                    + ag.irrigate
                    - soil.getEET()
                    - soil.getRPERC()
                    - soil.getSPERC();
 // if (longitude == -115.75 && latitude ==47.5 && veg.cmnt==5) cout <<" month: "<<pdm<<"pd_ISM: "<<pdstate[I_SM]<<" y_ism: "<<pstate[I_SM]<<" rain: "<<atms.getRAIN()<<" eet: "<<soil.getEET()<<" rperc: "<<soil.getRPERC()<<" sperc: "<<soil.getSPERC()<<" pet: "<<atms.getPET()<<endl;
  pdstate[I_RGRW] = soil.getRPERC() - soil.getRRUN();

  if (pdstate[I_RGRW] <0.0000001 && pdstate[I_RGRW] >-0.0000001 ) pdstate[I_RGRW] = 0.0;

  pdstate[I_SGRW] = soil.getSPERC() - soil.getSRUN();
  if (pdstate[I_SGRW] <0.0000001 && pdstate[I_SGRW]>-0.0000001) pdstate[I_SGRW] = 0.0;
  // Phenology
  //cout <<"I_RGRW: "<<pdstate[I_RGRW] <<" rperc: " << soil.getRPERC()<< " sperc: " <<soil.getSPERC()<<endl;
  pdstate[I_UNRMLF] = veg.getUNNORMLEAF();

  pdstate[I_LEAF] = veg.getLEAF();

  pdstate[I_LAI] = veg.getLAI();

  pdstate[I_FPC] = veg.getFPC();

  // Carbon fluxes in ecosystems

  pdstate[I_INGPP] = veg.getINGPP();

  pdstate[I_GPP] = veg.getGPP();

  pdstate[I_FOZONE] = veg.getFOZONE();

  pdstate[I_FINDOZONE] = veg.getFINDOZONE();

  pdstate[I_INNPP] = veg.getINNPP();

  pdstate[I_NPP] = veg.getNPP();

  pdstate[I_GPR] = veg.getGPR();

  pdstate[I_RVMNT] = veg.getRMAINT();

  pdstate[I_RVGRW] = veg.getRGRWTH();

  pdstate[I_ABVGPR] = veg.getABVGPR();

  pdstate[I_ROOTGPR] = veg.getROOTRESP();

//  pdstate[I_ISOPREN] = veg.voc.isoprene;

//  pdstate[I_TERPEN] = veg.voc.monoterpene;

//  pdstate[I_ORVOC] = veg.voc.otherReactive;

//  pdstate[I_OVOC] = veg.voc.other;

  //pdstate[I_LTRC] = veg.getLTRFALC();
  pdstate[I_LTRC] = totltfallc;

  pdstate[I_CDCMP] = microbe.getDECOMP();
  
  pdstate[I_RH] = microbe.getRH();

  //pdstate[I_NCE] = getNCE();
  //pdstate[I_NTCB] = getNTCB();
  pdstate[I_FFLC] = ag.getfflc();
  pdstate[I_FFDC] = ag.getffdc();
  pdstate[I_SWFC] = ag.getswfc();
  pdstate[I_DWFC] = ag.getdwfc();

  //pdstate[I_CNVRTC] = ag.getCONVRTFLXC();
  
  pdstate[I_DOCP] = microbe.getDOCPROD();

  pdstate[I_LCHDOC] = soil.getLEACHDOC();

  pdstate[I_ERDPOC] = soil.getERODEPOC();

//  pdstate[I_CO2DISS] = soil.getDISSCO2();

//  pdstate[I_LCHCO2] = soil.getLEACHCO2();

//  pdstate[I_HCO3P] = soil.getFORMHCO3();

//  pdstate[I_LCHHCO3] = soil.getLEACHHCO3();

//  pdstate[I_RHCO3P] = soil.getFORMRHCO3();


  // Nitrogen fluxes in ecosystems

  pdstate[I_AGFRTN] = ag.fertn;

  pdstate[I_BNFIX] = bnfix;

  pdstate[I_SNFIX] = veg.getNFIX();

  pdstate[I_ANFIX] = microbe.getNFIX();

  pdstate[I_INNUP] = veg.getINUPTAKE();
  pdstate[I_INNH4UP] = veg.getINH4UPTAKE();
  pdstate[I_INNO3UP] = veg.getINO3UPTAKE();

  pdstate[I_VNUP] = veg.getNUPTAKE();
  pdstate[I_VNH4UP] = veg.getNH4UPTAKE();
  pdstate[I_VNO3UP] = veg.getNO3UPTAKE();

  pdstate[I_VSUP] = veg.getSUPTAKE();
  pdstate[I_VLUP] = veg.getLUPTAKE();
  pdstate[I_VNMBL] = veg.getNMOBIL();
  pdstate[I_VNRSRB] = veg.getNRESORB();

  //pdstate[I_LTRN] = veg.getLTRFALN();
  pdstate[I_LTRN] = totltfalln;

  pdstate[I_NDCMP] = microbe.getNDECOMP();

  pdstate[I_DONP] = microbe.getDONPROD();

  pdstate[I_GMIN] = microbe.getGROSSNMIN();

  pdstate[I_NH4IMM] = microbe.getIMMNH4();

  pdstate[I_NIMM] = microbe.getIMMOB();

  pdstate[I_NMIN] = microbe.getNETNMIN();

  pdstate[I_AIMMNH4] = soil.getABIMMOB();

  pdstate[I_AMMN] = microbe.getAMMNVOL();

  pdstate[I_AIMMNO3] = ZERO;

  pdstate[I_NTRF] = microbe.getNITRIF();
  pdstate[I_NO3P] = microbe.getNO3PROD();
  pdstate[I_NOP] = microbe.getNOPROD();
  pdstate[I_N2OP] = microbe.getN2OPROD();
  pdstate[I_N2P] = microbe.getN2PROD();

  pdstate[I_DNTRF] = microbe.getDENITRIF();


  pdstate[I_NH3FLX] = soil.getNH3FLUX();

  pdstate[I_NOFLX] = soil.getNOFLUX();
  pdstate[I_N2OFLX] = soil.getN2OFLUX();
  pdstate[I_N2FLX] = soil.getN2FLUX();
  pdstate[I_LCHNO3] = soil.getLEACHNO3();
  pdstate[I_LCHDON] = soil.getLEACHDON();

  pdstate[I_ERDPON] = soil.getERODEPON();

  // Water fluxes

  pdstate[I_AGIRRIG] = ag.irrigate;
  
  pdstate[I_INEET] = soil.getINEET();
  pdstate[I_EET] = soil.getEET();
  pdstate[I_RPERC] = soil.getRPERC();
  pdstate[I_SPERC] = soil.getSPERC();
  pdstate[I_RRUN] = soil.getRRUN();
  pdstate[I_SRUN] = soil.getSRUN();


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void TTEM60::displayOptionalCflx( const scykey& s )
{
  switch( s )
  {
    case GET_LEAF:    cout << "   LEAF "; break;
    case GET_LAI:     cout << "    LAI "; break;
    case GET_FPC:     cout << "    FPC "; break;
    case GET_NSOLC:   cout << "  NSOLC "; break;
    case GET_TSOLC:   cout << "  TSOLC "; break;

    case GET_INGPP:    cout << "  INGPP "; break;
    case GET_GPP:      cout << "   GPP  "; break;
    case GET_INNPP:    cout << "  INNPP "; break;
    case GET_NPP:      cout << "   NPP  "; break;
    case GET_GPR:      cout << "    RA  "; break;
    case GET_RVMNT:    cout << "  RVMNT "; break;
    case GET_RVGRW:    cout << "  RVGRW "; break;
    case GET_LTRC:     cout << "   LTRC "; break;
    case GET_AGSTUBC:  cout << "AGSTUBC "; break;
    case GET_CDCMP:    cout << "  CDCMP "; break;
    case GET_RH:       cout << "    RH  "; break;
    case GET_RSOIL:    cout << "  RSOIL "; break;
    case GET_LCHDOC:   cout << " LCHDOC "; break;
    case GET_NEP:      cout << "   NEP  "; break;
    case GET_NTCB:     cout << "  NTCB  "; break;

    case GET_D40:      cout << " AOT40  "; break;
    case GET_FOZONE:   cout << " FOZONE "; break;

    case GET_CNVRTC:   cout << " CNVRTC "; break;
    case GET_VCNVRTC:  cout << "VCNVRTC "; break;
    case GET_SCNVRTC:  cout << "SCNVRTC "; break;
    case GET_SLASHC:   cout << " SLASHC "; break;
    case GET_CFLX:     cout << "  CFLUX "; break;
    case GET_NCE:     cout << "    NCE  "; break;

    case GET_AGPRDC:   cout << " AGPRODC "; break;
    case GET_PROD10C:  cout << " PROD10C "; break;
    case GET_PROD100C: cout << "PROD100C "; break;
    case GET_RESIDC:   cout << "  RESIDC "; break;

    case GET_AGFPRDC:  cout << " AGFPRDC "; break;
    case GET_PRDF10C:  cout << " PRDF10C "; break;
    case GET_PRDF100C: cout << "PRDF100C "; break;
    case GET_FRESIDC:  cout << " FRESIDC "; break;
    case GET_AGPRDFC:  cout << " AGPRDFC "; break;
    case GET_PRD10FC:  cout << " PRD10FC "; break;
    case GET_PRD100FC: cout << "PRD100FC "; break;
    case GET_TOTPRDFC: cout << "TOTPRDFC "; break;
    case GET_RESIDFC:  cout << " RESIDFC "; break;
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void TTEM60::displayOptionalNflx( const snykey& s )
{
  switch( s )
  {
    case GET_NH4:     cout << "    NH4 "; break;
    case GET_NO3:     cout << "    NO3 "; break;
    case GET_NSOLN:   cout << "  NSOLN "; break;
    case GET_TSOLN:   cout << "  TSOLN "; break;

    case GET_NINP:     cout << " NINPUT "; break;
    case GET_TNDEP:    cout << "TOTNDEP "; break;
    case GET_NH4DEP:   cout << " NH4DEP "; break;
    case GET_NO3DEP:   cout << " NO3DEP "; break;
    case GET_AGFRTN:   cout << "AGFERTN "; break;
    case GET_BNFIX:    cout << "  BNFIX "; break;
    case GET_SNFIX:    cout << "  SNFIX "; break;
    case GET_ANFIX:    cout << "  ANFIX "; break;
    case GET_INNUP:    cout << "  INNUP "; break;
    case GET_INNH4UP:  cout << "INNH4UP "; break;
    case GET_INNO3UP:  cout << "INNO3UP "; break;
    case GET_VNUP:     cout << " UPTAKE "; break;
    case GET_VNH4UP:   cout << "  NH4UP "; break;
    case GET_VNO3UP:   cout << "  NO3UP "; break;
    case GET_VSUP:     cout << " SUPTAK "; break;
    case GET_VLUP:     cout << " LUPTAK "; break;
    case GET_VNMBL:    cout << " NMOBIL "; break;
    case GET_VNRSRB:   cout << " NRTRAN "; break;
    case GET_LTRN:     cout << "   LTRN "; break;
    case GET_AGSTUBN:  cout << "AGSTUBN "; break;
    case GET_NDCMP:    cout << "  NDCMP "; break;
    case GET_DONP:     cout << "   DONP"; break;
    case GET_GMIN:     cout << "GRSNMIN "; break;
    case GET_NH4IMM:   cout << " NH4IMM "; break;
    case GET_NO3IMM:   cout << " NO3IMM "; break;
    case GET_NIMM:     cout << " NIMMOB "; break;
    case GET_NMIN:     cout << "NETNMIN "; break;
    case GET_LCHNO3:   cout << " LCHNO3"; break;
    case GET_LCHDON:   cout << " LCHDON"; break;
    case GET_NLST:     cout << "  NLOST "; break;
    case GET_NTNB:     cout << "  NTNB  "; break;

    case GET_CNVRTN:   cout << " CNVRTN "; break;
    case GET_VCNVRTN:  cout << "VCNVRTN "; break;
    case GET_SCNVRTN:  cout << "SCNVRTN "; break;
    case GET_SLASHN:   cout << " SLASHN "; break;
    case GET_NRETNT:   cout << " NRETNT "; break;
    case GET_NVRTNT:   cout << " NVRTNT "; break;
    case GET_NSRTNT:   cout << " NSRTNT "; break;

    case GET_AGPRDN:   cout << " AGPRODN "; break;
    case GET_PROD10N:  cout << " PROD10N "; break;
    case GET_PROD100N: cout << "PROD100N "; break;
    case GET_RESIDN:   cout << "  RESIDN "; break;

    case GET_AGFPRDN:  cout << " AGFPRDN "; break;
    case GET_PRDF10N:  cout << " PRDF10N "; break;
    case GET_PRDF100N: cout << "PRDF100N "; break;
    case GET_FRESIDN:  cout << " FRESIDN "; break;
    case GET_AGPRDFN:  cout << " AGPRDFN "; break;
    case GET_PRD10FN:  cout << " PRD10FN "; break;
    case GET_PRD100FN: cout << "PRD100FN "; break;
    case GET_TOTPRDFN: cout << "TOTPRDFN "; break;
    case GET_RESIDFN:  cout << " RESIDFN "; break;

    case GET_FIRENDEP: cout << "FIRENDEP "; break;
    
    case GET_L2SN:    cout << "   LCON "; break;
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void TTEM60::displayOptionalSoilTemp( const sstykey& s )
{
  switch( s )
  {
    case GET_FRONTD:   cout << " FRONTD"; break;
    case GET_ACTLAY:   cout << " ACTLAY"; break;
    case GET_THAWPCT:  cout << "THAWPCT"; break;
    case GET_THWBEG1:  cout << "THWBEG1"; break;
    case GET_THWEND1:  cout << "THWEND1"; break;
    case GET_THWBEG2:  cout << "THWBEG2"; break;
    case GET_THWEND2:  cout << "THWEND2"; break;
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void TTEM60::displayOptionalTraceGas( const sgykey& s )
{
  switch( s )
  {
    case GET_AIMMNH4: cout << "AIMMNH4"; break;
    case GET_AMMN:    cout << "AMMNVOL"; break;
    case GET_AIMMNO3: cout << "AIMMNO3"; break;
    case GET_NO3P:    cout << " NO3PRD"; break;
    case GET_NOP:     cout << "  NOPRD"; break;
    case GET_N2OP:    cout << " N2OPRD"; break;
    case GET_N2P:     cout << "  N2PRD"; break;
    case GET_GNLST:   cout << "  NLOST"; break;
    case GET_NH3FLX:  cout << " NH3FLX"; break;
    case GET_NOFLX:   cout << "  NOFLX"; break;
    case GET_N2OFLX:  cout << " N2OFLX"; break;
    case GET_N2FLX:   cout << "  N2FLX"; break;
    case GET_GNMIN:   cout << "NETNMIN"; break;
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void TTEM60::displayOptionalWflx( const swykey& s )
{
  switch( s )
  {
    case GET_SH2O:    cout << " SMOIST"; break;
    case GET_PCTP:    cout << "  PCTP "; break;
    case GET_VSM:     cout << "  VSM  "; break;

    case GET_RAIN:    cout << "  RAIN "; break;
    case GET_SNWFAL:  cout << " SNWFAL"; break;
    case GET_SNWINF:  cout << " SNWINF"; break;
    case GET_AGIRRIG: cout << " IRRIG "; break;
    case GET_PET:     cout << "  PET  "; break;
    case GET_INEET:   cout << "  INEET"; break;
    case GET_EET:     cout << "  EET  "; break;
    case GET_RPERC:   cout << "  RPERC"; break;
    case GET_SPERC:   cout << "  SPERC"; break;
    case GET_RRUN:    cout << "  RRUN "; break;
    case GET_SRUN:    cout << "  SRUN "; break;
    case GET_WYLD:    cout << "  WYLD "; break;
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TTEM60::ecdqc( const int& dcmnt )
{

  int qc = ACCEPT;

  if( vegcb[dcmnt] <= MISSING ) { return qc = 101; }
  if( vegcb[dcmnt] <= MISSING ) { return qc = 102; }
  // cout <<"vegca: " <<vegca[dcmnt]<<" dcmnt " <<dcmnt<<" missing: " <<MISSING <<" xx: " <<10 + vegca[dcmnt]<<endl;
  if( strna[dcmnt] <= MISSING ) { return qc = 103; }
  if( strnb[dcmnt] <= MISSING ) { return qc = 104; }

  if( stona[dcmnt] <= MISSING ) { return qc = 105; }
  if( stonb[dcmnt] <= MISSING ) { return qc = 106; }

  if( solca[dcmnt] <= MISSING ) { return qc = 107; }
  if( solcb[dcmnt] <= MISSING ) { return qc = 108; }

  if( solna[dcmnt] <= MISSING ) { return qc = 109; }
  if( solnb[dcmnt] <= MISSING ) { return qc = 110; }

  if( doccut[dcmnt] <= MISSING ) { return qc = 111; }
  if( doc1a[dcmnt] <= MISSING ) { return qc = 112; }
  if( doc1b[dcmnt] <= MISSING ) { return qc = 113; }
  if( doc2a[dcmnt] <= MISSING ) { return qc = 114; }
  if( doc2b[dcmnt] <= MISSING ) { return qc = 115; }

  if( doncut[dcmnt] <= MISSING ) { return qc = 116; }
  if( don1a[dcmnt] <= MISSING ) { return qc = 117; }
  if( don1b[dcmnt] <= MISSING ) { return qc = 118; }
  if( don2a[dcmnt] <= MISSING ) { return qc = 119; }
  if( don2b[dcmnt] <= MISSING ) { return qc = 120; }

  if( nh4a[dcmnt] <= MISSING ) { return qc = 121; }
  if( nh4b[dcmnt] <= MISSING ) { return qc = 122; }

  if( no3cut[dcmnt] <= MISSING ) { return qc = 123; }
  if( no31a[dcmnt] <= MISSING ) { return qc = 124; }
  if( no31b[dcmnt] <= MISSING ) { return qc = 125; }
  if( no32a[dcmnt] <= MISSING ) { return qc = 126; }
  if( no32b[dcmnt] <= MISSING ) { return qc = 127; }

  if( veg.getUNLEAF12( dcmnt ) <= -9.99 ) { return qc = 128; }
  if( veg.getINITLEAFMX( dcmnt ) <= -9.99 ) { return qc = 129; }

  if( veg.getCMAXCUT( dcmnt ) <= -99.99 ) { return qc = 130; }
  if( veg.getCMAX1A( dcmnt ) <= MISSING ) { return qc = 131; }
  if( veg.getCMAX1B( dcmnt ) <= MISSING ) { return qc = 132; }
  if( veg.getCMAX2A( dcmnt ) <= MISSING ) { return qc = 133; }
  if( veg.getCMAX2B( dcmnt ) <= MISSING ) { return qc = 134; }
  if( veg.getCFALL( dcmnt ) <= -99.99 ) { return qc = 135; }

  if( veg.getRMMAX( dcmnt ) <= -99.99 ) { return qc = 136; }
//  if( veg.getKRA( dcmnt ) <= -99.99 ) { return qc = 136; }
//  if( veg.getKRB( dcmnt ) <= -99.99 ) { return qc = 137; }

  if( veg.getRROOT( dcmnt ) <= -99.99 ) { return qc = 138; }

  if( microbe.getKDCUT( dcmnt ) <= -99.99 ) { return qc = 139; }
  if( microbe.getKD1A( dcmnt ) <= -99.99 ) { return qc = 140; }
  if( microbe.getKD1B( dcmnt ) <= -99.99 ) { return qc = 141; }
  if( microbe.getKD2A( dcmnt ) <= -99.99 ) { return qc = 142; }
  if( microbe.getKD2B( dcmnt ) <= -99.99 ) { return qc = 143; }

  if( microbe.getPROPFTOS( dcmnt ) <= -99.99 ) { return qc = 144; }

  if( soil.getPROPREACTA( dcmnt ) <= -99.99 ) 
  { 
    return qc = 145; 
  }

  if( soil.getPROPREACTB( dcmnt ) <= -99.99 ) 
  { 
    return qc = 146; 
  }

  if( soil.getNSOLPAR( dcmnt ) <= -99.99 ) 
  { 
    return qc = 147; 
  }
  
  if( microbe.getDOCPAR( dcmnt ) <= -99.99 ) { return qc = 148; }
  if( soil.getLCHDOMPAR( dcmnt ) <= -99.99 ) { return qc = 149; }

  if( veg.getNFIXPARA( dcmnt ) <= -99.99 ) { return qc = 150; }
  if( veg.getNFIXPARB( dcmnt ) <= -99.99 ) { return qc = 151; }

  if( veg.getNUPNH4CUT( dcmnt ) <= -99.99 ) { return qc = 152; }
  if( veg.getNUPNH41A( dcmnt ) <= MISSING ) { return qc = 153; }
  if( veg.getNUPNH41B( dcmnt ) <= MISSING ) { return qc = 154; }
  if( veg.getNUPNH42A( dcmnt ) <= MISSING ) { return qc = 155; }
  if( veg.getNUPNH42B( dcmnt ) <= MISSING ) { return qc = 156; }

  if( veg.getNUPNO3CUT( dcmnt ) <= -99.99 ) { return qc = 157; }
  if( veg.getNUPNO31A( dcmnt ) <= MISSING ) { return qc = 158; }
  if( veg.getNUPNO31B( dcmnt ) <= MISSING ) { return qc = 159; }
  if( veg.getNUPNO32A( dcmnt ) <= MISSING ) { return qc = 160; }
  if( veg.getNUPNO32B( dcmnt ) <= MISSING ) { return qc = 161; }

  if( veg.getNFALL( dcmnt ) <= -99.99 ) { return qc = 162; }
  if( veg.getLCCLNC( dcmnt ) <= -99.99 ) { return qc = 163; }

  if( microbe.getNFIXPAR( dcmnt ) <= -99.99 ) { return qc = 164; }

  if( microbe.getNH4IMMCUT( dcmnt ) <= -99.99 ) { return qc = 165; }
  if( microbe.getNH4IMM1A( dcmnt ) <= MISSING ) { return qc = 166; }
  if( microbe.getNH4IMM1B( dcmnt ) <= MISSING ) { return qc = 167; }
  if( microbe.getNH4IMM2A( dcmnt ) <= MISSING ) { return qc = 168; }
  if( microbe.getNH4IMM2B( dcmnt ) <= MISSING ) { return qc = 169; }

  if( microbe.getAMMNPAR( dcmnt ) <= MISSING ) { return qc = 170; }

  if( microbe.getNTRFPARCUT( dcmnt ) <= -99.99 ) { return qc = 171; }
  if( microbe.getNTRFPAR1A( dcmnt ) <= MISSING ) { return qc = 172; }
  if( microbe.getNTRFPAR1B( dcmnt ) <= MISSING ) { return qc = 173; }
  if( microbe.getNTRFPAR2A( dcmnt ) <= MISSING ) { return qc = 174; }
  if( microbe.getNTRFPAR2B( dcmnt ) <= MISSING ) { return qc = 175; }
  
  if( microbe.getINITNTRF( dcmnt ) <= -10000000000 ) { return qc = 176; }
  if( microbe.getALLNTRF( dcmnt ) <= MISSING ) { return qc = 177; }
  if( microbe.getTGMPAR( dcmnt ) <= -99.99 ) { return qc = 178; }

  if( soil.getLCHNO3PARCUT( dcmnt ) <= -99.99 ) { return qc = 179; }
  if( soil.getLCHNO3PAR1A( dcmnt ) <= MISSING ) { return qc = 180; }
  if( soil.getLCHNO3PAR1B( dcmnt ) <= MISSING ) { return qc = 181; }
  if( soil.getLCHNO3PAR2A( dcmnt ) <= MISSING ) { return qc = 182; }
  if( soil.getLCHNO3PAR2B( dcmnt ) <= MISSING ) { return qc = 183; }

  if( microbe.getDONPARCUT( dcmnt ) <= -99.99 ) { return qc = 184; }
  if( microbe.getDONPAR1A( dcmnt ) <= MISSING ) { return qc = 185; }
  if( microbe.getDONPAR1B( dcmnt ) <= MISSING ) { return qc = 186; }
  if( microbe.getDONPAR2A( dcmnt ) <= MISSING ) { return qc = 187; }
  if( microbe.getDONPAR2B( dcmnt ) <= MISSING ) { return qc = 188; }

  if( veg.getINITCNEVEN( dcmnt ) <= MISSING ) { return qc = 189; }
  if( veg.getCNMIN( dcmnt ) <= MISSING ) { return qc = 190; }
  if( veg.getC2NA( dcmnt ) <= MISSING ) { return qc = 191; }
  if( veg.getC2NB( dcmnt ) <= MISSING ) { return qc = 192; }
  if( veg.getC2NMIN( dcmnt ) <= MISSING ) { return qc = 193; }

  if( microbe.getCNSOIL( dcmnt ) <= MISSING ) { return qc = 194; }

  if( veg.getO3PARA( dcmnt ) <= -99.99 ) { return qc = 195; }
  if( veg.getO3PARB( dcmnt ) <= -99.99 ) { return qc = 196; }
  if( veg.getO3PARC( dcmnt ) <= -99.99 ) { return qc = 197; }


  return qc;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TTEM60::ECDsetODEstate( const int& pdcmnt,
                             const double& psiplusc )
{
  // Initialize the NUMEQ state variables used in the 
  //   ODE integrator from ECD and DAT files 
    
  y[I_VEGC] = vegca[pdcmnt] * psiplusc + vegcb[pdcmnt];
  
  if( y[I_VEGC] < ZERO ) { y[I_VEGC] = ZERO; }


  y[I_SOLC] = solca[pdcmnt] * psiplusc + solcb[pdcmnt];
  
  if( y[I_SOLC] < ZERO ) { y[I_SOLC] = ZERO; }

  if (y[I_VEGC]==ZERO)
  {
	  y[I_FOLC]=0.0;
	  y[I_STEMC]=0.0;
	  y[I_CROOTC]=0.0;
	  y[I_FROOTC]=0.0;

  }
  else if (veg.getifwoody(veg.cmnt)==1) //forest
	  {
		 y[I_FOLC] = 0.03 * y[I_VEGC];
		 y[I_FROOTC] = 0.03 * y[I_VEGC];
		 y[I_CROOTC] = 0.18 * y[I_VEGC];
		 y[I_STEMC]  = 0.76 * y[I_VEGC];
	  }
  else if (veg.getifwoody(veg.cmnt)==0) //grass
	  {
		  y[I_FOLC]=0.3* y[I_VEGC];
		  y[I_STEMC]=0.2* y[I_VEGC];
		  y[I_CROOTC]=0.2* y[I_VEGC];
		  y[I_FROOTC]=0.3* y[I_VEGC];
	  }

  else if (veg.getifwoody(veg.cmnt)==2) //shrub
	  {
		  y[I_FOLC] = 0.2 * y[I_VEGC];
		  y[I_FROOTC] = 0.2 * y[I_VEGC];
		  y[I_CROOTC] = 0.3 * y[I_VEGC];
		  y[I_STEMC]  = 0.3 * y[I_VEGC];
	  }

/*
  if( psiplusc <= doccut[pdcmnt] )
  { 
    y[I_DOC]  = doc1a[pdcmnt] * psiplusc + doc1b[pdcmnt];
  }
  else
  { 
    y[I_DOC]  = doc2a[pdcmnt] * psiplusc + doc2b[pdcmnt];
  }
  */
  y[I_DOC] = ZERO;
  y[I_DON] = ZERO;
  if( y[I_DOC] < ZERO ) { y[I_DOC] = ZERO; }

  y[I_DEADWOODC] = ZERO;
  y[I_DEADWOODN] = ZERO;
  y[I_AGE] = 0;
  y[I_CWD] = 0.0;
  y[I_AGR] = 0.0;
  y[I_BGR] = 0.0;
  y[I_AGL] = 0.0;
  y[I_BGL] = 0.0;
  y[I_SOC] = y[I_SOLC] - (y[I_AGR]+y[I_AGL]+y[I_BGR]+y[I_BGL]+y[I_CWD]) ;
  if (y[I_SOC]<0.0) y[I_SOC]=0.0;
  //cout <<"soc: "<<y[I_SOC]<<endl;

//  y[I_CO2G] = ZERO;
  
//  y[I_CO2W] = ZERO;
  
//  y[I_HCO3] = ZERO;
  
//  y[I_RHCO3] = ZERO;
  

  y[I_STRN] = strna[pdcmnt] * psiplusc + strnb[pdcmnt];
  
  if( y[I_STRN] < ZERO ) { y[I_STRN] = ZERO; }
    

  y[I_STON] = stona[pdcmnt] * psiplusc + stonb[pdcmnt];

  if( y[I_STON] < ZERO ) { y[I_STON] = ZERO; }
    

  y[I_SOLN] = solna[pdcmnt] * psiplusc + solnb[pdcmnt];
    
  if( y[I_SOLN] < ZERO ) { y[I_SOLN] = ZERO; }
        
 /*
  if( psiplusc <= doncut[pdcmnt] )
  { 
    y[I_DON]  = don1a[pdcmnt] * psiplusc + don1b[pdcmnt];
  }
  else
  { 
    y[I_DON]  = don2a[pdcmnt] * psiplusc + don2b[pdcmnt];
  }
      
  if( y[I_DON] < ZERO )
  {
    y[I_DON] = ZERO;
  }
  */

  y[I_NH4] = nh4a[pdcmnt] * psiplusc + nh4b[pdcmnt];
  
  if( y[I_NH4] < ZERO ) { y[I_NH4] = ZERO; }
   

  if( psiplusc <= no3cut[pdcmnt] )
  { 
    y[I_NO3] = no31a[pdcmnt] * psiplusc + no31b[pdcmnt];
  }
  else
  { 
    y[I_NO3] = no32a[pdcmnt] * psiplusc + no32b[pdcmnt];
  }
    
  if( y[I_NO3] < ZERO ) { y[I_NO3] = ZERO; }
    
  y[I_SM] = soil.getAWCAPMM() + soil.getWILTPT();

  y[I_RGRW] = ZERO;

  y[I_SGRW] =  ZERO;


  // Initialize all phenology and flux states to zero

  resetODEflux();


  // Reinitialize phenology state variables
  
  y[I_UNRMLF] = veg.getUNLEAF12( pdcmnt );

  y[I_LEAF] = veg.getINITLEAFMX( pdcmnt );  
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM60::getcropecd( const int& dv, const string& agecd )
{
  string agversion;
  string agsitename;
  string agdeveloper;
  string agsitecol;
  string agsiterow;
  string agupdated;

  string agdescription;

  fecd[dv].open( agecd.c_str(), ios::in );

  if( !fecd[dv] )
  {
    cerr << endl << "Cannot open " << agecd;
    cerr << " for ag siteECD input" << endl;
    
    exit( -1 );
  }

  veg.getXMLsiteRootNode( fecd[dv],
                          "siteECD",
                          agversion,
                          agsitename,
                          agsitecol,
                          agsiterow,
                          agdeveloper,
                          agupdated );

  ag.cmnt = veg.getXMLsiteCommunityNode( fecd[dv],
                                         "siteECD",
                                         agdescription );

  vegca[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegca",
                                              ag.cmnt );

  vegcb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegcb",
                                              ag.cmnt );

  strna[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "strna",
                                              ag.cmnt );

  strnb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "strnb",
                                              ag.cmnt );

  stona[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "stona",
                                              ag.cmnt );

  stonb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "stonb",
                                              ag.cmnt );

  solca[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "solca",
                                              ag.cmnt );

  solcb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "solcb",
                                              ag.cmnt );

  solna[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "solna",
                                              ag.cmnt );

  solnb[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "solnb",
                                              ag.cmnt );

  doccut[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "doccut",
                                               ag.cmnt );

  doc1a[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "doc1a",
                                              ag.cmnt );

  doc1b[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "doc1b",
                                              ag.cmnt );

  doc2a[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "doc2a",
                                              ag.cmnt );

  doc2b[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "doc2b",
                                              ag.cmnt );

  doncut[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "doncut",
                                               ag.cmnt );

  don1a[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "don1a",
                                              ag.cmnt );

  don1b[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "don1b",
                                              ag.cmnt );

  don2a[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "don2a",
                                              ag.cmnt );

  don2b[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "don2b",
                                              ag.cmnt );

  nh4a[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "nh4a",
                                             ag.cmnt );

  nh4b[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "nh4b",
                                             ag.cmnt );


  no3cut[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "no3cut",
                                               ag.cmnt );

  no31a[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "no31a",
                                              ag.cmnt );

  no31b[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "no31b",
                                              ag.cmnt );

  no32a[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "no32a",
                                              ag.cmnt );

  no32b[ag.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "no32b",
                                              ag.cmnt );

  veg.setUNLEAF12( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "unleaf12",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setINITLEAFMX( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "initleafmx",
                                                 ag.cmnt ),
                     ag.cmnt );


  veg.setCMAXCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "vegcmaxcut",
                                             ag.cmnt ),
                  ag.cmnt );

  veg.setCMAX1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax1a",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setCMAX1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax1b",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setCMAX2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax2a",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setCMAX2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax2b",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setCFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegcfall",
                                           ag.cmnt ),
                ag.cmnt );

//  veg.setKRA( veg.getXMLcmntArrayDouble( fecd[dv],
//                                         "siteECD",
//                                         "vegkra",
//                                         ag.cmnt ),
//              ag.cmnt );

//  veg.setKRB( veg.getXMLcmntArrayDouble( fecd[dv],
//                                         "siteECD",
//                                         "vegkrb",
//                                         ag.cmnt ),
//              ag.cmnt );

  veg.setRMMAX( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegrmmax",
                                           ag.cmnt ),
                ag.cmnt );

  veg.setRROOT( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegrroot",
                                           ag.cmnt ),
                ag.cmnt );


  microbe.setKDCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "microbekdcut",
                                               ag.cmnt ),
                    ag.cmnt );

  microbe.setKD1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbekd1a",
                                              ag.cmnt ),
                   ag.cmnt );

  microbe.setKD1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbekd1b",
                                              ag.cmnt ),
                   ag.cmnt );

  microbe.setKD2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbekd2a",
                                              ag.cmnt ),
                   ag.cmnt );

  microbe.setKD2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbekd2b",
                                              ag.cmnt ),
                   ag.cmnt );


  microbe.setPROPFTOS( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbepropftos",
                                                  ag.cmnt ),
                       ag.cmnt );

  soil.setPROPREACTA( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "soilpropReactA",
                                                 ag.cmnt ),
                      ag.cmnt );

  soil.setPROPREACTB( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "soilpropReactB",
                                                 ag.cmnt ),
                      ag.cmnt );

  soil.setNSOLPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "soilnonOMpar",
                                              ag.cmnt ),
                   ag.cmnt );

  microbe.setDOCPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbeDOCpar",
                                                ag.cmnt ),
                     ag.cmnt );

  soil.setLCHDOMPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "soillchDOMpar",
                                                ag.cmnt ),
                     ag.cmnt );

  veg.setNFIXPARA( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnfixpara",
                                               ag.cmnt ),
                   ag.cmnt );

  veg.setNFIXPARB( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnfixparb",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNUPNH4CUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "vegnupnh4cut",
                                               ag.cmnt ),
                    ag.cmnt );

  veg.setNUPNH41A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupnh41a",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNUPNH41B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupnh41b",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNUPNH42A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupnh42a",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNUPNH42B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupnh42b",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNUPNO3CUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "vegnupno3cut",
                                               ag.cmnt ),
                    ag.cmnt );

  veg.setNUPNO31A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupno31a",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNUPNO31B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupno31b",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNUPNO32A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupno32a",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNUPNO32B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupno32b",
                                              ag.cmnt ),
                   ag.cmnt );

  veg.setNFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegnfall",
                                           ag.cmnt ),
                ag.cmnt );

  veg.setLCCLNC( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "veglcclnc",
                                            ag.cmnt ),
                 ag.cmnt );

  microbe.setNFIXPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "microbenfixpar",
                                                 ag.cmnt ),
                      ag.cmnt );

  microbe.setNH4IMMCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbenh4immcut",
                                                   ag.cmnt ),
                        ag.cmnt );

  microbe.setNH4IMM1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenh4imm1a",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setNH4IMM1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenh4imm1b",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setNH4IMM2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenh4imm2a",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setNH4IMM2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenh4imm2b",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setAMMNPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "microbeammnpar",
                                                 ag.cmnt ),
                      ag.cmnt );

  microbe.setNTRFPARCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                                    "siteECD",
                                                    "microbentrfparcut",
                                                    ag.cmnt ),
                         ag.cmnt );

  microbe.setNTRFPAR1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbentrfpar1a",
                                                   ag.cmnt ),
                        ag.cmnt );

  microbe.setNTRFPAR1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbentrfpar1b",
                                                   ag.cmnt ),
                        ag.cmnt );

  microbe.setNTRFPAR2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbentrfpar2a",
                                                   ag.cmnt ),
                        ag.cmnt );

  microbe.setNTRFPAR2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbentrfpar2b",
                                                   ag.cmnt ),
                        ag.cmnt );

  microbe.setINITNTRF( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeinitntrf",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setALLNTRF( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "microbeallntrf",
                                                 ag.cmnt ),
                      ag.cmnt );

  microbe.setTGMPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbetgmpar",
                                                ag.cmnt ),
                     ag.cmnt );

  soil.setLCHNO3PARCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "soillchNO3parcut",
                                                   ag.cmnt ),
                        ag.cmnt );

  soil.setLCHNO3PAR1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "soillchNO3par1a",
                                                  ag.cmnt ),
                       ag.cmnt );

  soil.setLCHNO3PAR1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "soillchNO3par1b",
                                                  ag.cmnt ),
                       ag.cmnt );

  soil.setLCHNO3PAR2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "soillchNO3par2a",
                                                  ag.cmnt ), 
                       ag.cmnt );

  soil.setLCHNO3PAR2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "soillchNO3par2b",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setDONPARCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbeDONparcut",
                                                   ag.cmnt ),
                        ag.cmnt );

  microbe.setDONPAR1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeDONpar1a",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setDONPAR1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeDONpar1b",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setDONPAR2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeDONpar2a",
                                                  ag.cmnt ),
                       ag.cmnt );

  microbe.setDONPAR2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeDONpar2b",
                                                  ag.cmnt ),
                       ag.cmnt );
                                                       
  veg.setINITCNEVEN( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "veginitcneven",
                                                ag.cmnt ),
                     ag.cmnt );

  veg.setCNMIN( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegcnmin",
                                           ag.cmnt ),
                ag.cmnt );

  veg.setC2NA( veg.getXMLcmntArrayDouble( fecd[dv],
                                          "siteECD",
                                          "vegc2na",
                                          ag.cmnt ),
               ag.cmnt );

  veg.setC2NB( veg.getXMLcmntArrayDouble( fecd[dv],
                                          "siteECD",
                                          "vegc2nb",
                                          ag.cmnt ),
               ag.cmnt );

  veg.setC2NMIN( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegc2nmin",
                                            ag.cmnt ),
                 ag.cmnt );

  microbe.setCNSOIL( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbecnsoil",
                                                ag.cmnt ),
                     ag.cmnt );

  veg.setO3PARA( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3para",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setO3PARB( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3parb",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.setO3PARC( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3parc",
                                            ag.cmnt ),
                 ag.cmnt );

  veg.endXMLcommunityNode( fecd[dv] );


  fecd[dv].close();

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM60::getenviron( const int& pdyr,
                         const int& pdm,
                         ofstream& rflog1 )
{
  double nextsnowinf;
  double lastyrdst10;
  double diffdst10;

  // Determine monthly potential evapotranspiration

  atms.petjh( atms.getNIRR(), atms.getTAIR(), pdm );



  // Determine contribution of snowmelt to soil moisture

  soil.setSNOWINF( soil.snowmelt( elev, 
                                  atms.getTAIR(), 
                                  atms.getPREVTAIR(), 
                                  soil.getPREVSPACK() ) );

  // Determine new snow pack
  
  soil.setSNOWPACK( (soil.getPREVSPACK() 
                     + atms.getSNOWFALL() 
                     - soil.getSNOWINF()) );


  if( soil.getSNOWPACK() < ZERO )
  {
    soil.setSNOWPACK( ZERO );
  }

  soil.setACTLAYER( soil.getROOTZ() );

  if( 1 == soil.stmflg )
  {
    if ( 0 == pdyr && 0 == pdm ) 
    { 
      soil.stm.setKSWITCH( 0 ); 
    
      soil.stm.initSoilThermalRegime( atms.getAVETAIR() );
    }
    else { soil.stm.setKSWITCH( 1 ); }

    

    nextsnowinf = soil.snowmelt( elev,
                                 soil.stm.getNEXTTAIR(),
                                 atms.getTAIR(),
                                 soil.getSNOWPACK() );

    soil.stm.setNEXTSPACK( (soil.getSNOWPACK() 
                            + soil.stm.getNEXTSNOWFALL()
                            - nextsnowinf) );

    soil.stm.setSnowDensity( veg.cmnt, atms.getPREC(), pdm );

    // Obtain dst10 and tsoil value from last year
    
    lastyrdst10 = soil.getDST10();
    //cout <<"setmonthly soil temp run okay1: "<<endl;

    soil.stm.setMonthlySoilConditions( atms.getPREVTAIR(),
                                       atms.getTAIR(),
                                       soil.stm.getNEXTTAIR(),
                                       soil.getPREVSPACK(),
                                       soil.getSNOWPACK(),
                                       soil.stm.getNEXTSPACK(),
                                       veg.cmnt,
                                       pdm,
                                       rflog1 );

    //cout <<"setmonthly soil temp run okay2: "<<endl;
    // Determine year-to-year different in dst10 for 
    //   particular month
    
    diffdst10 = soil.getDST10() - lastyrdst10;
    
    // Impose active layer restrictions on rooting depth for 
    //   those areas underlain by permafrost 

    if( soil.stm.getTHAWEND1() > soil.getROOTZ() )
    {
      if( (soil.getROOTZ() - soil.stm.getTHAWBEGIN1()) < MINACTLAYERZ )
      {
        soil.setACTLAYER( MINACTLAYERZ ); 
      }
      else
      {
        soil.setACTLAYER( (soil.getROOTZ() - soil.stm.getTHAWBEGIN1()) );
      }
    }
    else
    {
      if( (soil.stm.getTHAWEND1() - soil.stm.getTHAWBEGIN1()) < MINACTLAYERZ )
      {
        soil.setACTLAYER( MINACTLAYERZ ); 
      }
      else
      {
      	if( soil.stm.getTHAWEND2() > soil.getROOTZ()
      	    && soil.stm.getTHAWBEGIN2() > soil.getROOTZ() )
        {
          soil.setACTLAYER( (soil.stm.getTHAWEND1() - soil.stm.getTHAWBEGIN1()) );
        }
      	
      	else if( soil.stm.getTHAWEND2() > soil.getROOTZ()
      	    && soil.stm.getTHAWBEGIN2() < soil.getROOTZ() )
        {
          soil.setACTLAYER( ((soil.stm.getTHAWEND1() - soil.stm.getTHAWBEGIN1()) 
                            + (soil.getROOTZ() - soil.stm.getTHAWBEGIN2())) );
        } 
      	
      	else
        {
          soil.setACTLAYER( ((soil.stm.getTHAWEND1() - soil.stm.getTHAWBEGIN1()) 
                            + (soil.stm.getTHAWEND2() - soil.stm.getTHAWBEGIN2())) );
        }
        
      }
    }

    // Update soil characteristics based on dynamic 
    //   active layer depth
              
    soil.updateActiveLayerRootZ();

    soil.setTSOIL( soil.stm.getTSOIL() );

    if( soil.stm.getDST10() > (atms.getTAIR() + 100.0) )
    {
      soil.setDST10( atms.getTAIR() );
    }
    else
    {  
      soil.setDST10( soil.stm.getDST10() );
    }
        
    // Determine the amount of soil (0 - 10 cm depth) that is 
    //   thawed

    veg.setThawPercent( soil.getPREVDST10(),
                        soil.getDST10(),
                        (soil.getNEXTDST10() + diffdst10) );
  }
  else { veg.setThawPercent( 1.0000000 ); }
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

double TTEM60::getOptionalCflx( const int& optflx )
{

  double outflux;

  switch( optflx )
  {
    case GET_LEAF:    outflux = y[I_LEAF]; break;
    case GET_LAI:     outflux = y[I_LAI]; break;
    case GET_FPC:     outflux = y[I_FPC]; break;
    case GET_NSOLC:   outflux = soil.getNSOLC(); break;
    case GET_TSOLC:   outflux = soil.getTSOLC(); break;

    case GET_INGPP:    outflux = y[I_INGPP]; break;
    case GET_GPP:      outflux = y[I_GPP]; break;
    case GET_INNPP:    outflux = y[I_INNPP]; break;
    case GET_NPP:      outflux = y[I_NPP]; break;
    case GET_GPR:      outflux = y[I_GPR]; break;
    case GET_RVMNT:    outflux = y[I_RVMNT]; break;
    case GET_RVGRW:    outflux = y[I_RVGRW]; break;
    case GET_LTRC:     outflux = y[I_LTRC]; break;
    case GET_CDCMP:    outflux = y[I_CDCMP]; break;
    case GET_RH:       outflux = y[I_RH]; break;
    case GET_RSOIL:    outflux = rsoil; break;
    case GET_LCHDOC:   outflux = y[I_LCHDOC]; break;
    case GET_NEP:      outflux = nep; break;
    case GET_NTCB:     outflux = ntcb; break;

    case GET_D40:      outflux = atms.getAOT40(); break;
    case GET_FOZONE:   outflux = y[I_FOZONE]; break;

    case GET_CNVRTC:   outflux = ag.getCONVRTFLXC();  break;
    case GET_VCNVRTC:  outflux = ag.getVCONVRTFLXC();  break;
    case GET_SCNVRTC:  outflux = ag.getSCONVRTFLXC();  break;
    case GET_SLASHC:   outflux = ag.getSLASHC();  break;
    case GET_CFLX:     outflux = ag.getCFLUX();  break;
    case GET_NCE:      outflux = nce;  break;

    case GET_AGSTUBC:  outflux = ag.getSTUBBLEC(); break;
    case GET_RESIDC:   outflux = ag.getCROPRESIDUEC();  break;

    case GET_AGPRDC:   outflux = ag.getPROD1C();  break;
    case GET_PROD10C:  outflux = ag.getPROD10C();  break;
    case GET_PROD100C: outflux = ag.getPROD100C();  break;

    case GET_FRESIDC:  outflux = ag.getFORMCROPRESIDUEC();  break;

    case GET_AGFPRDC:  outflux = ag.getCROPPRODC();  break;
    case GET_PRDF10C:  outflux = ag.getFORMPROD10C();  break;
    case GET_PRDF100C: outflux = ag.getFORMPROD100C();  break;
    case GET_AGPRDFC:  outflux = ag.getPROD1DECAYC();  break;
    case GET_PRD10FC:  outflux = ag.getPROD10DECAYC();  break;
    case GET_PRD100FC: outflux = ag.getPROD100DECAYC();  break;
    case GET_TOTPRDFC: outflux = ag.getTOTPRODDECAYC();  break;
    case GET_RESIDFC:  outflux = ag.getCROPRESIDUEFLXC();  break;

    default:           outflux = MISSING;
  }

  return outflux;

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

double TTEM60::getOptionalNflx( const int& optflx )
{

  double outflux;

  switch( optflx )
  {
    case GET_NH4:     outflux = y[I_NH4]; break;
    case GET_NO3:     outflux = y[I_NO3]; break;
    case GET_NSOLN:   outflux = soil.getNSOLN(); break;
    case GET_TSOLN:   outflux = soil.getTSOLN(); break;

    case GET_NINP:     outflux = soil.getNINPUT(); break;
    case GET_TNDEP:    outflux = atms.getTOTNDEP(); break;
    case GET_NH4DEP:   outflux = atms.getNH4DEP(); break;
    case GET_NO3DEP:   outflux = atms.getNO3DEP(); break;
    case GET_AGFRTN:   outflux = y[I_AGFRTN]; break;
    case GET_BNFIX:    outflux = y[I_BNFIX]; break;
    case GET_SNFIX:    outflux = y[I_SNFIX]; break;
    case GET_ANFIX:    outflux = y[I_ANFIX]; break;
    case GET_INNUP:    outflux = y[I_INNUP]; break;
    case GET_INNH4UP:  outflux = y[I_INNH4UP]; break;
    case GET_INNO3UP:  outflux = y[I_INNO3UP]; break;
    case GET_VNUP:     outflux = y[I_VNUP]; break;
    case GET_VNH4UP:   outflux = y[I_VNH4UP]; break;
    case GET_VNO3UP:   outflux = y[I_VNO3UP]; break;
    case GET_VSUP:     outflux = y[I_VSUP]; break;
    case GET_VLUP:     outflux = y[I_VLUP]; break;
    case GET_VNMBL:    outflux = y[I_VNMBL]; break;
    case GET_VNRSRB:   outflux = y[I_VNRSRB]; break;
    case GET_LTRN:     outflux = y[I_LTRN]; break;
    case GET_AGSTUBN:  outflux = ag.getSTUBBLEN(); break;
    case GET_NDCMP:    outflux = y[I_NDCMP]; break;
    case GET_DONP:     outflux = y[I_DONP]; break;
    case GET_GMIN:     outflux = y[I_GMIN]; break;
    case GET_NH4IMM:   outflux = y[I_NH4IMM]; break;
    case GET_NIMM:     outflux = y[I_NIMM]; break;
    case GET_NMIN:     outflux = y[I_NMIN]; break;
    case GET_LCHNO3:   outflux = y[I_LCHNO3]; break;
    case GET_LCHDON:   outflux = y[I_LCHDON]; break;
    case GET_NLST:     outflux = soil.getNLOST(); break;
    case GET_NTNB:     outflux = ntnb; break;

    case GET_CNVRTN:   outflux = ag.getCONVRTFLXN();  break;
    case GET_VCNVRTN:  outflux = ag.getVCONVRTFLXN();  break;
    case GET_SCNVRTN:  outflux = ag.getSCONVRTFLXN();  break;
    case GET_SLASHN:   outflux = ag.getSLASHN();  break;
    case GET_NRETNT:   outflux = ag.getNRETENT();  break;
    case GET_NVRTNT:   outflux = ag.getNVRETENT();  break;
    case GET_NSRTNT:   outflux = ag.getNSRETENT();  break;

    case GET_AGPRDN:   outflux = ag.getPROD1N();  break;
    case GET_PROD10N:  outflux = ag.getPROD10N();  break;
    case GET_PROD100N: outflux = ag.getPROD100N();  break;
    case GET_RESIDN:   outflux = ag.getCROPRESIDUEN();  break;

    case GET_AGFPRDN:  outflux = ag.getCROPPRODN();  break;
    case GET_PRDF10N:  outflux = ag.getFORMPROD10N();  break;
    case GET_PRDF100N: outflux = ag.getFORMPROD100N();  break;
    case GET_FRESIDN:  outflux = ag.getFORMCROPRESIDUEN();  break;
    case GET_AGPRDFN:  outflux = ag.getPROD1DECAYN();  break;
    case GET_PRD10FN:  outflux = ag.getPROD10DECAYN();  break;
    case GET_PRD100FN: outflux = ag.getPROD100DECAYN();  break;
    case GET_TOTPRDFN: outflux = ag.getTOTPRODDECAYN();  break;
    case GET_RESIDFN:  outflux = ag.getCROPRESIDUEFLXN();  break;

    case GET_FIRENDEP: outflux = ag.getFIRENDEP();  break;
    
    case GET_L2SN:     if ( y[I_VEGC] != ZERO )
                       {
                         outflux = y[I_STON]/y[I_STRN];
                       }
                       else { outflux = MISSING; } break;
    default:           outflux = MISSING;
  }

  return outflux;

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

double TTEM60::getOptionalSoilTemp( const int& optflx )
{

  double outflux;

  switch( optflx )
  {
    case GET_FRONTD:   outflux = soil.stm.getFRONTD() * 1000.0; break;
    case GET_ACTLAY:   outflux = soil.getACTLAYER() * 1000.0; break;
    case GET_THAWPCT:  outflux = veg.getTHAWPCT() * 100.0; break;
    case GET_THWBEG1:  outflux = soil.stm.getTHAWBEGIN1() * 1000.0; break;
    case GET_THWEND1:  outflux = soil.stm.getTHAWEND1() * 1000.0; break;
    case GET_THWBEG2:  outflux = soil.stm.getTHAWBEGIN2() * 1000.0; break;
    case GET_THWEND2:  outflux = soil.stm.getTHAWEND2() * 1000.0; break;
    default:          outflux = MISSING;
  }

  return outflux;

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

double TTEM60::getOptionalTraceGas( const int& optflx )
{

  double outflux;

  switch( optflx )
  {
    case GET_AIMMNH4: outflux = y[I_AIMMNH4]; break;
    case GET_AMMN:    outflux = y[I_AMMN]; break;
    case GET_AIMMNO3: outflux = y[I_AIMMNO3]; break;
    case GET_NO3P:    outflux = y[I_NO3P]; break;
    case GET_NOP:     outflux = y[I_NOP]; break;
    case GET_N2OP:    outflux = y[I_N2OP]; break;
    case GET_N2P:     outflux = y[I_N2P]; break;
    case GET_GNLST:   outflux = soil.getNLOST(); break;
    case GET_NH3FLX:  outflux = y[I_NH3FLX]; break;
    case GET_NOFLX:   outflux = y[I_NOFLX]; break;
    case GET_N2OFLX:  outflux = y[I_N2OFLX]; break;
    case GET_N2FLX:   outflux = y[I_N2FLX]; break;
    case GET_GNMIN:   outflux = y[I_NMIN]; break;
    default:          outflux = MISSING;
  }

  return outflux;

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

double TTEM60::getOptionalWflx( const int& optflx )
{

  double outflux;

  switch( optflx )
  {
    case GET_SH2O:    outflux = y[I_SM]; break;
    
    case GET_PCTP:    outflux = (y[I_SM]
                                * soil.getACTLAYER() / soil.getROOTZ())
                                / soil.getTOTPOR(); 
                      break;
    
    case GET_VSM:     outflux = 100.0*y[I_SM]
                                / (soil.getROOTZ()*1000.0); 
                      break;

    case GET_RAIN:    outflux = atms.getRAIN(); break;
    case GET_SNWFAL:  outflux = atms.getSNOWFALL(); break;
    case GET_SNWINF:  outflux = soil.getSNOWINF(); break;
    case GET_AGIRRIG: outflux = y[I_AGIRRIG]; break;
    case GET_PET:     outflux = atms.getPET(); break;
    case GET_INEET:   outflux = y[I_INEET]; break;
    case GET_EET:     outflux = y[I_EET]; break;
    case GET_RPERC:   outflux = y[I_RPERC]; break;
    case GET_SPERC:   outflux = y[I_SPERC]; break;
    case GET_RRUN:    outflux = y[I_RRUN]; break;
    case GET_SRUN:    outflux = y[I_SRUN]; break;
    case GET_WYLD:    outflux = soil.getH2OYLD(); break;
    default:          outflux = MISSING;
  }

  return outflux;

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM60::getsitecd( const int& numcmnt, ofstream& rflog1 )
{

  int dv;
  string ecd;

#ifndef PMODE
  cout << "Enter name of the site (.ECD) data file with the";
  cout << " parameter values  for cmnt" << endl;
#endif

  rflog1 << "Enter name of the site (.ECD) data file with the";
  rflog1 << " parameter values cmnt" << endl;

  for( dv = 0; dv < numcmnt; ++dv )
  {

#ifdef PMODE
	*fgo >> ecd;
#else
    cout << (dv+1) << ": ";
    cin >> ecd;
#endif
    rflog1 << (dv+1) << ": " << ecd << endl;

    getsitecd( dv, ecd );
  }
  rflog1 << endl;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM60::getsitecd( const int& dv, const string&  ecd )
{
  string kstring;

  fecd[dv].open( ecd.c_str(), ios::in );

  if( !fecd[dv] )
  {
    cerr << endl;
    cerr << "Cannot open " << ecd << " for site ECD input";
    cerr << endl;
    
    exit( -1 );
  }

  veg.getXMLsiteRootNode( fecd[dv],
                          "siteECD",
                          version,
                          sitename,
                          kstring,
                          kstring,
                          developer,
                          kstring );

  veg.cmnt = veg.getXMLsiteCommunityNode( fecd[dv],
                                          "siteECD",
                                          description );

  vegca[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "vegca",
                                               veg.cmnt );

  vegcb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "vegcb",
                                               veg.cmnt );

  strna[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "strna",
                                               veg.cmnt );

  strnb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "strnb",
                                               veg.cmnt );

  stona[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "stona",
                                               veg.cmnt );

  stonb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "stonb",
                                               veg.cmnt );

  solca[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "solca",
                                               veg.cmnt );

  solcb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "solcb",
                                               veg.cmnt );

  solna[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "solna",
                                               veg.cmnt );

  solnb[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "solnb",
                                               veg.cmnt );

  doccut[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "doccut",
                                                veg.cmnt );

  doc1a[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "doc1a",
                                               veg.cmnt );

  doc1b[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "doc1b",
                                               veg.cmnt );

  doc2a[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "doc2a",
                                               veg.cmnt );

  doc2b[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "doc2b",
                                               veg.cmnt );

  doncut[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "doncut",
                                                veg.cmnt );

  don1a[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "don1a",
                                               veg.cmnt );

  don1b[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "don1b",
                                               veg.cmnt );

  don2a[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "don2a",
                                               veg.cmnt );

  don2b[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "don2b",
                                               veg.cmnt );

  nh4a[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "nh4a",
                                              veg.cmnt );

  nh4b[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "nh4b",
                                             veg.cmnt );


  no3cut[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "no3cut",
                                                veg.cmnt );

  no31a[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "no31a",
                                               veg.cmnt );
  no31b[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "no31b",
                                               veg.cmnt );

  no32a[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "no32a",
                                               veg.cmnt );

  no32b[veg.cmnt] = veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "no32b",
                                               veg.cmnt );

  //cout <<"no31a: "<<no31a[veg.cmnt]<<" no32a: "<<no32a[veg.cmnt]<<" no31b: "<<no31b[veg.cmnt]<<" no32b: "<<no32b[veg.cmnt]<<endl;
  veg.setUNLEAF12( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "unleaf12",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setINITLEAFMX( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "initleafmx",
                                                veg.cmnt ),
                     veg.cmnt );

  veg.setCMAXCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                             "siteECD",
                                             "vegcmaxcut",
                                             veg.cmnt ),
                  veg.cmnt );

  veg.setCMAX1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax1a",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setCMAX1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax1b",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setCMAX2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax2a",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setCMAX2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegcmax2b",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setCFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegcfall",
                                           veg.cmnt ),
                veg.cmnt );
  veg.setLCFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "leafcfall",
                                           veg.cmnt ),
                veg.cmnt );
  veg.setSCFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "stemcfall",
                                           veg.cmnt ),
                veg.cmnt );
  veg.setFCFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "frootcfall",
                                           veg.cmnt ),
                veg.cmnt );
  veg.setCCFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "crootcfall",
                                           veg.cmnt ),
                veg.cmnt );

  //  veg.setKRA( veg.getXMLcmntArrayDouble( fecd[dv],
//                                         "siteECD",
//                                         "vegkra",
//                                         veg.cmnt ),
//              veg.cmnt );

//  veg.setKRB( veg.getXMLcmntArrayDouble( fecd[dv],
//                                         "siteECD",
//                                         "vegkrb",
//                                         veg.cmnt ),
//              veg.cmnt );

  veg.setRMMAX( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegrmmax",
                                           veg.cmnt ),
                veg.cmnt );

  veg.setRROOT( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegrroot",
                                           veg.cmnt ),
                veg.cmnt );

  microbe.setKDCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "microbekdcut",
                                               veg.cmnt ),
                    veg.cmnt );

  microbe.setKD1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbekd1a",
                                              veg.cmnt ),
                   veg.cmnt );

  microbe.setKD1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbekd1b",
                                              veg.cmnt ),
                   veg.cmnt );

  microbe.setKD2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbekd2a",
                                              veg.cmnt ),
                   veg.cmnt );

  microbe.setKD2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "microbekd2b",
                                              veg.cmnt ),
                   veg.cmnt );

  microbe.setPROPFTOS( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbepropftos",
                                                  veg.cmnt ),
                       veg.cmnt );

  soil.setPROPREACTA( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "soilpropReactA",
                                                 veg.cmnt ),
                      veg.cmnt );

  soil.setPROPREACTB( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "soilpropReactB",
                                                 veg.cmnt ),
                      veg.cmnt );

  soil.setNSOLPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "soilnonOMpar",
                                              veg.cmnt ),
                   veg.cmnt );

  microbe.setDOCPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbeDOCpar",
                                                veg.cmnt ),
                     veg.cmnt );

  soil.setLCHDOMPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "soillchDOMpar",
                                                veg.cmnt ),
                     veg.cmnt );

  veg.setNFIXPARA( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnfixpara",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNFIXPARB( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnfixparb",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNUPNH4CUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "vegnupnh4cut",
                                               veg.cmnt ),
                    veg.cmnt );

  veg.setNUPNH41A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupnh41a",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNUPNH41B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupnh41b",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNUPNH42A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupnh42a",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNUPNH42B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupnh42b",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNUPNO3CUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                               "siteECD",
                                               "vegnupno3cut",
                                               veg.cmnt ),
                    veg.cmnt );

  veg.setNUPNO31A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupno31a",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNUPNO31B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupno31b",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNUPNO32A( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupno32a",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNUPNO32B( veg.getXMLcmntArrayDouble( fecd[dv],
                                              "siteECD",
                                              "vegnupno32b",
                                              veg.cmnt ),
                   veg.cmnt );

  veg.setNFALL( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegnfall",
                                           veg.cmnt ),
                veg.cmnt );

  veg.setLCCLNC( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "veglcclnc",
                                            veg.cmnt ),
                 veg.cmnt );

  microbe.setNFIXPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "microbenfixpar",
                                                 veg.cmnt ),
                      veg.cmnt );

  microbe.setNH4IMMCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbenh4immcut",
                                                   veg.cmnt ),
                        veg.cmnt );

  microbe.setNH4IMM1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenh4imm1a",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setNH4IMM1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenh4imm1b",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setNH4IMM2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenh4imm2a",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setNH4IMM2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbenh4imm2b",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setAMMNPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "microbeammnpar",
                                                 veg.cmnt ),
                      veg.cmnt );

  microbe.setNTRFPARCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                                    "siteECD",
                                                    "microbentrfparcut",
                                                    veg.cmnt ),
                         veg.cmnt );

  microbe.setNTRFPAR1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbentrfpar1a",
                                                   veg.cmnt ),
                        veg.cmnt );

  microbe.setNTRFPAR1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbentrfpar1b",
                                                   veg.cmnt ),
                        veg.cmnt );

  microbe.setNTRFPAR2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbentrfpar2a",
                                                   veg.cmnt ),
                        veg.cmnt );

  microbe.setNTRFPAR2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbentrfpar2b",
                                                   veg.cmnt ), 
                        veg.cmnt );

  microbe.setINITNTRF( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeinitntrf",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setALLNTRF( veg.getXMLcmntArrayDouble( fecd[dv],
                                                 "siteECD",
                                                 "microbeallntrf",
                                                 veg.cmnt ),
                      veg.cmnt );

  microbe.setTGMPAR( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbetgmpar",
                                                veg.cmnt ),
                     veg.cmnt );

  soil.setLCHNO3PARCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "soillchNO3parcut",
                                                   veg.cmnt ),
                        veg.cmnt );

  soil.setLCHNO3PAR1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "soillchNO3par1a",
                                                  veg.cmnt ),
                       veg.cmnt );

  soil.setLCHNO3PAR1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "soillchNO3par1b",
                                                  veg.cmnt ),
                       veg.cmnt );

  soil.setLCHNO3PAR2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "soillchNO3par2a",
                                                  veg.cmnt ),
                       veg.cmnt );

  soil.setLCHNO3PAR2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "soillchNO3par2b",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setDONPARCUT( veg.getXMLcmntArrayDouble( fecd[dv],
                                                   "siteECD",
                                                   "microbeDONparcut",
                                                   veg.cmnt ),
                        veg.cmnt );

  microbe.setDONPAR1A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeDONpar1a",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setDONPAR1B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeDONpar1b",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setDONPAR2A( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeDONpar2a",
                                                  veg.cmnt ),
                       veg.cmnt );

  microbe.setDONPAR2B( veg.getXMLcmntArrayDouble( fecd[dv],
                                                  "siteECD",
                                                  "microbeDONpar2b",
                                                  veg.cmnt ),
                       veg.cmnt );

  veg.setINITCNEVEN( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "veginitcneven",
                                                veg.cmnt ),
                     veg.cmnt );

  veg.setCNMIN( veg.getXMLcmntArrayDouble( fecd[dv],
                                           "siteECD",
                                           "vegcnmin",
                                           veg.cmnt ),
                veg.cmnt );

  veg.setC2NA( veg.getXMLcmntArrayDouble( fecd[dv],
                                          "siteECD",
                                          "vegc2na",
                                          veg.cmnt ),
               veg.cmnt );

  veg.setC2NB( veg.getXMLcmntArrayDouble( fecd[dv],
                                          "siteECD",
                                          "vegc2nb",
                                          veg.cmnt ),
               veg.cmnt );

  veg.setC2NMIN( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "vegc2nmin",
                                            veg.cmnt ),
                 veg.cmnt );

  microbe.setCNSOIL( veg.getXMLcmntArrayDouble( fecd[dv],
                                                "siteECD",
                                                "microbecnsoil",
                                                veg.cmnt ),
                     veg.cmnt );

  veg.setO3PARA( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3para",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setO3PARB( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3parb",
                                            veg.cmnt ),
                 veg.cmnt );

  veg.setO3PARC( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "o3parc",
                                            veg.cmnt ),
                 veg.cmnt );
  microbe.setF0( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "factor0",
                                            veg.cmnt ),
                 veg.cmnt );
  microbe.setF1( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "factor1",
                                            veg.cmnt ),
                 veg.cmnt );
  microbe.setF2( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "factor2",
                                            veg.cmnt ),
                 veg.cmnt );
  microbe.setF3( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "factor3",
                                            veg.cmnt ),
                 veg.cmnt );
  microbe.setF4( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "factor4",
                                            veg.cmnt ),
                 veg.cmnt );
  microbe.setF5( veg.getXMLcmntArrayDouble( fecd[dv],
                                            "siteECD",
                                            "factor5",
                                            veg.cmnt ),
                 veg.cmnt );


  veg.endXMLcommunityNode( fecd[dv] );
  //cout <<"cmnt: " <<veg.cmnt<<" vegca: "<<vegca[veg.cmnt]<<endl;

  fecd[dv].close();

};

/* *************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void TTEM60::initrun( ofstream& rflog1 )
{
  
  avlnflag = nfeed = rheqflag = 0;

/* **************************************************************
		  Run Model with Nitrogen Limitation?
************************************************************** */


  #ifdef CALIBRATE_TEM
    cout << "Do you want to start calibration allowing available ";
    cout << "N to fluctuate?  ";
    
    if( 'Y' == toupper( getch() ) )
    {
      avlnflag = 1;
      cout << "Y" << endl;
    }
    else
    {
      avlnflag = 0;
      cout << "N" << endl;
    }
  #else

    #ifdef PMODE
		*fgo >> avlnflag;
	#else
		cout << endl << "Do you want to allow available N to fluctuate?";
		cout << endl;
		cout << "  Enter 0 for No" << endl;
		cout << "  Enter 1 for Yes: ";

		cin >> avlnflag;
	#endif

    rflog1 << endl << "Do you want to allow available N to fluctuate?";
    rflog1 << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "avlnflag = " << avlnflag << endl << endl;
  #endif
  
  baseline = initbase = 0;
  
  #ifdef CALIBRATE_TEM
    cout << "Do you want to start calibration with N feedback ";
    cout << "on GPP?  ";
    
    if( 'Y' == toupper( getch() ) )
    {
      nfeed = 1;
      cout << "Y" << endl;
//      cout << "Do you want to solve for baseline soil nitrogen?  ";
    
//      if ( 'Y' == toupper( getch() ) )
//      {
//         baseline = 1;
//         cout << "Y" << endl;
//      }
//      else { cout << "N" << endl; }
    }
    else
    {
       nfeed = 0;
       cout << "N" << endl;
    }
  #else

    #ifdef PMODE
		*fgo >> nfeed;
	#else
		cout << endl << "Do you want nitrogen feedback on GPP?" << endl;
		cout << "  Enter 0 for No" << endl;
		cout << "  Enter 1 for Yes: ";

		cin >> nfeed;
	#endif

	rflog1 << endl << "Do you want nitrogen feedback on GPP?";
    rflog1 << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "nfeed = " << nfeed << endl << endl;

    if( 1 == nfeed )
    {

     #ifdef PMODE
    	*fgo >> initbase;
	 #else

    	cout << endl;
    	cout << "Do you want to solve for baseline soil nitrogen?";
    	cout<< endl;
    	cout << "  Enter 0 for No" << endl;
    	cout << "  Enter 1 for Yes: ";

    	cin >> initbase;
	 #endif
      baseline = initbase;

      rflog1 << endl;
      rflog1 << "Do you want to solve for baseline soil nitrogen?";
      rflog1 << endl;
      rflog1 << "  Enter 0 for No" << endl;
      rflog1 << "  Enter 1 for Yes: " << endl;
      rflog1 << "baseline = " << baseline << endl << endl;
    }
  #endif
  
/* **************************************************************
			 Run Model with Moisture Limitation?
************************************************************** */

  moistlim = 0;
  
  #ifdef CALIBRATE_TEM
    cout << "Do you want to start calibration with moisture ";
    cout << "limitation?  ";
  
    if( 'Y' == toupper( getch() ) )
    {
      moistlim = 1;
      cout << "Y" << endl;
    }
    else
    {
      moistlim = 0;
      cout << "N" << endl;
    }
  #else

	#ifdef PMODE
		*fgo >> moistlim;
	#else

		cout << endl;
		cout << "Do you want to run the model with moisture limitation?";
		cout << endl;
		cout << "  Enter 0 for No" << endl;
		cout << "  Enter 1 for Yes: ";

		cin >> moistlim;
	#endif

    rflog1 << endl;
    rflog1 << "Do you want to run the model with moisture limitation?";
    rflog1 << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "moistlim = " << moistlim << endl << endl;
  #endif
  
/****************************************************************
			 Run Model with Ozone?
************************************************************** */

  o3flag = 0;
  
  #ifdef CALIBRATE_TEM
    cout << "Do you want to start calibration with ozone?  ";
    
    if( 'Y' == toupper(getch()) )
    {
      o3flag = 1;
      cout << "Y" << endl;
    }
    else
    {
      o3flag = 0;
      cout << "N" << endl;
    }  
  #else

	#ifdef PMODE
		*fgo >> o3flag;
	#else

    	cout << endl;
    	cout << "Do you want to run the model with ozone?" << endl;
    	cout << "  Enter 0 for No" << endl;
    	cout << "  Enter 1 for Yes: ";

    	cin >> o3flag;
	#endif

    rflog1 << endl;
    rflog1 << "Do you want to run the model with ozone?" << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "o3flag = " << o3flag << endl << endl;
  #endif

/* ***************************************************************
	       Details for Steady State Conditions
************************************************************** */

  #ifdef CALIBRATE_TEM
    double co2level;
    double dc2n;

    equil = 0;
    adapttol = 0;
    intbomb = 0;
    tolbomb = 0;
    cout << "Do you want the model to stop at steady state ";
    cout << "conditions?  ";
    
    if( 'Y' == toupper( getch() ) )
    {
      equil = 1;
      
      cout << "Y" << endl;
      cout << endl;
      cout << "How many years do you want to wait before checking ";
      cout << "equilibrium conditions?  ";

      cin >> strteq;

      strteq *= 12;

      if( 0 == nfeed )
      {
        cout << "Do you want decomposition to come into ";
        cout << "equilibrium?  ";
       
        if( 'Y' == toupper(getch()) )
        {
	  rheqflag = 1;
	  cout << "Y" << endl;
        }
        else { cout << "N" << endl; }
      }

      cout << "Enter the absolute tolerance for the water cycle?  ";

      cin >> wtol;

      cout << "Enter the absolute tolerance for the carbon cycle?  ";

      cin >> ctol;

      if( 1 == nfeed )
      {
        rheqflag = 1;
        
        cout << "Enter the absolute tolerance for the nitrogen ";
        cout << "cycle?  ";

        cin >> ntol;
      }

      cout << "Do you want to have the integrator tolerance ";
      cout << "adjusted automatically?  ";
      
      if( 'Y' == toupper( getch() ) )
      {
        cout << "Y" << endl;
        
        adapttol = 1;
        
        askODE( rflog1 );
        
        cout << "Enter the maximum number of years for the model ";
        cout << "to run:  ";

        cin >> maxyears;

        cout << "Enter the maximum number of attempts to reach a ";
        cout << "solution:  ";

        cin >> maxnrun;
      }
      else { cout << "N" << endl; }
    }
    else { cout << "N" << endl << endl; }

    if ( 0 == adapttol ) { askODE( rflog1 ); }

    cout << endl << "Enter the level of carbon dioxide in ppmv: ";
    
    cin >> co2level;
    
    atms.setCO2LEVEL( co2level );

    cout << endl << endl << endl;
    cout << "Enter the factor for changing C:N per ppmv of ";
    cout << "enhanced CO2: " << endl;
    cout << "               (Enter 0.0 for no change)" << endl;
    
    cin >> dc2n;
    
    veg.setDC2N( dc2n );
    
    cout << endl;
  #else
    maxyears = 0;
    maxnrun = 0;

	#ifdef PMODE
		*fgo >> strteq;
	#else
		cout << endl;
    	cout << "How many years do you want to wait before checking";
    	cout << " equilibrium conditions? ";

    	cin >> strteq;
	#endif
    rflog1 << endl;
    rflog1 << "How many years do you want to wait before checking";
    rflog1 << " equilibrium conditions? ";
    rflog1 << endl;
    rflog1 << "strteq = " << strteq << endl << endl;

	#ifdef PMODE
		*fgo >> maxyears;
	#else
		cout << endl;
		cout << "Enter the maximum number of years for the model to run:";

		cin >> maxyears;
	#endif

    rflog1 << endl;
    rflog1 << "Enter the maximum number of years for the model to run: ";
    rflog1 << endl;
    rflog1 << "maxyears = " << maxyears << endl << endl;

    runsize = maxyears;

	#ifdef PMODE
		*fgo >> maxnrun;
	#else
		cout << endl;
		cout << "Enter the maximum number of attempts to reach a solution: ";

		cin >> maxnrun;
	#endif
    rflog1 << endl;
    rflog1 << "Enter the maximum number of attempts to reach a solution: ";
    rflog1 << endl;
    rflog1 << "maxnrun = " << maxnrun << endl << endl;

    if( 0 == nfeed )
    {
	  #ifdef PMODE
    	*fgo >> rheqflag;
	  #else
    	cout << endl;
    	cout << "Do you want decomposition to come into equilibrium? ";
    	cout << "  Enter 0 for No" << endl;
    	cout << "  Enter 1 for Yes: ";

    	cin >> rheqflag;
	  #endif
      rflog1 << endl;
      rflog1 << "Do you want decomposition to come into equilibrium? ";
      rflog1 << endl;
      rflog1 << "  Enter 0 for No" << endl;
      rflog1 << "  Enter 1 for Yes: " << endl;
      rflog1 << "rheqflag = " << rheqflag << endl << endl;
    }
   
    wtol = 1000.0;

	#ifdef PMODE
		*fgo >> wtol;
	#else
		cout << endl;
		cout << "What absolute tolerance do you want to use for";
		cout << " checking equilibrium";
		cout << endl;
		cout << "of the water cycle? ";

		cin >> wtol;
	#endif
    rflog1 << endl;
    rflog1 << "What absolute tolerance do you want to use for";
    rflog1 << " checking equilibrium";
    rflog1 << endl;
    rflog1 << "of the water cycle? wtol = " << wtol;
    rflog1 << endl << endl;

    ctol = 1000.0;

	#ifdef PMODE
		*fgo >> ctol;
	#else
		cout << endl;
		cout << "What absolute tolerance do you want to use for";
		cout << " checking equilibrium";
		cout << endl;
		cout << "of the carbon cycle? ";

		cin >> ctol;
	#endif
    rflog1 << endl;
    rflog1 << "What absolute tolerance do you want to use for";
    rflog1 << " checking equilibrium";
    rflog1 << endl;
    rflog1 << "of the carbon cycle?" << endl;
    rflog1 << "ctol = " << ctol << endl << endl;

    ntol = 1000.0;

    if( 1 == nfeed )
    {
      rheqflag = 1;
	  #ifdef PMODE
      	  *fgo >> ntol;
	  #else
      	  cout << endl;
      	  cout << "What absolute tolerance do you want to use for";
      	  cout << " checking equilibrium";
      	  cout << endl;
      	  cout << "of the nitrogen cycle? ";

      	  cin >> ntol;
	  #endif
      rflog1 << endl;
      rflog1 << "What absolute tolerance do you want to use for";
      rflog1 << " checking equilibrium";
      rflog1 << endl;
      rflog1 << "of the nitrogen cycle?" << endl;
      rflog1 << "ntol = " << ntol << endl << endl;
    }

    if( 0 == equil )
    {

	 #ifdef PMODE
    	*fgo >> outputstartyr;
    	*fgo >> outputendyr;
    	*fgo >> diffyr;
	 #else
    	cout << endl << endl;
    	cout << "What year do you want to start collecting output data? ";

    	cin >> outputstartyr;

      	cout << endl << endl;
      	cout << "What year do you want to stop collecting output data? ";

      	cin >> outputendyr;

        cout << "How often (x years) should data be collected";
        cout << " after the initial year? ";

        cin >> diffyr;
	  #endif
      rflog1 << endl << endl;
      rflog1 << "What year do you want to start collecting output data? ";
      rflog1 << "outputstartyr = " << outputstartyr << endl;

      rflog1 << endl << endl;
      rflog1 << "What year do you want to stop collecting output data? ";
      rflog1 << "outputendyr = " << outputendyr << endl;

      rflog1 << "How often (x years) should data be collected";
      rflog1 << " after the initial year? ";
      rflog1 << "diffyr = " << diffyr << endl;
    }
  #endif

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM60::massbal( double y[NUMEQ], 
                      double prevy[MAXSTATE] )
{

//  double tmprootz;
  
//  if( 1 == dbugflg )
//  {
//    cout << endl << "Mass bal NH4 = " << y[I_NH4] << endl;
//    exit( -1 );
//  }

  if( (y[I_SM] - prevy[I_SM]) != (soil.getSNOWINF() + atms.getRAIN() 
      + y[I_AGIRRIG] - y[I_RPERC] - y[I_EET] - y[I_SPERC]) )
  {
    y[I_SPERC] = soil.getSNOWINF() 
                 + atms.getRAIN() 
                 + y[I_AGIRRIG] 
                 - y[I_RPERC] 
                 - y[I_EET]
                 - y[I_SM] 
                 + prevy[I_SM];
  }

  if( y[I_SPERC] < ZERO ) { y[I_SPERC] = ZERO; }


  if( (y[I_RGRW] - prevy[I_RGRW]) != (y[I_RPERC] - y[I_RRUN]) )
  {
    y[I_RRUN] = y[I_RPERC] - y[I_RGRW] + prevy[I_RGRW];
  }



  if( y[I_RRUN] < ZERO ) { y[I_RRUN] = ZERO; }


  if( (y[I_SGRW] - prevy[I_SGRW]) != (y[I_SPERC] - y[I_SRUN]) )
  {
    y[I_SRUN] = y[I_SPERC] - y[I_SGRW] + prevy[I_SGRW];
  }

  if( y[I_SRUN] < ZERO ) { y[I_SRUN] = ZERO; }

 
/************************* Carbon Cycle Balances **************************/

  if( y[I_INNPP] < y[I_NPP] ) { y[I_INNPP] = y[I_NPP]; }

  if( y[I_INGPP] < y[I_GPP] ) { y[I_INGPP] = y[I_GPP]; }

  if( y[I_GPR] != y[I_GPP] - y[I_NPP] )
  {
    y[I_GPR] = y[I_GPP] - y[I_NPP];
  }

  if( y[I_GPR] != y[I_RVMNT] + y[I_RVGRW] )
  {
    y[I_RVGRW] = y[I_GPR] - y[I_RVMNT];
  }


  if( y[I_VEGC] - prevy[I_VEGC]
      != y[I_NPP] - y[I_LTRC] )
  {
    //y[I_LTRC] = y[I_NPP]- y[I_VEGC] + prevy[I_VEGC];
	  //if (disturbflag ==0) cout <<"lon: "<<longitude<<" lat: "<<latitude<<" vegc1: "<<y[I_VEGC]<<" vegc0: "<<prevy[I_VEGC]<<" litc: "<<y[I_LTRC]<<" npp: "<<y[I_NPP]<<endl;
  }


//  if( y[I_ERDPOC] < 0.0 ) { y[I_ERDPOC] = 0.000000; } 
  y[I_ERDPOC] = ZERO; 
  /*//closed by cgs2014
  if( y[I_SOLC] - prevy[I_SOLC] != y[I_LTRC] + ag.getSLASHC() 
       - ag.getSCONVRTFLXC() - y[I_CDCMP] - y[I_ERDPOC]+veg.deadwoodltc)
  {
    y[I_CDCMP] = y[I_LTRC] 
                 + ag.getSLASHC()
                 + veg.deadwoodltc
                 - ag.getSCONVRTFLXC()
                 - y[I_ERDPOC] 
                 - y[I_SOLC] 
                 + prevy[I_SOLC];
  }
  */
  if( y[I_DOCP] < ZERO ) { y[I_DOCP] = ZERO; }

  if( y[I_LCHDOC] < ZERO ) { y[I_LCHDOC] = ZERO; }

  if( y[I_DOC] - prevy[I_DOC] != y[I_DOCP] - y[I_LCHDOC] )
  {
    if( y[I_DOCP] - y[I_DOC] + prevy[I_DOC] > ZERO ) 
    {
      y[I_LCHDOC] = y[I_DOCP] - y[I_DOC] + prevy[I_DOC];
    }
    else
    {
      y[I_DOCP] = y[I_DOC] - prevy[I_DOC] + y[I_LCHDOC];
    }    
  }
  
  /*
  if( y[I_CDCMP] != y[I_RH] + y[I_DOCP] )
  {
    y[I_RH] = y[I_CDCMP] - y[I_DOCP];
  }
  */


  /*********************Nitrogen Cycle Balances**********************/

  if( y[I_VNH4UP] < ZERO ) { y[I_VNH4UP] = ZERO; }
  if( y[I_VNO3UP] < ZERO ) { y[I_VNO3UP] = ZERO; }

  if( y[I_VNUP] != y[I_VNH4UP] + y[I_VNO3UP] )
  {
    y[I_VNUP] = y[I_VNH4UP] + y[I_VNO3UP];
  }

  if( y[I_VNUP] < ZERO ) { y[I_VNUP] = ZERO; }

  if( y[I_INNUP] < y[I_VNUP] ) { y[I_INNUP] = y[I_VNUP]; }

  if( y[I_VSUP] < ZERO ) { y[I_VSUP] = ZERO; }

  if( y[I_VSUP] > y[I_VNUP] ) { y[I_VSUP] = y[I_VNUP]; }

  if( y[I_VLUP] != y[I_VNUP] - y[I_VSUP] )
  {
    y[I_VLUP] = y[I_VNUP] - y[I_VSUP];
  }


  // DWK modified the following conditions for checking the mass
  //   balance on STON on 0020401

  if( y[I_STON] - prevy[I_STON]
       != y[I_VLUP] + y[I_VNRSRB] - y[I_VNMBL] )
  {
    y[I_VNRSRB] = y[I_STON] 
                  - prevy[I_STON] 
                  + y[I_VNMBL]
                  - y[I_VLUP];
  }

  if( y[I_BNFIX] < ZERO ) { y[I_BNFIX] = ZERO; }
  if ( y[I_SNFIX] < ZERO ) { y[I_SNFIX] = ZERO; }
  
  if( y[I_BNFIX] != y[I_SNFIX] + y[I_ANFIX] )
  {
    y[I_ANFIX] = y[I_BNFIX] - y[I_SNFIX];
  }    
  
  if( y[I_ANFIX] < ZERO ) { y[I_ANFIX] = ZERO; }
  /*
  if( y[I_STRN] - prevy[I_STRN] != y[I_VSUP] - y[I_LTRN]
       - y[I_VNRSRB] + y[I_VNMBL] + y[I_SNFIX] )
  {
    y[I_LTRN] = y[I_VSUP] 
                - y[I_VNRSRB] 
                + y[I_VNMBL]
                + y[I_SNFIX] 
                - y[I_STRN] 
                + prevy[I_STRN];
  }
  */
//  if( y[I_ERDPON] < ZERO ) { y[I_ERDPON] = ZERO; }
  y[I_ERDPON] = ZERO;
  /*
  if( y[I_SOLN] - prevy[I_SOLN] != y[I_LTRN] + ag.getSLASHN()
      - ag.getSCONVRTFLXN() - ag.getNSRETENT() - y[I_NDCMP] 
      - y[I_ERDPON] + y[I_ANFIX] + veg.deadwoodltn )
      //- y[I_ERDPON]  + veg.deadwoodltn )
  {
    y[I_NDCMP] = y[I_LTRN] 
                 + ag.getSLASHN() 
                 + veg.deadwoodltn
                 - ag.getSCONVRTFLXN()
                 - ag.getNSRETENT() 
                 + y[I_ANFIX]
                 - y[I_ERDPON]
                 - y[I_SOLN]  
                 + prevy[I_SOLN];
  }
  */
  if( y[I_DONP] < ZERO ) { y[I_DONP] = ZERO; }

  if( y[I_LCHDON] < ZERO ) { y[I_LCHDON] = ZERO; }

  if( y[I_DON] - prevy[I_DON] != y[I_DONP] - y[I_LCHDON] )
  {
    if( y[I_DONP] - y[I_DON] + prevy[I_DON] > ZERO ) 
    {
      y[I_LCHDON] = y[I_DONP] - y[I_DON] + prevy[I_DON];
    }
    else
    {
      y[I_DONP] = y[I_DON] - prevy[I_DON] + y[I_LCHDON];
    }    
  }
 
   
  //if( y[I_NDCMP] != y[I_NMIN] + y[I_DONP] )
  {
    //y[I_NMIN] = y[I_NDCMP] - y[I_DONP];
  }

  if( y[I_NMIN] != y[I_GMIN] - y[I_NH4IMM] )
  {
    if ( y[I_NMIN] < ZERO )
    {
      y[I_NH4IMM] = y[I_GMIN] - y[I_NMIN];
      
      if( y[I_NH4IMM] > y[I_NH4] + y[I_GMIN] + atms.getNH4DEP() 
                        + ag.getNRETENT() + ag.getFIRENDEP())
      {
        y[I_NH4IMM] = y[I_NH4] 
                      + y[I_GMIN] 
                      + atms.getNH4DEP()
                      + ag.getNRETENT()
                      + ag.getFIRENDEP();
      }
    }

    y[I_GMIN] = y[I_NMIN] + y[I_NH4IMM];
    y[I_NIMM] = y[I_NH4IMM];
    
  }


  if( y[I_AGFRTN] < ZERO ) { y[I_AGFRTN] = ZERO; }


  if( y[I_NTRF] < ZERO ) { y[I_NTRF] = ZERO; }

  if( y[I_NH4] - prevy[I_NH4]
      != atms.getNH4DEP() + ag.getNRETENT() + ag.getFIRENDEP() + y[I_GMIN] 
         - y[I_NH4IMM] - y[I_VNH4UP] - y[I_AMMN] - y[I_NTRF] )
  {
    y[I_AIMMNH4] =  atms.getNH4DEP() 
                    + ag.getNRETENT()
                    + ag.getFIRENDEP()
                    + y[I_GMIN] 
                    - y[I_NH4IMM]
                    - y[I_VNH4UP] 
                    - y[I_AMMN] 
                    - y[I_NTRF]
                    - y[I_NH4] 
                    + prevy[I_NH4];
  }

  if( y[I_NOP] < ZERO ) {  y[I_NOP] = ZERO; }

  if( y[I_NO3P] < ZERO ) {  y[I_NO3P] = ZERO; }

  if( y[I_NTRF] != y[I_NO3P] + y[I_NOP] )
  {
    y[I_NO3P] = y[I_NTRF] - y[I_NOP];

    if( y[I_NO3P] < ZERO )
    {
      y[I_NO3P] = ZERO;
      y[I_NOP] = y[I_NTRF];
    }
  }

  if( y[I_N2OP] < ZERO ) {  y[I_N2OP] = ZERO; }
  if( y[I_N2P] < ZERO ) {  y[I_N2P] = ZERO; }

  if( y[I_LCHNO3] < ZERO ) { y[I_LCHNO3] = ZERO; }

  if( y[I_NO3] - prevy[I_NO3]
       != atms.getNO3DEP() + y[I_AGFRTN] + y[I_NO3P] 
          - y[I_DNTRF] - y[I_VNO3UP] -  y[I_LCHNO3] )
  {
    if( soil.getLCHNO3PAR() > ZERO )
    {	
      if( atms.getNO3DEP() + y[I_AGFRTN] + y[I_NO3P] - y[I_DNTRF] 
           - y[I_VNO3UP] - y[I_LCHNO3] - y[I_NO3] 
           + prevy[I_NO3] > ZERO )
      {
        y[I_LCHNO3] =  atms.getNO3DEP() + y[I_AGFRTN] + y[I_NO3P] 
                       - y[I_DNTRF] - y[I_VNO3UP] 
                       - y[I_NO3] + prevy[I_NO3];
      }
      else
      {
        y[I_AIMMNO3] =  atms.getNO3DEP() + y[I_AGFRTN] + y[I_NO3P] 
                        - y[I_VNO3UP] - y[I_DNTRF] - y[I_LCHNO3] 
                        - y[I_NO3] + prevy[I_NO3];
      }
    }
    else
    {
      if( atms.getNO3DEP() + y[I_AGFRTN] + y[I_NO3P] - y[I_DNTRF] 
           - y[I_VNO3UP] - y[I_LCHNO3] - y[I_NO3] 
           + prevy[I_NO3] > ZERO )
      {
        y[I_DNTRF] = atms.getNO3DEP() + y[I_AGFRTN] + y[I_NO3P] 
                     - y[I_VNO3UP] - y[I_LCHNO3] - y[I_NO3] 
                     + prevy[I_NO3];      
      }
      else
      {
        y[I_AIMMNO3] =  atms.getNO3DEP() + y[I_AGFRTN] + y[I_NO3P] 
                        - y[I_VNO3UP] - y[I_LCHNO3] - y[I_DNTRF] 
                        - y[I_NO3] + prevy[I_NO3];
      }
    }
  }
  
  if( y[I_NH3FLX] < ZERO ) {  y[I_NH3FLX] = ZERO; }

  if( y[I_DNTRF] != y[I_N2OP] + y[I_N2P] )
  {
    y[I_N2P] = y[I_DNTRF] - y[I_N2OP];
 
    if( y[I_N2P] < ZERO )
    {
      y[I_N2P] = ZERO;
      y[I_N2OP] = y[I_DNTRF];
    }
  }   

  if( y[I_NOFLX] != y[I_NOP] ) { y[I_NOFLX] = y[I_NOP]; }  

//  if( y[I_NOFLX] < ZERO ) {  y[I_NOFLX] = ZERO; }

  if( y[I_N2OFLX] != y[I_N2OP] ) { y[I_N2OFLX] = y[I_N2OP]; }  
  
//  if( y[I_N2OFLX] < ZERO ) {  y[I_N2OFLX] = ZERO; }
  
  if( y[I_N2FLX] != y[I_N2P] ) { y[I_N2FLX] = y[I_N2P]; }  
  
//  if( y[I_N2FLX] < ZERO ) {  y[I_N2FLX] = ZERO; }


};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

int TTEM60::monthlyTransient( const int& pdyr,
                              const int& pdm,
                              const double& ptol,
                              ofstream& pflog1 )
{
  endeq = 0;
  initFlag = 1;

  // Reset to appropriate rooting depth and texture-dependent 
  //   vegetation parameters for each cohort during each month
  //   (Note: microbe and soil parameters are the same between
  //   disturbed [ag.cmnt] and natural [veg.cmnt] land cover 

  if( 0 == ag.state 
      || (ag.state > 0 && disturbflag > 0 && pdm < (disturbmonth-1)) )
  {
    soil.updateRootZ( veg.cmnt );

    veg.resetEcds( veg.cmnt, soil.getPSIPLUSC() );
  }
  else
  {
    soil.updateRootZ( ag.cmnt );
    
    veg.resetEcds( ag.cmnt, soil.getPSIPLUSC() );
  }
 

  // Reset texture-dependent parameters for microbes 
  //   and soil in current cohort

  microbe.resetEcds( veg.cmnt, soil.getPSIPLUSC() );

  soil.resetEcds( veg.cmnt );


  // Update monthly carbon, nitrogen and water dynamics of
  //   terrestrial ecosystem
  //if (veg.cmnt == 4 && longitude ==-80.75 && latitude == 51.25 && pdyr <200) cout <<"lat: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" ntcb0: "<<y[I_NTCB]<<endl;
  stepmonth( pdyr, pdm, intflag, ptol, pflog1 );
  //if (veg.cmnt == 4 && longitude ==-80.75 && latitude == 51.25 && pdyr <200) cout <<"lat: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" ntcb1: "<<y[I_NTCB]<<endl;
//  if( 11 == pdm ) { exit( -1 ); }


  // Update annual agricultural product pools and fluxes

//  if ( (CYCLE-1) == pdm ) { ag.updateyr( pdyr ); }

  if( totyr == modstartyr ) { wrtyr = 0;}

  if( (CYCLE-1) == pdm )
  {
    if( (totyr > modstartyr) ) { ++wrtyr; }
  }

  return wrtyr;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM60::natvegDynamics( const int& pdm, double pstate[] )
{
  int agstate = 0;
  int fertflag = 0;
  int perennial = 1;
  int tillflag = 0;
  int irrgflag = 0;

  double avlNH4;
  double avlNO3;

  double newpctp;
  double newsh2o;
  double newvsm;
  
  double propReact;


  soil.updateHydrology( elev,
                        atms.getTAIR(),
                        atms.getPREVTAIR(),
                        atms.getPREV2TAIR(),
                        atms.getRAIN(),
                        atms.getPET(),
                        pstate[I_SM],
                        pstate[I_RGRW],
                        pstate[I_SGRW],
                        irrgflag,
                        ag.irrigate,
                        pdm );
						
                                      
  // Assume no agricultural N fertilization (ag.fertn) and no
  //  nutrients resulting from agricultural conversion
  
  ag.fertn = ZERO;
  microbe.setAMMNVOL( ZERO );


  // Determine biological N fixation and partition it between
  //   symbiotic N fixation (veg.nfix) and nonsymbiotic N
  //   fixation (microbe.nfix)

  if( soil.getEET() > 0.1 )
  {
    bnfix = veg.getNFIXPARA( veg.cmnt ) * soil.getEET()
            + veg.getNFIXPARB( veg.cmnt );
  }
  else { bnfix = ZERO; }

  microbe.setNFIX( (bnfix * microbe.getNFIXPAR( veg.cmnt )) );
  
  if( microbe.getNFIX() < ZERO ) { microbe.setNFIX( ZERO ); }

  veg.setNFIX( (bnfix - microbe.getNFIX()) );

  if( veg.getNFIX() < ZERO ) { veg.setNFIX( ZERO ); }
  
  //if (veg.cmnt==7) cout <<"bnfix: "<<bnfix<<" vegnfix: "<<veg.getNFIX()<<" eet: "<<soil.getEET()<<endl;

  // Reduce EET if vegetation is not mature
  
  soil.setEVAPORATION( soil.getEET() 
                       * (1.0 - veg.getPROPTRANS( veg.cmnt)) );
					   
  veg.updateFoliage( veg.cmnt, pstate[I_FOLC],pstate[I_VEGC], soil.getEET() );

  soil.setEET( veg.getTRANSPIRATION() + soil.getEVAPORATION() );


  // Assume wetlands are wetter by the wfpsoff for determining
  //   moisture effects on vegetation and microbes

  newpctp = soil.getPCTP() + soil.getWFPSOFF();
  //if (veg.cmnt==7) cout <<"newpctp: "<<newpctp<<" getpctp: "<<soil.getPCTP()<<" wfpsoff: "<<soil.getWFPSOFF()<<endl;
  
  newsh2o = (newpctp * soil.getTOTPOR() * soil.getROOTZ()) 
		    / (100.0 * soil.getACTLAYER()); 

  newvsm = newsh2o / (soil.getROOTZ() * 1000.0);
  
  soil.setKH2O( newvsm, moistlim );

  // Get proportion of unfrozen organic matter in rooting zone
  
  propReact = soil.getThawedReactiveSOMProp( veg.cmnt );

  avlNH4 = pstate[I_NH4] * soil.getACTLAYER() / soil.getROOTZ();
  avlNO3 = pstate[I_NO3] * soil.getACTLAYER() / soil.getROOTZ();

  microbe.updateDynamics( veg.cmnt,
                          soil.getPCTFLDCAP(),
                          soil.getACTLAYER(),
                          (pstate[I_SOLC] * propReact),
                          pstate[I_AGR],
                          pstate[I_AGL],
                          pstate[I_BGR],
                          pstate[I_BGL],
                          pstate[I_CWD],
                          (pstate[I_SOLN] * propReact),
                          newsh2o,
                          newvsm,
                          avlNH4,
                          (atms.getNH4DEP()+ ag.getNRETENT() + ag.getFIRENDEP()),
                          moistlim,
                          tillflag,
                          ag.getTILLFACTOR( veg.cmnt ),
                          soil.getKH2O() );
  //closed by cgs to maintain vegetation flux during the month with disturbance
  /*if( disturbflag > 1 && pdm == (disturbmonth-1) )
  {
    // Set monthly vegetation fluxes to zero
    	
    veg.resetMonthlyFluxes();  	
  } 
  else*/
  {
    veg.updateDynamics( latitude,
    					longitude,
    					pdm,
    					veg.cmnt,
                        atms.getCO2(),
                        atms.getAOT40(),
						atms.getTAIR(),
                        atms.getNDEP(),
                        (ag.getNRETENT() + ag.getFIRENDEP()),
                        atms.getPAR(),
                        atms.getPET(),
                        atms.getPRVPETMX(),
                        soil.getEET(),
                        soil.getPRVEETMX(),
                        pstate[I_VEGC],
                        pstate[I_FOLC],
                        pstate[I_STEMC],
                        pstate[I_CROOTC],
                        pstate[I_FROOTC],
                        pstate[I_DEADWOODC],
                        pstate[I_STRN],
                        pstate[I_STON],
                        pstate[I_DEADWOODN],
                        newsh2o,
                        avlNH4,
                        avlNO3,
                        moistlim,
                        nfeed,
                        o3flag,
                        agstate,
                        perennial,
                        fertflag,
                        soil.getKH2O(),
                        microbe.getNETNMIN(),
                        microbe.getAMMNVOL(),
                        microbe.getNITRIF(),
                        microbe.getNO3PROD(),
                        ag.fertn,
                        pstate[I_AGE]);
  }                    

  //if (disturbflag > 1) cout <<"deadwoodc: "<<pstate[I_DEADWOODC] <<" woodfall rate: " <<veg.getwoodfall(veg.cmnt)<<" wooddecay: "<<veg.getdeadwoodltc()<<endl;
  // Determine carbon and nitrogen leaching losses
  //if (pstate[I_DEADWOODC] >0.0) cout <<"vegc: "<<pstate[I_VEGC]<<" deadwoodc: "<<pstate[I_DEADWOODC]<<" woodliter: "<<veg.deadwoodltc<<" woodfallp: "<<veg.getwoodfall(veg.cmnt)<<endl;
  soil.updateLeachingLosses( veg.cmnt,
                             (pstate[I_DOC] * propReact), 
                             (pstate[I_DON] * propReact), 
                             avlNO3, 
                             (pstate[I_SM] * soil.getACTLAYER() / soil.getROOTZ()) );
  

  if ( soil.getLEACHDOC() > (pstate[I_DOC]  
       + microbe.getDOCPROD()) )
  {
    soil.setLEACHDOC( (pstate[I_DOC] + microbe.getDOCPROD()) );
  }

  if ( soil.getLEACHDON() > (pstate[I_DON] 
         + microbe.getDONPROD()) )
  {
    soil.setLEACHDON( (pstate[I_DON] + microbe.getDONPROD()) );
  }

  
  // Determine loss of POC through erosion
  
  soil.setERODEPOC( ZERO );


  // Determine loss of PON through erosion
  
  soil.setERODEPON( ZERO );

               
  // Determine trace gas production based on microbe.no3prod,
  //   ag.fertn and water-filled pore space

  microbe.setTraceGasProduction( veg.cmnt,
                                 newpctp, 
                                 ag.fertn );

 
  if( soil.getLCHNO3PAR() > ZERO )
  {
    // Limit ecosystem N losses to total NO3 available for 
    //   leaching after denitrification has been substracted
    
    if ( soil.getLEACHNO3() > (avlNO3 + atms.getNO3DEP()
         + microbe.getNO3PROD() - microbe.getDENITRIF() 
         - veg.getNO3UPTAKE()) )
    {
      if ( (avlNO3 + atms.getNO3DEP() + microbe.getNO3PROD() 
           - microbe.getDENITRIF() - veg.getNO3UPTAKE()) > ZERO )
      {    
        soil.setLEACHNO3( (avlNO3
                           + atms.getNO3DEP()
                           + microbe.getNO3PROD()
                           - microbe.getDENITRIF()
                           - veg.getNO3UPTAKE()) );
      }
      else
      {
        soil.setLEACHNO3( ZERO );

        if( (avlNO3 + atms.getNO3DEP() + microbe.getNO3PROD() 
            - veg.getNO3UPTAKE()) > ZERO )
        {
          microbe.setDENITRIF( (avlNO3
                               + atms.getNO3DEP()
                               + microbe.getNO3PROD()
                               - veg.getNO3UPTAKE()) );
        }
        else
        {
          microbe.setDENITRIF( ZERO );
          microbe.setN2PROD( ZERO );
          microbe.setN2OPROD( ZERO );

          veg.setNO3UPTAKE( (avlNO3
                             + atms.getNO3DEP()
                             + microbe.getNO3PROD()) );
        }
      }
    }
  }
  else     
  {
    // If no leaching losses, limit ecosystem N losses to 
    //   total NO3 available for denitrification
    
    if( microbe.getDENITRIF() > (avlNO3 + atms.getNO3DEP() 
         + microbe.getNO3PROD() - veg.getNO3UPTAKE()) )
    {
      if( (avlNO3 + atms.getNO3DEP() + microbe.getNO3PROD() 
          - veg.getNO3UPTAKE()) > ZERO )
      {
        microbe.setDENITRIF( (avlNO3
                             + atms.getNO3DEP()
                             + microbe.getNO3PROD()
                             - veg.getNO3UPTAKE()) );
      }
      else
      {
        microbe.setDENITRIF( ZERO );
        microbe.setN2PROD( ZERO );
        microbe.setN2OPROD( ZERO );

        veg.setNO3UPTAKE( (avlNO3
                           + atms.getNO3DEP()
                           + microbe.getNO3PROD()) );
      }

                                 
      if( microbe.getN2OPROD() > ZERO )
      {
        microbe.setN2PROD( microbe.getDENITRIF() - microbe.getN2OPROD() );
      
        if( microbe.getN2PROD() < ZERO )
        {
          microbe.setN2PROD( ZERO );
          microbe.setN2OPROD( microbe.getDENITRIF() );
        }
      }
      else { microbe.setN2PROD( microbe.getDENITRIF() ); }           
    }    
  }


  // Determine trace gas fluxes based on ammonia volatilization
  //   (microbe.ammnvol), NO production (microbe.noprod),
  //   N2O production (microbe.n2oprod) and N2 production
  //   (microbe.n2prod)

  soil.setTraceGasFluxes( microbe.getAMMNVOL(),
                          microbe.getNOPROD(),
                          microbe.getN2OPROD(),
                          microbe.getN2PROD() );


  if( 0 == avlnflag )
  {
    // Do not allow changes in available nitrogen in soil
    //   (used during calibration process)

    // Equilibrate Ammonium (NH4) pool

    microbe.setAMMNVOL( (atms.getNH4DEP()
                        + microbe.getGROSSNMIN()
                        - microbe.getIMMNH4()
                        - veg.getNH4UPTAKE()
                        - microbe.getNITRIF()) );

    soil.setNH3FLUX( microbe.getAMMNVOL() );
    
    
    // Equilibrate nitrate (NO3) pool

    soil.setLEACHNO3( (atms.getNO3DEP()
                       + microbe.getNO3PROD()
                       - microbe.getDENITRIF()
                       - veg.getNO3UPTAKE()) );

    // Equilibrate DON pool

    soil.setLEACHDON( microbe.getDONPROD() );
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM60::pastureDynamics( const int& pdm, double pstate[] )
{
  int fertflag = 0;
  int irrgflag = 0;
  int perennial = 1;
  int tillflag = 0;

  double avlNH4;
  double avlNO3;

  double newpctp;
  double newsh2o;
  double newvsm;

  double propReact;
  
  soil.updateHydrology( elev,
                        atms.getTAIR(),
                        atms.getPREVTAIR(),
                        atms.getPREV2TAIR(),
                        atms.getRAIN(),
                        atms.getPET(),
                        pstate[I_SM],
                        pstate[I_RGRW],
                        pstate[I_SGRW],
                        irrgflag,
                        ag.irrigate,
                        pdm );

                                      
  // Assume no agricultural N fertilization (ag.fertn) and no
  //  nutrients resulting from agricultural conversion

  ag.fertn = ZERO;
  microbe.setAMMNVOL( ZERO );


  // Determine biological N fixation and partition it between
  //   symbiotic N fixation (veg.nfix) and nonsymbiotic N
  //   fixation (microbe.nfix)

  if( soil.getEET() > 0.1 )
  {
    bnfix = veg.getNFIXPARA( ag.cmnt ) * soil.getEET()
            + veg.getNFIXPARB( ag.cmnt );
  }
  else { bnfix = ZERO; }

  microbe.setNFIX( (bnfix * microbe.getNFIXPAR( veg.cmnt )) );
  
  if( microbe.getNFIX() < ZERO ) { microbe.setNFIX( ZERO ); }

  veg.setNFIX( (bnfix - microbe.getNFIX()) );

  if( veg.getNFIX() < ZERO ) { veg.setNFIX( ZERO ); }
  

  // Reduce EET if vegetation is not mature
  
  soil.setEVAPORATION( soil.getEET() 
                       * (1.0 - veg.getPROPTRANS( ag.cmnt)) );
					   
  veg.updateFoliage( ag.cmnt, pstate[I_FOLC], pstate[I_VEGC],soil.getEET() );

  soil.setEET( veg.getTRANSPIRATION() + soil.getEVAPORATION() );


  // Assume wetlands are wetter by the wfpsoff for determining
  //   moisture effects on vegetation and microbes

  newpctp = soil.getPCTP() + soil.getWFPSOFF();
  
  newsh2o = (newpctp * soil.getTOTPOR() * soil.getROOTZ()) 
		    / (100.0 * soil.getACTLAYER()); 

  newvsm = newsh2o / (soil.getROOTZ() * 1000.0);


  soil.setKH2O( newvsm, moistlim );

  // Get proportion of unfrozen organic matter in rooting zone
  
  propReact = soil.getThawedReactiveSOMProp( veg.cmnt );

  avlNH4 = pstate[I_NH4] * soil.getACTLAYER() / soil.getROOTZ();
  avlNO3 = pstate[I_NO3] * soil.getACTLAYER() / soil.getROOTZ();

  microbe.updateDynamics( veg.cmnt,
                          soil.getPCTFLDCAP(),
                          soil.getACTLAYER(),
                          (pstate[I_SOLC] * propReact),
                          pstate[I_AGR],
                          pstate[I_AGL],
                          pstate[I_BGR],
                          pstate[I_BGL],
                          pstate[I_CWD],
                          (pstate[I_SOLN] * propReact),
                          newsh2o,
                          newvsm,
                          avlNH4,
                          (atms.getNH4DEP()+ ag.getNRETENT() + ag.getFIRENDEP()),
                          moistlim,
                          tillflag,
                          ag.getTILLFACTOR( veg.cmnt ),
                          soil.getKH2O() );

  /*if( disturbflag > 1 && pdm == (disturbmonth-1) )
  {
    // Set monthly vegetation fluxes to zero
    	
    veg.resetMonthlyFluxes();  	
  } 
  else*/
  {
    veg.updateDynamics( latitude,
    					longitude,
    					pdm,
    					ag.cmnt,
                        atms.getCO2(),
                        atms.getAOT40(),
						atms.getTAIR(),
                        atms.getNDEP(),
                        (ag.getNRETENT() + ag.getFIRENDEP()),
                        atms.getPAR(),
                        atms.getPET(),
                        atms.getPRVPETMX(),
                        soil.getEET(),
                        soil.getPRVEETMX(),
                        pstate[I_VEGC],
                        pstate[I_FOLC],
                        pstate[I_STEMC],
                        pstate[I_CROOTC],
                        pstate[I_FROOTC],
                        pstate[I_DEADWOODC],
                        pstate[I_STRN],
                        pstate[I_STON],
                        pstate[I_DEADWOODN],
                        newsh2o,
                        avlNH4,
                        avlNO3,
                        moistlim,
                        nfeed,
                        o3flag,
                        ag.state,
                        perennial,
                        fertflag,
                        soil.getKH2O(),
                        microbe.getNETNMIN(),
                        microbe.getAMMNVOL(),
                        microbe.getNITRIF(),
                        microbe.getNO3PROD(),
                        ag.fertn,
                        pstate[I_AGE]);
  }                    

  // Determine carbon and nitrogen leaching losses
  
  soil.updateLeachingLosses( veg.cmnt,
                             (pstate[I_DOC] * propReact), 
                             (pstate[I_DON] * propReact), 
                             avlNO3, 
                             (pstate[I_SM] * soil.getACTLAYER() / soil.getROOTZ()) );
  

  if ( soil.getLEACHDOC() > (pstate[I_DOC]  
       + microbe.getDOCPROD()) )
  {
    soil.setLEACHDOC( (pstate[I_DOC] + microbe.getDOCPROD()) );
  }

  if ( soil.getLEACHDON() > (pstate[I_DON] 
         + microbe.getDONPROD()) )
  {
    soil.setLEACHDON( (pstate[I_DON] + microbe.getDONPROD()) );
  }

  
  // Determine loss of POC through erosion
  
  soil.setERODEPOC( ZERO );


  // Determine loss of PON through erosion
  
  soil.setERODEPON( ZERO );

               
  // Determine trace gas production based on microbe.no3prod,
  //   ag.fertn and water-filled pore space

  microbe.setTraceGasProduction( veg.cmnt,
                                 newpctp, 
                                 ag.fertn );
  
  
  if( soil.getLCHNO3PAR() > ZERO )
  {
    // Limit ecosystem N losses to total NO3 available for 
    //   leaching after denitrification has been substracted
    
    if ( soil.getLEACHNO3() > (avlNO3 + atms.getNO3DEP()
         + microbe.getNO3PROD() - microbe.getDENITRIF() 
         - veg.getNO3UPTAKE()) )
    {
      if ( (avlNO3 + atms.getNO3DEP() + microbe.getNO3PROD() 
           - microbe.getDENITRIF() - veg.getNO3UPTAKE()) > ZERO )
      {    
        soil.setLEACHNO3( (avlNO3
                           + atms.getNO3DEP()
                           + microbe.getNO3PROD()
                           - microbe.getDENITRIF()
                           - veg.getNO3UPTAKE()) );
      }
      else
      {
        soil.setLEACHNO3( ZERO );

        if( (avlNO3 + atms.getNO3DEP() + microbe.getNO3PROD() 
            - veg.getNO3UPTAKE()) > ZERO )
        {
          microbe.setDENITRIF( (avlNO3
                               + atms.getNO3DEP()
                               + microbe.getNO3PROD()
                               - veg.getNO3UPTAKE()) );
        }
        else
        {
          microbe.setDENITRIF( ZERO );
          microbe.setN2PROD( ZERO );
          microbe.setN2OPROD( ZERO );

          veg.setNO3UPTAKE( (avlNO3
                             + atms.getNO3DEP()
                             + microbe.getNO3PROD()) );
        }
      }
    }
  }
  else     
  {
    // If no leaching losses, limit ecosystem N losses to 
    //   total NO3 available for denitrification
    
    if( microbe.getDENITRIF() > (avlNO3 + atms.getNO3DEP() 
         + microbe.getNO3PROD() - veg.getNO3UPTAKE()) )
    {

      if( (avlNO3 + atms.getNO3DEP() + microbe.getNO3PROD() 
          - veg.getNO3UPTAKE()) > ZERO )
      {
        microbe.setDENITRIF( (avlNO3
                             + atms.getNO3DEP()
                             + microbe.getNO3PROD()
                             - veg.getNO3UPTAKE()) );
      }
      else
      {
        microbe.setDENITRIF( ZERO );
        microbe.setN2PROD( ZERO );
        microbe.setN2OPROD( ZERO );

        veg.setNO3UPTAKE( (avlNO3
                           + atms.getNO3DEP()
                           + microbe.getNO3PROD()) );
      }

                                 
      if( microbe.getN2OPROD() > ZERO )
      {
        microbe.setN2PROD( microbe.getDENITRIF() - microbe.getN2OPROD() );
      
        if( microbe.getN2PROD() < ZERO )
        {
          microbe.setN2PROD( ZERO );
          microbe.setN2OPROD( microbe.getDENITRIF() );
        }
      }
      else { microbe.setN2PROD( microbe.getDENITRIF() ); }           
    }    
  }


  // Determine trace gas fluxes based on ammonia volatilization
  //   (microbe.ammnvol), NO production (microbe.noprod),
  //   N2O production (microbe.n2oprod) and N2 production
  //   (microbe.n2prod)

  soil.setTraceGasFluxes( microbe.getAMMNVOL(),
                          microbe.getNOPROD(),
                          microbe.getN2OPROD(),
                          microbe.getN2PROD() );

  if( 0 == avlnflag )
  {
    // Do not allow changes in available nitrogen in soil
    //   (used during calibration process)

    // Equilibrate Ammonium (NH4) pool

    microbe.setAMMNVOL( (atms.getNH4DEP()
                        + microbe.getGROSSNMIN()
                        - microbe.getIMMNH4()
                        - veg.getNH4UPTAKE()
                        - microbe.getNITRIF()) );

    soil.setNH3FLUX( microbe.getAMMNVOL() );
    
    
    // Equilibrate nitrate (NO3) pool

    soil.setLEACHNO3( (atms.getNO3DEP()
                       + microbe.getNO3PROD()
                       - microbe.getDENITRIF()
                       - veg.getNO3UPTAKE()) );

    // Equilibrate DON pool

    soil.setLEACHDON( microbe.getDONPROD() );
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void TTEM60::pcdisplayDT( const double& tottime,
                          const double& deltat )
{

  window( 1,15,39,15 );
  gotoxy( 1,1 );
  printf( "TIME = %10.8lf   DT = %10.8lf    ",
          tottime,
          deltat );

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void TTEM60::pcdisplayMonth( const int& dyr, const int& pdm )
{
  char limitC[2];

  double outflux1;
  double outflux2;
  double outflux3;
  double outflux4;
  double outflux5;
  double outflux6;
  double outflux7;
  double outflux8;
  double outflux9;
  double outflux10;
  double outflux11;
  double outflux12;
  double outflux13;
  double outflux14;
  double outflux15;
  double outflux16;
  double outflux17;
  double outflux18;
  double outflux19;
  double outflux20;

  if ( 0 == topwind )
  {
    window( 1,1,80,24 );
    clrscr();
    switch( calwind )
    {

      case 1:  window( 1,1,43,1 );
               cout << "    TIME   AVAILW   RGRNDW  SNOWPCK SRGRNDW";
               break;

      case 2:  window( 1,1,32,1 );
               cout << " TIME   VEG C    SOIL C    DOC ";
               break;

      case 3:  window( 1,1,48,1 );
               cout << " TIME     STRN   STON  SOIL N    NH4     NO3     DON";
               break;

      case 4:  window( 1,1,57,1 );
               cout << "  TIME   TSOIL   DST0   DST5  DST10  DST20  DST50 DST100 DST200";
               break;

      case 5:  window( 1,1,57,1 );
               cout << "  TIME     NH4      NO3    NMIN   NITRIF   DNITRF  NLOST";
               break;

    }

    topwind = 1;
  }

  // Assign values of optional variables to outflux? for later 
  //   screen display

  outflux1 = getOptionalCflx( scy[0] );

  outflux2 = getOptionalCflx( scy[1] );

  outflux3 = getOptionalCflx( scy[2] );

  outflux4 = getOptionalCflx( scy[3] );

  outflux5 = getOptionalCflx( scy[4] );

  outflux6 = getOptionalCflx( scy[5] );

  outflux7 = getOptionalNflx( sny[0] );

  outflux8 = getOptionalNflx( sny[1] );

  outflux9 = getOptionalNflx( sny[2] );

  outflux10 = getOptionalWflx( swy[0] );

  outflux11 = getOptionalWflx( swy[1] );

  outflux12 = getOptionalWflx( swy[2] );

  outflux13 = getOptionalWflx( swy[3] );

  outflux14 = getOptionalWflx( swy[4] );

  outflux15 = getOptionalSoilTemp( ssty[0] );

  outflux16 = getOptionalSoilTemp( ssty[1] );

  outflux17 = getOptionalTraceGas( sgy[0] );

  outflux18 = getOptionalTraceGas( sgy[1] );

  outflux19 = getOptionalTraceGas( sgy[2] );

//  outflux20 = getOptionalTraceGas( sgy[3] );


  window( 1,2,80,13 );
  gotoxy( 1,1 );
  delline();
  gotoxy( 1,12 );

  // Display monthly values for selected C and N pools and fluxes

  switch( calwind )
  {
    case 1: printf( "%4d-%2d %9.2lf %8.2lf %8.2lf %7.2lf %6.1lf %6.1lf %6.1lf %6.1lf %7.2lf",
                    dyr,
                    (pdm+1),
                    soil.getAVLH2O(),
                    y[I_RGRW],
                    soil.getSNOWPACK(),
                    y[I_SGRW],
                    outflux10,
                    outflux11,
                    outflux12,
                    outflux13,
                    outflux14 );
            break;
            
    case 2: // Productivity is nitrogen limited

            if ( (y[I_INNPP] > y[I_NPP])
                 && (y[I_INNUP] == y[I_VNUP]) )
            {
              strcpy( limitC, "N" );
            }

            // Productivity is either carbon or nitrogen limited

            if ( (y[I_INNPP] == y[I_NPP])
                 && (y[I_INNUP] == y[I_VNUP]) )
            {
              strcpy( limitC, "E" );
            }

            // Productivity is carbon limited (climate)

            if ( (y[I_INNPP] == y[I_NPP])
                 && (y[I_INNUP] > y[I_VNUP]) )
            {
              strcpy( limitC, "C" );
            }

            // Productivity is limited by both carbon and nitrogen

            if ( (y[I_INNPP] > y[I_NPP])
                 && (y[I_INNUP] > y[I_VNUP]) )
            {
              strcpy( limitC, "B" );
            }

            // Unknown limits on productivity

            if ( y[I_INNPP] < y[I_NPP]
                 || y[I_INNUP] < y[I_VNUP] )
            {
              strcpy( limitC, "?" );
            }

            printf( "%3d-%2d %8.2lf %8.2lf %7.2lf %s %6.3lf %7.3lf %7.3lf %7.3lf %6.3lf %7.3lf",
                    dyr,
                    (pdm+1),
                    y[I_VEGC],
                    y[I_SOLC],
                    y[I_DOC],
                    limitC,
                    outflux1,
                    outflux2,
                    outflux3,
                    outflux4,
                    outflux5,
                    outflux6 );

            break;
    
    case 3: printf( "%3d-%2d %7.2lf %6.3lf %7.2lf %7.4lf %7.4lf %7.4lf %8.4lf %8.4lf %8.4lf",
                    dyr,
                    (pdm+1),
                    y[I_STRN],
                    y[I_STON],
                    y[I_SOLN],
                    y[I_NH4],
                    y[I_NO3],
                    y[I_DON],
                    outflux7,
                    outflux8,
                    outflux9 );

    case 4: printf( "%4d-%2d %6.1lf %6.1lf %6.1lf %6.1lf %6.1lf %6.1lf %6.1lf %6.1lf %6.1lf %6.1lf",
                    dyr,
                    (pdm+1),
                    soil.stm.getTSOIL(),
                    soil.stm.getDST0(),
                    soil.stm.getDST5(),
                    soil.stm.getDST10(),
                    soil.stm.getDST20(),
                    soil.stm.getDST50(),
                    soil.stm.getDST100(),
                    soil.stm.getDST200(),
                    outflux15,
                    outflux16 );
            break;
    case 5: printf("%4d-%2d %7.1lf %7.1lf %7.1lf %7.1lf %7.1lf %7.1lf %7.1lf %7.1lf %7.1lf",
                    dyr,
                    (pdm+1),
                    (y[I_NH4] * 1000.0),
                    (y[I_NO3] * 1000.0),
                    (y[I_NMIN] * 1000.0),
                    (y[I_NTRF] * 1000.0),
                    (y[I_DNTRF] * 1000.0),
                    (soil.getNLOST() * 1000.0),
                    (outflux17 * 1000.0),
                    (outflux18 * 1000.0),
                    (outflux19 * 1000.0) );
            break;

  }

  window( 1,14,80,14 );
  gotoxy( 1,1 );
  delline();
//  printf("                                                                         ");

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

#ifdef CALIBRATE_TEM

void TTEM60::pcdisplayODEerr( const int& test, double pstate[] )
{
  int errtest;

  window( 40,15,80,15 );
  gotoxy( 1,1 );
  if ( test != ACCEPT )
  {
    errtest = test - 1;
    switch( errtest )
    {
      case   I_VEGC: printf( "   VEGC = %8.2lf  Error = %11.8lf  ",
                             pstate[I_VEGC],
                             error[I_VEGC] );
                     break;

      case   I_SOLC: printf( "  SOILC = %8.2lf  Error = %11.8lf  ",
                             pstate[I_SOLC],
                             error[I_SOLC] );
                     break;

      case   I_STRN: printf( "   STRN = %8.2lf  Error = %11.8lf  ",
                             pstate[I_STRN],
                             error[I_STRN] );
                     break;

      case   I_STON: printf( "LABILEN = %8.3lf  Error = %11.8lf  ",
                             pstate[I_STON],
                             error[I_STON] );
                     break;

      case   I_SOLN: printf( "  SOILN = %8.2lf  Error = %11.8lf  ",
                             pstate[I_SOLN],
                             error[I_SOLN] );
                     break;

      case   I_NH4: printf( "    NH4 = %8.5lf  Error = %11.8lf  ",
                             pstate[I_NH4],
                             error[I_NH4] );
                     break;

      case   I_NO3: printf( "    NO3 = %8.5lf  Error = %11.8lf  ",
                             pstate[I_NO3],
                             error[I_NO3] );

      case  I_INGPP: printf( "INITGPP = %8.1lf  Error = %11.8lf  ",
                             pstate[I_INGPP],
                             error[I_INGPP] );
                     break;

      case    I_GPP: printf( "    GPP = %8.1lf  Error = %11.8lf  ",
                             pstate[I_GPP],
                             error[I_GPP] );
                     break;

      case  I_INNPP: printf( "INITNPP = %8.1lf  Error = %11.8lf  ",
                             pstate[I_INNPP],
                             error[I_INNPP] );
                     break;

      case    I_NPP: printf( "    NPP = %8.1lf  Error = %11.8lf  ",
                             pstate[I_NPP],
                             error[I_NPP] );
                     break;

      case    I_GPR: printf( "     RA = %8.1lf  Error = %11.8lf  ",
                             pstate[I_GPR],
                             error[I_GPR] );
                     break;

      case  I_RVMNT: printf( "     RM = %8.1lf  Error = %11.8lf  ",
                             pstate[I_RVMNT],
                             error[I_RVMNT] );
                     break;

      case  I_RVGRW: printf( "     RG = %8.1lf  Error = %11.8lf  ",
                             pstate[I_RVGRW],
                             error[I_RVGRW] );
                     break;

      case   I_LTRC: printf( "   LTRC = %8.1lf  Error = %11.8lf  ",
                             pstate[I_LTRC],
                             error[I_LTRC] );
                     break;

      case     I_RH: printf( "     RH = %8.1lf  Error = %11.8lf  ",
                             pstate[I_RH],
                             error[I_RH] );
                     break;

      case I_LCHDOC: printf( " LCHDOC = %8.1lf  Error = %11.8lf  ",
                             pstate[I_LCHDOC],
                             error[I_LCHDOC] );
                     break;


      case I_AGFRTN: printf( "AGFERTN = %8.3lf  Error = %11.8lf  ",
                              pstate[I_AGFRTN],
                              error[I_AGFRTN] );
                     break;

      case  I_SNFIX: printf( "  SNFIX = %8.3lf  Error = %11.8lf  ",
                             pstate[I_SNFIX],
                             error[I_SNFIX] );
                     break;

      case  I_ANFIX: printf( "  ANFIX = %8.3lf  Error = %11.8lf  ",
                             pstate[I_ANFIX],
                             error[I_ANFIX] );
                     break;

      case  I_INNUP: printf( "INUPTAK = %8.3lf  Error = %11.8lf  ",
                             pstate[I_INNUP],
                             error[I_INNUP] );
                     break;

      case I_INNH4UP: printf( "INNH4UP = %8.3lf  Error = %11.8lf  ",
                              pstate[I_INNH4UP],
                              error[I_INNH4UP] );
                      break;

      case I_INNO3UP: printf( "INNO3UP = %8.3lf  Error = %11.8lf  ",
                              pstate[I_INNO3UP],
                              error[I_INNO3UP] );
                      break;

      case   I_VNUP: printf( "NUPTAKE = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VNH4UP],
                             error[I_VNH4UP] );
                     break;

      case I_VNH4UP: printf( "  NH4UP = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VNH4UP],
                             error[I_VNH4UP] );
                     break;

      case I_VNO3UP: printf( "  NO3UP = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VNO3UP],
                             error[I_VNO3UP] );
                     break;

      case   I_VSUP: printf( "SUPTAKE = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VSUP],
                             error[I_VSUP] );
                     break;

      case   I_VLUP: printf( "LUPTAKE = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VLUP],
                             error[I_VLUP] );
                     break;

      case  I_VNMBL: printf( " NMOBIL = %8.3lf  Error = %11.8lf  ",
                             pstate[I_VNMBL],
                             error[I_VNMBL] );
                     break;

      case   I_LTRN: printf( "   LTRN = %8.3lf  Error = %11.8lf  ",
                             pstate[I_LTRN],
                             error[I_LTRN] );
                     break;

      case   I_GMIN: printf( "   GMIN = %8.3lf  Error = %11.8lf  ",
                             pstate[I_GMIN],
                             error[I_GMIN] );
                     break;

      case   I_NIMM: printf( "MCRONUP = %8.3lf  Error = %11.8lf  ",
                             pstate[I_NIMM],
                             error[I_NIMM] );
                     break;

      case   I_NMIN: printf( "   NMIN = %8.3lf  Error = %11.8lf  ",
                             pstate[I_NMIN],
                             error[I_NMIN] );
                     break;

      case   I_AMMN: printf( "AMMNVOL = %8.3lf  Error = %11.8lf  ",
                             pstate[I_AMMN],
                             error[I_AMMN] );
                     break;

      case   I_NTRF: printf( " NITRIF = %8.3lf  Error = %11.8lf  ",
                             pstate[I_NTRF],
                             error[I_NTRF] );
                     break;

      case   I_NO3P: printf( " NO3PRD = %8.3lf  Error = %11.8lf  ",
                             pstate[I_NO3P],
                             error[I_NO3P] );
                     break;

      case   I_NOP: printf( "  NOPRD = %8.3lf  Error = %11.8lf  ",
                             pstate[I_NOP],
                             error[I_NOP] );
                     break;

      case   I_N2OP: printf( " N2OPRD = %8.3lf  Error = %11.8lf  ",
                             pstate[I_N2OP],
                             error[I_N2OP] );
                     break;

      case    I_N2P: printf( "  N2PRD = %8.3lf  Error = %11.8lf  ",
                             pstate[I_N2P],
                             error[I_N2P] );
                     break;

      case  I_DNTRF: printf( "DNITRIF = %8.3lf  Error = %11.8lf  ",
                             pstate[I_DNTRF],
                             error[I_DNTRF] );
                     break;

      case I_NH3FLX: printf( " NH3FLX = %8.3lf  Error = %11.8lf  ",
                             pstate[I_NH3FLX],
                             error[I_NH3FLX] );
                     break;

      case  I_NOFLX: printf( "  NOFLX = %8.3lf  Error = %11.8lf  ",
                             pstate[I_NOFLX],
                             error[I_NOFLX] );
                     break;

      case I_N2OFLX: printf( " N2OFLX = %8.3lf  Error = %11.8lf  ",
                             pstate[I_N2OFLX],
                             error[I_N2OFLX] );
                     break;

      case  I_N2FLX: printf( "  N2FLX = %8.3lf  Error = %11.8lf  ",
                             pstate[I_N2FLX],
                             error[I_N2FLX] );
                     break;

      case I_LCHNO3: printf( " LCHNO3 = %8.1lf  Error = %11.8lf  ",
                             pstate[I_LCHNO3],
                             error[I_LCHNO3] );
                     break;

      case I_LCHDON: printf( " LCHDON = %8.1lf  Error = %11.8lf  ",
                             pstate[I_LCHDON],
                             error[I_LCHDON] );
                     break;

      case I_UNRMLF: printf( " UNLEAF = %8.3lf  Error = %11.8lf  ",
                             pstate[I_UNRMLF],
                             error[I_UNRMLF] );
                     break;

      case   I_LEAF: printf( "   LEAF = %8.3lf  Error = %11.8lf  ",
                             pstate[I_LEAF],
                             error[I_LEAF] );
                     break;

      case   I_RGRW: printf( " RGRNDW = %8.2lf  Error = %11.8lf  ",
                             pstate[I_RGRW],
                             error[I_RGRW] );
                     break;

      case   I_SGRW: printf( " SGRNDW = %8.2lf  Error = %11.8lf  ",
                             pstate[I_SGRW],
                             error[I_SGRW] );
                     break;

      case     I_SM: printf( " SMOIST = %8.5lf  Error = %11.8lf  ",
                             pstate[I_SM],
                             error[I_SM] );
                     break;

      case  I_RPERC: printf( "  RPERC = %8.3lf  Error = %11.8lf  ",
                             pstate[I_RPERC],
                             error[I_RPERC] );
                     break;

      case   I_RRUN: printf( "   RRUN = %8.1lf  Error = %11.8lf  ",
                             pstate[I_RRUN],
                             error[I_RRUN] );
                     break;

      case  I_SPERC: printf( "  SPERC = %8.1lf  Error = %11.8lf  ",
                             pstate[I_SPERC],
                             error[I_SPERC] );
                     break;

      case   I_SRUN: printf( "   SRUN = %8.1lf  Error = %11.8lf  ",
                             pstate[I_SRUN],
                             error[I_SRUN] );
                     break;

      case    I_EET: printf( "    EET = %8.3lf  Error = %11.8lf  ",
                             pstate[I_EET],
                             error[I_EET] );
                     break;
    }
  }

};

#endif

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TTEM60::resetMonthlyELMNTFluxes( void )
{

  // Reset all monthly fluxes to zero

  atms.resetMonthlyFluxes();
  
  veg.resetMonthlyFluxes();
  
  microbe.resetMonthlyFluxes();
  
  soil.resetMonthlyFluxes();
  
  ag.resetMonthlyFluxes();
  
  // Carbon fluxes 

  rsoil = ZERO;
  nep = ZERO;
  nce = ZERO;
  ntcb = ZERO;

  // Nitrogen fluxes 

  bnfix = ZERO;
  ntnb = ZERO;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM60::resetODEflux( void )
{
  int i;

  for( i = MAXSTATE; i < NUMEQ; ++i ) { y[i] = ZERO; }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM60::resetYrFluxes( void )
{

  atms.resetYrFluxes();
  
  veg.resetYrFluxes();
  
  microbe.resetYrFluxes();
  
  soil.resetYrFluxes();
  
  ag.resetYrFluxes();
  
  yrtotalc = ZERO;

  // Annual carbon fluxes

  yrrsoil = ZERO;
  yrnep = ZERO;
  yrnce = ZERO;
  yrntcb = ZERO;

  // Annual nitrogen fluxes
 
  yrbnfix = ZERO;
  yrntnb = ZERO;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM60::rkf( const int& numeq, 
                  double pstate[], 
                  double& pdt, 
                  const int& pdm )
{

  int i;
  double ptdt = ZERO;

  for( i = 0; i < numeq; ++i )
  {
    dum4[i] = dum5[i] = pstate[i];
    yprime[i] = rk45[i] = error[i] = ZERO;
  }

  ptdt = pdt * 0.25;

  delta( pdm, dum4, f11 );
  //if (f11[1]<0.0) cout <<"pstate1: "<< dum4[1]<<" pdstate1: "<<f11[1]<<endl;
  step( numeq, yprime, f11, yprime, a1 );

  step( numeq, rk45, f11, rk45, b1 );

  step( numeq, dum4, f11, ydum, ptdt );

  delta( pdm, ydum, f2 );
  //cout <<"MDC1: "<<pstate[I_MDC]<<" dum4: "<<dum4[I_MDC]<<" errmdc: "<<error[I_MDC]<<endl;
  
  for( i = 0; i < numeq; ++i ) 
  {
    f13[i] = a31*f11[i] + a32*f2[i];
  }
  
  step( numeq, dum4, f13, ydum, pdt );

  delta( pdm, ydum, f3 );

  step( numeq, yprime, f3, yprime, a3 );

  step( numeq, rk45, f3, rk45, b3 );
  
  for( i = 0; i < numeq; ++i )
  {
    f14[i] = a41*f11[i] + a42*f2[i] + a43*f3[i];
  }
  
  step( numeq, dum4, f14, ydum, pdt );

  delta( pdm, ydum, f4 );

  step( numeq, yprime, f4, yprime, a4 );

  step( numeq, rk45, f4, rk45, b4 );
  
  for( i = 0; i < numeq; ++i )
  {
    f15[i] = a51*f11[i] + a52*f2[i] + a53*f3[i] + a54*f4[i];
  }
  
  step( numeq, dum4, f15, ydum, pdt );

  delta( pdm, ydum, f5 );

  step( numeq, yprime, f5, yprime, a5 );

  step( numeq, rk45, f5, rk45, b5 );
  
  for( i = 0; i < numeq; ++i )
  {
    f16[i] = b61*f11[i] + b62*f2[i] + b63*f3[i] + b64*f4[i] + b65*f5[i];
  }
  //cout <<"MDC1.5: "<<pstate[I_MDC]<<" dum4: "<<dum4[I_MDC]<<" errmdc: "<<error[I_MDC]<<" pdt: "<<pdt<<endl;
  step( numeq, dum4, f16, ydum, pdt );

  delta( pdm, ydum, f6 );

  step( numeq, rk45, f6, rk45, b6 );

  step( numeq, dum4, yprime, dum4, pdt );

  step( numeq, dum5, rk45, dum5, pdt );
  
  for ( i = 0; i < numeq; ++i )
  {
    error[i] = fabs( dum4[i] - dum5[i] );
  }

};

/***************************************************************
 ***************************************************************/


/* *************************************************************
************************************************************** */

void TTEM60::setELMNTecd( const int& pdcmnt,
                          const double& psiplusc )
{
  // Initialize TEM parameters dependent upon a grid cell's
  //   soil texture

  veg.resetEcds( pdcmnt, psiplusc );
  
  microbe.resetEcds( pdcmnt, psiplusc );
  	 
  soil.resetEcds( pdcmnt );  
  
};

/* *************************************************************
************************************************************** */

/* *************************************************************
************************************************************** */

void TTEM60::setEquilC2N( const int& pdcmnt, const double& co2 )
{

  atms.yrpet = 1.0;
  soil.yreet = 1.0;

  // Determine vegetation C/N parameter as a function of
  //   vegetation type, annual PET, and annual EET (annual EET
  //   initially equal to yrpet)

  veg.updateC2N( pdcmnt,
                 soil.yreet,
                 atms.yrpet,
                 co2,
                 atms.getINITCO2() );

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM60::setEquilEvap( const double& nirr,
                           const double& tair,
                           const int& pdm )
{

  // Determine initial values for atms.prvpetmx, atms.prveetmx,
  //   and veg.topt

  atms.petjh( nirr, tair, pdm );

  if( 0 == pdm )
  {
    atms.setPRVPETMX( atms.getPET() );
    soil.setPRVEETMX( atms.getPET() );
    veg.setTOPT( tair );
  }
  else
  {
    if( atms.getPET() > atms.getPRVPETMX() )
    {
      atms.setPRVPETMX( atms.getPET() );
      soil.setPRVEETMX( atms.getPET() );
      veg.setTOPT( tair );
    }
  }
  

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM60::setPrevState( void )
{
  for( int i = 0; i < MAXSTATE; ++i ) { prevy[i] = y[i]; }

};

/* *************************************************************
************************************************************* */


/***************************************************************
 ***************************************************************/

void TTEM60::step( const int& numeq, 
                   double pstate[], 
                   double pdstate[], 
                   double ptstate[],
	           double& pdt )
{
  for( int i = 0; i < numeq; ++i )
  {
    ptstate[i] = pstate[i] + (pdt * pdstate[i]);
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TTEM60::stepmonth( const int& pdyr,
                       const int& pdm,
                       int& intflag,
                       const double& ptol,
                       ofstream& rflog1 )
{

//  int i;
  int mintflag;

  // Reset all monthly fluxes to zero
  
  resetMonthlyELMNTFluxes();
  
  // Reset ODE state variables representing monthly
  //   fluxes to zero
  
  resetODEflux();
  for( int i = 0; i < NUMEQ; i++ )
  {
	  f111[i]=0.0;
  }


  //calculate new standage
  if (disturbflag == 0) y[I_AGE] = ag.getage();
  else if (disturbflag > 0 && pdm <(disturbmonth-1) ) y[I_AGE] = prevy[I_AGE];
  else if (disturbflag >0 && pdm >= (disturbmonth-1) ) y[I_AGE] = ag.getage();

  //if (veg.cmnt == 4 && longitude ==-80.75 && latitude == 51.25 && pdyr <200) cout <<"lat: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" I_NCE: "<<I_NCE<<" I_NTCB: "<<I_NTCB<<" NCE: "<<y[I_NCE]<<" NTCB: "<<y[I_NTCB]<<endl;

  // If 12 months have passed since disturbance, reset
  //   all immediate fluxes associated with disturbance 
  //   (i.e., slash, convrtflx, nretent) to zero
  ag.resetMonthlyDisturbFluxes();
  /*
  if( CYCLE == distmnthcnt )
  {
    distmnthcnt = 0;
    //ag.resetMonthlyDisturbFluxes(); //closed by cgs2014. gas emissions will only occur in the disturbance month rather than 12 months after disturbance
  } 

  // Count the number of months since disturbance
  
  if( distmnthcnt > 0 ) { ++distmnthcnt; }


  // If (FRI times 12 months have passed since a fire
  //   disturbance, reset ag.firendep to zero
  
  if( (ag.getFRI() * 12) == firemnthcnt )
  {
    firemnthcnt = 0;
    ag.setFIRENDEP( ZERO );
  }

  // Count the number of months since a fire disturbance
  
  if( firemnthcnt > 0 ) { ++firemnthcnt; };
  
  if( 0 == ag.state && 2 == distmnthcnt )
  {
     // Establish natural vegetation the month
     //   following disturbance if cohort is not
     //   in agriculture

      prevy[I_VEGC] = y[I_VEGC] = ag.getNATSEEDC();
      prevy[I_STRN] = y[I_STRN] = ag.getNATSEEDSTRN();
      prevy[I_STON] = y[I_STON] = ag.getNATSEEDSTON();  
  }
  */
  if( 0 == pdm )
  {
    if( 0 == pdyr )
    {
      microbe.setKD( microbe.getKDC() );
      ag.setKD( microbe.getKD() );
      ag.setNATSOIL( y[I_SOLC] );
    }
    else
    {
      if( 0 == ag.state && 0 == ag.prvstate )
      {
        microbe.setKD( microbe.yrkd( nfeed,
                                     veg.yrltrc,
                                     veg.yrltrn,
                                     veg.getLCCLNC( veg.cmnt ) ) );
        ag.setKD( microbe.getKD() );
        ag.setNATSOIL( y[I_SOLC] );
      }
      else
      {
        if( y[I_SOLC] < ag.getNATSOIL() )
        {
          microbe.setKD( (ag.getKD() 
                         * y[I_SOLC] 
                         / ag.getNATSOIL()) );
        }
        else { microbe.setKD( ag.getKD() ); }
      }
    }


    // Reset annual fluxes to zero

    resetYrFluxes();
  }

  //calculated month drought code here. Source: Girardin and Wotton. 2009. Journal of Applied Meteorology and climatology.

  int monday[12]={31,28,31,30,31,30,31,31,30,31,30,31};
  double lfn[12]={-1.6,-1.6,-1.6,0.9,3.8,5.8,6.4,5.0,2.4,0.4,-1.6,-1.6};//north hemis
  double lfs[12]={6.4,5.0,2.4,0.4,-1.6,-1.6,-1.6,-1.6,-1.6,0.9,3.8,5.8}; //south hemisphere
  //calculate daylength factor (LF)
  double LF;
  if (latitude>15.0) LF = lfn[pdm];
  else if (latitude <= 15.0 && latitude > -15.0) LF=1.39;
  else if (latitude <=-15.0) LF=lfs[pdm];
  //calculate potential ET
  double temptair = atms.getTAIR();
  double tempmaxtair=temptair + 2.8; //maximum temperature is estimated to be noon observation + 2.8 oC
  if (tempmaxtair <0.0) tempmaxtair = 0.0;

  double potet = monday[pdm] * (0.36 * atms.getMXTAIR() + LF);
  if (potet <0.0) potet =0.0;
  //calculate first half month drought code
  double MDChalf = prevy[I_MDC] + 0.25 * potet;
  //calculated end month drought code (MDCm)
  double RMeff = 0.83 *atms.getPREC();
  double Qmr = 800.0 * exp(MDChalf/(-400.0)) + 3.937 * RMeff;
  if (Qmr >800.0) Qmr = 800.0;
  if (Qmr <=0.0) Qmr = 0.001;
  double MDCm = 400 * log(800/Qmr) + 0.25 * potet;
  y[I_MDC] = (prevy[I_MDC] + MDCm) / 2.0;
  if (y[I_MDC]>800.0) y[I_MDC]=800.0;
  if (y[I_MDC]<=0.0) y[I_MDC]=1.00;


  //if (longitude == -127.75 && latitude ==55.25 && disturbflag ==3 && (pdm==7||pdm==8) ) cout <<" pdm-1: "<<pdm<<" totvegc: "<<totvegc<<" totvegc2: "<< y[I_CROOTC] + y[I_FROOTC]+y[I_STEMC]+y[I_FOLC]<<" yvegc: "<<y[I_VEGC]<<" prevyvegc: "<<prevy[I_VEGC]<<" soilc: "<<y[I_SOLC]<<endl;


  //cout <<" lat: "<<latitude<<" lon: "<<longitude<<" pdm: "<<pdm<<" prec: "<<atms.getPREC()<<" MDC: "<<y[I_MDC]<<" fuel: "<<y[I_AGL]+y[I_AGL]<< " FFFC: "<<FFFC<<" FFFC2: "<<FFFC2<<endl;
  //cout <<" lat: "<<latitude<<" lon: "<<longitude<<" pdm: "<<pdm<<" prec: "<<atms.getPREC()<<" maxtair: "<<tempmaxtair<<" LF: "<<LF<<" potet: "<<potet<<" dchalf: "<<MDChalf<<" Qmr: "<<Qmr<<" MDCm: "<<MDCm<<" MDC: " <<y[I_MDC]<<" prevMDC: "<<prevy[I_MDC]<<endl;

  // Implement disturbance effects

  if( disturbflag > 0 && pdm == (disturbmonth-1) )
  {
    distmnthcnt = 1;

    // Save proportion of vegetation for "seed" for regrowth
    // after disturbance

    double leftcarbon = y[I_FOLC] * (1-ag.leafmortpar) + y[I_STEMC] * (1-ag.stemmortpar) + (y[I_FROOTC] + y[I_CROOTC]) * (1-ag.rootmortpar);
    double leftstrn = leftcarbon * y[I_STRN]/omitzero(y[I_VEGC]); //added by cgs. need to improve in the later
    double leftston = leftcarbon * y[I_STON]/omitzero(y[I_VEGC]);

    //ag.setNATSEEDC( (y[I_VEGC] * ag.getVRESPAR()) );
    //ag.setNATSEEDSTRN( (y[I_STRN] * ag.getVRESPAR()) );
    //ag.setNATSEEDSTON( (y[I_STON] * ag.getVRESPAR()) );
    ag.setNATSEEDC(leftcarbon);
    ag.setNATSEEDSTRN(leftstrn);
    ag.setNATSEEDSTON(leftston);

      
    // Determine carbon and nitrogen lost during the first year
    //   after disturbance

    int ifwetland =0;
    int cveg = veg.getCURRENTVEG();
    if (cveg== 3 || cveg== 5 || cveg== 7|| cveg== 11 || cveg== 17 || cveg== 20 || cveg== 22|| cveg== 23 || cveg== 24 || cveg== 25 || cveg== 26|| cveg== 28|| cveg== 29|| cveg== 30) ifwetland = 1;


    ag.conversion( veg.cmnt,
                   y[I_VEGC],
                   y[I_FOLC],
                   y[I_STEMC],
                   y[I_CROOTC] + y[I_FROOTC],
                   y[I_DEADWOODC],
                   y[I_AGL] + y[I_AGR],
                   y[I_CWD],
                   y[I_MDC],
                   y[I_DEADWOODN],
                   y[I_STRN],
                   y[I_STON],
                   y[I_SOLC],
                   y[I_SOLN],
                   disturbflag,
                   ifwetland,
                   randomnum);
    //if (veg.cmnt ==7) cout <<"deadwoodc: "<< y[I_DEADWOODC]<<endl;

    // Create potential 10-year and 100-year products
    //   from human disturbance

    ag.createWoodProducts( pdyr,
                           y[I_VEGC],
                           y[I_STEMC],
                           y[I_STRN],
                           y[I_STON] );

    //if (ag.sconvert >0.0) cout <<"sconvert: "<<ag.sconvert<<endl;
    //prevy[I_VEGC] = y[I_VEGC] = ZERO; //assign initial vegetation carbon pools for new cohort after disturbance
    //prevy[I_STRN] = y[I_STRN] = ZERO;
    //prevy[I_STON] = y[I_STON] = ZERO;

    //added by cgs2014 to update the carbon/nitrogen pools after disturbance

    y[I_FOLC] = y[I_FOLC] * (1-ag.leafmortpar);
    y[I_STEMC] = y[I_STEMC] * (1-ag.stemmortpar);
	y[I_CROOTC] = y[I_CROOTC] * (1-ag.rootmortpar);
	y[I_FROOTC] = y[I_FROOTC] * (1-ag.rootmortpar);
	y[I_VEGC]  = leftcarbon;
	y[I_STRN]  = leftstrn;
	y[I_STON]  = leftston;
	y[I_DEADWOODC] = y[I_DEADWOODC]
	                   +ag.getdeadwoodc()
	                   -ag.getswfc()
			 	 	   -ag.getdeadslashc();
	y[I_DEADWOODN] = y[I_DEADWOODN]
	                   +ag.getdeadwoodn()
	                   - ag.getswfc() * y[I_DEADWOODN]/omitzero(y[I_DEADWOODC])
				 	   - ag.getdeadslashn();
	y[I_SOLC] = y[I_SOLC]
	            + ag.getSLASHC()
	            - ag.getSCONVRTFLXC();


	  if (veg.getifwoody(veg.cmnt)==1) //forest
	  {
	  y[I_CWD] = y[I_CWD]
			  	  	   +ag.getstem2slash() * 0.7 //disturbed stem cwd
			  	  	   +ag.getdeadslashc() * 0.8
	                   -ag.getdwfc();
	  y[I_AGR] = y[I_AGR]
			  	  	   + ag.getstem2slash()*0.3//disturbed stem litter
			  	  	   + ag.getdeadslashc()*0.2
			  	  	   -y[I_AGR]* ag.getfflc() / omitzero(y[I_AGR] + y[I_AGL]) ;
	  y[I_AGL] = y[I_AGL]
	                   +ag.getleaf2slash()//disturbed leaf litter
	                   -ag.getfflc() * y[I_AGL] / omitzero(y[I_AGR] + y[I_AGL]) ;
	  y[I_BGR] = y[I_BGR]
	                   +ag.getroot2slash()*0.8; //disturbed root

	  y[I_BGL] = y[I_BGL]
	                   +ag.getroot2slash()*0.2; //disturbed root
	  }

	  if (veg.getifwoody(veg.cmnt)==2) //shrub
	  {
		  y[I_CWD] = y[I_CWD]
	                   + ag.getstem2slash() * 0.5//disturbed stem litter
	                   -ag.getdwfc();
		  y[I_AGR] = y[I_AGR]
				       + ag.getstem2slash()*0.5//disturbed stem litter
				       -ag.getfflc() * y[I_AGR] / omitzero(y[I_AGR] + y[I_AGL]);
		  y[I_AGL] = y[I_AGL]
		                   +ag.getleaf2slash()//disturbed leaf litter
		                   -ag.getfflc() * y[I_AGL] / omitzero(y[I_AGR] + y[I_AGL]) ;
		  y[I_BGR] = y[I_BGR]
		                   +ag.getroot2slash()*0.8; //disturbed root

		  y[I_BGL] = y[I_BGL]
		                   +ag.getroot2slash()*0.2; //disturbed root
	  }
	  if (veg.getifwoody(veg.cmnt)==0) //herbaceous
	  {
		  y[I_CWD] = y[I_CWD]
	                   + ag.getstem2slash() * 0.0//disturbed stem litter
	                   -ag.getdwfc();
		  y[I_AGR] = y[I_AGR]
				       + ag.getstem2slash()*1.0//disturbed stem litter
				       -ag.getfflc() * y[I_AGR] / omitzero(y[I_AGR] + y[I_AGL]);
		  y[I_AGL] = y[I_AGL]
		                   +ag.getleaf2slash()//disturbed leaf litter
		                   -ag.getfflc() * y[I_AGL] / omitzero(y[I_AGR] + y[I_AGL]) ;
		  y[I_BGR] = y[I_BGR]
			  	  	   +ag.getroot2slash()*0.3;
		  y[I_BGL] = y[I_BGL]
	                   +ag.getroot2slash()*0.7;
	  }


	  y[I_SOLN] = y[I_SOLN]
	                    + ag.getSLASHN()
	                    + ag.getdeadslashn()
	                    - ag.getSCONVRTFLXN()
	                    - ag.getNSRETENT();

	//if (longitude == -127.75 && latitude ==55.25 && disturbflag ==3 && (pdm==7||pdm==8) ) cout <<" pdm0: "<<pdm<<" totvegc: "<< y[I_CROOTC] + y[I_FROOTC]+y[I_STEMC]+y[I_FOLC]<<" yvegc: "<<y[I_VEGC]<<" soilc: "<<y[I_SOLC]<<" deadwoodc: "<<y[I_DEADWOODC]<<" litrc: "<<veg.getLTRFALC()<<endl;

	if (disturbflag==1)//no standing deadwood and CWD after land use change
    {
    //reset deadwoodc and deadwoodN as 0
      y[I_SOLC]+=y[I_DEADWOODC];//the left deadwood will be added to soil
      y[I_SOLN]+=y[I_DEADWOODN];
      y[I_DEADWOODC] = ZERO;
      y[I_DEADWOODN] = ZERO;

      prevy[I_DEADWOODC] = ZERO;
      prevy[I_DEADWOODN] = ZERO;

    }
    if( 3 == disturbflag )   // fire disturbance = 3
    {
//      ag.setFireNDEP();
      ag.setFIRENDEP( ZERO );
      firemnthcnt = 1;
    }
      
	  
    if( 1 == disturbflag )
    {

      if( 1 == ag.state)
      {
        // Update soil characteristics for rooting depth that
        //   is appropriate to crops
        soil.updateRootZ( ag.cmnt );


        // Establish crops

//        prevy[I_VEGC] = y[I_VEGC] = ag.cropseedC[ag.cmnt];
//        prevy[I_STRN] = y[I_STRN] = ag.cropseedSTRN[ag.cmnt];
//        prevy[I_STON] = y[I_STON] = ag.cropseedSTON[ag.cmnt];


        // Update soil texture-dependent vegetation parameters
        //   for crops

        veg.resetEcds( ag.cmnt, soil.getPSIPLUSC() );

        // Update other adaptive parameters

        atms.yrpet = 1.0;
        soil.yreet = 1.0;
      
        veg.updateC2N( ag.cmnt,
                       soil.yreet,
                       atms.yrpet,
                       atms.getPREVCO2(),
                       atms.getINITCO2() );

        veg.setPRVLEAFMX( ag.getCROPPRVLEAFMX() );
     
        veg.setTOPT( ag.getCROPTOPT() );
      
        atms.setPRVPETMX( ag.getCROPPRVPETMX() );
      
        soil.setPRVEETMX( ag.getCROPPRVEETMX() );
      
        ag.setPRVCROPNPP( ZERO );
      }
      
      if( 2 == ag.state )
      {
        // Update soil characteristics for rooting depth that
        //   is appropriate to pasture

        soil.updateRootZ( ag.cmnt );

        // Establish pasture

//        prevy[I_VEGC] = y[I_VEGC] = ag.cropseedC[ag.cmnt];
//        prevy[I_STRN] = y[I_STRN] = ag.cropseedSTRN[ag.cmnt];
//        prevy[I_STON] = y[I_STON] = ag.cropseedSTON[ag.cmnt];


        // Update soil texture-dependent vegetation parameters
        //   for crops

        veg.resetEcds( ag.cmnt, soil.getPSIPLUSC() );

        // Update other adaptive parameters

        atms.yrpet = 1.0;
        soil.yreet = 1.0;
      
        veg.updateC2N( ag.cmnt,
                       soil.yreet,
                       atms.yrpet,
                       atms.getPREVCO2(),
                       atms.getINITCO2() );

        veg.setPRVLEAFMX( ag.getCROPPRVLEAFMX() );
     
        veg.setTOPT( ag.getCROPTOPT() );
      
        atms.setPRVPETMX( ag.getCROPPRVPETMX() );
      
        soil.setPRVEETMX( ag.getCROPPRVEETMX() );
      
        ag.setPRVCROPNPP( ZERO );
      }      

      if( 3 == ag.state )
      {
        // Update rooting depth to be appropriate to urban areas

        soil.updateRootZ( ag.cmnt );

        // Update soil texture-dependent vegetation parameters
        //   for urban areas

        veg.resetEcds( ag.cmnt, soil.getPSIPLUSC() );

        // Update other adaptive parameters

        atms.yrpet = 1.0;
        soil.yreet = 1.0;

        veg.updateC2N( ag.cmnt,
                       soil.yreet,
                       atms.yrpet,
                       atms.getPREVCO2(),
                       atms.getINITCO2() );

        veg.setPRVLEAFMX( ag.getCROPPRVLEAFMX() );

        veg.setTOPT( ag.getCROPTOPT() );

        atms.setPRVPETMX( ag.getCROPPRVPETMX() );

        soil.setPRVEETMX( ag.getCROPPRVEETMX() );

        ag.setPRVCROPNPP( ZERO );
      }
    }
  }
  //end of disturbflag > 0

//  else { ag.setNoWoodProducts( pdyr ); }

  // Revert to natural vegetation after cropland abandonment

  if( 0 == ag.state && ag.prvstate > 0 )
  {
    // Update soil characteristics to rooting depth that is 
    //   appropriate to natural vegetation

    soil.updateRootZ( veg.cmnt );


    // Establish natural vegetation

    prevy[I_VEGC] = y[I_VEGC] = ag.getNATSEEDC();
    prevy[I_STRN] = y[I_STRN] = ag.getNATSEEDSTRN();
    prevy[I_STON] = y[I_STON] = ag.getNATSEEDSTON();

    // Update soil texture-dependent vegetation parameters
    //   for natural vegetation

    veg.resetEcds( veg.cmnt, soil.getPSIPLUSC() );


    // Update other adaptive parameters

    atms.yrpet = ag.getNATYRPET();
    soil.yreet = ag.getNATYREET();

    veg.updateC2N( veg.cmnt,
                   soil.yreet,
                   atms.yrpet,
                   atms.getPREVCO2(),
                   atms.getINITCO2() );

    veg.setPRVLEAFMX( ag.getNATPRVLEAFMX() );
    veg.setTOPT( ag.getNATTOPT() );
    atms.setPRVPETMX( ag.getNATPRVPETMX() );
    soil.setPRVEETMX( ag.getNATPRVEETMX() );
    ag.setPRVCROPNPP( ZERO );
  }
  

  // set ozone efect from previous month for initial month

  if( 0 == pdyr && 0 == pdm ) { veg.setFPREVOZONE( 1.0 ); }


  // Get environmental conditions for month "dm"
  //if ((pdyr == 3 || pdyr == 4) && pdm ==6) cout <<"pdyr: "<<pdyr<<" getSNOWPACK1: "<<soil.getSNOWPACK()<<"avtive1: " <<soil.getACTLAYER()<<endl;

  getenviron( pdyr, pdm, rflog1 );
  //if (pdm ==2) cout <<"pdyr: "<<pdyr<<" getSNOWPACK1: "<<soil.getSNOWPACK()<<" snowpack2: "<<y[I_SNWPCK]<<" avtive1: " <<soil.getACTLAYER()<<" active2: " <<y[I_ACTLAYER] <<" I_snowpack: " <<I_SNWPCK<<" I_VEGC:  "<<I_VEGC<<endl;

  // Determine effect of air temperature on GPP (i.e. temp)
  //   and GPR (i.e. respq10)
  
  if( 0 == ag.state )
  {
    veg.setTEMP( veg.cmnt, atms.getTAIR() );
	veg.setNewRESPQ10( veg.cmnt, atms.getTAIR() );
  }
  else
  {
	veg.setTEMP( ag.cmnt, atms.getTAIR() );
	veg.setNewRESPQ10( ag.cmnt, atms.getTAIR() );
  }

  //if (ag.state>0) cout <<"pdyr: "<<pdyr<<" agstate: "<<ag.state<<" ag.cmnt: "<<ag.cmnt<<" veg.cmnt: "<<veg.cmnt<<endl;
	// Determine effect of temperature on decomposition 
    //   (i.e. dq10)

	microbe.setNewDQ10( veg.cmnt,
						atms.getTAIR(),
						soil.getTSOIL(),
						soil.getSNOWPACK(),
                        soil.stmflg );

  // Update growing degree days

  if( atms.getTAIR() >= GDDMIN )
  {
    ag.setGROWDD( ag.getGROWDD() 
                  + ((atms.getTAIR() - GDDMIN) 
                  * atms.getNDAYS( pdm )) );
    //if (veg.getCURRENTVEG()==50 && pdyr >3000) cout <<" dyr: "<<pdyr<<" mon: "<<pdm<<" gddcum: "<< ag.getGROWDD()<<endl;
  }
  else
  {
    // If "cold snap" (i.e. TAIR < GDDMIN) hits after crops 
    //   begin to grow, crops are assumed to die and resulting
    //   detritus is added to soil organic carbon and nitrogen

    if( 1 == ag.state && ag.getGROWDD() >= GDDSEED  )
    {       	
      ag.frostDamage( y[I_VEGC], y[I_STRN], y[I_STON] );
      y[I_VEGC] = ZERO;
      prevy[I_VEGC] = ZERO;
      y[I_FOLC] = ZERO;
      prevy[I_FOLC] = ZERO;
      y[I_STEMC] = ZERO;
      prevy[I_STEMC] = ZERO;
      y[I_FROOTC] = ZERO;
      prevy[I_FROOTC] = ZERO;
      y[I_CROOTC] = ZERO;
      prevy[I_CROOTC] = ZERO;

      y[I_STRN] = ZERO;
      prevy[I_STRN] = ZERO;
      y[I_STON] = ZERO;
      prevy[I_STON] = ZERO;
      y[I_SOLC] += ag.getSTUBBLEC();
      prevy[I_SOLC] = y[I_SOLC];
      y[I_SOLN] += ag.getSTUBBLEN();
      prevy[I_SOLN] = y[I_SOLN];
      //if (veg.getCURRENTVEG()==50) cout <<" dyr: "<<pdyr<<" mon: "<<pdm<<endl;
    }

    ag.setGROWDD( ZERO );
  }

//  if ( 191 == pdyr ) { dbugflg = 1; }

  // Run TEM for a monthly time step

  //if (veg.cmnt == 4 && longitude ==-80.75 && latitude == 51.25 && pdyr <200) cout <<"lat0: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" disturbmon: "<<disturbmonth<<" currentveg: "<<veg.getCURRENTVEG()<<" vegc: "<<y[I_VEGC] <<" soilc: " <<y[I_SOLC]<<" standdead: "<<y[I_DEADWOODC]<<" stemc: "<<y[I_STEMC]<<" rootc: "<<y[I_CROOTC] + y[I_FROOTC]<<" abovel: "<<y[I_AGR] + y[I_AGL] <<" cwd: "<<y[I_CWD]<<" ntcb: "<<y[I_NTCB]<<" convert0: "<<ag.getCONVRTFLXC()<<" convertc: "<<y[I_CNVRTC]<<" rh: "<< microbe.getRH()<<endl;
  //if (veg.cmnt == 7 && longitude ==-93.5 && latitude == 41.25 && veg.getCURRENTVEG()==50 && pdyr <10) cout <<"lat1: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" disturbmon: "<<disturbmonth<<" currentveg: "<<veg.getCURRENTVEG()<<" vegc: "<<y[I_VEGC] <<" soilc: " <<y[I_SOLC]<<" stemc: "<<y[I_STEMC]<<" rootc: "<<y[I_CROOTC] + y[I_FROOTC]<<" ntcb: "<<y[I_NTCB]<<" npp: "<<y[I_NPP]<<" nce: "<<y[I_NCE]<<" rh: "<<y[I_RH]<<" strn: "<<y[I_STRN]<<" ston: "<<y[I_STON]<<" litn: "<<veg.getLTRFALN()<<" sup: "<<veg.getSUPTAKE() <<" resorbn: "<<veg.getNRESORB()<<" nmobil: "<<veg.getNMOBIL()<<" nfix: "<<veg.getNFIX()<<endl;
  for( int i = 0; i < MAXSTATE; i++ )
  {
	  if (y[i]<0.0) y[i]=0.0;
	  if (prevy[i]<0.0) prevy[i]=0.0;
  }
  y[I_SOC]=y[I_SOLC]-y[I_AGR]-y[I_AGL]-y[I_BGR]-y[I_BGL]-y[I_CWD];
  /*//functions: adapt() was closed by cgs to avoid the problem for carbon imbalance

  mintflag = adapt( NUMEQ, y, ptol, pdm );

  if( 1 == mintflag ) { intflag = 1; }

  if ( blackhol != 0 )
  {
    if( 0 == initFlag ) { qualcon[0] = 10; }
    else { qualcon[pdyr] += 10; }
  }

  else*/
  {
		  delta( pdm, y, f111);
		  //if (veg.cmnt==4&&pdyr<10&& longitude ==-89.75 && latitude == 51.5) cout <<" soilc1: "<<y[I_SOC]<<" prevsoilc: "<<prevy[I_SOC]<<" f111: "<<f111[I_SOC]<<" SOLC: "<<f111[I_SOLC]<<microbe.getSOIL_LITTER()-microbe.getRH_SOIL()<<endl;
		  for( int i = 0; i < NUMEQ; ++i )
		  {
		  	  y[i]=f111[i]+y[i];
		  }
		  //step( NUMEQ, y, y, f111, tempvara1);
  }

  for( int i = 0; i < MAXSTATE; i++ )
  {
	  if (y[i]<0.0) y[i]=0.0;
	  if (prevy[i]<0.0) prevy[i]=0.0;
  }
 // if (veg.cmnt==4&&pdyr<10&& longitude ==-89.75 && latitude == 51.5) cout <<"equil: "<<equil<<" soilc2: "<<y[I_SOC]<<endl;
  //if (veg.cmnt == 7 && longitude ==-93.5 && latitude == 41.25 && veg.getCURRENTVEG()==50 && pdyr <10) cout <<"lat2: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" disturbmon: "<<disturbmonth<<" currentveg: "<<veg.getCURRENTVEG()<<" vegc: "<<y[I_VEGC] <<" soilc: " <<y[I_SOLC]<<" stemc: "<<y[I_STEMC]<<" rootc: "<<y[I_CROOTC] + y[I_FROOTC]<<" ntcb: "<<y[I_NTCB]<<" npp: "<<y[I_NPP]<<" nce: "<<y[I_NCE]<<" rh: "<<y[I_RH]<<" strn: "<<y[I_STRN]<<" ston: "<<y[I_STON]<<" litn: "<<veg.getLTRFALN()<<" sup: "<<veg.getSUPTAKE() <<" resorbn: "<<veg.getNRESORB()<<" nmobil: "<<veg.getNMOBIL()<<" nfix: "<<veg.getNFIX()<<endl;
   //if (veg.cmnt == 7 && veg.getCURRENTVEG()==50) cout <<"lat2: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" no3: "<< y[I_NO3]<<" nh4: "<<y[I_NH4]<<" strn: "<<y[I_STRN]<<" ston: "<<y[I_STON]<<" litn: "<<veg.getLTRFALN()<<" sup: "<<veg.getSUPTAKE() <<" resorbn: "<<veg.getNRESORB()<<" nmobil: "<<veg.getNMOBIL()<<" nfix: "<<veg.getNFIX()<<endl;

  //if (ag.getCONVRTFLXC() > 0.0) cout << "dm: " <<pdm<<"soil c2: " <<y[I_SOLC]<<" sconvertc: " <<y[I_SCNVRTC]<<" mintflag: "<<mintflag<<endl;

 //if (disturbflag ==3 && (pdm==7||pdm==8) && longitude == -127.75 && latitude ==55.25) cout <<" equil: "<<equil<<" pdm1: "<<pdm<<" totvegc: "<< y[I_CROOTC] + y[I_FROOTC]+y[I_STEMC]+y[I_FOLC]<<" yvegc: "<<y[I_VEGC]<<" soilc: "<<y[I_SOLC]<<" deadwoodc: "<<y[I_DEADWOODC]<< " litrc: "<<veg.getLTRFALC()<<endl;
 //
  // Check mass balance

  //massbal( y, prevy );


  // Determine veg.plant.nitrogen
  
  veg.setVEGN( (y[I_STRN] + y[I_STON]) );
  
  // Determine soil.availn.total
  
  soil.setAVLN( (y[I_NH4] + y[I_NO3]) );

  // Determine water yield (soil.h2oyld)

  soil.setH2OYLD( (y[I_RRUN] + y[I_SRUN]) );


  // Determine total monthly N inputs to ecosystem

  soil.setNINPUT( (ag.getNRETENT()
                   + ag.getFIRENDEP() 
                   + atms.getTOTNDEP()
                   + y[I_AGFRTN] 
                   + y[I_BNFIX]) );


  // Determine monthly Soil Respiration (rsoil)

//  rsoil = y[I_RH] + y[I_ROOTGPR] - y[I_CO2DISS];
  rsoil = y[I_RH] + y[I_ROOTGPR];


  // Determine Net Ecosystem Production (nep)

  nep = y[I_NPP] - y[I_RH];

//  cout << nep << " " << endl;

  // Determine fluxes from crop residues

  ag.updateCropResidueFluxes();


  // Determine fluxes from decay of products

  ag.decayProducts();


  // Reset growing degree days to zero if crops were
  // harvested this month

  if( 1 == ag.state && ag.getGROWDD() >= GDDHARVST )
  {
    ag.harvest( pdm, y[I_VEGC], y[I_STRN], y[I_STON] );
    y[I_VEGC] = ZERO;
    y[I_STRN] = ZERO;
    y[I_STON] = ZERO;
    y[I_SOLC] += ag.getSTUBBLEC();
    y[I_SOLN] += ag.getSTUBBLEN();
    y[I_FOLC] = ZERO;
    y[I_STEMC] = ZERO;
    y[I_FROOTC] = ZERO;
    y[I_CROOTC] = ZERO;
    y[I_AGR]+=0.1 * ag.getSTUBBLEC() * 0.5; //the litter carbon will be added to different pools after harvest. added by cgs2014
    y[I_AGL]+=0.9 * ag.getSTUBBLEC() * 0.5;
    y[I_BGR]+=0.1 * ag.getSTUBBLEC() * 0.5;
    y[I_BGL]+=0.9 * ag.getSTUBBLEC() * 0.5;
//if (veg.getCURRENTVEG()==50) cout <<" dyr: "<<pdyr<<" mon: "<<pdm<<" gddharvest: "<<GDDHARVST<<" gddcum: "<< ag.getGROWDD()<<endl;
    ag.setGROWDD( ZERO );

  }
  else { ag.setNoCrops( pdm ); }

  // Determine crop residue
  
  ag.updateCropResidue();


  // Determine standing stock of products

  ag.updateProducts();

  ag.updateTotalProductFormation();

//  cout << nep << " " << endl;

  // Graze veg biomass every month if cohort is in pasture

  if( 2 == ag.state )
  {
   ag.grazing(  y[I_VEGC], y[I_STRN] );
    y[I_VEGC] -= ag.getFORAGEC();
    y[I_STRN] -= ag.getFORAGEN();
    y[I_SOLC] += ag.getMANUREC();
    y[I_SOLN] += ag.getMANUREN();
    y[I_NH4] += ag.getURINE();
  }
  else { ag.setNoGrazing(); }


  // Determine CFLUX from ecosystem from NEP plus
  //   fluxes from burning associated with agricultural
  //   conversion or

  ag.setCFLUX( (nep
                - ag.getCONVRTFLXC()
                - ag.getCROPRESIDUEFLXC()
                - ag.getANIMALRESP()) );

  // Determine Net Carbon Exchange (NCE) with atmosphere
  //   (CFLUX plus decay of agricultural and wood products)
  
  nce = ag.getCFLUX() - ag.getTOTPRODDECAYC();



  // Determine Net Terrestrial Carbon Balance (NTCB)
  //   of ecosystem
  
//  y[I_ERDPOC] = ZERO;
//  y[I_LCHCO2] = ZERO;
//  y[I_LCHHCO3] = ZERO;
  
//  ntcb = nce 
//         - y[I_LCHDOC]
//         - y[I_ERDPOC]
//         - y[I_LCHCO2]
//         - y[I_LCHHCO3];

  ntcb = nce - y[I_LCHDOC];


  y[I_NTCB] = ntcb;
  y[I_NCE] = nce;
  y[I_NEP] =nep;
  y[I_CNVRTC]=ag.getCONVRTFLXC();
  if (y[I_NCE] <0.0001 && y[I_NCE] >-0.0001) y[I_NCE] = 0.0;
  if ( y[I_NTCB] <0.0001 &&  y[I_NTCB] >-0.0001)  y[I_NTCB] = 0.0;
  if (y[I_NEP] <0.0001 && y[I_NEP] >-0.0001) y[I_NEP] = 0.0;

  //if (veg.cmnt == 4 && longitude ==-80.75 && latitude == 51.25 && pdyr <200) cout <<"lat: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" disturbmon: "<<disturbmonth<<" currentveg: "<<veg.getCURRENTVEG()<<" vegc: "<<y[I_VEGC] <<" soilc: " <<y[I_SOLC]<<" standdead: "<<y[I_DEADWOODC]<<" stemc: "<<y[I_STEMC]<<" rootc: "<<y[I_CROOTC] + y[I_FROOTC]<<" abovel: "<<y[I_AGR] + y[I_AGL] <<" cwd: "<<y[I_CWD]<<" ntcb: "<<y[I_NTCB]<<" npp: "<<y[I_NPP]<<" convertc: "<<y[I_CNVRTC]<<" nce: "<<y[I_NCE]<<" rh: "<<y[I_RH]<<" MDC: "<<y[I_MDC]<<endl;
  //if (veg.cmnt == 7 && longitude ==-93.5 && latitude == 41.25 && veg.getCURRENTVEG()==50 && pdyr <200) cout <<"lat3: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" disturbmon: "<<disturbmonth<<" currentveg: "<<veg.getCURRENTVEG()<<" vegc: "<<y[I_VEGC] <<" soilc: " <<y[I_SOLC]<<" stemc: "<<y[I_STEMC]<<" rootc: "<<y[I_CROOTC] + y[I_FROOTC]<<" ntcb: "<<y[I_NTCB]<<" npp: "<<y[I_NPP]<<" nce: "<<y[I_NCE]<<" rh: "<<y[I_RH]<<" strn: "<<y[I_STRN]<<" ston: "<<y[I_STON]<<" litn: "<<veg.getLTRFALN()<<" sup: "<<veg.getSUPTAKE() <<" nmobil: "<<veg.getNMOBIL()<<endl;

  //if (pdm >=10) cout <<"mon1: " <<pdm<<" totsolc: "<<soil.getTSOLC()<<" I_SOLC: " <<y[I_SOLC]<<endl;
  // Determine carbon storage in ecosystem

  ag.setTOTEC( (y[I_VEGC]
             + y[I_SOLC]
             + soil.getNSOLC()
             + y[I_DOC]
             + y[I_DEADWOODC] ) );


  // Determine total carbon in ecosystems plus
  //   agricultural and wood products
                 
  totalc = ag.getTOTEC() + ag.getTOTPRODC(); 
 

  // Add nitrogen released from products back into nitrate pool
  //closed by cgs2014. No product nitrogen goes back to the ecosystem. all product nitrogen will be lost to atmosphere or the landfill
  //y[I_NO3] +=  ag.getTOTPRODDECAYN();


  // Determine total loss of nitrogen from ecosystem
    
  soil.setNLOST( (y[I_NH3FLX]
                  + y[I_NOFLX]
                  + y[I_N2OFLX]
                  + y[I_N2FLX]
                  + y[I_LCHNO3]
                  + y[I_LCHDON]
                  + y[I_ERDPON]
                  + ag.getCONVRTFLXN()
                  + ag.getCROPRESIDUEFLXN()
                  + ag.getTOTPRODDECAYN()) );// added by cgs2014

  // Determine net terrestrial nitrogen balance (NTNB) 
  //   in ecosystem
      
  ntnb = soil.getNINPUT() - soil.getNLOST();


  // Update ANNUAL carbon, nitrogen and water pools and fluxes 
  //   from integrator results

  updateYearSummary();


  #ifdef CALIBRATE_TEM
    // Display monthly results to DOS screen
    pcdisplayMonth( pdyr, pdm );
  #endif


  if( 1 == ag.state && ag.getPRVCROPNPP() < y[I_NPP] )
  {
     ag.setPRVCROPNPP( y[I_NPP] );
  }
  else { ag.setPRVCROPNPP( ZERO ); }

  if( 1 == ag.state && ag.getPRVCROPNPP() > y[I_NPP] )
  {
    y[I_UNRMLF] = ZERO;
  }

  // Reset growing degree days to zero if crops were
  // harvested this month

  if( 1 == ag.state && ag.getGROWDD() >= GDDHARVST )
  {
    ag.setGROWDD( ZERO );
  }

  if( atms.getTAIR() < GDDMIN ) { ag.setGROWDD( ZERO ); }


//  if ( dyr > 0 || dm > 0 ) { veg.fprevozone = veg.fozone[dm]; }
  if( pdyr > 0 || pdm > 0 )
  {
    veg.setFPREVOZONE( (1.0 - 0.5 * (1.0 - y[I_FOZONE])) );
  }

  if( y[I_NPP] <= ZERO ) { veg.setFPREVOZONE( 1.0 ); }


  // Update atms.prevco2 for next month

  atms.setPREVCO2( atms.getCO2() );
  
  // Update atms.prevtair and atms.prev2tair for next month

  atms.setPREV2TAIR( atms.getPREVTAIR() );
  atms.setPREVTAIR( atms.getTAIR() );

  atms.setPREVRAIN( atms.getTAIR() );
  // Update previous snowpack and previous soil temperature
  // at 10 cm for next month

  soil.setPREVSPACK( soil.getSNOWPACK() );
  soil.setPREVDST10( soil.getDST10() );

  // Update previous unnormalized relative leaf area for
  //   next month
  
  veg.setPREVUNRMLEAF( y[I_UNRMLF] );
  
  // Update ag.prevPROD1, ag.prevPROD10 and ag.prevPROD100
  // for next month

  ag.setPREVPROD1C( ag.getPROD1C() );
  ag.setPREVPROD1N( ag.getPROD1N() );

  ag.setPREVPROD10C( ag.getPROD10C() );
  ag.setPREVPROD10N( ag.getPROD10N() );

  ag.setPREVPROD100C( ag.getPROD100C() );
  ag.setPREVPROD100N( ag.getPROD100N() );

  // Update ag.prevCropResidue for next month

  ag.setPREVCROPRESIDUEC( ag.getCROPRESIDUEC() );
  ag.setPREVCROPRESIDUEN( ag.getCROPRESIDUEN() );

  //  Update maximum EET, maximum PET, GPP optimum temperature
  //    (veg.topt), and maximum leaf cover (veg.prvleafmx) for
  //    the current year

  if( 0 == pdm )
  {
    soil.setEETMX( y[I_EET] );
    atms.setPETMX( atms.getPET() );
    veg.setNEWTOPT( atms.getTAIR() );
    veg.setNEWLEAFMX( y[I_UNRMLF] );
  }
  else
  {
    if( y[I_EET] > soil.getEETMX() )
    {
      soil.setEETMX( y[I_EET] );
    }

    if( atms.getPET() > atms.getPETMX() )
    {
      atms.setPETMX( atms.getPET() );
    }

    if( 0 == ag.state )
    {
      veg.resetNEWTOPT( veg.cmnt, 
                        atms.getTAIR(), 
                        y[I_UNRMLF] );
    }
    else 
    {
      veg.resetNEWTOPT( ag.cmnt, 
                        atms.getTAIR(), 
                        y[I_UNRMLF] );
    }    
  }
  /*
  //total soil carbon: TOTSOLC or I_TSOLC. CGS move below codes from the bottom to here
  if( 0 == initFlag )
  {
    soil.setNSOLC( y[I_SOLC]
                   / (1.0 - soil.getNSOLPAR( veg.cmnt ))
                   - y[I_SOLC] );

    soil.setNSOLN( ((y[I_SOLN]
                   / (1.0 - soil.getNSOLPAR( veg.cmnt )))
                   - y[I_SOLN]) );
  }
  soil.setTSOLC( (y[I_SOLC] + soil.getNSOLC() + y[I_DOC]) );
  soil.setTSOLN( (y[I_SOLN] + soil.getNSOLN() + y[I_DON]) );
  if (pdm ==11 && pdyr <10) cout <<"mon2: " <<pdm<<" totsolc: "<<soil.getTSOLC()<<" I_SOLC: " <<y[I_SOLC]<<endl;
  */

  // Save state of all the ODE state variables 
  //   representing pools to allow checking
  //   of mass balance
  
  setPrevState();


  // Update annual parameters for next year
  
  if( (CYCLE-1) == pdm )
  {
    soil.setPRVEETMX( soil.getEETMX() );
    atms.setPRVPETMX( atms.getPETMX() );
    veg.setTOPT( veg.getNEWTOPT() );
    veg.setPRVLEAFMX( veg.getNEWLEAFMX() );

    // Update optimum temperature parameters for GPP

    if( 0 == ag.state )
    {
      veg.boundTOPT( veg.cmnt );
      
    // Update adaptive parameters

      ag.setNATYRPET( atms.yrpet );
      ag.setNATYREET( soil.yreet );
      ag.setNATPRVPETMX( atms.getPRVPETMX() );
      ag.setNATPRVEETMX( soil.getPRVEETMX() );
      ag.setNATTOPT( veg.getTOPT() );
      ag.setNATPRVLEAFMX( veg.getPRVLEAFMX() );

      // Determine vegetation C/N parameter as a function
      // of vegetation type, annual PET, annual EET,
      // CO2 concentration

      veg.updateC2N( veg.cmnt,
                     soil.yreet,
                     atms.yrpet,
                     atms.getPREVCO2(),
                     atms.getINITCO2() );


    }
    else
    {
      veg.boundTOPT( ag.cmnt );

      ag.setCROPPRVPETMX( atms.getPRVPETMX() );
      ag.setCROPPRVEETMX( soil.getPRVEETMX() );
      ag.setCROPTOPT( veg.getTOPT() );
      ag.setCROPPRVLEAFMX( veg.getPRVLEAFMX() );

     // Determine vegetation C/N parameter as a function of
     //   vegetation type, annual PET, annual EET, CO2
     //   concentration

      veg.updateC2N( ag.cmnt,
                     soil.yreet,
                     atms.yrpet,
                     atms.getPREVCO2(),
                     atms.getINITCO2() );

    }

    // Update next year ag.prvstate with current year ag.state
    //if (longitude==-98.5&&latitude==29.5&&pdyr>4990&&veg.cmnt==5) cout <<" pdyr: "<<pdyr<<" tair: "<<atms.getTAIR()<<" temp: "<<veg.getTEMP()<<" q10: "<<veg.getQ10()<<endl;

    ag.prvstate = ag.state;

    veg.yrcarbon  /= 12.0;
    soil.yrorgc /= 12.0;
    soil.yrnonorgc /= 12.0;
    soil.yrtotorgc /= 12.0;
    soil.yrDOM.carbon /= 12.0;

    soil.yrgasCO2 /= 12.0;
    soil.yraqCO2 /= 12.0;
    soil.yrHCO3 /= 12.0;
    soil.yrRHCO3 /= 12.0;
    soil.yralk /= 12.0;

    yrtotalc /= 12.0;
    veg.yrdeadwoodc /= 12.0;
    veg.yrdeadwoodn /= 12.0;

    veg.yrnitrogen  /= 12.0;
    veg.yrstructn /= 12.0;

    if( veg.yrstructn != ZERO )
    {
      veg.yrc2n  = veg.yrcarbon / veg.yrstructn;
    }
//  else { veg.yrc2n = MISSING; }

    veg.yrstoren /= 12.0;

    soil.yrorgn /= 12.0;

    if( soil.yrorgn != ZERO )
    {
      soil.yrc2n = soil.yrorgc / soil.yrorgn;
    }
//  else { soil.yrc2n = MISSING; }

    soil.yrnonorgn /= 12.0;
    soil.yrtotorgn /= 12.0;

    soil.yrDOM.nitrogen /= 12.0;

    soil.yravln.total  /= 12.0;
    soil.yravln.nh4  /= 12.0;
    soil.yravln.no3  /= 12.0;

    soil.yravlh2o /= 12.0;
    soil.yrsmoist /= 12.0;
    soil.yrvsm /= 12.0;
    soil.yrpctp /= 12.0;
    soil.yrsnowpack /= 12.0;
    soil.yrrgrndh2o /= 12.0;
    soil.yrsgrndh2o /= 12.0;

    veg.yrunleaf /= 12.0;
    veg.yrleaf /= 12.0;
    veg.yrlai /= 12.0;
    veg.yrfpc /= 12.0;

    soil.yrtsoil /= 12.0;
    soil.stm.yrdst0 /= 12.0;
    soil.stm.yrdst5 /= 12.0;
    soil.yrdst10 /= 12.0;
    soil.stm.yrdst20 /= 12.0;
    soil.stm.yrdst50 /= 12.0;
    soil.stm.yrdst100 /= 12.0;
    soil.stm.yrdst200 /= 12.0;
    soil.stm.yrfrontd /= 12.0;
    soil.stm.yrfrontd2 /= 12.0;
    soil.stm.yrthawbegin1 /= 12.0;
    soil.stm.yrthawbegin2 /= 12.0;
    soil.stm.yrthawend1 /= 12.0;
    soil.stm.yrthawend2 /= 12.0;
    veg.yrthawpct /= 12.0;

//    if( 1 == baseline )
//    {
//      y[I_SOLN] = y[I_SOLC] / microbe.getCNSOIL( veg.cmnt );
//    }

    if( 0 == initFlag )
    {
      soil.setNSOLC( ((y[I_SOLC]
                     / (1.0 - soil.getNSOLPAR( veg.cmnt )))
                     - y[I_SOLC]) );

      soil.setNSOLN( ((y[I_SOLN]
                     / (1.0 - soil.getNSOLPAR( veg.cmnt )))
                     - y[I_SOLN]) );
    }

    soil.setTSOLC( (y[I_SOLC] + soil.getNSOLC() + y[I_DOC]) );

    soil.setTSOLN( (y[I_SOLN] + soil.getNSOLN() + y[I_DON]) );
    //if (pdm ==11 && pdyr <10) cout <<"mon2: " <<pdm<<" totsolc: "<<soil.getTSOLC()<<" I_SOLC: " <<y[I_SOLC]<<" nonsolcpar: "<<soil.getNSOLPAR( veg.cmnt )<<endl;

    if( endeq > 0 ) 
    { 
      ++endeq; 
    }
  }
  //if (veg.cmnt == 4 && longitude ==-80.75 && latitude == 51.25 && pdyr <200) cout <<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" vegc: "<<y[I_VEGC] <<" soilc: " <<y[I_SOLC]<<" standdead: "<<y[I_DEADWOODC]<<" stemc: "<<y[I_STEMC]<<" rootc: "<<y[I_CROOTC] + y[I_FROOTC]<<" cwd: "<<y[I_CWD]<<" abovel: "<<y[I_AGR] + y[I_AGL] <<" vconvert: "<<ag.getVCONVRTFLXC() <<" sconvert: "<<ag.getSCONVRTFLXC()<<" fflc: "<<ag.getfflc()<<" ffdc: "<<ag.getffdc()<<" swfc: "<<ag.getswfc()<<" dwfc: "<<ag.getdwfc()<<" deadconvc: "<<ag.getdeadconvc()<<" MDC: "<<y[I_MDC]<<" GPP: "<<y[I_GPP] <<" NPP: "<<y[I_NPP] <<" RH: "<<y[I_RH]<<" NTCB: "<<y[I_NTCB]<<endl;
  //if (veg.cmnt == 4 && longitude ==-80.75 && latitude == 51.25 && pdyr <1000 && pdyr >499) cout <<" pdyr: " <<pdyr<<" pdm2: "<<pdm <<" vegc: "<<y[I_VEGC] <<" soilc: " <<y[I_SOLC]<<" standdead: "<<y[I_DEADWOODC]<<" stemc: "<<y[I_STEMC]<<" rootc: "<<y[I_CROOTC] + y[I_FROOTC]<<" cwd: "<<y[I_CWD]<<" abovel: "<<y[I_AGR] + y[I_AGL] <<" vconvert: "<<ag.getVCONVRTFLXC() <<" sconvert: "<<ag.getSCONVRTFLXC()<<" fflc: "<<ag.getfflc()<<" ffdc: "<<ag.getffdc()<<" swfc: "<<ag.getswfc()<<" dwfc: "<<ag.getdwfc()<<" deadconvc: "<<ag.getdeadconvc()<<" MDC: "<<y[I_MDC]<<" GPP: "<<y[I_GPP] <<" NPP: "<<y[I_NPP] <<" RH: "<<y[I_RH]<<" NTCB: "<<y[I_NTCB]<<endl;
  //if (veg.cmnt == 7 && veg.getCURRENTVEG()==50) cout <<"lat3: "<<latitude<<" lon: "<<longitude<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" no3: "<< y[I_NO3]<<" nh4: "<<y[I_NH4]<<" strn: "<<y[I_STRN]<<" ston: "<<y[I_STON]<<" litn: "<<veg.getLTRFALN()<<" sup: "<<veg.getSUPTAKE() <<" resorbn: "<<veg.getNRESORB()<<" nmobil: "<<veg.getNMOBIL()<<" nfix: "<<veg.getNFIX()<<" decayn: "<<ag.getTOTPRODDECAYN()<<endl;

  return endeq;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TTEM60::testEquilibrium( void )
{

  if( 0 == nfeed && 0 == rheqflag
      && wtol >= fabs( atms.yrrain + atms.yrsnowfall + ag.yrirrig //added by cgs to test equilibrium for cropland
      - soil.yreet - soil.yrrrun - soil.yrsrun )
      //&& (ctol >= fabs( veg.yrnpp - veg.yrltrc )) )
      && (ctol >= fabs( yrntcb )) )//changed by cgs for reach equilibrium in cropland

  {
    return 1;
  }

  if( 0 == nfeed && 1 == rheqflag
      && wtol >= fabs( atms.yrrain 
                       + atms.yrsnowfall + ag.yrirrig
                       - soil.yreet 
                       - soil.yrrrun 
                       - soil.yrsrun )
                        
      && (ctol >= fabs( yrntcb ))
      //&& (ctol >= fabs( veg.yrnpp - veg.yrltrc ))//closed by cgs for reach equilibrium in cropland
      && (ctol >= fabs( veg.yrltrc + ag.yrstubC- microbe.yrdecomp )) //veg.deadwoodltc = 0, ag.getSCONVRTFLXC()=0, ag.getslashhc =0
       
      && (ctol >= fabs( microbe.yrdecomp 
                        - microbe.yrrh 
                        - microbe.yrDOMprod.carbon )) )
  {
    return 1;
  }

  if( 1 == nfeed && 1 == rheqflag
      && wtol >= fabs( atms.yrrain 
                       + atms.yrsnowfall + ag.yrirrig
                       - soil.yreet 
                       - soil.yrrrun 
                       - soil.yrsrun )
              
      && (ntol >= fabs( yrntnb - soil.yrabNimmob ))

      && (ntol >= fabs( veg.yrnup.total 
                        + veg.yrnfix 
                        - veg.yrltrn - ag.yrstubN ))
                         
      && (ntol >= fabs( veg.yrltrn + ag.yrstubN
                        + microbe.yrnfix 
                        - microbe.yrndecomp ))
                         
      && (ntol >= fabs( microbe.yrndecomp 
                        - microbe.yrnmin 
                        - microbe.yrDOMprod.nitrogen ))
        
      && (ntol >= fabs( atms.yrndep.nh4 
                        + microbe.yrnmin 
                        - soil.yrabNimmob 
                        - veg.yrnup.nh4
                        - microbe.yrnitrif 
                        - microbe.yrammnvol ))
                          
      && (ntol >= fabs( microbe.yrnitrif 
                        - microbe.yrno3prd 
                        - microbe.yrnoprd ))
                           
      && (ntol >= fabs( microbe.yrno3prd 
                        + atms.yrndep.no3 
                        - veg.yrnup.no3 
                        - soil.yrlchNO3
                        - microbe.yrdenitrf ))
                         
      && (ntol >= fabs( microbe.yrdenitrf 
                        - microbe.yrn2oprd 
                        - microbe.yrn2prd ))
                         
      && (ntol >= fabs( microbe.yrDOMprod.nitrogen 
                        - soil.yrlchDOM.nitrogen ))

      && (ntol >= fabs( microbe.yrammnvol - soil.yrnh3flx ))
      && (ntol >= fabs( microbe.yrnoprd - soil.yrnoflx ))
      && (ntol >= fabs( microbe.yrn2oprd - soil.yrn2oflx ))
      && (ntol >= fabs( microbe.yrn2prd - soil.yrn2flx ))

      && (ctol >= fabs( yrntcb ))
      //&& (ctol >= fabs( veg.yrnpp - veg.yrltrc )) //closed by cgs for reach equilibrium in cropland
      && (ctol >= fabs( veg.yrltrc + ag.yrstubC- microbe.yrdecomp ))
       
      && (ctol >= fabs( microbe.yrdecomp 
                        - microbe.yrrh 
                        - microbe.yrDOMprod.carbon )) )
  {

//    if( (soil.yrabNimmob < 0.0 
//        && (ntol >= fabs( yrntnb - soil.yrabNimmob )))
//        || (soil.yrabNimmob >= 0.0
//        && (ntol >= fabs( yrntnb ))) )
//    {   
      return 1;
//    }
  }

//  cout << "NTCB = " << yrntcb << endl;
  return 0;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM60::updateYearSummary( void )
{

  // Update sum of annual carbon storage in ecosystems

  veg.yrcarbon  += y[I_VEGC];
  soil.yrorgc += y[I_SOLC];
  soil.yrDOM.carbon += y[I_DOC];
  veg.yrdeadwoodc+=y[I_DEADWOODC];
  veg.yrdeadwoodn+=y[I_DEADWOODN];
//  soil.yrgasCO2 += y[I_CO2G];
//  soil.yraqCO2 += y[I_CO2W];
//  soil.yrHCO3 += y[I_HCO3];
//  soil.yrRHCO3 += y[I_RHCO3];

  soil.yrnonorgc += soil.getNSOLC();
  soil.yrtotorgc += soil.getTSOLC();
//  soil.yralk += soil.alkalinity;
  yrtotalc += totalc;

  // Update sum of annual nitrogen storage in ecosystems

  veg.yrstructn += y[I_STRN];
  veg.yrstoren += y[I_STON];
  soil.yrorgn += y[I_SOLN];
  soil.yrDOM.nitrogen += y[I_DON];
  soil.yravln.nh4  += y[I_NH4];
  soil.yravln.no3  += y[I_NO3];

  veg.yrnitrogen  += veg.getVEGN();
  soil.yrnonorgn += soil.getNSOLN();
  soil.yrtotorgn += soil.getTSOLN();
  soil.yravln.total  += soil.getAVLN();

  // Update sum of annual water storage in ecosystems

  soil.yravlh2o += ((y[I_SM] * soil.getACTLAYER() / soil.getROOTZ())
                   - soil.getWILTPT());
                   
  soil.yrsmoist += y[I_SM];
  
  soil.yrvsm += (y[I_SM] / (soil.getROOTZ() * 1000.0)) ;
  
  soil.yrpctp += (100.0 * (y[I_SM] * soil.getACTLAYER() / soil.getROOTZ())
                 / soil.getTOTPOR());
  
  soil.yrsnowpack += soil.getSNOWPACK();
  soil.yrrgrndh2o += y[I_RGRW];
  soil.yrsgrndh2o += y[I_SGRW];

  // Update sum of annual phenology in natural ecosystems

  veg.yrunleaf += y[I_UNRMLF];
  veg.yrleaf += y[I_LEAF];
  veg.yrlai += y[I_LAI];
  veg.yrfpc += y[I_FPC];

  // Update sum of annual carbon fluxes in ecosystems

  veg.yringpp += y[I_INGPP];
  veg.yrgpp   += y[I_GPP];
  veg.yrinnpp += y[I_INNPP];
  veg.yrnpp   += y[I_NPP];
  veg.yrgpr   += y[I_GPR];
  veg.yrrmaint += y[I_RVMNT];
  veg.yrrgrowth += y[I_RVGRW];
  veg.yrabvgrndResp += y[I_ABVGPR];
  veg.yrrootResp += y[I_ROOTGPR];

//  veg.yrvoc.isoprene += y[I_ISOPREN];
//  veg.yrvoc.monoterpene += y[I_TERPEN];
//  veg.yrvoc.otherReactive += y[I_ORVOC];
//  veg.yrvoc.other += y[I_OVOC];

  veg.yrltrc  += y[I_LTRC];
  microbe.yrdecomp += y[I_CDCMP];
  microbe.yrrh += y[I_RH];

  microbe.yrDOMprod.carbon += y[I_DOCP];
  soil.yrlchDOM.carbon += y[I_LCHDOC];
  soil.yrerodePOM.carbon += y[I_ERDPOC];

//  soil.yrdissCO2 += y[I_CO2DISS];
//  soil.yrlchCO2 += y[I_LCHCO2];
//  soil.yrformHCO3 += y[I_HCO3P];
//  soil.yrformRHCO3 += y[I_RHCO3P];
//  soil.yrlchHCO3 += y[I_LCHHCO3];

  yrrsoil += rsoil;
  yrnep += nep;
//  veg.yrvoc.total += veg.voc.total;
  yrnce += nce;
//  soil.yrlchALK += soil.leachALK;  
  yrntcb += ntcb;


 // Update sum of annual nitrogen fluxes in ecosystems

  ag.yrfertn += y[I_AGFRTN];

  yrbnfix += y[I_BNFIX];
  veg.yrnfix += y[I_SNFIX];
  microbe.yrnfix += y[I_ANFIX];

  veg.yrinnup.total += y[I_INNUP];
  veg.yrinnup.nh4 += y[I_INNH4UP];
  veg.yrinnup.no3 += y[I_INNO3UP];

  veg.yrnup.total   += y[I_VNUP];
  veg.yrnup.nh4   += y[I_VNH4UP];
  veg.yrnup.no3   += y[I_VNO3UP];

  veg.yrsup    += y[I_VSUP];
  veg.yrlup    += y[I_VLUP];
  veg.yrnmobil += y[I_VNMBL];
  veg.yrnrsorb += y[I_VNRSRB];

  veg.yrltrn  += y[I_LTRN];

  microbe.yrndecomp += y[I_NDCMP];
  microbe.yrDOMprod.nitrogen += y[I_DONP];

  microbe.yrgmin += y[I_GMIN];
  microbe.yrimmnh4 += y[I_NH4IMM];
  microbe.yrimmb += y[I_NIMM];
  microbe.yrnmin  += y[I_NMIN];

  soil.yrabNimmob += y[I_AIMMNH4];

  microbe.yrammnvol += y[I_AMMN];

  microbe.yrnitrif += y[I_NTRF];
  microbe.yrno3prd += y[I_NO3P];
  microbe.yrnoprd += y[I_NOP];
  microbe.yrn2oprd += y[I_N2OP];
  microbe.yrn2prd += y[I_N2P];
  microbe.yrdenitrf += y[I_DNTRF];

  soil.yrnh3flx += y[I_NH3FLX];
  soil.yrnoflx += y[I_NOFLX];
  soil.yrn2oflx += y[I_N2OFLX];
  soil.yrn2flx += y[I_N2FLX];

  soil.yrlchNO3 += y[I_LCHNO3];
  soil.yrlchDOM.nitrogen += y[I_LCHDON];
  soil.yrerodePOM.nitrogen += y[I_ERDPON];

  atms.yrndep.total += atms.getTOTNDEP();
  atms.yrndep.nh4 += atms.getNH4DEP();
  atms.yrndep.no3 += atms.getNO3DEP();

  soil.yrnin   += soil.getNINPUT();
  soil.yrnlost += soil.getNLOST();
  yrntnb += ntnb;


   // Update sum of annual water fluxes in ecosystems

  ag.yrirrig += y[I_AGIRRIG];
  soil.yrineet += y[I_INEET];
  soil.yreet += y[I_EET];
  soil.yrrperc += y[I_RPERC];
  soil.yrsperc += y[I_SPERC];
  soil.yrrrun += y[I_RRUN];
  soil.yrsrun += y[I_SRUN];

  atms.yrrain += atms.getRAIN();
  atms.yrsnowfall += atms.getSNOWFALL();
  atms.yrpet += atms.getPET();
  soil.yrsnowinf += soil.getSNOWINF();
  soil.yrh2oyld += soil.getH2OYLD();


  ag.yrstubC += ag.getSTUBBLEC();
  ag.yrstubN += ag.getSTUBBLEN();

 // Update sum of annual carbon and nitrogen fluxes from
 //   agricultural conversion

  ag.yrconvrtC += ag.getCONVRTFLXC();
  ag.yrvconvrtC += ag.getVCONVRTFLXC();
  ag.yrsconvrtC += ag.getSCONVRTFLXC();
  ag.yrslashC += ag.getSLASHC();
  ag.yrcflux += ag.getCFLUX();
  
  ag.yrconvrtN += ag.getCONVRTFLXN();
  ag.yrvconvrtN += ag.getVCONVRTFLXN();
  ag.yrsconvrtN += ag.getSCONVRTFLXN();
  ag.yrslashN += ag.getSLASHN();

  ag.yrnrent += ag.getNRETENT();
  ag.yrnvrent += ag.getNVRETENT();
  ag.yrnsrent += ag.getNSRETENT();

  
 // Update sum of annual carbon and nitrogen fluxes from
 //   products

  ag.yrformPROD1C += ag.getCROPPRODC();
  ag.yrformPROD1N += ag.getCROPPRODN();

  ag.yrformResidueC += ag.getFORMCROPRESIDUEC();
  ag.yrformResidueN += ag.getFORMCROPRESIDUEN();

  ag.yrdecayPROD1C += ag.getPROD1DECAYC();
  ag.yrdecayPROD1N += ag.getPROD1DECAYN();

  ag.yrfluxResidueC += ag.getCROPRESIDUEFLXC();
  ag.yrfluxResidueN += ag.getCROPRESIDUEFLXN();

  ag.yrformPROD10C += ag.getFORMPROD10C();
  ag.yrformPROD10N += ag.getFORMPROD10N();

  ag.yrdecayPROD10C += ag.getPROD10DECAYC();
  ag.yrdecayPROD10N += ag.getPROD10DECAYN();

  ag.yrformPROD100C += ag.getFORMPROD100C();
  ag.yrformPROD100N += ag.getFORMPROD100N();

  ag.yrdecayPROD100C += ag.getPROD100DECAYC();
  ag.yrdecayPROD100N += ag.getPROD100DECAYN();

  ag.yrformTOTPRODC += ag.getFORMTOTPRODC();
  ag.yrformTOTPRODN += ag.getFORMTOTPRODN();

  ag.yrdecayTOTPRODC += ag.getTOTPRODDECAYC();
  ag.yrdecayTOTPRODN += ag.getTOTPRODDECAYN();

  // Update sum of soil thermal dynamic variables

  soil.yrtsoil += soil.getTSOIL();
  soil.stm.yrdst0 += soil.stm.getDST0();
  soil.stm.yrdst5 += soil.stm.getDST5();
  soil.yrdst10 += soil.getDST10();
  soil.stm.yrdst20 += soil.stm.getDST20();
  soil.stm.yrdst50 += soil.stm.getDST50();
  soil.stm.yrdst100 += soil.stm.getDST100();
  soil.stm.yrdst200 += soil.stm.getDST200();
  soil.stm.yrfrontd += soil.stm.getFRONTD();
  soil.stm.yrfrontd2 += soil.stm.getFRONTD2();
  soil.stm.yrthawbegin1 += soil.stm.getTHAWBEGIN1();
  soil.stm.yrthawbegin2 += soil.stm.getTHAWBEGIN2();
  soil.stm.yrthawend1 += soil.stm.getTHAWEND1();
  soil.stm.yrthawend2 += soil.stm.getTHAWEND2();
  veg.yrthawpct += veg.getTHAWPCT();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM60::urbanDynamics( const int& pdm, double pstate[] )
{
  int fertflag = 0;
  int irrgflag = 0;
  int perennial = 1;
  int tillflag = 0;

  double newpctp;
  double newsh2o;
  double newvsm;

  double propReact;
  
  soil.updateHydrology( elev,
                        atms.getTAIR(),
                        atms.getPREVTAIR(),
                        atms.getPREV2TAIR(),
                        atms.getRAIN(),
                        atms.getPET(),
                        pstate[I_SM],
                        pstate[I_RGRW],
                        pstate[I_SGRW],
                        irrgflag,
                        ag.irrigate,
                        pdm );
                       
                                      
  // Assume no agricultural N fertilization (ag.fertn) and no
  //  nutrients resulting from agricultural conversion

  ag.fertn = ZERO;
  microbe.setAMMNVOL( ZERO );

  // Add fireNDEP to ecosystem if fire occurred within the past 
  //   FRI years.  Otherwise do not add any nitrogen to ecosyste
 
//  if( firemnthcnt > 0 ) 
//  { 
//    soil.setNINPUT( ag.getNRETENT() + ag.getFIRENDEP() ); 
//  }
//  else { soil.setNINPUT( ZERO ); } 


  // Determine biological N fixation and partition it between
  //   symbiotic N fixation (veg.nfix) and nonsymbiotic N
  //   fixation (microbe.nfix)

  if( soil.getEET() > 0.1 )
  {
    bnfix = veg.getNFIXPARA( ag.cmnt ) * soil.getEET()
            + veg.getNFIXPARB( ag.cmnt );
  }
  else { bnfix = ZERO; }

  microbe.setNFIX( (bnfix * microbe.getNFIXPAR( veg.cmnt )) );
  
  if( microbe.getNFIX() < ZERO ) { microbe.setNFIX( ZERO ); }

  veg.setNFIX( (bnfix - microbe.getNFIX()) );

  if( veg.getNFIX() < ZERO ) { veg.setNFIX( ZERO ); }
  

  // Reduce EET if vegetation is not mature
  
  soil.setEVAPORATION( soil.getEET() 
                       * (1.0 - veg.getPROPTRANS( ag.cmnt)) );
					   
  veg.updateFoliage( ag.cmnt, pstate[I_FOLC], pstate[I_VEGC], soil.getEET() );

  soil.setEET( veg.getTRANSPIRATION() + soil.getEVAPORATION() );


  // Assume wetlands are wetter by the wfpsoff for determining
  //   moisture effects on vegetation and microbes

  newpctp = soil.getPCTP() + soil.getWFPSOFF();
  
  newsh2o = (newpctp * soil.getTOTPOR() * soil.getROOTZ()) 
		    / (100.0 * soil.getACTLAYER()); 

  newvsm = newsh2o / (soil.getROOTZ() * 1000.0);


  soil.setKH2O( newvsm, moistlim );

  // Get proportion of unfrozen organic matter in rooting zone
  
  propReact = soil.getThawedReactiveSOMProp( veg.cmnt );

  microbe.updateDynamics( veg.cmnt,
                          soil.getPCTFLDCAP(),
                          soil.getACTLAYER(),
                          (pstate[I_SOLC] * propReact),
                          pstate[I_AGR],
                          pstate[I_AGL],
                          pstate[I_BGR],
                          pstate[I_BGL],
                          pstate[I_CWD],
                          (pstate[I_SOLN] * propReact),
                          newsh2o,
                          newvsm,
                          (pstate[I_NH4] * soil.getACTLAYER() / soil.getROOTZ()),
                          (atms.getNH4DEP()+ ag.getNRETENT() + ag.getFIRENDEP()),
                          moistlim,
                          tillflag,
                          ag.getTILLFACTOR( veg.cmnt ),
                          soil.getKH2O() );
  /*
  if( disturbflag > 1 && pdm == (disturbmonth-1) )
  {
    // Set monthly vegetation fluxes to zero
    	
    veg.resetMonthlyFluxes();  	
  } 
  else*/
  {
    veg.updateDynamics( latitude,
    					longitude,
    					pdm,
    					ag.cmnt,
                        atms.getCO2(),
                        atms.getAOT40(),
						atms.getTAIR(),
                        atms.getNDEP(),
                        (ag.getNRETENT() + ag.getFIRENDEP()),
                        atms.getPAR(),
                        atms.getPET(),
                        atms.getPRVPETMX(),
                        soil.getEET(),
                        soil.getPRVEETMX(),
                        pstate[I_VEGC],
                        pstate[I_FOLC],
                        pstate[I_STEMC],
                        pstate[I_CROOTC],
                        pstate[I_FROOTC],
                        pstate[I_DEADWOODC],
                        pstate[I_STRN],
                        pstate[I_STON],
                        pstate[I_DEADWOODN],
                        newsh2o,
                        (pstate[I_NH4] * soil.getACTLAYER() / soil.getROOTZ()),
                        (pstate[I_NO3] * soil.getACTLAYER() / soil.getROOTZ()),
                        moistlim,
                        nfeed,
                        o3flag,
                        ag.state,
                        perennial,
                        fertflag,
                        soil.getKH2O(),
                        microbe.getNETNMIN(),
                        microbe.getAMMNVOL(),
                        microbe.getNITRIF(),
                        microbe.getNO3PROD(),
                        ag.fertn,
                        pstate[I_AGE]);
  }                    

  // Determine carbon and nitrogen leaching losses
  
  soil.updateLeachingLosses( veg.cmnt,
                             (pstate[I_DOC] * propReact), 
                             (pstate[I_DON] * propReact), 
                             (pstate[I_NO3] * soil.getACTLAYER() / soil.getROOTZ()), 
                             (pstate[I_SM] * soil.getACTLAYER() / soil.getROOTZ()) );
  

  if ( soil.getLEACHDOC() > (pstate[I_DOC]  
       + microbe.getDOCPROD()) )
  {
    soil.setLEACHDOC( (pstate[I_DOC] + microbe.getDOCPROD()) );
  }

  if ( soil.getLEACHDON() > (pstate[I_DON] 
         + microbe.getDONPROD()) )
  {
    soil.setLEACHDON( (pstate[I_DON] + microbe.getDONPROD()) );
  }

  
  // Determine loss of POC through erosion
  
  soil.setERODEPOC( ZERO );


  // Determine loss of PON through erosion
  
  soil.setERODEPON( ZERO );

               
  // Determine trace gas production based on microbe.no3prod,
  //   ag.fertn and water-filled pore space

  microbe.setTraceGasProduction( veg.cmnt,
                                 newpctp, 
                                 ag.fertn );
 
 
  
  if( soil.getLCHNO3PAR() > ZERO )
  {
    // Limit ecosystem N losses to total NO3 available for 
    //   leaching after denitrification has been substracted
    
    if ( soil.getLEACHNO3() > (pstate[I_NO3] + atms.getNO3DEP()
         + microbe.getNO3PROD() - microbe.getDENITRIF() 
         - veg.getNO3UPTAKE()) )
    {
      soil.setLEACHNO3( (pstate[I_NO3]
                         + atms.getNO3DEP()
                         + microbe.getNO3PROD()
                         - microbe.getDENITRIF()
                         - veg.getNO3UPTAKE()) );
    }
  }
  else     
  {
    // If no leaching losses, limit ecosystem N losses to 
    //   total NO3 available for denitrification
    
    if ( microbe.getDENITRIF() > (pstate[I_NO3] + atms.getNO3DEP()
         + microbe.getNO3PROD() - veg.getNO3UPTAKE()) )
    {
      microbe.setDENITRIF( (pstate[I_NO3]
                           + atms.getNO3DEP()
                           + microbe.getNO3PROD()
                           - veg.getNO3UPTAKE()) );
                           
      if( microbe.getN2OPROD() > ZERO )
      {
        microbe.setN2PROD( microbe.getDENITRIF() - microbe.getN2OPROD() );
      
        if( microbe.getN2PROD() < ZERO )
        {
          microbe.setN2PROD( ZERO );
          microbe.setN2OPROD( microbe.getDENITRIF() );
        }
      }
      else { microbe.setN2PROD( microbe.getDENITRIF() ); }           
    }    
  }

  // Determine trace gas fluxes based on ammonia volatilization
  //   (microbe.ammnvol), NO production (microbe.noprod),
  //   N2O production (microbe.n2oprod) and N2 production
  //   (microbe.n2prod)

  soil.setTraceGasFluxes( microbe.getAMMNVOL(),
                          microbe.getNOPROD(),
                          microbe.getN2OPROD(),
                          microbe.getN2PROD() );

  if( 0 == avlnflag )
  {
    // Do not allow changes in available nitrogen in soil
    //   (used during calibration process)

    // Equilibrate Ammonium (NH4) pool

    microbe.setAMMNVOL( (atms.getNH4DEP()
                        + microbe.getGROSSNMIN()
                        - microbe.getIMMNH4()
                        - veg.getNH4UPTAKE()
                        - microbe.getNITRIF()) );

    soil.setNH3FLUX( microbe.getAMMNVOL() );
    
    
    // Equilibrate nitrate (NO3) pool

    soil.setLEACHNO3( (atms.getNO3DEP()
                       + microbe.getNO3PROD()
                       - microbe.getDENITRIF()
                       - veg.getNO3UPTAKE()) );

    // Equilibrate DON pool

    soil.setLEACHDON( microbe.getDONPROD() );
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM60::writesitecd( ofstream& fout, const int& dcmnt )
{
  fout << "<?xml version = \"1.0\"?>" << endl << endl;

  fout << "<siteECD version = \"" << version << "\"" << endl;
  fout << "  site = \"" << sitename << "\"" << endl;
  fout << "  longitude = \"" << sitecol << "\"" << endl;
  fout << "  latitude = \"" << siterow << "\"" << endl;
  fout << "  developedBy = \"" << developer << "\"" << endl;
  fout << "  updated = \"" << updated << "\">" << endl;
  fout << endl;
  fout << "  <community type = \"" << veg.cmnt << "\"" << endl;
  fout << "    description = \"" << description << "\">" << endl;
  fout << endl;
  fout << "    <vegca>" << vegca[dcmnt] << "</vegca>" << endl;
  fout << "    <vegcb>" << vegcb[dcmnt] << "</vegcb>" << endl;
  fout << "    <strna>" << strna[dcmnt] << "</strna>" << endl;
  fout << "    <strnb>" << strnb[dcmnt] << "</strnb>" << endl;
  fout << "    <stona>" << stona[dcmnt] << "</stona>" << endl;
  fout << "    <stonb>" << stonb[dcmnt] << "</stonb>" << endl;
  fout << "    <solca>" << solca[dcmnt] << "</solca>" << endl;
  fout << "    <solcb>" << solcb[dcmnt] << "</solcb>" << endl;
  fout << "    <solna>" << solna[dcmnt] << "</solna>" << endl;
  fout << "    <solnb>" << solnb[dcmnt] << "</solnb>" << endl;
  fout << "    <doccut>" << doccut[dcmnt] << "</doccut>" << endl;
  fout << "    <doc1a>" << doc1a[dcmnt] << "</doc1a>" << endl;
  fout << "    <doc1b>" << doc1b[dcmnt] << "</doc1b>" << endl;
  fout << "    <doc2a>" << doc2a[dcmnt] << "</doc2a>" << endl;
  fout << "    <doc2b>" << doc2b[dcmnt] << "</doc2b>" << endl;
  fout << "    <doncut>" << doncut[dcmnt] << "</doncut>" << endl;
  fout << "    <don1a>" << don1a[dcmnt] << "</don1a>" << endl;
  fout << "    <don1b>" << don1b[dcmnt] << "</don1b>" << endl;
  fout << "    <don2a>" << don2a[dcmnt] << "</don2a>" << endl;
  fout << "    <don2b>" << don2b[dcmnt] << "</don2b>" << endl;
  fout << "    <nh4a>" << nh4a[dcmnt] << "</nh4a>" << endl;
  fout << "    <nh4b>" << nh4b[dcmnt] << "</nh4b>" << endl;
  fout << "    <no3cut>" << no3cut[dcmnt] << "</no3cut>" << endl;
  fout << "    <no31a>" << no31a[dcmnt] << "</no31a>" << endl;
  fout << "    <no31b>" << no31b[dcmnt] << "</no31b>" << endl;
  fout << "    <no32a>" << no32a[dcmnt] << "</no32a>" << endl;
  fout << "    <no32b>" << no32b[dcmnt] << "</no32b>" << endl;
  fout << endl;
  fout << "    <unleaf12>" << veg.getUNLEAF12( dcmnt );
  fout << "</unleaf12>" << endl;
  fout << "    <initleafmx>" << veg.getINITLEAFMX( dcmnt );
  fout << "</initleafmx>" << endl;
  fout << endl;
  fout << "    <vegcmaxcut>" << veg.getCMAXCUT( dcmnt );
  fout << "</vegcmaxcut>" << endl;
  fout << "    <vegcmax1a>" << veg.getCMAX1A( dcmnt );
  fout << "</vegcmax1a>" << endl;
  fout << "    <vegcmax1b>" << veg.getCMAX1B( dcmnt );
  fout << "</vegcmax1b>" << endl;
  fout << "    <vegcmax2a>" << veg.getCMAX2A( dcmnt );
  fout << "</vegcmax2a>" << endl;
  fout << "    <vegcmax2b>" << veg.getCMAX2B( dcmnt );
  fout << "</vegcmax2b>" << endl;
  fout << "    <vegcfall>" << veg.getCFALL( dcmnt );
  fout << "</vegcfall>" << endl;
  //fout << "    <vegkra>" << veg.getKRA( dcmnt );
  //fout << "</vegkra>" << endl;
  //fout << "    <vegkrb>" << veg.getKRB( dcmnt );
  //fout << "</vegkrb>" << endl;
  fout << "    <vegrmmax>" << veg.getRMMAX( dcmnt );
  fout << "</vegrmmax>" << endl;
  fout << "    <vegrroot>" << veg.getRROOT( dcmnt );
  fout << "</vegrroot>" << endl;
  fout << "    <microbekdcut>" << microbe.getKDCUT( dcmnt );
  fout << "</microbekdcut>" << endl;
  fout << "    <microbekd1a>" << microbe.getKD1A( dcmnt );
  fout << "</microbekd1a>" << endl;
  fout << "    <microbekd1b>" << microbe.getKD1B( dcmnt );
  fout << "</microbekd1b>" << endl;
  fout << "    <microbekd2a>" << microbe.getKD2A( dcmnt );
  fout << "</microbek2da>" << endl;
  fout << "    <microbekd2b>" << microbe.getKD2B( dcmnt );
  fout << "</microbekd2b>" << endl;
  fout << "    <microbepropftos>" << microbe.getPROPFTOS( dcmnt );
  fout << "</microbepropftos>" << endl;
  fout << "    <soilnonOMpar>" << soil.getNSOLPAR( dcmnt );
  fout << "</soilnonOMpar>" << endl;
  fout << "    <microbeDOCpar>" << microbe.getDOCPAR( dcmnt );
  fout << "</microbeDOCpar>" << endl;
  fout << "    <soillchDOMpar>" << soil.getLCHDOMPAR( dcmnt );
  fout << "</soillchDOMpar>" << endl;  
  fout << endl;
  fout << "    <vegnfixpara>" << veg.getNFIXPARA( dcmnt );
  fout << "</vegnfixpara>" << endl;
  fout << "    <vegnfixparb>" << veg.getNFIXPARB( dcmnt );
  fout << "</vegnfixparb>" << endl;
  fout << "    <vegnupnh4cut>" << veg.getNUPNH4CUT( dcmnt );
  fout << "</vegnupnh4cut>" << endl;
  fout << "    <vegnupnh41a>" << veg.getNUPNH41A( dcmnt );
  fout << "</vegnupnh41a>" << endl;
  fout << "    <vegnupnh41b>" << veg.getNUPNH41B( dcmnt );
  fout << "</vegnupnh41b>" << endl;
  fout << "    <vegnupnh42a>" << veg.getNUPNH42A( dcmnt );
  fout << "</vegnupnh42a>" << endl;
  fout << "    <vegnupnh42b>" << veg.getNUPNH42B( dcmnt );
  fout << "</vegnupnh42b>" << endl;
  fout << "    <vegnupno3cut>" << veg.getNUPNO3CUT( dcmnt );
  fout << "</vegnupno3cut>" << endl;
  fout << "    <vegnupno31a>" << veg.getNUPNO31A( dcmnt );
  fout << "</vegnupno31a>" << endl;
  fout << "    <vegnupno31b>" << veg.getNUPNO31B( dcmnt );
  fout << "</vegnupno31b>" << endl;
  fout << "    <vegnupno32a>" << veg.getNUPNO32A( dcmnt );
  fout << "</vegnupno32a>" << endl;
  fout << "    <vegnupno32b>" << veg.getNUPNO32B( dcmnt );
  fout << "</vegnupno32b>" << endl;
  fout << "    <vegnfall>" << veg.getNFALL( dcmnt );
  fout << "</vegnfall>" << endl;
  fout << "    <veglcclnc>" << veg.getLCCLNC( dcmnt );
  fout << "</veglcclnc>" << endl;
  fout << "    <microbenfixpar>" << microbe.getNFIXPAR( dcmnt );
  fout << "</microbenfixpar>" << endl;
  fout << "    <microbenh4immcut>" << microbe.getNH4IMMCUT( dcmnt );
  fout << "</microbenh4immcut>" << endl;
  fout << "    <microbenh4imm1a>" << microbe.getNH4IMM1A( dcmnt );
  fout << "</microbenh4imm1a>" << endl;
  fout << "    <microbenh4imm1b>" << microbe.getNH4IMM1B( dcmnt );
  fout << "</microbenh4imm1b>" << endl;
  fout << "    <microbenh4imm2a>" << microbe.getNH4IMM2A( dcmnt );
  fout << "</microbenh4imm2a>" << endl;
  fout << "    <microbenh4imm2b>" << microbe.getNH4IMM2B( dcmnt );
  fout << "</microbenh4imm2b>" << endl;
  fout << "    <microbeammnpar>" << microbe.getAMMNPAR( dcmnt );
  fout << "</microbeammnpar>" << endl;
  fout << "    <microbentrfparcut>" << microbe.getNTRFPARCUT( dcmnt );
  fout << "</microbentrfparcut>" << endl;
  fout << "    <microbentrfpar1a>" << microbe.getNTRFPAR1A( dcmnt );
  fout << "</microbentrfpar1a>" << endl;
  fout << "    <microbentrfpar1b>" << microbe.getNTRFPAR1B( dcmnt );
  fout << "</microbentrfpar1b>" << endl;
  fout << "    <microbentrfpar2a>" << microbe.getNTRFPAR2A( dcmnt );
  fout << "</microbentrfpar2a>" << endl;
  fout << "    <microbentrfpar2b>" << microbe.getNTRFPAR2B( dcmnt );
  fout << "</microbentrfpar2b>" << endl;
  fout << "    <microbeinitntrf>" << microbe.getINITNTRF( dcmnt );
  fout << "</microbeinitntrf>" << endl;
  fout << "    <microbeallntrf>" << microbe.getALLNTRF( dcmnt );
  fout << "</microbeallntrf>" << endl;
  fout << "    <soillchNO3parcut>" << soil.getLCHNO3PARCUT( dcmnt );
  fout << "</soillchNO3parcut>" << endl;
  fout << "    <soillchNO3par1a>" << soil.getLCHNO3PAR1A( dcmnt );
  fout << "</soillchNO3par1a>" << endl;
  fout << "    <soillchNO3par1b>" << soil.getLCHNO3PAR1B( dcmnt );
  fout << "</soillchNO3par1b>" << endl;
  fout << "    <soillchNO3par2a>" << soil.getLCHNO3PAR2A( dcmnt );
  fout << "</soillchNO3par2a>" << endl;
  fout << "    <soillchNO3par2b>" << soil.getLCHNO3PAR2B( dcmnt );
  fout << "</soillchNO3par2b>" << endl;
  fout << "    <microbeDONparcut>" << microbe.getDONPARCUT( dcmnt );
  fout << "</microbeDONparcut>" << endl;
  fout << "    <microbeDONpar1a>" << microbe.getDONPAR1A( dcmnt );
  fout << "</microbeDONpar1a>" << endl;
  fout << "    <microbeDONpar1b>" << microbe.getDONPAR1B( dcmnt );
  fout << "</microbeDONpar1b>" << endl;
  fout << "    <microbeDONpar2a>" << microbe.getDONPAR2A( dcmnt );
  fout << "</microbeDONpar2a>" << endl;
  fout << "    <microbeDONpar2b>" << microbe.getDONPAR2B( dcmnt );
  fout << "</microbeDONpar2b>" << endl;
  fout << "    <veginitcneven>" << veg.getINITCNEVEN( dcmnt );
  fout << "</veginitcneven>" << endl;
  fout << "    <vegcnmin>" << veg.getCNMIN( dcmnt );
  fout << "</vegcnmin>" << endl;
  fout << "    <vegc2na>" << veg.getC2NA( dcmnt );
  fout << "</vegc2na>" << endl;
  fout << "    <vegc2nb>" << veg.getC2NB( dcmnt );
  fout << "</vegc2nb>" << endl;
  fout << "    <vegc2nmin>" << veg.getC2NMIN( dcmnt );
  fout << "</vegc2nmin>" << endl;
  fout << "    <microbecnsoil>" << microbe.getCNSOIL( dcmnt );
  fout << "</microbecnsoil>" << endl;
  fout << endl;
  fout << "    <o3para>" << veg.getO3PARA( dcmnt );
  fout << "</o3para>" << endl;
  fout << "    <o3parb>" << veg.getO3PARB( dcmnt );
  fout << "</o3parb>" << endl;
  fout << "    <o3parc>" << veg.getO3PARC( dcmnt );
  fout << "</o3parc>" << endl;
  fout << "  </community>" << endl;
  fout << "</siteECD>" << endl;


};



double TTEM60::omitzero( const double& data1)
{
	double tempx;
	if (data1 ==0.0) tempx = mpe;
	else tempx = data1;
	return tempx;

};







