/* *************************************************************
****************************************************************
QSOILTEMP602.H - Soil temperature model developed by Qianlai
                  Zhuang for the Terrestrial Ecosystem Model.
                  Model based on Goodrich one-dimensional
                  soil thermal model (Goodrich 1976, 1978a, 
                  1978b; Romanovsky et al. 1997)
****************************************************************

Modifications:

20040420 - DWK created by modifying qsoiltemp50b1.h
20040420 - DWK changed const int MXTSDATA = 1500 to
           MAXTSDATA = 3
20040420 - DWK changed const int MXTSFITTED = 3000 to
           MXTSFITTED = 60
20040420 - DWK changed include from qsoiltemp50b1.cpp to
           qsoiltemp50b2.cpp at bottom of file
20040421 - DWK added const int MAXSPDATA = 25
20040421 - DWK changed double DEPTH[MAXCMNT][25] to 
           double DEPTH[MAXCMNT][MAXSPDATA]
20040421 - DWK changed double TEMP[MAXCMNT][25] to 
           double TEMP[MAXCMNT][MAXSPDATA]
20040421 - DWK changed double spare[6] to 
           double spare[MAXSLAYERS]
20040421 - DWK deleted const double& xtair and 
           const int& xtairflg from function call to 
           initMonthlySoilConditions() and 
           setMonthlySoilConditions()
20040421 - DWK changed x{MXTSDATA} to x[], y[MXTSDATA] to y[],
           xx[MXTSFITTED] to xx[], and yy[MXTSFITTED] to yy[]
           in function call to interpolate()
20040421 - DWK renamed updateSoilTemperature() to be 
           interpolateAirHeat() and eliminated integer *index
           from function call
20040421 - DWK renamed makeSnowPack() to be setSnowPack()
20040421 - DWK renamed setHeatFlux() to be setHeatSources()
20040421 - DWK deleted private double caltim[MXTSFITTED],
           double final, and double per
20040421 - DWK added private double total
20040421 - DWK changed void setSnowDensity() to 
           double setSnowDensity()
20040423 - DWK renamed setFreezeThawIndex() to be 
           setAirFreezeThawIndex()
20040423 - DWK renamed updateTempDownsweep() to be
           crankNicholsonBackward()
20040423 - DWK renamed updateTempUpsweep() to be
           crankNicholsonForward()
20040423 - DWK added private function setSnowCharSturm()
20040424 - DWK added private functions setFullSnowNode() and 
           interpolateSnow() to replace makeSnowPack()
20040424 - DWK changed STMmonth.previous to be STMmonth.beginning
20040424 - DWK changed STMmonth.current to be STMmonth.middle
20040424 - DWK changed STMmonth.next to be STMmonth.end
20040424 - DWK added const double& thresholdFreezeT and 
           const double& latentHeat to resetXFA()
20040424 - DWK added const int& spaceTimeIndex to interpolate()
20040424 - DWK added private enum interpIndex to provide indexes
           for interpolate() to act over space or time
20040424 - DWK renamed setSoilType13() as 
           setSoilCharKerstenGravelSand()                                                                    
20040424 - DWK renamed setSoilType14() as 
           setSoilCharKerstenClaySilt()                                                                    
20040424 - DWK renamed setSoilType21() as setSnowCharWoodside()                                                                    
20040424 - DWK renamed setSoilType22() as setSnowCharDevaux()                                                                    
20040424 - DWK renamed setSoilType23() as 
           setSoilCharDeVriesGravelSand()                                                                    
20040424 - DWK renamed setSoilType24() as 
           setSoilCharDeVriesClay()                                                                    
20040424 - DWK renamed setSoilType25() as 
           setPeatCharDeVries()                                                                    
20040424 - DWK renamed setSoilType26() as 
           setSoilCharDeVriesSand()                                                                    
20040424 - DWK renamed setSoilType27() as 
           setSoilCharKerstenSand()                                                                    
20040424 - DWK renamed setTempDownsweep() with 
           gaussianEliminForward()
20040424 - DWK renamed setTempUpsweep() with 
           gaussianEliminBackward()
20040424 - DWK deleted double sigma
20040425 - DWK added private functions resetTimeStep(),
           createSinglePhasePlane() and createMorePhasePlanes()
20040425 - DWK renamed setSinglePhasePlane() as 
           updateSinglePhasePlaneTemps()
20040425 - DWK renamed setMultiplePhasePlanes as
           updateMultiplePhasePlaneTemps()
20040425 - DWK added global const double SECSPERDAY
20040425 - DWK changed private double tair to double tstair
20040426 - DWK added private functions equilibriumCheck() and
           updateActiveLayer()
20040426 - DWK deleted public double diffsoilt[10]
20040427 - DWK added private function setSnowPackGrid()
20040428 - DWK changed double *acum to const double& acum in 
           function call to updateSnowPackGrid()
20040428 - DWK added const integer& cmnt to the function calls
           of setFullSnowNode() and setSnowPackGrid()
20040428 - DWK deleted const int& outmon from function call of 
           setSoilProfileGrid(() 
20040501 - DWK changed global const int MXTSFITTED = 60 to
           MXTSFITTED = 120
20040501 - DWK deleted private function updateSurfaceHeatFlux()
20040502 - DWK changed private double capx[MAXNODES] to 
           double capx[MAXNODES+2]
20040502 - DWK changed private double conx[MAXNODES] to 
           double conx[MAXNODES+2]
20040502 - DWK changed private double ddry[MAXNODES] to 
           double ddry[MAXNODES+2]                                                        
20040502 - DWK changed private double dx[MAXNODES] to 
           double dx[MAXNODES+2]
20040502 - DWK changed private double dx9[MAXNODES] to 
           double dx9[MAXNODES+2]
20040502 - DWK changed private double e[MAXNODES] to 
           double e[MAXNODES+2]
20040502 - DWK changed private double ht[MAXNODES] to 
           double ht[MAXNODES+2]
20040502 - DWK changed private double htold[MAXNODES] to 
           double htold[MAXNODES+2]
20040502 - DWK changed private integer mater[MAXNODES] to 
           integer mater[MAXNODES+2]
20040502 - DWK changed private double s[MAXNODES] to 
           double s[MAXNODES+2]
20040502 - DWK changed private double t[MAXNODES] to 
           double t[MAXNODES+2]
20040502 - DWK changed private double t9[MAXNODES] to 
           double t9[MAXNODES+2]
20040502 - DWK changed private double tanee[MAXNODES] to 
           double tanee[MAXNODES+2]
20040502 - DWK changed private double thalt[MAXNODES] to 
           double thalt[MAXNODES+2]
20040502 - DWK changed private double told[MAXNODES] to 
           double told[MAXNODES+2]
20040502 - DWK changed private double water[MAXNODES] to 
           double water[MAXNODES+2]
20040502 - DWK changed private double water9[MAXNODES] to 
           double water9[MAXNODES+2]
20040502 - DWK changed private double x[MAXNODES] to 
           double x[MAXNODES+2]
20040502 - DWK changed private double x9[MAXNODES] to 
           double x9[MAXNODES+2]
20040502 - DWK changed private double xfa[MAXNODES] to 
           double xfa[MAXNODES+2]
20040502 - DWK changed private double xfa9[MAXNODES] to 
           double xfa9[MAXNODES+2]
20040502 - DWK changed private double xfb[MAXNODES] to 
           double xfb[MAXNODES+2]
20040502 - DWK changed private double xfb9[MAXNODES] to 
           double xfb9[MAXNODES+2]
20040502 - DWK changed private double xhet[MAXNODES] to 
           double xhet[MAXNODES+2]
20040502 - DWK changed private double ONOFF[MAXCMNT][MAXNODES] 
           to double ONOFF[MAXCMNT][MAXNODES+2]
20040502 - DWK deleted private integer LISO[MAXCMNT], 
           double TISO[MAXCMNT], integer LMAX[MAXCMNT], 
           double VDEPTH[MAXCMNT], double DEPTEM1[MAXCMNT], 
           double DEPTEM2[MAXCMNT], double DEPTEM3[MAXCMNT],
           double DEPTEM4[MAXCMNT], double DEPTEM5[MAXCMNT],
           integer NDEPF[MAXCMNT], double VEP[MAXCMNT], 
           double DEPFLX1[MAXCMNT], double DEPFLX2[MAXCMNT],
           double DEPFLX3[MAXCMNT], double DEPFLX4[MAXCMNT],
           double DEPFLX5[MAXCMNT], double HLAT[MAXCMNT],
           double initFINAL[MAXCMNT], double initPER[MAXCMNT],
           double initTHETA[MAXCMNT], double TOP{MAXCMNT],
           integer vegIG[MAXCMNT], double EPSMIN[MAXCMNT],
           double VSPACE[MAXCMNT], double VDEN[MAXCMNT]
           double VDEP1[MAXCMNT], double DEPHET[MAXCMNT],
           integer SNOFAL[MAXCMNT], double EPSSNO[MAXCMNT],
           double CONVRT[MAXCMNT], double initETAO[MAXCMNT],
           double DENFAC[MAXCMNT], initFCMELT[MAXCMNT],
           double DENMAX[MAXCMNT], integer FIRST[MAXCMNT], 
           and integer KINTI[MAXCMNT]
20040504 - DWK deleted private integer INDEX{MAXCMNT] and 
           double VDEPP[MAXCMNT]
20040504 - DWK added public function initSoilThermalRegime()
20040504 - DWK changed global const integer MAXMAT = 30 to
           MAXMAT = 14
20040504 - DWK deleted integer *ktop from function call to
           private function setUpperBoundaryConditions()
20040504 - DWK deleted private double emiss, integer ktop, 
           integer nan, integer nst, double zz and 
           double spare[MAXSLAYERS]
20040515 - DWK changed include from qsoiltemp50b2.cpp to
           qsoiltemp50b3.cpp at bottom of file
20040515 - DWK deleted const int& cmnt from setFullSnowNode()
20040515 - DWK deleted const integer& cmnt from 
           interpolateAirHeat()
20040515 - DWK deleted const integer& cmnt from 
           setUpperBoundaryConditions()
20040515 - DWK deleted const int& cmnt from setSnowPackGrid()
20040515 - DWK deleted const int& cmnt from setHeatSources()
20040515 - DWK deleted private functions equilibriumCheck(),
           setPeatCharDeVries(), setSnowCharDevaux(),
           setSnowCharWoodside(), setSoilCharDeVriesClay(),
           setSoilCharDeVriesGravelSand(), 
           setSoilCharDeVriesSand(), 
           setSoilCharKerstenClaySilt(), 
           setSoilCharKerstenGravelSand(),
           setSoilCharKerstenSand()
20040515 - DWK changed global const MAXMAT = 14 to MAXMAT = 5
20040515 - DWK deleted struct STMsoiltype
20040515 - DWK changed global const MAXSNODES = 10 to 
           MAXSNODES = 9
20040515 - DWK replaced array size MAXNODES+2 with MAXNODES+1
20040517 - DWK changed access from private to public for 
           int is9, double smass9, 
           double weight9[MAXSNODES], double x9[MAXNODES], 
           double dx9[MAXNODES], double xfa9[MAXNODES],
           double xfb9[MAXNODES], double water9[MAXNODES],
20040517 - DWK changed access from private to public for 
           function setSnowPackGrid(), and renamed function as 
           setInitSnowPackGrid()
20040518 - DWK deleted const int& cmnt and ofstream& rflog1 
           from arguments to initSoilThermalRegime() and 
           changed access from public to private
20040518 - DWK changed access of setInitSnowPackGrid() from 
           public to private  
20040707 - DWK changed include from tprocessXML431a.h to 
           tprocessXML60.h
20040707 - DWK changed class Soilthermal to Soilthermal60
20040707 - DWK changed inheritance of ProcessXML to inheritance
           of ProcessXML60
20040707 - DWK changed include from qsoiltemp50b3.cpp to
           qsoiltemp60.cpp at bottom of file
20040707 - DWK changed void getsnowecd( char ecd[MAXFNAME] ) to
           void getsnowecd( string ecd )
20040707 - DWK changed void getsoillecd( char ecd[MAXFNAME] ) to
           void getsoillecd( string ecd )
20040707 - DWK changed void getsoiltecd( char ecd[MAXFNAME] ) to
           void getsoiltecd( string ecd )
20051117 - DWK changed tprocessXML60.h to tprocessXML602.h
20051117 - DWK deleted values.h standard include
20051117 - DWK deleted qsoiltem60.cpp from bottom of file
20051125 - DWK changed public double activeLayer[CYCLE] to 
           double activeLayer
20051125 - DWK changed public double dst0[CYCLE] to double dst0 
20051125 - DWK changed public double dst5[CYCLE] to double dst5                     
20051125 - DWK changed public double dst10[CYCLE] to 
           double dst10
20051125 - DWK changed public double dst20[CYCLE] to 
           double dst20           
20051125 - DWK changed public double dst50[CYCLE] to 
           double dst50
20051125 - DWK changed public double dst100[CYCLE] to 
           double dst100
20051125 - DWK changed public double dst200[CYCLE] to 
           double dst200
20051125 - DWK changed public double frontd[CYCLE] to 
           double frontd 
20051125 - DWK changed public double thawbegin[CYCLE] to 
           double thawbegin                                                      
20051125 - DWK changed public double thawend[CYCLE] to 
           double thawend
20051125 - DWK changed public double tsoil[CYCLE] to 
           double tsoil
20051125 - DWK deleted const int& outmon from function call to
           setMonthlySoilTemp() and updateActiveLayer()
20051213 - DWK added public function resetYrFluxes()
20051218 - DWK changed private double thawbegin to 
           double thawbegin1
20051218 - DWK changed private double thawend to 
           double thawend1
20051218 - DWK added private double thawbegin2 and 
           double thawend2
20051218 - DWK added private double frontd2
20051219 - DWK deleted const int& pcmnt from function call to
           setMonthlySoilTemp()
20060105 - DWK added public function InterpolatST()
20060109 - DWK changed private function initSoilThermalRegime() 
           to be public function initSoilThermalRegime()
20060109 - DWK added const double& tstart to function call of 
           initSoilThermalRegime()
                                                        
****************************************************************

References:

Goodrich LE (1976) A numerical model for assessing the influence
  of snow cover on the ground thermal regime.  PhD thesis, 
  410 pp., McGill University, Montreal, Quebec.
  
Goodrich LE (1978a) Some results of a numerical study of ground 
  thermal regimes. pp. 29-34. IN: Proceedings of the 3rd 
  International Conference on Permafrost, Vol. 1, Nat. Res.
  Counc. of Canada, Ottawa.
  
Goodrich LE (1978b) Efficient numerical technique for one-
  dimensional thermal problems with phase change. International
  Journal of Heat and Mass Transfer 21, 615-621.
  
Liston, G. E. and R. A. Pielke (2000) A climate version of the 
  regional atmospheric modeling system.  Theoretical and 
  Applied Climatology 66, 29-47.
    
Romanovsky VE, Osterkamp TE, Duxbury NS (1997) An evaluation of
  three numerical models used in simulation of the active layer
  and permafrost temperature regimes.  Cold Reg. Scie. Technol.
  26, 195-203.
  
Sturm, M., J. Holmgren, M. Konig and K. Morris (1997) The
  thermal conductivity of seasonal snow.  Journal of Glaciology
  43(143), 26-41.
  
Sturm, M., J. Holmgren and G. E. Liston (1995) A seasonal snow
  cover classification system for local to global applications.
  Journal of Climate 8(5), 1261-1283.
    
Zhuang Q, Romanovsky VE, McGuire AD (2001) Incorporation of a 
  permafrost model into a large-scale ecosystem model: Evaluation
  of temporal and spatial scaling issues in simulating soil
  thermal dynamics. Journal of Geophysical Research 106(D24),
  33,648-33,670.

Zhuang Q, McGuire AD, O'Neill KP, Harden JW, Romanovsky VE, Yarie J
  (2002) Modeling the soil thermal and carbon dynamics of a fire
  chronosequence in interior Alaska. Journal of Geophysical Research,
  108, D1, 8147, doi:10.1029/2001JD001244.

Zhuang Q, McGuire AD, Melillo JM, Clein JS, Dargaville RJ,
  Kicklighter DW, Myneni RB, Dong J, Romanovsky VE, Harden J,
  Hobbie JE (2003) Carbon cycling in extratropical terrestrial
  ecosystems of the Northern Hemisphere during the 20th century:
  a modeling analysis of the influences of soil thermal dynamics,
  Tellus 55B: 751-776.

******************************************************************
*************************************************************** */


#ifndef QSOILTEMP602_H
#define QSOILTEMP602_H  

//#include<values.h>

#include "tprocessXML602.h"

typedef long int integer;
typedef char *address;
typedef short int shortint;
//typedef float double;
typedef double doubledouble;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

// Maximum number of nodes in the snowpack 
 
#ifndef MAXSNODES_CONST
#define MAXSNODES_CONST
  const int MAXSNODES = 9;
#endif

// Maximum number of soil layers with different thermal 
//   characteristics

const int MAXSLAYERS = 6;

// Maximum number of nodes in snowpack and soil combined

#ifndef MAXNODES_CONST
#define MAXNODES_CONST
  const int MAXNODES = 210;
#endif 

// Maximum number of nodes considered in the output of the soil
//  temperature profile

const int MAXDFST = 30;

// Maximum number of materials for which soil thermal properties 
//   are considered

const int MAXMAT = 5;

// Maximum number of soil depths with "raw" temperature data

const int MAXSPDATA = 25;

// Maximum number of data points per month 

// const int MXTSDATA = 1500;
const int MXTSDATA = 3;

// Maximum number of time steps per month of the interpolated
//   fitted data (assumes 0.5 day time steps)

// const int MXTSFITTED = 3000;
const int MXTSFITTED = 120; // 2x 0.5-day time steps in a month

// Seconds per day

const double SECSPERDAY = 86400.0;

const double PI = 3.1415926536;

const double ADJUST = 0.0001;  

const doubledouble TEN = 10.0; 
const long THOUSANDX100 = 100000;


// Absolute values less than ZEROTOL are considered "zero"

const double ZEROTOL = 0.00001;

struct STMsoilwater 
{
  double layer1;  // moss/litter layer water 
  double layer2;  // organic layer water
  double layer3;

};


struct  STMmonth 
{
  double beginning;
  double middle;
  double end;

};


class Soilthermal60 : public ProcessXML60 
{

  public:

    Soilthermal60();

/* **************************************************************
		 Public Functions
************************************************************** */

     void   getsnowecd( ofstream& rflog1 );
     void   getsnowecd( const string& ecd );
     void   getsoillecd( ofstream& rflog1 );
     void   getsoillecd( const string& ecd );
     void   getsoiltecd( ofstream& rflog1 );
     void   getsoiltecd( const string& ecd );

     void initSoilThermalRegime( const double& tstart );

     double InterpolatST( const double& dst1, 
                          const double& dst2,
                          const double& dst3,
                          const double& dst4, 
                          const double& dst5, 
                          const double& dst6,  
                          double x );

     void resetYr( void );

    
    int    setMonthlySoilConditions( const double& prevtair,
                                     const double& tair,
                                     const double& nexttair,
                                     const double& prevspack,
                                     const double& snowpack,
                                     const double& nextspack,
                                     const int& cmnt,
                                     const int& outmon,
                                     ofstream& rflog1 );

             // setSnowDensity() returns snow density (g cm^-3)
             
    double   setSnowDensity( const int& cmnt, 
                             const double& prec,
                             const int& outmon );

    int setInitSnowPackGrid( const double& prevspack, 
                             const double& spack,  
                             const double& sdensity );
                           
    void   showsnowecd( int& pdcmnt );
    void   showsoillecd( int& pdcmnt );
    void   showsoiltecd( int& pdcmnt );
    void   updateyrCDM( double tair[CYCLE] );


     // "Get" and "Set" functions for private variables and 
     //   parameters
     
     
     // activeLayer ********************************************
      
     inline double getACTLAYER( void ) { return activeLayer; };

     inline void setACTLAYER( const double& pactlayer ) 
     { 
       activeLayer = pactlayer; 
     };


    // cdsnow **************************************************

     inline double getCDSNOW( const int& pcmnt ) 
     { 
       return cdsnow[pcmnt]; 
     };

     inline void setCDSNOW( const double& pcdsnow, 
                            const int& pcmnt ) 
     { 
       cdsnow[pcmnt] = pcdsnow; 
     };


    // dst0 ****************************************************
      
     inline double getDST0( void ) { return dst0; };


    // dst5 ****************************************************
      
     inline double getDST5( void ) { return dst5; };


    // dst10 ***************************************************
      
     inline double getDST10( void ) { return dst10; };


    // dst20 ***************************************************
      
     inline double getDST20( void ) { return dst20; };


    // dst50 ***************************************************
      
     inline double getDST50( void ) { return dst50; };


    // dst100 **************************************************
      
     inline double getDST100( void ) { return dst100; };


    // dst200 **************************************************
      
     inline double getDST200( void ) { return dst200; };

    // dst300 **************************************************

     inline double getDST300( void ) { return dst300; };
    // dx9 *****************************************************

     inline double getDX9( const int& pnode ) 
     { 
       return dx9[pnode]; 
     };

     inline void setDX9( const double& pdx9, 
                        const int& pnode ) 
     { 
       dx9[pnode] = pdx9; 
     };


    // frontd **************************************************
      
     inline double getFRONTD( void ) { return frontd; };


    // frontd2 *************************************************
      
     inline double getFRONTD2( void ) { return frontd2; };


    // gflux ***************************************************

     inline double getGFLUX( const int& pcmnt ) 
     { 
       return gflux[pcmnt]; 
     };

     inline void setGFLUX( const double& pgflux, 
                           const int& pcmnt ) 
     { 
       gflux[pcmnt] = pgflux; 
     };


    // initcondf *************************************************

     inline double getINITCONDF( const int& pcmnt,
                                 const int& pslayr ) 
     { 
       return initCONDF[pcmnt][pslayr]; 
     };

     inline void setINITCONDF( const double& pinitcondf, 
                               const int& pcmnt,
                               const int& pslayr ) 
     { 
       initCONDF[pcmnt][pslayr] = pinitcondf; 
     };


    // initcondt *************************************************

     inline double getINITCONDT( const int& pcmnt,
                                 const int& pslayr ) 
     { 
       return initCONDT[pcmnt][pslayr]; 
     };

     inline void setINITCONDT( const double& pinitcondt, 
                               const int& pcmnt,
                               const int& pslayr ) 
     { 
       initCONDT[pcmnt][pslayr] = pinitcondt; 
     };


    // initdense *************************************************

     inline double getINITDENSE( const int& pcmnt,
                                 const int& pslayr ) 
     { 
       return initDENSE[pcmnt][pslayr]; 
     };

     inline void setINITDENSE( const double& pinitdense, 
                               const int& pcmnt,
                               const int& pslayr ) 
     { 
       initDENSE[pcmnt][pslayr] = pinitdense; 
     };


    // initdtday ***********************************************

     inline double getINITDTDAY( const int& pcmnt ) 
     { 
       return initDTDAY[pcmnt]; 
     };

     inline void setINITDTDAY( const double& pinitdtday, 
                               const int& pcmnt ) 
     { 
       initDTDAY[pcmnt] = pinitdtday; 
     };


    // initdxa *************************************************

     inline double getINITDXA( const int& pcmnt,
                               const int& pslayr ) 
     { 
       return initDXA[pcmnt][pslayr]; 
     };

     inline void setINITDXA( const double& pinitdxa, 
                             const int& pcmnt,
                             const int& pslayr ) 
     { 
       initDXA[pcmnt][pslayr] = pinitdxa; 
     };


    // initdxb *************************************************

     inline double getINITDXB( const int& pcmnt,
                               const int& pslayr ) 
     { 
       return initDXB[pcmnt][pslayr]; 
     };

     inline void setINITDXB( const double& pinitdxb, 
                             const int& pcmnt,
                             const int& pslayr ) 
     { 
       initDXB[pcmnt][pslayr] = pinitdxb; 
     };


    // initsphf *************************************************

     inline double getINITSPHF( const int& pcmnt,
                                 const int& pslayr ) 
     { 
       return initSPHF[pcmnt][pslayr]; 
     };

     inline void setINITSPHF( const double& pinitsphf, 
                               const int& pcmnt,
                               const int& pslayr ) 
     { 
       initSPHF[pcmnt][pslayr] = pinitsphf; 
     };


    // initspht *************************************************

     inline double getINITSPHT( const int& pcmnt,
                                 const int& pslayr ) 
     { 
       return initSPHT[pcmnt][pslayr]; 
     };

     inline void setINITSPHT( const double& pinitspht, 
                               const int& pcmnt,
                               const int& pslayr ) 
     { 
       initSPHT[pcmnt][pslayr] = pinitspht; 
     };


    // initthick *************************************************

     inline double getINITTHICK( const int& pcmnt,
                                 const int& pslayr ) 
     { 
       return initTHICK[pcmnt][pslayr]; 
     };

     inline void setINITTHICK( const double& pinitthick, 
                               const int& pcmnt,
                               const int& pslayr ) 
     { 
       initTHICK[pcmnt][pslayr] = pinitthick; 
     };


    // initwater *************************************************

     inline double getINITWATER( const int& pcmnt,
                                 const int& pslayr ) 
     { 
       return initWATER[pcmnt][pslayr]; 
     };

     inline void setINITWATER( const double& pinitwater, 
                               const int& pcmnt,
                               const int& pslayr ) 
     { 
       initWATER[pcmnt][pslayr] = pinitwater; 
     };


    // is9 *****************************************************
      
     inline int getIS9( void ) { return is9; };

     inline void setIS9( const int& pis9 ) 
     { 
       is9 = pis9; 
     };


    // ism19 *****************************************************
      
     inline int getISM19( void ) { return ism19; };

     inline void setISM19( const int& pism19 ) 
     { 
       ism19 = pism19; 
     };


     // kswitch ************************************************
     
     inline void setKSWITCH( const int& pkswitch ) 
     { 
       kswitch = pkswitch; 
     };


     // nextsnowfall *******************************************
      
     inline double getNEXTSNOWFALL( void ) 
     { 
       return nextsnowfall; 
     };

     inline void setNEXTSNOWFALL( const double& pnextsnwfal ) 
     { 
       nextsnowfall = pnextsnwfal; 
     };


     // nextspack **********************************************
      
     inline double getNEXTSPACK( void ) 
     { 
       return nextspack; 
     };

     inline void setNEXTSPACK( const double& pnextspack ) 
     { 
       nextspack = pnextspack; 
     };


     // nexttair ********************************************
      
     inline double getNEXTTAIR( void ) { return nexttair; }

     inline void setNEXTTAIR( const double& pnexttair ) 
     { 
       nexttair = pnexttair; 
     };


     // prevspack **********************************************
     
     inline double getPREVSPACK( void ) { return prevspack; }

     inline void setPREVSPACK( const double& pprvspack ) 
     { 
       prevspack = pprvspack; 
     };


    // smass9 **************************************************
      
     inline double getSMASS9( void ) { return smass9; };

     inline void setSMASS9( const double& psmass9 ) 
     { 
       smass9 = psmass9; 
     };


    // t9 ******************************************************

     inline double getT9( const int& pnode ) 
     { 
       return t9[pnode]; 
     };

     inline void setT9( const double& pt9, 
                        const int& pnode ) 
     { 
       t9[pnode] = pt9; 
     };


    // thawbegin1 **********************************************
      
     inline double getTHAWBEGIN1( void ) { return thawbegin1; };


    // thawbegin2 **********************************************
      
     inline double getTHAWBEGIN2( void ) { return thawbegin2; };


    // thawend1 ************************************************
      
     inline double getTHAWEND1( void ) { return thawend1; };


    // thawend2 ************************************************
      
     inline double getTHAWEND2( void ) { return thawend2; };


    // tsoil ***************************************************
      
     inline double getTSOIL( void ) { return tsoil; };


    // vcond ***************************************************

     inline double getVCOND( const int& pcmnt,
                             const int& pslayr ) 
     { 
       return vcond[pcmnt][pslayr]; 
     };

     inline void setVCOND( const double& pvcond, 
                           const int& pcmnt,
                           const int& pslayr ) 
     { 
       vcond[pcmnt][pslayr] = pvcond; 
     };


    // vegmat **************************************************

     inline integer getVEGMAT( const int& pcmnt,
                               const int& pslayr ) 
     { 
       return vegMAT[pcmnt][pslayr]; 
     };

     inline void setVEGMAT( const integer& pvegmat, 
                            const int& pcmnt,
                            const int& pslayr ) 
     { 
       vegMAT[pcmnt][pslayr] = pvegmat; 
     };


    // vegwindsp ***********************************************

     inline integer getVEGWINDSP( const int& pcmnt ) 
     { 
       return vegWindSp[pcmnt]; 
     };

     inline void setVEGWINDSP( const integer& pvegwindsp, 
                               const int& pcmnt ) 
     { 
       vegWindSp[pcmnt] = pvegwindsp; 
     };


    // vsph ***************************************************

     inline double getVSPH( const int& pcmnt,
                            const int& pslayr ) 
     { 
       return vsph[pcmnt][pslayr]; 
     };

     inline void setVSPH( const double& pvsph, 
                          const int& pcmnt,
                          const int& pslayr ) 
     { 
       vsph[pcmnt][pslayr] = pvsph; 
     };


    // water9 **************************************************

     inline double getWATER9( const int& pnode ) 
     { 
       return water9[pnode]; 
     };

     inline void setWATER9( const double& pwater9, 
                            const int& pnode ) 
     { 
       water9[pnode] = pwater9; 
     };


    // weight9 *************************************************

     inline double getWEIGHT9( const int& pnode ) 
     { 
       return weight9[pnode]; 
     };

     inline void setWEIGHT9( const double& pweight9, 
                             const int& pnode ) 
     { 
       weight9[pnode] = pweight9; 
     };


    // x9 ******************************************************

     inline double getX9( const int& pnode ) 
     { 
       return x9[pnode]; 
     };

     inline void setX9( const double& px9, 
                        const int& pnode ) 
     { 
       x9[pnode] = px9; 
     };
    

    // xfa9 ****************************************************

     inline double getXFA9( const int& pnode ) 
     { 
       return xfa9[pnode]; 
     };

     inline void setXFA9( const double& pxfa9, 
                          const int& pnode ) 
     { 
       xfa9[pnode] = pxfa9; 
     };


    // xfa9 ****************************************************

     inline double getXFB9( const int& pnode ) 
     { 
       return xfb9[pnode]; 
     };

     inline void setXFB9( const double& pxfb9, 
                          const int& pnode ) 
     { 
       xfb9[pnode] = pxfb9; 
     };


 /* ************************************************************
		 Public Variables
************************************************************* */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif


    // Annual sum of dst0
    double yrdst0; 
      
    // Annual sum of dst5
    double yrdst5;    

    // Annual sum of dst10
    double yrdst10;

    // Annual sum of dst20
    double yrdst20;   

    // Annual sum of dst50
    double yrdst50;   

    // Annual sum of dst100
    double yrdst100;  

    // Annual sum of dst200
    double yrdst200;  

    // Annual sum of frontd
    double yrfrontd;

    // Annual sum of frontd2
    double yrfrontd2;

    // Annual sum of thawbegin1
    double yrthawbegin1;

    // Annual sum of thawbegin2
    double yrthawbegin2;

    // Annual sum of thawend1
    double yrthawend1;

    // Annual sum of thawend2
    double yrthawend2;

    // Annual sum of tsoil
    double yrtsoil;


  private:

     enum interpIndex { SPACEINDX, TIMEINDX };

/* **************************************************************
		 Private Functions
************************************************************** */

    inline integer abs( integer x ) 
                   {
                     if ( x >= 0 ) { return x; }
                     else { return -(x); }
                   };

        // crankNicholsonBackward() ormerly TABV()

    int crankNicholsonBackward( integer *i1 );


        // crankNicholsonForward() formerly TBLO()

    int crankNicholsonForward( integer *i1 );


    void createMorePhasePlanes( integer& itzero,
                                integer& ibzero,
                                double& r1, 
                                double& topold,
                                double& topnew );


    void createSinglePhasePlane(  const integer& cmnt,
                                  const integer& itzero, 
                                  double& r1, 
                                  double& topold, 
                                  const double& topnew );
    
    
    inline doubledouble dabs( doubledouble x ) 
                        { 
                          if ( x >= 0 ) { return x; }
                          else { return -(x); }
                        };

    inline doubledouble dmax( doubledouble a, doubledouble b ) 
                        { 
                          if ( a >= b ) { return a; }
                          else { return b; } 
                        };
    inline doubledouble dmin( doubledouble a, doubledouble b ) 
                        { 
                          if ( a <= b ) { return a; }
                          else { return b; } 
                        };


        // gaussianEliminBackward() formerly ASMBLO()

    int gaussianEliminBackward( integer *i1, 
                                const integer& cmnt );


        // gaussianElimin() formerly ASMABV()

    int gaussianEliminForward( integer *i1 );


    void initMonthlySoilConditions( const int& outmon,
                                    const double& prevtair,
                                    const double& tair,
                                    const double& nexttair, 
                                    const double& prevspack, 
                                    const double& snowpack, 
                                    const double& nextspack,
                                    const int& cmnt,
                                    ofstream& rflog1 );

 
        // initSoilTemp() formerly TINITL
 
    int initSoilTemp( const int& cmnt, ofstream& rflog1 );


        // interpolate() formerly INTRPL()

    int interpolate( double x[], 
                     double y[], 
                     integer *nmax,
  	             double xx[], 
  	             double yy[], 
  	             integer *nnmax, 
                     const int& spaceTimeIndex,
  	             ofstream& rflog1 );


        // interpolateAirHeat(), formerly DATA(),  

    int interpolateAirHeat( double fitted[MXTSFITTED],
                            ofstream& rflog1 );

        // interpolateSnow(), formerly part of SNOFAL(), 
        //   interpolates the change in snow cover during each 
        //   time step from monthly changes in snow pack
        
    int interpolateSnow( const integer& maxTimeSteps,
                         const double& snowDensity ); 


    inline integer max( integer a, integer b )
                   {
                     if ( a >= b ) { return a; }
                     else { return b; }
                   };

    inline integer min( integer a, integer b )
                   {
                     if ( a <= b ) { return a; }
                     else { return b; }
                   };

    int resetNodeHeat( void );

    void resetTimeStep( const integer& cmnt, 
                        const int& outmon,
                        ofstream& rflog1 );
    
    int resetXFA( const double& thresholdFreezeT,  
                  const double& latentHeat );

    int setAirFreezeThawIndex( void );


        // setFullSnowNode(), formerly part of SNOFAL(),
        //   determines the maximum mass of snow that a "full"
        //   snow node can hold (NOTE: During simulations, the
        //   mass of a "full" snow node remains constant whereas 
        //   the thickness and density of snow within the node 
        //   are allowed to vary

    double setFullSnowNode( const double& timeStepSecs ); 


        // setHeatSources() formerly SETHET()

    int setHeatSources( integer *mmax, 
                        integer kint[2], 
                        double heat[MXTSFITTED][2],
                        double source[MXTSFITTED], 
                        ofstream& kflog1 );

    int setInitSnowPackGrid( const double& initSnowPack );


        // setLowerBoundaryConditions() formerly BCBOT()

    int setLowerBoundaryConditions( double *tber );


    int setMonthlySoilTemp( void );

    int setSnowCharSturm( double *ta, 
                          double *sph, 
                          double *cnd );


         // setSoilProfileGrid() formerly GRID()

    int setSoilProfileGrid( const int& cmnt, 
                            ofstream& rflog1 );  


        // setUpperBoundaryConditions() formerly BCTOP()

    int setUpperBoundaryConditions( double taer[MXTSFITTED], 
                                    ofstream& rflog1 ); 

    void updateActiveLayer( void );

        // updateMultiplePhasePlaneTemps() formerly ASMTWO()

    int updateMultiplePhasePlaneTemps( integer *itzero, 
                                       integer *ibzero,
                                       const integer& cmnt,
                                       ofstream& rflog1 );


        //  updatePhasePlane() formerly PHASE()

    int updatePhasePlane( integer *dnode, 
                          double *rlw, 
                          double *cu, 
                          double *rsu, 
                          double *cd, 
                          double *rsd, 
                          double *tt, 
                          double *td, 
                          double *dxx, 
                          double *gold, 
                          double *g, 
                          double *t1, 
                          double *t2, 
                          ofstream& rflog1 );


        // updateSinglePhasePlaneTemps() formerly ASMONE()

    int updateSinglePhasePlaneTemps( integer *dnode, 
                                     const integer& cmnt, 
                                     ofstream& rflog1 );


        // updateSnowPackGrid() formerly NEIGE()

    int updateSnowPackGrid( const double& acum, 
                            double *snoden, 
                            double *weight, 
                            double *comp, 
                            double *fcmelt, 
                            double *tdrive );


        // updateSoilProfileHeat() formerly INHET()

    int updateSoilProfileHeat( const int& m, 
                               integer *kint, 
                               double *source );


    int updateSoilThermalProperties( integer *mat, 
                                     double *ta, 
	                             double *sph, 
                                     double *cnd, 
                                     ofstream& rflog1 );


    int updateTimeStep( const int& m,
                        const integer& cmnt,
                        double& r1, 
                        const integer& imaxp,
                        ofstream& rflog1 );


/* **************************************************************
		 Private Variables
************************************************************** */

    // Width ( m ) of the active layer
    double activeLayer;

    // Variable used to calculate snow thermal conductivity 
    double calcu_cdsnow; 

    // Monthly soil temperature ( degrees C ) at the surface
    double dst0;

    // Monthly soil temperature ( degrees C ) at 5 cm
    double dst5;

    // Monthly soil temperature ( degrees C ) at 10 cm
    double dst10;

    // Monthly soil temperature ( degrees C ) at 20 cm
    double dst20;     

    // Monthly soil temperature ( degrees C ) at 50 cm
    double dst50;     

    // Monthly soil temperature ( degrees C ) at 100 cm
    double dst100;    

    // Monthly soil temperature ( degrees C ) at 200 cm
    double dst200;    

    // Monthly soil temperature ( degrees C ) at 300 cm
    double dst300;
    // Initial or previous thickness of each node ( m ) in the
    //   snowpack and soil
    double dx9[MAXNODES+1];

    // Monthly depth ( m ) of the node containing the upper 
    //   phase change plane (frozen front from ground surface)
    double frontd;

    double frontd2;

    // Initial or previous index of node representing the 
    //   surface of the snowpack    
    int is9;

    // Debugging flag
    int kdbgflg;

    // Flag to indicate that the soil thermal model needs to
    //  be initialized:
    //  kswitch = 0 during initialization
    //  kswitch = 1 when soil thermal model has already been 
    //              initialized
    int kswitch; 

    double nextdst10; // Next month's dst10

    double nextsnowfall;     // Next month's snowfall
    
    double nextspack;        // Next month's snowpack

    double nexttair;         // Next month's air temperature

    double prevdst10; // Previous month's dst10
    
    double prevspack;        // Previous month's snowpack

    // Initial or previous mass of snowpack ( kg m^-2 )
    double smass9;

    // Monthly snow density ( g cm^-3 )
    double snow_dens;
    
    // Initial or previous temperature at each node  
    double t9[MAXNODES+1];

    // Monthly depth ( m ) of thawed soil just below 
    //   the upper phase change plane associated with
    //   the freezing front from top of soil profile
    double thawbegin1;

    // Monthly depth ( m ) of thawed soil just below 
    //   lower phase change plane associated with 
    //   frozen soil lens in middle of profile
    //   (if it exists)
    double thawbegin2;

    // Monthly depth ( m ) of phase change plane 
    //   associated with frozen soil lens in middle
    //   of profile (if it exists)
    double thawend1;

    // Monthly depth ( m ) of the node just above the lower 
    //   phase change plane represented by permafrost 
    //   (if it exists)
    double thawend2;

    // Mean monthly soil temperature ( degrees C ) of the top 
    //   20 cm of the soil profile
    double tsoil;      

    // Initial or previous water content at each node
    double water9[MAXNODES+1];
    
    // Initial or previous mass of snow within each node 
    //   of the snowpack
    double weight9[MAXSNODES];

    // Initial or previous depth of the top of each node
    double x9[MAXNODES+1];   

    // Initial or previous depth of phase change plane from the 
    //   top of the node that contains the phase change plane 
    double xfa9[MAXNODES+1];

    // Initial or previous heighth of the phase change plane 
    //   above the bottom of the node that contains the phase 
    //   change plane
    double xfb9[MAXNODES+1];




    // Interpolated air temperatures for the beginning, 
    //   middle and end of a month   
    STMmonth airheat;  

    // The term "CAPi" in Figure 5.2 of Goodrich [1976] and 
    //   Figure 1 in Figure 1 of Goodrich [1978b] for each
    //   node in the snowpack and soil
    double capx[MAXNODES+1];

    // Cooling degree month ( degrees C )
    //   see p.1274 in Sturm et al. (1995)
    double CDM; 

    // Maximum Mass ( g m^-2? )of snow in a "filled" snow node 
    double comp;
     
    // Frozen thermal conductivity ( W m^-1 K^-1 ) for each
    //   material type  
    double condf[MAXMAT];
    
    // Thawed thermal conductivity ( W m^-1 K^-1 ) for each
    //   material type   
    double condt[MAXMAT];
    
    // The term "CNi" in Figure 5.2 of Goodrich [1976] and 
    //   Figure 1 in Figure 1 of Goodrich [1978b] for each
    //   node in the snowpack and soil
    double conx[MAXNODES+1];

    // Number of days per month
    double daze[CYCLE];
    
    // Dry density ( kg m^-3 ) of each node in the snowpack
    //   and soil 
    double ddry[MAXNODES+1];  
    
    // Actual time step used ( seconds )
    double dt;       
    
    // Actual time step used ( days )
    double dtday;
    
    // Amount of time ( seconds) the phase change plane existed 
    //   in the node (NOTE: will not be larger than dt)
    double dtfaz;
    
    // Thickness of each node ( m ) in the snowpack and soil 
    double dx[MAXNODES+1];    
        
    // Thickness ( m ) of the current node
    double dxx;

    // Coefficient in Crank-Nicholson finite difference method
    //   (see Equation 5.13 in Goodrich [1976] or Equations [7a]
    //   and [9a] in Goodrich [1978b])
    double e[MAXNODES+1];

    // Interpolated time ( days ) for a time step within a 
    //   particular month
    double eltim[MXTSFITTED];
//    double eltim[150000];
    
    // Mass of snow (kg m^-3) that melts per day per degree C
    double fcmelt;
    
    // final day of time period
    double final;   
    
    // Air freezing index 
    double frza;

    // Heat at the bottom of the soil profile
    double hbot;

    // Interpolated heat from up to two heat sources during
    //   a 0.5-day time step 
    double heat[MXTSFITTED][2];
    
    // Source of heat
    double heatt[2];

    // Latent heat
    double hlat;    

    // Heat at each node in the soil profile
    double ht[MAXNODES+1];    
    
    // Previous heat at each node in the soil profile
    double htold[MAXNODES+1]; 
    
    // Heat above the surface of the snowpack or soil 
    //   (if snowpack does not exist)
    double htop;   

    // Index of node representing the ground surface 
    //   (i.e., maximum number of nodes in snowpack)
    integer ig;      
    
    // Index of node just above the ground surface
    //   (i.e., bottom node of snowpack if it exists)
    integer igm1;

    // Index of the node that would represent the bottom of the
    //   bottom node of the soil profile (i.e., the bottom of 
    //   the soil profile
    integer imax;
    
    // Index of the node that represents the top of the bottom
    //   node of of the soil profile
    integer imax1;

    // Index of node representing the surface of the
    //   snowpack if it exists
    integer is;   
    
    // Index of node just above the surface of the 
    //   snowpack (i.e. air)
    integer ism1;
    
    // Initial or previous index of node just above the surface 
    //   of the snowpack 
    integer ism19;

    integer kint[2];
    
    // Material index of each node:
    //    0 - Missing data
    //    1 - Moss/litter conductivity data read in from a file
    //    2 - Organic soil conductivity data read in from a file
    //    3 - Mineral soil conductivity data read in from a file
    //    4 - Sturm snow 
    integer mater[MAXNODES+1];
    
    // Maximum number of time steps in a month
    integer mmax;            
   
    // Period of calculation cycle (days)
    double per;     

    // Coefficient in Crank-Nicholson finite difference method
    //   (see Equation 5.13 in Goodrich [1976] or Equations 
    //   [7b] and [9b] in Goodrich [1978b])
    double s[MAXNODES+1];
    
    // Interpolated snow density ( kg m^-3 ) that occurs within 
    //   a 0.5 day time step    
    double sden[MXTSFITTED];
    
    // Heat flux from surface of snowpack or soil (if snowpack
    //   does not exist)
    double sflux;    
       
    // Mass of snowpack ( kg m^-2 )
    double smass;
        
    // Interpolated mass of snow ( kg m^-2 ) that accumulates
    //   or declines within a 0.5 day time step 
    double snow[MXTSFITTED];

    // Interpolated snow pack depth ( m ) for the beginning, 
    //   middle and end of a month
    STMmonth snowdepth; 
    
    // Water content of 3 hydrologic soil layers of the HM
    STMsoilwater soilwater; 

//    double spare[MAXSLAYERS];
    
    // Frozen volumetric heat capacity ( MJ m^-1 K^-1 x 1000 )
    double sphf[MAXMAT];
    
    // Thawed volumetric heat capacity ( MJ m^-1 K^-1 x 1000 )
    double spht[MAXMAT];

    // Temperature (degrees C ) at each node
    double t[MAXNODES+1];   
    
    // Interpolated air temperature ( degrees C ) for a 
    //   0.5 day time step
    double taer[MXTSFITTED];
    
    // Air temperature ( degrees C ) during current 0.5-day 
    //   time step
    double tstair;
    
    // Highest air temperature ( degrees C ) for a time step 
    //   during a month
    double tairhi;
    
    // Lowest air temperature ( degrees C ) for a time step 
    //   during a month
    double tairlo;
    
    // Mean monthly temperature ( degrees C ) at each node
    //   in the snowpack or soil
    double tanee[MAXNODES+1];  
    
    // Interpolated temperature ( degrees C ) at the bottom
    //   of the soil profile for a 0.5 day time step
    double tber[MXTSFITTED];
    
    // Temperature ( degrees C ) at the bottom of the 
    //   soil profile
    double tbot;
    
    // Temperature ( degrees C ) of the node below the 
    //   current node
    double td;

    // Threshold temperature ( degrees C ) for phase change 
    //   of water/ice
    double tf;       

    // Previous temperature ( degrees C ) at each node in
    //   the snowpack and soil
    double thalt[MAXNODES+1];
    
    double theta;           // theta parameter
    double theta1;
    
    // Air thawing index
    double thwa;
    
    double time;    // cumulative time of calculation
    double told[MAXNODES+1];  // previous temperature at each node

    // Variable used to indicate an assumed constant 30 days 
    //   per month when calculating mmax and looking at 
    //   within month thermal and snow dynamics  
    double total;
    
    // Water content ( volumetric % ) at each node
    double water[MAXNODES+1];

    double water1; // moss/litter layer water content
    double water2; 
    double water3;

    // Actual mass ( kg m^-2 ) of snow within each node of the 
    //   snowpack
    double weight[MAXSNODES];

    // Depth ( m ) of the top of each node 
    //   (negative if in snowpack; positive if in soil)
    double x[MAXNODES+1];    
    
    
    // Depth ( m ) of the phase change plane from the top of 
    //   the node that contains the phase change plane 
    //   (NOTE: any node may contain a phase change plane;
    //   xfa[] must be less than the thickness of the node
    double xfa[MAXNODES+1];
    
    
    // Heighth ( m ) of the phase change plane above the bottom 
    //   of the node that contains the phase change plane 
    //   (NOTE: any node may contain a phase change plane;
    //   xfb must be less than the thickness of the node
    double xfb[MAXNODES+1];
    
    
    // Depth ( m ) of heat source from the top of the node that
    //   contains the heat source (NOTE: All nodes within the 
    //   soil may contain a heat source
    double xhet[MAXNODES+1];

    // ********************* Parameters ************************

    // Vegetation-specific
    double cdsnow[MAXCMNT];

    // Vegetation-specific depths of initial soil temperatures
    double DEPTH[MAXCMNT][MAXSPDATA];

    // Vegetation-specific
    double gflux[MAXCMNT];

    // Vegetation-specific Frozen thermal conductivity of each 
    //   soil layer    
    double initCONDF[MAXCMNT][MAXSLAYERS];

    // Vegetation-specific thawed thermal conductivity of 
    //   soil layer    
    double initCONDT[MAXCMNT][MAXSLAYERS];

    // Vegetation-specific dry density of soil layer ( kg m^-3 )
    double initDENSE[MAXCMNT][MAXSLAYERS];
    
    // Vegetation-specific time step as proportion of a day
    double initDTDAY[MAXCMNT];   
    
    // Vegetation-specific depth step ( m ) of heat flux from 
    //   top within soil layer
    double initDXA[MAXCMNT][MAXSLAYERS];
    
    // Vegetation-specific depth step ( m ) of heat flux from 
    //   bottom within soil layer
    double initDXB[MAXCMNT][MAXSLAYERS];
    
    // Vegetation-specific frozen volumetric heat capacity of 
    //   soil layer
    double initSPHF[MAXCMNT][MAXSLAYERS];
    
    // Vegetation-specific thawed volumetric heat capacity 
    //   of soil layer
    double initSPHT[MAXCMNT][MAXSLAYERS];

   
    // Vegetation-specific thickness ( m ) of soil layers:
    //   initTHICK[MAXCMNT][0] - moss/litter layer
    //   initTHICK[MAXCMNT][1] - upper organic soil layer
    //   initTHICK[MAXCMNT][2] - lower organic soil layer
    //   initTHICK[MAXCMNT][3] - upper mineral soil layer
    //   initTHICK[MAXCMNT][4] - lower mineral soil layer
    
    double initTHICK[MAXCMNT][MAXSLAYERS];
    
    // Vegetation-specific initial water content of each 
    //   soil layer
    double initWATER[MAXCMNT][MAXSLAYERS];

    // Vegetation-specific initial temperature data
    double TEMP[MAXCMNT][MAXSPDATA];
    
    // Vegetation-specific
    double vcond[MAXCMNT][MAXSLAYERS];
   
    // Vegetation-specific
    double VDEPTH1[MAXCMNT];
   
    // Vegetation-specific material index of each soil layer
    integer vegMAT[MAXCMNT][MAXSLAYERS];
    
    // Vegetation-specific wind speed index: 
    //   0 for low
    //   1 for high 
    integer vegWindSp[MAXCMNT];

    // Vegetation-specific
    double vsph[MAXCMNT][MAXSLAYERS]; 
    
};

#endif

