/* *************************************************************
****************************************************************
QSOILTEMP602.CPP - Soil temperature model developed by Qianlai
                    Zhuang for the Terrestrial Ecosystem Model
****************************************************************

Modifications:

20040421 - DWK created by modifying qsoiltemp50b1.cpp
20040421 - DWK deleted const double& xtair and 
           const int& xtairflg from initMonthlySoilConditions() 
           and setMonthlySoilConditions()
20040421 - DWK changed x{MXTSDATA} to x[], y[MXTSDATA] to y[],
           xx[MXTSFITTED] to xx[], and yy[MXTSFITTED] to yy[]
           in function call to interpolate()
20040421 - DWK renamed updateSoilTemperature() to be 
           interpolateAirHeat() and eliminated integer *index
           from function call
20040421 - DWK renamed makeSnowPack() to be setSnowPack()
20040421 - DWK renamed setHeatFlux() to be setHeatSources()
20040421 - DWK changed void setSnowDensity() to 
           double setSnowDensity()
20040423 - DWK renamed setFreezeThawIndex() to be 
           setAirFreezeThawIndex()
20040423 - DWK renamed updateTempDownsweep() to be
           crankNicholsonBackward()
20040423 - DWK renamed updateTempUpsweep() to be
           crankNicholsonForward()
20040423 - DWK added private function setSnowCharSturm()
20040423 - DWK changed double tstart[MAXNODES] to 
           double tstart[MAXSPDATA] in initSoilTemp() 
20040423 - DWK changed double xstart[MAXNODES] to 
           double xstart[MAXSPDATA] in initSoilTemp() 
20040424 - DWK added setFullSnowNode() and interpolateSnow() 
           to replace makeSnowPack()
20040424 - DWK changed airheat.previous to be airheat.beginning
20040424 - DWK changed airheat.current to be airheat.middle
20040424 - DWK changed airheat.next to be airheat.end
20040424 - DWK changed snowdepth.previous to be 
           snowdepth.beginning
20040424 - DWK changed snowdepth.current to be snowdepth.middle
20040424 - DWK changed snowdepth.next to be snowdepth.end
20040424 - DWK added const double& thresholdFreezeT and 
           const double& latentHeat to resetXFA()
20040424 - DWK added const int& spaceTimeIndex to interpolate()
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
20040425 - DWK added resetTimeStep(), createSinglePhasePlane() 
           and createMorePhasePlanes()
20040425 - DWK renamed setSinglePhasePlane() as 
           updateSinglePhasePlaneTemps()
20040425 - DWK renamed setMultiplePhasePlanes as
           updateMultiplePhasePlaneTemps()
20040425 - DWK changed tair to tstair in setAirFreezeThawIndex(),
           updateSnowPackGrid(), updateSurfaceHeatFlux() and
           updateTimeStep()
20040426 - DWK added private functions equilibriumCheck() and
           updateActiveLayer()
20040426 - DWK deleted diffsoilt[] from setMonthlySoilTemp()
20040427 - DWK added setSnowPackGrid()
20040428 - DWK deleted const int& outmon from function call of 
           setSoilProfileGrid(() 
20040501 - DWK deleted updateSurfaceHeatFlux()
20040502 - DWK multiplied snow_dens by 1000.0 to change the 
           units from g cm^-3 to kg m^-3 in setSnowDensity().  
           In addition, DWK changed (1000.0 * snow_dens) to
           snow_dens in denominator of algorithms calculating
           snowdepth.beginning, snowdepth.middle and 
           snowdepth.end
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
           and integer KINTI[MAXCMNT] from functions
20040504 - DWK deleted private integer INDEX{MAXCMNT] and 
           double VDEPP[MAXCMNT] from functions
20040504 - DWK added initSoilThermalRegime()
20040504 - DWK deleted integer *ktop from 
           setUpperBoundaryConditions()
20040505 - DWK moved initSoilTemp() and resetXFA() from
           initMonthlySoilConditions() to initSoilThermalRegime()
20040505 - DWK moved initialization of tf and hlat from
           initMonthlySoilConditions() to initSoilThermalRegime()
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
20040515 - DWK replaced calcu_cdsnow with the Sturm et al. (1997)
           equation in setSnowCharSturm()
20040515 - DWK replaced conditional statement "dnode > lmimax"
           with "dnode == lmimax" in setSoilProfileGrid()
20040515 - DWK replaced array size MAXNODES+2 with MAXNODES+1
20040517 - DWK deleted call to setSnowPackGrid() from
           initMonthlySoilConditions()  
20040517 - DWK deleted const double& initSnowPack from from  
           setSnowPackGrid() and added const double& prevspack, 
           const double& snowpack and const double sdensity to  
           setSnowPackGrid() and renamed it as 
           setInitSnowPackGrid()
20040707 - DWK changed Soilthermal:: to Soilthermal60::
20040707 - DWK changed ProcessXML() to ProcessXML51() in
           Soilthermal60()
20040707 - DWK changed void getsnowecd( char ecd[MAXFNAME] ) to
           void getsnowecd( string ecd )
20040707 - DWK changed void getsoillecd( char ecd[MAXFNAME] ) to
           void getsoillecd( string ecd )
20040707 - DWK changed void getsoiltecd( char ecd[MAXFNAME] ) to
           void getsoiltecd( string ecd )
20051117 - DWK added include qsoiltem602.cpp and standard includes
20051213 - DWK added public function resetYrFluxes()
20060105 - DWK added public function InterpolatST()
20060109 - DWK added const double& tstart to 
           initSoilThermalRegime()
                                 
**************************************************************** 
************************************************************* */

#include<cstdio>

  using std::printf;
  using std::sprintf;

#include<iostream>

  using std::cin;
  using std::cout;
  using std::ios;
  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std:: atoi;
  
#include<string>
  
  using std::string;

#include<cmath>

  using std::fabs;
  using std::log;
  using std::pow;
  using std::sqrt;
  
//#include<values.h>

#include "qsoiltemp603.h"

/* *************************************************************
************************************************************* */

Soilthermal60::Soilthermal60() : ProcessXML60() 
{ 

// Initialize number of days per month for each month of year

  daze[0] = daze[2] = daze[4] = daze[6] = daze[7] = 31.0;
  daze[9] = daze[11] = 31.0;
  daze[3] = daze[5] = daze[8] = daze[10] = 30.0;
  daze[1] = 28.0;

  is9 = 9;
  ism19 = is9 - 1;
  smass9 = ZERO;

  kswitch = 0;

  kdbgflg = 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::crankNicholsonBackward( integer *i1 )
{

  // formerly TABV(I1)
  // Evaluates new node temperatures between is and i1
  // Backward substitution for t (see Equation [8] in Goodrich
  // [1978b])
  
  integer dnode;
  integer j;
  integer jj;

  if( *i1 < is ) { return 0; }
  jj = is + *i1;

  for( j = is; j <= *i1; ++j ) 
  {
    dnode = jj - j;
    t[dnode] = s[dnode] * t[dnode+1] + e[dnode];
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::crankNicholsonForward( integer *i1 )
{

  // formerly TBLO(I1)
  // Evaluates new node temperatures between i1 and imax
  // Forward substitution for t (see Equation [10] in Goodrich
  // [1978b]  
  
  integer dnode;

  if( *i1 > imax ) { return 0; }

  for( dnode = *i1; dnode <= imax; ++dnode ) 
  {
    t[dnode] = s[dnode] * t[dnode-1] + e[dnode];
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::createMorePhasePlanes( integer& itzero,
                                           integer& ibzero,
                                           double& r1, 
                                           double& topold,
                                           double& topnew )
{
  integer dnode;
  double tset;
  
      r1 = topnew - tf;

      if( r1 >= ZERO ) { tset = tf - ADJUST; }
      else { tset = tf + ADJUST; }

      for( dnode = is; dnode <= itzero; ++dnode ) 
      {
        e[dnode] = ZERO;
        s[dnode] = ZERO;
        t[dnode] = tset;  // ***** should be careful 
      }
   
      // to protect

      if( ZERO == (topnew - topold) ) { topnew = topold + 0.01; }

      dtfaz = dt * (topnew - tf) / (topnew - topold);
      ibzero = itzero;

      // TO ENSURE 2ND CALL TO updateMultiplePhasePlaneTemps()
      //   IGNORES LOWER PHASE PLANE -- IN CASE PREVIOUSLY HAD
      //   MULTIPLE PHASE PLANES */

      if( ibzero > itzero ) { ibzero = 1; }
      itzero = is;
      r1 = topnew - tf;
      
      if( r1 >= ZERO ) { topold = tf + ADJUST; }
      else { topold = tf - ADJUST; }

      t[is] = topold;
      
      if( xfa[is] > ZERO ) 
      {
        xfb[is] = xfa[is];
      }

      xfa[is] = ZERO;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::createSinglePhasePlane( const integer& cmnt,
                                            const integer& itzero, 
                                            double& r1, 
                                            double& topold, 
                                            const double& topnew )
{

  integer itz1;
  
  xfa[itzero] = ZERO;

  if( per > ZERO ) 
  { 
    // to protect
    if ( ZERO == (topold - topnew) ) 
    { 
      topold =  topnew + 0.01; 
    } 

    // Interpolate the amount of time that the phase change 
    //   plane existed within the time step
     
    dtfaz = dt * (topold - tf) / (topold - topnew);
    
    itz1 = itzero + 1;
    
    // Set temperature of the surface node (i.e., is) to be
    //   the threshold freezing temperature
    
    t[itzero] = tf;

    // Determine Crank-Nicholson coefficients e[] and s[] for 
    //   all nodes below the surface node in the snowpack 
    //   and soil

    gaussianEliminBackward( &itz1, cmnt );

    // Recalculate temperatures for all nodes below the surface 
    //   node (i.e. is) in the snowpack and soil

    crankNicholsonForward( &itz1 );

    r1 = topnew - tf;
    
    if( r1 >= 0 ) { topold = tf + ADJUST; }
    else { topold = tf - ADJUST; }

    t[itzero] = topold;
    dtfaz = dt - dtfaz;
  }
  else
  {
    // TRANSIENT PROBLEMS ONLY 

    topold = topnew;
    t[itzero] = topold;
  }
  

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::gaussianEliminBackward( integer *i1, 
                                           const integer& cmnt )
{

// formerly ASMBLO(Il)

// Assembles e[i], s[i] for ordinary nodes (imax to i1) for 
//   upsweep e, s (see Equations [9a] and [9b] in Goodrich 
//   [1978b]) to be used in the Crank-Nicholson finite 
//   difference approach
   

  double con;
  double cond;
  double condt1;
  double condth;

  double conu;
  double conut1;
  double conuth;

  double denm;

  integer dnode;
  integer i2;
  integer im1;
  integer ip1;

  integer j;
  integer jj;

  double r;
  double rc;
  double rhs;

  // Bottom node ordinary implicit
  
  con = conx[imax1];

  // to protect
  if( ZERO == dt ) { dt = 0.01; } 

  rc = capx[imax1] * 0.5 / dt;
  denm = rc + con + hbot;

  // to protected   debugged Q. Z. 12/Dec/2000; 
  if( ZERO == denm ) { denm = 30.0; } 

  s[imax] = con / denm;
  e[imax] = (rc * t[imax] - gflux[cmnt] + (hbot * tbot) 
              + ht[imax]) / denm;
  if( *i1 >= imax ) { return 0; }

  // Theta method
  
  i2 = *i1;
  if( *i1 == is ) { i2 = is + 1; }
  jj = i2 + imax1;

  for( j = i2; j <= imax1; ++j ) 
  {
    dnode = jj - j;
    im1 = dnode - 1;
    ip1 = dnode + 1;

    if( ZERO == dt ) { dt = 0.01; } // to protect

    rc = (capx[dnode] + capx[im1]) / dt;
    conut1 = conx[im1] * theta1;
    condt1 = conx[dnode] * theta1;
    rhs = (rc - conut1 - condt1) * t[dnode] + conut1 
          * t[im1] + condt1 * t[ip1];
    rhs = rhs + theta * ht[dnode] + theta1 * htold[dnode];
    conuth = conx[im1] * theta;
    condth = conx[dnode] * theta;
    r = conuth + rc + condth;
    denm = r - condth * s[ip1];

    // to protected   debugged Q. Z. 12/Dec/2000;
    if( denm == 0.0 )  { denm = 30.0; } 
    
    s[dnode] = conuth / denm;
    e[dnode] = (rhs + condth * e[ip1]) / denm;
  }
  
  if( *i1 > is ) { return 0; }

  // Surface node ordinary implicit - prescribed flux permitted
  
  // to protect
  if( ZERO == dt ) { dt = 0.01; }; 
  
  rc = capx[is] * 0.5 / dt;
  conu = conx[ism1];
  cond = conx[is];
  rhs = rc * t[is] + sflux + ht[dnode];
  r = conu + rc + cond;
  denm = r - cond * s[i2];

  // to protected   debugged Q. Z. 12/Dec/2000;
  if( ZERO == denm ) { denm = 30.0; } 
  
  s[is] = conu / denm;
  e[is] = (rhs + cond * e[i2]) / denm;

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::gaussianEliminForward( integer *i1 )
{

// formerly ASMABV(Il)

// Assembles e[i], s[i] for ordinary nodes (is to i1) for 
//   downsweep e,s (see Equations [7a] and [7b] in Goodrich  
//   [1978b]) to be used in the Crank-Nicholson finite 
//   difference approach

  double cond;
  double condt1;
  double condth;

  double conu;
  double conut1;
  double conuth;

  double denm;

  integer dnode;
  integer im1;
  integer ip1;
  integer isp1;

  double r;
  double rc;
  double rhs;

  if( *i1 < is ) { return 0; }

  // Surface node ordinary implicit - prescribed flux permitted

  // to protect
  if( ZERO == dt ) { dt =0.01; }; 

  rc = capx[is] * 0.5 / dt;

  conu = conx[ism1];
  cond = conx[is];
  rhs = rc * t[is] + sflux + ht[is];
  r = conu + rc + cond;
  denm = r - conu * s[ism1];
 
 // to protected - debugged Q. Z. 12/Dec/2000;  
  if( denm == 0.0 ) { denm = 30.0; } 

  s[is] = cond / denm;

  e[is] = (rhs + conu * e[ism1]) / denm;
  if( *i1 == is ) { return 0; }

  // Theta method
  
  isp1 = is + 1;

  for( dnode = isp1; dnode <= *i1; ++dnode ) 
  {
    im1 = dnode - 1;
    ip1 = dnode + 1;

    if( dt == 0.0 ) { dt = 0.01; } // to protect
    
    rc = (capx[dnode] + capx[im1]) / dt;
    conu = conx[im1];
    cond = conx[dnode];
    conut1 = theta1 * conu;
    condt1 = theta1 * cond;
    rhs = (rc - conut1 - condt1) * t[dnode] + conut1 *
   	  t[im1] + condt1 * t[ip1];
    rhs = rhs + theta * ht[dnode] + theta1 
          * htold[dnode];
    conuth = theta * conu;
    condth = theta * cond;
    r = conuth + rc + condth;
    denm = r - conuth * s[im1];

    if( denm == 0.0 )  { denm = 30.0; } // to protect debugged Q. Z. 12/Dec/2000;
    
    s[dnode] = condth / denm;
    e[dnode] = (rhs + conuth * e[im1]) / denm;
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::getsnowecd( ofstream& rflog1 ) 
{

  string ecd;


#ifdef PMODE
  *fgo >> ecd;
#else
  cout << "Enter name of the snow (.ECD) data file with the ";
  cout << "parameter values:" << endl;
  cin >> ecd;
#endif

  rflog1 << "Enter name of the snow (.ECD) data file with the ";
  rflog1 << "parameter values: ";
  rflog1 << ecd << endl << endl;

  getsnowecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::getsnowecd( const string& ecd ) 
{

  ifstream infile;
//  int dnode;
  int dcmnt;
  int comtype;

//  char string[80];

  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for SNOWPACK ECD input" << endl;
    
    exit( -1 );
  }

  getXMLrootNode( infile, "snowPackECD" ); 


  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "snowPackECD" );
    
    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in snowPackECD" << endl;
      
      exit( -1 );
    }

    vegWindSp[comtype] = getXMLcmntArrayInt( infile, 
                                             "snowPackECD", 
                                             "vegWindSp", 
                                             comtype ); 
    
    gflux[comtype] = getXMLcmntArrayDouble( infile, 
                                         "snowPackECD", 
                                         "gflux",  
                                         comtype ); 

    cdsnow[comtype] = getXMLcmntArrayDouble( infile, 
                                         "snowPackECD", 
                                         "cdsnow",  
                                         comtype ); 

    initDTDAY[comtype] = getXMLcmntArrayDouble( infile, 
                                                "snowPackECD", 
                                                "DTDAY",  
                                                comtype ); 

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in snowPackECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::getsoillecd( ofstream& rflog1 )
{

  string ecd;

#ifdef PMODE
  *fgo >> ecd;
#else
  cout << "Enter name of soil layer (.ECD) data file with ";
  cout << "parameter values:" << endl;
  cin >> ecd;
#endif

  rflog1 << "Enter name of soil layer (.ECD) data file with ";
  rflog1 << "parameter values:" << ecd << endl << endl;

  getsoillecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::getsoillecd( const string& ecd )
{

  ifstream infile;
  int dlayer;
  int dcmnt;
  int comtype;

  char string[80];


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << "\nCannot open " << ecd;
    cerr << " for SOIL LAYER ECD input" << endl;
    
    exit( -1 );
  }

  getXMLrootNode( infile, "soilLayerECD" ); 

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "soilLayerECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in soilLayerECD" << endl;
      
      exit( -1 );
    }
    
    for( dlayer = 0; dlayer < MAXSLAYERS; ++dlayer )
    {
      sprintf( string, "%s%d", "THICK", (dlayer+1) );	
      initTHICK[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                          "soilLayerECD", 
                                                          string,  
                                                          comtype ); 

      sprintf( string, "%s%d", "DXA", (dlayer+1) );	
      initDXA[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                        "soilLayerECD", 
                                                        string,  
                                                        comtype ); 

      sprintf( string, "%s%d", "DXB", (dlayer+1) );	
      initDXB[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                        "soilLayerECD", 
                                                        string,  
                                                        comtype ); 

      sprintf( string, "%s%d", "MAT", (dlayer+1) );	
      vegMAT[comtype][dlayer] = getXMLcmntArrayLong( infile, 
                                                     "soilLayerECD", 
                                                     string,  
                                                     comtype ); 

      sprintf( string, "%s%d", "DENSE", (dlayer+1) );	
      initDENSE[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                          "soilLayerECD", 
                                                          string,  
                                                          comtype ); 

      sprintf( string, "%s%d", "WATER", (dlayer+1) );	
      initWATER[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                          "soilLayerECD", 
                                                          string,  
                                                          comtype ); 

      sprintf( string, "%s%d", "vcond", (dlayer+1) );	
      vcond[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                      "soilLayerECD", 
                                                      string,  
                                                      comtype ); 

      sprintf( string, "%s%d", "vsph", (dlayer+1) );	
      vsph[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                     "soilLayerECD", 
                                                     string,  
                                                     comtype ); 

      sprintf( string, "%s%d", "cond", (dlayer+1) );	
      initCONDT[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                          "soilLayerECD", 
                                                          string,  
                                                          comtype ); 

      sprintf( string, "%s%d", "spht", (dlayer+1) );	
      initSPHT[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                         "soilLayerECD", 
                                                         string,  
                                                         comtype ); 

      sprintf( string, "%s%d", "condf", (dlayer+1) );	
      initCONDF[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                          "soilLayerECD", 
                                                          string,  
                                                          comtype ); 

      sprintf( string, "%s%d", "sphf", (dlayer+1) );	
      initSPHF[comtype][dlayer] = getXMLcmntArrayDouble( infile, 
                                                         "soilLayerECD", 
                                                         string,  
                                                         comtype ); 
    }

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in soilLayerECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::getsoiltecd( ofstream& rflog1 )
{

  string ecd;

#ifdef PMODE
  *fgo >> ecd;
#else
  cout << "Enter name of the soil temperature initial (.ECD) ";
  cout << "data file with the parameter values:" << endl;
  cin >> ecd;
#endif

  rflog1 << "Enter name of the soil temperature initial (.ECD) ";
  rflog1 << "data file with the parameter values:" << ecd;
  rflog1 << endl << endl;

  getsoiltecd( ecd );


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::getsoiltecd( const string& ecd )
{
  ifstream infile;
  int dnode;
  int dcmnt;
  int comtype;

  char string[80];


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for soil temperature ECD input" << endl;
    
    exit( -1 );
  }

  getXMLrootNode( infile, "soilTempECD" ); 


  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "soilTempECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in soilTempECD" << endl;
      
      exit( -1 );
    }

    for( dnode = 0; dnode < 25; ++dnode )
    {
      sprintf( string, "%s%d", "DEPTH", (dnode+1) );	
      DEPTH[comtype][dnode] = getXMLcmntArrayDouble( infile, 
                                                     "soilTempECD", 
                                                     string,  
                                                     comtype ); 

      sprintf( string, "%s%d", "TEMP", (dnode+1) );	
      TEMP[comtype][dnode] = getXMLcmntArrayDouble( infile, 
                                                    "soilTempECD", 
                                                    string,  
                                                    comtype ); 
    }

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in soilTempECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::initMonthlySoilConditions( const int& outmon,
                                               const double& prevtair,
                                               const double& tair, 
                                               const double& nexttair, 
                                               const double& prevspack, 
                                               const double& snowpack, 
                                               const double& nextspack,
                                               const int& cmnt,
                                               ofstream& rflog1 )
{

  int dnode;

//  double ksnow_dens;

  // Initial snow pack depth at beginning of current month
  //   (added by DWK on 20040427)
  double initSnowPack;
  
  // Determine air temperatures ( degrees C ) at the beginning,
  //   middle and end of current month based on air temperatures 
  //   for the previous month, the current month, and the 
  //   next month

  airheat.beginning = (prevtair + tair)/ 2.0; 
  airheat.middle = tair;
  airheat.end = (tair + nexttair) / 2.0;

  // Determine snow depths ( m ) at the beginning, middle and 
  //   end of the current month based on the snowpack for the 
  //   previous month, the current month and the next month 
  //   (see Equation [2] in Zhuang et al., 2003)

  snowdepth.beginning = (prevspack + snowpack) 
                        / (2.0 * snow_dens); 
  
  snowdepth.middle = snowpack / snow_dens;
  
  snowdepth.end = (snowpack + nextspack) 
                  / (2.0 * snow_dens);

  
  // linkage of STM with HYM
  // Initialize soil water layers with current month PCTP

//  soilwater.layer1 = pctp / 100.0; 
//  soilwater.layer2 = pctp / 100.0;  
//  soilwater.layer3 = *water3;     


//  if( 0 == kswitch )
//  {
//    initSoilThermalRegime();
//  }

  
  for( dnode = 0; dnode  < MAXNODES; ++dnode ) 
  {
    x[dnode] = x9[dnode];
    dx[dnode] = dx9[dnode];
    xfa[dnode] = xfa9[dnode];
    xfb[dnode] = xfb9[dnode];
    xhet[dnode] = ZERO;

    water[dnode] = water9[dnode];
    ddry[dnode] = ZERO;

    conx[dnode] = ZERO;
    capx[dnode] = ZERO;
    e[dnode] = ZERO;
    s[dnode] = ZERO;
    t[dnode] = t9[dnode];
    ht[dnode] = ZERO;
    htold[dnode] = ht[dnode];
  }
  
  resetTimeStep( cmnt, outmon, rflog1 );

  // Determine soil profile properties based on gridded 
  //   depth steps (NOTE: Time step ( dtday, dt ) may change 
  //   within setSoilProfileGrid())

  setSoilProfileGrid( cmnt, rflog1 );


  // Determine upper boundary conditions for soil profile

  setUpperBoundaryConditions( taer, rflog1 );


  // Determine lower boundary conditions for soil profile

  setLowerBoundaryConditions( tber );

  
  // Determine heat sources within the soil profile

  setHeatSources( &mmax, kint, heat, heatt, rflog1 );


  
 
//  if( 0 == kswitch )
//  {
    // Determine initial gridded snowpack characteristics 

    initSnowPack = snowdepth.beginning * snow_dens;
    setInitSnowPackGrid( initSnowPack );
//  }  

  // Determine gridded snowpack characteristics at beginning
  //   of month 

  is = is9;
  ism1 = is - 1;
  smass = smass9;
 
  // Determine snow accumulation or decline rate 
  //   (i.e., snow[]) and density of snow for each time step

  // ksnow_dens = 250.0;
  interpolateSnow( mmax, snow_dens );



  if( 0 == kswitch )
  {
    // Initialize temperatures in soil profile

//    initSoilTemp( cmnt, rflog1 );

    // SET XFA VALUE FOR INITIAL PHASE BOUNDARY 

    resetXFA( tf, hlat );
  }

  resetNodeHeat();

/* --- CALCULATE AIR FREEZING & THAWING INDICES */

  setAirFreezeThawIndex();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::initSoilTemp( const int& cmnt, 
                                 ofstream& rflog1 )
{

// formerly TINITL
// Creates initial temperatures for nodes on calculation grid
// Used for nodes between ig and imax only


  integer dnode;
  double dtemp;
  double dxstr;

  integer j;
  integer jmax;

  integer n;
  integer nmax;
  integer nn;
  integer nn1;
  integer num;

  double temp;
  double tt[MAXNODES+1];

 
  // Initial temperatures from raw data
  
  double tstart[MAXSPDATA];  

  double x0;
  double xbot;
  
  // raw depth (user units)
  // NOTE: last xstart[] has value less than or equal to zero
  
  double xstart[MAXSPDATA];  
  double xx[MAXNODES+1];

  // Assumes INDEX is equal to 0 rather than 
  //   2 as previously read in from ecd file
  
//  return 0;
  
  // This part of the code is currently skipped over
  
  n = 0;

  do {
    ++n;

    xstart[n-1] = DEPTH[cmnt][n-1];
    tstart[n-1] = TEMP[cmnt][n-1];

  } while ( xstart[n-1] > -1.0 ); 

  nmax = n - 1;

  if( nmax <= 2 ) 
  { 
    // DWK added conditional statement immediately below 20040504 
    if( nmax <= 1 )
    {
      for( dnode = ig; dnode <= imax; ++dnode ) 
      {
        t[dnode] = tstart[0];
      }
    }
    else
    {  
      // Linear initial distribution
    
      if( ZERO == (xstart[1] - xstart[0]) ) 
      { 
        xstart[1] = xstart[0] + 0.01; 
      } // to protect

      dtemp = (tstart[1] - tstart[0]) / (xstart[1] - xstart[0]);
      temp = tstart[0];
      x0 = xstart[0];
 
      for( dnode = ig; dnode <= imax; ++dnode ) 
      {
        t[dnode] = temp + dtemp * (x[dnode] - x0);
      }
    }
  }
  else
  {
    nn = n - 1;
    xbot = x[imax];

    if( xstart[nn-1] < xbot ) 
    { 
      // Linear interpolation from table of real values
      
      dxstr = xstart[nn-1] - xstart[nn-2];

      if( ZERO == dxstr ) { dxstr = 0.01; } // to protect
      
      num = integer ((xbot - xstart[nn-1]) / dxstr);
     
      if( num <= 1 ) { num = 1; }
      
      dxstr = (xbot - xstart[nn - 1]) / num;
      dtemp = (tstart[nn] - tstart[nn-1]) / num;
      nn1 = nn + 1;

      nmax = num + nn1;
      
      for( n = nn1; n <= nmax; ++n ) 
      {
        xstart[n-1] = xstart[n-2] + dxstr;
        tstart[n-1] = tstart[n-2] + dtemp;
      }
    }

    jmax = imax - igm1;
    for( j = 1; j <= jmax; ++j ) 
    {
      dnode = igm1 + j;
      xx[j-1] = x[dnode];
    }

    interpolate( xstart, 
                 tstart, 
                 &nmax, 
                 xx, 
                 tt, 
                 &jmax,
                 SPACEINDX, 
                 rflog1 );


    for( j = 1; j <= jmax; ++j ) 
    {
      dnode = igm1 + j;
      t[dnode] = tt[j-1];
    }
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::initSoilThermalRegime( const double& tstart )
{

  integer dnode;
  
  ig = 9;
  
  for( dnode = 0; dnode < (MAXNODES+1); ++dnode ) 
  {
    x9[dnode] = ZERO;
    t9[dnode] = tstart;
    xfa9[dnode] = -1e10f;
    xfb9[dnode] = ZERO;
    water9[dnode] = 1.00000;
    dx9[dnode] = ZERO;
  }

  // Assume no snow cover
  
  is = ig;
  is9 = is;
  ism1 = is - 1;
  ism19 = ism1;
  smass = ZERO;
  smass9 = smass;
  
  for( dnode = 0; dnode < MAXSNODES; ++dnode ) 
  {
    weight9[dnode] = ZERO;
  }

  hlat = 334000.0;     // Initialize latent heat
  
  // Initialize freezing temperature threshold to 0 degrees C
  tf = ZERO;       

}; 

/* *************************************************************
************************************************************* */

int Soilthermal60::interpolate( double x[], 
                                double y[], 
                                integer *nmax, 
                                double xx[], 
                                double yy[], 
                                integer *nnmax,
                                const int& spaceTimeIndex,
                                ofstream& rflog1 )
{

// former INTRPL(X,Y,NMAX,XX,YY,NNMAX)
// Linear interpolation subroutine

  integer m;
  integer n;
  double slope;
  double x1;
  double x2;
  double y1;
  double y2;

  if( (SPACEINDX == spaceTimeIndex   
      && (*nmax > MAXSPDATA || *nnmax > (MAXNODES+1)))
      || (TIMEINDX == spaceTimeIndex 
      && (*nmax > MXTSDATA || *nnmax > MXTSFITTED)) )
  { 
    rflog1 << "Index in interpolate is out of range" << endl;
    rflog1 << " nmax = " << *nmax;
    rflog1 << " nnmax = " << *nnmax << endl;
    
    exit( -1 ); 
  }
  
  xx[*nnmax + 1] = 1e10f;

  m = 0;
  
  for( n = 1; n < *nmax; ++n ) 
  {
    y1 = y[n];
    y2 = y[n+1];
    x1 = x[n];
    x2 = x[n+1];
    
    if( x2 == x1 )  { x2 = x1 + 0.001; } // added

    slope = (y2 - y1) / (x2 - x1);

    while( m < *nnmax && xx[m] <= x2 )
    {
      yy[m] = y1 + slope * (xx[m] - x1);
      ++m;
    }
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::interpolateAirHeat( double fitted[MXTSFITTED],
                                       ofstream& rflog1 )
{

  // former DATA(INDEX,FITTED)
  // Creates temperature values for each time step m.
  // Allowed for 4 different functions of temperature vs time
  // based on INDEX.  We always assume INDEX = 3. 

  double elraw[MXTSDATA];

  integer m;

  integer nmax;

  double raw[MXTSDATA];


  for( m = 0; m < mmax; ++m ) 
  {
    eltim[m] = (m+1) * dtday;
  }

  elraw[0] = ZERO;
  raw[0] = airheat.beginning;
  elraw[1] = total * 0.5;
  raw[1] = airheat.middle;
  elraw[2] = total;
  raw[2] = airheat.end;
  nmax = 3;

  interpolate( elraw, 
               raw, 
               &nmax, 
               eltim, 
               fitted, 
               &mmax,
               TIMEINDX, 
               rflog1 );


  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::interpolateSnow( const integer& maxTimeSteps,
                                    const double& snowDensity ) 
{

  // formerly part of SNOFAL(SNOW,SDEN,COMP,ETAO,DENFAC,FCMELT)
  // Determines snow accumulation rate for each full time step
   
  integer m;   // Time step
  integer mhmax;  // Number of time steps in 0.5 month

  if ( 0 == maxTimeSteps ) 
  { 
    cerr << "maxTimeSteps = 0 in interpolateSnow()" << endl;
    exit( -1 );
  } 

  mhmax = maxTimeSteps / 2;

  // Interpolate snowfall for each 0.5 day time step
  
  for( m = 0; m < mhmax; ++m ) 
  {
    sden[m] = snowDensity;
    snow[m] = (snowdepth.middle - snowdepth.beginning) 
              * sden[m] * 2.0 / (double) maxTimeSteps;
  }

  for( m = mhmax; m < maxTimeSteps; ++m ) 
  {
    sden[m] = snowDensity;
    snow[m] = (snowdepth.end - snowdepth.middle) 
              * sden[m] * 2.0 / (double) maxTimeSteps;
  }

  return 0;

}

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Soilthermal60::InterpolatST( const double& dst1, 
                                    const double& dst2,
                                    const double& dst3,
                                    const double& dst4, 
                                    const double& dst5, 
                                    const double& dst6,  
                                    double x )
{

//Interpolating soil temperature to obtain every 1 cm depth soil 
//   temperatures with Lagrange method

// Numerical method of computation and FORTRAN lanuage, 
//   Xianyi Zheng, 1986

  double f;
  double p;
  int i;
  int j;
  double x0[6];
  double y0[6];

  x0[0] = 5.0;
  x0[1] = 10.0;
  x0[2] = 20.0;
  x0[3] = 50.0;
  x0[4] = 100.0;
  x0[5] = 200.0;

  y0[0] = dst1;
  y0[1] = dst2;
  y0[2] = dst3;
  y0[3] = dst4;
  y0[4] = dst5;
  y0[5] = dst6;


  f = 0.0;
  for( i = 0; i < 6; i++ )
  {
    p = 1.0;
    for( j =0 ; j < 6; j++ )
    {
      if( i != j ) 
      { 
        p = p * (1.0/x - x0[j]) / (x0[i] - x0[j]); 
      }
    }

    f = f + p * y0[i];
  }


  return f;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::resetNodeHeat( void )
{

  int dnode;

  for( dnode = 0; dnode <= imax; ++dnode ) 
  {
    ht[dnode] = ZERO;
    thalt[dnode] = -999.9;
    htold[dnode] = ZERO;
  }


  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::resetTimeStep( const integer& cmnt,
                                   const int& outmon,
                                   ofstream& rflog1 )
{
//  double dtgiv;
//  double FIRST = 0.000000;
//  final = 30.000000; 
  per = ZERO;

//  if( per <= 0.0 && final <= FIRST ) 
//  {
//    rflog1 << "Final time is before starting time (FIRST)" << endl;
//    exit( -1 );
//  }

  total = daze[outmon];
  total = 30.0;
  
//  if( per > 0.0 && final <= 0.0 ) 
//  {
//    final = FIRST + per;
//  }

//  total = final - FIRST;

     
  dtday = initDTDAY[cmnt];

  if( ZERO == dtday ) 
  { 
    rflog1 << "dtday = 0.0 in setSoilProfileGrid()" << endl;
    rflog1 << "Check ecd file!" << endl;
    
    exit( -1 ); 
  } 

  dt = dtday * SECSPERDAY;
  
  mmax = integer (total / dtday);

  if( 0 == mmax ) 
  {  
    rflog1 << "mmax = 0 in setSoilProfileGrid()" << endl;
    rflog1 << "Adjust dtday in ecd files!" << endl;
    
    exit( -1 ); 	
  } 


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::resetXFA( const double& thresholdFreezeT,  
                             const double& latentHeat )
{

  // Set xfa value for initial phase boundary
  
  int dnode;
  double g;
  double tu;

  dx[imax] = ZERO;

  for( dnode = ig; dnode <= imax; ++dnode ) 
  {
    dxx = dx[dnode];
    tu = t[dnode] - thresholdFreezeT;

    if( dnode < imax ) 
    {
      td = t[dnode+1] - thresholdFreezeT;
    }

    xfa[dnode] = -1e10f;
    xfb[dnode] = ZERO;

    if( (tu * td) < ZERO ) 
    { 
      // to protect
      if( ZERO == (tu - td) ) { tu = td + 0.01; } 

      g = dxx * tu / (tu - td);
      
      if( dabs( latentHeat ) >= ZEROTOL ) 
      {
        xfa[dnode] = g;
      }
    }
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::resetYr( void )
{
  // Reset annual summary variables to zero 

  yrtsoil = ZERO;
  yrdst0 = ZERO;
  yrdst5 = ZERO;
  yrdst10 = ZERO;
  yrdst20 = ZERO;
  yrdst50 = ZERO;
  yrdst100 = ZERO;
  yrdst200 = ZERO;
  yrfrontd = ZERO;
  yrfrontd2 = ZERO;
  yrthawbegin1 = ZERO;
  yrthawbegin2 = ZERO;
  yrthawend1 = ZERO;
  yrthawend2 = ZERO;
  	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::setAirFreezeThawIndex( void )
{

  int dts;

  frza = ZERO;
  thwa = ZERO;
  tairlo = 1e3f;
  tairhi = -1e3f;

  for( dts = 0; dts < mmax; ++dts ) 
  {
    tstair = taer[dts];
    
    // Determine maximum air temperature during a time step for a 
    //   particular month
    if( tstair >= tairhi ) 
    {
      tairhi = tstair;
    }
      
    // Determine minimum air temperature during a time step for a 
    //   particular month
    if( tstair <= tairlo ) 
    {
      tairlo = tstair;
    }

    // Determine cumulative number of thawing or freezing days
    //   during a particular month
    if( tstair >= ZERO ) { thwa += tstair; }
    else { frza -= tstair; }
  }

  thwa *= dtday;  // Air thawing index
  frza *= dtday;  // Air freezing index


  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Soilthermal60::setFullSnowNode( const double& timeStepSecs )
{

  // formerly part of SNOFAL(SNOW,SDEN,COMP,ETAO,DENFAC,FCMELT)
  // Determines the maximum mass (kg m^-2) of snow that a "full" 
  //   snow node can hold
  
  // Abel's constant see Table 1.1 and Equation 6.2 in 
  //   Goodrich (1976)
//  const double ABEL = 2.9e-6f;

  // Maximum density of snow in a full snow node (kg m^-3)
//  const double DENMAX = 1400.0;   // used with Abel equation
  
  // DENMAX for Sturm et al. (1997) equation that provides the 
  //   same snow thermal conductivity as the Abel's equation
  //   with a DENMAX = 1400.0  
  const double DENMAX = 908.2;
   
  // (Minimum?) fourier number for snow?
  const double EPSSNO = 0.40;
  
  // Volumetric heat capacity of snow used with Abel's 
  //   conductivity
  const double SPHSNO = 2000.0; 

  // Volumetric heat capacity of snow used with Sturm et al.'s 
  //   conductivity
//  const double SPHSNO = 2093.0; 
  
  // Thermal conductivity of snow
  double cnd;  

  // Maximum mass of snow within a "full" snow node
  double compactSnow;  
  
    
  // from Equation 6.2 in Goodrich (1976)
 //  cnd = ABEL * DENMAX * DENMAX;

  // Replace Abel's equation with Sturm et al. (1997)
  
  cnd = pow( 10, (0.002650 * DENMAX - 1.652) );

  compactSnow = sqrt( cnd * DENMAX * timeStepSecs 
                      * EPSSNO / SPHSNO );
  
  if( EPSSNO < 0.1 ) { compactSnow = 1e10f; }

  return compactSnow;

}

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::setHeatSources( integer *mmax, 
                                   integer kint[2], 
                                   double heat[MXTSFITTED][2], 
                                   double source[MXTSFITTED], 
                                   ofstream& rflog1 )
{

  // formerly SETHET(MAXM,KINT,HEAT,SOURCE)
  // Sets internal point heat source functions
  // Sets kint = 0 if source is inactive
  // Individual sources can be turned off independently: 
  //   kint[iheat] = 0

  integer dnode;
  integer iheat;
//  integer index;

  integer KINTI;
  
  integer m;
 
  double sh_dephet;
  double sh_vdep;


  for( iheat = 0; iheat < 2; ++iheat ) 
  {
    kint[iheat] = 1;
  }

  iheat = 0;
//  ++iheat;

  sh_vdep = ZERO;
  sh_dephet = ZERO; // probably problem here ??

  KINTI = 99999;
  if( 99999 == KINTI ) { return 0; }

  sh_dephet = sh_vdep * sh_dephet + x[ig];
  kint[iheat] = KINTI;

  // Locate corresponding element and set all other
  //   xhet[] to 999.9
  
  for( dnode = ig; dnode <= imax1; ++dnode ) 
  {
    xhet[dnode] = 999.9;
    if( x[dnode] <= sh_dephet && sh_dephet <= x[dnode+1] ) 
    {
      // Heat source number "iheat" located within element
      //   dnode
      
      xhet[dnode] = sh_dephet;
    }
  }

  // Determine source strength and store in heat[][]
  
  interpolateAirHeat( source, rflog1 );

  for( m = 0; m < *mmax; ++m ) 
  {
    heat[m][iheat] = source[m];
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::setInitSnowPackGrid( const double& initSnowPack )
{
  integer dnode;  // node index
  
  // Assume air temperature is equal to -1.0 if snow pack exists
  //   for updateSnowPackGrid()
  double initTair = -1.0;

  // Initialize fcmelt
  fcmelt = 5.0;
 
  // Determine the maximum mass (kg m^-2) of snow within 
  //   a "full" snow node
    
  comp = setFullSnowNode( dt );

  // Build a gridded snowpack based on snowdepth at the
  //   beginning of the month (i.e., set snow pack surface to 
  //   ground surface)
  is = ig;
  ism1 = is - 1;
  smass = ZERO;
  
  // Initialize snow characteristics in all nodes assuming
  //   no snow pack exists
     
  for( dnode = 0; dnode <= igm1; ++dnode ) 
  {
    xfa[dnode] = -1e10f;
    xfb[dnode] = ZERO;
    water[dnode] = 1.0;
    dx[dnode] = ZERO;
    weight[dnode] = ZERO;
  }


  if( initSnowPack <= 0.000001 )
  {
    // No snowpack exists this month       

    is9 = is;
    ism19 = ism1;
    smass9 = smass;
            
    return 0;
  }
  else
  {
    // Create snow pack grid based on initSnowPack and 
    //   snow_dens (NOTE: initTair = -1.0 keeps the snow 
    //   pack from melting during this initial call to 
    //   updateSnowPackGrid()

	//if (is <=1) cout<<"check2: "<<"is: " <<is<<endl;
    updateSnowPackGrid( initSnowPack, 
                        &snow_dens, 
                        weight, 
                        &comp,
                        &fcmelt, 
                        &initTair );
    //if (is <=1) cout<<"check2: "<<"is: " <<is<<endl;
    
    is9 = is;
                          
    // Initialize index of node just above snow pack surface

    ism1 = is - 1;
    ism19 = ism1;

    smass9 = smass;
  }

  return 0;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::setLowerBoundaryConditions( double *tber )
{

  // formerly BCBOT(MAXM,TBER)
  // Determines type of boundary condition at level x[imax-1]
  // Constant heat flux (positive downward - KBOT = 2)
  
  integer m;

  hbot = ZERO;
  for( m = 0; m < mmax; ++m ) { tber[m] = ZERO; }

  return 0;

};

/* *************************************************************
************************************************************* */


/* ************************************************************* 
************************************************************* */

int Soilthermal60::setMonthlySoilConditions( const double& prevtair,
                                             const double& tair, 
                                             const double& nexttair, 
                                             const double& prevspack, 
                                             const double& snowpack, 
                                             const double& nextspack, 
                                             const int& cmnt,
                                             const int& outmon,
                                             ofstream& rflog1 )
{

  // modified from the former SOILTEMP()

  integer dnode; // node index

  // Index of the top of the node that would exist two nodes
  //   below the bottom node of the soil profile 
  integer imaxp;
  
  integer m;  // time step index
  
  double r1;


  // Initialize data each month for soil temperature calculations

  initMonthlySoilConditions( outmon, 
                             prevtair, 
                             tair, 
                             nexttair,
                             prevspack, 
                             snowpack, 
                             nextspack,
                             cmnt, 
                             rflog1 );


  imaxp = imax + 1;

  for( dnode = ig; dnode <= imax; ++dnode )
  {
    tanee[dnode] = ZERO;
  }

  /* ************* BEGIN TIME STEPPING ********************** */

  time = ZERO;

  if( 0 == mmax ) 
  { 
    rflog1 << " mmax = 0 in setMonthlySoilConditions()" << endl;
    
    exit( -1 ); 
  } 

  /* --- TIME = ELAPSED TIME IN DAYS AT END OF MONTH */
  //cout <<"updateTimeStep runs okay1: "<<endl;
  for( m = 0; m < mmax; ++m )
  {
    time += dtday;

    updateTimeStep( m, cmnt, r1, imaxp, rflog1 );


    for( dnode = ig; dnode <= imax; ++dnode )
    {
      tanee[dnode] += t[dnode]/ mmax;
    }
  }
  //cout <<"updateTimeStep runs okay2: "<<endl;
  setMonthlySoilTemp();

  updateActiveLayer(); 

  return 0;

}; 

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::setMonthlySoilTemp( void )
{
  int dnode;

  int dst5flag = 0;
  int dst10flag = 0;
  int dst20flag = 0;
  int dst50flag = 0;
  int dst100flag = 0;
  int dst200flag = 0;
  int dst300flag = 0;

  int tsoilflag = 0;

 
  // define the top 20cm is from the top of the moss layer

  tsoil = ZERO;

  // Surface
  
  dst0 = tanee[9];

    
  for( dnode = 10; dnode < MAXNODES; ++dnode ) 
  {
    // Collect summary data for temperatures in
    //   the top 20 cm
    
    if( x[dnode] < 0.20 && 0 == tsoilflag )
    {
      tsoil += ((tanee[dnode] + tanee[dnode-1]) * 0.5)
               * ((x[dnode] - x[dnode-1]) / 0.20);

    }
    else { tsoilflag = 1; }
    
    
    // depth 5 cm

    if( x[dnode] >= 0.05 && 0 == dst5flag )
    {
      if( 0.05 == x[dnode] )
      {
      	dst5 = tanee[dnode];
      }
      else
      {
      	dst5 = tanee[dnode-1] 
      	       + (tanee[dnode] - tanee[dnode-1])
      	       * ((0.05 - x[dnode-1]) / (x[dnode] - x[dnode-1]));
      }
      
      dst5flag = 1;
    } 	
    
  // depth 10 cm

    if( x[dnode] >= 0.10 && 0 == dst10flag )
    {
      if( 0.10 == x[dnode] )
      {
      	dst10 = tanee[dnode];
      }
      else
      {
      	dst10 = tanee[dnode-1] 
      	        + (tanee[dnode] - tanee[dnode-1])
      	        * ((0.10 - x[dnode-1]) / (x[dnode] - x[dnode-1]));
      }
      
      dst10flag = 1;
    } 	

  // depth 20 cm

    if( x[dnode] >= 0.20 && 0 == dst20flag )
    {
      if( 0.20 == x[dnode] )
      {
      	dst20 = tanee[dnode];

        tsoil += ((tanee[dnode] + tanee[dnode-1]) * 0.5)
                 * ((x[dnode] - x[dnode-1]) / 0.20);

      }
      else
      {
      	dst20 = tanee[dnode-1] 
      	        + (tanee[dnode] - tanee[dnode-1])
      	        * ((0.20 - x[dnode-1]) / (x[dnode] - x[dnode-1]));

        tsoil += ((tanee[dnode] + tanee[dnode-1]) * 0.5)
                 * ((0.20 - x[dnode-1]) / 0.20);

      }
      
      dst20flag = 1;
    } 	

  // depth 50 cm

    if( x[dnode] >= 0.50 && 0 == dst50flag )
    {
      if( 0.50 == x[dnode] )
      {
      	dst50 = tanee[dnode];
      }
      else
      {
      	dst50 = tanee[dnode-1] 
      	        + (tanee[dnode] - tanee[dnode-1])
      	        * ((0.50 - x[dnode-1]) / (x[dnode] - x[dnode-1]));
      }
      
      dst50flag = 1;
    } 	

    // depth 100 cm

    if( x[dnode] >= 1.00 && 0 == dst100flag )
    {
      if( 1.00 == x[dnode] )
      {
      	dst100 = tanee[dnode];
      }
      else
      {
      	dst100 = tanee[dnode-1] 
      	        + (tanee[dnode] - tanee[dnode-1])
      	        * ((1.00 - x[dnode-1]) / (x[dnode] - x[dnode-1]));
      }
      
      dst100flag = 1;
    } 	

    // depth 200 cm

    if( x[dnode] >= 2.00 && 0 == dst200flag )
    {
      if( 2.00 == x[dnode] )
      {
      	dst200 = tanee[dnode];
      }
      else
      {
      	dst200 = tanee[dnode-1] 
      	        + (tanee[dnode] - tanee[dnode-1])
      	        * ((2.00 - x[dnode-1]) / (x[dnode] - x[dnode-1]));
      }
      
      dst200flag = 1;
    }
    if( x[dnode] >= 3.00 && 0 == dst300flag )
    {
      if( 3.00 == x[dnode] )
      {
      	dst300 = tanee[dnode];
      }
      else
      {
      	dst300 = tanee[dnode-1]
      	        + (tanee[dnode] - tanee[dnode-1])
      	        * ((3.00 - x[dnode-1]) / (x[dnode] - x[dnode-1]));
      }

      dst300flag = 1;
    }


  }
  
     
  // Save snow and soil temperature variables for starting 
  //   conditions for the next month

  is9 = is;
  ism19 = ism1;
  smass9 = smass;
  for( dnode = 0; dnode < MAXNODES; ++dnode ) 
  {
    x9[dnode] = x[dnode];
    dx9[dnode] = dx[dnode];
    xfa9[dnode] = xfa[dnode];
    xfb9[dnode] = xfb[dnode];
    water9[dnode] = water[dnode];
    t9[dnode] = t[dnode];
  }

  for( dnode = 0; dnode < MAXSNODES; ++dnode ) 
  {
    weight9[dnode] = weight[dnode];
  }
  
  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::setSnowCharSturm( double *ta, 
                                     double *sph, 
                                     double *cnd )
{ 

// Snow conductivity - Sturm et al. (1997)  formula

  if ( *ta >= 0.0 ) 
  {
    // Thawed
    
    *sph = 4086.0;
    *cnd = 10.0;
  }
  else
  {
    // Frozen
    
    *sph = 2093.0;
    *cnd = pow( 10, (0.002650 * snow_dens - 1.652) );
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Soilthermal60::setSnowDensity( const int& cmnt, 
                                      const double& prec, 
                                      const int& outmon )
{

  // snow classification determine the snow density, and therefore 
  // the snow depth, and snow thermal conductivity
  // NOTE: CDM is cooling degree month
  
  /* ***********************************************************
   Snow classes following Sturm et al. (1995):
   snow_class = 1: Lowland tundra 
   snow_class = 2: Alpine tundra
   snow_class = 3: Taiga forest
   snow_class = 4: Maritime forest
   snow_class = 5: Prairie grasslands
   snow_class = 6: Ephemeral snow cover
  *********************************************************** */

  int snow_class;
  
  // wind_speed is a relative index:
  // 0 = low  ( <= 0.5 m s^-1)
  // 1 = high ( >= 2 m s^-1)
  int wind_speed; 

  // Use prescribed wind speed levels

  wind_speed = vegWindSp[cmnt];
  
  // Following classification based on Table 7 in Sturm et 
  //   al. (1995)
  
  if( CDM < 50.0 ) { snow_class = 6; } 
  else
  {
    if ( CDM > 125.0 ) 
    {
      if ( (prec / daze[outmon]) > 2.0 )
      {
        if ( wind_speed == 1 ) { snow_class = 1; }
        else { snow_class = 2; }
      }
      else
      {
        if ( wind_speed == 0 ) { snow_class = 2; }
        else { snow_class = 1; }
      }
    }
    else 
    {
      if ( (prec / daze[outmon]) < 2.0 )
      {
        if ( wind_speed == 1 ) { snow_class = 4; }
        else { snow_class = 3; }
      }
      else
      {
        if ( wind_speed == 0 ) { snow_class = 5; }
        else { snow_class = 5; }
      }
    }
  }

//  Determine snow density (snow_dens) in kg m^-3 based on
//    Table 1 in Liston and Pielke (2000)

  switch ( snow_class ) 
  {
    case 1: snow_dens = 280.0;  break;
    case 2: snow_dens = 225.0;  break;
    case 3: snow_dens = 250.0;  break;
    case 4: snow_dens = 250.0;  break;
    case 5: snow_dens = 300.0;  break;
    case 6: snow_dens = 350.0;  break;
  }

  return snow_dens;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::setSoilProfileGrid( const int& cmnt,
                                       ofstream& rflog1 )
{

  // formerly GRID(MAXM)
  // Sets grid spacing, node depths, material type, density, water 
  //   content and explicit thermal properties if mat <= 10
  
  double cnd;

  doubledouble d__1;
  doubledouble d__2;

  double dense;
  double df;
  integer dnode;
  double dtbig;  // biggest usable time step (days)
  double dtmax;
  double dxa;
  double dxagiv;
  double dxb;
  double dxmin;
  double dxx;

  double eps;
  double err;

  double f;
  double fwas;

  double grdepsmin;
  double grdtop;
  double grdvcond;
//  double grdvdens;
//  double grdvspace;
  double grdvsph;

  integer j;
  integer k;

  integer layer;
  integer lmimax;

  integer mat;       // material index

  integer n;
  integer nn;
  integer nodes;

  double power;

  double r;

  double ratio;
  double rnew;
  double rwas;

  double sph;
  double sph_dense;

  double temp;
  double thick;
  double thmdxa;
//  double total;

//  double usedxa;

  double wat;
 
  double xx;

  lmimax = MAXNODES - 2;

  theta = 0.50;     // Initialize theta
  
  ig = 9;          // Initialize index for ground surface node
  
  grdtop = 0.000000; // Initialize grdtop
  
  grdepsmin = 0.10;  // Initialize epsmin

  if( theta <= 0.5 ) { theta = 0.5; }
  
  theta *= 2.0;
  theta1 = 2.0 - theta;
  igm1 = ig - 1;
  xx = grdtop;
  x[ig] = xx;
  dnode = igm1;
  layer = 0;
  dtbig = 1e10f;

  do {
    // Assign initial soil layer characteristics  	
  	
    ++layer;

    thick = initTHICK[cmnt][layer-1]; 
    dxa = initDXA[cmnt][layer-1]; 
    dxb = initDXB[cmnt][layer-1];
    mat = vegMAT[cmnt][layer-1];
    dense = initDENSE[cmnt][layer-1]; 
    wat = initWATER[cmnt][layer-1];

//    if( kswitch != 0 )
//    {
//      // adding for using soil water from TEM, 29/March/2001
//      if ( layer == 3 )
//      {
//        wat = soilwater.layer2;
//      }
//      else { wat = initWATER[cmnt][layer-1]; }
 
//      wat = initWATER[cmnt][layer-1];
//    }

    if( thick <= ZERO ) { break; }

    if( thick < dxa ) 
    { 
      rflog1 << "thick is less than dxa" << endl;
      exit( -1 ); 
    }

    if( mat < 4 ) 
    { 
      // Use constant thermal properties if mat < 4
      
      grdvcond = vcond[cmnt][layer-1]; 
      grdvsph = vsph[cmnt][layer-1]; 
      condt[mat] = initCONDT[cmnt][layer-1];
      spht[mat] = initSPHT[cmnt][layer-1]; 
      condf[mat] = initCONDF[cmnt][layer-1];
      sphf[mat] = initSPHF[cmnt][layer-1];

      if( ZERO == grdvcond ) { grdvcond = 1.0; }
      
      if( ZERO == grdvsph ) { grdvsph = 1.0; }
      
      condt[mat] *= grdvcond;
      spht[mat] *= grdvsph;
      condf[mat] *= grdvcond;
      sphf[mat] *= grdvsph; 

      cnd = condf[mat];
      sph = sphf[mat];
    }
    else
    {
      if( mat > MAXMAT ) 
      { 
        rflog1 << "Material index (mat) = " << mat; 
        rflog1 << " is undefined" << endl;	
        
        exit( -1 ); 
      }
        
      temp = tf + 1.0;
      for( n = 0; n < 10; ++n ) 
      {
        updateSoilThermalProperties( &mat, 
                                     &temp,   
                                     &sph, 
                                     &cnd, 
                                     rflog1 );

        temp += -1.0;
      }
    }

    /* --- CALCULATE MINIMUM ALLOWED GRID SPACING */

    if( ZERO == (sph * dense) )  {  sph_dense = 0.01; }
    else  { sph_dense = sph * dense; }

//    eps = epsmin * 86400.0 * cnd / (sph * dense);

    eps = grdepsmin * SECSPERDAY * cnd / sph_dense;   // to protect
    dxagiv = dxa;

    do {
      // to protect
      if ( (eps * dtday) < ZERO ) { eps = 0.01; }; 

      dxmin = sqrt( eps * dtday );

      if( dxa < dxmin ) 
      { 
        dxa = dxmin;

        if( dxa >= thick ) 
        { 
          dxa = dxagiv;
          dtday *= 0.5;
          mmax *= 2;
        }
      }
    } while( (dxa < dxmin) && (dxa >= thick) );  

    if( ZERO == eps ) { eps = 0.01;} // to protect

    dtmax = dxa * dxa / eps;
    if( dtmax < dtbig ) { dtbig = dtmax; }

    /* --- SET VALUES AT CALCULATION NODES FOR LAYER N */

    if( dxa >= thick * 0.5 || dxb <= dxa ) 
    {   
      /* --- LINEAR GRID SPACING WITHIN LAYER */

      if( ZERO == dxa ) { dxa = 0.01; } //to protect
    
      nodes = integer (thick / dxa + 0.01);
      if( 0 == nodes ) { nodes = 1; }

      dxa = thick / (double) nodes;
//      usedxa = dxa / grdvspace;

      for( j = 0; j < nodes; ++j ) 
      {
        ++dnode;

        if( dnode == lmimax ) { break; }
     
        x[dnode] = xx;
        dx[dnode] = dxa;
        xx += dxa;
        mater[dnode] = mat;
        ddry[dnode] = dense;
        water[dnode] = wat;
      }
    }
    else
    {
      if( ZERO == (thick - dxb) ) { thick = dxb + 0.01; } // to protect
      
      r = (thick - dxa) / (thick - dxb);

      if( ZERO == dxa ) { dxa = 0.01; } //to protect
      
      ratio = dxb / dxa;

      if( 1.0 == r ) { r = 0.01; }  //to protect
      
      if( r < ZERO ) { r = 0.01; }
      
      if( ZERO == dxa ) { dxa = 0.01; } //to protect
      
      nn = integer (log( fabs( dxb / dxa ) ) * 10.0 / log( r ) + 5.0);
      nn /= 10;
      nodes = nn + 1;

      /* --- REVISE VALUE OF R BY SECANT METHOD */

      thmdxa = thick - dxa;
      power = 1.0 / nn;
      d__1 = (doubledouble) ratio;
      d__2 = (doubledouble) power;
      rwas = pow( d__1,d__2 );

      fwas = dxa * pow( rwas, (double) nodes ) - rwas * thick + thmdxa;
      err = 1e-6f;
      for( k = 0; k < 20; ++k ) 
      {
        f = dxa * pow( r, (double) nodes ) - r * thick + thmdxa;

        if( fabs( f ) <= err ) { break; }
        
        df = f - fwas;

        if( dabs( df ) < ZEROTOL ) { break; }
        
        rnew = r - (r - rwas) * f / df;
        rwas = r;
        fwas = f;
        r = rnew;
      }

      r = rnew;
      dxx = dxa;
      for( j = 0; j < nodes; ++j ) 
      {
        ++dnode;
        if( dnode == lmimax ) { break; }  

        x[dnode] = xx;
        dx[dnode] = dxx;
        xx += dxx;
        dxx *= r;
        mater[dnode] = mat;
        ddry[dnode] = dense;
        water[dnode] = wat;
      }
    }
  } while( dnode <= lmimax );

  imax1 = dnode;
  imax = dnode + 1;
  dt = dtday * SECSPERDAY;
  
  return 0;

}; 

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::setUpperBoundaryConditions( double taer[MXTSFITTED],
                                               ofstream& rflog1 )
{

  // formerly BCTOP(MAXM,TAER,KTOP,SFLX1,SFLX2,SFLX3)
  // Determines type of boundary condition at top surface, level 
  //   x[is]

  double fac;
  double fznfac;   // freezing n-factor

//  integer index;

  integer m;

  double thnfac;    // thawing n-factor
  double tint;

  htop = 1e20f;
  sflux = ZERO;

  interpolateAirHeat( taer, rflog1 );

  thnfac = 1.0;
  fznfac = 1.0;

  for( m = 0; m < mmax; ++m ) 
  {
    tint = taer[m];
    fac = thnfac;
    if ( tint <= ZERO ) { fac = fznfac; }
    taer[m] = fac * tint;
  }
  
  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::showsnowecd( int& pdcmnt ) 
{

 // clrscr();
  cout << endl;
  cout << "             SNOW PARAMETERS FOR VEGETATION TYPE STANDS ";
  cout << endl << endl << endl << endl;
                 
  printf( "GFLUX  = %5.1lf   CDSNOW  = %5.1lf DTDAY  = %5.1lf\n", 
          gflux[pdcmnt],           
          cdsnow[pdcmnt],   
          initDTDAY[pdcmnt] );
                    
  cout << endl << "Press any key to continue . . ." << endl;
 // while (getch() == '\0');

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::showsoillecd( int& pdcmnt ) 
{
  // clrscr();
  cout << endl;
  cout << "         PARAMETERS FOR SOIL LAYERS " << endl << endl;
  printf( "THICK1 =%5.1f  DXA1 =%5.1f DXB1 =%5.1f MAT1 =%d    DENSE1=%5.1f WATER1= %5.3f\n",
          initTHICK[pdcmnt][0], 
          initDXA[pdcmnt][0], 
          initDXB[pdcmnt][0], 
          vegMAT[pdcmnt][0],
          initDENSE[pdcmnt][0], 
          initWATER[pdcmnt][0] );
          
  printf( "VCOND1 =%5.1f  VSPH1=%5.1f CONDT1=%5.1f SPHT1=%5.1f CONDF1=%5.1f SPHF1 = %5.1f\n",
          vcond[pdcmnt][0], 
          vsph[pdcmnt][0], 
          initCONDT[pdcmnt][0], 
          initSPHT[pdcmnt][0], 
          initCONDF[pdcmnt][0], 
          initSPHF[pdcmnt][0] );
          
  printf( "THICK2 =%5.1f  DXA2 =%5.1f DXB2 =%5.1f MAT2 =%d    DENSE2=%5.1f WATER2= %5.3f\n",
          initTHICK[pdcmnt][1], 
          initDXA[pdcmnt][1], 
          initDXB[pdcmnt][1], 
          vegMAT[pdcmnt][1],
          initDENSE[pdcmnt][1], 
          initWATER[pdcmnt][1] );
          
  printf( "VCOND2 =%5.1f  VSPH2=%5.1f CONDT2=%5.1f SPHT2=%5.1f CONDF2=%5.1f SPHF2 = %5.1f\n",
          vcond[pdcmnt][1], 
          vsph[pdcmnt][1], 
          initCONDT[pdcmnt][1], 
          initSPHT[pdcmnt][1],
          initCONDF[pdcmnt][1], 
          initSPHF[pdcmnt][1] );
          
  printf( "THICK3 =%5.1f  DXA3 =%5.1f DXB3 =%5.1f MAT3 =%d    DENSE3=%5.1f WATER3= %5.3f\n",
          initTHICK[pdcmnt][2], 
          initDXA[pdcmnt][2], 
          initDXB[pdcmnt][2], 
          vegMAT[pdcmnt][2],
          initDENSE[pdcmnt][2], 
          initWATER[pdcmnt][2]);
          
  printf( "VCOND3 =%5.1f  VSPH3=%5.1f CONDT3=%5.1f SPHT3=%5.1f CONDF3=%5.1f SPHF3 = %5.1f\n",
          vcond[pdcmnt][2], 
          vsph[pdcmnt][2], 
          initCONDT[pdcmnt][2], 
          initSPHT[pdcmnt][2], 
          initCONDF[pdcmnt][2], 
          initSPHF[pdcmnt][2] );
          
  printf( "THICK4 =%5.1f  DXA4 =%5.1f DXB4 =%5.1f MAT4 =%d    DENSE4=%5.1f WATER4= %5.3f\n",
          initTHICK[pdcmnt][3], 
          initDXA[pdcmnt][3], 
          initDXB[pdcmnt][3], 
          vegMAT[pdcmnt][3], 
          initDENSE[pdcmnt][3], 
          initWATER[pdcmnt][3]);
          
  printf( "VCOND4 =%5.1f  VSPH4=%5.1f CONDT4=%5.1f SPHT4=%5.1f CONDF4=%5.1f SPHF4 = %5.1f\n",
          vcond[pdcmnt][3], 
          vsph[pdcmnt][3], 
          initCONDT[pdcmnt][3], 
          initSPHT[pdcmnt][3], 
          initCONDF[pdcmnt][3], 
          initSPHF[pdcmnt][3] );
          
  printf( "THICK5 =%5.1f  DXA5 =%5.1f DXB5 =%5.1f MAT5 =%d    DENSE5=%5.1f WATER5= %5.3f\n",
          initTHICK[pdcmnt][4], 
          initDXA[pdcmnt][4], 
          initDXB[pdcmnt][4], 
          vegMAT[pdcmnt][4],
          initDENSE[pdcmnt][4], 
          initWATER[pdcmnt][4] );
          
  printf( "VCOND5 =%5.1f  VSPH5=%5.1f CONDT5=%5.1f SPHT5=%5.1f CONDF5=%5.1f SPHF5 = %5.1f\n",
          vcond[pdcmnt][4], 
          vsph[pdcmnt][4], 
          initCONDT[pdcmnt][4], 
          initSPHT[pdcmnt][4], 
          initCONDF[pdcmnt][4], 
          initSPHF[pdcmnt][4] );
          
  printf( "THICK6 =%5.1f  DXA6 =%5.1f DXB6 =%5.1f MAT6 =%d    DENSE6=%5.1f WATER6= %5.3f\n",
          initTHICK[pdcmnt][5], 
          initDXA[pdcmnt][5], 
          initDXB[pdcmnt][5], 
          vegMAT[pdcmnt][5],
          initDENSE[pdcmnt][5], 
          initWATER[pdcmnt][5]);
          
  cout << endl << "Press any key to continue . . ." << endl;
//  while (getch() == '\0');

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::showsoiltecd( int& pdcmnt ) 
{
  int j;
//  clrscr();
  cout << "INITIALIZATION FOR SOIL TEMPERATURES:" << endl;
  printf( "DEPTH     TEMP \n" );
  
  for( j = 0; j < 15; ++j )
  { 
    cout << DEPTH[pdcmnt][j] << "    " << TEMP[pdcmnt][j] << endl;
  }
  
  cout << endl << "Press any key to continue . . ." << endl;
//  while ( getch() == '\0' );

  printf( "DEPTH     TEMP \n" );
  
  for ( j = 15; j < 25; ++j )
  { 
    cout << DEPTH[pdcmnt][j] << "    " << TEMP[pdcmnt][j]<< endl;
  }
  cout << endl << "Press any key to continue . . ." << endl;
//  while (getch() == '\0');


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::updateActiveLayer( void )
{

/* *************************************************************
Determine the active layer, including the depth of the frozen 
  front from the ground surface (frontd) and the depth to the 
  frozen front of permafrost (deepd) if it exists
************************************************************* */

  integer dnode; // node index
  
  // Index of node containing the first phase change plane
  //   (must be above node 37)
  integer first_phase;

  // Index of the node just above the node containing the second
  //   phase change plane (must be above node 37)  
  integer second_phase;

  // Index of the node just above the node containing the third 
  //   phase change plane (must be above node 37)  
  integer third_phase;

  // Index of the node just above the node containing the fourth 
  //   phase change plane (must be above node 37)  
  integer fourth_phase;  
  
  
  frontd = ZERO; 
  first_phase = 8;
  for( dnode = 9; dnode < 37; ++dnode )
  {
    if( tanee[dnode] < tf )
    {
      first_phase = dnode;
    }
    else { dnode = 37; } // break from loop
  }


  // Determine the depth of the freezing front from soil surface
  
  if( first_phase > 8 && first_phase < 36 )
  {
//      frontd = (x[first_phase] 
//               + ((x[first_phase+1] - x[first_phase]) / 2.0));

    frontd = (x[first_phase] 
             + ((x[first_phase+1] - x[first_phase])
             *((tf - tanee[first_phase])
             /(tanee[first_phase+1] - tanee[first_phase]))));

    if( frontd < 0.05 ) { frontd = ZERO; }
  }
  else 
  { 
    if( first_phase <= 8 ) { frontd = ZERO; }
    else { frontd = x[first_phase]; }
  } 

  thawbegin1 = frontd;
  
  second_phase = first_phase;
  if( first_phase < 36 )
  {
    for( dnode = (first_phase+1); dnode < 37; ++dnode )
    {
      if( tanee[dnode] > tf )
      {
        second_phase = dnode;
      }
      else { dnode = 37; } // break from loop
    }
  }
  
  
  // Determine top of permafrost at bottom of soil profile
  
  if( second_phase < 36 )
  {
//      thawend = (x[second_phase] 
//                + ((x[second_phase+1] - x[second_phase]) / 2.0 ));

      thawend1 = (x[second_phase] 
                 + ((x[second_phase+1] - x[second_phase])
                 * (tanee[second_phase] - tf)
                 /(tanee[second_phase] - tanee[second_phase+1])));
  }
  else { thawend1 = x[36]; } 

  frontd2 = thawend1;

  third_phase = second_phase;
  if( third_phase < 36 )
  {
    for( dnode = (second_phase+1); dnode < 37; ++dnode )
    {
      if( tanee[dnode] < tf )
      {
        third_phase = dnode;
      }
      else { dnode = 37; } // break from loop
    }
 
    if( third_phase < 36 )
    {
      frontd2 = (x[third_phase] 
                + ((x[third_phase+1] - x[third_phase])
                *((tf - tanee[third_phase])
                /(tanee[third_phase+1] - tanee[third_phase]))));    
    }
    else { frontd2 = x[third_phase]; }  
  }
  
  thawbegin2 = frontd2;

  fourth_phase = third_phase;
  if( fourth_phase < 36 )
  {
    for( dnode = (third_phase+1); dnode < 37; ++dnode )
    {
      if( tanee[dnode] > tf )
      {
        fourth_phase = dnode;
      }
      else { dnode = 37; } // break from loop
    }
  }
  
  
  // Determine top of permafrost at bottom of soil profile
  
  if( fourth_phase < 36 )
  {
      thawend2 = (x[fourth_phase] 
                 + ((x[fourth_phase+1] - x[fourth_phase])
                 * (tanee[fourth_phase] - tf)
                 /(tanee[fourth_phase] - tanee[fourth_phase+1])));
  }
  else { thawend2 = x[36]; } 


  activeLayer = (thawend1 - thawbegin1) 
                + (thawend2 - thawbegin2);
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::updateMultiplePhasePlaneTemps( integer *itzero, 
                                                  integer *ibzero,
                                                  const integer& cmnt,
                                                  ofstream& rflog1 )
{

// formerly ASMTWO(ITZERO,IBZERO)

// Calculates new positions of more than one phase plane
// Calculates new temperatures at all nodes
// Phase planes are uncoupled
// Temperature = 0.0 between phase planes
// Only exterior phase planes are mobile

  double cd;
  double cnd;
  integer counter;
  double cu;

  double dens;
  double dtwas;
  double dxx;

  double ga;
  double gaold;

  double gb;
  double gbold;

  double gdn;

  double gold;
  double gup;
  double gv;

  integer dnode;
  integer im1;
  integer ip1;
  integer ip2;

  double r1;
  double rlw;
  double rsd;
  double rsu;

  double sph;

  double t1;
  double t2;

  double td;
  double tt;

  double v;
  double vdn;
  double vup;

  double wat;

  double zero1 = ZERO;

  gaold = xfa[*itzero];
  gbold = xfa[*ibzero] + xfb[*ibzero];
  
/* --- SECTION FOR UPPER PHASE PLANE */
  
  dnode = *itzero;
  dtwas = dtfaz;
  v = xfb[dnode];
  im1 = dnode - 1;
  ip1 = dnode + 1;
  ip2 = dnode + 2;

  counter = 0;
  while( counter < THOUSANDX100 ) 
  {
    gaussianEliminForward( &im1 );

    dxx = dx[dnode];
    wat = water[dnode];
    dens = ddry[dnode];
    tt = t[dnode];
    td = t[ip1];

    updateSoilThermalProperties( &mater[dnode],
                                 &tt,
                                 &sph,
                                 &cnd, 
                                 rflog1 );
    cu = cnd;
    rsu = dens * sph;
    rlw = dens * hlat * wat;

    updatePhasePlane( &dnode, 
                      &rlw, 
                      &cu, 
                      &rsu, 
                      &zero1, 
                      &zero1, 
                      &tt, 
                      &tf,
                      &dxx, 
                      &gaold, 
                      &ga, 
                      &t1, 
                      &t2, 
                      rflog1 );

    counter++;
    
    /* --- EVALUATE WHETHER GA GOES OUT OF ELEMENT */

    if( dabs( v ) >= ZEROTOL )
    {
      // Second phase plane within element exists
      
      t[dnode] = t1;
      
      crankNicholsonBackward( &im1 );
      
      gv = v + gaold;
      xfb[dnode] = gv-ga;
      xfa[dnode] = ga;
      if ( gv > ga ) { break; }
      xfa[dnode] = -1e10f;
      xfb[dnode] = ZERO;
      if ( dnode < *ibzero ) { break; }
      xfa[dnode] = ga;
      xfb[dnode] = gv-ga;
      break;
    }
    else
    {
      if( ga <= dxx )
      {
        // ga remains within element "dnode"
        
        t[dnode] = t1;
        
        crankNicholsonBackward( &im1 );
        
        xfa[dnode] = ga;
        break;
      }      
      else
      {
        /* --- OUT BELOW */

        if( dxx < 0.01 ) { dxx = 0.01; } // to protect
        
        gdn = xfa[ip1];
        capx[dnode] = rsu * dxx;
        conx[dnode] = cu / dxx;
        xfa[dnode] = -1e10f;
        
        if( dnode == imax1 ) { return 0; }
        
        ++dnode;
        im1 = dnode - 1;
        ip1 = dnode + 1;
        ip2 = dnode + 2;
        dtfaz = dtwas - dtfaz;
        dtwas = dtfaz;
        r1 = td - tf;
        
        if( r1 >= ZERO )
        {
          t[dnode] = tf - ADJUST;
        }
        else { t[dnode] = tf + ADJUST; }

        gaussianEliminForward( &im1 );

        crankNicholsonBackward( &im1 );

        gaold = ZERO;

        if( gdn != -1e10f )
        {
          vdn = xfb[dnode];

          if( dabs( vdn ) >= ZEROTOL )
          {
            gaold = vdn + gdn;
            v = ZERO;
            xfb[dnode] = ZERO;
          }
          else { v = gdn; }
        }
      }
    } 
  }
  
  if( THOUSANDX100 == counter )
  {
    rflog1 << "First phase plane of multiple phase planes not ";
    rflog1 << "resolved (counter reached 100,000)" << endl;
    
    exit( -1 );
  }
 
  if( *ibzero < *itzero ) { return 0; }

/* --- SECTION FOR LOWER PHASE PLANE */

  dnode = *ibzero;
  dtwas = dt;
  dtfaz = dt;
  gold = xfa[dnode];
  v = xfb[dnode];
  im1 = dnode - 1;
  ip1 = dnode + 1;
  ip2 = dnode + 2;

  counter = 0;
  while( counter < THOUSANDX100 ) 
  {
    gaussianEliminBackward( &ip2, cmnt );
    dxx = dx[dnode];
    wat = water[dnode];
    dens = ddry[dnode];
    tt = t[dnode];
    td = t[ip1];
 
    updateSoilThermalProperties( &mater[dnode], 
                                 &td, 
                                 &sph, 
                                 &cnd, 
                                 rflog1 );

    cd = cnd;
    rsd = dens * sph;
    rlw = dens * hlat * wat;

    updatePhasePlane( &dnode, 
                      &rlw, 
                      &zero1, 
                      &zero1, 
                      &cd, 
                      &rsd, 
                      &tf, 
                      &td, 
                      &dxx, 
                      &gbold, 
                      &gb, 
                      &t1, 
                      &t2, 
                      rflog1 );


    /* --- EVALUATE WHETHER G GOES OUT OF ELEMENT */

    if( dabs( v ) >= ZEROTOL ) 
    { 
      // Second phase plane within element exits
      
      t[ip1] = t2;

      crankNicholsonForward( &ip2 );

      xfa[dnode] = gold;
      xfb[dnode] = gb - gold;
      
      if(gb > gold) { return 0; }

      xfa[dnode] = -1e10f;
      xfb[dnode] = ZERO;

      return 0;
    }
    else
    {
      if( gb >= ZERO )
      { 
      	// gb remains within element "dnode"
      	
        t[ip1] = t2;

        crankNicholsonForward( &ip2 );

        xfa[dnode] = gb;
        return 0;
      }
      else
      { 
        /* --- OUT ABOVE */

        gup = xfa[im1];
        capx[dnode] = rsd * dxx;

        if( dxx < 0.01 ) { dxx = 0.01; } // to protect
        
        conx[dnode] = cd / dxx;
        xfa[dnode] = -1e10f;
        
        if( dnode == is ) { return 0; }
        
        --dnode;
        im1 = dnode - 1;
        ip1 = dnode + 1;
        ip2 = dnode + 2;
        dtfaz = dtwas - dtfaz;
        dtwas = dtfaz;
        r1 = tt - tf;
        
        if( r1 >= ZERO ) 
        { 
          t[ip1] = tf - ADJUST; 
        }
        else { t[ip1] = tf + ADJUST; }

        gaussianEliminBackward( &ip2, cmnt );

        crankNicholsonForward( &ip2 );

        gbold = dx[dnode];
        v = ZERO;

        if( gup != -1e10f )
        { 
          gold = gup;
          vup = xfb[dnode];

          if( dabs( vup ) >= ZEROTOL ) 
          { 
            gbold = gup;
            v = ZERO;
            xfb[dnode] = ZERO;
          }
          else { v = gbold - gup; }
        }
      }
    }
  }

  rflog1 << "Second phase plane of multiple phase planes not ";
  rflog1 << "resolved (counter reached 100,000)" << endl;
  
  exit( -1 );  

  return 0;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

 int Soilthermal60::updatePhasePlane( integer *dnode, 
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
                                      ofstream& rflog1 )
{

// formerly PHASE(I,RLW,CU,RSU,CD,RSD,TT,TD,DXX,GOLD,G,T1,T2)
// Secant solution for phase plane or dry node equation
// Uses time averaged heat fluxes for phase element
// Formulation is ordinary implicit if flux = 0 on one side of 
//   phase plane; also, if g is closer than star to element boundary

  double a1;
  double a1p;
  double a2;
  double a2p;

  double alf1;
  double alf2;

  double b1;
  double b1p;
  double b2;
  double b2p;

  double bet1;
  double bet2;

  double c1; 
  double c2;

  double capd;
  double capu;
  double cdtf;
  double cond;
  double conu;
  double cutf;

  double d1;
  double d2;

  double df;
  double dg;
  double dgdgdg; 

  double err;

  double f;
  double fld;
  double flu;
  double ft;
  double fwas;

  double gam1;
  double gam2;
  double gnew;
  double gwas;

  double h;
  double hold;

  integer im1;
  integer ip1;
  integer ip2;

  integer n;

  double r1;
  double r2;
  double rlwt;
  double rsd4;
  double rsu4;

  double star;

  double ttt1;
  double ttt2;


  star = *dxx * 0.05;
  hold = *dxx - *gold;
  im1 = *dnode - 1;
  ip1 = *dnode + 1;
  ip2 = *dnode + 2;


  if( ZERO == dtfaz )  
  {  
    dtfaz = dtday * SECSPERDAY; // debugging by Q.Z. Nov/28/2000
  }

  conu = conx[im1];
  capu = capx[im1] * 0.5 / dtfaz;
  cond = conx[ip1];
  capd = capx[ip1] * 0.5 / dtfaz;
  rsu4 = *rsu * 0.25;
  d1 = rsu4 / dtfaz;
  c1 = conu * (1.000000 - s[im1]) + capu + d1 * *gold;
  b1 = d1 * *tt;
  a1 = (capu + d1 * *gold) * *tt + conu * 1.00001 * e[im1] 
       + ht[*dnode];

  if( *dnode == is ) 
  {
    a1 += sflux;
  }

  cutf = *cu * tf;
  a1p = a1 - tf * c1;
  b1p = b1 - tf * d1;
  rsd4 = *rsd * 0.25;
  d2 = rsd4 / dtfaz;
  c2 = cond * (1.0 - s[ip2]) + capd + d2 * hold;
  b2 = d2 * *td;
  a2 = (capd + d2 * hold) * *td + cond * e[ip2] 
       + ht[ip1];
  cdtf = *cd * tf;
  a2p = a2 - tf * c2;
  b2p = b2 - tf * d2;

  r1 = a1p * *cu;
  if( r1 < ZERO ) { *rlw = -(*rlw); } 
  
  if( dabs( r1 ) < ZEROTOL )
  { 
    r2 = a2p * *cd;
    
    if( dabs( r2 ) < ZEROTOL ) { return 0; } 
    
    if( r2 > ZERO ) { *rlw = -(*rlw); }
  } 
 
  if ( ZERO == dtfaz ) 
  {    
    dtfaz = dtday * SECSPERDAY; // debugging by Q.Z. Nov/28/2000
  }

  rlwt = *rlw / dtfaz;
  alf1 = *cu * (a1p + b1p * *gold);
  gam1 = c1 + d1 * *gold;
  bet1 = *cu + gam1 * *gold;

  alf2 = *cd * (a2p + b2p * hold);

  gam2 = c2 + d2 * hold;
  bet2 = *cd + gam2 * hold;

  // to protect
  if( ZERO == (rlwt * bet1 * bet2 + alf1 * gam2 - alf2 * gam1) )  
  {  
    dgdgdg = 0.01;
  }
  else { dgdgdg = rlwt * bet1 * bet2 + alf1 * gam2 - alf2 * gam1; }

//    dg = (alf2 * bet1 + alf1 * bet2) / (rlwt * bet1 * bet2 
//         + alf1 * gam2 - alf2 * gam1);
  dg = (alf2 * bet1 + alf1 * bet2) / dgdgdg;

  if( dg >= *dxx ) { dg = *dxx * 0.99; }
  
  if( dg <= ZERO ) { dg = *dxx * 0.10; }
  
  *g = *gold + dg * 0.50;
  h = *dxx - *g;
  ft = rlwt * (*g - *gold);

  if( ZERO == (*cu + (c1 + d1 * *g) * *g) ) 
  {   
    *cu = (c1 + d1 * *g) * *g + 0.01; // to protect Q. Z. Dec/12/2000
  }

  if( ZERO == (*cd + (c2 + d2 * h) * h) ) 
  { 
    *cd = -(c2 + d2 * h) * h + 0.01; // to protect
  }

  flu = *cu * (a1p + b1p * *g) / (*g * (c1 + d1 * *g) + *cu);
  fld = *cd * (a2p + b2p * h) / (h * (c2 + d2 * h) + *cd);

  fwas = ft - flu - fld;
  gwas = *g;
  *g = *gold + dg;

  if( dabs( *rlw ) >=  ZEROTOL ) 
  { 
    err = 1e-5f;

    for( n = 0; n < 50; ++n ) 
    {
      h = *dxx - *g;
      ft = rlwt * (*g - *gold);

    // the following modification is made on 28/NOv/2000 QZ
      ttt1 = *g * (c1 + d1 * *g) + *cu;
      ttt2 = h * (c2 + d2 * h) + *cd ;
      
      if( ZERO == ttt1 ) { ttt1 = 3.01; }
      
      if( ZERO == ttt2 ) { ttt2  = 3.001; }

//	flu = *cu * (a1p + b1p * *g) / (*g * (c1 + d1 * *g) + *cu);
//     fld = *cd * (a2p + b2p * h) / (h * (c2 + d2 * h) + *cd);

      flu = *cu * (a1p + b1p * *g) / ttt1;
      fld = *cd * (a2p + b2p * h) / ttt2;

      if( *gold > star && hold > star ) 
      { 
        if( dabs( *gold ) < ZEROTOL ) 
        {   
          //printf(" = 0.0 ");
          *gold = 0.01; 
        } // to protect

        if( ZERO == hold ) { hold = 0.001; } // 04/20/2001

        flu = (theta * flu + theta1 * *cu 
              * (*tt - tf) / *gold) * 0.5;
        fld = (theta * fld + theta1 * *cd 
              * (*td - tf) / hold) * 0.5;
      }

      f = ft - flu - fld;
      if( dabs( f ) <= err ) { break; }

      df = f - fwas;

      if( dabs( df ) < ZEROTOL ) { break; }

      gnew = *g - (*g - gwas) * f / df;
      gwas = *g;
      fwas = f;
      *g = gnew;
      
    }

    if( n == 50 ) 
    { 
//      rflog1 << "Error in phase plane" << endl;
      return 0;
//      exit( -1 ); 
    }

  }
  else
  {
    /* --- NO LATENT HEAT - uses ordinary implicit formulation */

    err = 1e-5f;
    
    for( n = 0; n < 50; ++n ) 
    {
      h = *dxx - *g;
      flu = ZERO;
      if ( *cu > ZERO ) 
      {
        // added 11/Dec/2000

        if( ZERO == (*cu + (c1 + d1 * *g) * *g) ) 
        {   
          *cu = -(c1 + d1 * *g) * *g + 0.01; // to protect
        }
 
        flu = *cu * (a1p + b1p * *g) / (*cu + (c1 + d1 * *g) * *g);
      }

      fld = ZERO;
      
      if( *cd > ZERO ) 
      {
        if( ZERO == (*cd + (c2 + d2 * h) * h) ) 
        {  
          *cd = -(c2 + d2 * h) * h + 0.001; // to protect 04/20/2001
        }

        fld = *cd * (a2p + b2p * h) / (*cd + (c2 + d2 * h) * h);
      }
	
      f = -flu - fld;

      if( dabs( f ) <= err ) { break; }

      df = f - fwas;

      if( dabs( df ) < ZEROTOL ) { break; }

      gnew = *g - (*g - gwas) * f / df;
      gwas = *g;
      fwas = f;
      *g = gnew;
    }

    if( 50 == n ) 
    { 
      rflog1 << "Failed to converge ";
      rflog1 << "- no latent heat in phase plane" << endl;
      
      exit( -1 ); 
//      return 0;
    }
  }

  if( ZERO == (*g * (c1 + d1 * *g) + *cu) ) 
  { 
    *cu = -*g * (c1 + d1 * *g) + 0.01; // to proctect
  }

  if( ZERO == (h * (c2 + d2 * h) + *cd) ) 
  {  
    *cd = -h * (c2 + d2 * h) + 0.001; // to protect 04/20/2001
  }

  *t1 = (*g * (a1 + b1 * *g) + cutf) / (*g * (c1 + d1 * *g) + *cu);
  *t2 = (h * (a2 + b2 * h) + cdtf) / (h * (c2 + d2 * h) + *cd);

  if( *g < ZERO ) 
  { 
    // Out at upper side
    
    if( dabs(*gold - *g) < ZEROTOL )  
    {
      dtfaz = dtday * SECSPERDAY; // debugging by Q.Z. 12/Dec/2000
    }
    else  { dtfaz = *gold * dtfaz / (*gold - *g); }

//     if ( dabs(*gold - *g) < ZEROTOL )  { *gold = *g +0.01; printf(" bad"); }  // 04/20/2001

//     dtfaz = *gold * dtfaz / (*gold - *g);

    *t1 = *tt;
  }
  else
  {
    if( *g <= *dxx ) { return 0; }

    // Out at lower side
    
//    if( (*g - *gold) > 10.0)   { *g = *gold +1.0;}

    if( dabs(*g - *gold) < ZEROTOL ) 
    {
    //     *g = *gold + 0.01;
      dtfaz = dtday * SECSPERDAY; // debugging by Q.Z. 12/Dec/2000
    } // to protect
    else
    {
      dtfaz = hold * dtfaz / (*g - *gold);
    }

    *t2 = *td;
  }

  return 0;

 };

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::updateSinglePhasePlaneTemps( integer *dnode,
                                                const integer& cmnt,
                                                ofstream& rflog1 )
{

// formerly ASMONE(I)

// Calculates new position of single phase plane
// Calculates new temperatures at all nodes

  double cd;
  double cnd;

  integer counter;

  double cu;

  double dens;
  double dtwas;
  double dxx;

  double g;
  double gold;

  integer im1;
  integer ip1;
  integer ip2;

  double r1;
  double rlw;
  double rsd;
  double rsu;

  double sph;

  double t1;
  double t2;

  double td;
  double tt;

  double wat;

 
  dtwas = dtfaz;
  gold = xfa[*dnode];
  im1 = *dnode - 1;
  ip1 = *dnode + 1;
  ip2 = *dnode + 2;

  counter = 0;
  while( counter < THOUSANDX100 )
  {
    // Determines Crank-Nicholson coefficients e and s for nodes
    //   is (i.e., surface) down to iml ( i.e., the node just
    //   above the node with the phase change plane) 	
    gaussianEliminForward( &im1 );

    // Determines Crank-Nicholson coefficients e and s for nodes
    //   ip2 (i.e., the second node below the node with the 
    //   change plane) down to imax1 ( i.e., the bottom node of 
    //   the soil profile)	
    
    gaussianEliminBackward( &ip2, cmnt );

    dxx = dx[*dnode];
    wat = water[*dnode];
    dens = ddry[*dnode];
    tt = t[*dnode];

    updateSoilThermalProperties( &mater[*dnode], 
                                 &tt,
                                 &sph, 
                                 &cnd, 
                                 rflog1 );

    cu = cnd;
    rsu = dens * sph;
    td = t[ip1];

    updateSoilThermalProperties( &mater[*dnode],
                                 &tf,
                                 &sph,
                                 &cnd,
                                 rflog1 );

    cd = cnd;
    rsd = dens * sph;
    rlw = dens * hlat * wat;

    updatePhasePlane( dnode, 
                      &rlw, 
                      &cu, 
                      &rsu, 
                      &cd, 
                      &rsd, 
                      &tt, 
                      &td, 
                      &dxx, 
                      &gold, 
                      &g, 
                      &t1, 
                      &t2, 
                      rflog1 );

    counter++;

    if( g < ZERO ) 
    { 

      // Out above
      
      capx[*dnode] = rsd * dxx;

      if( dxx < 0.01 ) { dxx = 0.01; } // to protect
      
      conx[*dnode] = cd / dxx;
      xfa[*dnode] = -1e10f;
      
      if( *dnode == is ) { return 0; }
      
      --(*dnode);
      im1 = *dnode - 1;
      ip1 = *dnode + 1;
      ip2 = *dnode + 2;
      r1 = tt - tf;

      if( r1 >= ZERO ) 
      { 
        t[ip1] = tf - ADJUST; 
      }
      else { t[ip1] = tf + ADJUST; }
 
      // Determines Crank-Nicholson coefficients e and s for 
      //   nodes "is" (i.e., surface) down to dnode 	

      gaussianEliminForward( dnode );

      // Calculates temperatures for nodes "is" down to dnode
       
      crankNicholsonBackward( dnode );

      // Determines Crank-Nicholson coefficients e and s for 
      //   nodes ip2 (i.e., the second node below dnode) 
      //   down to imax1 ( i.e., the bottom node of the soil 
      //   profile)	

      gaussianEliminBackward( &ip2, cmnt );

      // Calculates temperatures for nodes "ip2" (i.e., the 
      //   second node below dnode) down to imax

      crankNicholsonForward( &ip2 );

      dtfaz = dtwas - dtfaz;
      dtwas = dtfaz;
      gold = dx[*dnode];
    }
    else
    {
      if( g <= dxx ) 
      { 
        // g remains within element i
        
        t[*dnode] = t1;

        crankNicholsonBackward( &im1 );

        t[ip1] = t2;

        crankNicholsonForward( &ip2 );

        xfa[*dnode] = g;

        return 0;
      }
      else
      {
      	// Out below
      	
        if( dxx < 0.01 ) { dxx = 0.01; } // to protect

        capx[*dnode] = rsu * dxx;
        conx[*dnode] = cu / dxx;
        xfa[*dnode] = -1e10f;
        
        if( *dnode == imax1 ) { return 0; }
        
        ++(*dnode);
        im1 = *dnode - 1;
        ip1 = *dnode + 1;
        ip2 = *dnode + 2;
        r1 = td - tf;

        if( r1 >= ZERO ) 
        { 
          t[*dnode] = tf - ADJUST; 
        }
        else { t[*dnode] = tf + ADJUST; }

        // Determines Crank-Nicholson coefficients e and s for 
        //   nodes "is" (i.e., surface) down to iml ( i.e., the 
        //   node just above the node with the phase change 
        //   plane) 	

        gaussianEliminForward( &im1 );

        // Calculates temperatures for nodes "is" down to iml
  
        crankNicholsonBackward( &im1 );
        
        // Determines Crank-Nicholson coefficients e and s for 
        //   nodes ip1 (i.e., the node just below node with the
        //   phase change plane) down to imax1 ( i.e., the 
        //   bottom node of the soil profile)	

        gaussianEliminBackward( &ip1, cmnt );

        // Calculates temperatures for nodes "ip1" (i.e., the 
        //   node just below the node with the phase change 
        //   plane) down to imax

        crankNicholsonForward( &ip1 );

        dtfaz = dtwas - dtfaz;
        dtwas = dtfaz;
        gold = ZERO;
      }
    }
  }

  rflog1 << "Single phase plane not resolved ";
  rflog1 << "(counter reached 100,000)" << endl;
//  exit( -1 );  

  return 0;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::updateSnowPackGrid( const double& acum, 
                                       double *snoden, 
                                       double *weight, 
                                       double *comp, 
                                       double *fcmelt, 
                                       double *tdrive )
{

  // formerly NEIGE(ACUM,SNODEN,WEIGHT,COMP,ETAO,DENFAC,FCMELT,TDRIVE)

  // Calculates density and thickness at time m of snow layers is,
  //   is+1, .. igm1 where "is" is the top of the top snow layer and 
  //   "igml" is the top of the bottom snow layer

  // Also calculates compaction of surface layer "is"
  // Number of layers changes only when sufficient snow has fallen
  //   (i.e., a layer may correspond to several time steps)
  
  // comp = initial mass of a full layer
  
  // Includes accumulation, mechanical ablation and melting
  
  // Densification by successive iteration centered difference
  
  // All units are S.I.U. except acum and fcmelt (days)
  
  double ablat;  // ablation
  double dsmass; // Change in snowpack mass ( kg m^-2 )
  double dxold;
  double fonte;

  integer dnode;  // node index
  integer j;
  integer jj;

  double oldremain;

  double ratio;
  double remain;
  double rit;
  double tsno;
  double wtis;
  double wtisd;
  double xx;

  // Determine new snow mass from snowfall during time step
  
  dsmass = acum * dtday;

  // No snowpack or snowfall during time step?
  
  if( smass <= ZERO && dsmass <= ZERO ) 
  {
    return 0;
  }

  // Add new snow mass to snowpack mass
  
  smass += dsmass;
  //if (is <=0) is =0; //added by cgs to avoid segmentation fault for is <0;

  if( dsmass > ZERO ) 
  {  
    // Accumulation
     	
    tsno = t[is];
    
    // If no snowpack, create first snowpack node

    if( is == ig ) 
    {
      is = igm1;
    }

    wtis = weight[is] + dsmass;
    if( wtis <= *comp ) 
    { 
      // Add new snow mass to fill only a fraction of a node
      if( ZERO == *snoden ) { *snoden = 0.01; } // to protect

      if( ZERO == (dx[is] + dsmass / *snoden) ) 
      { 
        dx[is]= -(dsmass / *snoden)+ 0.01;
      } // to protect

      ddry[is] = wtis / (dx[is] + dsmass / *snoden);
      if( weight[is] <= ZERO ) 
      {
        xfa[is] = -1e10f;
      }

      weight[is] = wtis;
      t[is] = tsno;
    }
    else
    {
      // Add new snowpack nodes
      
      wtisd = *comp - weight[is];

      if( ZERO == *snoden ) { *snoden = 0.01; } // to protect
      if( ZERO == (dx[is] +  wtisd / *snoden) ) 
      { 
        dx[is] =  -( wtisd / *snoden) + 0.01; // to protect
      }

      ddry[is] = *comp / (dx[is] +  wtisd / *snoden);
      weight[is] = *comp;
      remain = dsmass - wtisd;
      t[is] = tsno;

      do {
        --is;
        if (is <0) is =0;//added by cgs to avoid segmentation fault for is <0;
        weight[is] = *comp;
        ddry[is] = *snoden;
        t[is] = tsno;
        xfa[is] = -1e10f;
        oldremain = remain;
        if( remain <= *comp ) 
        { 
          weight[is] = remain;
        }
        else 
        { 
          remain -= *comp; 
        }
      } while( oldremain > *comp );
    }
  }
  else
  {
    // Mechanical ablation
    
    ablat = -dsmass;

    while( ablat >= weight[is] )  
    {
      if( is == ig ) 
      { 
        // Complete disappearance of snow cover
        
        smass = ZERO;
        weight[igm1] = ZERO;
        xfa[igm1] = -1e10f;
        
        return 0;
      }

      ablat -= weight[is];
      weight[is] = ZERO;
      dx[is] = ZERO;
      xfa[is] = -1e10f;
      ++is;
    }
 
     weight[is] -= ablat;
  }

  if( *tdrive > ZERO ) 
  { 
    // Melt ablation
    
    fonte = *fcmelt * dtday * tstair;
    smass -= fonte;
    
    if( smass > ZERO ) { *tdrive = 0.01; }

    while( fonte >= weight[is] )  
    {
      if( is == ig ) 
      { 
        // Complete disappearance of snow cover
        
        smass = ZERO;
        weight[igm1] = ZERO;
        xfa[igm1] = -1e10f;
        
        return 0;
      }

      fonte -= weight[is];
      weight[is] = ZERO;
      dx[is] = ZERO;
      xfa[is] = -1e10f;
      ++is;
    }

     weight[is] -= fonte;
  }

  // Compaction of all layers including surface layer (is)
  // Assumes constant mass of each layer
  // sload includes one-half weight of node 1

  for( dnode = is; dnode <= igm1; ++dnode ) 
  {
    // Calculate layer thickness (assume constant mass)
    
    // DWK changed mater[dnode-1] = 21 to mater[dnode] = 4
    //   on 20040504
    	   	
    mater[dnode] = 4;  
    rit = ddry[dnode];
    dxold = dx[dnode];

    if( ZERO == rit ) { rit = 0.01; } // to protect
    dx[dnode] = weight[dnode] / rit;

    if( xfa[dnode] >= ZERO ) 
    { 
      // Update phase boundaries for changed dx[dnode-1]
      
      if( ZERO == dxold ) { dxold = 0.01; } // to protect
      
      ratio = dx[dnode] / dxold;
      
      if ( dabs( xfb[dnode] ) >= ZEROTOL ) 
      { 
        xfb[dnode] += (ratio - 1.0) * (dx[dnode] - xfa[dnode]);
        
        if( xfb[dnode] < ZERO ) 
        {
          xfb[dnode] = dx[dnode] * 0.01;
        }
      }

      xfa[dnode] = ratio * xfa[dnode];
    }
  }

  // Update heights of nodes
  
  xx = x[ig];
  jj = is + igm1;

  for( j = is; j <= igm1; ++j ) 
  {
    dnode = jj - j;
    xx -= dx[dnode];
    x[dnode] = xx;
  }

  // Restrain snow temperatures and ground surface temperatures
  
  for( dnode = is; dnode <= ig; ++dnode ) 
  {
    if( t[dnode] > ZERO ) 
    {
      t[dnode] = 0.01;
    }
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::updateSoilProfileHeat( const int& m,
                                          integer *kint, 
                                          double *source )
{

  // formerly INHET(KINT,HEATT)
  // Updates internal heat source vectors
  // Formulation allows for point sources at depths xhet[]
  // Elements having value for xhet[] other than 999.9 have heat
  //   sources numbered 1 ... nheat in increasing order of depth

  integer dnode;
//  int dnode1;

  double hti;
  double htip1;

  integer iheat;
  integer ip1;

  integer nheat;


  for( iheat = 0; iheat < 2; ++iheat ) 
  {
    heatt[iheat] = heat[m][iheat];
  }

  // Initialize ht[] and update htold[]

  for( dnode = is; dnode <= imax; ++dnode ) 
  {
    htold[dnode] = ht[dnode];
    ht[dnode] = ZERO;
  }

  // ASSEMBLE NEW NODAL HEATING VECTOR 

  nheat = 0;
 
  for( dnode = ig; dnode <= imax1; ++dnode ) 
  {
    ip1 = dnode + 1;
    if ( xhet[dnode] != 999.9 ) 
    {
      // ELEMENT CONTAINS A POINT HEAT SOURCE 

      ++nheat;
      if( nheat > 2 ) { return 0; }

      if( 1 == kint[nheat-1] ) 
      { 
        if( t[dnode] < -10.0 && t[ip1] < -10.0 ) 
        {
          source[nheat-1] = ZERO;
        }

        if( ZERO == dx[dnode] ) 
        { 
          dx[dnode] = 0.01; // to protect
        }

        hti = source[nheat-1] * (x[ip1] 
              - xhet[dnode]) / dx[dnode];
        htip1 = source[nheat-1] * (xhet[dnode] - x[dnode])
	        / dx[dnode];
        ht[dnode] += hti;
        ht[ip1] += htip1;
      }
    }
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::updateSoilThermalProperties( integer *mat, 
                                                double *ta,  
                                                double *sph, 
                                                double *cnd, 
                                                ofstream& rflog1 )
{
 
  // formerly PRAM(MAT,TA,TB,DENSE,WAT,SPH,CND)
  
  if( *mat > 3 ) 
  { 
    if( *mat > (MAXMAT-1) ) 
    { 
      rflog1 << "Material index (mat) = " << *mat; 
      rflog1 << " is undefined" << endl;	
      
      exit( -1 ); 
    }

    setSnowCharSturm( ta, sph, cnd );
    
    return 0; 
  }

  if( *ta <= tf ) 
  { 
    *cnd = condf[*mat];
    *sph = sphf[*mat];
  }
  else
  {
    *cnd = condt[*mat];
    *sph = spht[*mat];
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soilthermal60::updateTimeStep( const int& m,
                                   const integer& cmnt,
                                   double& r1,
                                   const integer& imaxp,
                                   ofstream& rflog1 )
{

  double acum;
  double cnd;

  integer dnode;
  integer dnode1;
  integer ibzero;
  integer itzero;

  integer j;
  integer jj;

//  double sflx1;
//  double sflx2;
//  double sflx3;

  double snoden;
  double sph;

  double tdrive;
  double topnew;
  double topold;
  double tt;


  // EVALUATE NEW AIR TEMPERATURE (TSTAIR) 

  tstair = taer[m];
  tdrive = tstair;

  dtfaz = dt;

  // SNOW COVER GRID

  acum = snow[m];
  if( tstair > ZERO ) {  acum = ZERO; }
 
  snoden = sden[m];

  //if (is <=1) cout <<"check1: "<<" is: " <<is<<endl;
  updateSnowPackGrid( acum, 
                      &snoden, 
                      weight, 
                      &comp,
                      &fcmelt, 
                      &tdrive );

  //if (is <=1) cout<<"check2: "<<"is: " <<is<<endl;

  topold = t[is];
  ism1 = is - 1;
  t[ism1] = tdrive;
  e[ism1] = tdrive;
  s[ism1] = ZERO;
 

  // UPDATE INTERNAL HEAT SOURCES 

  updateSoilProfileHeat( m, kint, heatt );


  // EVALUATE NEW BOTTOM CONDITIONS 

  tbot = tber[m];
  t[imaxp] = tbot;


  // EVALUATE THERMAL PROPERTIES 

  tt = t[is];

  for( dnode = is; dnode <= imax1; ++dnode ) 
  {
    dxx = dx[dnode];
    td = t[dnode+1];

    // Determine new thermal conductivity and volumetric 
    //   heat capacity of each node in the snowpack and
    //   soil based on current thermal and moisture
    //   conditions
    

    updateSoilThermalProperties( &mater[dnode], 
                                 &tt, 
                                 &sph, 
                                 &cnd, 
                                 rflog1 );

    tt = td;
    if( dxx < 0.01 ) { dxx = 0.01; }; // to protect

    // Determine CNi for each node (see Figure 1 in 
    //   Goodrich [1978b])
    
    conx[dnode] = cnd / dxx;
    
    // Determine CAPi for each node (see Figure 1 in 
    //   Goodrich [1978b])
    
    capx[dnode] = ddry[dnode] * sph * dxx;
  }

  // Determine CNi for node just above surface of snowpack
  //   or soil (if snowpack does not exist)

  conx[ism1] = htop;

  // Determine CAPi for node just above surface of snowpack
  //   or soil (if snowpack does not exist)

  capx[ism1] = ZERO;

  if ( hlat <= ZERO ) 
  { 
    // NO LATENT HEAT FOR ANY ELEMENT 

    // Determine Crank-Nicholson coefficients e[] and s[] for 
    //   all nodes in the snowpack and soil
    
    gaussianEliminBackward( &is, cmnt );
    
    // Determine new temperatures for each node in the 
    //   snowpack and soil
    
    crankNicholsonForward( &is );
  }
  else
  {
    // LOCATE UPPER & LOWER PHASE PLANES 

    for( dnode = is; dnode <= imax; ++dnode ) 
    {
      // itzero indicates the node location of the upper
      //   phase change plane
      	
      itzero = dnode; 
      if( xfa[dnode] >= ZERO ) 
      {
        break;
      }
    }

    jj = ism1 + imax1;
    for( j = ism1; j <= imax1; ++j ) 
    {
      dnode = jj - j;

      // ibzero indicates the node location of the lower
      //   phase change plane

      ibzero = dnode;
      if( xfa[dnode] >= ZERO ) 
      {
        break;
      }
    }

    // dnode1 is an index to indicate the number of 
    //   phase change planes in the snowpack and soil:
    //   - negative values indicate no phase change plane
    //   - zero indicates a single phase change plane
    //   - positive values indicate multiple phase change planes
    
    dnode1 = ibzero - itzero;

    if( dnode1 < 0 ) 
    { 
      // No phase change plane exists in the snowpack and soil
      
      // Determine Crank-Nicholson coefficients e[] and s[] for 
      //   all nodes in the snowpack and soil

      gaussianEliminBackward( &is, cmnt );

      topnew = tdrive * s[is] + e[is];
      
      if( (topnew - tf) * (topold - tf) >= ZERO ) 
      {
        // Update temperatures for all nodes in the snowpack
        //   and soil 
        
        crankNicholsonForward( &is );
        
        return 0;
      
      }
      else
      {
        // Create single phase change plane at the surface of 
        //   the snowpack or soil (if no snowpack exists)
        
        itzero = is;
        ibzero = is;
        createSinglePhasePlane( cmnt, 
                                itzero, 
                                r1, 
                                topold, 
                                topnew );
      }
    } 

    if( dnode1 <= 0 && xfb[itzero] <= ZERO )   
    { 
      // SINGLE PHASE PLANE 

      // Calculate temperatures for all nodes and determine
      //   new position of the phase change plane

      updateSinglePhasePlaneTemps( &itzero, cmnt, rflog1 );
    }
    else
    {
      // MORE THAN ONE PHASE PLANE

      // Calculate temperatures for all nodes and determine
      //   new positions of the phase change planes
      
      updateMultiplePhasePlaneTemps( &itzero, 
                                     &ibzero, 
                                     cmnt, 
                                     rflog1 );
    }

    topnew = t[is];


    /* --- START OF SECOND PHASE PLANE AT SURFACE */
    /* --- TEMPERATURES SET TO ZERO BETWEEN PHASE PLANES */

    while((topnew - tf) * (topold - tf) < ZERO ) 
    {
      createMorePhasePlanes( itzero, 
                             ibzero, 
                             r1, 
                             topold, 
                             topnew );    

      updateMultiplePhasePlaneTemps( &itzero, 
                                     &ibzero, 
                                     cmnt, 
                                     rflog1 );
                                     
      topnew = t[is];
    }
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soilthermal60::updateyrCDM( double tair[CYCLE] )
{
  // Updates "cooling degree month" (CDM) each year
  
  int dm;

  CDM = ZERO;
  for ( dm = 0; dm < CYCLE; ++dm )
  {
    CDM += (10.0 - tair[dm]);
  }

};


