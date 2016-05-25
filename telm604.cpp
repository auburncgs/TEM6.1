/* **************************************************************
*****************************************************************
TELM604.CPP - Runs TEM for a single grid cell

Modifications:

20030712 - DWK created by modifying telm50b1.cpp
20030712 - DWK changed TEMelmnt:: to TEMelmnt51::
20030712 - DWK added clm.d40[], clm.o3[] and cln.ndep[] to
           atmswritepred()
20030712 - DWK added double d40[CYCLE] and double ndep[CYCLE] to
           temgisqc()
20030712 - DWK added tem.atms.d40[CYCLE] and tem.atms.ndep[CYCLE]
           to temgisqc() in setTEMequilState()
20030712 - DWK added tem.veg.fprevozone to gtem.restoreTEMState()
           and saveTEMState() in setTEMequilState(),
           updateTEMmonth() and updateTEMyear()
20030712 - DWK added tem.I_FOZONE and tem.I_FINDOZONE to
           temwritemiss() and temwritepred()
20030712 - DWK changed tem.microbe.nup to tem.microbe.nimm in
           restoreTEMState() and saveTEMState() function calls
           in setTEMequilState(), updateTEMmonth() and
           updateTEMyear()
20030712 - DWK added tem.atms.ndep[], tem.veg.nfix[],
           tem.microbe.grossmin[], tem.microbe.n2oprod[],
           tem.microbe.no3prod[], tem.microbe.nitrif[],
           tem.microbe.ammonif[], tem.microbe.denitrif[] and
           tem.soil.leach[] to temwritepred()
20030712 - DWK renamed microbe.nuptake[] to microbe.immob[] in
           temwritepred()
20030712 - DWK added tem.I_TOTVNFX and tem.I_AGVNFX to
           temwritepred()
20030712 - DWK changed ofstream ftemout[MAXPRED] to ofstream
           ftemout[NUMTEM} in setTEMequilState(), setTEMmiss(),
           temwritemiss(), temwritepred(), updateTEMmonth() and
           updateTEMyear()
20030712 - DWK changed TEMList50 gtem to TEMList51 gtem in
           function calls to setTEMequilState(), updateTEMmonth()
           and updateTEMyear()
20030712 - DWK added spinoutfg = 2 condition to setTEMmiss(),
           updateTEMmonth() and updateTEMyear()
20031106 - DWK added tem.qualcon[outyr] = 0; to temwritepred()
20031231 - DWK added initialization of veg.prvleafmx to 
           setTEMequilState()
20031231 - DWK deleted tem.veg.cmax, tem.microbe.kdc, 
           tem.microbe.decay, tem.veg.nmax and tem.microbe.nimm 
           from functions calls of restoreTEMState() and 
           saveTEMState() in setTEMequilState(), updateTEMmonth()  
           and updateTEMyear() 
20040229 - DWK changed TEMelmnt51:: to TEMelmnt60::
20040714 - DWK changed char predmap[MAXPRED][9] to 
           string predmap[MAXPRED] in atmswritemiss(),
           atmswritepred(), setTEMequilState(), setTEMmiss(),
           temwritemiss(), temwritepred(), updateTEMmonth() and
           updateTEMyear()
20040716 - DWK changed Clmdata atmspred to Clmdata50 atmspred
           in atmswritemiss() and atmswritepred()
20040716 - DWK changed Temdata tempred to Temdata43 tempred
           in temwritemiss() and temwritepred()
20040716 - DWK changed char varname1[9] to const string& varname1  
           in coregerr()
20040716 - DWK changed char varname2[9] to const string& varname2  
           in coregerr()
20040716 - DWK changed Soildata fao to Soildata43 fao in 
           setGIStopography()
20040716 - DWK changed Elevdata elv to Elevdata43 elv in 
           setGIStopography()
20040716 - DWK added subarea to function call of 
           gtem.restoreTEMState() and gtem.saveTEMState() in 
           setTEMequilState(), updateTEMmonth() and
           updateTEMyear()           
20040716 - DWK added subarea to function call of 
           gtem.setTEMStateConst() in setTEMmiss()
20040716 - DWK added subarea to function call of 
           tempred.poutdel() and tempred.outdel() in 
           temwritemiss() and temwritepred()           
20040716 - DWK changed TEMList51& gtem to TEMList60& gtem in 
           setTEMequilState(), setTEMmiss(), updateTEMmonth() 
           and updateTEMyear()
20040716 - DWK changed InorgN ndep[CYCLE] to InorgN60 ndep[CYCLE]
           in temgisqc()
20041003 - DWK added tem.ntns[] to temwritepred()
20041003 - DWK added spinoutyrs to setTEMequilState(),
           setTEMmiss(), updateTEMmonth(), and updateTEMyear()
20050409 - DWK added tem.microbe.decomp[] and 
           tem.microbe.ndecomp[] to temwritepred()
20051117 - DWK added inlcud telm602.h and standard includes
20051117 - DWK changed Soildata43 to Soildata60 in 
           setGIStopography()
20051117 - DWK changed Elevdata43 to Elevdata60 in 
           setGIStopography()
20051118 - DWK added vector<string> predstr to TEMelmnt60()
20051119 - DWK deleted TEMList60& gtem from setTEMequilState(),
           updateTEMmonth() and updateTEMyear()
20051119 - DWK added public functions restoreCohortState() and
           saveCohortState()
20051122 - DWK added public functions readCohortState(),
           getTEMCohortState(), saveTEMCohortState()
           and writeCohortState()
20051123 - DWK added public function TEMequilibrium()
20051123 - DWK added public function outputTEMmonth()
20051128 - DWK added public function initializeTEMCohortState()
20060607 - DWK changed float col to double col in coregerr()
20060607 - DWK changed float row to double row in coregerr()
20060909 - DWK changed include from telm602.h to telm603.h
20070514 - DWK changed include from telm603.h to telm603a.h
20070514 - DWK added int standage to readCohortState(), 
           temwritepred() and writeCohortState()
20070830 - DWK changed include from telm603b.h to telm603c.h
20090127 - DWK changed include from telm603c.h to telm604.h
                                   
****************************************************************
************************************************************* */

#include<cstdio>

  using std::fscanf;
  using std::FILE;
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

#include<iomanip>

  using std::setprecision;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std::atoi;

#include<cmath>

  using std::exp;
  using std::fabs;
  using std::pow;
#include<math.h>
#include <float.h>
#include <limits>

#include<vector>

  using std::vector;
      
#include<string>
  
  using std::string;
using namespace std;

#include "telm604.h"

/* *************************************************************
************************************************************* */
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
TEMelmnt60::TEMelmnt60()
{

  col = MISSING;
  row = MISSING;
  carea = -999;
  subarea = -999;
  fatalerr = 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMelmnt60::atmswritemiss( ofstream fout[NUMATMS],
                                const vector<string>& predname,
                                const int& pdyr,
                                const int& natmspred,
                                const double value )
{

  int i;
  int dm;
  Clmdata60 atmspred;

  for( i = 0; i < natmspred; ++i )
  {
    for( dm = 0; dm <= CYCLE; ++dm )
    {
      atmspred.mon[dm] = value;
    }

	 atmspred.outdel( fout[i],
                          col,
                          row,
                          predname[i],
                          carea,
	                  atmstotyr[pdyr],
                          atmspred.mon,
                          contnent );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMelmnt60::atmswritepred( ofstream fout[NUMATMS],
                                const vector<string>& predname,
                                const int& natmspred )
{

  int i;
  int dm;
  Clmdata60 atmspred;

  // Covert cal/cm2/day to W/m2 (4.186 Joules / calorie)

  const double  cal2Watts = 0.4845;


  for( i = 0; i < natmspred; ++i )
  {
    if( predname.at( i ) == clm.predstr.at( clm.I_GIRR ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm )
      {
        atmspred.mon[dm] = climate[clm.I_GIRR][dm] * cal2Watts;
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_NIRR ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm )
      {
        atmspred.mon[dm] = climate[clm.I_NIRR][dm] * cal2Watts;
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_PAR ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm )
      {
        atmspred.mon[dm] = climate[clm.I_PAR][dm] * cal2Watts;
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_CLDS ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm ) 
      { 
        atmspred.mon[dm] = climate[clm.I_CLDS][dm]; 
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_TAIR ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm ) 
      { 
      	atmspred.mon[dm] = climate[clm.I_TAIR][dm]; 
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_PREC ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm ) 
      { 
      	atmspred.mon[dm] = climate[clm.I_PREC][dm]; 
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_RAIN ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm ) 
      { 
      	atmspred.mon[dm] = climate[clm.I_RAIN][dm]; 
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_SNWFAL ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm ) 
      { 
      	atmspred.mon[dm] = climate[clm.I_SNWFAL][dm]; 
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_AOT40 ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm ) 
      { 
      	atmspred.mon[dm] = climate[clm.I_AOT40][dm]; 
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_TNDEP ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm )
      {
        atmspred.mon[dm] = climate[clm.I_TNDEP][dm] * 1000.0;
      }
    }

    else if( predname.at( i ) == clm.predstr.at( clm.I_NH4DEP ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm )
      {
        atmspred.mon[dm] = climate[clm.I_NH4DEP][dm] * 1000.0;
      }
    }
    else if( predname.at( i ) == clm.predstr.at( clm.I_NO3DEP ) )
    {
      for( dm = 0; dm <= CYCLE; ++dm )
      {
        atmspred.mon[dm] = climate[clm.I_NO3DEP][dm] * 1000.0;
      }
    }
    else
    {
      for( dm = 0; dm <= CYCLE; ++dm )
      {
        atmspred.mon[dm] = MISSING;
      }
    }

    atmspred.outdel( fout[i],
                     col,
                     row,
                     predname[i],
                     carea,
                     year,
                     atmspred.mon,
                     contnent );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMelmnt60::coregerr( ofstream& rflog1,
                          const string& varname1,
                          const double& col1,
                          const double& row1,
                          const string& varname2,
                          const double& col2,
                          const double& row2 )
{

  int fatalerr = 0;

  if( col1 != col2 || row1 != row2 )
  {
    fatalerr = 1;

    cout << "ERROR:  " << varname1 << " data and ";
    cout << varname2 << "data are not coregistered." << endl;
    cout << "COL = " << col1 << " and ROW = " << row1;
    cout << " in " << varname1 << " data" << endl;
    cout << "COL = " << col2 << " and ROW = " << row2;
    cout << " in " << varname2 << " data" << endl;

    rflog1 << "ERROR:  " << varname1 << " data and ";
    rflog1 << varname2 << "data are not coregistered." << endl;
    rflog1 << "COL = " << col1 << " and ROW = " << row1;
    rflog1 << " in " << varname1 << " data" << endl;
    rflog1 << "COL = " << col2 << " and ROW = " << row2;
    rflog1 << " in " << varname2 << " data" << endl;
  }

  return fatalerr;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

int TEMelmnt60::equilibrateTEM( const int& pichrt, 
                                const double& ptol, 
                                ofstream& rflog1 )
{
  // Run TEM until steady state conditions occur (equilibrium)

  int dyr = 0;
  int dm;

  // Initialize standing stocks of carbon and nitrogen from 
  //   calibration ("ECD") data and fluxes for integrator

  tem.ECDsetODEstate( tem.veg.cmnt, tem.soil.getPSIPLUSC() );
   
  // Set previous value of TEM ODE state variables to the 
  //   current values of TEM ODE state variables for initial
  //   conditions
 
  tem.setPrevState();


  // Initialize all disturbance fluxes to zero

  tem.ag.resetMonthlyDisturbFluxes();

  tem.ag.setFIRENDEP( ZERO );


  // Initialize all agricultural and wood product pools to 
  //   zero

  tem.ag.resetPROD();


  tem.totyr = 0;
  tem.endeq = 0;
  tem.intflag = 0;
  tem.initFlag = 0;
  
  // Determine CDM for long-term mean climate data for 
  //   soil thermal calculations (STM)

  tem.soil.stm.updateyrCDM( climate[clm.I_TAIR] ); 


  // Initialize tem.atms.prevtair for water balance (WBM) and 
  //   soil thermal (STM) calculations
  
  tem.atms.setPREVTAIR( climate[clm.I_TAIR][CYCLE-1] );
  tem.atms.setPREVRAIN( climate[clm.I_RAIN][CYCLE-1] );

  // Initialize tem.atms.prev2tair for water balance  
  //   calculations (WBM)

  tem.atms.setPREV2TAIR( climate[clm.I_TAIR][CYCLE-2] );

  
  // Initialize tem.atms.prevtair for water balance (WBM) and 
  //   soil thermal (STM) calculations
  
  tem.soil.setPREVSPACK( ZERO );
  
  
  // Initialize tem.soil.prevdst10 for GPP calculations
  
  tem.soil.setPREVDST10( ZERO );
  
  
  // Initialize tem.veg.prvleafmx and tem.veg.prevunrmleaf
  //   for phenology calculations
  
  tem.veg.setPRVLEAFMX( tem.veg.getINITLEAFMX( tem.veg.cmnt ) );
  
  tem.veg.setPREVUNRMLEAF( tem.veg.getUNLEAF12( tem.veg.cmnt ) );  

  tem.setLON(col); //added by cgs to pass the latitude and longitude information to tem
  tem.setLAT(row);

  while( (dyr < tem.runsize) && (tem.endeq < 2) )
  {
	  double tempxx=0.0;
	  double tempyy=0.0;
	  double tempzz=0.0;
	  double tempmm=0.0;
	  double tempnn=0.0;
	  double tempoo=0.0;
	  double temppp=0.0;
	for( dm = 0; dm < CYCLE; ++dm )
    {
      // Assign telmnt[0].climate to TEM.atms
      
      tem.atms.setGIRR( climate[clm.I_GIRR][dm] );
      tem.atms.setCLDS( climate[clm.I_CLDS][dm] );
      tem.atms.setNIRR( climate[clm.I_NIRR][dm] );
      tem.atms.setPAR(  climate[clm.I_PAR][dm] );
      tem.atms.setTAIR( climate[clm.I_TAIR][dm] );
      tem.atms.setRAIN( climate[clm.I_RAIN][dm] );
      tem.atms.setPREC( climate[clm.I_PREC][dm] );
      tem.atms.setSNOWFALL( climate[clm.I_SNWFAL][dm] );
      tem.atms.setCO2( climate[clm.I_CO2][dm] );
      tem.atms.setAOT40( climate[clm.I_AOT40][dm] );
      tem.atms.setTOTNDEP( climate[clm.I_TNDEP][dm] );
      tem.atms.setNH4DEP( climate[clm.I_NH4DEP][dm] );
      tem.atms.setNO3DEP( climate[clm.I_NO3DEP][dm] );
      tem.atms.setMXTAIR( mxtair );

      tem.soil.stm.setNEXTTAIR( climate[clm.I_TAIR][dm+1] );
      tem.soil.stm.setNEXTSNOWFALL( climate[clm.I_SNWFAL][dm+1] );

	  //cout <<"dyr: " <<dyr <<"dm: " <<dm<<"tair: " <<climate[clm.I_TAIR][dm]<<" NIRR: " <<climate[clm.I_NIRR][dm]<<" atmtair: " <<tem.atms.getTAIR()<<endl;

//      if( (CYCLE-1) == dm )
//      {
//        tem.soil.setNEXTDST10( cohort[pichrt].dst10[0] );
//      }
//      else
//      {
        tem.soil.setNEXTDST10( cohort[pichrt].dst10[dm+1] );
//      }
        //cout <<"cmnt: "<< tem.veg.cmnt <<" leafmax: "<<tem.veg.getINITLEAFMX( tem.veg.cmnt )<< " tair: " <<climate[clm.I_TAIR][dm]<<" co2: " <<climate[clm.I_CO2][dm]<< " ndep: " <<climate[clm.I_TNDEP][dm]<< endl;

       tem.endeq = tem.stepmonth( dyr,
                                 dm, 
                                 tem.intflag, 
                                 ptol, 
                                 rflog1 );
     
		
      // Update telmnt[0].cohort.dst10 for future use 
      //double xx= tem.getY(I_VEGC);
      //cout <<" xx: "<<xx<<endl;
      //if (xx != xx) cout <<" I_soilC1:" <<xx<<" col: " <<col<< "row: "<<row <<"value of soilc1: "<<is_valid(xx)<< " or: "<<is_nan(xx) <<" or: " <<is_infinite(xx)<<endl;
      //cohort[pichrt].dst10[dm] = tem.soil.getDST10();
      
      // Save TEM output to telmnt[0].output	  
      
      outputTEMmonth( dm );
      //tempxx+=tem.getNTCB();
      //tempyy+=tem.veg.getNFIX();
      //tempzz+=tem.soil.getEET();
      //tempmm+=tem.veg.getNUPTAKE();
      //tempnn+=tem.veg.getGPP();
      //tempoo+=tem.veg.getNPP();
      //temppp+=tem.microbe.getRH();
      //if (col ==-83.5 && row == 36.5 && tem.veg.cmnt==5) cout << "year: " <<dyr<<" cmnt: "<<tem.veg.cmnt<<" vegc: "<<tem.getY(I_VEGC)<<" TSOILC: "<<tem.soil.getTSOLC()<<" gpp: " <<tem.veg.getGPP() << " gpr: "<<tem.veg.getGPR() <<" litterc: "<<tem.veg.getLTRFALC()<<endl;
      //if (col == -80.75 && row== 51.25 && tem.veg.cmnt==4 && dm==11) cout << "year: " <<dyr<<" mon: "<<dm<<" cmnt: "<<tem.veg.cmnt<<" vegc: "<<tem.getY(I_VEGC)<<" GPP: "<< tempnn <<" NPP: "<<tempoo<<" RH: "<<temppp<<" ntcb: "<<tempxx<<endl;
      //if (col ==-83.75 && row == 36.75 && dyr > 1000) exit(-1);
      //if (tem.veg.cmnt==7 && dm==11) cout <<"dyr: "<<dyr<<" dm: "<<dm<<" totalntcb: "<<tempxx<<" nfix: "<<tempyy<<" gpp: "<<tempnn<<" avn: "<<tem.soil.getAVLN()<<" nup: "<<tempmm<<" vegc: "<<tem.getY(I_VEGC)<<" TSOILC: "<<tem.soil.getTSOLC()<<" NEP: "<<tem.getNEP()<<" nce: "<<tem.getNCE()<<" ntcb: "<<tem.getNTCB()<<endl;
      //cout <<"pdyr: "<<dyr<<" vegtype: "<<tem.veg.cmnt<<" vegc: "<<tem.getY(I_VEGC)<<" stemc: "<<output[tem.I_STEMC][dm]<<" rootc :"<< output[tem.I_ROOTC][dm]<<" storageN: "<<tem.getY(I_STON) << " structN: "<<tem.getY(I_STRN)<<endl;
      //if (tem.veg.cmnt==4 && col == -127.75 && row ==55.25&& dyr<10) cout <<" col: "<<col<<" row: "<<row<<" ttotyr: "<<dyr<<" pdm: "<<dm<<" cohort: "<<pichrt<<" cmnt: "<<tem.veg.cmnt<<" vegc: "<< output[tem.I_VEGC][dm]<<" vegc2 :"<<output[tem.I_FOLC][dm]+ output[tem.I_STEMC][dm]+output[tem.I_CROOTC][dm]+output[tem.I_FROOTC][dm] <<" vegn: "<<output[tem.I_STRN][dm]<<" ston: "<<output[tem.I_STON][dm]<<" stemc: "<<output[tem.I_STEMC][dm]<<" leafc: "<<output[tem.I_FOLC][dm]<<" crootc: "<<output[tem.I_CROOTC][dm]<<" frootc: "<<output[tem.I_FROOTC][dm]<<" soc: "<<output[tem.I_SOLC][dm]<<" agr: "<<output[tem.I_AGR][dm]<<" agl: "<<output[tem.I_AGL][dm]<<" bgr: "<<output[tem.I_BGR][dm]<<" bgl: "<<output[tem.I_BGL][dm]<<" CWD: "<<output[tem.I_CWD][dm]<<" ntcb: "<<output[tem.I_NTCB][dm]<<" GPP: "<<output[tem.I_GPP][dm] <<" NPP: "<<output[tem.I_NPP][dm]<<" RH: "<<output[tem.I_RH][dm]<<endl;
      tempxx += output[tem.I_NPP][dm];
      tempyy += output[tem.I_NTCB][dm];
      tempzz += tem.veg.getGPP();
      tempmm += output[tem.I_NEP][dm];
      tempnn += output[tem.I_NCE][dm];
      tempoo += output[tem.I_RH][dm];
      temppp = output[tem.I_LCHDOC][11];

      //if (tem.veg.cmnt==7&&dyr>5990){
    	  //cout <<"year: "<<dyr<<" mon: "<<dm<<" tsoil: "<<tem.soil.getTSOIL()<<" dst10: "<<tem.soil.getDST10()<<" tair: "<<tem.atms.getTAIR()<<" moistlim: "<<tem.veg.getGV()<<" dq10: "<<tem.microbe.getDQ10()<<" temp: "<<tem.veg.getTEMP()<<" repq10: "<<tem.veg.getNewRESPQ10()<<" npp: "<<output[tem.I_NPP][dm]<<" gpp: "<<output[tem.I_GPP][dm]<<" foliage: "<<tem.veg.getFOLIAGE()<<" leaf: "<<tem.veg.getLEAF()<<" thawp: "<<tem.veg.getTHAWPCT()<<" par: "<<tem.atms.getPAR()<<" predst10: "<<tem.soil.getPREVDST10()<<" nextdst10: "<<tem.soil.getNEXTDST10()<<endl;
    	  //cout <<"year: "<<dyr<<" mon: "<<dm<<" unnormalleaf: "<<tem.veg.getUNNORMLEAF()<<"prvmaxleaf: "<< tem.veg.getPRVLEAFMX()<<" leaf: "<<tem.veg.getLEAF()<<" pet: "<<tem.atms.getPET()<<" eet: "<<tem.soil.getEET()<<" prveetmx: "<<tem.soil.getPRVEETMX()<<endl;
      //}
      //int dddm = dm-1;
      //if (dddm <=0) dddm =0;
      //if (tem.veg.cmnt==7) cout <<"year: "<<dyr<<" dm: "<<dm<<" dq10: "<<tem.microbe.getdq10()<<" soc: "<< output[tem.I_SOLC][dm]<<" moisture: "<<tem.veg.getGV()<<" vegc: "<<output[tem.I_VEGC][dm]<<" gpp: "<<tem.veg.getGPP()<<" npp: "<<tem.veg.getNPP()<<" gpr: "<<tem.veg.getGPR()<<" litc: "<<tem.veg.getLTRFALC()<<" prev_VEGC: "<<output[tem.I_VEGC][dddm]<<endl;
    }
    //if (tem.veg.cmnt==7 && cohort[pichrt].currentveg==50)&& col == -143.25 && row ==64.25
	//&& col == -89.75 && row ==51.5
	//if (tem.veg.cmnt==7&&dyr>5990){
    //cout <<" col0: "<<col<<" row: "<<row<<" year: "<<dyr<<" dm: "<<dm<<" cohort: "<<pichrt<<" currentveg: "<< cohort[pichrt].currentveg<<" cmnt: "<<tem.veg.cmnt<<" vegc11: "<< output[tem.I_VEGC][11]<<" leafc11: "<<output[tem.I_FOLC][11] <<" stemc11: "<<output[tem.I_STEMC][11] <<" SOLC11: "<< output[tem.I_SOLC][11]<<" soilc: "<<output[tem.I_SOC][11]<< " NPP: "<<tempxx<<" GPP: "<<tempzz<<" rh: "<<tempoo<<" docleach: "<<temppp<<" ntcb: "<<tempyy<<" nep: "<<tempmm<<" nce: "<<tempnn<<endl;
    //cout <<" year: "<<dyr<<" dm: "<<dm<<" coho: "<<pichrt<<" cmnt: "<<tem.veg.cmnt<<" vegc: "<< output[tem.I_VEGC][7]<<" soc: "<< output[tem.I_SOLC][6]<<" vegc2 :"<<output[tem.I_FOLC][7]+ output[tem.I_STEMC][7]+output[tem.I_CROOTC][7]+output[tem.I_FROOTC][7] <<" stemc: "<<output[tem.I_STEMC][7]<<" leafc: "<<output[tem.I_FOLC][7]<<" crootc: "<<output[tem.I_CROOTC][7]<<" frootc: "<<output[tem.I_FROOTC][7]<<" agr: "<<output[tem.I_AGR][7]<<" agl: "<<output[tem.I_AGL][7]<<" bgr: "<<output[tem.I_BGR][7]<<" bgl: "<<output[tem.I_BGL][7]<<" CWD: "<<output[tem.I_CWD][7]<<" standead: "<<output[tem.I_DEADWOODC][7]<<" forestage: "<< output[tem.I_AGE][11]<<endl;
	//if (dyr>501) exit(-1);
	//}

	++dyr;
    ++tem.totyr;


    if( dyr >= tem.strteq && 0 == tem.endeq )
    {
      tem.endeq = tem.testEquilibrium();
      //if (tem.veg.cmnt==7 && cohort[pichrt].currentveg==50 && dyr>500) cout <<" endeq: "<<tem.endeq<<" ctrh: "<<tem.veg.yrltrc - tem.microbe.yrdecomp <<" ctdecom: "<<tem.microbe.yrdecomp - tem.microbe.yrrh - tem.microbe.yrDOMprod.carbon<<endl;
    }
  }
  //cout <<" totyr: "<<tem.totyr<<" endeq: " <<tem.endeq<<" initflag: " <<tem.initFlag<<endl;

  if( tem.totyr >= tem.runsize && tem.endeq < 2 ) 
  { 
    tem.nattempt += 1; 
    tem.initFlag = 0;
  }
  else { tem.initFlag = 1; }
  return tem.nattempt;

};
/* *************************************************************
************************************************************* */



/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt60::initializeCohortTEMState( const int& pichrt )
{
  int dm;
  int dnode;
  int i;

  for( i = 0; i < MAXSTATE; ++i )
  {
    cohort[pichrt].y[i] = ZERO;
    cohort[pichrt].prevy[i] = ZERO;
  }

  cohort[pichrt].aggrowdd = ZERO;

  cohort[pichrt].agkd = ZERO;

  cohort[pichrt].c2n = ZERO;

  cohort[pichrt].cneven = ZERO;
  
  cohort[pichrt].convrtflx.carbon = ZERO;
  cohort[pichrt].convrtflx.nitrogen = ZERO;

  cohort[pichrt].cropprveetmx = ZERO;

  cohort[pichrt].cropprvleafmx = ZERO;

  cohort[pichrt].cropprvpetmx = ZERO;

  cohort[pichrt].cropResidue.carbon = ZERO;
  cohort[pichrt].cropResidue.nitrogen = ZERO;

  cohort[pichrt].croptopt = ZERO;

  cohort[pichrt].distmnthcnt = 0;

  cohort[pichrt].disturbflag = 0;
  
  cohort[pichrt].disturbmonth = 0;
  
  for( dm = 0; dm < CYCLE; ++dm )
  {
    cohort[pichrt].dst10[dm] = ZERO;
  }

  cohort[pichrt].eetmx = ZERO;

  cohort[pichrt].firemnthcnt = 0;
  
  cohort[pichrt].firendep = ZERO;

  cohort[pichrt].formPROD10.carbon = ZERO;
  cohort[pichrt].formPROD10.nitrogen = ZERO;

  cohort[pichrt].formPROD100.carbon = ZERO;
  cohort[pichrt].formPROD100.nitrogen = ZERO;
  
  cohort[pichrt].fprevozone = ZERO;

  for( dm = 0; dm < CYCLE; ++dm )
  {
    cohort[pichrt].initPROD1[dm].carbon = ZERO;
    cohort[pichrt].initPROD1[dm].nitrogen = ZERO;
  }

  for( i = 0; i < 10; ++i )
  {
    cohort[pichrt].initPROD10[i].carbon = ZERO;
    cohort[pichrt].initPROD10[i].nitrogen = ZERO;
  }
    
  for( i = 0; i < 100; ++i )
  {
    cohort[pichrt].initPROD100[i].carbon = ZERO;
    cohort[pichrt].initPROD100[i].nitrogen = ZERO;
  }
  
  cohort[pichrt].kd = ZERO;

  cohort[pichrt].natprveetmx = ZERO;

  cohort[pichrt].natprvleafmx = ZERO;

  cohort[pichrt].natprvpetmx = ZERO;

  cohort[pichrt].natseedC = ZERO;

  cohort[pichrt].natseedSTRN = ZERO;

  cohort[pichrt].natseedSTON = ZERO;

  cohort[pichrt].natsoil = ZERO;

  cohort[pichrt].nattopt = ZERO;

  cohort[pichrt].natyreet = ZERO;

  cohort[pichrt].natyrpet= ZERO;

  cohort[pichrt].newleafmx = ZERO;

  cohort[pichrt].newtopt = ZERO;

  cohort[pichrt].nonsolc = ZERO;

  cohort[pichrt].nonsoln = ZERO;

  cohort[pichrt].nretent = ZERO;

  cohort[pichrt].nsretent = ZERO;

  cohort[pichrt].nvretent = ZERO;

  cohort[pichrt].petmx = ZERO;

  cohort[pichrt].prev2tair = ZERO;

  cohort[pichrt].prevco2 = ZERO;

  cohort[pichrt].prevCropResidue.carbon = ZERO;
  cohort[pichrt].prevCropResidue.nitrogen = ZERO;

  cohort[pichrt].prevdst10 = ZERO;

  cohort[pichrt].prevPROD1.carbon = ZERO;
  cohort[pichrt].prevPROD1.nitrogen = ZERO;

  cohort[pichrt].prevPROD10.carbon = ZERO;
  cohort[pichrt].prevPROD10.nitrogen = ZERO;

  cohort[pichrt].prevPROD100.carbon = ZERO;
  cohort[pichrt].prevPROD100.nitrogen = ZERO;

  cohort[pichrt].prevspack = ZERO;

  cohort[pichrt].prevtair = ZERO;
  cohort[pichrt].prevrain = ZERO;
  cohort[pichrt].prevunrmleaf = ZERO;

  cohort[pichrt].productYear = 0;
  
  cohort[pichrt].prvcropnpp = ZERO;

  cohort[pichrt].prveetmx = ZERO;

  cohort[pichrt].prvleafmx = ZERO;

  cohort[pichrt].prvpetmx = ZERO;

  cohort[pichrt].qc = 0;

  cohort[pichrt].sconvrtflx.carbon = ZERO;
  cohort[pichrt].sconvrtflx.nitrogen = ZERO;

  cohort[pichrt].slash.carbon = ZERO;
  cohort[pichrt].slash.nitrogen = ZERO;
  cohort[pichrt].deadwood.carbon = ZERO;
  cohort[pichrt].deadwood.nitrogen = ZERO;

  cohort[pichrt].STMis9 = 0;
  
  cohort[pichrt].STMsmass9 = ZERO;
  
  for( dnode = 0; dnode < MAXNODES; ++dnode )
  {
    cohort[pichrt].STMdx9[dnode] = ZERO;
     
    cohort[pichrt].STMt9[dnode] = ZERO;
    
    cohort[pichrt].STMwater9[dnode] = ZERO;
    
    cohort[pichrt].STMx9[dnode] = ZERO;

    cohort[pichrt].STMxfa9[dnode] = ZERO;

    cohort[pichrt].STMxfb9[dnode] = ZERO;
  }

  for( dnode = 0; dnode < MAXSNODES; ++dnode )
  {
    cohort[pichrt].STMweight9[dnode] = MISSING;
  }

  cohort[pichrt].topt = ZERO;

  cohort[pichrt].tqc = 0;

  cohort[pichrt].vconvrtflx.carbon = ZERO;
  cohort[pichrt].vconvrtflx.nitrogen = ZERO;

  cohort[pichrt].yrltrc = ZERO;
  cohort[pichrt].yrltrn = ZERO;
  cohort[pichrt].leafmortpar = ZERO;
  cohort[pichrt].rootmortpar = ZERO;
  cohort[pichrt].stemmortpar = ZERO;
  cohort[pichrt].leafslash = ZERO;
  cohort[pichrt].leafconv = ZERO;
  cohort[pichrt].rootslash = ZERO;
  cohort[pichrt].rootconv = ZERO;
  cohort[pichrt].stemslash = ZERO;
  cohort[pichrt].stemconv = ZERO;
  cohort[pichrt].standdead = ZERO;
  cohort[pichrt].deadconv = ZERO;
  cohort[pichrt].deadslash = ZERO;
  cohort[pichrt].sconvert = ZERO;
  cohort[pichrt].prod10par = ZERO;
  cohort[pichrt].prod100par = ZERO;


};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TEMelmnt60::outputTEMmonth( const int& pdm )
{

  // Ecosystem carbon pools determined in integrator

  output[tem.I_VEGC][pdm] = tem.getY( tem.I_VEGC );

  output[tem.I_SOLC][pdm] = tem.getY( tem.I_SOLC );

  output[tem.I_DOC][pdm] = tem.getY( tem.I_DOC );

//  output[tem.I_CO2G][pdm] = tem.getY( tem.I_CO2G );

//  output[tem.I_CO2W][pdm] = tem.getY( tem.I_CO2W );

//  output[tem.I_HCO3][pdm] = tem.getY( tem.I_HCO3 );

//  output[tem.I_RHCO3][pdm] = tem.getY( tem.I_RHCO3 );

//  output[tem.I_ALK][pdm] = tem.soil.alkalinity;

  // Ecosystem nitrogen pools determined in integrator

  output[tem.I_STRN][pdm] = tem.getY( tem.I_STRN );

  output[tem.I_STON][pdm] = tem.getY( tem.I_STON );

  output[tem.I_SOLN][pdm] = tem.getY( tem.I_SOLN );

  output[tem.I_DON][pdm] = tem.getY( tem.I_DON );

  output[tem.I_NH4][pdm] = tem.getY( tem.I_NH4 );

  output[tem.I_NO3][pdm] = tem.getY( tem.I_NO3 );
  output[tem.I_DEADWOODC][pdm] = tem.getY( tem.I_DEADWOODC );
  output[tem.I_DEADWOODN][pdm] = tem.getY( tem.I_DEADWOODN );
  output[tem.I_FOLC][pdm] = tem.getY( tem.I_FOLC);
  output[tem.I_STEMC][pdm] = tem.getY( tem.I_STEMC);
  output[tem.I_CROOTC][pdm] = tem.getY( tem.I_CROOTC);
  output[tem.I_FROOTC][pdm] = tem.getY( tem.I_FROOTC);
  output[tem.I_AGR][pdm] = tem.getY( tem.I_AGR);
  output[tem.I_AGL][pdm] = tem.getY( tem.I_AGL);
  output[tem.I_BGR][pdm] = tem.getY( tem.I_BGR);
  output[tem.I_BGL][pdm] = tem.getY( tem.I_BGL);
  output[tem.I_CWD][pdm] = tem.getY( tem.I_CWD);
  output[tem.I_MDC][pdm] = tem.getY( tem.I_MDC);
  //output[tem.I_AGE][pdm] = tem.ag.getage();
  output[tem.I_AGE][pdm] = tem.getY(tem.I_AGE);
  output[tem.I_SOC][pdm] = tem.getY(tem.I_SOC);
  // Ecosystem water pools determined in integrator

  output[tem.I_AVLW][pdm] = (tem.getY( tem.I_SM )
                            * tem.soil.getACTLAYER() / tem.soil.getROOTZ())
                            - tem.soil.getWILTPT();

  if( output[tem.I_AVLW][pdm] < ZERO )
  {
    output[tem.I_AVLW][pdm] = ZERO;
  }

  output[tem.I_SM][pdm] = tem.getY( tem.I_SM );

  output[tem.I_VSM][pdm] = tem.getY( tem.I_SM )
                           / (tem.soil.getROOTZ() * 1000.0);

  output[tem.I_PCTP][pdm] = (tem.getY( tem.I_SM )
                            * tem.soil.getACTLAYER() / tem.soil.getROOTZ())
                            / tem.soil.getTOTPOR();

  output[tem.I_RGRW][pdm] = tem.getY( tem.I_RGRW );

  output[tem.I_SGRW][pdm] = tem.getY(tem.I_SGRW );

  // Monthly phenology determined in integrator

  output[tem.I_UNRMLF][pdm] = tem.getY( tem.I_UNRMLF );

  output[tem.I_LEAF][pdm] = tem.getY( tem.I_LEAF );

  output[tem.I_LAI][pdm] = tem.getY( tem.I_LAI );

  output[tem.I_FPC][pdm] = tem.getY( tem.I_FPC );

  // Monthly carbon fluxes in ecosystems determined in integrator

  output[tem.I_INGPP][pdm] = tem.getY( tem.I_INGPP );

  output[tem.I_GPP][pdm] = tem.getY( tem.I_GPP );
  //if (output[tem.I_GPP][pdm]>0.0) cout <<"mon: "<<pdm<<" outputgpp: "<<output[tem.I_GPP][pdm]<<endl;

  output[tem.I_FOZONE][pdm] = tem.getY( tem.I_FOZONE );

  output[tem.I_FINDOZONE][pdm] = tem.getY( tem.I_FINDOZONE );

  output[tem.I_INNPP][pdm] = tem.getY( tem.I_INNPP );

  output[tem.I_NPP][pdm] = tem.getY( tem.I_NPP );

  output[tem.I_GPR][pdm] = tem.getY( tem.I_GPR );

  output[tem.I_RVMNT][pdm] = tem.getY( tem.I_RVMNT );

  output[tem.I_RVGRW][pdm] = tem.getY( tem.I_RVGRW );

  output[tem.I_ABVGPR][pdm] = tem.getY( tem.I_ABVGPR );

  output[tem.I_ROOTGPR][pdm] = tem.getY( tem.I_ROOTGPR );

//  output[tem.I_ISOPREN][pdm] = tem.getY( tem.I_ISOPREN );

//  output[tem.I_TERPEN][pdm] = tem.getY( tem.I_TERPEN );

//  output[tem.I_ORVOC][pdm] = tem.getY( tem.I_ORVOC );

//  output[tem.I_OVOC][pdm] = tem.getY( tem.I_OVOC );

//  output[tem.I_VOC][pdm] = tem.veg.voc.total;

  output[tem.I_LTRC][pdm] = tem.getY( tem.I_LTRC );

  output[tem.I_CDCMP][pdm] = tem.getY( tem.I_CDCMP );

  output[tem.I_RH][pdm] = tem.getY( tem.I_RH );

  output[tem.I_DOCP][pdm] = tem.getY( tem.I_DOCP );

  output[tem.I_LCHDOC][pdm] = tem.getY( tem.I_LCHDOC );

  output[tem.I_ERDPOC][pdm] = tem.getY( tem.I_ERDPOC );

  // Monthly nitrogen fluxes in ecosystems determined in
  //   integrator

  output[tem.I_AGFRTN][pdm] = tem.getY( tem.I_AGFRTN );

  output[tem.I_BNFIX][pdm] = tem.getY( tem.I_BNFIX );

  output[tem.I_SNFIX][pdm] = tem.getY( tem.I_SNFIX );

  output[tem.I_ANFIX][pdm] = tem.getY( tem.I_ANFIX );

  output[tem.I_INNUP][pdm] = tem.getY( tem.I_INNUP );

  output[tem.I_INNH4UP][pdm] = tem.getY( tem.I_INNH4UP );

  output[tem.I_INNO3UP][pdm] = tem.getY( tem.I_INNO3UP );

  output[tem.I_VNUP][pdm] = tem.getY( tem.I_VNUP );

  output[tem.I_VNH4UP][pdm] = tem.getY( tem.I_VNH4UP );

  output[tem.I_VNO3UP][pdm] = tem.getY( tem.I_VNO3UP );

  output[tem.I_VSUP][pdm] = tem.getY( tem.I_VSUP );

  output[tem.I_VLUP][pdm] = tem.getY( tem.I_VLUP );

  output[tem.I_VNMBL][pdm] = tem.getY( tem.I_VNMBL );

  output[tem.I_VNRSRB][pdm] = tem.getY( tem.I_VNRSRB );

  output[tem.I_LTRN][pdm] = tem.getY( tem.I_LTRN );

  output[tem.I_NDCMP][pdm] = tem.getY( tem.I_NDCMP );

  output[tem.I_DONP][pdm] = tem.getY( tem.I_DONP );

  output[tem.I_GMIN][pdm] = tem.getY( tem.I_GMIN );

  output[tem.I_NH4IMM][pdm] = tem.getY( tem.I_NH4IMM );

  output[tem.I_NIMM][pdm] = tem.getY( tem.I_NIMM );

  output[tem.I_NMIN][pdm] = tem.getY( tem.I_NMIN );

  output[tem.I_AIMMNH4][pdm] = tem.getY( tem.I_AIMMNH4 );

  output[tem.I_AIMMNO3][pdm] = tem.getY( tem.I_AIMMNO3 );

  output[tem.I_AMMN][pdm] = tem.getY( tem.I_AMMN );

  output[tem.I_NTRF][pdm] = tem.getY( tem.I_NTRF );

  output[tem.I_NO3P][pdm] = tem.getY( tem.I_NO3P );

  output[tem.I_NOP][pdm] = tem.getY( tem.I_NOP );

  output[tem.I_N2OP][pdm] = tem.getY( tem.I_N2OP );

  output[tem.I_N2P][pdm] = tem.getY( tem.I_N2P );

  output[tem.I_DNTRF][pdm] = tem.getY( tem.I_DNTRF );

  output[tem.I_NH3FLX][pdm] = tem.getY( tem.I_NH3FLX );

  output[tem.I_NOFLX][pdm] = tem.getY( tem.I_NOFLX );

  output[tem.I_N2OFLX][pdm] = tem.getY( tem.I_N2OFLX );

  output[tem.I_N2FLX][pdm] = tem.getY( tem.I_N2FLX );

  output[tem.I_LCHNO3][pdm] = tem.getY( tem.I_LCHNO3 );

  output[tem.I_LCHDON][pdm] = tem.getY( tem.I_LCHDON );

  output[tem.I_ERDPON][pdm] = tem.getY( tem.I_ERDPON );

  // Monthly water fluxes in ecosystems

  output[tem.I_AGIRRIG][pdm] = tem.getY( tem.I_AGIRRIG );

  output[tem.I_INEET][pdm] = tem.getY( tem.I_INEET );

  output[tem.I_EET][pdm] = tem.getY( tem.I_EET );

  output[tem.I_RPERC][pdm] = tem.getY( tem.I_RPERC );

  output[tem.I_SPERC][pdm] = tem.getY( tem.I_SPERC );

  output[tem.I_RRUN][pdm] = tem.getY( tem.I_RRUN );

  output[tem.I_SRUN][pdm] = tem.getY( tem.I_SRUN );

  // Other ecosystem carbon pools

  output[tem.I_NSOLC][pdm] = tem.soil.getNSOLC();

  output[tem.I_TSOLC][pdm] = tem.soil.getTSOLC();

  output[tem.I_TOTEC][pdm] = tem.ag.getTOTEC();

  output[tem.I_TOTC][pdm] = tem.getTOTALC();

  // Other ecosystem nitrogen pools

  output[tem.I_VEGN][pdm] = tem.veg.getVEGN();

  output[tem.I_NSOLN][pdm] = tem.soil.getNSOLN();

  output[tem.I_TSOLN][pdm] = tem.soil.getTSOLN();

  output[tem.I_AVLN][pdm] = tem.soil.getAVLN();


  // Other ecosystem water pools
  
  output[tem.I_SNWPCK][pdm] = tem.soil.getSNOWPACK();

  // Other monthly carbon fluxes in ecosystems

  output[tem.I_RSOIL][pdm] = tem.getRSOIL();

  output[tem.I_NEP][pdm] = tem.getNEP();

  output[tem.I_NCE][pdm] = tem.getNCE();

  output[tem.I_NTCB][pdm] = tem.getNTCB();

  // Other monthly nitrogen fluxes in ecosystems

  output[tem.I_NINP][pdm] = tem.soil.getNINPUT();

  output[tem.I_NLST][pdm] = tem.soil.getNLOST();

  output[tem.I_NTNB][pdm] = tem.getNTNB();

  // Other monthly water fluxes in ecosystems

  output[tem.I_PET][pdm] = tem.atms.getPET();

  output[tem.I_SNWINF][pdm] = tem.soil.getSNOWINF();

  output[tem.I_WYLD][pdm] = tem.soil.getH2OYLD();

  // Carbon in Human product pools

  output[tem.I_AGPRDC][pdm] = tem.ag.getPROD1C();

  output[tem.I_PROD10C][pdm] = tem.ag.getPROD10C();

  output[tem.I_PROD100C][pdm] = tem.ag.getPROD100C();

  output[tem.I_TOTPRDC][pdm] = tem.ag.getTOTPRODC();

  // Carbon in crop residue pool

  output[tem.I_RESIDC][pdm] = tem.ag.getCROPRESIDUEC();

  output[tem.I_AGSTUBC][pdm] = tem.ag.getSTUBBLEC();

  // Nitrogen in Human product pools

  output[tem.I_AGPRDN][pdm] = tem.ag.getPROD1N();

  output[tem.I_PROD10N][pdm] = tem.ag.getPROD10N();

  output[tem.I_PROD100N][pdm] = tem.ag.getPROD100N();

  output[tem.I_TOTPRDN][pdm] = tem.ag.getTOTPRODN();

  // Nitrogen in crop residue pool

  output[tem.I_RESIDN][pdm] = tem.ag.getCROPRESIDUEN();

  output[tem.I_AGSTUBN][pdm] = tem.ag.getSTUBBLEN();

  // Monthly carbon fluxes associated with
  //  agricultural conversion

  output[tem.I_CNVRTC][pdm] = tem.ag.getCONVRTFLXC();

  output[tem.I_VCNVRTC][pdm] = tem.ag.getVCONVRTFLXC();

  output[tem.I_SCNVRTC][pdm] = tem.ag.getSCONVRTFLXC();
  //if (tem.ag.getSCONVRTFLXC()>0.0) cout <<"sconvertc: "<<output[tem.I_SCNVRTC][pdm]<<endl;

  output[tem.I_SLASHC][pdm] = tem.ag.getSLASHC();

  output[tem.I_CFLX][pdm] = tem.ag.getCFLUX();
  output[tem.I_FFLC][pdm] = tem.ag.getfflc();
  output[tem.I_FFDC][pdm] = tem.ag.getffdc();
  output[tem.I_SWFC][pdm] = tem.ag.getswfc();
  output[tem.I_DWFC][pdm] = tem.ag.getdwfc();

  // Monthly nitrogen fluxes associated with
  //  agricultural conversion

  output[tem.I_CNVRTN][pdm] = tem.ag.getCONVRTFLXN();

  output[tem.I_VCNVRTN][pdm] = tem.ag.getVCONVRTFLXN();

  output[tem.I_SCNVRTN][pdm] = tem.ag.getSCONVRTFLXN();

  output[tem.I_SLASHN][pdm] = tem.ag.getSLASHN();

  output[tem.I_NRETNT][pdm] = tem.ag.getNRETENT();

  output[tem.I_NVRTNT][pdm] = tem.ag.getNVRETENT();

  output[tem.I_NSRTNT][pdm] = tem.ag.getNSRETENT();

  // Monthly carbon and nitrogen fluxes from agricultural
  //   ecosystems

  output[tem.I_AGFPRDC][pdm] = tem.ag.getCROPPRODC();
  output[tem.I_AGFPRDN][pdm] = tem.ag.getCROPPRODN();

  output[tem.I_FRESIDC][pdm] = tem.ag.getFORMCROPRESIDUEC();
  output[tem.I_FRESIDN][pdm] = tem.ag.getFORMCROPRESIDUEN();

  output[tem.I_AGPRDFC][pdm] = tem.ag.getPROD1DECAYC();
  output[tem.I_AGPRDFN][pdm] = tem.ag.getPROD1DECAYN();

  output[tem.I_RESIDFC][pdm] = tem.ag.getCROPRESIDUEFLXC();
  output[tem.I_RESIDFN][pdm] = tem.ag.getCROPRESIDUEFLXN();


  // Monthly carbon and nitrogen fluxes from products

  output[tem.I_PRDF10C][pdm] = tem.ag.getFORMPROD10C();
  output[tem.I_PRDF10N][pdm] = tem.ag.getFORMPROD10N();

  output[tem.I_PRD10FC][pdm] = tem.ag.getPROD10DECAYC();
  output[tem.I_PRD10FN][pdm] = tem.ag.getPROD10DECAYN();

  output[tem.I_PRDF100C][pdm] = tem.ag.getFORMPROD100C();
  output[tem.I_PRDF100N][pdm] = tem.ag.getFORMPROD100N();

  output[tem.I_PRD100FC][pdm] = tem.ag.getPROD100DECAYC();
  output[tem.I_PRD100FN][pdm] = tem.ag.getPROD100DECAYN();

  output[tem.I_TOTFPRDC][pdm] = tem.ag.getFORMTOTPRODC();
  output[tem.I_TOTFPRDN][pdm] = tem.ag.getFORMTOTPRODN();

  output[tem.I_TOTPRDFC][pdm] = tem.ag.getTOTPRODDECAYC();
  output[tem.I_TOTPRDFN][pdm] = tem.ag.getTOTPRODDECAYN();

  //  Output agricultural area-specific vs natural area-specific
  //    results

  if( 1 == tem.ag.state )
  {
    output[tem.I_CROPC][pdm] = tem.getY( I_VEGC );
    output[tem.I_NATVEGC][pdm] = ZERO;

    output[tem.I_CROPN][pdm] = tem.veg.getVEGN();
    output[tem.I_NATVEGN][pdm] = ZERO;

    output[tem.I_CSTRN][pdm] = tem.getY( I_STRN );
    output[tem.I_NATSTRN][pdm] = ZERO;

    output[tem.I_CSTON][pdm] = tem.getY( I_STON );
    output[tem.I_NATSTON][pdm] = ZERO;

    output[tem.I_CROPULF][pdm] = tem.getY( I_UNRMLF );
    output[tem.I_NATULF][pdm] = ZERO;

    output[tem.I_CROPLEAF][pdm] = tem.getY( I_LEAF );
    output[tem.I_NATLEAF][pdm] = ZERO;

    output[tem.I_CROPLAI][pdm] = tem.getY( I_LAI );
    output[tem.I_NATLAI][pdm] = ZERO;

    output[tem.I_CROPFPC][pdm] = tem.getY( I_FPC );
    output[tem.I_NATFPC][pdm] = ZERO;

    output[tem.I_AGINGPP][pdm] = tem.getY( I_INGPP );
    output[tem.I_NATINGPP][pdm] = ZERO;

    output[tem.I_AGGPP][pdm] = tem.getY( I_GPP );
    output[tem.I_NATGPP][pdm] = ZERO;

    output[tem.I_AGINNPP][pdm] = tem.getY( I_INNPP );
    output[tem.I_NATINNPP][pdm] = ZERO;

    output[tem.I_AGNPP][pdm] = tem.getY( I_NPP );
    output[tem.I_NATNPP][pdm] = ZERO;

    output[tem.I_AGGPR][pdm] = tem.getY( I_GPR );
    output[tem.I_NATGPR][pdm] = ZERO;

    output[tem.I_AGRVMNT][pdm] = tem.getY( I_RVMNT );
    output[tem.I_NATRVMNT][pdm] = ZERO;

    output[tem.I_AGRVGRW][pdm] = tem.getY( I_RVGRW );
    output[tem.I_NATRVGRW][pdm] = ZERO;

    output[tem.I_AGLTRC][pdm] = tem.getY( I_LTRC );
    output[tem.I_NATLTRC][pdm] = ZERO;

    output[tem.I_AGSNFX][pdm] = tem.getY( I_SNFIX );
    output[tem.I_NATSNFX][pdm] = ZERO;

    output[tem.I_AGINNUP][pdm] = tem.getY( I_INNUP );
    output[tem.I_NATINNUP][pdm] = ZERO;

    output[tem.I_AINNH4UP][pdm] = tem.getY( I_INNH4UP );
    output[tem.I_NINNH4UP][pdm] = ZERO;

    output[tem.I_AINNO3UP][pdm] = tem.getY( I_INNO3UP );
    output[tem.I_NINNO3UP][pdm] = ZERO;

    output[tem.I_AGVNUP][pdm] = tem.getY( I_VNUP );
    output[tem.I_NATVNUP][pdm] = ZERO;

    output[tem.I_AGVNH4UP][pdm] = tem.getY( I_VNH4UP );
    output[tem.I_NVNH4UP][pdm] = ZERO;

    output[tem.I_AGVNO3UP][pdm] = tem.getY( I_VNO3UP );
    output[tem.I_NVNO3UP][pdm] = ZERO;

    output[tem.I_AGVSUP][pdm] = tem.getY( I_VSUP );
    output[tem.I_NATVSUP][pdm] = ZERO;

    output[tem.I_AGVLUP][pdm] = tem.getY( I_VLUP );
    output[tem.I_NATVLUP][pdm] = ZERO;

    output[tem.I_AGVNMBL][pdm] = tem.getY( I_VNMBL );
    output[tem.I_NATVNMBL][pdm] = ZERO;

    output[tem.I_AGVNRSRB][pdm] = tem.getY( I_VNRSRB );
    output[tem.I_NVNRSRB][pdm] = ZERO;

    output[tem.I_AGLTRN][pdm] = tem.getY( I_LTRN );
    output[tem.I_NATLTRN][pdm] = ZERO;
  }
  else
  {
    output[tem.I_CROPC][pdm] = ZERO;
    output[tem.I_NATVEGC][pdm] = tem.getY( I_VEGC );

    output[tem.I_CROPN][pdm] = ZERO;
    output[tem.I_NATVEGN][pdm] = tem.veg.getVEGN();

    output[tem.I_CSTRN][pdm] = ZERO;
    output[tem.I_NATSTRN][pdm] = tem.getY( I_STRN );

    output[tem.I_CSTON][pdm] = ZERO;
    output[tem.I_NATSTON][pdm] = tem.getY( I_STON );

    output[tem.I_CROPULF][pdm] = ZERO;
    output[tem.I_NATULF][pdm] = tem.getY( I_UNRMLF );

    output[tem.I_CROPLEAF][pdm] = ZERO;
    output[tem.I_NATLEAF][pdm] = tem.getY( I_LEAF );

    output[tem.I_CROPLAI][pdm] = ZERO;
    output[tem.I_NATLAI][pdm] = tem.getY( I_LAI );

    output[tem.I_CROPFPC][pdm] = ZERO;
    output[tem.I_NATFPC][pdm] = tem.getY( I_FPC );

    output[tem.I_AGINGPP][pdm] = ZERO;
    output[tem.I_NATINGPP][pdm] = tem.getY( I_INGPP );

    output[tem.I_AGGPP][pdm] = ZERO;
    output[tem.I_NATGPP][pdm] = tem.getY( I_GPP );

    output[tem.I_AGINNPP][pdm] = ZERO;
    output[tem.I_NATINNPP][pdm] = tem.getY( I_INNPP );

    output[tem.I_AGNPP][pdm] = ZERO;
    output[tem.I_NATNPP][pdm] = tem.getY( I_NPP );

    output[tem.I_AGGPR][pdm] = ZERO;
    output[tem.I_NATGPR][pdm] = tem.getY( I_GPR );

    output[tem.I_AGRVMNT][pdm] = ZERO;
    output[tem.I_NATRVMNT][pdm] = tem.getY( I_RVMNT );

    output[tem.I_AGRVGRW][pdm] = ZERO;
    output[tem.I_NATRVGRW][pdm] = tem.getY( I_RVGRW );

    output[tem.I_AGLTRC][pdm] = ZERO;
    output[tem.I_NATLTRC][pdm] = tem.getY( I_RVGRW );

    output[tem.I_AGSNFX][pdm] = ZERO;
    output[tem.I_NATSNFX][pdm] = tem.getY( I_SNFIX );

    output[tem.I_AGINNUP][pdm] = ZERO;
    output[tem.I_NATINNUP][pdm] = tem.getY( I_INNUP );

    output[tem.I_AINNH4UP][pdm] = ZERO;
    output[tem.I_NINNH4UP][pdm] = tem.getY( I_INNH4UP );

    output[tem.I_AINNO3UP][pdm] = ZERO;
    output[tem.I_NINNO3UP][pdm] = tem.getY( I_INNO3UP );

    output[tem.I_AGVNUP][pdm] = ZERO;
    output[tem.I_NATVNUP][pdm] = tem.getY( I_VNUP );

    output[tem.I_AGVNH4UP][pdm] = ZERO;
    output[tem.I_NVNH4UP][pdm] = tem.getY( I_VNH4UP );

    output[tem.I_AGVNO3UP][pdm] = ZERO;
    output[tem.I_NVNO3UP][pdm] = tem.getY( I_VNO3UP );

    output[tem.I_AGVSUP][pdm] = ZERO;
    output[tem.I_NATVSUP][pdm] = tem.getY( I_VSUP );

    output[tem.I_AGVLUP][pdm] = ZERO;
    output[tem.I_NATVLUP][pdm] = tem.getY( I_VLUP );

    output[tem.I_AGVNMBL][pdm] = ZERO;
    output[tem.I_NATVNMBL][pdm] = tem.getY( I_VNMBL );

    output[tem.I_AGVNRSRB][pdm] = ZERO;
    output[tem.I_NVNRSRB][pdm] = tem.getY( I_VNRSRB );

    output[tem.I_AGLTRN][pdm] = ZERO;
    output[tem.I_NATLTRN][pdm] = tem.getY( I_LTRN );
  }


  // Monthly soil thermal dynamics

  output[tem.I_TSOIL][pdm] = tem.soil.getTSOIL();

  output[tem.I_DST0][pdm] = tem.soil.stm.getDST0();

  output[tem.I_DST5][pdm] = tem.soil.stm.getDST5();

  output[tem.I_DST10][pdm] = tem.soil.getDST10();

  output[tem.I_DST20][pdm] = tem.soil.stm.getDST20();

  output[tem.I_DST50][pdm] = tem.soil.stm.getDST50();

  output[tem.I_DST100][pdm] = tem.soil.stm.getDST100();

  output[tem.I_DST200][pdm] = tem.soil.stm.getDST200();
  output[tem.I_DST300][pdm] = tem.soil.stm.getDST300();

  output[tem.I_FRONTD][pdm] = tem.soil.stm.getFRONTD();

  output[tem.I_THAWBE][pdm] = tem.soil.stm.getTHAWBEGIN1();

  output[tem.I_THAWEND][pdm] = tem.soil.stm.getTHAWEND1();

  output[tem.I_THAWPCT][pdm] = tem.veg.getTHAWPCT();

  output[tem.I_ACTLAYER][pdm] = tem.soil.getACTLAYER();

};



/* **************************************************************
************************************************************** */

int TEMelmnt60::setGIStopography( ofstream& rflog1,
                                  int& ftlerr,
                                  FILE* fstxt,
                                  FILE* felev )
{

  int gisend;

  //Soildata60 fao;  //fmy: as public so that can be called from main - Aug 2013
  //Elevdata60 elv;  //fmy: as public so that can be called from main - Aug 2013

  gisend = fao.getdel( fstxt );

  if( -1 == gisend )
  {
    rflog1 << "Ran out of Soil texture data" << endl << endl;

    exit( -1 );
  }

  ftlerr = coregerr( rflog1,
                     "Climate",
                     col,
                     row,
                     "TEXTURE",
                     fao.col,
                     fao.row );

  tem.soil.setPCTSILT( fao.pctsilt );
  tem.soil.setPCTCLAY( fao.pctclay );

  gisend = elv.getdel( felev );

  if( gisend == -1 )
  {
    rflog1 << "Ran out of Elevation data" << endl << endl;

    exit( -1 );
  }

  ftlerr = coregerr( rflog1,
                     "Climate",
                     col,
                     row,
                     "ELEV",
                     elv.col,
                     elv.row );

  tem.elev = elv.elev;

  return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TEMelmnt60::setTEMequilState( ofstream& rflog1,
                                   const int& equil,
                                   const int& totsptime,
                                   const int& pichrt )
{

  int dm;

  const int dyr = 0;


  // Set all TEM state-related variables in cohort to MISSING
  //   (i.e. start with "clean slate" in cohort)

  initializeCohortTEMState( pichrt );

  // Assign cohort data to TEM (i.e. transfer information from
  //   the land cover/land use module to TEM and start with
  //   "clean slate" for TEM cohort data)

  dm = 0;
  getTEMCohortState( pichrt, dm );
  //cout <<" y[3]1:" <<cohort[pichrt].y[3]<<endl;

  cohort[pichrt].qc = ACCEPT;
  cohort[pichrt].tqc = ACCEPT;

  tem.qualcon[dyr] = 0;

  tem.totyr = 0;

//  tem.ag.fert1950flag = 0;

  //tem.atms.setMXTAIR( mxtair ); //moved by cgs to below to get the monthly maxtemp
  tem.atms.yrprec = yrprec;

  // Check TEM climate input for valid data

  for( dm = 0; dm < CYCLE; ++dm )
  {
    tem.atms.setGIRR( climate[clm.I_GIRR][dm] );
    tem.atms.setCLDS( climate[clm.I_CLDS][dm] );
    tem.atms.setNIRR( climate[clm.I_NIRR][dm] );
    tem.atms.setPAR(  climate[clm.I_PAR][dm] );
    tem.atms.setTAIR( climate[clm.I_TAIR][dm] );
    tem.atms.setRAIN( climate[clm.I_RAIN][dm] );
    tem.atms.setPREC( climate[clm.I_PREC][dm] );
    tem.atms.setSNOWFALL( climate[clm.I_SNWFAL][dm] );
    tem.atms.setCO2( climate[clm.I_CO2][dm] );
    tem.atms.setAOT40( climate[clm.I_AOT40][dm] );
    tem.atms.setTOTNDEP( climate[clm.I_TNDEP][dm] );
    tem.atms.setNH4DEP( climate[clm.I_NH4DEP][dm] );
    tem.atms.setNO3DEP( climate[clm.I_NO3DEP][dm] );
    tem.atms.setMXTAIR( mxtair );

    cohort[pichrt].qc = temgisqc( cohort[pichrt].chrtarea,
                                  tem.soil.getPCTSILT(),
                                  tem.soil.getPCTCLAY(),
                                  tem.veg.cmnt,
                                  tem.elev,
                                  tem.atms.getNIRR(),
                                  tem.atms.getPAR(),
                                  tem.atms.getTAIR(),
                                  tem.atms.getMXTAIR(),
                                  avetair,
                                  tem.atms.yrprec,
                                  tem.atms.getRAIN(),
                                  tem.atms.getSNOWFALL(),
                                  tem.atms.getCO2(),
                                  tem.atms.getAOT40(),
                                  tem.atms.getNDEP() );

    //if (cohort[pichrt].qc >0) cout <<"qc: "<<cohort[pichrt].qc<<endl;
    if( cohort[pichrt].qc != ACCEPT )
    {
      rflog1 << "temgisqc = " << cohort[pichrt].qc;
      rflog1 << " during month " << (dm+1) << endl;
      //break;
      //cout <<" input data have errors: "<<" qc: "<<cohort[pichrt].qc<<" during month " << (dm+1)<<endl;
      //exit(-1);
    }

	//cout <<" cmnt: " <<tem.veg.cmnt<<"tair: " << tem.atms.getTAIR() <<" co2: " <<tem.atms.getCO2()<<" ndep: " <<tem.atms.getNDEP()<<" texture:  "<<tem.soil.getPCTCLAY()<<endl;

	// Determine initial values for tem.atms.prvpetmx,
    //   tem.atms.prveetmx and and tem.veg.topt based on
    //   long-term mean climate

    tem.setEquilEvap( tem.atms.getNIRR(),
                      tem.atms.getTAIR(),
                      dm );


    // Set initial value of soil temperature at 10 cm for cohort
    //   to zero for all months (value will later be updated
    //   monthly) for use as tem.soil.nextdst10

    cohort[pichrt].dst10[dm] = avetair;
  }
  //cout <<"qc0: "<<cohort[pichrt].qc<<" nupnh41a: " << tem.veg.getNUPNH41A(tem.veg.cmnt)<<" cmnt0: "<< tem.veg.cmnt<<endl;

  // Check TEM parameters for specific vegetation types

  if( ACCEPT == cohort[pichrt].qc )
  {

	cohort[pichrt].qc = tem.ecdqc( tem.veg.cmnt );
	if (cohort[pichrt].chrtarea <0.05) cohort[pichrt].qc=198; //added by cgs to assign no-data to skip the grid cells with land area ~0.0 km2

    if( cohort[pichrt].qc != ACCEPT )
    {
      // Note: If a TEM parameter is invalid,
      //   cohort[pichrt].qc will have a value greater than
      //   100
      //cout <<"col: " <<col <<" row: " <<row<<" ecdqc: " <<cohort[pichrt].qc	<<endl;
      rflog1 << "temecdqc = " << cohort[pichrt].qc << endl;
    }
  }
  //cout <<"qc1.0: "<<cohort[pichrt].qc<<" nupnh41a: " << tem.veg.getNUPNH41A(tem.veg.cmnt)<<endl;

  if( cohort[pichrt].qc != ACCEPT )
  {
    // If environmental conditions are too extreme for the
    //   existence of vegetation (e.g., no precipitation or
    //   constant freezing air temperatures), assign zero to
    //   all TEM variables if the plant community is anything
    //   besides ice and open water; and all TEM parameters
    //   are valid (i.e. cohort[pichrt].qc < 100 )

	//modified by cgs. if input data is wrong, output values are MISSING

	if( cohort[pichrt].qc < 100
         && tem.veg.cmnt >= 0
         && (mxtair <= -5.0 || yrprec == ZERO) )
	//if( cohort[pichrt].qc == 9
	        // && tem.veg.cmnt >= 0)

    {
      // Set tqc flag to assign zero to all TEM variables
      //   during simulation

      cohort[pichrt].tqc = TQCZEROFLAG;
    }
    else { cohort[pichrt].tqc = REJECT;}


	//cohort[pichrt].tqc = REJECT;

    // Set missing values to telmnt[0].output

    setTEMmiss( dyr,
                equil,
                totsptime,
                pichrt  );
  }
  else // "cohort[pichrt].qc == ACCEPT"
  {

/* *************************************************************
                   Start Equilibrium Conditions
************************************************************* */

    // Determine soil properties of element based on
    //   soil texture

    tem.soil.xtext( tem.veg.cmnt,
                    tem.soil.getPCTSILT(),
                    tem.soil.getPCTCLAY() );

    //if (tem.veg.cmnt==5) cout <<"silt: "<<tem.soil.getPCTSILT()<<" clay: "<<tem.soil.getPCTCLAY()<<" psi: "<<tem.soil.getPSIPLUSC()<<endl;
    // Initialize tem.atms.prevco2

    tem.atms.setPREVCO2( climate[clm.I_CO2][CYCLE-1] );


    // Initialize TEM parameters based on element's
    //   (i.e. grid cell) vegetation type, soil texture
    //   and atmospheric CO2 concentration

    tem.setELMNTecd( tem.veg.cmnt, tem.soil.getPSIPLUSC() );

    tem.setEquilC2N( tem.veg.cmnt,
                     tem.atms.getPREVCO2() );


    // Assume potential vegetation when determining
    //   equilibrium conditions
    //closed by CGS2014. Open crop and management.
    //tem.ag.state = 0;
    //tem.ag.prvstate = 0;

    //tem.ag.tillflag = 0;
    //tem.ag.fertflag = 0;
    //tem.ag.irrgflag = 0;

    tem.disturbflag = 0;
    tem.distmnthcnt = 0;


   // Initialize agricultural growing degree days to zero

    tem.ag.setGROWDD( ZERO );


    // "While" loop to allow adaptive integrator tolerance
    //   (i.e. tem.tol) to be reduced if chaotic behavior
    //   occurs.  Try up to "tem.maxnrun" times to equilibrate
    //   TEM.  If TEM does not equilibrate within "tem.runsize"
    //   iterations, decrease tem.tol by an order of magnitude
    //   and try again

    tem.nattempt = 0;
    tem.tol = tem.inittol;
    tem.baseline = tem.initbase;
    tem.initFlag = 0;

    while( tem.nattempt < tem.maxnrun
           && 0 == tem.initFlag )
    {
      tem.nattempt = equilibrateTEM( pichrt,
	                             tem.tol,
	                             rflog1 );

      if( tem.nattempt < tem.maxnrun
          && 0 == tem.initFlag )
      {
      	tem.tol /= 10.0;
      }
      //double xx= tem.getY(I_VEGC);

      //if (xx!=xx) cout <<" I_vegC1:" <<xx<<" col: " <<col<< "row: "<<row <<"value of soilc1: "<<is_valid(xx)<< " or: "<<is_nan(xx) <<" or: " <<is_infinite(xx)<<endl;
      //if (xx >0.0 && xx < 0.001) cout <<"vegc: "<<xx<<endl;
    }
    //cout <<"model check2: " <<" initFlag: "<< tem.initFlag<<" totyr2: " << tem.totyr<<"nattempt; "<< tem.nattempt<<endl;
    // Update summary variables for initial agricultural
    //   state of cohort at end of equilibrium portion
    //   of the TEM simulation
    //if (col > -109 && col <-76 && row == 54) cout <<" pichrt2: " <<pichrt<<" maxcohort: " << maxcohorts <<endl;
    tem.ag.setNATSEEDC( ZERO );
    tem.ag.setNATSEEDSTRN( ZERO );
    tem.ag.setNATSEEDSTON( ZERO );
    tem.ag.setCROPPRVLEAFMX( 1.0 );
    tem.ag.setCROPTOPT( tem.veg.getTOPT() );
    tem.ag.setCROPPRVPETMX( tem.atms.getPRVPETMX() );
    tem.ag.setCROPPRVEETMX( tem.soil.getPRVEETMX() );
    tem.ag.setPRVCROPNPP( ZERO );


    // Save quality control information about the simulation
    //   conditions when the equilibrium portion ended
    //   (i.e. did the carbon and nitrogen fluxes really come
    //         to equilibrium or was the run terminated after
    //         running chaotically for a specified maximum
    //         number of years?)

    tem.qualcon[dyr] += (tem.nattempt + 1);


    // If simulation is part of a transient simulation, reset
    //   tem.totyr to represent an actual year rather than
    //   the number of iterations required to reach equilibrum

    if( 0 == equil )
    {
      tem.totyr = clm.modstartyr - totsptime - 1;
      ttotyr[dyr] = tem.totyr;
     //closed by cgs2014. to deal with vegetation carbon = 0 but soil carbon is not 0
     cohort[pichrt].tqc = transqc( tem.maxyears,
	                            tem.totyr,
	                            output[tem.I_VEGC] );
    }
    else { ttotyr[dyr] = tem.totyr; }
  } // End of "cohort.qc == ACCEPT"

  //cout <<"model check3: "<<" totyr: "<<tem.totyr<<endl;
  // Save TEM state of cohort to telmnt[0].cohort
  //cout <<"qc2: "<<cohort[pichrt].qc<<"tqc: " <<cohort[pichrt].tqc<<" missing: " << MISSING<< " pichrt: "<<pichrt<<" cmnt: "<<tem.veg.cmnt<< endl;

  saveTEMCohortState( pichrt );
  //if (col > -109 && col <-76 && row == 54) cout <<" pichrt3: " <<pichrt<<" maxcohort: " << maxcohorts <<" col: " <<col<<" row: " <<row<<endl;


};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt60::setTEMmiss( const int& pdyr,
                             const int& equil,
                             const int& totsptime,
                             const int& pichrt )
{
  int dm;
  int i;

  if( 0 == equil )
  {
    ttotyr[pdyr] = clm.modstartyr
                   - totsptime - 1
                   + (pdyr * tem.diffyr);
  }
  else
  {
    ttotyr[pdyr] = -999;
  }

  tem.totyr = ttotyr[pdyr];
  /*
  if( TQCZEROFLAG == cohort[pichrt].tqc )
  {

    if( 1 == equil ) { ttotyr[pdyr] = 1; }

    // Assign zero to all TEM state variables

    for( i = 0; i < MAXSTATE; ++i )
    {
      tem.setY( ZERO, i );
      tem.setPREVY(ZERO, i );
    }

    for( i = MAXSTATE; i < NUMEQ; ++i )
    {
      tem.setY( ZERO, i );
    }

    // Assign zero to all TEM ouput variables

    for( i = 0; i < NUMTEM; ++i )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        output[i][dm] = ZERO;
      }
    }
  }
  else*/
  {
    // Assign missing values to grid cells that are covered by ice or open
    // water, or where TEM did not converge on a solution

    for( i = 0; i < MAXSTATE; ++i )
    {
      tem.setY( MISSING, i );
      tem.setPREVY( MISSING, i );
    }
    for( i = MAXSTATE; i < NUMEQ; ++i )
    {
      tem.setY( MISSING, i );
    }

    for( i = 0; i < NUMTEM; ++i )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        output[i][dm] = MISSING;
      }
    }
  }

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

int TEMelmnt60::temgisqc( const double& subarea,
                          const double& pctsilt,
                          const double& pctclay,
                          const int& cmnt,
                          const double& elev,
                          const double& nirr,
                          const double& par,
                          const double& tair,
                          const double& mxtair,
                          const double& avtair,
                          const double& yrprec,
                          const double& rain,
                          const double& snowfall,
                          const double& co2,
                          const double& aot40,
                          const InorgN60& ndep )


{
  int qc;

  qc = ACCEPT;

  //if( subarea < 1 ) { return qc = 1; } //closed by cgs to allow subarea =0;
  if( pctsilt < ZERO ) { return qc = 2; }
  if( pctclay < ZERO ) { return qc = 3; }
  if( cmnt < 0 || cmnt > NUMVEG ) { return qc = 4; }
  if( elev <= -999.0 ) { return qc = 5;}

  if( nirr <= -1.0 ) { return qc = 6; }
  if( par <= -1.0 ) { return qc = 7; }
  if( tair <= -99.0 ) { return qc = 8; }
  //modified by cgs. to run the model at grid cells with exteme climate in some years.
  //avoid assigning zero or missing to these extreme years
  if( mxtair <= -50.0 ) { return qc = 9; }
  //if( mxtair <= -5.0 ) { return qc = 9; }
  if( avtair <= -99.0 ) { return qc = 10; }
  if( yrprec < ZERO ) { return qc = 11; }
 // if( yrprec < ZERO ) { return qc = 11; }//modified by cgs to run the model at grid cells with exteme climate in some years
  if( rain <= -1.0 ) { return qc = 12; }
  if( snowfall <= -1.0 ) { return qc = 13; }
  if( co2 <= -1.0 ) { return qc = 14; }
  if( aot40 <= -1.0 ) { return qc = 15; }
  if( ndep.nh4 <= -1.0 ) { return qc = 16; }
  if( ndep.no3 <= -1.0 ) { return qc = 17; }

  return qc;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt60::temwritepred( ofstream fout[NUMTEM],
                               const vector<string>& predname,
                               const int& pdyr,
                               const int& pichrt,
                               const int& ntempred )
{
  // Units conversion from grams to milligrams
  const double GRAMS2MG = 1000.0;

  // Units conversion from meters to millimeters
  const double METERS2MM = 1000.0;

  // Units conversion from proportion to percent
  const double PROP2PCT = 100.0;

  int i;
  int dm;
  Temdata60 tempred;

  //if (tem.getY(tem.I_GPP) > 0.0) cout <<"pdyr: "<<pdyr<<" outputgpp2: "<<output[tem.I_GPP][dm]<< " getgpp: "<<tem.getY(tem.I_GPP)<<" I_GPP: "<<tem.I_GPP<<"i_deadwoodc: "<<tem.I_DEADWOODC<<endl;
  for( i = 0; i < ntempred; ++i )
  {
    // ************** Carbon stocks in ecosystems  *************

    //cout <<"ntempred: "<<ntempred<<" i: "<<i<<"predname.at: "<<predname.at( i )<<" tem.predstr.at: "<<tem.predstr.at(tem.I_GPP)<<endl;
    if( predname.at( i ) == tem.predstr.at( tem.I_VEGC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VEGC][dm];
        //if (dm==2) cout <<" I_VEGC: "<<tempred.mon[dm]<<" "<<"outputnumber: " <<ntempred<<" predstr.I_snwpck: "<< tem.predstr.at( tem.I_SNWPCK )<<" ";
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SOLC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SOLC][dm];
        //double xx = tempred.mon[dm];
        //if (col == -137.625 && row == 59.625) cout <<"I_solc: " <<output[tem.I_SOLC][dm] <<endl;
        //if (output[tem.I_SOLC][dm]!=output[tem.I_SOLC][dm] && dm == 7 && pdyr<=1) cout <<" I_soilC1:" <<tempred.mon[dm]<<" col: " <<col<< "row: "<<row <<"value of soilc1: "<<is_valid(xx)<< " or: "<<is_nan(xx) <<" or: " <<is_infinite(xx)<<endl;

      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DOC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DOC][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_DEADWOODC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DEADWOODC][dm];
        //if (output[tem.I_DEADWOODC][dm]>0.01) cout <<"col: "<< col<<" row: "<<row<<" pdyr: "<<pdyr+1800<<" disturbflag: "<<tem.disturbflag<<" deadwoodc: "<< output[tem.I_DEADWOODC][dm]<<endl;
        //if (col == -128.125 && row == 59.875 && pdyr <= 10) cout <<"pdyr: " <<pdyr + 1800 <<" month: " <<dm<<" cohort: " <<pichrt <<" cmnt: " <<tem.veg.cmnt<<" curveg: " <<tem.veg.getCURRENTVEG() << " I_deadwoodc: "<<output[tem.I_DEADWOODC][dm]<<endl;

      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_CO2G ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CO2G][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CO2W ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CO2W][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_HCO3) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_HCO3][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RHCO3 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RHCO3][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NSOLC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NSOLC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TSOLC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TSOLC][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_ALK ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ALK][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTEC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTEC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTC][dm];
      }
    }


    // *************** Nitrogen stocks in ecosystems ***********

    else if( predname.at( i ) == tem.predstr.at( tem.I_STRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_STRN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_STON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_STON][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SOLN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SOLN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DON][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NH4 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NH4][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NO3 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NO3][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VEGN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VEGN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NSOLN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NSOLN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TSOLN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TSOLN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AVLN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AVLN][dm] * GRAMS2MG;
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_DEADWOODN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DEADWOODN][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_FOLC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FOLC][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_STEMC ) )
        {
          for( dm = 0; dm < CYCLE; ++dm )
          {
            tempred.mon[dm] = output[tem.I_STEMC][dm];
          }
        }
    else if( predname.at( i ) == tem.predstr.at( tem.I_CROOTC ) )
        {
          for( dm = 0; dm < CYCLE; ++dm )
          {
            tempred.mon[dm] = output[tem.I_CROOTC][dm];
          }
        }
    else if( predname.at( i ) == tem.predstr.at( tem.I_FROOTC ) )
        {
          for( dm = 0; dm < CYCLE; ++dm )
          {
            tempred.mon[dm] = output[tem.I_FROOTC][dm];
          }
        }
    else if( predname.at( i ) == tem.predstr.at( tem.I_AGR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGR][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_AGL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGL][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_BGR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_BGR][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_BGL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_BGL][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_CWD ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CWD][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_MDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_MDC][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_AGE ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGE][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_SOC) )
    {
       for( dm = 0; dm < CYCLE; ++dm )
       {
         tempred.mon[dm] = output[tem.I_SOC][dm];
       }
    }
    // *****************Water stocks in ecosystems *************

    else if( predname.at( i ) == tem.predstr.at( tem.I_AVLW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AVLW][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SM ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SM][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VSM ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VSM][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PCTP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PCTP][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SNWPCK ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SNWPCK][dm];
        //if (dm==2) cout <<tempred.mon[dm]<<" ";
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RGRW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RGRW][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SGRW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SGRW][dm];
      }
    }


   // ******************** Phenology ***************************


    else if( predname.at( i ) == tem.predstr.at( tem.I_UNRMLF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_UNRMLF][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LEAF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LEAF][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LAI ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LAI][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_FPC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FPC][dm] * PROP2PCT;
      }
    }


    // *************** Carbon fluxes in ecosystems *************


    else if( predname.at( i ) == tem.predstr.at( tem.I_INGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INGPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_GPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_GPP][dm];

        //if (output[tem.I_GPP][dm]>0) cout <<"pdyr: "<<pdyr + 1800<<" mon: "<<dm <<" gppoutput3: "<<tempred.mon[dm]<<endl;
      }
    }

    // *********************** Ozone Effects *******************

    else if( predname.at( i ) == tem.predstr.at( tem.I_FOZONE ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FOZONE][dm] * PROP2PCT;

      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_FINDOZONE ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FINDOZONE][dm] * PROP2PCT;

      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_INNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INNPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_GPR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_GPR][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RVMNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RVMNT][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RVGRW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RVGRW][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_ABVGPR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ABVGPR][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_ROOTGPR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ROOTGPR][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_ISOPREN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ISOPREN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TERPEN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TERPEN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_ORVOC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ORVOC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_OVOC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_OVOC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VOC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VOC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LTRC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LTRC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CDCMP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CDCMP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RH ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RH][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RSOIL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RSOIL][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DOCP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DOCP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LCHDOC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LCHDOC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_ERDPOC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ERDPOC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NEP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NEP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NCE ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NCE][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CO2DISS ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CO2DISS][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LCHCO2 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LCHCO2][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_HCO3P ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_HCO3P][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LCHHCO3 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LCHHCO3][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RHCO3P ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RHCO3P][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LCHALK ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LCHALK][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NTCB ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NTCB][dm];
      }
    }

    // ************** Nitrogen fluxes in ecosystems ************


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGFRTN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGFRTN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_BNFIX ) )
    {
      for ( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_BNFIX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SNFIX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SNFIX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_ANFIX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ANFIX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_INNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INNUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_INNH4UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INNH4UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_INNO3UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INNO3UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VNUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VNH4UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VNH4UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VNO3UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VNO3UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VSUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VSUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VLUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VLUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VNMBL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VNMBL][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VNRSRB ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VNRSRB][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LTRN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NDCMP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NDCMP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DONP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DONP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_GMIN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_GMIN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NH4IMM ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NH4IMM][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NIMM ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NIMM][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NMIN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NMIN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AMMN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AMMN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NTRF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NTRF][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NO3P ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NO3P][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NOP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NOP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_N2OP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_N2OP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_N2P ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_N2P][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DNTRF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DNTRF][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NH3FLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NH3FLX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NOFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NOFLX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_N2OFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_N2OFLX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_N2FLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_N2FLX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LCHNO3 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LCHNO3][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_LCHDON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_LCHDON][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_ERDPON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ERDPON][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NINP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NINP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NLST ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NLST][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NTNB ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NTNB][dm];
      }
    }

    // *****************Water fluxes in ecosystems *************

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGIRRIG ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGIRRIG][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_INEET ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_INEET][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_EET ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_EET][dm];
      }
    }

   else if( predname.at( i ) == tem.predstr.at( tem.I_RPERC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RPERC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SPERC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SPERC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RRUN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RRUN][dm];
      }
    }


    else if( predname.at( i ) == tem.predstr.at( tem.I_SRUN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SRUN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PET ) )
    {
      for ( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PET][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SNWINF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SNWINF][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_WYLD ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_WYLD][dm];
      }
    }


// ************** Carbon stocks in products ********************


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGPRDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGPRDC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PROD10C ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PROD10C][dm];      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PROD100C ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PROD100C][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTPRDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTPRDC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RESIDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RESIDC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGSTUBC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGSTUBC][dm];
      }
    }

    // ************** Nitrogen stocks in products **************


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGPRDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGPRDN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PROD10N ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PROD10N][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PROD100N ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PROD100N][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTPRDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTPRDN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RESIDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RESIDN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGSTUBN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGSTUBN][dm];
      }
    }


    // *** Carbon fluxes during agricultural conversion ********


    else if( predname.at( i ) == tem.predstr.at( tem.I_CNVRTC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CNVRTC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VCNVRTC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VCNVRTC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SCNVRTC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SCNVRTC][dm];
        //if (output[tem.I_SCNVRTC][dm] >0.0) cout <<"sconvertc2: "<<tempred.mon[dm]<<" I_SCONV: "<<tem.I_SCNVRTC<<endl;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SLASHC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SLASHC][dm];
        //if (col == -128.125 && row == 59.875 && pdyr <=120) cout <<"pdyr: " <<pdyr + 1800 <<" month: " <<dm<<" cohort: " <<pichrt <<" cmnt: " <<tem.veg.cmnt<<" curveg: " <<tem.veg.getCURRENTVEG() << " I_SLASHC: "<<output[tem.I_SLASHC][dm]<<endl;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CFLX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CFLX][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_FFLC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FFLC][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_FFDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FFDC][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_SWFC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SWFC][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_DWFC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DWFC][dm];
      }
    }


    // *** Nitrogen fluxes during agricultural conversion ******


    else if( predname.at( i ) == tem.predstr.at( tem.I_CNVRTN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CNVRTN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_VCNVRTN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_VCNVRTN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SCNVRTN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_SCNVRTN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_SLASHN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
      	 tempred.mon[dm] = output[tem.I_SLASHN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NRETNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NRETNT][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NVRTNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NVRTNT][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NSRTNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NSRTNT][dm];
      }
    }


    // ************** Carbon fluxes to/from products ***********

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGFPRDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGFPRDC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRDF10C ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRDF10C][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRDF100C ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRDF100C][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTFPRDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTFPRDC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_FRESIDC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FRESIDC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGPRDFC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGPRDFC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRD10FC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRD10FC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRD100FC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRD100FC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTPRDFC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTPRDFC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RESIDFC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RESIDFC][dm];
      }
    }

    // ************** Nitrogen fluxes to/from products *********


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGFPRDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGFPRDN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRDF10N ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRDF10N][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRDF100N ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRDF100N][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTFPRDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTFPRDN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_FRESIDN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FRESIDN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGPRDFN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGPRDFN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRD10FN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRD10FN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_PRD100FN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_PRD100FN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_TOTPRDFN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TOTPRDFN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_RESIDFN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_RESIDFN][dm];
      }
    }

    // ************** Carbon stocks in crops   *****************


    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVEGC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVEGC][dm];
      }
    }


    // ************** Nitrogen stocks in crops *****************


    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVEGN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVEGN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CSTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CSTRN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATSTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATSTRN][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CSTON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CSTON][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATSTON ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATSTON][dm];
      }
    }

    // ******************** Crop Phenology *********************


    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPULF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPULF][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATULF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATULF][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPLEAF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPLEAF][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATLEAF ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATLEAF][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPLAI ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPLAI][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATLAI ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATLAI][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_CROPFPC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_CROPFPC][dm] * PROP2PCT;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATFPC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATFPC][dm] * PROP2PCT;
      }
    }

    // ************** Carbon fluxes in croplands ***************

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGINGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGINGPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATINGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATINGPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGGPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATGPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATGPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGINNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGINNPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATINNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATINNPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGNPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATNPP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATNPP][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGGPR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGGPR][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATGPR ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATGPR][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGRVMNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGRVMNT][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATRVMNT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATRVMNT][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGRVGRW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGRVGRW][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATRVGRW ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATRVGRW][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGLTRC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGLTRC][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATLTRC ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATLTRC][dm];
      }
    }

    // ************** Nitrogen fluxes in croplands *************


    else if( predname.at( i ) == tem.predstr.at( tem.I_AGSNFX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGSNFX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATSNFX ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATSNFX][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGINNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGINNUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATINNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATINNUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AINNH4UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AINNH4UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NINNH4UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NINNH4UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AINNO3UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AINNO3UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NINNO3UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NINNO3UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVNUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVNUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVNUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVNH4UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVNH4UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NVNH4UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NVNH4UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVNO3UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVNO3UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NVNO3UP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NVNO3UP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVSUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVSUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVSUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVSUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVLUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVLUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVLUP ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVLUP][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVNMBL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVNMBL][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATVNMBL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATVNMBL][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGVNRSRB ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGVNRSRB][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NVNRSRB ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NVNRSRB][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_AGLTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_AGLTRN][dm] * GRAMS2MG;
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_NATLTRN ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_NATLTRN][dm] * GRAMS2MG;
      }
    }

    // ************** Soil temperature data ********************


    else if( predname.at( i ) == tem.predstr.at( tem.I_TSOIL ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_TSOIL][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DST0 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DST0][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DST5 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DST5][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DST10 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DST10][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DST20 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DST20][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DST50 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DST50][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DST100 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DST100][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_DST200 ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_DST200][dm];
      }
    }
    else if( predname.at( i ) == tem.predstr.at( tem.I_DST300 ) )
    {
       for( dm = 0; dm < CYCLE; ++dm )
       {
         tempred.mon[dm] = output[tem.I_DST300][dm];
       }
     }
    else if( predname.at( i ) == tem.predstr.at( tem.I_FRONTD ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_FRONTD][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_THAWBE ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_THAWBE][dm];
      }
    }

    else if( predname.at( i ) == tem.predstr.at( tem.I_THAWEND ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_THAWEND][dm];
      }
    }

   // ************** Thaw percent (i.e., f(FT)) ****************

   else if( predname.at( i ) == tem.predstr.at( tem.I_THAWPCT ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_THAWPCT][dm] * PROP2PCT;
      }
    }


    else if( predname.at( i ) == tem.predstr.at( tem.I_ACTLAYER ) )
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = output[tem.I_ACTLAYER][dm] * METERS2MM;
        //if (dm==2) cout <<"actlayer: "<< tempred.mon[dm]<<endl;
      }
    }

    else
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        tempred.mon[dm] = MISSING;
      }
    }

    // Write output data to files
    //if (i==12 && output[115][6]>0.0) cout <<"mondata: "<<tempred.mon[6]<<" output: "<<output[115][6]<<endl;

    if( predname.at( i ) == tem.predstr.at( tem.I_VSM )
        || predname.at( i ) == tem.predstr.at( tem.I_PCTP )
        || predname.at( i ) == tem.predstr.at( tem.I_LEAF )
        || predname.at( i ) == tem.predstr.at( tem.I_THAWPCT )
        || predname.at( i ) == tem.predstr.at( tem.I_FOZONE )
        || predname.at( i ) == tem.predstr.at( tem.I_FINDOZONE ) )
    {
      tempred.poutdel( fout[i],
                       col,
                       row,
                       predname.at( i ),
                       (pichrt+1),
                       cohort[pichrt].standage,
                       tem.veg.getPOTVEG(),
                       tem.veg.getCURRENTVEG(),
                       tem.veg.getSUBTYPE(),
                       tem.veg.cmnt,
                       (PROP2PCT * tem.soil.getPSIPLUSC()),
                       tem.qualcon[pdyr],
                       carea,
                       cohort[pichrt].chrtarea,
                       ttotyr[pdyr],
                       tempred.mon,
                       region );
    }
    else
    {
      tempred.outdel( fout[i],
                      col,
                      row,
                      predname.at( i ),
                      (pichrt+1),
                      cohort[pichrt].standage,
                      tem.veg.getPOTVEG(),
                      tem.veg.getCURRENTVEG(),
                      tem.veg.getSUBTYPE(),
                      tem.veg.cmnt,
                      (PROP2PCT * tem.soil.getPSIPLUSC()),
                      tem.qualcon[pdyr],
                      carea,
                      cohort[pichrt].chrtarea,
                      ttotyr[pdyr],
                      tempred.mon,
                      region );


      	  	  	  	  //below code is added to output site specific outputs

       	   	   	   	  double xx=0.0;
       	   	   	      //if ( (col==-98.25 && row==55.5 && ttotyr[pdyr]>=1960) || (col==-99.75 && row==59.5 && ttotyr[pdyr]>=1993) || (col==-99.75 && row==59.5 && ttotyr[pdyr]>=1993) || (col==-99.25 && row==55.75 && ttotyr[pdyr]>=1988) || (col==-98.75 && row==55.75 && ttotyr[pdyr]>=1980))//disturb1961, disturb1998/1994, disturb1989, disturb1981
       	   	   	   	  //if (col==-99.25 && row==55.75 && ttotyr[pdyr]>=1988)//disturb1989
       	   	   	      //if (col == -76.50 && row ==70.0 && ttotyr[pdyr] > 499)
       	   	   	      {
       	   	   	   	  for (int t=0; t<12;t++)
       	   	   	   	  {
       	   	   	   		  //if (predname.at(i)=="NPP"||predname.at(i)== "SCONVRTC"|| predname.at(i)=="VCONVRTC"||predname.at(i)=="NTCB") xx+=tempred.mon[t];
       	   	   	   		  //else if (predname.at(i)=="TOTSOLC" ||predname.at(i)=="VEGC"||predname.at(i)=="DEADWOODC") xx= tempred.mon[11];
       	   	   	   		  //if (predname.at(i)=="NTCB") cout <<"col: "<<col<<" row: "<<row<<" year: "<<ttotyr[pdyr]<<" mon: "<<t<<" NTCB: "<<tempred.mon[t]<<endl;

       	   	   	   	  }/*
      	  	  	  	  //if (predname.at(i)=="NPP" || predname.at(i)=="TOTSOLC" ||predname.at(i)=="VEGC"||predname.at(i)=="DEADWOODC"||predname.at(i)== "SCONVRTC"|| predname.at(i)=="VCONVRTC"||predname.at(i)=="NTCB")
       	   	   	   		  //cout <<"col: "<<col<<" row: "<<row<<" var: "<<predname.at(i)<<" cohort: "<<pichrt+1<<" age: "<<cohort[pichrt].standage<<" cmnt: "<<tem.veg.cmnt<<" year: " <<ttotyr[pdyr]<<" value: "<<xx<<endl;

                      if (predname.at(i)=="VCONVRTC" && output[tem.I_VCNVRTC][6]>0.0)
                      {
                    	  cout <<"col: "<<col<<" row: "<<row<<" name: "<<predname.at( i )<<" currentveg: "<< tem.veg.getCURRENTVEG()<<" mon6: "<< tempred.mon[6]<<" outputvconvert: "<<output[tem.I_VCNVRTC][6]<<endl;
                    	  cout <<"xx: "<<tempred.mon[6]<<endl;
                    	  exit(-1);
                      }*/

       	   	   	      }

      	  	  	  	  //if (output[tem.I_VEGC][dm]>0.0) cout <<" col: "<<col<<" row: "<<row<<" ttotyr: "<<ttotyr[pdyr]<<" pdm: "<<dm<<" cohort: "<<pichrt<<" cmnt: "<<tem.veg.cmnt<<" vegc: "<< output[tem.I_VEGC][dm]<<" stemc: "<<output[tem.I_STEMC][dm]<<" leafc: "<<output[tem.I_FOLC][dm]<<" crootc: "<<output[tem.I_CROOTC][dm]<<" frootc: "<<output[tem.I_FROOTC][dm]<<" soc: "<<output[tem.I_SOLC][dm]<<" agr: "<<output[tem.I_AGR][dm]<<" agl: "<<output[tem.I_AGL][dm]<<" bgr: "<<output[tem.I_BGR][dm]<<" bgl: "<<output[tem.I_BGL][dm]<<" CWD: "<<output[tem.I_CWD][dm]<<endl;

    }
  }


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMelmnt60::transqc( int& maxyears,
                         int& totyr,
                         double plantc[CYCLE] )
{

  int dm;
  int ccc;
  double sumcarbon = ZERO;
  ccc = ACCEPT;

  if( totyr < 0 || totyr >= maxyears ) { return ccc = 30; }

  for( dm = 0; dm < CYCLE; ++dm ) { sumcarbon += plantc[dm]; }

  if( sumcarbon <= 0.1) { return ccc = TQCZEROFLAG; }

  return ccc;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMelmnt60::updateTEMmonth( ofstream& rflog1,
                                 const int& equil,
                                 const int& totsptime,
                                 const int& pdyr,
                                 const int& pdm,
                                 const int& pichrt )
{

  // Pass cohort characteristics information to TEM
  getTEMCohortState( pichrt, pdm );

  //if (pdyr == 1800 && pdm == 6 && cohort[pichrt].y[I_SOLC] !=cohort[pichrt].y[I_SOLC]) cout <<"year1: " <<pdyr<<" qc: " <<cohort[pichrt].qc <<" tqc: " <<cohort[pichrt].tqc<<" cmnt1: " <<tem.veg.cmnt<<" cmnt2: "<<cohort[pichrt].cmnt<<" pichrt: " <<pichrt <<" vegc: "<<cohort[pichrt].y[I_VEGC]<<" soilc: " <<cohort[pichrt].y[I_SOLC]<<endl;
  //if (pdm ==6) cout <<"col: " <<col <<" row: " <<row <<" year1: " <<pdyr<<" qc: " <<cohort[pichrt].qc <<" tqc: " <<cohort[pichrt].tqc<<" pichrt: " <<pichrt<<" sconvert: " <<cohort[pichrt].sconvert<<endl;
  //modified by cgs col == -75.5 && row == 82.75 &&
  //if( ACCEPT == cohort[pichrt].qc )
  {
	//cohort[pichrt].qc = tem.ecdqc( tem.veg.cmnt );
	//if (cohort[pichrt].chrtarea <0.05) cohort[pichrt].qc=198; //added by cgs to assign no-data to the grid cells with land area ~0.0 km2

	//if ((cohort[pichrt].qc!=0 ||cohort[pichrt].tqc!=0) && tem.veg.cmnt !=1) cout <<"qc: "<<cohort[pichrt].qc<<" tqc: "<<cohort[pichrt].tqc <<" cmnt: "<<tem.veg.cmnt<<endl;
  }
  //if (pdyr == 1800 && pdm == 6 && cohort[pichrt].y[I_SOLC] !=cohort[pichrt].y[I_SOLC]) cout <<"qc2: "<<cohort[pichrt].qc<<endl;
  if( ACCEPT == cohort[pichrt].qc
      && ACCEPT == cohort[pichrt].tqc )
  //if(cohort[pichrt].qc<100) //modified by cgs
  {
    // Pass monthly climate information to TEM

    tem.atms.setGIRR( climate[clm.I_GIRR][pdm] );
    tem.atms.setCLDS( climate[clm.I_CLDS][pdm] );
    tem.atms.setNIRR( climate[clm.I_NIRR][pdm] );
    tem.atms.setPAR(  climate[clm.I_PAR][pdm] );
    tem.atms.setTAIR( climate[clm.I_TAIR][pdm] );
    tem.atms.setRAIN( climate[clm.I_RAIN][pdm] );
    tem.atms.setPREC( climate[clm.I_PREC][pdm] );
    tem.atms.setSNOWFALL( climate[clm.I_SNWFAL][pdm] );
    tem.atms.setCO2( climate[clm.I_CO2][pdm] );
    tem.atms.setAOT40( climate[clm.I_AOT40][pdm] );
    tem.atms.setTOTNDEP( climate[clm.I_TNDEP][pdm] );
    tem.atms.setNH4DEP( climate[clm.I_NH4DEP][pdm] );
    tem.atms.setNO3DEP( climate[clm.I_NO3DEP][pdm] );
    tem.atms.setMXTAIR( mxtair );
    tem.atms.yrprec = yrprec;

    tem.soil.stm.setNEXTTAIR( climate[clm.I_TAIR][pdm+1] );
    tem.soil.stm.setNEXTSNOWFALL( climate[clm.I_SNWFAL][pdm+1] );
    //if (pdyr ==3 || pdyr ==4) cout <<"girr: "<<climate[clm.I_GIRR][pdm]<<" yrprec: "<<yrprec<<" snowfall: "<< climate[clm.I_SNWFAL][pdm]<<" rain: "<<climate[clm.I_RAIN][pdm]<<" Tair: "<< climate[clm.I_TAIR][pdm]<<endl;
    tem.setLON(col); //added by cgs to pass the latitude and longitude information to tem
    tem.setLAT(row);

    // Determine CDM for current year for soil thermal model

    if( 0 == pdm )
    {
      tem.soil.stm.updateyrCDM( climate[clm.I_TAIR] );
    }

    tem.baseline = 0;
    tem.wrtyr = -99;

    tem.totyr = clm.modstartyr
                - totsptime - 1
                + (pdyr * tem.diffyr);

//    if( 1 == tem.ag.state &&  tem.totyr >= 1950 )
//    {
//      tem.ag.fertflag = 1;
//    }
//    else
//    {
//      tem.ag.fertflag = 0;
//    }


    // Run the Terrestrial Ecosystem Model (TEM) under
    //   transient conditions
    //if (col > -109 && col <-76 && row == 54 && pdyr <=1) cout <<"pdyr2: "<<pdyr<<" pichrt: " <<pichrt<<endl;
    //if (pdyr <=1) cout <<"pdyr: " <<pdyr<<" cohort0: " <<pichrt<<" veg c0: "<<tem.veg.getVEGC() <<" soil c0: "<<tem.soil.getTSOLC()<<" y1: " <<cohort[pichrt].y[1]<<endl;
    wrtyr = tem.monthlyTransient( pdyr,
                                  pdm,
                                  tem.tol,
                                  rflog1 );

    // Update telmnt[0].cohort.dst10 for future use

    cohort[pichrt].dst10[pdm] = tem.soil.getDST10();

    // Save TEM output to telmnt[0].output

    outputTEMmonth( pdm );
    //if (tem.veg.cmnt==4 && col == -98.0 && row== 55.50 && pichrt==0 && pdyr<10) cout <<" pdyr2: "<<pdyr<<" vegc: "<< output[tem.I_VEGC][pdm]<<" soc: "<<output[tem.I_SOLC][pdm] <<" gpp: "<<output[tem.I_GPP][pdm]<<endl;
   //if (tem.veg.cmnt == 4 && col ==-80.75 && row == 51.25 && pdyr <200) cout <<"lat: "<<row<<" lon: "<<col<<" pdyr: " <<pdyr<<" pdm: "<<pdm <<" ntcb2: "<<output[tem.I_NTCB][pdm]<<endl;

    //if (col ==-83.5 && row == 36.5 && tem.veg.cmnt==5) cout << "year: " <<dyr<<" cmnt: "<<tem.veg.cmnt<<" vegc: "<<tem.getY(I_VEGC)<<" TSOILC: "<<tem.soil.getTSOLC()<<" gpp: " <<tem.veg.getGPP() << " gpr: "<<tem.veg.getGPR() <<" litterc: "<<tem.veg.getLTRFALC()<<endl;

    //if (output[tem.I_SOLC][pdm]!=output[tem.I_SOLC][pdm] && pdm == 7 && pdyr<=1) cout <<" I_soilC0:" <<xx<<" col: " <<col<< "row: "<<row <<"value of soilc1: "<<is_valid(xx)<< " or: "<<is_nan(xx) <<" or: " <<is_infinite(xx)<<endl;
    //if (col == -123.0 && row== 44.25 && tem.veg.cmnt==7)
    //if (tem.veg.cmnt==7) cout <<"pdyr: "<<pdyr<<" pdm: "<<pdm<<" stemc: "<<output[tem.I_STEMC][pdm]<<" GPP :"<< output[tem.I_GPP][pdm]<<" vegc : "<<output[tem.I_VEGC][pdm]<< endl;
    ttotyr[pdyr] = tem.totyr;
  } // End of qc == ACCEPT and tqc = ACCEPT
  else
  {
    if( (CYCLE-1) == pdm )
    {
      // Set missing values to telmnt[0].output

      setTEMmiss( pdyr,
                  equil,
                  totsptime,
                  pichrt );
    }
  }

//  if( 12 == pdyr && 1 == pichrt && (CYCLE-1) == pdm ) { exit( -1 ); }

  // Save TEM state for cohort
  //cout <<"qc2: "<<cohort[pichrt].qc<<"tqc: " <<cohort[pichrt].tqc<<endl;
  saveTEMCohortState( pichrt );
  //if (pdm == 6 && cohort[pichrt].cmnt == 14) cout <<" year2: " <<pdyr<<" cmnt1: " <<tem.veg.cmnt<<" cmnt2: "<<cohort[pichrt].cmnt<<" pichrt: " <<pichrt <<" vegc: "<<cohort[pichrt].y[I_VEGC]<<" soilc: " <<cohort[pichrt].y[I_SOLC]<<endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMelmnt60::writeCohortState( ofstream& ofstate,
                                   const int& pichrt )
{
  int dm;
  int dnode;
  int i;

  ofstate << col << " ";
  ofstate << row << " ";
  ofstate << (pichrt+1) << " ";

  ofstate << cohort[pichrt].srcCohort << " ";
  ofstate << cohort[pichrt].standage << " ";
  ofstate << cohort[pichrt].chrtarea << " ";
  ofstate << cohort[pichrt].potveg << " ";
  ofstate << cohort[pichrt].currentveg << " ";
  ofstate << cohort[pichrt].subtype << " ";
  ofstate << cohort[pichrt].cmnt << " ";


  for( i = 0; i < MAXSTATE; ++i )
  {
    ofstate << cohort[pichrt].y[i] << " ";
    ofstate << cohort[pichrt].prevy[i] << " ";
    //cout <<"i: " <<2<<" y[2]: " <<cohort[pichrt].y[2]<<endl;
  }

  ofstate << cohort[pichrt].agcmnt << " ";

  ofstate << cohort[pichrt].aggrowdd << " ";

  ofstate << cohort[pichrt].agkd << " ";

  ofstate << cohort[pichrt].agprvstate << " ";

  ofstate << cohort[pichrt].agstate << " ";

  ofstate << cohort[pichrt].c2n << " ";

  ofstate << cohort[pichrt].cneven << " ";

  ofstate << cohort[pichrt].convrtflx.carbon << " ";
  ofstate << cohort[pichrt].convrtflx.nitrogen << " ";

  ofstate << cohort[pichrt].cropprveetmx << " ";

  ofstate << cohort[pichrt].cropprvleafmx << " ";

  ofstate << cohort[pichrt].cropprvpetmx << " ";

  ofstate << cohort[pichrt].cropResidue.carbon << " ";
  ofstate << cohort[pichrt].cropResidue.nitrogen << " ";

  ofstate << cohort[pichrt].croptopt << " ";

  ofstate << cohort[pichrt].distmnthcnt << " ";

  ofstate << cohort[pichrt].disturbflag << " ";

  ofstate << cohort[pichrt].disturbmonth << " ";

  for( dm = 0; dm < CYCLE; ++dm )
  {
    ofstate << cohort[pichrt].dst10[dm] << " ";
  }

  ofstate << cohort[pichrt].eetmx << " ";

  ofstate << cohort[pichrt].fertflag << " ";

  ofstate << cohort[pichrt].firemnthcnt << " ";

  ofstate << cohort[pichrt].firendep << " ";

  ofstate << cohort[pichrt].formPROD10.carbon << " ";
  ofstate << cohort[pichrt].formPROD10.nitrogen << " ";

  ofstate << cohort[pichrt].formPROD100.carbon << " ";
  ofstate << cohort[pichrt].formPROD100.nitrogen << " ";

  ofstate << cohort[pichrt].fprevozone << " ";

  ofstate << cohort[pichrt].FRI << " ";

  for( dm = 0; dm < CYCLE; ++dm )
  {
    ofstate << cohort[pichrt].initPROD1[dm].carbon << " ";
    ofstate << cohort[pichrt].initPROD1[dm].nitrogen << " ";
  }

  for( i = 0; i < 10; ++i )
  {
    ofstate << cohort[pichrt].initPROD10[i].carbon << " ";
    ofstate << cohort[pichrt].initPROD10[i].nitrogen << " ";
  }

  for( i = 0; i < 100; ++i )
  {
    ofstate << cohort[pichrt].initPROD100[i].carbon << " ";
    ofstate << cohort[pichrt].initPROD100[i].nitrogen << " ";
  }

  ofstate << cohort[pichrt].irrgflag << " ";

  ofstate << cohort[pichrt].kd << " ";

  ofstate << cohort[pichrt].natprveetmx << " ";

  ofstate << cohort[pichrt].natprvleafmx << " ";

  ofstate << cohort[pichrt].natprvpetmx << " ";

  ofstate << cohort[pichrt].natseedC << " ";

  ofstate << cohort[pichrt].natseedSTRN << " ";

  ofstate << cohort[pichrt].natseedSTON << " ";

  ofstate << cohort[pichrt].natsoil << " ";

  ofstate << cohort[pichrt].nattopt << " ";

  ofstate << cohort[pichrt].natyreet << " ";

  ofstate << cohort[pichrt].natyrpet << " ";

  ofstate << cohort[pichrt].newleafmx << " ";

  ofstate << cohort[pichrt].newtopt << " ";

  ofstate << cohort[pichrt].nonsolc << " ";

  ofstate << cohort[pichrt].nonsoln << " ";

  ofstate << cohort[pichrt].nretent << " ";

  ofstate << cohort[pichrt].nsretent << " ";

  ofstate << cohort[pichrt].nvretent << " ";

  ofstate << cohort[pichrt].petmx << " ";

  ofstate << cohort[pichrt].prev2tair << " ";

  ofstate << cohort[pichrt].prevco2 << " ";

  ofstate << cohort[pichrt].prevCropResidue.carbon << " ";
  ofstate << cohort[pichrt].prevCropResidue.nitrogen << " ";

  ofstate << cohort[pichrt].prevdst10 << " ";

  ofstate << cohort[pichrt].prevPROD1.carbon << " ";
  ofstate << cohort[pichrt].prevPROD1.nitrogen << " ";

  ofstate << cohort[pichrt].prevPROD10.carbon << " ";
  ofstate << cohort[pichrt].prevPROD10.nitrogen << " ";

  ofstate << cohort[pichrt].prevPROD100.carbon << " ";
  ofstate << cohort[pichrt].prevPROD100.nitrogen << " ";

  ofstate << cohort[pichrt].prevspack << " ";

  ofstate << cohort[pichrt].prevtair << " ";
  ofstate << cohort[pichrt].prevrain << " ";

  ofstate << cohort[pichrt].prevunrmleaf << " ";

  //ofstate << cohort[pichrt].prod10par << " ";

  //ofstate << cohort[pichrt].prod100par << " ";

  ofstate << cohort[pichrt].productYear << " ";

  ofstate << cohort[pichrt].prvchrtarea << " ";

  ofstate << cohort[pichrt].prvcropnpp << " ";

  ofstate << cohort[pichrt].prveetmx << " ";

  ofstate << cohort[pichrt].prvleafmx << " ";

  ofstate << cohort[pichrt].prvpetmx << " ";

  ofstate << cohort[pichrt].qc << " ";

  //ofstate << cohort[pichrt].sconvert << " ";

  ofstate << cohort[pichrt].sconvrtflx.carbon << " ";
  ofstate << cohort[pichrt].sconvrtflx.nitrogen << " ";

  ofstate << cohort[pichrt].slash.carbon << " ";
  ofstate << cohort[pichrt].slash.nitrogen << " ";

  //ofstate << cohort[pichrt].slashpar << " ";

  ofstate << cohort[pichrt].STMis9 << " ";

  ofstate << cohort[pichrt].STMsmass9 << " ";

  for( dnode = 0; dnode < MAXNODES; ++dnode )
  {
    ofstate << cohort[pichrt].STMdx9[dnode] << " ";
    ofstate << cohort[pichrt].STMt9[dnode] << " ";
    ofstate << cohort[pichrt].STMwater9[dnode] << " ";
    ofstate << cohort[pichrt].STMx9[dnode] << " ";
    ofstate << cohort[pichrt].STMxfa9[dnode] << " ";
    ofstate << cohort[pichrt].STMxfb9[dnode] << " ";
  }

  for( dnode = 0; dnode < MAXSNODES; ++dnode )
  {
    ofstate << cohort[pichrt].STMweight9[dnode] << " ";
  }

  ofstate << cohort[pichrt].tillflag << " ";

  ofstate << cohort[pichrt].topt << " ";

  ofstate << cohort[pichrt].tqc << " ";

  //ofstate << cohort[pichrt].vconvert << " ";

  ofstate << cohort[pichrt].vconvrtflx.carbon << " ";
  ofstate << cohort[pichrt].vconvrtflx.nitrogen << " ";

  ofstate << cohort[pichrt].vrespar << " ";

  ofstate << cohort[pichrt].wfpsoff << " ";

  ofstate << cohort[pichrt].yrltrc << " ";
  ofstate << cohort[pichrt].yrltrn << " ";

  ofstate <<cohort[pichrt].leafmortpar<< " ";
  ofstate <<cohort[pichrt].rootmortpar<< " ";
  ofstate <<cohort[pichrt].stemmortpar<< " ";
  ofstate <<cohort[pichrt].leafslash<< " ";
  ofstate <<cohort[pichrt].leafconv<< " ";
  ofstate <<cohort[pichrt].rootslash<< " ";
  ofstate <<cohort[pichrt].rootconv<< " ";
  ofstate <<cohort[pichrt].stemslash<< " ";
  ofstate <<cohort[pichrt].stemconv<< " ";
  ofstate <<cohort[pichrt].standdead<< " ";
  ofstate <<cohort[pichrt].deadconv<< " ";
  ofstate <<cohort[pichrt].deadslash<< " ";
  ofstate <<cohort[pichrt].sconvert<< " ";
  ofstate <<cohort[pichrt].prod10par<< " ";
  ofstate <<cohort[pichrt].prod100par<< " ";
  ofstate << endl;
  //if (col == -137.625 && row == 59.625) cout <<"yrltrn: " <<cohort[pichrt].yrltrn <<" yrltrc: " <<cohort[pichrt].yrltrc<<endl;

};

void TEMelmnt60::readCohortState0( ifstream& ifstate,
                                  const int& pichrt )
{
  int i;
  int dm;
  int dnode;
  double dumflt1;
  double dumflt2;
  int dumint;
  string dumstr;

ifstate >> dumflt1;    // Longitude of element
ifstate >> dumflt2;    // Latitude of element

ifstate >> dumint;   // ichrt+1
//cout <<"lon: "<<dumflt1<<" lat: "<<dumflt2<<" col: "<<col<<"row: "<<row<<endl;
ifstate >> cohort[pichrt].srcCohort;
ifstate >> cohort[pichrt].standage;
ifstate >> cohort[pichrt].chrtarea;
ifstate >> cohort[pichrt].potveg;
ifstate >> cohort[pichrt].currentveg;
ifstate >> cohort[pichrt].subtype;
ifstate >> cohort[pichrt].cmnt;

//changed by cgs2014 to read different format initial cohort state data
//for( i = 0; i < MAXSTATE; ++i )
for( i = 0; i < 12; ++i )
{
  ifstate >> cohort[pichrt].y[i];
  ifstate >> cohort[pichrt].prevy[i];
  //if (i<=14) cout <<" y[" <<i <<"]: " <<cohort[pichrt].y[i]<<" ";

}
cohort[pichrt].y[12]=0.0;
cohort[pichrt].y[13]=0.0;

//if (cohort[pichrt].y[13]!=cohort[pichrt].y[13]) cout <<"y13 has no data " <<" cmnt: " <<cohort[pichrt].cmnt<<" icohort: " <<dumint<<" y[1]: " <<cohort[pichrt].y[1]<<" y0: " <<cohort[pichrt].y[0]<<endl;

//if (col == -129.5 && (row ==55))
//cout <<"col: " <<col<<" row: " <<row<<" ichrt: " <<pichrt<<" src: " <<cohort[pichrt].srcCohort <<" standage: " <<cohort[pichrt].standage<<" chrtarea: " <<cohort[pichrt].chrtarea<<" cmnt: " <<cohort[pichrt].cmnt<<" y0: " <<cohort[pichrt].y[0]<< " y1: " <<cohort[pichrt].y[1]<<endl;
ifstate >> cohort[pichrt].agcmnt;

ifstate >> cohort[pichrt].aggrowdd;

ifstate >> cohort[pichrt].agkd;

ifstate >> cohort[pichrt].agprvstate;

ifstate >> cohort[pichrt].agstate;

ifstate >> cohort[pichrt].c2n;

ifstate >> cohort[pichrt].cneven;

ifstate >> cohort[pichrt].convrtflx.carbon;
ifstate >> cohort[pichrt].convrtflx.nitrogen;

ifstate >> cohort[pichrt].cropprveetmx;

ifstate >> cohort[pichrt].cropprvleafmx;

ifstate >> cohort[pichrt].cropprvpetmx;

ifstate >> cohort[pichrt].cropResidue.carbon;
ifstate >> cohort[pichrt].cropResidue.nitrogen;

ifstate >> cohort[pichrt].croptopt;

ifstate >> cohort[pichrt].distmnthcnt;

ifstate >> cohort[pichrt].disturbflag;

ifstate >> cohort[pichrt].disturbmonth;

for( dm = 0; dm < CYCLE; ++dm )
{
  ifstate >> cohort[pichrt].dst10[dm];
}

ifstate >> cohort[pichrt].eetmx;

ifstate >> cohort[pichrt].fertflag;

ifstate >> cohort[pichrt].firemnthcnt;

ifstate >> cohort[pichrt].firendep;

ifstate >> cohort[pichrt].formPROD10.carbon;
ifstate >> cohort[pichrt].formPROD10.nitrogen;

ifstate >> cohort[pichrt].formPROD100.carbon;
ifstate >> cohort[pichrt].formPROD100.nitrogen;

ifstate >> cohort[pichrt].fprevozone;

ifstate >> cohort[pichrt].FRI;

for( dm = 0; dm < CYCLE; ++dm )
{
  ifstate >> cohort[pichrt].initPROD1[dm].carbon;
  ifstate >> cohort[pichrt].initPROD1[dm].nitrogen;
}

for( i = 0; i < 10; ++i )
{
  ifstate >> cohort[pichrt].initPROD10[i].carbon;
  ifstate >> cohort[pichrt].initPROD10[i].nitrogen;
}

for( i = 0; i < 100; ++i )
{
  ifstate >> cohort[pichrt].initPROD100[i].carbon;
  ifstate >> cohort[pichrt].initPROD100[i].nitrogen;
}

ifstate >> cohort[pichrt].irrgflag;

ifstate >> cohort[pichrt].kd;

//if (cohort[pichrt].kd!=cohort[pichrt].kd) cout <<"kd has no data: " <<endl;

ifstate >> cohort[pichrt].natprveetmx;
//cout <<"natprveetmx: "<<cohort[pichrt].natprveetmx<<endl;

ifstate >> cohort[pichrt].natprvleafmx;

ifstate >> cohort[pichrt].natprvpetmx;

ifstate >> cohort[pichrt].natseedC;

ifstate >> cohort[pichrt].natseedSTRN;

ifstate >> cohort[pichrt].natseedSTON;

ifstate >> cohort[pichrt].natsoil;

ifstate >> cohort[pichrt].nattopt;

ifstate >> cohort[pichrt].natyreet;

ifstate >> cohort[pichrt].natyrpet;

ifstate >> cohort[pichrt].newleafmx;

ifstate >> cohort[pichrt].newtopt;

ifstate >> cohort[pichrt].nonsolc;

ifstate >> cohort[pichrt].nonsoln;

ifstate >> cohort[pichrt].nretent;

ifstate >> cohort[pichrt].nsretent;

ifstate >> cohort[pichrt].nvretent;

ifstate >> cohort[pichrt].petmx;

ifstate >> cohort[pichrt].prev2tair;

ifstate >> cohort[pichrt].prevco2;

ifstate >> cohort[pichrt].prevCropResidue.carbon;
ifstate >> cohort[pichrt].prevCropResidue.nitrogen;

ifstate >> cohort[pichrt].prevdst10;

ifstate >> cohort[pichrt].prevPROD1.carbon;
ifstate >> cohort[pichrt].prevPROD1.nitrogen;

ifstate >> cohort[pichrt].prevPROD10.carbon;
ifstate >> cohort[pichrt].prevPROD10.nitrogen;

ifstate >> cohort[pichrt].prevPROD100.carbon;
ifstate >> cohort[pichrt].prevPROD100.nitrogen;

ifstate >> cohort[pichrt].prevspack;

ifstate >> cohort[pichrt].prevtair;
ifstate >> cohort[pichrt].prevrain;

ifstate >> cohort[pichrt].prevunrmleaf;

ifstate >> cohort[pichrt].prod10par;

ifstate >> cohort[pichrt].prod100par;

ifstate >> cohort[pichrt].productYear;

ifstate >> cohort[pichrt].prvchrtarea;

ifstate >> cohort[pichrt].prvcropnpp;

ifstate >> cohort[pichrt].prveetmx;

ifstate >> cohort[pichrt].prvleafmx;

ifstate >> cohort[pichrt].prvpetmx;

ifstate >> cohort[pichrt].qc;

ifstate >> cohort[pichrt].sconvert;

ifstate >> cohort[pichrt].sconvrtflx.carbon;
ifstate >> cohort[pichrt].sconvrtflx.nitrogen;

ifstate >> cohort[pichrt].slash.carbon;
ifstate >> cohort[pichrt].slash.nitrogen;

ifstate >> cohort[pichrt].slashpar;

ifstate >> cohort[pichrt].STMis9;

ifstate >> cohort[pichrt].STMsmass9;

for( dnode = 0; dnode < MAXNODES; ++dnode )
{
  ifstate >> cohort[pichrt].STMdx9[dnode];
  ifstate >> cohort[pichrt].STMt9[dnode];
  ifstate >> cohort[pichrt].STMwater9[dnode];
  ifstate >> cohort[pichrt].STMx9[dnode];
  ifstate >> cohort[pichrt].STMxfa9[dnode];
  ifstate >> cohort[pichrt].STMxfb9[dnode];
}

for( dnode = 0; dnode < MAXSNODES; ++dnode )
{
  ifstate >> cohort[pichrt].STMweight9[dnode];
}

ifstate >> cohort[pichrt].tillflag;

ifstate >> cohort[pichrt].topt;

ifstate >> cohort[pichrt].tqc;

ifstate >> cohort[pichrt].vconvert;

ifstate >> cohort[pichrt].vconvrtflx.carbon;
ifstate >> cohort[pichrt].vconvrtflx.nitrogen;

ifstate >> cohort[pichrt].vrespar;

ifstate >> cohort[pichrt].wfpsoff;

ifstate >> cohort[pichrt].yrltrc;
ifstate >> cohort[pichrt].yrltrn;
  /*
  else
  {
	  string tempstring;
	  for (int i=0;i<1646;i++) //total: 1648.
	  {
		ifstate>>tempstring;
		//cout<<tempstring<<" ";
	  }
	  //cout <<"end: "<<endl;

  }
  */
  //cout <<" yrtrn0: "<<cohort[pichrt].yrltrn<<" yrltrc: "<<cohort[pichrt].yrltrc <<endl;

//ifstate >> cohort[pichrt].standdeadpar;
//ifstate >> cohort[pichrt].leafmortpar;
//ifstate >> cohort[pichrt].stemmortpar;
//ifstate >> cohort[pichrt].rootmortpar;
//ifstate >> cohort[pichrt].deadconvert;
//ifstate >> cohort[pichrt].deadslash;

 ifstate.seekg( 0, ios::cur );

};

/* *************************************************************
************************************************************** */

void TEMelmnt60::getTEMCohortState( const int& pichrt,
                                    const int& pdm )
{
  int dnode;
  int i;
  int kdm;

  tem.veg.setPOTVEG( cohort[pichrt].potveg );

  tem.veg.setCURRENTVEG( cohort[pichrt].currentveg );

  tem.veg.setSUBTYPE( cohort[pichrt].subtype );

  tem.veg.cmnt = cohort[pichrt].cmnt;

  for( i = 0; i < MAXSTATE; ++i )
  {
    tem.setY( cohort[pichrt].y[i], i );
    tem.setPREVY( cohort[pichrt].prevy[i], i );
  }
  //if (col == -80.75 && row== 51.25 && tem.veg.cmnt==4) cout <<"pdm: "<<pdm<<" vegc: "<<cohort[pichrt].y[0]<<" prevegc: "<<cohort[pichrt].prevy[0]<<endl;
  tem.ag.cmnt = cohort[pichrt].agcmnt;
  //cout <<"agcmnt: "<<cohort[pichrt].agcmnt<<endl;

  tem.ag.setGROWDD( cohort[pichrt].aggrowdd );

  tem.ag.setKD( cohort[pichrt].agkd );

  tem.ag.prvstate = cohort[pichrt].agprvstate;

  tem.ag.state = cohort[pichrt].agstate;

  tem.veg.setC2N( cohort[pichrt].c2n );

  tem.veg.setCNEVEN( cohort[pichrt].cneven );

  tem.ag.setCONVRTFLXC( cohort[pichrt].convrtflx.carbon );
  tem.ag.setCONVRTFLXN( cohort[pichrt].convrtflx.nitrogen );

  tem.ag.setCROPPRVEETMX( cohort[pichrt].cropprveetmx );

  tem.ag.setCROPPRVLEAFMX( cohort[pichrt].cropprvleafmx );

  tem.ag.setCROPPRVPETMX( cohort[pichrt].cropprvpetmx );

  tem.ag.setCROPRESIDUEC( cohort[pichrt].cropResidue.carbon );
  tem.ag.setCROPRESIDUEN( cohort[pichrt].cropResidue.nitrogen );

  tem.ag.setCROPTOPT( cohort[pichrt].croptopt );

  tem.distmnthcnt = cohort[pichrt].distmnthcnt;

  tem.disturbflag = cohort[pichrt].disturbflag;

  tem.disturbmonth = cohort[pichrt].disturbmonth;

  if( (CYCLE-1) == pdm )
  {
    tem.soil.setNEXTDST10( cohort[pichrt].dst10[0] );
  }
  else
  {
    tem.soil.setNEXTDST10( cohort[pichrt].dst10[pdm+1] );
  }

  tem.soil.setEETMX( cohort[pichrt].eetmx );

  tem.ag.fertflag = cohort[pichrt].fertflag;

  tem.firemnthcnt = cohort[pichrt].firemnthcnt;

  tem.ag.setFIRENDEP( cohort[pichrt].firendep );

  tem.ag.setFORMPROD10C( cohort[pichrt].formPROD10.carbon );
  tem.ag.setFORMPROD10N( cohort[pichrt].formPROD10.nitrogen );

  tem.ag.setFORMPROD100C( cohort[pichrt].formPROD100.carbon );
  tem.ag.setFORMPROD100N( cohort[pichrt].formPROD100.nitrogen );

  tem.veg.setFPREVOZONE( cohort[pichrt].fprevozone );

  tem.ag.setFRI( cohort[pichrt].FRI );

  for( kdm = 0; kdm < CYCLE; ++kdm )
  {
    tem.ag.setINITPROD1C( cohort[pichrt].initPROD1[kdm].carbon,
                          kdm );

    tem.ag.setINITPROD1N( cohort[pichrt].initPROD1[kdm].nitrogen,
                          kdm );
  }

  for( i = 0; i < 10; ++i )
  {
    tem.ag.setINITPROD10C( cohort[pichrt].initPROD10[i].carbon, i );
    tem.ag.setINITPROD10N( cohort[pichrt].initPROD10[i].nitrogen, i );
  }

  for( i = 0; i < 100; ++i )
  {
    tem.ag.setINITPROD100C( cohort[pichrt].initPROD100[i].carbon, i );
    tem.ag.setINITPROD100N( cohort[pichrt].initPROD100[i].nitrogen, i );
  }

  tem.ag.irrgflag = cohort[pichrt].irrgflag;

  tem.microbe.setKD( cohort[pichrt].kd );

  tem.ag.setNATPRVEETMX( cohort[pichrt].natprveetmx );

  tem.ag.setNATPRVLEAFMX( cohort[pichrt].natprvleafmx );

  tem.ag.setNATPRVPETMX( cohort[pichrt].natprvpetmx );

  tem.ag.setNATSEEDC( cohort[pichrt].natseedC );

  tem.ag.setNATSEEDSTRN( cohort[pichrt].natseedSTRN );

  tem.ag.setNATSEEDSTON( cohort[pichrt].natseedSTON );

  tem.ag.setNATSOIL( cohort[pichrt].natsoil );

  tem.ag.setNATTOPT( cohort[pichrt].nattopt );

  tem.ag.setNATYREET( cohort[pichrt].natyreet );

  tem.ag.setNATYRPET( cohort[pichrt].natyrpet );

  tem.veg.setNEWLEAFMX( cohort[pichrt].newleafmx );

  tem.veg.setNEWTOPT( cohort[pichrt].newtopt );

  tem.soil.setNSOLC( cohort[pichrt].nonsolc );

  tem.soil.setNSOLN( cohort[pichrt].nonsoln );

  tem.ag.setNRETENT( cohort[pichrt].nretent );

  tem.ag.setNSRETENT( cohort[pichrt].nsretent );

  tem.ag.setNVRETENT( cohort[pichrt].nvretent );

  tem.atms.setPETMX( cohort[pichrt].petmx );

  tem.atms.setPREV2TAIR( cohort[pichrt].prev2tair );

  tem.atms.setPREVCO2( cohort[pichrt].prevco2 );

  tem.ag.setPREVCROPRESIDUEC( cohort[pichrt].prevCropResidue.carbon );
  tem.ag.setPREVCROPRESIDUEN( cohort[pichrt].prevCropResidue.nitrogen );

  tem.soil.setPREVDST10( cohort[pichrt].prevdst10 );

  tem.ag.setPREVPROD1C( cohort[pichrt].prevPROD1.carbon );
  tem.ag.setPREVPROD1N( cohort[pichrt].prevPROD1.nitrogen );

  tem.ag.setPREVPROD10C( cohort[pichrt].prevPROD10.carbon );
  tem.ag.setPREVPROD10N( cohort[pichrt].prevPROD10.nitrogen );

  tem.ag.setPREVPROD100C( cohort[pichrt].prevPROD100.carbon );
  tem.ag.setPREVPROD100N( cohort[pichrt].prevPROD100.nitrogen );

  tem.soil.setPREVSPACK( cohort[pichrt].prevspack );

  tem.atms.setPREVTAIR( cohort[pichrt].prevtair );

  tem.veg.setPREVUNRMLEAF( cohort[pichrt].prevunrmleaf );

  //tem.ag.setPROD10PAR( cohort[pichrt].prod10par );

  //tem.ag.setPROD100PAR( cohort[pichrt].prod100par );

  tem.ag.setPRODUCTYEAR( cohort[pichrt].productYear );

  tem.ag.setPRVCROPNPP( cohort[pichrt].prvcropnpp );

  tem.soil.setPRVEETMX( cohort[pichrt].prveetmx );

  tem.veg.setPRVLEAFMX( cohort[pichrt].prvleafmx );

  tem.atms.setPRVPETMX( cohort[pichrt].prvpetmx );

  //tem.ag.setSCONVERT( cohort[pichrt].sconvert );

  tem.ag.setSCONVRTFLXC( cohort[pichrt].sconvrtflx.carbon );
  tem.ag.setSCONVRTFLXN( cohort[pichrt].sconvrtflx.nitrogen );

  tem.ag.setSLASHC( cohort[pichrt].slash.carbon );
  tem.ag.setSLASHN( cohort[pichrt].slash.nitrogen );

  tem.ag.setdeadwoodc( cohort[pichrt].deadwood.carbon );
  tem.ag.setdeadwoodn( cohort[pichrt].deadwood.nitrogen );

  //tem.ag.setSLASHPAR( cohort[pichrt].slashpar );

  tem.soil.stm.setIS9( cohort[pichrt].STMis9 );

  tem.soil.stm.setSMASS9( cohort[pichrt].STMsmass9 );

  for( dnode = 0; dnode < MAXNODES; ++dnode )
  {
    tem.soil.stm.setDX9( cohort[pichrt].STMdx9[dnode], dnode );

    tem.soil.stm.setT9( cohort[pichrt].STMt9[dnode], dnode );

    tem.soil.stm.setWATER9( cohort[pichrt].STMwater9[dnode],
                            dnode );

    tem.soil.stm.setX9( cohort[pichrt].STMx9[dnode], dnode );

    tem.soil.stm.setXFA9( cohort[pichrt].STMxfa9[dnode],
                          dnode );

    tem.soil.stm.setXFB9( cohort[pichrt].STMxfb9[dnode],
                          dnode );
  }

  for( dnode = 0; dnode < MAXSNODES; ++dnode )
  {
    tem.soil.stm.setWEIGHT9( cohort[pichrt].STMweight9[dnode],
                            dnode );
  }

  tem.ag.tillflag = cohort[pichrt].tillflag;

  tem.veg.setTOPT( cohort[pichrt].topt );

  //tem.ag.setVCONVERT( cohort[pichrt].vconvert );
  tem.ag.setVCONVRTFLXC( cohort[pichrt].vconvrtflx.carbon );
  tem.ag.setVCONVRTFLXN( cohort[pichrt].vconvrtflx.nitrogen );
  //modified by cgs
  tem.ag.setVRESPAR( 1-cohort[pichrt].stemmortpar);

  tem.soil.setWFPSOFF( cohort[pichrt].wfpsoff );

  tem.veg.yrltrc = cohort[pichrt].yrltrc;
  tem.veg.yrltrn = cohort[pichrt].yrltrn;


  tem.ag.leafmortpar = cohort[pichrt].leafmortpar;
  tem.ag.rootmortpar = cohort[pichrt].rootmortpar;
  tem.ag.stemmortpar = cohort[pichrt].stemmortpar;
  tem.ag.leafslash = cohort[pichrt].leafslash;
  tem.ag.leafconv = cohort[pichrt].leafconv;
  tem.ag.rootslash = cohort[pichrt].rootslash;
  tem.ag.rootconv = cohort[pichrt].rootconv;
  tem.ag.stemslash = cohort[pichrt].stemslash;
  tem.ag.stemconv = cohort[pichrt].stemconv;
  tem.ag.standdead = cohort[pichrt].standdead;
  tem.ag.deadconv = cohort[pichrt].deadconv;
  tem.ag.deadslash = cohort[pichrt].deadslash;
  tem.ag.sconvert = cohort[pichrt].sconvert;
  tem.ag.prod10par = cohort[pichrt].prod10par;
  tem.ag.prod100par = cohort[pichrt].prod100par;
  tem.ag.forestage = cohort[pichrt].forestage;
  //cout <<"forestage: "<<tem.ag.forestage<<endl;
  //cout <<"col: "<<col <<" row: "<<row<<" cohort: "<<pichrt<<" chrtarea: "<<cohort[pichrt].chrtarea<< " qc: "<<cohort[pichrt].qc<< " tqc: "<<cohort[pichrt].tqc<<endl;
};
void TEMelmnt60::getTEMCohortState0( const int& pichrt,
                                    const int& pdm )
{
  int dnode;
  int i;
  int kdm;

  tem.veg.setPOTVEG( cohort[pichrt].potveg );

  tem.veg.setCURRENTVEG( cohort[pichrt].currentveg );

  tem.veg.setSUBTYPE( cohort[pichrt].subtype );

  tem.veg.cmnt = cohort[pichrt].cmnt;

  for( i = 0; i < MAXSTATE; ++i )
  {
    tem.setY( cohort[pichrt].y[i], i );
    tem.setPREVY( cohort[pichrt].prevy[i], i );
    //cout <<"i: "<<i<<" yi2: "<<cohort[pichrt].y[i]<<endl;
  }

  tem.ag.cmnt = cohort[pichrt].agcmnt;

  tem.ag.setGROWDD( cohort[pichrt].aggrowdd );

  tem.ag.setKD( cohort[pichrt].agkd );

  tem.ag.prvstate = cohort[pichrt].agprvstate;

  tem.ag.state = cohort[pichrt].agstate;

  tem.veg.setC2N( cohort[pichrt].c2n );

  tem.veg.setCNEVEN( cohort[pichrt].cneven );

  tem.ag.setCONVRTFLXC( cohort[pichrt].convrtflx.carbon );
  tem.ag.setCONVRTFLXN( cohort[pichrt].convrtflx.nitrogen );

  tem.ag.setCROPPRVEETMX( cohort[pichrt].cropprveetmx );

  tem.ag.setCROPPRVLEAFMX( cohort[pichrt].cropprvleafmx );

  tem.ag.setCROPPRVPETMX( cohort[pichrt].cropprvpetmx );

  tem.ag.setCROPRESIDUEC( cohort[pichrt].cropResidue.carbon );
  tem.ag.setCROPRESIDUEN( cohort[pichrt].cropResidue.nitrogen );

  tem.ag.setCROPTOPT( cohort[pichrt].croptopt );

  tem.distmnthcnt = cohort[pichrt].distmnthcnt;

  tem.disturbflag = cohort[pichrt].disturbflag;

  tem.disturbmonth = cohort[pichrt].disturbmonth;

  if( (CYCLE-1) == pdm )
  {
    tem.soil.setNEXTDST10( cohort[pichrt].dst10[0] );
  }
  else
  {
    tem.soil.setNEXTDST10( cohort[pichrt].dst10[pdm+1] );
  }

  tem.soil.setEETMX( cohort[pichrt].eetmx );

  tem.ag.fertflag = cohort[pichrt].fertflag;

  tem.firemnthcnt = cohort[pichrt].firemnthcnt;

  tem.ag.setFIRENDEP( cohort[pichrt].firendep );

  tem.ag.setFORMPROD10C( cohort[pichrt].formPROD10.carbon );
  tem.ag.setFORMPROD10N( cohort[pichrt].formPROD10.nitrogen );

  tem.ag.setFORMPROD100C( cohort[pichrt].formPROD100.carbon );
  tem.ag.setFORMPROD100N( cohort[pichrt].formPROD100.nitrogen );

  tem.veg.setFPREVOZONE( cohort[pichrt].fprevozone );

  tem.ag.setFRI( cohort[pichrt].FRI );

  for( kdm = 0; kdm < CYCLE; ++kdm )
  {
    tem.ag.setINITPROD1C( cohort[pichrt].initPROD1[kdm].carbon,
                          kdm );

    tem.ag.setINITPROD1N( cohort[pichrt].initPROD1[kdm].nitrogen,
                          kdm );
  }

  for( i = 0; i < 10; ++i )
  {
    tem.ag.setINITPROD10C( cohort[pichrt].initPROD10[i].carbon, i );
    tem.ag.setINITPROD10N( cohort[pichrt].initPROD10[i].nitrogen, i );
  }

  for( i = 0; i < 100; ++i )
  {
    tem.ag.setINITPROD100C( cohort[pichrt].initPROD100[i].carbon, i );
    tem.ag.setINITPROD100N( cohort[pichrt].initPROD100[i].nitrogen, i );
  }

  tem.ag.irrgflag = cohort[pichrt].irrgflag;

  tem.microbe.setKD( cohort[pichrt].kd );

  tem.ag.setNATPRVEETMX( cohort[pichrt].natprveetmx );

  tem.ag.setNATPRVLEAFMX( cohort[pichrt].natprvleafmx );

  tem.ag.setNATPRVPETMX( cohort[pichrt].natprvpetmx );

  tem.ag.setNATSEEDC( cohort[pichrt].natseedC );

  tem.ag.setNATSEEDSTRN( cohort[pichrt].natseedSTRN );

  tem.ag.setNATSEEDSTON( cohort[pichrt].natseedSTON );

  tem.ag.setNATSOIL( cohort[pichrt].natsoil );

  tem.ag.setNATTOPT( cohort[pichrt].nattopt );

  tem.ag.setNATYREET( cohort[pichrt].natyreet );

  tem.ag.setNATYRPET( cohort[pichrt].natyrpet );

  tem.veg.setNEWLEAFMX( cohort[pichrt].newleafmx );

  tem.veg.setNEWTOPT( cohort[pichrt].newtopt );

  tem.soil.setNSOLC( cohort[pichrt].nonsolc );

  tem.soil.setNSOLN( cohort[pichrt].nonsoln );

  tem.ag.setNRETENT( cohort[pichrt].nretent );

  tem.ag.setNSRETENT( cohort[pichrt].nsretent );

  tem.ag.setNVRETENT( cohort[pichrt].nvretent );

  tem.atms.setPETMX( cohort[pichrt].petmx );

  tem.atms.setPREV2TAIR( cohort[pichrt].prev2tair );

  tem.atms.setPREVCO2( cohort[pichrt].prevco2 );

  tem.ag.setPREVCROPRESIDUEC( cohort[pichrt].prevCropResidue.carbon );
  tem.ag.setPREVCROPRESIDUEN( cohort[pichrt].prevCropResidue.nitrogen );

  tem.soil.setPREVDST10( cohort[pichrt].prevdst10 );

  tem.ag.setPREVPROD1C( cohort[pichrt].prevPROD1.carbon );
  tem.ag.setPREVPROD1N( cohort[pichrt].prevPROD1.nitrogen );

  tem.ag.setPREVPROD10C( cohort[pichrt].prevPROD10.carbon );
  tem.ag.setPREVPROD10N( cohort[pichrt].prevPROD10.nitrogen );

  tem.ag.setPREVPROD100C( cohort[pichrt].prevPROD100.carbon );
  tem.ag.setPREVPROD100N( cohort[pichrt].prevPROD100.nitrogen );

  tem.soil.setPREVSPACK( cohort[pichrt].prevspack );

  tem.atms.setPREVTAIR( cohort[pichrt].prevtair );
  tem.atms.setPREVRAIN( cohort[pichrt].prevrain );

  tem.veg.setPREVUNRMLEAF( cohort[pichrt].prevunrmleaf );

  //tem.ag.setPROD10PAR( cohort[pichrt].prod10par );

  //tem.ag.setPROD100PAR( cohort[pichrt].prod100par );

  tem.ag.setPRODUCTYEAR( cohort[pichrt].productYear );

  tem.ag.setPRVCROPNPP( cohort[pichrt].prvcropnpp );

  tem.soil.setPRVEETMX( cohort[pichrt].prveetmx );

  tem.veg.setPRVLEAFMX( cohort[pichrt].prvleafmx );

  tem.atms.setPRVPETMX( cohort[pichrt].prvpetmx );

  //tem.ag.setSCONVERT( cohort[pichrt].sconvert );

  tem.ag.setSCONVRTFLXC( cohort[pichrt].sconvrtflx.carbon );
  tem.ag.setSCONVRTFLXN( cohort[pichrt].sconvrtflx.nitrogen );

  tem.ag.setSLASHC( cohort[pichrt].slash.carbon );
  tem.ag.setSLASHN( cohort[pichrt].slash.nitrogen );
  tem.ag.setdeadwoodc( cohort[pichrt].deadwood.carbon );
  tem.ag.setdeadwoodn( cohort[pichrt].deadwood.nitrogen );

  //tem.ag.setSLASHPAR( cohort[pichrt].slashpar );

  tem.soil.stm.setIS9( cohort[pichrt].STMis9 );

  tem.soil.stm.setSMASS9( cohort[pichrt].STMsmass9 );

  for( dnode = 0; dnode < MAXNODES; ++dnode )
  {
    tem.soil.stm.setDX9( cohort[pichrt].STMdx9[dnode], dnode );

    tem.soil.stm.setT9( cohort[pichrt].STMt9[dnode], dnode );

    tem.soil.stm.setWATER9( cohort[pichrt].STMwater9[dnode],
                            dnode );

    tem.soil.stm.setX9( cohort[pichrt].STMx9[dnode], dnode );

    tem.soil.stm.setXFA9( cohort[pichrt].STMxfa9[dnode],
                          dnode );

    tem.soil.stm.setXFB9( cohort[pichrt].STMxfb9[dnode],
                          dnode );
  }

  for( dnode = 0; dnode < MAXSNODES; ++dnode )
  {
    tem.soil.stm.setWEIGHT9( cohort[pichrt].STMweight9[dnode],
                            dnode );
  }

  tem.ag.tillflag = cohort[pichrt].tillflag;

  tem.veg.setTOPT( cohort[pichrt].topt );

  //tem.ag.setVCONVERT( cohort[pichrt].vconvert );

  tem.ag.setVCONVRTFLXC( cohort[pichrt].vconvrtflx.carbon );
  tem.ag.setVCONVRTFLXN( cohort[pichrt].vconvrtflx.nitrogen );

  tem.ag.setVRESPAR( 1-cohort[pichrt].stemmortpar );
  //tem.ag.setVRESPAR( cohort[pichrt].vrespar );

  tem.soil.setWFPSOFF( cohort[pichrt].wfpsoff );

  tem.veg.yrltrc = cohort[pichrt].yrltrc;
  tem.veg.yrltrn = cohort[pichrt].yrltrn;
  tem.ag.leafmortpar = cohort[pichrt].leafmortpar;
  tem.ag.rootmortpar = cohort[pichrt].rootmortpar;
  tem.ag.stemmortpar = cohort[pichrt].stemmortpar;

  tem.ag.leafslash = cohort[pichrt].leafslash;
  tem.ag.leafconv = cohort[pichrt].leafconv;
  tem.ag.rootslash = cohort[pichrt].rootslash;
  tem.ag.rootconv = cohort[pichrt].rootconv;
  tem.ag.stemslash = cohort[pichrt].stemslash;
  tem.ag.stemconv = cohort[pichrt].stemconv;
  tem.ag.standdead = cohort[pichrt].standdead;
  tem.ag.deadconv = cohort[pichrt].deadconv;
  tem.ag.deadslash = cohort[pichrt].deadslash;
  tem.ag.sconvert = cohort[pichrt].sconvert;
  tem.ag.prod10par = cohort[pichrt].prod10par;
  tem.ag.prod100par = cohort[pichrt].prod100par;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMelmnt60::readCohortState( ifstream& ifstate,
                                  const int& pichrt )
{
  int i;
  int dm;
  int dnode;
  double dumflt1;
  double dumflt2;
  int dumint;
  string dumstr;

ifstate >> dumflt1;    // Longitude of element
ifstate >> dumflt2;    // Latitude of element

ifstate >> dumint;   // ichrt+1
//cout <<"lon: "<<dumflt1<<" lat: "<<dumflt2<<" col: "<<col<<"row: "<<row<<endl;
ifstate >> cohort[pichrt].srcCohort;
ifstate >> cohort[pichrt].standage;
ifstate >> cohort[pichrt].chrtarea;
ifstate >> cohort[pichrt].potveg;
ifstate >> cohort[pichrt].currentveg;
ifstate >> cohort[pichrt].subtype;
ifstate >> cohort[pichrt].cmnt;

//changed by cgs2014 to read different format initial cohort state data
for( i = 0; i < MAXSTATE; ++i )
//for( i = 0; i < 12; ++i )
{
  ifstate >> cohort[pichrt].y[i];
  ifstate >> cohort[pichrt].prevy[i];
  //if (i<=14) cout <<" y[" <<i <<"]: " <<cohort[pichrt].y[i]<<" ";
}

//if (cohort[pichrt].cmnt==7 && cohort[pichrt].currentveg==50) cout <<" col: "<<dumflt1<<" row: "<<dumflt2<<" nh4: "<<cohort[pichrt].y[7]<<" no3: "<<cohort[pichrt].y[8]<<endl;
//if (col == -123.0 && row== 44.25 && cohort[pichrt].currentveg==13) cout <<"cohort: "<< dumint<<" vegc: "<<cohort[pichrt].y[0]<<" soilc: "<<cohort[pichrt].y[1]<<endl;
//cohort[pichrt].y[24]=1000; //initial forestage
//cohort[pichrt].prevy[24]=1000;

//if (cohort[pichrt].y[13]!=cohort[pichrt].y[13]) cout <<"y13 has no data " <<" cmnt: " <<cohort[pichrt].cmnt<<" icohort: " <<dumint<<" y[1]: " <<cohort[pichrt].y[1]<<" y0: " <<cohort[pichrt].y[0]<<endl;

//if (col == -129.5 && (row ==55))
//cout <<"col: " <<col<<" row: " <<row<<" ichrt: " <<pichrt<<" src: " <<cohort[pichrt].srcCohort <<" standage: " <<cohort[pichrt].standage<<" chrtarea: " <<cohort[pichrt].chrtarea<<" cmnt: " <<cohort[pichrt].cmnt<<" y0: " <<cohort[pichrt].y[0]<< " y1: " <<cohort[pichrt].y[1]<<endl;
ifstate >> cohort[pichrt].agcmnt;

ifstate >> cohort[pichrt].aggrowdd;

ifstate >> cohort[pichrt].agkd;

ifstate >> cohort[pichrt].agprvstate;

ifstate >> cohort[pichrt].agstate;

ifstate >> cohort[pichrt].c2n;

ifstate >> cohort[pichrt].cneven;

ifstate >> cohort[pichrt].convrtflx.carbon;
ifstate >> cohort[pichrt].convrtflx.nitrogen;

ifstate >> cohort[pichrt].cropprveetmx;

ifstate >> cohort[pichrt].cropprvleafmx;

ifstate >> cohort[pichrt].cropprvpetmx;

ifstate >> cohort[pichrt].cropResidue.carbon;
ifstate >> cohort[pichrt].cropResidue.nitrogen;

ifstate >> cohort[pichrt].croptopt;

ifstate >> cohort[pichrt].distmnthcnt;

ifstate >> cohort[pichrt].disturbflag;

ifstate >> cohort[pichrt].disturbmonth;

for( dm = 0; dm < CYCLE; ++dm )
{
  ifstate >> cohort[pichrt].dst10[dm];
}

ifstate >> cohort[pichrt].eetmx;

ifstate >> cohort[pichrt].fertflag;

ifstate >> cohort[pichrt].firemnthcnt;

ifstate >> cohort[pichrt].firendep;

ifstate >> cohort[pichrt].formPROD10.carbon;
ifstate >> cohort[pichrt].formPROD10.nitrogen;

ifstate >> cohort[pichrt].formPROD100.carbon;
ifstate >> cohort[pichrt].formPROD100.nitrogen;

ifstate >> cohort[pichrt].fprevozone;

ifstate >> cohort[pichrt].FRI;

for( dm = 0; dm < CYCLE; ++dm )
{
  ifstate >> cohort[pichrt].initPROD1[dm].carbon;
  ifstate >> cohort[pichrt].initPROD1[dm].nitrogen;
}

for( i = 0; i < 10; ++i )
{
  ifstate >> cohort[pichrt].initPROD10[i].carbon;
  ifstate >> cohort[pichrt].initPROD10[i].nitrogen;
}

for( i = 0; i < 100; ++i )
{
  ifstate >> cohort[pichrt].initPROD100[i].carbon;
  ifstate >> cohort[pichrt].initPROD100[i].nitrogen;
}

ifstate >> cohort[pichrt].irrgflag;

ifstate >> cohort[pichrt].kd;

//if (cohort[pichrt].kd!=cohort[pichrt].kd) cout <<"kd has no data: " <<endl;

ifstate >> cohort[pichrt].natprveetmx;
//cout <<"natprveetmx: "<<cohort[pichrt].natprveetmx<<endl;

ifstate >> cohort[pichrt].natprvleafmx;

ifstate >> cohort[pichrt].natprvpetmx;

ifstate >> cohort[pichrt].natseedC;

ifstate >> cohort[pichrt].natseedSTRN;

ifstate >> cohort[pichrt].natseedSTON;

ifstate >> cohort[pichrt].natsoil;

ifstate >> cohort[pichrt].nattopt;

ifstate >> cohort[pichrt].natyreet;

ifstate >> cohort[pichrt].natyrpet;

ifstate >> cohort[pichrt].newleafmx;

ifstate >> cohort[pichrt].newtopt;

ifstate >> cohort[pichrt].nonsolc;

ifstate >> cohort[pichrt].nonsoln;

ifstate >> cohort[pichrt].nretent;

ifstate >> cohort[pichrt].nsretent;

ifstate >> cohort[pichrt].nvretent;

ifstate >> cohort[pichrt].petmx;

ifstate >> cohort[pichrt].prev2tair;

ifstate >> cohort[pichrt].prevco2;

ifstate >> cohort[pichrt].prevCropResidue.carbon;
ifstate >> cohort[pichrt].prevCropResidue.nitrogen;

ifstate >> cohort[pichrt].prevdst10;

ifstate >> cohort[pichrt].prevPROD1.carbon;
ifstate >> cohort[pichrt].prevPROD1.nitrogen;

ifstate >> cohort[pichrt].prevPROD10.carbon;
ifstate >> cohort[pichrt].prevPROD10.nitrogen;

ifstate >> cohort[pichrt].prevPROD100.carbon;
ifstate >> cohort[pichrt].prevPROD100.nitrogen;

ifstate >> cohort[pichrt].prevspack;

ifstate >> cohort[pichrt].prevtair;
ifstate >> cohort[pichrt].prevrain;

ifstate >> cohort[pichrt].prevunrmleaf;

//ifstate >> cohort[pichrt].prod10par;

//ifstate >> cohort[pichrt].prod100par;

ifstate >> cohort[pichrt].productYear;

ifstate >> cohort[pichrt].prvchrtarea;

ifstate >> cohort[pichrt].prvcropnpp;

ifstate >> cohort[pichrt].prveetmx;

ifstate >> cohort[pichrt].prvleafmx;

ifstate >> cohort[pichrt].prvpetmx;

ifstate >> cohort[pichrt].qc;

//ifstate >> cohort[pichrt].sconvert;

ifstate >> cohort[pichrt].sconvrtflx.carbon;
ifstate >> cohort[pichrt].sconvrtflx.nitrogen;

ifstate >> cohort[pichrt].slash.carbon;
ifstate >> cohort[pichrt].slash.nitrogen;

//ifstate >> cohort[pichrt].slashpar;

ifstate >> cohort[pichrt].STMis9;

ifstate >> cohort[pichrt].STMsmass9;

for( dnode = 0; dnode < MAXNODES; ++dnode )
{
  ifstate >> cohort[pichrt].STMdx9[dnode];
  ifstate >> cohort[pichrt].STMt9[dnode];
  ifstate >> cohort[pichrt].STMwater9[dnode];
  ifstate >> cohort[pichrt].STMx9[dnode];
  ifstate >> cohort[pichrt].STMxfa9[dnode];
  ifstate >> cohort[pichrt].STMxfb9[dnode];
}

for( dnode = 0; dnode < MAXSNODES; ++dnode )
{
  ifstate >> cohort[pichrt].STMweight9[dnode];
}

ifstate >> cohort[pichrt].tillflag;

ifstate >> cohort[pichrt].topt;

ifstate >> cohort[pichrt].tqc;

//ifstate >> cohort[pichrt].vconvert;

ifstate >> cohort[pichrt].vconvrtflx.carbon;
ifstate >> cohort[pichrt].vconvrtflx.nitrogen;

ifstate >> cohort[pichrt].vrespar;

ifstate >> cohort[pichrt].wfpsoff;

ifstate >> cohort[pichrt].yrltrc;
ifstate >> cohort[pichrt].yrltrn;
  /*
  else
  {
	  string tempstring;
	  for (int i=0;i<1646;i++) //total: 1648.
	  {
		ifstate>>tempstring;
		//cout<<tempstring<<" ";
	  }
	  //cout <<"end: "<<endl;

  }
  */
  //cout <<" yrtrn0: "<<cohort[pichrt].yrltrn<<" yrltrc: "<<cohort[pichrt].yrltrc <<endl;


ifstate >> cohort[pichrt].leafmortpar;
ifstate >> cohort[pichrt].rootmortpar;
ifstate >> cohort[pichrt].stemmortpar;
ifstate >> cohort[pichrt].leafslash;
ifstate >> cohort[pichrt].leafconv;
ifstate >> cohort[pichrt].rootslash;
ifstate >> cohort[pichrt].rootconv;
ifstate >> cohort[pichrt].stemslash;
ifstate >> cohort[pichrt].stemconv;
ifstate >> cohort[pichrt].standdead;
ifstate >> cohort[pichrt].deadconv;
ifstate >> cohort[pichrt].deadslash;
ifstate >> cohort[pichrt].sconvert;
ifstate >> cohort[pichrt].prod10par;
ifstate >> cohort[pichrt].prod100par;
//if (col==-127.75 && row== 69.25) cout <<"col2: "<<dumflt2 <<" row: "<<dumflt1<<" cohort: "<<pichrt<<" chrtarea: "<<cohort[pichrt].chrtarea<< " qc: "<<cohort[pichrt].qc<< " tqc: "<<cohort[pichrt].tqc<<endl;


ifstate.seekg( 0, ios::cur );

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt60::saveTEMCohortState( const int& pichrt )
{
  int dm;
  int dnode;
  int i;


  cohort[pichrt].potveg = tem.veg.getPOTVEG();

  cohort[pichrt].currentveg = tem.veg.getCURRENTVEG();

  cohort[pichrt].subtype = tem.veg.getSUBTYPE();

  cohort[pichrt].cmnt = tem.veg.cmnt;

  for( i = 0; i < MAXSTATE; ++i )
  {
    cohort[pichrt].y[i] = tem.getY( i );
    cohort[pichrt].prevy[i] = tem.getPREVY( i );
    if (cohort[pichrt].y[i] >= -0.00001 && cohort[pichrt].y[i] <= 0.00001) cohort[pichrt].y[i] = 0.0;
    if (cohort[pichrt].prevy[i]>= -0.00001 && cohort[pichrt].y[i] <= 0.00001) cohort[pichrt].prevy[i] = 0.0;
  }
  //if (col == -80.75 && row== 51.25 && tem.veg.cmnt==4) cout <<"vegc_save: "<<cohort[pichrt].y[0]<<" prevegc: "<<cohort[pichrt].prevy[0]<<endl;
  cohort[pichrt].agcmnt = tem.ag.cmnt;

  cohort[pichrt].aggrowdd = tem.ag.getGROWDD();

  cohort[pichrt].agkd = tem.ag.getKD();

  cohort[pichrt].agprvstate = tem.ag.prvstate;

  cohort[pichrt].agstate = tem.ag.state;

  cohort[pichrt].c2n = tem.veg.getC2N();

  cohort[pichrt].cneven = tem.veg.getCNEVEN();

  cohort[pichrt].convrtflx.carbon = tem.ag.getCONVRTFLXC();
  cohort[pichrt].convrtflx.nitrogen = tem.ag.getCONVRTFLXN();

  cohort[pichrt].cropprveetmx = tem.ag.getCROPPRVEETMX();

  cohort[pichrt].cropprvleafmx = tem.ag.getCROPPRVLEAFMX();

  cohort[pichrt].cropprvpetmx = tem.ag.getCROPPRVPETMX();

  cohort[pichrt].cropResidue.carbon = tem.ag.getCROPRESIDUEC();
  cohort[pichrt].cropResidue.nitrogen = tem.ag.getCROPRESIDUEN();

  cohort[pichrt].croptopt = tem.ag.getCROPTOPT();

  cohort[pichrt].distmnthcnt = tem.distmnthcnt;

  cohort[pichrt].disturbflag = tem.disturbflag;

  cohort[pichrt].disturbmonth = tem.disturbmonth;

  cohort[pichrt].eetmx = tem.soil.getEETMX();

  cohort[pichrt].fertflag = tem.ag.fertflag;

  cohort[pichrt].firemnthcnt = tem.firemnthcnt;

  cohort[pichrt].firendep = tem.ag.getFIRENDEP();

  cohort[pichrt].formPROD10.carbon = tem.ag.getFORMPROD10C();
  cohort[pichrt].formPROD10.nitrogen = tem.ag.getFORMPROD10N();

  cohort[pichrt].formPROD100.carbon = tem.ag.getFORMPROD100C();
  cohort[pichrt].formPROD100.nitrogen = tem.ag.getFORMPROD100N();

  cohort[pichrt].fprevozone = tem.veg.getFPREVOZONE();

  cohort[pichrt].FRI = tem.ag.getFRI();

  for( dm = 0; dm < CYCLE; ++dm )
  {
    cohort[pichrt].initPROD1[dm].carbon = tem.ag.getINITPROD1C( dm );
    cohort[pichrt].initPROD1[dm].nitrogen = tem.ag.getINITPROD1N( dm );
  }

  for( i = 0; i < 10; ++i )
  {
    cohort[pichrt].initPROD10[i].carbon = tem.ag.getINITPROD10C( i );
    cohort[pichrt].initPROD10[i].nitrogen = tem.ag.getINITPROD10N( i );
  }

  for( i = 0; i < 100; ++i )
  {
    cohort[pichrt].initPROD100[i].carbon = tem.ag.getINITPROD100C( i );
    cohort[pichrt].initPROD100[i].nitrogen = tem.ag.getINITPROD100N( i );
  }

  cohort[pichrt].irrgflag = tem.ag.irrgflag;

  cohort[pichrt].kd = tem.microbe.getKD();

  cohort[pichrt].natprveetmx = tem.ag.getNATPRVEETMX();

  cohort[pichrt].natprvleafmx = tem.ag.getNATPRVLEAFMX();

  cohort[pichrt].natprvpetmx = tem.ag.getNATPRVPETMX();

  cohort[pichrt].natseedC = tem.ag.getNATSEEDC();

  cohort[pichrt].natseedSTRN = tem.ag.getNATSEEDSTRN();

  cohort[pichrt].natseedSTON = tem.ag.getNATSEEDSTON();

  cohort[pichrt].natsoil = tem.ag.getNATSOIL();

  cohort[pichrt].nattopt = tem.ag.getNATTOPT();

  cohort[pichrt].natyreet = tem.ag.getNATYREET();

  cohort[pichrt].natyrpet = tem.ag.getNATYRPET();

  cohort[pichrt].newleafmx = tem.veg.getNEWLEAFMX();

  cohort[pichrt].newtopt = tem.veg.getNEWTOPT();

  cohort[pichrt].nonsolc = tem.soil.getNSOLC();

  cohort[pichrt].nonsoln = tem.soil.getNSOLN();

  cohort[pichrt].nretent = tem.ag.getNRETENT();

  cohort[pichrt].nsretent = tem.ag.getNSRETENT();

  cohort[pichrt].nvretent = tem.ag.getNVRETENT();

  cohort[pichrt].petmx = tem.atms.getPETMX();

  cohort[pichrt].prev2tair = tem.atms.getPREV2TAIR();

  cohort[pichrt].prevco2 = tem.atms.getPREVCO2();

  cohort[pichrt].prevCropResidue.carbon = tem.ag.getPREVCROPRESIDUEC();
  cohort[pichrt].prevCropResidue.nitrogen = tem.ag.getPREVCROPRESIDUEN();

  cohort[pichrt].prevdst10 = tem.soil.getPREVDST10();

  cohort[pichrt].prevPROD1.carbon = tem.ag.getPREVPROD1C();
  cohort[pichrt].prevPROD1.nitrogen = tem.ag.getPREVPROD1N();

  cohort[pichrt].prevPROD10.carbon = tem.ag.getPREVPROD10C();
  cohort[pichrt].prevPROD10.nitrogen = tem.ag.getPREVPROD10N();

  cohort[pichrt].prevPROD100.carbon = tem.ag.getPREVPROD100C();
  cohort[pichrt].prevPROD100.nitrogen = tem.ag.getPREVPROD100N();

  cohort[pichrt].prevspack = tem.soil.getPREVSPACK();

  cohort[pichrt].prevtair = tem.atms.getPREVTAIR();
  cohort[pichrt].prevrain = tem.atms.getPREVRAIN();

  cohort[pichrt].prevunrmleaf = tem.veg.getPREVUNRMLEAF();

  //cohort[pichrt].prod10par = tem.ag.getPROD10PAR();

  //cohort[pichrt].prod100par = tem.ag.getPROD100PAR();

  cohort[pichrt].productYear = tem.ag.getPRODUCTYEAR();

  cohort[pichrt].prvcropnpp = tem.ag.getPRVCROPNPP();

  cohort[pichrt].prveetmx = tem.soil.getPRVEETMX();

  cohort[pichrt].prvleafmx = tem.veg.getPRVLEAFMX();

  cohort[pichrt].prvpetmx = tem.atms.getPRVPETMX();

  //cohort[pichrt].sconvert = tem.ag.getSCONVERT();

  cohort[pichrt].sconvrtflx.carbon = tem.ag.getSCONVRTFLXC();
  cohort[pichrt].sconvrtflx.nitrogen = tem.ag.getSCONVRTFLXN();

  cohort[pichrt].slash.carbon = tem.ag.getSLASHC();
  cohort[pichrt].slash.nitrogen = tem.ag.getSLASHN();
  cohort[pichrt].deadwood.carbon = tem.ag.getdeadwoodc();
  cohort[pichrt].deadwood.nitrogen = tem.ag.getdeadwoodn();

  //cohort[pichrt].slashpar = tem.ag.getSLASHPAR();

  cohort[pichrt].STMis9 = tem.soil.stm.getIS9();

  cohort[pichrt].STMsmass9 = tem.soil.stm.getSMASS9();

  for( dnode = 0; dnode < MAXNODES; ++dnode )
  {
    cohort[pichrt].STMdx9[dnode] = tem.soil.stm.getDX9( dnode );

    cohort[pichrt].STMt9[dnode] = tem.soil.stm.getT9( dnode );

    cohort[pichrt].STMwater9[dnode] = tem.soil.stm.getWATER9( dnode );

    cohort[pichrt].STMx9[dnode] = tem.soil.stm.getX9( dnode );

    cohort[pichrt].STMxfa9[dnode] = tem.soil.stm.getXFA9( dnode );

    cohort[pichrt].STMxfb9[dnode] = tem.soil.stm.getXFB9( dnode );
  }

  for( dnode = 0; dnode < MAXSNODES; ++dnode )
  {
    cohort[pichrt].STMweight9[dnode] = tem.soil.stm.getWEIGHT9( dnode );
  }

  cohort[pichrt].tillflag = tem.ag.tillflag;

  cohort[pichrt].topt = tem.veg.getTOPT();

  cohort[pichrt].vconvert = tem.ag.getVCONVERT();

  cohort[pichrt].vconvrtflx.carbon = tem.ag.getVCONVRTFLXC();
  cohort[pichrt].vconvrtflx.nitrogen = tem.ag.getVCONVRTFLXN();

  cohort[pichrt].vrespar = tem.ag.getVRESPAR();

  cohort[pichrt].wfpsoff = tem.soil.getWFPSOFF();

  cohort[pichrt].yrltrc = tem.veg.yrltrc;
  cohort[pichrt].yrltrn = tem.veg.yrltrn;

  cohort[pichrt].leafmortpar =  tem.ag.leafmortpar;
  cohort[pichrt].rootmortpar = tem.ag.rootmortpar;
  cohort[pichrt].stemmortpar = tem.ag.stemmortpar;
  cohort[pichrt].leafslash =  tem.ag.leafslash;
  cohort[pichrt].leafconv =  tem.ag.leafconv;
  cohort[pichrt].rootslash =  tem.ag.rootslash;
  cohort[pichrt].rootconv =  tem.ag.rootconv;
  cohort[pichrt].stemslash =  tem.ag.stemslash;
  cohort[pichrt].stemconv =  tem.ag.stemconv;
  cohort[pichrt].standdead =  tem.ag.standdead;
  cohort[pichrt].deadconv =  tem.ag.deadconv;
  cohort[pichrt].deadslash =  tem.ag.deadslash;
  cohort[pichrt].sconvert =  tem.ag.sconvert;
  cohort[pichrt].prod10par =  tem.ag.prod10par;
  cohort[pichrt].prod100par =  tem.ag.prod100par;
};

/* *************************************************************
************************************************************** */

/* *************************************************************
************************************************************** */

void TEMelmnt60::setCohortTEMState( const ElmntCohort60& firstchrt,
                                    ElmntCohort60& targetchrt )
{
  int dm;
  int dnode;
  int i;


  for( i = 0; i < MAXSTATE; ++i )
  {
    targetchrt.y[i] = firstchrt.y[i];
    targetchrt.prevy[i] = firstchrt.prevy[i];
  }

  targetchrt.aggrowdd = firstchrt.aggrowdd;

  targetchrt.agkd = firstchrt.agkd;

  targetchrt.c2n = firstchrt.c2n;

  targetchrt.cneven = firstchrt.cneven;

  targetchrt.convrtflx.carbon = firstchrt.convrtflx.carbon;
  targetchrt.convrtflx.nitrogen = firstchrt.convrtflx.nitrogen;

  targetchrt.cropprveetmx = firstchrt.cropprveetmx;

  targetchrt.cropprvleafmx = firstchrt.cropprvleafmx;

  targetchrt.cropprvpetmx = firstchrt.cropprvpetmx;

  targetchrt.cropResidue.carbon = firstchrt.cropResidue.carbon;
  targetchrt.cropResidue.nitrogen = firstchrt.cropResidue.nitrogen;

  targetchrt.croptopt = firstchrt.croptopt;

  targetchrt.distmnthcnt = firstchrt.distmnthcnt;

  for( dm = 0; dm < CYCLE; ++dm )
  {
    targetchrt.dst10[dm] = firstchrt.dst10[dm];
  }

  targetchrt.eetmx = firstchrt.eetmx;

  targetchrt.firemnthcnt = firstchrt.firemnthcnt;

  targetchrt.firendep = firstchrt.firendep;

  targetchrt.formPROD10.carbon = firstchrt.formPROD10.carbon;
  targetchrt.formPROD10.nitrogen = firstchrt.formPROD10.nitrogen;

  targetchrt.formPROD100.carbon = firstchrt.formPROD100.carbon;
  targetchrt.formPROD100.nitrogen = firstchrt.formPROD100.nitrogen;

  targetchrt.fprevozone = firstchrt.fprevozone;

  for( dm = 0; dm < CYCLE; ++dm )
  {
    targetchrt.initPROD1[dm].carbon = firstchrt.initPROD1[dm].carbon;
    targetchrt.initPROD1[dm].nitrogen = firstchrt.initPROD1[dm].nitrogen;
  }

  for( i = 0; i < 10; ++i )
  {
    targetchrt.initPROD10[i].carbon = firstchrt.initPROD10[i].carbon;
    targetchrt.initPROD10[i].nitrogen = firstchrt.initPROD10[i].nitrogen;
  }

  for( i = 0; i < 100; ++i )
  {
    targetchrt.initPROD100[i].carbon = firstchrt.initPROD100[i].carbon;
    targetchrt.initPROD100[i].nitrogen = firstchrt.initPROD100[i].nitrogen;
  }

  targetchrt.kd = firstchrt.kd;

  targetchrt.natprveetmx = firstchrt.natprveetmx;

  targetchrt.natprvleafmx = firstchrt.natprvleafmx;

  targetchrt.natprvpetmx = firstchrt.natprvpetmx;

  targetchrt.natseedC = firstchrt.natseedC;

  targetchrt.natseedSTRN = firstchrt.natseedSTRN;

  targetchrt.natseedSTON = firstchrt.natseedSTON;

  targetchrt.natsoil = firstchrt.natsoil;

  targetchrt.nattopt = firstchrt.nattopt;

  targetchrt.natyreet = firstchrt.natyreet;

  targetchrt.natyrpet = firstchrt.natyrpet;

  targetchrt.newleafmx =  firstchrt.newleafmx;

  targetchrt.newtopt = firstchrt.newtopt;

  targetchrt.nonsolc = firstchrt.nonsolc;

  targetchrt.nonsoln = firstchrt.nonsoln;

  targetchrt.nretent = firstchrt.nretent;

  targetchrt.nsretent = firstchrt.nsretent;

  targetchrt.nvretent = firstchrt.nvretent;

  targetchrt.petmx = firstchrt.petmx;

  targetchrt.prev2tair = firstchrt.prev2tair;

  targetchrt.prevco2 = firstchrt.prevco2;

  targetchrt.prevCropResidue.carbon = firstchrt.prevCropResidue.carbon;
  targetchrt.prevCropResidue.nitrogen = firstchrt.prevCropResidue.nitrogen;

  targetchrt.prevdst10 = firstchrt.prevdst10;

  targetchrt.prevPROD1.carbon = firstchrt.prevPROD1.carbon;
  targetchrt.prevPROD1.nitrogen = firstchrt.prevPROD1.nitrogen;

  targetchrt.prevPROD10.carbon = firstchrt.prevPROD10.carbon;
  targetchrt.prevPROD10.nitrogen = firstchrt.prevPROD10.nitrogen;

  targetchrt.prevPROD100.carbon = firstchrt.prevPROD100.carbon;
  targetchrt.prevPROD100.nitrogen = firstchrt.prevPROD100.nitrogen;

  targetchrt.prevspack = firstchrt.prevspack;

  targetchrt.prevtair = firstchrt.prevtair;
  targetchrt.prevrain = firstchrt.prevrain;

  targetchrt.prevunrmleaf = firstchrt.prevunrmleaf;

  targetchrt.productYear = firstchrt.productYear;

  targetchrt.prvcropnpp = firstchrt.prvcropnpp;

  targetchrt.prveetmx = firstchrt.prveetmx;

  targetchrt.prvleafmx = firstchrt.prvleafmx;

  targetchrt.prvpetmx = firstchrt.prvpetmx;

  targetchrt.qc = firstchrt.qc;

  //targetchrt.sconvrtflx.carbon = firstchrt.sconvrtflx.carbon;
  //targetchrt.sconvrtflx.nitrogen = firstchrt.sconvrtflx.nitrogen;

  targetchrt.slash.carbon = firstchrt.slash.carbon;
  targetchrt.slash.nitrogen = firstchrt.slash.nitrogen;
  //targetchrt.deadwood.carbon = firstchrt.deadwood.carbon;
  //targetchrt.deadwood.nitrogen = firstchrt.deadwood.nitrogen;

  targetchrt.STMis9 = firstchrt.STMis9;

  targetchrt.STMsmass9 = firstchrt.STMsmass9;

  for( dnode = 0; dnode < MAXNODES; ++dnode )
  {
    targetchrt.STMdx9[dnode] = firstchrt.STMdx9[dnode];

    targetchrt.STMt9[dnode] = firstchrt.STMt9[dnode];

    targetchrt.STMwater9[dnode] = firstchrt.STMwater9[dnode];

    targetchrt.STMx9[dnode] = firstchrt.STMx9[dnode];

    targetchrt.STMxfa9[dnode] = firstchrt.STMxfa9[dnode];

    targetchrt.STMxfb9[dnode] = firstchrt.STMxfb9[dnode];
  }

  for( dnode = 0; dnode < MAXSNODES; ++dnode )
  {
    targetchrt.STMweight9[dnode] = firstchrt.STMweight9[dnode];
  }

  targetchrt.topt = firstchrt.topt;

  targetchrt.tqc = firstchrt.tqc;

  //targetchrt.vconvrtflx.carbon = firstchrt.vconvrtflx.carbon;
  //targetchrt.vconvrtflx.nitrogen = firstchrt.vconvrtflx.nitrogen;

  targetchrt.wfpsoff = firstchrt.wfpsoff;

  targetchrt.yrltrc = firstchrt.yrltrc;
  targetchrt.yrltrn = firstchrt.yrltrn;

  /*
  targetchrt.leafmortpar =  firstchrt.leafmortpar;
  targetchrt.stemmortpar = firstchrt.stemmortpar;
  targetchrt.rootmortpar = firstchrt.rootmortpar;
  targetchrt.leafslash =  firstchrt.leafslash;
  targetchrt.leafconv =  firstchrt.leafconv;
  targetchrt.rootslash =  firstchrt.rootslash;
  targetchrt.rootconv =  firstchrt.rootconv;
  targetchrt.stemslash =  firstchrt.stemslash;
  targetchrt.stemconv =  firstchrt.stemconv;
  targetchrt.standdead =  firstchrt.standdead;
  targetchrt.deadconv =  firstchrt.deadconv;
  targetchrt.deadslash =  firstchrt.deadslash;
  targetchrt.sconvert =  firstchrt.sconvert;
  targetchrt.prod10par =  firstchrt.prod10par;
  targetchrt.prod100par =  firstchrt.prod100par;
  */


};

/* *************************************************************
************************************************************** */


