/* *************************************************************
****************************************************************
TCLM602.CPP - describes climate module used by TEM

20030513 - DWK created by modifying tclm425.cpp to incorporate
           open-N cycle dynamics into TEM
20030513 - DWK changed TEMclm43:: to TEMclm51::
20030513 - DWK added initNdep() and setNdepFlags()
20030513 - DWK added NDEP to constructor
20040714 - DWK changed TEMclm51:: to TEMclm60::
20040716 - DWK changed inheritance of Atmosphere51 to
           inheritance of Atmosphere60
20051117 - DWK added include tclm602.h and standard includes
20051117 - DWK changed CO2data co2dat to CO2data60 co2dat in 
           loadyrCO2()
20051121 - DWK added switch() to initCO2()
20051121 - DWK added I_RAIN and I_SNWFAL to TEMclm60()

****************************************************************
************************************************************* */

#include<cstdio>

  using std::FILE;
  using std::fscanf;

#include<iostream>

  using std::cin;
  using std::cout;
  using std::ios;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setprecision;

#include<cstdlib>
  
  using std::exit;
  
#include<cmath>

  using std::cos;
  using std::pow;
  using std::sin;

#include<vector>
  using std::vector;
    
#include<string>

  using std::string;


#include "tclm602.h"

/* *************************************************************
************************************************************* */

TEMclm60::TEMclm60() : Atmosphere60(), predstr( NUMATMS )
{

// Initialize predstr array

  predstr.at( I_GIRR ) = "GIRR";
  predstr.at( I_NIRR ) = "NIRR";
  predstr.at( I_PAR ) = "PAR";
  predstr.at( I_CLDS ) = "CLDINESS";
  predstr.at( I_TAIR ) = "TAIR";
  predstr.at( I_PREC ) = "PREC";
  predstr.at( I_RAIN ) = "RAIN";
  predstr.at( I_SNWFAL ) = "SNOWFALL";
  predstr.at( I_CO2 ) = "CO2";
  predstr.at( I_AOT40 ) = "AOT40";
  predstr.at( I_TNDEP ) = "TOTNDEP";
  predstr.at( I_NH4DEP ) = "NH4DEP";
  predstr.at( I_NO3DEP ) = "NO3DEP";

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMclm60::initCO2( ofstream& rflog1 )
{
  switch( tco2flag )
  {
    case 1:
		#ifdef PMODE
    		*fgo >> ico2fname;
    	#else
    		cout << "Please enter the file name containing ";
            cout << "the annual transient atmospheric CO2 ";
            cout << "data:" << endl;
            cout << "               (e.g., CO2) " << endl;

            cin >> ico2fname;
		#endif
            rflog1 << "Please enter the file name containing ";
            rflog1 << "the annual transient atmospheric CO2 ";
            rflog1 << "data: " << endl;
            rflog1 << "               (e.g., CO2) " << endl;
            rflog1 << ico2fname << endl << endl;
            break;

    case 2:
		#ifdef PMODE
    		*fgo >> ico2fname;
    		*fgo >> ico2end;
		#else

    		cout << "Please enter the first part of the file ";
            cout << "name containing the monthly transient ";
            cout << "atmospheric CO2 data: " << endl;
            cout << "               (e.g., CO2) " << endl;

            cin >> ico2fname;

            cout << "Please enter the file extension ";
            cout << "(include " "the '.'): " << endl;

            cin >> ico2end;
		#endif
            rflog1 << "Please enter the file name containing the monthly";
            rflog1 << " transient atmospheric CO2 data: " << endl;
            rflog1 << "               (e.g., CO2) " << endl;
            rflog1 << ico2fname << endl << endl;
            
            rflog1 << "Please enter the file extension ";
            rflog1 << "(include the '.'): " << endl;
            rflog1 << ico2end;
            
            break;
            
  }

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void TEMclm60::initNdep( ofstream& rflog1 )
{

  if( 1 == tndepflag )
  {
#ifdef PMODE
	*fgo >> indepfname;
	*fgo >> independ;
#else
    cout << "Please enter the first part of the file name";
    cout << " containing the monthly N deposition data: " << endl;
    cout << "               (e.g., NDEP) " << endl;

    cin >> indepfname;

    cout << "Please enter the file extension (include the '.'): ";

    cin >> independ;
#endif
    rflog1 << "Please enter the first part of the file name";
    rflog1 << " containing the monthly N deposition data: " << endl;
    rflog1 << "               (e.g., NDEP) " << endl;
    rflog1 << indepfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << independ << endl << endl;
  }
  else
  {

#ifdef PMODE
	*fgo >> indepfname;

#else

	cout << "Please enter the file name containing the";
    cout << " long-term mean monthly N deposition data: " << endl;
    cout << "               (e.g., NDEP) " << endl;

    cin >> indepfname;
#endif

    rflog1 << "Please enter the file name containing the";
    rflog1 << " long-term mean monthly N deposition data: " << endl;
    rflog1 << "               (e.g., NDEP) " << endl;
    rflog1 << indepfname << endl << endl;
  }

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void TEMclm60::initO3( ofstream& rflog1 )
{

  if( 1 == to3flag )
  {
#ifdef PMODE
	*fgo >> io3fname;
	*fgo >> io3end;
#else

    cout << "Please enter the first part of the file name";
    cout << " containing the monthly O3 data: " << endl;
    cout << "               (e.g., D40) " << endl;

    cin >> io3fname;

    cout << "Please enter the file extension (include the '.'): ";

    cin >> io3end;
#endif

    rflog1 << "Please enter the first part of the file name";
    rflog1 << " containing the monthly O3 data: " << endl;
    rflog1 << "               (e.g., D40) " << endl;
    rflog1 << io3fname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << io3end << endl << endl;
  }
  else
  {
#ifdef PMODE
	*fgo >> io3fname;
#else
    cout << "Please enter the file name containing the";
    cout << " long-term mean monthly O3 data: " << endl;
    cout << "               (e.g., D40) " << endl;

    cin >> io3fname;
#endif
    rflog1 << "Please enter the file name containing the";
    rflog1 << " long-term mean monthly O3 data: " << endl;
    rflog1 << "               (e.g., D40) " << endl;
    rflog1 << io3fname << endl << endl;
  }

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void TEMclm60::initPrec( ofstream& rflog1 )
{

  if( 1 == tprecflag )
  {
#ifdef PMODE
	*fgo >> iprecfname;
	*fgo >> iprecend;
#else
    cout << "Please enter the first part of the file name";
    cout << " containing the monthly precipitation data: " << endl;
    cout << "               (e.g., PREC) " << endl;

    cin >> iprecfname;

    cout << "Please enter the file extension (include the '.'): ";
    
    cin >> iprecend;
#endif

    rflog1 << "Please enter the first part of the file name";
    rflog1 << " containing the monthly precipitation data: " << endl;
    rflog1 << "               (e.g., PREC) " << endl;
    rflog1 << iprecfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << iprecend << endl << endl;
  }
  else
  {
#ifdef PMODE
	*fgo >> iprecfname;
#else
    cout << "Please enter the file name containing the long-term";
    cout << " mean monthly precipitation data: " << endl;
    cout << "               (e.g., PREC) " << endl;

    cin >> iprecfname;
#endif
    rflog1 << "Please enter the file name containing the long-term";
    rflog1 << " mean monthly precipitation data: " << endl;
    rflog1 << "               (e.g., PREC) " << endl;
    rflog1 << iprecfname << endl << endl;
  }

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void TEMclm60::initSolarRad( ofstream& rflog1 )
{

  if( 1 == sradflag )
  {
    if( 1 == tcldsflag )
    {
      if( 1 == cldflag )
      {
#ifdef PMODE
    	*fgo >> icldsfname;
	    *fgo >> icldsend;
#else
        cout << "Please enter the first part of the file name";
        cout << " containing the monthly cloudiness data: " << endl;
        cout << "               (e.g., CLDS) " << endl;
        rflog1 << "Please enter the first part of the file name";
        rflog1 << " containing the monthly cloudiness data: " << endl;
        rflog1 << "               (e.g., CLDS) " << endl;

        cin >> icldsfname;

        rflog1 << icldsfname << endl << endl;

        cout << "Please enter the file extension (include the '.'): ";
        
        cin >> icldsend;
#endif
        rflog1 << "Please enter the file extension (include the '.'): ";
        rflog1 << icldsend << endl << endl;
      }
      else
      {

#ifdef PMODE
	    *fgo >> inirrfname;
	    *fgo >> inirrend;
#else
        cout << "Please enter the first part of the file name";
        cout << " containing the monthly surface solar";
        cout << " radiation data: " << endl;
        cout << "               (e.g., NIRR) " << endl;

        cin >> inirrfname;

        cout << "Please enter the file extension (include the '.'): ";

        cin >> inirrend;
#endif
        rflog1 << "Please enter the first part of the file name";
        rflog1 << " containing the monthly surface solar";
        rflog1 << " radiation data: " << endl;
        rflog1 << "               (e.g., NIRR) " << endl;
        rflog1 << inirrfname << endl << endl;

        rflog1 << "Please enter the file extension (include the '.'): ";
        rflog1 << inirrend << endl << endl;
      }
    }
    else
    {
      if( 1 == cldflag )
      {
#ifdef PMODE
	    *fgo >> icldsfname;

#else
        cout << "Please enter the file name containing the";
        cout << " long-term mean monthly cloudiness data: ";
        cout << endl;
        cout << "               (e.g., CLDS) " << endl;

        cin >> icldsfname;
#endif
        rflog1 << "Please enter the file name containing the";
        rflog1 << " long-term mean monthly cloudiness data: ";
        rflog1 << endl;
        rflog1 << "               (e.g., CLDS) " << endl;
        rflog1 << icldsfname << endl << endl;
      }
      else
      {

#ifdef PMODE
	    *fgo >> inirrfname;

#else
        cout << "Please enter the file name containing the";
        cout << " long-term mean monthly surface solar";
        cout << " radiation data: " << endl;
        cout << "               (e.g., NIRR) " << endl;

        cin >> inirrfname;
#endif
        rflog1 << inirrfname << endl << endl;
        rflog1 << "Please enter the file name containing the";
        rflog1 << " long-term monthly surface solar";
        rflog1 << " radiation data: " << endl;
        rflog1 << "               (e.g., NIRR) " << endl;

      }
    }
  }
  else
  {
    if( 1 == tcldsflag )
    {
#ifdef PMODE
	  *fgo >> igirrfname;
	  *fgo >> igirrend;

	  *fgo >> inirrfname;
	  *fgo >> inirrend;
#else

      cout << "Please enter the first part of the file name";
      cout << " containing the monthly solar radiation at the";
      cout << " top of the atmosphere data: " << endl;
      cout << "               (e.g., GIRR) " << endl;
      cin >> igirrfname;

      cout << "Please enter the file extension (include the '.'): ";
      cin >> igirrend;

      cout << "Please enter the first part of the file name";
      cout << " containing the monthly surface solar radiation";
      cout << " data: " << endl;
      cout << "               (e.g., NIRR) " << endl;
      cin >> inirrfname;

      cout << "Please enter the file extension (include the '.'): ";
      cin >> inirrend;
#endif
      rflog1 << "Please enter the first part of the file name";
      rflog1 << " containing the monthly solar radiation at the";
      rflog1 << " top of the atmosphere data: " << endl;
      rflog1 << "               (e.g., GIRR) " << endl;

      rflog1 << igirrfname << endl << endl;

      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << igirrend << endl << endl;

      rflog1 << "Please enter the first part of the file name";
      rflog1 << " containing the monthly surface solar radiation";
      rflog1 << " data: " << endl;
      rflog1 << "               (e.g., NIRR) " << endl;
      rflog1 << inirrfname << endl << endl;

      rflog1 << "Please enter the file extension (include the '.'): ";
      rflog1 << inirrend << endl << endl;

      if( 1 == parflag )
      {
#ifdef PMODE
	    *fgo >> iparfname;
	    *fgo >> iparend;
#else
        cout << "Please enter the first part of the file name";
        cout << " containing the monthly photosynthetically";
        cout << " active radiation data: " << endl;
        cout << "               (e.g., PAR) " << endl;

        cin >> iparfname;

        cout << "Please enter the file extension (include the '.'): ";

        cin >> iparend;
#endif
        rflog1 << "Please enter the first part of the file name";
        rflog1 << " containing the monthly photosynthetically";
        rflog1 << " active radiation data: " << endl;
        rflog1 << "               (e.g., PAR) " << endl;
        rflog1 << iparfname << endl << endl;

        rflog1 << "Please enter the file extension (include the '.'): ";
        rflog1 << iparend << endl << endl;
      }
    }
    else
    {
#ifdef PMODE
	  *fgo >> igirrfname;
	  *fgo >> inirrfname;
#else
      cout << "Please enter the file name containing the";
      cout << " long-term mean monthly solar radiation at the";
      cout << " top of the atmosphere data: " << endl;
      cout << "               (e.g., GIRR) " << endl;

      cin >> igirrfname;

      cout << "Please enter the file name containing the";
      cout << " long-term mean monthly surface solar";
      cout << " radiation data: " << endl;
      cout << "               (e.g., NIRR) " << endl;

      cin >> inirrfname;
#endif
      rflog1 << "Please enter the file name containing the";
      rflog1 << " long-term mean monthly solar radiation at the";
      rflog1 << " top of the atmosphere data: " << endl;
      rflog1 << "               (e.g., GIRR) " << endl;

      rflog1 << igirrfname << endl << endl;

      rflog1 << "Please enter the file name containing the";
      rflog1 << " long-term mean monthly surface solar";
      rflog1 << " radiation data: " << endl;
      rflog1 << "               (e.g., NIRR) " << endl;

      rflog1 << inirrfname << endl << endl;

      if( 1 == parflag )
      {
#ifdef PMODE
	    *fgo >> iparfname;

#else
        cout << "Please enter the file name containing the";
        cout << " long-term mean monthly photosynthetically";
        cout << " active radiation data: " << endl;
        cout << "               (e.g., PAR) " << endl;

        cin >> iparfname;
#endif
        rflog1 << "Please enter the file name containing the";
        rflog1 << " long-term mean monthly photosynthetically";
        rflog1 << " active radiation data: " << endl;
        rflog1 << "               (e.g., PAR) " << endl;

        rflog1 << iparfname << endl << endl;
      }
    }
  }

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void TEMclm60::initTair( ofstream& rflog1 )
{

  if( 1 == ttairflag )
  {
#ifdef PMODE
	*fgo >> itairfname;
	*fgo >> itairend;
#else
    cout << "Please enter the first part of the file name";
    cout << " containing the monthly air temperature data: ";
    cout << endl;
    cout << "               (e.g., TAIR) " << endl;

    cin >> itairfname;

    cout << "Please enter the file extension (include the '.'): ";

    cin >> itairend;
#endif

    rflog1 << "Please enter the first part of the file name";
    rflog1 << " containing the monthly air temperature data: ";
    rflog1 << endl;
    rflog1 << "               (e.g., TAIR) " << endl;
    rflog1 << itairfname << endl << endl;
    rflog1 << "Please enter the file extension (include the '.'): ";
    rflog1 << itairend << endl << endl;
  }
  else
  {
#ifdef PMODE
	*fgo >> itairfname;

#else
    cout << "Please enter the file name containing the";
    cout << " long-term mean monthly air temperature data: ";
    cout << endl;
    cout << "               (e.g., TAIR) " << endl;

    cin >> itairfname;
#endif
    rflog1 << "Please enter the file name containing the";
    rflog1 << " long-term mean monthly air temperature data: ";
    rflog1 << endl;
    rflog1 << "               (e.g., TAIR) " << endl;
    rflog1 << itairfname << endl << endl;
  }

};

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

double TEMclm60::mkclds( const double& girr,
                         const double& nirr )
{

  double clouds;

  if ( nirr >= (0.76 * girr) ) { return clouds = ZERO; }
  else
  {
    clouds = 1.0 - (((nirr/girr) - 0.251)/0.509);
    clouds *= 100.0;
  }
  if ( clouds > 100.0 ) { clouds = 100.0; }

  return clouds;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double TEMclm60::mkd40( const double& lon,
                        const double& lat,
		        const string& contnent,
                        const double& o3,
		        const int& pdm )
{
  int k;
  int itmp;
  double d40;
  double tozone;
  double tmin;
  double tmax;
  double y;
  double tdum1;
  double tdum2;
  double slope1;
  double slope2;

  if ( lon > -90.0
       && lon < -55.0
       && lat > 10.0
       && lat < 90.0
       && "NAMERICA" == contnent )
  {
    //  Eastern North America
    tmin = 6.;
    tmax = 15.;
    y = 12.8;
  }
  else if ( lon > -124.0
            && lon < -112.0
            && lat > 47.0
            && lat < 55.0
            && "NAMERICA" == contnent )
  {
    //  Northwest North America
    tmin = 5.;
    tmax = 15.;
    y = 11.0;
  }
  else if ( lon > -180.0
            && lon < -123.0
            && lat > 53.0
            && lat < 90.0
            && "NAMERICA" == contnent )
  {
    //  Alaska and North
    tmin = 6.;
    tmax = 17.;
    y = 4.3;
  }
  else if ( lon > -180.0
            && lon < -90.0
            && lat > 10.0
            && lat < 90.0
            && "NAMERICA" == contnent )
  {
    //  West North America
    tmin = 6.;
    tmax = 16.;
    y = 7.9;
  }
  else if ( lon > -15.0
            && lon < 180.0
            && lat > 10.0
            && lat < 90.0
            && ("EUROPE" == contnent
            || "ASIA" == contnent) )
  {
    //  Eur-Asia
    tmin = 5.;
    tmax = 15.;
    y = -1.04 * lat + 75.72;
    if ( y < ZERO ) { y = ZERO; }
  }
  else
  {
    itmp=-999;
  }

  //  Next loop

  if ( 0 == pdm
       || 1 == pdm 
       || 2 == pdm 
       || 3 == pdm 
       || 10 == pdm 
       || 11 == pdm 
       || -999 == itmp )
  {
    d40 = ZERO;
  }
  else
  {
    tdum1 = tmin+(tmax-tmin)/2.0;
    tdum2 = ((24.-tmax)+tmin)/2.0;
    slope1 = (2*y)/(tmax-tmin);
    slope2 = (-1*2*y)/((24-tmax)+tmin);

    for ( k = 0; k < 24; ++k )
    {
      if ( k <= (tmin-1) )
      {
        tozone = ((slope2)*k+(-1*slope2*(tmin-tdum2)))+ o3;
      }

      if ( k > (tmin-1) && k <= (tmax-1) )
      {
        tozone = ((slope1)*k+(-1*slope1*tdum1))+ o3;
      }

      if ( k > (tmax-1) )
      {
        tozone = ((slope2)*k+(-1*slope2*(tmax+tdum2)))+ o3;
      }

      if ( tozone >= 40.0 && k >= 6 && k <= 18 )
      {
        d40 += (tozone-40.0);
      }
    }

    d40 *= 31.0;
  }

  return d40;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMclm60::setCldsFlags( ofstream& rflog1,
                             const int& requil )
{

#ifdef PMODE
	*fgo >> sradflag;

#else
  cout << endl;
  cout << "Do you want to run the TEMclm model for solar";
  cout << " radiation variables?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes" << endl;

  cin >> sradflag;
#endif
  rflog1 << endl;
  rflog1 << "Do you want to run the TEMclm model for solar";
  rflog1 << " radiation variables?" << endl;
  rflog1 << "  Enter 0 for No" << endl;
  rflog1 << "  Enter 1 for Yes" << endl;
  rflog1 << " sradflag = " << sradflag << endl << endl;

#ifdef PMODE
	*fgo >> cldflag;

#else
  cout << "Do you have spatially explicit surface solar";
  cout << " radiation data or cloudiness data?:" << endl;
  cout << "Enter 0 for surface solar radiation data (W/ sq. m):";
  cout << endl;
  cout << "Enter 1 for cloudiness data (percent cloudiness): ";
  cout << endl;

  cin >> cldflag;
#endif
  rflog1 << "Do you have spatially explicit surface solar";
  rflog1 << " radiation data or cloudiness data?:" << endl;
  rflog1 << "Enter 0 for surface solar radiation data (W/ sq. m):";
  rflog1 << endl;
  rflog1 << "Enter 1 for cloudiness data (percent cloudiness):";
  rflog1 << endl;
  rflog1 << "cldflag = " << cldflag << endl;

  parflag = 0;
  
  if( 0 == sradflag )
  {
#ifdef PMODE
	*fgo >> parflag;

#else
    cout << "Do you have spatially explicit PAR data? ";
    
    cin >> parflag;
#endif
    rflog1 << "Do you have spatially explicit PAR data? ";
    rflog1 << parflag << endl << endl;
  }

  tcldsflag = 0;
  
  if( 0 == requil )
  {
    if( 1 == cldflag )
    {
#ifndef PMODE
      cout << "Do you have transient cloudiness data?:" << endl;
#endif
      rflog1 << endl;
      rflog1 << "Do you have transient cloudiness data?:";
      rflog1 << endl;
    }
    else
    {
#ifndef PMODE
      cout << "Do you have transient solar radiation data?:";
      cout << endl;
#endif
      rflog1 << endl;
      rflog1 << "Do you have transient solar radiation data?:";
      rflog1 << endl;
    }

#ifdef PMODE
    *fgo >>tcldsflag;
#else
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes: ";

    cin >> tcldsflag;
#endif
    rflog1 << "Enter 0 for no:" << endl;
    rflog1 << "Enter 1 for yes: " << endl;
    rflog1 << "tcldsflag = " << tcldsflag << endl << endl;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMclm60::setCO2Flags( ofstream& rflog1, const int& requil )
{

  tco2flag = 0;
  
  if( 0 == requil )
  {
#ifdef PMODE
    *fgo >>tco2flag;
#else
    cout << "Do you have transient CO2 data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for annual transient CO2 data: ";
    cout << "Enter 2 for monthly transient CO2 data: ";

    cin >> tco2flag;
#endif
    rflog1 << "Do you have transient CO2 data?:" << endl;
    rflog1 << "Enter 0 for No:" << endl;
    rflog1 << "Enter 1 for annual transient CO2 data: ";
    rflog1 << "Enter 2 for monthly transient CO2 data: ";
    rflog1 << "tco2flag = " << tco2flag << endl << endl;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMclm60::setNdepFlags( ofstream& rflog1, const int& requil )
{

  tndepflag = 0;
  
  if( 0 == requil )
  {
#ifdef PMODE
    *fgo >>tndepflag;
#else
    cout << "Do you have transient N deposition data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";

    cin >> tndepflag;
#endif
    rflog1 << endl;
    rflog1 << "Do you have transient N deposition data?:";
    rflog1 << endl;
    rflog1 << "Enter 0 for No:" << endl;
    rflog1 << "Enter 1 for Yes: " << endl;
    rflog1 << "tndepflag = " << tndepflag << endl << endl;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMclm60::setO3Flags( ofstream& rflog1, const int& requil )
{

  to3flag = 0;
  
  if( 0 == requil )
  {
#ifdef PMODE
    *fgo >>to3flag;
#else
    cout << "Do you have transient O3 data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    
    cin >> to3flag;
#endif
    rflog1 << "Do you have transient O3 data?:" << endl;
    rflog1 << "Enter 0 for No:" << endl;
    rflog1 << "Enter 1 for Yes: " << endl;
    rflog1 << "to3flag = " << to3flag << endl << endl;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMclm60::setPrecFlags( ofstream& rflog1, 
                             const int& requil )
{

  tprecflag = 0;
  
  if( 0 == requil )
  {
#ifdef PMODE
    *fgo >>tprecflag;
#else
    cout << "Do you have transient precipitation data?:";
    cout << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    
    cin >> tprecflag;
#endif
    rflog1 << endl;
    rflog1 << "Do you have transient precipitation data?:";
    rflog1 << endl;
    rflog1 << "Enter 0 for No:" << endl;
    rflog1 << "Enter 1 for Yes: " << endl;
    rflog1 << "tprecflag = " << tprecflag << endl << endl;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMclm60::setTairFlags( ofstream& rflog1,
                             const int& requil )
{

  ttairflag = 0;
  
  if( 0 == requil )
  {
#ifdef PMODE
    *fgo >>ttairflag;
#else
    cout << "Do you have transient air temperature data?:";
    cout << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    
    cin >> ttairflag;
#endif
    rflog1 << endl;
    rflog1 << "Do you have transient air temperature data?:";
    rflog1 << endl;
    rflog1 << "Enter 0 for No:" << endl;
    rflog1 << "Enter 1 for Yes: " << endl;
    rflog1 << "ttairflag = " << ttairflag << endl << endl;
  }

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

double TEMclm60::xgirr( const double& plat,
                        const int& pdm,
                        double& psumday )
{

  const double pi = 3.141592654;                // Greek "pi"
  const double sp = 1368.0 * 3600.0 / 41860.0;  // solar constant

  double lambda;
  double sumd;
  double sig;
  double eta;
  double sinbeta;
  double sb;
  double sotd;
  int day;
  int hour;
  double gross;

  ndays[0] = 31;
  ndays[1] = 28;
  ndays[2] = 31;
  ndays[3] = 30;
  ndays[4] = 31;
  ndays[5] = 30;
  ndays[6] = 31;
  ndays[7] = 31;
  ndays[8] = 30;
  ndays[9] = 31;
  ndays[10] = 30;
  ndays[11] = 31;

  lambda = plat * pi / 180.0;

  gross = ZERO;
  for( day = 0; day < ndays[pdm]; ++day )
  {
    ++psumday;
    sumd = 0;
    sig = -23.4856*cos(2 * pi * (psumday + 10.0)/365.25);
    sig *= pi / 180.0;

    for( hour = 0; hour < 24; ++hour )
    {
      eta = (double) ((hour+1) - 12) * pi / 12.0;
      sinbeta = sin( lambda ) * sin( sig )
                + cos( lambda ) * cos( sig ) * cos( eta );

      sotd = 1 - (0.016729 
             * cos( 0.9856 * (psumday - 4.0) * pi / 180.0) );

      sb = sp * sinbeta / pow( sotd, 2.0 );
      if( sb >= ZERO ) { sumd += sb; }
    }

    gross += sumd;
  }

  gross /= (double) ndays[pdm];

  return gross;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double TEMclm60::xnirr( const double& clds, const double& girr )
{

  double nirr;

  if( clds >= ZERO )
  {
    nirr = girr * (0.251 + (0.509*(1.0 - clds/100.0)));
  }
  else { nirr = MISSING; }

  return nirr;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double TEMclm60::xpar( const double& clds, const double& nirr )
{

  double par;

  if( clds >= ZERO )
  {
      par = nirr * ((0.2 * clds / 100.0) + 0.45);
  }
  else { par = MISSING; }

  return par;

};


