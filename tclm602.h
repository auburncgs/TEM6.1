/* *************************************************************
****************************************************************
TCLM602.H - describes climate module used by TEM

20030513 - DWK created by modifying tclm43.h to incorporate
           open-N cycle dynamics into TEM
20030513 - DWK changed include from atms43.h to atms51.h
20030513 - DWK changed Class TEMclm43 to Class TEMclm51
20030513 - DWK added setNdepFlags() and initNdep()
20030513 - DWK added char indepfname[] and char independ[]
20030513 - DWK changed include from tclm43.cpp to tclm51.cpp at
           bottom of file
20030718 - DWK changed long co2year[MAXRTIME] to
           int co2year[MAXRTIME]
20030718 - DWK added I_CO2, I_TNDEP, I_NH4DEP and I_NO3DEP to
           clmkey
20030718 - DWK changed global const int NUMATMS = 9 to
           const int NUMATMS = 12
20030718 - DWK changed public char predstr[NUMATMS+4][9] to
           char predstr[NUMATMS+3][9]
20031019 - DWK changed char igirrfname[MAXFNAME],
           char igirrend[25], char icldsfname[MAXFNAME],
           char icldsend, char inirrfname[MAXFNAME],
           char inirrend[25], char iparfname[MAXFNAME],
           char iparend[25], char itairfname[MAXFNAME],
           char itairend[25], char iprecfname[MAXFNAME],
           char iprecend[25], char ico2fname[MAXFNAME],
           char ico2end[25], char io3fname[MAXFNAME],
           char io3end[25], char indepfname[MAXFNAME],
           and char independ[25] to strings
20040714 - DWK changed class TEMclm51 to TEMclm60
20040714 - DWK changed public char predstr[NUMATMS+3][9] to
           string predstr[NUMATMS];
20040714 - DWK changed include from tclm51.cpp to tclm60.cpp at
           bottom of file
20040716 - DWK changed include from atms51.h to atms60.h
20040716 - DWK changed inheritance of Atmosphere51 to
           inheritance of Atmosphere60
20051117 - DWK added include temconstants.hpp
20051117 - DWK changed include from atms60.h to atms602.h
20051117 - DWK changed include from tco2dat425.h to tco2dat602.h
20051117 - DWK deleted tclm60.cpp from bottom of file
20051121 - DWK added I_RAIN and I_SNWFAL to enum clmkey
20051121 - DWK deleted I_O3 from enum clmkey
20051121 - DWK changed I_D40 to I_AOT40 in enum clmkey
20051121 - DWK changed public string predstr[NUMATMS] to
           vector<string> predstr
20051121 - DWK deleted public function loadyrCO2()
20051121 - DWK deleted public double tco2[MAXRTIME][CYCLE]
20051121 - DWK deleted public int co2year[MAXRTIME]
20051121 - DWK deleted include tco2dat602.h
20051124 - DWK added public flags for transient data: 
           int tco2flag, int to3flag, int tndepflag, 
           int ttairflag, int tprecflag and int tcldsflag
20051124 - DWK changed char contnent[9] to string contnent
           in function call to mkd40()
                                         
****************************************************************
************************************************************* */

#ifndef TCLM602_H
#define TCLM602_H



// TEMclm uses the global constants CYCLE, MISSING, NUMATMS 
//   and NUMATMS
#include "temconsts602.hpp"

// TEMclm60 inherits class Atmosphere60
#include "atms602.h"        


class TEMclm60 : public Atmosphere60
{

  public:

     TEMclm60();

     enum clmkey { I_GIRR,   I_NIRR,   I_PAR,    I_CLDS, I_TAIR,
                   I_PREC,   I_RAIN,   I_SNWFAL, I_CO2,  I_AOT40,  
                   I_TNDEP,  I_NH4DEP, I_NO3DEP };

/* *************************************************************
		 Public Functions
************************************************************* */

     // Get file names of input data sets

     void initCO2( ofstream& rflog1 );
     void initNdep( ofstream& rflog1 );
     void initO3 ( ofstream& rflog1 );
     void initPrec( ofstream& rflog1 );
     void initSolarRad( ofstream& rflog1 );
     void initTair( ofstream& rflog1 );


     /* Determine cloudiness based on solar radiation at the top
        of the atmosphere and solar radiation at the top of the
        vegetation canopy */

     double mkclds( const double& girr, const double& nirr );

     double mkd40( const double& lon,
                   const double& lat,
	  	   const string& contnent,
                   const double& o3,
	  	   const int& pdm );

     void setCldsFlags( ofstream& rflog1, const int& requil );
     void setCO2Flags( ofstream& rflog1, const int& requil );
     void setNdepFlags( ofstream& rflog1, const int& requil );
     void setO3Flags( ofstream& rflog1, const int& requil );
     void setPrecFlags( ofstream& rflog1, const int& requil );
     void setTairFlags( ofstream& rflog1, const int& requil );

     /* Determine solar radiation at the top of the atmosphere (i.e.,
        "gross irradiance" or GIRR) based on the solar constant, time
        of year and latitude as described by S. M. Turton (1986)
        Solar radiation under cloudless skies.  Weatherwise 39:
        223-224.  */

     double xgirr( const double& plat, 
                   const int& pdm, 
                   double& psumday );

     /* Determine solar radiation at the top of the vegetation canopy
        (i.e., "net irradiance" or NIRR) based on GIRR and cloudiness */

     double xnirr( const double& clds, const double& girr );

     /* Determine phototsynthetically active radiation (PAR) based on
        NIRR and cloudiness.  Algorithm determined by Jim Raich from
        a variety of studies [e.g., McCree (1966) Agricul. Meteorol.
        3: 353-366; Stigter and Musabilha (1982) Journal of Applied
        Ecology 19: 853-858] that indicate that cloud cover increases
        the proportion of PAR. */

     double xpar( const double& clds, const double& nirr );


/* *************************************************************
		 Public Variables
************************************************************* */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif

     // "Do data sets have cloudiness data or NIRR data?" flag
     int cldflag;     

     ifstream fco2;

     // Name of file extension and beginning of file name
     //  for cloudiness data

     string icldsend;
     string icldsfname;

     // Name of file extension and beginning of file name
     //  for atmospheric CO2 concentration data

     string ico2end;
     string ico2fname;

     // Name of file extension and beginning of file name
     //  for gross irradiance (i.e., solar radiation at the
     //  top of the atmosphere) data

     string igirrend;
     string igirrfname;

     // Name of file extension and beginning of file name
     //  for atmospheric nitrogen deposition data

     string independ;
     string indepfname;

     // Name of file extension and beginning of file name
     //  for net irradiance (i.e. solar radiation at the 
     //  top of the vegetation canopy) data

     string inirrend;
     string inirrfname;

     // Name of file extension and beginning of file name
     //  for ozone (i.e. AOT40) data

     string io3end;
     string io3fname;

     // Name of file extension and beginning of file name
     //  for photosynthetically active radiation data

     string iparend;
     string iparfname;

     // Name of file extension and beginning of file name
     //  for precipitation data

     string iprecend;
     string iprecfname;

     // Name of file extension and beginning of file name
     //  for air temperature data
     
     string itairend;
     string itairfname;

     int parflag;     // Read in PAR from spatially explicit data sets?

     int predflag;    // Write climate module output to files?

     // Names of climate output variables
     
     vector<string> predstr;

     int sradflag;    // Run solar radiation module?

     // Initial year of transient portion of simulation
     
     int modstartyr;

     // Flag for transient cloudiness data 
     //   ( =0 for long-term mean data
     //      1 for transient data)

     int tcldsflag;


     // Flag for transient atmospheric CO2 data 
     //   ( =0 for long-term mean data
     //      1 for transient data)

     int tco2flag;

     // Flag for transient atmospheric N deposition data 
     //   ( =0 for long-term mean data
     //      1 for transient data)

     int tndepflag;

     // Flag for transient ozone (AOT40) data 
     //   ( =0 for long-term mean data
     //      1 for transient data)

     int to3flag;

     // Flag for transient precipitation data 
     //   ( =0 for long-term mean data
     //      1 for transient data)

     int tprecflag;


     // Flag for transient air temperature data 
     //   ( =0 for long-term mean data
     //      1 for transient data)

     int ttairflag;

     // Year represented by the climate data 
     
     int year;

     // Number of days since beginning of year
     
     double yrsumday;


};

#endif
