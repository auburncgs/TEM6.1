/* **************************************************************
*****************************************************************
ATMS602.H - object describes physical characteristics of the
	       atmosphere

Modifications:

20030513 - DWK created by modifying atms43.h to incorporate
           open-N cycle dynamics into TEM
20030513 - DWK changed Atmosphere43 to Atmosphere51
20030513 - DWK added int tndepflag
20030513 - DWK added double ndep[CYCLE], double yrtotndep,
           double yrnh4dep, double yrno3dep;
20030513 - DWK changed include from atms43.cpp to atms51.cpp at
           bottom of file
20030717 - DWK changed public double no2[CYCLE] to nox[CYCLE]
20030718 - DWK added include "inorgn51.h"
20030718 - DWK changed public double ndep[CYCLE] to
           InorgN ndep[CYCLE]
20040716 - DWK changed include from inorgn51.hpp to inorgn60.hpp
20040716 - DWK changed class Atmosphere51 to class Atmosphere60
20040716 - DWK changed public InorgN ndep[CYCLE] to
           InorgN60 ndep[CYCLE]
20040716 - DWK changed public InorgN yrndep to InorgN60 yrndep
20040716 - DWK changed size of girr[], nirr[], par[], tair[] and 
           prec[] arrays from CYCLE to CYCLE+1
20040716 - DWK changed include from atms51.cpp to atms60.cpp at
           bottom of file
20051117 - DWK deleted include "atms.cpp" at bottom of file
20051117 - DWK added include "temconsts602.hpp"
20051122 - DWK changed public double d40[CYCLE] to 
           double aot40[CYCLE]
20051124 - DWK changed public double tair[CYCLE+1] to
           double tair
20051124 - DWK changed public double co2[CYCLE+1] to 
           double co2
20051124 - DWK changed public InorgN60 ndep[CYCLE+1]
           to InorgN60 ndep
20051124 - DWK changed public double aot40[CYCLE+1] to
           double aot40
20051124 - DWK changed public double eet[CYCLE] to
           double eet
20051124 - DWK changed public double ineet[CYCLE] to
           double ineet
20051124 - DWK changed public double pet[CYCLE] to
           double pet
20051124 - DWK changed public double  prec[CYCLE+1] to
           double prec
20051124 - DWK changed public double  snowfall[CYCLE+1] to
           double snowfall
20051124 - DWK changed public double  rain[CYCLE+1] to 
           double rain
20051124 - DWK changed public double clds[CYCLE+1] to
           double clds
20051124 - DWK changed public double girr[CYCLE+1] to
           double girr
20051124 - DWK changed public double nirr[CYCLE+1] to
           double nirr
20051124 - DWK changed public double par[CYCLE+1] to
           double par
20051124 - DWK deleted public flags for transient data: 
           int tco2flag, int to3flag, int tndepflag, 
           int ttairflag, int tprecflag and int tcldsflag
20051124 - DWK added public double prevco2
20051130 - DWK deleted public double eet, double eetmx, 
           double ineet, double prveetmx, double yreet and 
           double yrineet
20051130 - DWK deleted public function xeet()
20051201 - DWK added public function resetMonthlyFluxes()
20051201 - DWK added public function resetYrFluxes()
                                                                                                                                            
*****************************************************************
************************************************************** */

// Class uses global constants CYCLE and MAXRTIME

#ifndef ATMS602_H
#define ATMS602_H

#include "temconsts602.hpp"
#include "inorgn60.hpp"      // Atmospere60 uses InorgN60 class


class Atmosphere60
{

  public:

     Atmosphere60();

/* *************************************************************
		 Public Functions
************************************************************* */

     /* Determine Potential Evapotranspiration based on 
        algorithms from Jensen M. E. and H. R. Haise (1963) 
        Estimating evapotranspiration from solar radiation.  
        Journal of the Irrigation and Drainage Division 4: 
        14-41.  */

     void petjh( const double& nirr,
                 const double& tair,
                 const int& pdm );

     // Determine in precipitation occurs as rain or snow

     void precsplt( const double& prec,
                    const double& tair,
                    double& rain,
                    double& snowfall );

     void resetMonthlyFluxes( void );
     
     void resetYrFluxes( void );


     // "Get" and "Set" private variables

     // aot40 **************************************************
     
     inline double getAOT40( void ) { return aot40; }

     inline void setAOT40( const double& paot40 ) 
     { 
       aot40 = paot40; 
     }

     // avetair ************************************************
     
     inline double getAVETAIR( void ) { return avetair; }

     inline void setAVETAIR( const double& pavetair ) 
     { 
       avetair = pavetair; 
     }

     // clds ***************************************************
     
     inline double getCLDS( void ) { return clds; }

     inline void setCLDS( const double& pclds ) 
     { 
       clds = pclds; 
     }


     // co2 ****************************************************
     
     inline double getCO2( void ) { return co2; }

     inline void setCO2( const double& pco2 ) { co2 = pco2; }

     
     // co2level ***********************************************
     
     inline double getCO2LEVEL( void ) { return co2level; }

     inline void setCO2LEVEL( const double& pco2level ) 
     { 
       co2level = pco2level; 
     }

 
     // girr ***************************************************
     
     inline double getGIRR( void ) { return girr; }

     inline void setGIRR( const double& pgirr ) 
     { 
       girr = pgirr; 
     }


     // initco2 ************************************************
     
     inline double getINITCO2( void ) { return initco2; }

     inline void setINITCO2( const double& pinitco2 ) 
     { 
       initco2 = pinitco2; 
     }


     // mxtair *************************************************
     
     inline double getMXTAIR( void ) { return mxtair; }

     inline void setMXTAIR( const double& pmxtair ) 
     { 
       mxtair = pmxtair; 
     }


     // ndays[] ************************************************
     
     inline double getNDAYS( const int& pdm ) 
     { 
       return ndays[pdm]; 
     }


     // ndep ***************************************************
     
     inline InorgN60 getNDEP( void ) { return ndep; }


     // ndep.nh4 ***********************************************
     
     inline double getNH4DEP( void ) { return ndep.nh4; }

     inline void setNH4DEP( const double& pnh4dep ) 
     { 
       ndep.nh4 = pnh4dep; 
     }

     // ndep.no3 ***********************************************
     
     inline double getNO3DEP( void ) { return ndep.no3; }

     inline void setNO3DEP( const double& pno3dep ) 
     { 
       ndep.no3 = pno3dep; 
     }

     // ndep.total *********************************************
     inline double getTOTNDEP( void ) { return ndep.total; }

     inline void setTOTNDEP( const double& ptotndep ) 
     { 
       ndep.total = ptotndep; 
     }

     // nirr ***************************************************
     
     inline double getNIRR( void ) { return nirr; }

     inline void setNIRR( const double& pnirr ) 
     { 
       nirr = pnirr; 
     }


     // par ****************************************************
     
     inline double getPAR( void ) { return par; }

     inline void setPAR( const double& ppar ) { par = ppar; }


     // pet ****************************************************
     
     inline double getPET( void ) { return pet; }

     inline void setPET( const double& ppet ) { pet = ppet; }


     // petmx **************************************************
     
     inline double getPETMX( void ) { return petmx; }

     inline void setPETMX( const double& ppetmx ) 
     { 
       petmx = ppetmx; 
     }


     // prec ***************************************************
     
     inline double getPREC( void ) { return prec; }

     inline void setPREC( const double& pprec ) 
     { 
       prec = pprec; 
     }


     // prev2tair **********************************************
     
     inline double getPREV2TAIR( void ) { return prev2tair; }

     inline void setPREV2TAIR( const double& pprev2tair ) 
     { 
       prev2tair = pprev2tair; 
     }


     // prevco2 ************************************************
     
     inline double getPREVCO2( void ) { return prevco2; }

     inline void setPREVCO2( const double& pprevco2 ) 
     { 
       prevco2 = pprevco2; 
     }


     // prevtair ***********************************************
     
     inline double getPREVTAIR( void ) { return prevtair; }

     inline void setPREVTAIR( const double& pprevtair ) 
     { 
       prevtair = pprevtair; 
     }

     // prevrain ***********************************************

     inline double getPREVRAIN( void ) { return prevrain; }

     inline void setPREVRAIN( const double& pprevrain )
     {
       prevrain = pprevrain;
     }
     //prevmdc
     inline double getPREVMDC( void ) { return prevmdc; }

     inline void setPREVMDC( const double& pprevmdc )
     {
       prevmdc = pprevmdc;
     }
     // prvpetmx ***********************************************
     
     inline double getPRVPETMX( void ) { return prvpetmx; }

     inline void setPRVPETMX( const double& pprvpetmx ) 
     { 
       prvpetmx = pprvpetmx; 
     }


     // rain ***************************************************
     
     inline double getRAIN( void ) { return rain; }

     inline void setRAIN( const double& prain ) 
     { 
       rain = prain; 
     }


     // snowfall ***********************************************
     
     inline double getSNOWFALL( void ) { return snowfall; }

     inline void setSNOWFALL( const double& psnowfall ) 
     { 
       snowfall = psnowfall; 
     }


     // tair ***************************************************
     
     inline double getTAIR( void ) { return tair; }
        
     inline void setTAIR( const double& ptair ) 
     { 
       tair = ptair; 
     }

  



/* **************************************************************
		 Public Variables
************************************************************** */

     // Number of days per month     
     double ndays[CYCLE];   

    // Annual total N deposition (mg N / (sq. meter * year))
     InorgN60 yrndep;

     // Annual NH4 deposition (mg N / (sq. meter * year))
     double yrnh4dep;

     // Annual NO3 deposition (mg N / (sq. meter * year))
     double yrno3dep;

     // Annual potential evapotranspiration (mm / year)
     double yrpet;

     // Annual total precipitation (mm / year)
     double yrprec;

     // Annual sum of rainfall (mm / year)
     double yrrain;

     // Annual snow ( mm / year)    
     double yrsnowfall;         


   private:
   
/* **************************************************************
		 Private Variables
************************************************************** */

     // Atmospheric ozone AOT 40 index (i.e. ozone 
     //   concentrations above 40 ppb-hr)
     double aot40;

     // Mean annual air temperature (degrees C)
     double avetair;

     // Cloudiness (%)    
     double clds;   
     
     // Atmospheric carbon dioxide concentration (ppmv)
     double co2; 

     // Constant atmospheric CO2 concentration (ppmv) used to 
     //   calibrate TEM or equilibrate TEM at beginning of 
     //   extrapolation     
     double co2level;

     // Gross Irradiance (cal/(sq. cm * day))
     double girr;


     // initial CO2 concentration (ppmv)     
     double initco2;     

     // Maximum monthly air temperature (degrees C)
     double mxtair;


     // Monthly N deposition (mg N / (sq. meter * month)
     InorgN60 ndep;

     // Net Irradiance   (cal/(sq. cm * day))
     double nirr;

     // Photosynthetically Active Radiation  (cal/(sq.cm * day))
     double par;

     // Monthly potential evapotranspiration (mm / month)
     double pet;

     // Maximum PET of current year (mm / month)
     double petmx;

     // Total Precipitation (mm / month)
     double  prec;
     
     // Previous Month's Atmospheric CO2 Concentration (ppmv)
     double prevco2;

     // Previous Month's Air Temperature (degrees C)
     double prevtair;
     // Previous Month's rain (mm)
     double prevrain;
     //previous monthly drought code (MDC)
     double prevmdc;

     // Previous 2 Month's Air Temperature (degrees C)
     double prev2tair;

     // Maximum PET of previous year (mm / month)
     double prvpetmx;

     // Rainfall (mm / month)     
     double  rain;       

     // Snowfall (mm / month)     
     double  snowfall;   

     // Surface Air Temperature (degrees C)
     double tair;



};

#endif
