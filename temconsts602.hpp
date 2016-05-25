#ifndef ACCEPT_CONST
#define ACCEPT_CONST
  const int ACCEPT = 0;
#endif

// TEMVEG index of crops
#ifndef CROPVEG_CONST
#define CROPVEG_CONST
  const int CROPVEG = 50;
#endif

// Number of months in an annual cycle
#ifndef CYCLE_CONST    
#define CYCLE_CONST
  const int CYCLE = 12;
#endif

// Maximum number of cohorts in a grid cell
#ifndef MAXCHRTS_CONST    
#define MAXCHRTS_CONST
  const int MAXCHRTS = 3000;
//  const int MAXCHRTS = 11;
#endif

// Maximum number of community types for TEM parameters
#ifndef MAXCMNT_CONST
#define MAXCMNT_CONST
  const int MAXCMNT = 15;
#endif

// Maximum number of ecosystem (i.e. carbon and nitrogen)
//   pool variables in TEM
#ifndef MAXESTAT_CONST
#define MAXESTAT_CONST
  //const int MAXESTAT = 9;
  //const int MAXESTAT = 11; //modified by cgs to add two pools: I_deadwoodc, I_deadwoodn
  //const int MAXESTAT = 15; //modified by cgs to add three pools: I_FOLC, I_STEMC, I_CROOTC, I_FROOTC
  //const int MAXESTAT = 20; //modified by cgs to add five pools: I_CWD, I_AGR, I_AGL, I_BGR,I_BGL
  //const int MAXESTAT = 21; //modified by cgs to add monthly drought code (I_MDC)
  //const int MAXESTAT = 22; //modified by cgs to add forest age output (I_AGE)
  const int MAXESTAT = 23; //modified by cgs to add pure soil organic carbon output (I_SOC)
#endif

#ifndef MAXFRI_CONST
#define MAXFRI_CONST
  const double MAXFRI = 2000;
#endif

#ifndef MAXGRID_CONST    
#define MAXGRID_CONST
  const int MAXGRID = 1;
#endif

// Maximum number of snow and soil nodes for STM
#ifndef MAXNODES_CONST    
#define MAXNODES_CONST
  const int MAXNODES = 210;
#endif

// Maximum number of years in a simulation
#ifndef MAXRTIME_CONST    
#define MAXRTIME_CONST
  const int MAXRTIME = 1200;
//  const int MAXRTIME = 20100;
#endif

// Maximum number of nodes in snowpack for STM
#ifndef MAXSNODES_CONST
#define MAXSNODES_CONST
  const int MAXSNODES = 9;
#endif

// Maximum number of water pool variables
#ifndef MAXWSTAT_CONST
#define MAXWSTAT_CONST
  const int MAXWSTAT = 3;
#endif

// Minimum active layer depth is 0.0001 m (used to avoid 
//   dividing by zero when calculating volumeteric soil 
//   moisture)
#ifndef MINACTLAYERZ_CONST
#define MINACTLAYERZ_CONST
  const double MINACTLAYERZ = 0.0001;
#endif

// Default value for missing data
#ifndef MISSING_CONST    
#define MISSING_CONST
  const double MISSING = -999999.9;
#endif

// Total number of TEMCLM output variables
#ifndef NUMATMS_CONST
#define NUMATMS_CONST
const int NUMATMS = 13;
#endif

// Maximum number of vegetation subtypes in a mosaic
#ifndef NUMMSAC_CONST
#define NUMMSAC_CONST
  const int NUMMSAC = 5;
#endif

// Total number of STM output variables
#ifndef NUMSTM_CONST
#define NUMSTM_CONST
const int NUMSTM = 13;
#endif

// Maximum number of vegetation types
#ifndef NUMVEG_CONST
#define NUMVEG_CONST
  const int NUMVEG = 56;
#endif

#ifndef REJECT_CONST
#define REJECT_CONST
  const int REJECT = 1;
#endif

#ifndef ZERO_CONST    
#define ZERO_CONST
  const double ZERO = 0.00;
#endif

/* *************************************************************
     Constants that are combinations of other constants
************************************************************* */

// Maximum number of pool variables in simulation
#ifndef MAXSTATE_CONST
#define MAXSTATE_CONST
  const int MAXSTATE = MAXWSTAT + MAXESTAT;
  //maxwstat: 3 + maxestat: 15
  //maxwstat: 3 + maxestat: 23
#endif

// Maximum number of state (pools plus fluxes) variables
//  used to simulate hydrology in TEM 
#ifndef NUMWEQ_CONST
#define NUMWEQ_CONST
const int NUMWEQ = MAXWSTAT + 7; //10
#endif

// Maximum number of state (pools plus fluxes) variables
//  used to simulate carbon and nitrogen dynamics in TEM 
#ifndef NUMEEQ_CONST
#define NUMEEQ_CONST
//cgs2014 added a variable "DST300", change 58 to 59
  //const int NUMEEQ = MAXESTAT + 59; //83
  const int NUMEEQ = MAXESTAT + 63; //cgs2014 added four fluxes: ffdc, fflc, swfc, dwfc


#endif

// Total number of state (pools plus fluxes variables
//  (carbon, nitrogen and water) used in TEM 
#ifndef NUMEQ_CONST
#define NUMEQ_CONST
  const int NUMEQ = NUMWEQ + NUMEEQ; //93
  //const int NUMEQ = NUMWEQ + NUMEEQ+151; //93 //modified by cgs2014 to assign initial values in resetODEflux ()
#endif

// Total number of output variables from TEM 
#ifndef NUMTEM_CONST
#define NUMTEM_CONST
  const int NUMTEM = NUMEQ + 151; //93 + 151
  //const int NUMTEM = NUMEQ; //93 + 151
#endif
                      
