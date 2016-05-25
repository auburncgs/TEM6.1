/* **************************************************************
*****************************************************************
TELM604.H - Determines the interface among CLM, LCLUC and TEM
              modules

Modifications:

20030712 - DWK created by modifying telm50b1.h to incorporate
           open-N cycle dynamics into TEM
20030712 - DWK changed include from tclm425.h to tclm51.h
20030712 - DWK changed include from ttem50b1.h to ttem51.h
20030712 - DWK changed Class TEMelmnt to Class TEMelmnt51
20030712 - DWK changed TEMclm clm to TEMclm51 clm
20030712 - DWK changed TTEM50 tem to TTEM51 tem
20030712 - DWK added double d40[CYCLE] to temgisqc() function call
20030712 - DWK added double ndep[CYCLE] to temgisqc() function call
20030712 - DWK changed ofstream ftemout[MAXPRED] to ofstream
           ftemout[NUMTEM] in setTEMequilState(), setTEMmiss(),
           temwritemiss(), and temwritepred()
20030712 - DWK changed TEMList50 gtem to TEMList51 gtem in
           function calls to setTEMequilState(), updateTEMmonth()
           and updateTEMyear()
20030712 - DWK changed include from telm50b1.cpp to telm51.cpp
           at bottom of file
20031021 - DWK changed include from tlcluc431a.h to tlcluc51.h
20031021 - DWK added public double subarea
20040203 - DWK changed include from ttem51.h to ttem511.h
20040203 - DWK changed include from ttem51.cpp to ttem511.cpp 
           at bottom of file
20040203 - DWK changed include from ttem511.h to ttem512.h
20040203 - DWK changed include from ttem511.cpp to ttem512.cpp 
           at bottom of file
20040220 - DWK changed include from ttem512.h to ttem513.h
20040220 - DWK changed include from ttem512.cpp to ttem513.cpp 
           at bottom of file
20040229 - DWK changed include from ttem513.h to ttem60.h
20040229 - DWK changed class TEMelmnt51 to class TEMelmnt60
20040229 - DWK changed TTEM51 tem to TTEM60 tem
20040229 - DWK changed include from telm513.cpp to telm60.cpp 
           at bottom of file
20040714 - DWK changed public char temstateAfname[MAXFNAME]
           to string temstateAfname
20040714 - DWK changed public char temstateBfname[MAXFNAME]
           to string temstateBfname
20040714 - DWK changed include from tlcluc51.h to tlcluc60.h
20040714 - DWK changed public TEMlcluc51 lcluc to TEMlcluc60 lcluc
20040714 - DWK changed char predmap[MAXPRED][9] to 
           string predmap[MAXPRED] in atmswritemiss(),
           atmswritepred(), setTEMequilState(), setTEMmiss(),
           temwritemiss(), temwritepred(), updateTEMmonth() and
           updateTEMyear()
20040714 - DWK changed include from tclm51.h to tclm60.h
20040714 - DWK changed TEMclm51 clm to TEMclm60 clm
20040716 - DWK changed include from ttemdat425.h to ttemdat433.h
20040716 - DWK changed public int carea to long carea
20040716 - DWK changed public double subarea to long subarea
20040716 - DWK changed public char contnent[9] to string contnent
20040716 - DWK changed public long atmstotyr[MAXRTIME] to 
           int atmstotyr[MAXRTIME]
20040716 - DWK changed public long ttotyr[MAXRTIME] to 
           int ttotyr[MAXRTIME]
20040716 - DWK changed char varname1[9] to const string& varname1  
           in public coregerr()
20040716 - DWK changed char varname2[9] to const string& varname2  
           in public coregerr()
20040716 - DWK changed include from tsoldat425.h to tsoldat433.h
20040716 - DWK changed include from telvdat425.h to telvdat433.h
20040716 - DWK changed TEMList51& gtem to TEMList60& gtem in 
           setTEMequilState(), setTEMmiss(), updateTEMmonth() 
           and updateTEMyear()
20040716 - DWK changed InorgN ndep[CYCLE] to InorgN60 ndep[CYCLE]
           in temgisqc()
20040830 - DWK added public string region
20040830 - DWK changed include from ttemdat433.h to ttemdat60.h
20041003 - DWK added spinoutyrs to setTEMequilState(),
           setTEMmiss(), updateTEMmonth(), and updateTEMyear()
20050409 - DWK changed include from ttem60.h to ttem601.h
20050409 - DWK changed include from telm60.cpp tp telm601.cpp at
           bottom of the file
20050622 - DWK changed include from tlcluc60.h to tlcluc601.h           
20050622 - DWK changed include from telm601.cpp tp telm601a.cpp at
           bottom of the file
20051117 - DWK changed include from tclm60.h to tclm602.h
20051117 - DWK changed include from tlcluc601.h to tlcluc602.h
20051117 - DWK changed include from ttem601.h to ttem602.h
20051117 - DWK changed include from tsoldat433.h to tsoldat602.h
20051117 - DWK changed include from telvdat433.h to telvdat602.h
20051117 - DWK changed include from ttemdat60.h to ttemdat602.h
20051117 - DWK deleted telm601a.cpp from bottom of file
20051118 - DWK added public vector<string> predstr
20051118 - DWK added enum outkey
20051118 - DWK added public climate[NUMATMS][CYCLE]
20051119 - DWK added public double mxtair
20051119 - DWK added public double yrprec
20051119 - DWK added int year
20051119 - DWK added public double 
           output[NUMTEM][CYCLE]
20051119 - DWK added public ElmntCohort60 cohort[MAXCHRTS] and 
           associated includes
20051119 - DWK deleted TEMList60& gtem from setTEMequilState(),
           setTEMmiss(), updateTEMmonth() and updateTEMyear()
20051119 - DWK added const int* pichrt to setTEMequilState(),
           setTEMmiss(), updateTEMmonth() and updateTEMyear()
20051122 - DWK added const int& pichrt to public 
           setTEMequilState(), setTEMmiss(), temwritepred(),
           updateTEMmonth() and updateTEMyear()
20051122 - DWK deleted public temwritemiss()
20051122 - DWK changed double nirr[CYCLE] to const double& nirr,
           double par[CYCLE] to const double& par,
	   double tair[CYCLE] to const double& tair,
           double co2[CYCLE] to const double& co2,
           double d40[CYCLE] to const double& aot40,
           InorgN60 ndep[CYCLE] to const InorgN60& ndep in
           private temgisqc() 
20051122 - DWK deleted double prec[CYCLE] from private temgisqc()
20051122 - DWK added const double& rain and 
           const double& snowfall to private temgisqc()
20051122 - DWK added public functions readCohortState(),
           getTEMCohortState(), saveTEMCohortState()
           and writeCohortState()
20051123 - DWK added public function TEMequilibrium()
20051123 - DWK added public function  setMonth()
20051126 - DWK changed Biomass plant[CYCLE] to 
           double plantc[CYCLE] in function call to tranqc()
20051128 - DWK added public function initializeCohortTEMState()
20051129 - DWK deleted const vector<string>& predmap,
           ofstream ftemout[NUMTEM], const int& istateflag, 
           const int& ostateflag and ofstream& ofstate from 
           function call to setTEMequilState()                                                                                   
20051129 - DWK added public function setCohortState()
20051129 - DWK deleted public function  updateTEMyear()           
20051129 - DWK added public int wrtyr
20051129 - DWK deleted const vector<string>& predmap,
           ofstream ftemout[NUMTEM], const int& ostateflag,
           ofstream& ofstate and const int& stateyear from 
           function call to updateTEMmonth()
20051130 - DWK deleted const int& spinoutfg and
           const int& spinoutyrs from function call to 
           setTEMmiss(), setTEMequilState() and updateTEMmonth()
20060124 - DWK added public double avetair
20060607 - DWK changed public float col to double col
20060607 - DWK changed public float row to double row
20060909 - DWK changed include from ttem602.h to ttem603.h
20070514 - DWK changed include tlcluc602.h to tlcluc603a.h  
20070514 - DWK changed include ttemdat602.h to ttemdat603a.h  
20070514 - DWK changed include telmntcohort602.h to 
           telmntcohort603a.h
20070829 - DWK changed include from ttem603.h to ttem603b.h
20070830 - DWK changed include from ttem603b.h to ttem603c.h
20070830 - DWK changed include from telmntcohort603a.hpp to
           telmntcohort603c.hpp
20070830 - DWK changed include from tlcluc603a.h to tlcluc603c.h
20090127 - DWK changed include from ttem603c.h to ttem604.h
		    
*****************************************************************
************************************************************** */

#ifndef TELM604_H
#define TELM604_H

const int TQCZEROFLAG = 31;


//Modules representing climate and TEM

#include "tclm602.h"   // TEMelmnt60 uses the TEMclm60 class

// TEMelmnt60 uses the Elmntcohort60 
#include "telmntcohort603c.hpp"  

#include "tlcluc603c.h" // TEMelmnt60 uses the TEMlcluc60 class
#include "ttem604.h"   // TEMelmnt60 uses the TTEM60 class

// Modules describing the interface with spatially explicit 
//   data sets

#include "tclmdat602.h"  //TEMelmnt60 uses the Clmdata60 class
#include "tsoldat602.h"  //TEMelmnt60 uses the Soildata60 class
#include "telvdat602.h"  //TEMelmnt60 uses the Elevdata60 class
#include "ttemdat603a.h"  //TEMelmnt60 uses the Temdata60 class

class TEMelmnt60
{

  public:

     TEMelmnt60();

     enum outkey
     {
       I_VEGC,     I_SOLC,     I_NSOLC,     I_DOC,      I_TSOLC,       
       I_CO2G,     I_CO2W,     I_HCO3,      I_RHCO3,    I_ALK,      
       I_CROPC,    I_NATVEGC,  I_AGPRDC,    I_PROD10C,  I_PROD100C,  
       I_TOTPRDC,  I_RESIDC,   I_TOTEC,     I_TOTGC,    I_TOTC,     
       
       I_INGPP,    I_GPP,      I_FOZONE,    I_FINDOZONE, I_INNPP,     
       I_NPP,      I_GPR,      I_RVMNT,     I_RVGRW,     I_ABVGPR,    
       I_ROOTGPR,  I_ISOPREN,  I_TERPEN,    I_ORVOC,     I_OVOC,      
       I_VOC,      I_LTRC,     I_CDCMP,     I_RH,        I_RSOIL,     
       I_DOCP,      I_LCHDOC,  I_ERDPOC,    I_NEP,       I_CO2DISS,   
       I_LCHCO2,    I_HCO3P,   I_RHCO3P,    I_LCHHCO3,   I_LCHALK,   
       
       I_CNVRTC,   I_VCNVRTC,  I_SCNVRTC,   I_SLASHC,   I_CFLX,     
       I_AGSTUBC,  I_AGFPRDC,  I_PRDF10C,   I_PRDF100C, I_TOTFPRDC, 
       I_FRESIDC,  I_AGPRDFC,  I_PRD10FC,   I_PRD100FC, I_TOTPRDFC, 
       I_RESIDFC,  I_NCE,      I_NTCS,

       I_AGINGPP,  I_NATINGPP, I_AGGPP,     I_NATGPP,   I_AGINNPP,
       I_NATINNPP, I_AGNPP,    I_NATNPP,    I_AGGPR,    I_NATGPR,
       I_AGRVMNT,  I_NATRVMNT, I_AGRVGRW,   I_NATRVGRW, I_AGLTRC,
       I_NATLTRC,

       I_VEGN,     I_STRN,     I_STON,      I_SOLN,     I_NSOLN,    
       I_DON,      I_TSOLN,    I_AVLN,      I_NH4,      I_NO3,
       I_CROPN,    I_NATVEGN,  I_CSTRN,     I_NATSTRN,  I_CSTON,    
       I_NATSTON,  I_AGPRDN,   I_PROD10N,   I_PROD100N, I_TOTPRDN,  
       I_RESIDN,

       I_NINP,     I_BNFIX,    I_SNFIX,     I_ANFIX,    I_INNUP,    
       I_INNH4UP,  I_INNO3UP,  I_VNUP,      I_VNH4UP,   I_VNO3UP,   
       I_VSUP,     I_VLUP,     I_VNMBL,     I_VNRSRB,   I_LTRN,     
       I_NDCMP,    I_DONP,     I_GMIN,      I_NH4IMM,   I_NO3IMM,   
       I_NIMM,     I_NMIN,     I_AIMMNH4,   I_AMMN,     I_AIMMNO3,  
       I_NTRF,     I_NO3P,     I_NOP,       I_N2OP,     I_N2P,      
       I_DNTRF,    I_NLST,     I_NH3FLX,    I_NOFLX,    I_N2OFLX,   
       I_N2FLX,    I_LCHNO3,   I_LCHDON,    I_ERDPON,   I_NTNS,
      
       I_CNVRTN,   I_VCNVRTN,  I_SCNVRTN,   I_SLASHN,   I_NRETNT,   
       I_NVRTNT,   I_NSRTNT,   I_AGFRTN,    I_AGSTUBN,  I_AGFPRDN,  
       I_PRDF10N,  I_PRDF100N, I_TOTFPRDN,  I_FRESIDN,  I_AGPRDFN,  
       I_PRD10FN,  I_PRD100FN, I_TOTPRDFN,  I_RESIDFN,
      
       I_AGSNFX,   I_NATSNFX,  I_AGINNUP,   I_NATINNUP, I_AINNH4UP,
       I_NINNH4UP, I_AINNO3UP, I_NINNO3UP,  I_AGVNUP,   I_NATVNUP,
       I_AGVNH4UP, I_NVNH4UP,  I_AGVNO3UP,  I_NVNO3UP,  I_AGVSUP,
       I_NATVSUP,  I_AGVLUP,   I_NATVLUP,   I_AGVNMBL,  I_NATVNMBL,
       I_AGVNRSRB, I_NVNRSRB,  I_AGLTRN,    I_NATLTRN,

       I_UNRMLF,   I_LEAF,     I_LAI,       I_FPC,

       I_CROPULF,  I_NATULF,   I_CROPLEAF,  I_NATLEAF, I_CROPLAI,
       I_NATLAI,   I_CROPFPC,  I_NATFPC,

       I_AVLW,     I_SM,       I_VSM,       I_PCTP,     I_SNWPCK,
       I_RGRW,     I_SGRW,

       I_SNWINF,   I_AGIRRIG,  I_PET,       I_INEET,    I_EET,      
       I_RPERC,     I_SPERC,   I_RRUN,      I_SRUN,     I_WYLD,

       I_TSOIL,    I_DST0,     I_DST5,      I_DST10,    I_DST20,
       I_DST50,    I_DST100,   I_DST200,    I_DST300, I_FRONTD,   I_THAWBE,
       I_THAWEND,  I_THAWPCT,  I_ACTLAYER, I_DEADWOODC, I_DEADWOODN
     }; 

/* *************************************************************
		 Public Function Declarations
************************************************************* */

     void atmswritemiss( ofstream fout[NUMATMS],
                         const vector<string>& predname,
                         const int& pdyr,
                         const int& atmspred,
                         const double value );

     void atmswritepred( ofstream fout[NUMATMS],
                         const vector<string>& predname,
                         const int& atmspred );

     int coregerr( ofstream& rflog1,
                   const string& varname1,
                   const double& col1,
		   const double& row1,
                   const string& varname2,
                   const double& col2,
		   const double& row2 );

     int equilibrateTEM( const int& pichrt,
                         const double& ptol, 
                         ofstream& rflog1 );

     void getTEMCohortState( const int& pichrt,
                             const int& pdm );
     void getTEMCohortState0( const int& pichrt,
                             const int& pdm );

     void initializeCohortTEMState( const int& pichrt );

     void outputTEMmonth( const int& pdm );

     void readCohortState( ifstream& ifstate,
                           const int& pichrt );
     void readCohortState0( ifstream& ifstate,
                           const int& pichrt );
     void saveTEMCohortState( const int& pichrt );

     void setCohortTEMState( const ElmntCohort60& firstchrt,
                             ElmntCohort60& targetchrt );
     
     int setGIStopography( ofstream& flog1,
                           int& ftlerr,
                           FILE* fstxt,
                           FILE*felev );
          
     void setTEMequilState( ofstream& rflog1,
                            const int& equil,
                            const int& totsptime,
                            const int& pichrt );

     void setTEMmiss( const int& pdyr,
                      const int& equil,
                      const int& totsptime,
                      const int& pichrt );

     void temwritepred( ofstream fout[NUMTEM],
                        const vector<string>& predname,
                        const int& pdyr,
                        const int& pichrt,
                        const int& ntempred );

     void updateTEMmonth( ofstream& rflog1,
                          const int& equil,
                          const int& totsptime,
                          const int& outyr,
                          const int& pdm,
                          const int& pichrt );

     void writeCohortState( ofstream& ofstate,
                            const int& pichrt );

/* *************************************************************
		 Public Variables
************************************************************* */

     int atmstotyr[MAXRTIME];
     
     // Mean annual air temperature
     double avetair;

     double carea;

     double climate[NUMATMS][CYCLE+1];
     
     TEMclm60 clm;

     ElmntCohort60 cohort[MAXCHRTS];

     Soildata60 fao;

     Elevdata60 elv;
     
     double col;

     string contnent;
     
     int fatalerr;

     double lat;

     TEMlcluc60 lcluc;

     double lon;

     int lonlatflag;
     
     // Maximum number of cohorts in an element
     //   during the current year
     int maxcohorts;
     
     // Maximum monthly air temperature during the year
     double mxtair;

     // Maximum number of "natural" cohorts in an element
     int natcohorts;

     // Number of climate predicted variables
     int natmspred;
     
     // Number of TEM predicted variables
     int ntempred;

     double output[NUMTEM][CYCLE];
     
//     ofstream ofstateA;
//     ofstream ofstateB;

     int outyr;

     // Maximum number of cohorts in an element
     //   during the previous year
     int prvmxcohrts;

     string region;
     
     double row;

     double subarea;

     TTEM60 tem;

//     string temstateAfname;
//     string temstateBfname;

     int totpred;

     int ttotyr[MAXRTIME];

     int wrtyr;

     // Annual precipitation during the year
     double yrprec;
     
     int year;
     //double dumflt1;
     //double dumflt2;

/* *************************************************************
		 Private Function Declarations
************************************************************* */

  private:

     int temgisqc( const double& subarea,
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
                   const InorgN60& ndep );

     int transqc( int& maxyears,
                  int& totyr,
                  double plantc[CYCLE] );

};

#endif
