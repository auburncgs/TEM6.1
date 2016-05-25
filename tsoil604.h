/* **************************************************************
*****************************************************************
TSOIL604.H - object describing characteristics of soil used by
	          the Terrestrial Ecosystem Model (TEM)

Modifications:

20030718 - DWK created by modifying tsoil50b1.h to incorporate
           open-N cycle dynamics into TEM
20030718 - DWK added include inorgn51.hpp
20030718 - DWK changed class Tsoil50 to Tsoil60
20030718 - DWK changed include from tsoil50b1.cpp to tsoil60.cpp
           at bottom of file
20030718 - DWK changed double availn[CYCLE] to InorgN availn[CYCLE]
20030718 - DWK changed double yravln to InorgN yravln
20030718 - DWK added parameters lchNO3par, lchNO3parcut[MAXCMNT],
           lchNO3par1a[MAXCMNT], lchNO3par1b[MAXVMNT],
           lchNO3par2a[MAXCMNT], lchNO3par2b[MAXCMNT],
           DOCpar[MAXCMNT], DONpar[MAXCMNT] and lchDOMpar[MAXCMNT]
20030718 - DWK deleted public int tndepflag
20030718 - DWK deleted public nintot[MAXCMNT] and nloss[MAXCMNT]
20030720 - DWK added public double leachNO3[CYCLE], double
           yrlchNO3, Biomass leachDOM[CYCLE] and Biomass yrlchDOM
20030721 - DWK added public double nh3flux[CYCLE], noflux[CYCLE],
           n2oflux[CYCLE], double n2flux[CYCLE], double yrnh3flx,
           double yrnoflx, double yrn2oflx, and  double yrn2flx
20030721 - DWK added public function setTraceGasFluxes()
20030731 - DWK added public double abioticNimmob[CYCLE] and
           yrabNimmob
20030731 - DWK added public function setDOMleaching()
20030731 - DWK added public double maxdocCN[MAXCMNT] and
           double mindocCN[MAXCMNT]
20030826 - DWK added public Biomass DOM[CYCLE] and Biomass yrDOM
20031015 - DWK added double gasCO2[CYCLE], double yrgasCO2,
           double aqueousCO2[CYCLE] double yraqCO2,
           double HCO3[CYCLE], double yrHCO3,
           double CO3[CYCLE], double yrCO3,
           double alkalinity[CYCLE], and double yralk
20031015 - double leachALK[CYCLE], double yrlchALK,
           double dissolveCO2[CYCLE], double yrdissCO2,
           double leachCO2[CYCLE], double yrlchCO2,
           double formCO3[CYCLE], double yrformCO3,
           double leachCO3[CYCLE], double yrlchCO3,
           double formHCO3[CYCLE], double yrformHCO3,
           double leachHCO3[CYCLE], double yrlchHCO3,
           Biomass erodePOM[CYCLE] and Biomass yrerodePOM
20031015 - DWK replaced char ecd[MAXFNAME] with string ecd in
           function calls to getecd() and getrootz()
20031016 - DWK changed include from qsoiltemp50b1.h to
           qsoiltemp51.h
20031016 - DWK changed inheritance of Soilthermal to 
           Soilthermal51
20040707 - DWK changed include from qsoiltemp51.h 
           to qsoiltemp60.h
20040707 - DWK changed inheritance of Soilthermal51 to 
           Soilthermal60
20040716 - DWK changed include from inorgn51.hpp to inorgn60.hpp
20040716 - DWK changed public InorgN availn[CYCLE] to
           InorgN60 availn[CYCLE]
20040716 - DWK changed public InorgN yravln to InorgN60 yravln 
20040919 - DWK changed public double DONpar[MAXCMNT] to
           double DONpar and added double DONpara[MAXCMNT]
           and double DONparb[MAXCMNT]
20051117 - DWK changed include from qsoiltemp60.h to 
           qsoiltemp602.h 
20051117 - DWK deleted tsoil60.cpp from bottom of file 
20051118 - DWK changed public double CO3[] to RHCO3[]
20051118 - DWK deleted public double leachCO3[]
20051125 - DWK changed public double avlh2o[CYCLE] to 
           double avlh2o
20051125 - DWK changed public double h2oyld[CYCLE] to 
           double h2oyld
20051125 - DWK changed public double moist[CYCLE] to 
           double moist
20051125 - DWK changed public double abioticNimmob[CYCLE] to 
           double abioticNimmob
20051125 - DWK changed public double alkalinity[CYCLE] to 
           double alkalinity
20051125 - DWK changed public double aqueousCO2[CYCLE] to 
           double aqueousCO2
20051125 - DWK changed public InorgN60 availn[CYCLE] to 
           InorgN60 availn
20051125 - DWK changed public double dissolveCO2[CYCLE] to 
           double dissolveCO2
20051125 - DWK changed public Biomass DOM[CYCLE] to 
           Biomass DOM 
20051125 - DWK changed public Biomass erodePOM[CYCLE] to 
           Biomass erodePOM
20051125 - DWK changed public double formHCO3[CYCLE] to 
           double formHCO3
20051125 - DWK changed public double formRHCO3[CYCLE] to 
           double formRHCO3
20051125 - DWK changed public double gasCO2[CYCLE] to 
           double gasCO2
20051125 - DWK changed public double HCO3[CYCLE] to 
           double HCO3
20051125 - DWK changed public double leachALK[CYCLE] to 
           double leachALK 
20051125 - DWK changed public double leachCO2[CYCLE] to 
           double leachCO2 
20051125 - DWK changed public Biomass leachDOM[CYCLE] to 
           Biomass leachDOM
20051125 - DWK changed public double leachHCO3[CYCLE] to 
           double leachHCO3
20051125 - DWK changed public double leachNO3[CYCLE] to 
           double leachNO3
20051125 - DWK changed public double n2flux[CYCLE] to 
           double n2flux
20051125 - DWK changed public double n2oflux[CYCLE] to 
           double n2oflux 
20051125 - DWK changed public double nh3flux[CYCLE] to 
           double nh3flux
20051125 - DWK changed public double ninput[CYCLE] to 
           double ninput
20051125 - DWK changed public double nlost[CYCLE] to 
           double nlost
20051125 - DWK changed public double noflux[CYCLE] to 
           double noflux
20051125 - DWK changed public Biomass nonReactiveOrg[CYCLE] to 
           Biomass nonReactiveOrg
20051125 - DWK changed public Biomass org[CYCLE] to Biomass org
20051125 - DWK changed public double pcfc[CYCLE] to double pcfc          
20051125 - DWK changed public double pctp[CYCLE] to double pctp
20051125 - DWK changed public double rgrndh2o[CYCLE] to 
           double rgrndh2o
20051125 - DWK changed public double RHCO3[CYCLE] to 
           double RHCO3           
20051125 - DWK changed public double rperc[CYCLE] to 
           double rperc
20051125 - DWK changed public double rrun[CYCLE] to double rrun                      
20051125 - DWK changed public double sgrndh2o[CYCLE] to 
           double sgrndh2o
20051125 - DWK changed public double snowinf[CYCLE] to 
           double snowinf
20051125 - DWK changed public double snowpack[CYCLE] to 
           double snowpack
20051125 - DWK changed public double sperc[CYCLE] to 
           double sperc
20051125 - DWK changed public double srun[CYCLE] to double srun                                           
20051125 - DWK changed public Biomass totalOrg[CYCLE] to 
           Biomass totalOrg
20051125 - DWK changed public double vsm[CYCLE] to double vsm           
20051125 - DWK deleted const int& dm from function call to
           percol(), setDOMleaching(), setTraceGasFluxes() 
           and updateActiveLayerRootZ()
20051130 - DWK added public double eet, double eetmx, 
           double ineet, double prveetmx, double yreet and 
           double yrineet
20051130 - DWK added public function xeet()
20051130 - DWK added public function updateHydrology()
20051130 - DWK added private double ndays[CYCLE]
20051130 - DWK deleted rperc from rrunoff()
20051130 - DWK deleted sperc from srunoff()
20051130 - DWK deleted const double& snowinf and
           const double& eet from percol()
20051130 - DWK deleted public parameters double DOCpar[MAXCMNT],
           double DONpar, double DONparcut[MAXCMNT],
           double DONpar1a[MAXCMNT], double DONpar1b[MAXCMNT],
           double DONpar2a[MAXCMNT], and 
           double DONpar2b[MAXCMNT]
20051201 - DWK added public double kh2o
20051201 - DWK added public function setKH2O()
20051201 - DWK added public function updateLeachingLosses()
20051202 - DWK added public function resetMonthlyFluxes()
20051202 - DWK added public function resetYrFluxes()
20051208 - DWK added public function resetEcds()
20051213 - DWK deleted inheritance of Soilthermal60 class
20051213 - DWK added public Soilthermal60 stm
20051213 - DWK added private double activeLayer
20051213 - DWK added include tprocessXML602.h
20051213 - DWK inherits class ProcessXML60
20060608 - DWK added private double drainage and double surfrun
20060909 - DWK added private parameters 
           double propReactA[MAXCMNT] and 
           double propReactB[MAXCMNT]
20060909 - DWK added public function 
           getThawedReactiveSOMProp()
20060909 - DWK added getPROPREACTA() and setPROPREACTA()
20060909 - DWK added getPROPREACTB() and setPROPREACTB()
20070411 - DWK changed include from qsoiltemp602.h to 
           qsoiltemp603.h
20070829 - DWK added private double evaporation
20070830 - DWK added private double wfpsoff
                                            
****************************************************************
************************************************************* */

#ifndef TSOIL604_H
#define TSOIL604_H

// Tsoil60 also uses the global constants CYCLE, MAXRTIME, MISSING,
//   NUMVEG and MAXFNAME
#include "temconsts602.hpp"

#include "tprocessXML602.h"

#include "bioms423.hpp"      // Tsoil60 uses Biomass class
#include "inorgn60.hpp"      // Tsoil60 used InorgN60 class
#include "qsoiltemp603.h"    // Tsoil60 inherits Soilthermal60 class


class Tsoil60 : public ProcessXML60 
{

  public:
    
     Tsoil60( void );

/* **************************************************************
		 Public Functions
************************************************************** */

     void getecd ( ofstream& rflog1 );

     void getecd ( const string& ecd );
     
     void getrootz( ofstream& rflog1 );

     void getrootz( const string& ecd );

     double getThawedReactiveSOMProp( const int& pdcmnt );

     void lake( const double& tair,
                const double& prec,
                double& rain,
                double& snowfall,
                const double& pet,
                double& eet );

     void resetEcds( const int& pcmnt );
     
     void resetMonthlyFluxes( void );
     
     void resetYrFluxes( void );

     void setDOMleaching( const int& pdcmnt,
                          const double& soilc,
                          const double& soiln,
                          const double& soilh2o,
                          const double& rain,
                          const double& eet );

     void setKH2O( const double& vsm,
                   const int& moistlim );

     void setTraceGasFluxes( const double& nh3vol,
                             const double& noprd,
                             const double& n2oprd,
                             const double& n2prd  );

     void showecd( void );

     double snowmelt( const double& elev,
                      const double& tair,
                      const double& prevtair,
                      const double& snowpack );

     void updateActiveLayerRootZ( void );
     
     void updateHydrology( const double& elev,
                           const double& tair,
                           const double& prevtair,
                           const double& prev2tair,
                           const double& rain,
                           const double& pet,
                           const double& avlh2o,
                           const double& rgrndh2o,
                           const double& sgrndh2o,
                           const int& irrgflag,
                           double& irrigate,
                           const int& pdm );

     void updateLeachingLosses( const int& pdcmnt,
                                const double& doc,
                                const double& don,
                                const double& soilno3, 
                                const double& soilh2o );
                                    
     void updateRootZ( const int& pdcmnt );

     void xtext( const int& pdcmnt, 
                 const double& pctsilt, 
                 const double& pctclay );


     // "Get" and "Set" private variables and parameters
     
     // abioticNimmob ******************************************
     
     inline double getABIMMOB( void ) { return abioticNimmob; }


     // activeLayer ********************************************
      
     inline double getACTLAYER( void ) { return activeLayer; }

     inline void setACTLAYER( const double& pactlayer ) 
     { 
       activeLayer = pactlayer; 
     }


     // alkalinity *********************************************
     
     inline double getALK( void ) { return alkalinity; }


     // aqueousCO2 *********************************************
     
     inline double getAQCO2( void ) { return aqueousCO2; }
     
     
     // availn.nh4 *********************************************
     
     inline double getNH4( void ) { return availn.nh4; }


     // availn.no3 *********************************************
     
     inline double getNO3( void ) { return availn.no3; }


     // availn.total *******************************************
     
     inline double getAVLN( void ) { return availn.total; }

     inline void setAVLN( const double& pavln ) 
     { 
       availn.total = pavln; 
     }


     // avlh2o *************************************************
     
     inline double getAVLH2O( void ) { return avlh2o; }

     inline void setAVLH2O( const double& pavlh2o ) 
     { 
       avlh2o = pavlh2o; 
     }

     // awcapmm ************************************************
     
     inline double getAWCAPMM( void ) { return awcapmm; }

     inline void setAWCAPMM( const double& pawcapmm ) 
     { 
       awcapmm = pawcapmm; 
     }


     // dissolveCO2 ********************************************
     
     inline double getDISSCO2( void ) { return dissolveCO2; }


     // DOM.carbon *********************************************
     
     inline double getDOC( void ) { return DOM.carbon; }
     
     
     // DOM.nitrogen *******************************************
     
     inline double getDON( void ) { return DOM.nitrogen; }
    

     // drainage ***********************************************
     
     inline double getDRAINAGE( void ) { return drainage; }

     inline void setDRAINAGE( const double& pdrainage ) 
     { 
       drainage = pdrainage; 
     }


     // dst10 **************************************************
     
     inline double getDST10( void ) { return dst10; }

     inline void setDST10( const double& pdst10 ) 
     { 
       dst10 = pdst10; 
     }
    
    
     // eet ****************************************************
     
     inline double getEET( void ) { return eet; }

     inline void setEET( const double& peet ) 
     { 
       eet = peet; 
     }


     // eetmx **************************************************
     
     inline double getEETMX( void ) { return eetmx; }

     inline void setEETMX( const double& peetmx ) 
     { 
       eetmx = peetmx; 
     }


     // erodePOM.carbon ****************************************
     
     inline double getERODEPOC( void ) 
     { 
       return erodePOM.carbon; 
     }

     inline void setERODEPOC( const double& perodepoc ) 
     { 
       erodePOM.carbon = perodepoc; 
     }


     // erodePOM.nitrogen **************************************
     
     inline double getERODEPON( void ) 
     { 
       return erodePOM.nitrogen; 
     }

     inline void setERODEPON( const double& perodepon ) 
     { 
       erodePOM.nitrogen = perodepon; 
     }


     // evaporation **************************************
     
     inline double getEVAPORATION( void ) 
     { 
       return evaporation; 
     }

     inline void setEVAPORATION( const double& pevap ) 
     { 
       evaporation = pevap; 
     }


     // fldcap *************************************************
     
     inline double getFLDCAP( void ) { return fldcap; }


     // fldcapa ************************************************
     
     inline double getFLDCAPA( void ) { return fldcapa; }

     inline void setFLDCAPA( const double& pfldcapa ) 
     { 
       fldcapa = pfldcapa; 
     }


     // fldcapb ************************************************
     
     inline double getFLDCAPB( void ) { return fldcapb; }

     inline void setFLDCAPB( const double& pfldcapb ) 
     { 
       fldcapb = pfldcapb; 
     }


     // formHCO3 ***********************************************
     
     inline double getFORMHCO3( void ) { return formHCO3; }


     // formRHCO3 **********************************************
     
     inline double getFORMRHCO3( void ) { return formRHCO3; }


     // gasCO2 *************************************************
     
     inline double getGASCO2( void ) { return gasCO2; }

     
     // h2oyld *************************************************
     
     inline double getH2OYLD( void ) { return h2oyld; }

     inline void setH2OYLD( const double& ph2oyld ) 
     { 
       h2oyld = ph2oyld; 
     }

     
     // HCO3 ***************************************************
     
     inline double getHCO3( void ) { return HCO3; }


     // ineet **************************************************
     
     inline double getINEET( void ) { return ineet; }

     inline void setINEET( const double& pineet ) 
     { 
       ineet = pineet; 
     }


     // kh2o ***************************************************
     
     inline double getKH2O( void ) { return kh2o; }


     // lchDOMpar **********************************************
     
     inline double getLCHDOMPAR( const int& pcmnt ) 
     { 
       return lchDOMpar[pcmnt]; 
     }
     
     inline void setLCHDOMPAR( const double& plchdompar, 
                               const int& pcmnt ) 
     { 
       lchDOMpar[pcmnt] = plchdompar; 
     }


     // lchNO3par **********************************************
     
     inline double getLCHNO3PAR( void ) { return lchNO3par; }

     inline void setLCHNO3PAR( const double& plchno3par ) 
     { 
       lchNO3par = plchno3par; 
     }


     // lchNO3parcut *******************************************
     
     inline double getLCHNO3PARCUT( const int& pcmnt ) 
     { 
       return lchNO3parcut[pcmnt]; 
     }
     
     inline void setLCHNO3PARCUT( const double& plchno3parcut, 
                                  const int& pcmnt ) 
     { 
       lchNO3parcut[pcmnt] = plchno3parcut; 
     }

     // lchNO3par1a ********************************************
     
     inline double getLCHNO3PAR1A( const int& pcmnt ) 
     { 
       return lchNO3par1a[pcmnt]; 
     }
     
     inline void setLCHNO3PAR1A( const double& plchno3par1a, 
                                 const int& pcmnt ) 
     { 
       lchNO3par1a[pcmnt] = plchno3par1a; 
     }

     // lchNO3par1b ********************************************

     inline double getLCHNO3PAR1B( const int& pcmnt ) 
     { 
       return lchNO3par1b[pcmnt]; 
     }
     
     inline void setLCHNO3PAR1B( const double& plchno3par1b, 
                                 const int& pcmnt ) 
     { 
       lchNO3par1b[pcmnt] = plchno3par1b; 
     }

     // lchNO3par2a ********************************************

     inline double getLCHNO3PAR2A( const int& pcmnt ) 
     { 
       return lchNO3par2a[pcmnt]; 
     }
     
     inline void setLCHNO3PAR2A( const double& plchno3par2a, 
                                 const int& pcmnt ) 
     { 
       lchNO3par2a[pcmnt] = plchno3par2a; 
     }

     // lchNO3par2b ********************************************

     inline double getLCHNO3PAR2B( const int& pcmnt ) 
     { 
       return lchNO3par2b[pcmnt]; 
     }
     
     inline void setLCHNO3PAR2B( const double& plchno3par2b, 
                                 const int& pcmnt ) 
     { 
       lchNO3par2b[pcmnt] = plchno3par2b; 
     }

     
     // leachALK ***********************************************
     
     inline double getLEACHALK( void ) { return leachALK; }
     
     
     // leachCO2 ***********************************************

     inline double getLEACHCO2( void ) { return leachCO2; }
     
     
     // leachDOM.carbon ****************************************

     inline double getLEACHDOC( void ) 
     { 
       return leachDOM.carbon; 
     }

     inline void setLEACHDOC( const double& pleachdoc )
     { 
       leachDOM.carbon = pleachdoc; 
     }


     // leachDOM.nitrogen **************************************
     
     inline double getLEACHDON( void ) 
     { 
       return leachDOM.nitrogen; 
     }

     inline void setLEACHDON( const double& pleachdon ) 
     { 
       leachDOM.nitrogen = pleachdon; 
     }


     // leachHCO3 **********************************************
     
     inline double getLEACHHCO3( void ) { return leachHCO3; }


     // leachNO3 ***********************************************
     
     inline double getLEACHNO3( void ) { return leachNO3; }

     inline void setLEACHNO3( const double& pleachno3 ) 
     { 
       leachNO3 = pleachno3; 
     }


     // minrootz ***********************************************
     
     inline double getMINROOTZ( const int& pcmnt ) 
     { 
       return minrootz[pcmnt]; 
     }

     inline void setMINROOTZ( const double& pminrootz,
                              const int& pcmnt ) 
     { 
       minrootz[pcmnt] = pminrootz; 
     }


     // moist **************************************************
     
     inline double getMOIST( void ) { return moist; }

     inline void setMOIST( const double& psh2o ) 
     { 
       moist = psh2o; 
     }


     // n2flux *************************************************
     
     inline double getN2FLUX( void ) { return n2flux; }


     // n2oflux ************************************************
     
     inline double getN2OFLUX( void ) { return n2oflux; }


     // nextdst10 **********************************************
     
     inline double getNEXTDST10( void ) { return nextdst10; }

     inline void setNEXTDST10( const double& pnextdst10 ) 
     { 
       nextdst10 = pnextdst10; 
     }


     // nh3flux ************************************************
     
     inline double getNH3FLUX( void ) { return nh3flux; }

     inline void setNH3FLUX( const double& pnh3flx ) 
     { 
       nh3flux = pnh3flx; 
     }


     // ninput *************************************************
     
     inline double getNINPUT( void ) { return ninput; }

     inline void setNINPUT( const double& pninput ) 
     { 
       ninput = pninput; 
     }


     // nlost **************************************************
     
     inline double getNLOST( void ) { return nlost; }

     inline void setNLOST( const double& pnlst ) 
     { 
       nlost = pnlst; 
     }


     // noflux *************************************************
     
     inline double getNOFLUX( void ) { return noflux; }


     // nonReactiveOMpar ***************************************
     
     inline double getNSOLPAR( const int& pcmnt ) 
     { 
       return nonReactiveOMpar[pcmnt]; 
     }
     
     inline void setNSOLPAR( const double& pnsolcpar, 
                             const int& pcmnt ) 
     { 
       nonReactiveOMpar[pcmnt] = pnsolcpar; 
     }


     // nonReactiveOrg.carbon **********************************
     
     inline double getNSOLC( void ) 
     { 
       return nonReactiveOrg.carbon; 
     }

     inline void setNSOLC( const double& pnsolc ) 
     { 
       nonReactiveOrg.carbon = pnsolc; 
     }


     // nonReactiveOrg.nitrogen ********************************

     inline double getNSOLN( void ) 
     { 
       return nonReactiveOrg.nitrogen; 
     }

     inline void setNSOLN( const double& pnsoln ) 
     { 
       nonReactiveOrg.nitrogen = pnsoln; 
     }


     // pctclay ************************************************
     
     inline double getPCTCLAY( void ) { return pctclay; }

     inline void setPCTCLAY( const double& ppctclay ) 
     { 
       pctclay = ppctclay; 
     }


     // pcfldcap ***********************************************
     
     inline double getPCTFLDCAP( void ) { return pcfldcap; }

     inline void setPCTFLDCAP( const double& ppctfldcap ) 
     { 
       pcfldcap = ppctfldcap; 
     }

     
     // pctp ***************************************************
     
     inline double getPCTP( void ) { return pctp; }

     inline void setPCTP( const double& ppctp ) 
     { 
       pctp = ppctp; 
     }


     // pctpora ************************************************
     
     inline double getPCTPORA( void ) { return pctpora; }

     inline void setPCTPORA( const double& ppctpora ) 
     { 
       pctpora = ppctpora; 
     }


     // pctporb ************************************************
     
     inline double getPCTPORB( void ) { return pctporb; }

     inline void setPCTPORB( const double& ppctporb ) 
     { 
       pctporb = ppctporb; 
     }


     // pctsand ************************************************
     
     inline double getPCTSAND( void ) { return pctsand; }

     inline void setPCTSAND( const double& ppctsand ) 
     { 
       pctsand = ppctsand; 
     }


     // pctsilt ************************************************
     
     inline double getPCTSILT( void ) { return pctsilt; }

     inline void setPCTSILT( const double& ppctsilt ) 
     { 
       pctsilt = ppctsilt; 
     }


     // prevdst10 **********************************************
     
     inline double getPREVDST10( void ) { return prevdst10; }

     inline void setPREVDST10( const double& pprevdst10 ) 
     { 
       prevdst10 = pprevdst10; 
     }


     // prevspack **********************************************
     
     inline double getPREVSPACK( void ) { return prevspack; }

     inline void setPREVSPACK( const double& pprvspack ) 
     { 
       prevspack = pprvspack; 
     }


     // propReactA *************************************************
     
     inline double getPROPREACTA( const int& pcmnt ) 
     { 
       return propReactA[pcmnt]; 
     }

     inline void setPROPREACTA( const double& ppropreacta,
                                const int& pcmnt ) 
     { 
       propReactA[pcmnt] = ppropreacta; 
     }


     // propReactB *************************************************
     
     inline double getPROPREACTB( const int& pcmnt ) 
     { 
       return propReactB[pcmnt]; 
     }

     inline void setPROPREACTB( const double& ppropreactb,
                                const int& pcmnt ) 
     { 
       propReactB[pcmnt] = ppropreactb; 
     }


     // prveetmx ***********************************************
     
     inline double getPRVEETMX( void ) { return prveetmx; }

     inline void setPRVEETMX( const double& pprveetmx ) 
     { 
       prveetmx = pprveetmx; 
     }


     // psiplusc ***********************************************
     
     inline double getPSIPLUSC( void ) { return psiplusc; }

     inline void setPSIPLUSC( const double& ppsiplusc ) 
     { 
       psiplusc = ppsiplusc; 
     }


     // rootz **************************************************
     
     inline double getROOTZ( void ) { return rootz; }

     inline void setROOTZ( const double& prootz ) 
     { 
       rootz = prootz; 
     }


     // rootza *************************************************
     
     inline double getROOTZA( const int& pcmnt ) 
     { 
       return rootza[pcmnt]; 
     }

     inline void setROOTZA( const double& prootza,
                            const int& pcmnt ) 
     { 
       rootza[pcmnt] = prootza; 
     }


     // rootzb *************************************************
     
     inline double getROOTZB( const int& pcmnt ) 
     { 
       return rootzb[pcmnt]; 
     }

     inline void setROOTZB( const double& prootzb,
                            const int& pcmnt ) 
     { 
       rootzb[pcmnt] = prootzb; 
     }


     // rootzc *************************************************
     
     inline double getROOTZC( const int& pcmnt ) 
     { 
       return rootzc[pcmnt]; 
     }

     inline void setROOTZC( const double& prootzc,
                            const int& pcmnt ) 
     { 
       rootzc[pcmnt] = prootzc; 
     }


     // rperc **************************************************
     
     inline double getRPERC( void ) { return rperc; }


     // rrun ***************************************************
     
     inline double getRRUN( void ) { return rrun; }
 
 
     // snowinf ************************************************
     
     inline double getSNOWINF( void ) { return snowinf; }

     inline void setSNOWINF( const double& psnwinf ) 
     { 
       snowinf = psnwinf; 
     }


     // snowpack ***********************************************
     
     inline double getSNOWPACK( void ) { return snowpack; }

     inline void setSNOWPACK( const double& psnwpck ) 
     { 
       snowpack = psnwpck; 
     }


     // sperc **************************************************
     
     inline double getSPERC( void ) { return sperc; }


     // srun ***************************************************
     
     inline double getSRUN( void ) { return srun; }


     // surfrun ************************************************
     
     inline double getSURFRUN( void ) { return surfrun; }

     inline void setSURFRUN( const double& psurfrun ) 
     { 
       surfrun = psurfrun; 
     }


     // totpor *************************************************
     
     inline double getTOTPOR( void ) { return totpor; }

     inline void setTOTPOR( const double& ptotpor ) 
     { 
       totpor = ptotpor; 
     }


     // tsoil **************************************************
     
     inline double getTSOIL( void ) { return tsoil; };

     inline void setTSOIL( const double& ptsoil ) 
     { 
       tsoil = ptsoil; 
     }


     // totalOrg.carbon ****************************************
     
     inline double getTSOLC( void ) { return totalOrg.carbon; }

     inline void setTSOLC( const double& ptsolc ) 
     { 
       totalOrg.carbon = ptsolc; 
     }


     // totalOrg.nitrogen **************************************
     
     inline double getTSOLN( void ) 
     { 
       return totalOrg.nitrogen; 
     }

     inline void setTSOLN( const double& ptsoln ) 
     { 
       totalOrg.nitrogen = ptsoln; 
     }


     // vsm *************************************************
     
     inline double getVSM( void ) { return vsm; }
     
     inline void setVSM( const double& pvsm ) 
     { 
       vsm = pvsm; 
     }


     // wfpsoff *************************************************
     
     inline double getWFPSOFF( void ) 
     { 
       return wfpsoff; 
     }

     inline void setWFPSOFF( const double& pwfpsoff ) 
     { 
       wfpsoff = pwfpsoff; 
     }


     // wiltpt *************************************************
     
     inline double getWILTPT( void ) { return wiltpt; }
     
     inline void setWILTPT( const double& pwiltpt ) 
     { 
       wiltpt = pwiltpt; 
     }


     // wiltpta ************************************************
     
     inline double getWILTPTA( void ) { return wiltpta; }

     inline void setWILTPTA( const double& pwiltpta ) 
     { 
       wiltpta = pwiltpta; 
     }


     // wiltptb ************************************************
     
     inline double getWILTPTB( void ) { return wiltptb; }

     inline void setWILTPTB( const double& pwiltptb ) 
     { 
       wiltptb = pwiltptb; 
     }



     // wsoil **************************************************
     
     inline double getWSOIL( void ) { return wsoil; }

     
/* *************************************************************
		 Public Variables
************************************************************* */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif

     // Soil thermal module
     Soilthermal60 stm;

    // Flag to indicate that soil thermal model should be run:
    //   stmflg = 0 if soil thermal model is not run
    //   stmflg = 1 if soil thermal model is run
    int stmflg; 

     
     // Annual sum of abioticNimmob
     double yrabNimmob;  // (g N / (sq. meter * year))

     // Annual sum of alkalinity
     double yralk;

     // Annual sum of aqueousCO2
     double yraqCO2;

     // Annual sum of avlh2o
     double yravlh2o;

     // Annual sum of availn
     InorgN60 yravln;

     // Ratio of soil reactive organic carbon to 
     //   soil reactive organic nitrogen
     double yrc2n;

     // Annual sum of dissolveCO2
     double yrdissCO2;

     // Annual sum of DOM
     Biomass yrDOM;   // (g C or g N / sq. meter)

     // Annual sum of dst10
     double yrdst10;

     // Annual estimated actual evapotranspiration (mm / year)
     double yreet;

     // Annual sum of erodePOM
     Biomass yrerodePOM;  // (g C or g N / sq. meter * year)

     // Annual sum of formHCO3
     double yrformHCO3;   // (g C / (sq. meter * year))

     // Annual sum of formRHCO3
     double yrformRHCO3;   // (g C / (sq. meter * year))

     // Annual sum of gasCO2
     double yrgasCO2;

     // Annual sum of HCO3
     double yrHCO3;

     // Annual sum of h2oyld (mm / year)
     double yrh2oyld;

     // Annual initial estimated actual evapotranspiration
     //   (mm / year)
     double yrineet;

     // Annual sum of leachALK
     double yrlchALK;   // (g C / (sq. meter * year))

     // Annual sum of leachCO2
     double yrlchCO2;   // (g C / (sq. meter * year))

     // Annual sum of leachDOM
     Biomass yrlchDOM;  // (g C or g N / (sq. meter * year))

     // Annual sum of leachHCO3
     double yrlchHCO3;  // (g C / (sq. meter * year))

     // Annual sum of leachNO3
     double yrlchNO3;   // (g N / (sq. meter * year))

     // Annual sum of n2flux
     double yrn2flx;    // (g N / (sq. meter * year))

     // Annual sum of n2oflux
     double yrn2oflx;  // (g N / (sq. meter * year))

     // Annual sum of nh3flux
     double yrnh3flx;  // (g N / (sq. meter * year))

     // Annual sum of ninput
     double yrnin;     // (g N / (sq. meter * year))

     // Annual sum of nlost
     double yrnlost;   // (g N / (sq. meter * year))

     // Annual sum of noflux
     double yrnoflx;   // (g N / (sq. meter * year))

     // Annual sum of nonReactiveOrg.carbon
     double yrnonorgc;

     // Annual sum of nonReactiveOrg.nitrogen
     double yrnonorgn;

     // Annual sum of org.carbon 
     double yrorgc;
     
     // Annual sum of org.nitrogen 
     double yrorgn;

     // Annual sum of pctp
     double yrpctp;

     // Annual sum of rgrdnh2o
     double yrrgrndh2o;

     // Annual sum of RHCO3
     double yrRHCO3;

     // Annual sum of rperc
     double yrrperc;   // (mm / year)

     // Annual sum of rrun
     double yrrrun;  // (mm / year)

     // Annual sum of sgrndh2o
     double yrsgrndh2o;

     // Annual sum of moist
     double yrsmoist;

     // Annual sum of snowinf
     double yrsnowinf;      // (mm / year)

     // Annual sum of snowpack
     double yrsnowpack;

     // Annual sum of sperc
     double yrsperc; // (mm / year)

     // Annual sum of srun
     double yrsrun;  // (mm / year)

     // Annual sum of totalOrg.carbon
     double yrtotorgc;

     // Annual sum of totalOrg.nitrogen     
     double yrtotorgn;

     // Annual sum of tsoil
     double yrtsoil;

     // Annual sum of vsm
     double yrvsm;
   
   
   private:
   
/* **************************************************************
		 Private Functions
************************************************************** */

     void percol( const double& rain,                  
                  const double& avlh2o );

     double rrunoff( const double& rgrndh2o );

     double srunoff( const double& elev,
                     const double& tair,
                     const double& prevtair,
                     const double& prev2tair,
		     const double& sgrndh2o );

     /* Estimated "actual" evapotransipiration (i.e. EET) as 
        described in Vorosmarty et al. (1989) Global Biogeochemical 
        Cycles 3: 241-265.  */

     double xeet( const double& rain,
                  const double& pet,
                  const double& avlh2o,
                  const int& pdm );


/* **************************************************************
		 Private Variables
************************************************************** */
     
     // Monthly abiotic immobilization
     double abioticNimmob; // (g N / (sq. meter * month))

     // Active Layer Depth in mm
     double activeLayer;
     
     // Soil Alkalinity (g C / sq. meter)
     double alkalinity;

     // Aqueous carbon dioxide (g C / sq. meter)
     double aqueousCO2;

     // Available inorganic nitrogen (g N / sq. meter)
     InorgN60 availn;

     // Available water (mm)
     double avlh2o;  

     // Available water capacity (mm)
     double awcapmm;        

     // Dissolution of soil carbon dioxide 
     double dissolveCO2; // ( g C / (sq. meter * month))

     // Dissolved organic matter (DOM)
     Biomass DOM;  // (g C or g N / sq. meter )

     // Monthly drainage (mm / month)
     double drainage;

    // Monthly soil temperature ( degrees C ) at 10 cm
    double dst10;

     // Monthly estimated actual Evapotranspiration (mm / month)
     double eet;

     // Maximum EET of current year (mm / month)
     double eetmx;

     // Erosion of particulate organic matter (POM)
     Biomass erodePOM;  // (g C or g N / (sq. meter * month)

     // Evaporation (mm)
	 double evaporation;

     // Volume of water at field capacity (mm)
     double fldcap;         

     //  Formation of dissolved bicarbonate (HCO3) in the soil 
     //    column
     double formHCO3; // (g C / (sq. meter * month))

     //  Formation of dissolved bicarbonate (CO3) from rocks in 
     //    soil column
     double formRHCO3;     // (g C / (sq. meter * month))

     // Gasous carbon dioxide in soil column
     double gasCO2;  // g C / sq. meter

      // Water yield (mm / month)
     double h2oyld;  

     // Dissolved bicarbonate
     double HCO3;    // g C / sq. meter

     // Initial Estimated Actual Evapotranspiration (mm / month)
     double ineet;

     // Relative hydraulic conductivity through soil profile
     double kh2o;
     
     // Leaching of alkalinity (ALK) from soil column
     double leachALK;   // (g C / (sq. meter * month))

     // Leaching of dissolved carbon dioxide (CO2) from soil column
     double leachCO2;   // (g C / (sq. meter * month))

     // Leaching of dissolved organic matter (DOM) from
     //   soil column
     Biomass leachDOM;  // (g C or g N / (sq. meter * month))


     // Leaching of dissolved bicarbonate (HCO3) from soil column
     double leachHCO3;  // (g C / (sq. meter * month))

     // Leaching of nitrate (NO3) from soil column
     double leachNO3;   // (g N / (sq. meter * month))

     // Mean annual volumetric soil moisture (%)
     double meanvsm;

     // Soil moisture (mm)
     double moist;   

     // N2 flux (g N / (sq. meter * month))
     double n2flux;

     // N2O flux
     double n2oflux;   // (g N / (sq. meter * month))

     double ndays[CYCLE];

     // Next month's soil temperature at 10 cm
     double nextdst10; 

     // NH3 flux
     double nh3flux;   // (g N / (sq. meter * month))

     // Total nitrogen input to soils
     double ninput;    // (g N / (sq. meter * month)) 

     // Total nitrogen lost from soils
     double nlost;     // (g N / (sq. meter * month))

     // NO flux
     double noflux;    // (g N / (sq. meter * month))

     // Nonreactive soil organic matter
     Biomass nonReactiveOrg; //  (g C or g N / sq. meter)

      // Reactive soil organic matter
     Biomass org;      //  (g C or g N / sq. meter)

     // Soil moisture as %field capacity
     double pcfc;    

     // Percent clay in soil texture
     double pctclay;        

     // Soil moisture as %total pore space
     double pctp;    

     // Percent sand in soil texture
     double pctsand;        
     
     // Percent silt in soil texture
     double pctsilt;

     // Previous month's soil temperature at 10 cm
     double prevdst10; 

     // Previous month's snow pack
     double prevspack;
     
     // Maximum EET of previous year (mm / month)
     double prveetmx;
             
     // Proportion silt and clay in soil texture
     double psiplusc;       

     // Rain runoff storage (mm / month)
     double rgrndh2o;

     // Dissolved bicarbonate from rocks
     double RHCO3;

     // Rain percolation (excess)
     double rperc;   // (mm / month)

     // Rain Runoff (mm / month)
     double rrun;

     // Snowmelt runoff storage (mm)
     double sgrndh2o;

     // Snow melt infiltration (mm / month)
     double snowinf; 

     // Snowpack (mm)
     double snowpack;

     // Snow melt percolation (excess)
     double sperc;   // (mm / month)

     // Snow runoff
     double srun;    // (mm / month)

     // Surface runoff (mm / month)
     double surfrun;

     // Soil texture (categorical data)
     int text;              

      // Total soil organic matter
     Biomass totalOrg;   // (g C or g N / sq. meter)

     // volume of total pore space (mm)
     double totpor;         

     // Mean soil temperature (degrees C )of top 20 cm of soil 
     //   profile
     double tsoil;
     
     // Volumetric soil moisture (as %rooting depth)
     double vsm;     
 
     // Water-filled pore space offset in wetlands
	 double wfpsoff;

     // Volume of water at wilting point (mm)
     double wiltpt;         

     // wetland soil type designation (categorical data)
     int wsoil;


/* *************************************************************
		 Private Parameters
************************************************************* */

     // DOM

     double DOCpar[MAXCMNT];

     double DONpar;
     double DONparcut[MAXCMNT];
     double DONpar1a[MAXCMNT];
     double DONpar1b[MAXCMNT];
     double DONpar2a[MAXCMNT];
     double DONpar2b[MAXCMNT];

     // Field capacity (%soil volume)

     double fldcapa;
     double fldcapb;
     double pcfldcap;

     // Leaching parameters

     double lchDOMpar[MAXCMNT];

     double lchNO3par;
     double lchNO3parcut[MAXCMNT];
     double lchNO3par1a[MAXCMNT];
     double lchNO3par1b[MAXCMNT];
     double lchNO3par2a[MAXCMNT];
     double lchNO3par2b[MAXCMNT];

     double maxdocCN[MAXCMNT];
     double mindocCN[MAXCMNT];


     //Proportion of total soil organic matter that is 
     //   nonreactive

     double nonReactiveOMpar[MAXCMNT];


     // Porosity of soil (%soil volume)

     double pctpor;
     double pctpora;
     double pctporb;

     //Proportion of reactive soil organic matter that is 
     //   unfrozen within rooting zone

     double propReactA[MAXCMNT];
     double propReactB[MAXCMNT];


     // Effective rooting depth (m)

     double minrootz[MAXCMNT];
     double rootz;
     double rootza[MAXCMNT];
     double rootzb[MAXCMNT];
     double rootzc[MAXCMNT];


     // Wilting point (%soil volume)

     double pcwiltpt;
     double wiltpta;
     double wiltptb;
     
};

#endif
