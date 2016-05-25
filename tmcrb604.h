/* **************************************************************
*****************************************************************
TMCRB604.H - object describing characteristics of soil microbes
	     used in the Terrestrial Ecosystem Model (TEM)

Modifications:

20030718 - DWK created by modifying tmcrb50b1.h
20030718 - DWK changed class Tmicrobe50 to class Tmicrobe51
20030718 - DWK added public functions gminFastxClm(),
           gminSlowxClm(), nimmFastxClm() and
           setTraceGasProduction()
20030718 - DWK renamed nuptake[] to be immob[] and yrnuptake
           to be yrimmb
20030718 - DWK added public double grossnmin[], double nfix[],
           double ammnvol[], double nitrif[], double noprod[],
           double noflux[], double no3prod[], double denitrif[],
           double n2oprod[], double n2oflux[], double n2prod[],
           double n2flux[], double yrgmin, double yrnfix,
           double yrammnvol, double yrnitrif, double yrnoprd,
           double yrno3prd, double yrdenitrif,
           double yrn2oprd,  double yrn2prd,

20030718 - DWK added parameters kdcut[], kd1a[], kd1b[], kd2a[],
           and kd2b[]; and deleted parameters kda[] and kdb[]
20030718 - DWK added public double nfixpara[MAXCMNT],
           double nfixparb[MAXCMNT], double ntrfpar[MAXCMNT],
           double ammnpar[MAXCMNT], and double ntrfpar[MAXCMNT]
20030718 - DWK changed parameters "nup", "nupa", "nupb" to
           "nimm", "nimma" and "nimmb", respectively
20030718 - DWK changed include from tmcrb50b1.cpp to tmcrb51.cpp
           at bottom of file
20030722 - DWK added public function setNitrification()
20030722 - DWK changed public double ntrfpar[MAXCMNT] to double
           maxntrfCN[MAXCMNT] and added double minntrfCN[MAXCMNT]
20030825 - DWK added public double relntrf
20030826 - DWK added public Biomass DOMprod[CYCLE] and Biomass
           yrDOMprd
20030902 - DWK added public parameters double no3imm, double
           no3imma[MAXCMNT] and double no3immb[MAXCMNT]; and
           changed public parameter double nimm to double nh4imm,
           double nimma[MAXCMNT] to nh4imma[MAXCMNT] and
           double nimmb[MAXCMNT] to nh4immb[MAXCMNT]
20030902 - DWK modified function call to public function
           nimmFastxClm()
20031016 - DWK replaced char ecd[MAXFNAME] with string ecd in
           function call to getvegecd()
20031021 - DWK changed include from tprocessXML431a.h
           to tprocessXML51.h
20031021 - DWK changed inheritance from ProcessXML to ProcessXML51
20031204 - DWK added public parameter double tgmpar[MAXCMNT]
20031204 - DWK added const int& dcmnt to function call of 
           setTraceGasProduction()
20040203 - DWK rename setNitrification() as setNitrifStatus()
           and changed the determination of pntrf to be based on
           monthly net N mineralization instead of soil C:N 
20040203 - DWK added public function nitrfxclm()
20040203 - DWK renamed double maxntrfCN[MAXCMNT] to be 
           double initntrf[MAXCMNT]
20040203 - DWK renamed double minntrfCN[MAXCMNT] to be
           double allntrf[MAXCMNT]
20040203 - DWK changed include from tmcrb51.cpp to tmcrb512.cpp
           at bottom of file
20040212 - DWK added public parameter double ntrfpar[MAXCMNT]
20040228 - DWK changed class Tmicrobe51 to class Tmicrobe60
20040228 - DWK changed double nh4imma[MAXCMNT] to 
           nh4imm1a[MAXCMNT] 
20040228 - DWK changed double nh4immb[MAXCMNT] to 
           nh4imm1b[MAXCMNT]
20040228 - DWK added double nh4immcut[MAXCMNT], 
           double nh4imm2a[MAXCMNT] and double nh4imm2b[MAXCMNT]
20040228 - DWK deleted double immno3[CYCLE], double yrimmno3, 
           double no3imm, double no3imma[MAXCMNT] 
           and double no3immb[MAXCMNT]
20040228 - DWK changed include from tmcrb51.cpp to tmcrb60.cpp
           at bottom of file           
20040707 - DWK changed include from tprocessXML51.h to
           tprocessXML60.h
20040707 - DWK changed inheritance of ProcessXML51 to 
           inheritance of ProcessXML60
20040919 - DWK changed public double ntrfpar[MAXCMNT] to 
           double ntrfpar and added double ntrfparcut[MAXCMNT],
           double ntrfpar1a[MAXCMNT], double ntrfpar1b[MAXCMNT],
           double ntrfpar2a[MAXCMNT], and 
           double ntrfpar2b[MAXCMNT]
20050409 - DWK added public double decomp[CYCLE] and 
           double ndecomp[CYCLE]
20050409 - DWK changed include from tmcrb60.cpp to tmcrb601.cpp
           at bottom of file
20051117 - DWK added include temconsts602.hpp
20051117 - DWK changed include from tprocessXML60.h to
           tprocessXML602.h
20051117 - DWK deleted tmcrb601.cpp from bottom of file
20051125 - DWK changed public double ammnvol[CYCLE] to 
           double ammnvol
20051125 - DWK changed public double decomp[CYCLE] to 
           double decomp
20051125 - DWK changed public double denitrif[CYCLE] to 
           double denitrif                                                       
20051125 - DWK changed public Biomass DOMprod[CYCLE] to 
           Biomass DOMprod
20051125 - DWK changed public double grossnmin[CYCLE] to 
           double grossnmin
20051125 - DWK changed public double immnh4[CYCLE] to 
           double immnh4
20051125 - DWK changed public double immob[CYCLE] to double immob 
20051125 - DWK changed public double n2oprod[CYCLE] to 
           double n2oprod
20051125 - DWK changed public double n2prod[CYCLE] to 
           double n2prod
20051125 - DWK changed public double ndecomp[CYCLE] to 
           double ndecomp
20051125 - DWK changed public double netnmin[CYCLE] to 
           double netnmin
20051125 - DWK changed public double nfix[CYCLE] to double nfix                                                                            
20051125 - DWK changed public double nitrif[CYCLE] to 
           double nitrif
20051125 - DWK changed public double noprod[CYCLE] to 
           double noprod
20051125 - DWK changed public double no3prod[CYCLE] to 
           double no3prod                      
20051125 - DWK changed public double rh[CYCLE] to double rh           
20051125 - DWK deleted const int& outmon from function call to
           setNitrifStatus() and setTraceGasProduction()
20051130 - DWK added public function updateDynamics()
20051130 - DWK added public parameters double DOCpar[MAXCMNT],
           double DONpar, double DONparcut[MAXCMNT],
           double DONpar1a[MAXCMNT], double DONpar1b[MAXCMNT],
           double DONpar2a[MAXCMNT], and 
           double DONpar2b[MAXCMNT]
20051201 - DWK added public function setRHMOIST()
20051202 - DWK added public function setDQ10()
20051202 - DWK added private double dq10
20051202 - DWK added public function resetMonthlyFluxes()
20051202 - DWK added public function resetYrFluxes()
20051207 - DWK added public function resetEcds()
           
*****************************************************************
************************************************************** */

// Tmicrobe51 uses the global constants CYCLE, NUMMSAC and NUMVEG

#ifndef TMCRB604_H
#define TMCRB604_H

#include "temconsts602.hpp"
#include "tprocessXML602.h"
#include "bioms423.hpp"

class Tmicrobe60: public ProcessXML60
{

   public:

     Tmicrobe60();

/* **************************************************************
		 Public Functions
************************************************************** */

     void getvegecd( ofstream& rflog1 );

     void getvegecd( const string& ecd );

     void resetEcds( const int& pcmnt, const double& psiplusc );
      
     void resetMonthlyFluxes( void );

     void resetYrFluxes( void );

     void setDQ10( const int& pdcmnt, 
                   const double& tair, 
                   const double& tsoil,
                   const double& snowpack, 
                   const int& tsoilflg );

     void setNewDQ10( const int& pdcmnt, 
                      const double& tair, 
                      const double& tsoil,
                      const double& snowpack, 
                      const int& tsoilflg );

     inline double getDQ10(void)
     {
       return dq10;
     };

     double setNitrifStatus( const int& dcmnt );

     double setRHMOIST( const int& pdcmnt,
                        const double& pcfldcap,
                        const double& vsm,
                        const int& moistlim );

     void setTraceGasProduction( const int& pdcmnt,
                                 const double& pctp,
                                 double& agfertn );

     void showecd( const int& pdcmnt );

     void updateDynamics( const int& pcmnt,
                          const double& pcfldcap,
                          const double& actlayer,
                          const double& soilorgc,
                          const double& agr,
                          const double& agl,
                          const double& bgr,
                          const double& bgl,
                          const double& cwd,
                          const double& soilorgn,
                          const double& soilh2o,
                          const double& vsm,
                          const double& nh4,
                          const double& nh4dep,
                          const int& moistlim,
                          const int& tillflag,
                          const double& tillfactor,
                          const double& ksoil );
     
     double yrkd( const int& nfeed,
                  const double& yrltrc,
                  const double& yrltrn,
                  const double& lcclnc );

     // "Get" and "Set" private variables and parameters
     
     // alpha *****************************************************

     inline double getALPHA( const int& pcmnt )
     {
       return alpha[pcmnt];
     };

     inline void setALPHA( const double& palpha, const int& pcmnt )
     {
       alpha[pcmnt] = palpha;
     };


     // allntrf ************************************************
     
     inline double getALLNTRF( const int& pcmnt ) 
     { 
       return allntrf[pcmnt]; 
     }

     inline void setALLNTRF( const double& pallntrf, 
                             const int& pcmnt ) 
     { 
       allntrf[pcmnt] = pallntrf; 
     }


     // ammnpar ************************************************
     
     inline double getAMMNPAR( const int& pcmnt ) 
     { 
       return ammnpar[pcmnt]; 
     }

     inline void setAMMNPAR( const double& pammnpar, 
                             const int& pcmnt ) 
     { 
       ammnpar[pcmnt] = pammnpar; 
     }

     
     // ammnvol ************************************************
     
     inline double getAMMNVOL( void ) { return ammnvol; }

     inline void setAMMNVOL( const double& pammnvol ) 
     { 
       ammnvol = pammnvol; 
     }


     // beta *****************************************************

     inline double getBETA( const int& pcmnt )
     {
       return beta[pcmnt];
     };

     inline void setBETA( const double& pbeta, const int& pcmnt )
     {
       beta[pcmnt] = pbeta;
     };


     // cnsoil *************************************************
     
     inline double getCNSOIL( const int& pcmnt ) 
     { 
       return cnsoil[pcmnt]; 
     }

     inline void setCNSOIL( const double& pcnsoil, 
                            const int& pcmnt ) 
     { 
       cnsoil[pcmnt] = pcnsoil; 
     }


     // decay **************************************************
     
     inline double getDECAY( void ) { return decay; }


     // decomp *************************************************
     
     inline double getDECOMP( void ) { return decomp; }


     // denitrif ***********************************************
     
     inline double getDENITRIF( void ) { return denitrif; }

     inline void setDENITRIF( const double& pdenitrif ) 
     { 
       denitrif = pdenitrif; 
     }


     // DOCpar *************************************************
     
     inline double getDOCPAR( const int& pcmnt ) 
     { 
       return DOCpar[pcmnt]; 
     }

     inline void setDOCPAR( const double& pdocpar, 
                            const int& pcmnt ) 
     { 
       DOCpar[pcmnt] = pdocpar; 
     }


     // DOMprod.carbon *****************************************
     
     inline double getDOCPROD( void ) { return DOMprod.carbon; }


     // DONpar *************************************************
     
     inline double getDONPAR( void ) { return DONpar; }
     
     inline void setDONPAR( const double& pdonpar ) 
     { 
       DONpar = pdonpar; 
     }
     
     
     // DONparcut **********************************************
     
     inline double getDONPARCUT( const int& pcmnt ) 
     { 
       return DONparcut[pcmnt]; 
     }

     inline void setDONPARCUT( const double& pdonparcut, 
                               const int& pcmnt ) 
     { 
       DONparcut[pcmnt] = pdonparcut; 
     }


     // DONpar1a ***********************************************
     
     inline double getDONPAR1A( const int& pcmnt ) 
     { 
       return DONpar1a[pcmnt]; 
     }

     inline void setDONPAR1A( const double& pdonpar1a, 
                              const int& pcmnt ) 
     { 
       DONpar1a[pcmnt] = pdonpar1a; 
     }


     // DONpar1b ***********************************************
     
     inline double getDONPAR1B( const int& pcmnt ) 
     { 
       return DONpar1b[pcmnt]; 
     }

     inline void setDONPAR1B( const double& pdonpar1b, 
                              const int& pcmnt ) 
     { 
       DONpar1b[pcmnt] = pdonpar1b; 
     }


     // DONpar2a ***********************************************
     
     inline double getDONPAR2A( const int& pcmnt ) 
     { 
       return DONpar2a[pcmnt]; 
     }

     inline void setDONPAR2A( const double& pdonpar2a, 
                              const int& pcmnt ) 
     { 
       DONpar2a[pcmnt] = pdonpar2a; 
     }


     // DONpar2b ***********************************************
     
     inline double getDONPAR2B( const int& pcmnt ) 
     { 
       return DONpar2b[pcmnt]; 
     }

     inline void setDONPAR2B( const double& pdonpar2b, 
                              const int& pcmnt ) 
     { 
       DONpar2b[pcmnt] = pdonpar2b; 
     }


     // DOMprod.nitrogen ***************************************
     
     inline double getDONPROD( void ) 
     { 
       return DOMprod.nitrogen; 
     }


     // gamma *****************************************************

     inline double getGAMMA( const int& pcmnt )
     {
       return gamma[pcmnt];
     };

     inline void setGAMMA( const double& pgamma, const int& pcmnt )
     {
       gamma[pcmnt] = pgamma;
     };


     // grossnmin **********************************************
     
     inline double getGROSSNMIN( void ) { return grossnmin; }


     // immnh4 *************************************************
     
     inline double getIMMNH4( void ) { return immnh4; }


     // immob **************************************************
     
     inline double getIMMOB( void ) { return immob; }


     // initntrf ***********************************************
     
     inline double getINITNTRF( const int& pcmnt ) 
     { 
       return initntrf[pcmnt]; 
     }

     inline void setINITNTRF( const double& pinntrf, 
                              const int& pcmnt ) 
     { 
       initntrf[pcmnt] = pinntrf; 
     }


     // kd *****************************************************
     
     inline double getKD( void ) { return kd; }

     inline void setKD( const double& pkd ) { kd = pkd; }


     // kdcut **************************************************
     
     inline double getKDCUT( const int& pcmnt ) 
     { 
       return kdcut[pcmnt]; 
     }

     inline void setKDCUT( const double& pkdcut, 
                           const int& pcmnt ) 
     { 
       kdcut[pcmnt] = pkdcut; 
     }


     // kd1a ***************************************************
     
     inline double getKD1A( const int& pcmnt ) 
     { 
       return kd1a[pcmnt]; 
     }

     inline void setKD1A( const double& pkd1a, 
                          const int& pcmnt ) 
     { 
       kd1a[pcmnt] = pkd1a; 
     }


     // kd1b ***************************************************
     
     inline double getKD1B( const int& pcmnt ) 
     { 
       return kd1b[pcmnt]; 
     }

     inline void setKD1B( const double& pkd1b, 
                          const int& pcmnt ) 
     { 
       kd1b[pcmnt] = pkd1b; 
     }


     // kd2a ***************************************************
     
     inline double getKD2A( const int& pcmnt ) 
     { 
       return kd2a[pcmnt]; 
     }

     inline void setKD2A( const double& pkd2a, 
                          const int& pcmnt ) 
     { 
       kd2a[pcmnt] = pkd2a; 
     }


     // kd2b ***************************************************
     
     inline double getKD2B( const int& pcmnt ) 
     { 
       return kd2b[pcmnt]; 
     }

     inline void setKD2B( const double& pkd2b, 
                          const int& pcmnt ) 
     { 
       kd2b[pcmnt] = pkd2b; 
     }


     // kdc ****************************************************
     
     inline double getKDC( void ) { return kdc; }

     inline void setKDC( const double& pkdc ) { kdc = pkdc; }


     // kn2 ****************************************************
     
     inline double getKN2( const int& pcmnt ) 
     { 
       return kn2[pcmnt]; 
     }

     inline void setKN2( const double& pkn2, 
                         const int& pcmnt ) 
     { 
       kn2[pcmnt] = pkn2; 
     }


//     // lcclnc *************************************************
     
//     inline double getLCCLNC( const int& pcmnt ) 
//     { 
//       return lcclnc[pcmnt]; 
//     }

//     inline void setLCCLNC( const double& plcclnc, 
//                            const int& pcmnt ) 
//     { 
//       lcclnc[pcmnt] = plcclnc; 
//     }


     // moistmax ***********************************************
     
     inline double getMOISTMAX( const int& pcmnt ) 
     { 
       return moistmax[pcmnt]; 
     }

     inline void setMOISTMAX( const double& pmoistmax, 
                              const int& pcmnt ) 
     { 
       moistmax[pcmnt] = pmoistmax; 
     }


     // moistmin ***********************************************
     
     inline double getMOISTMIN( const int& pcmnt ) 
     { 
       return moistmin[pcmnt]; 
     }

     inline void setMOISTMIN( const double& pmoistmin, 
                              const int& pcmnt ) 
     { 
       moistmin[pcmnt] = pmoistmin; 
     }


     // moistopt ***********************************************
     
     inline double getMOISTOPT( const int& pcmnt ) 
     { 
       return moistopt[pcmnt]; 
     }

     inline void setMOISTOPT( const double& pmoistopt, 
                              const int& pcmnt ) 
     { 
       moistopt[pcmnt] = pmoistopt; 
     }


     // n2oprod ************************************************
     
     inline double getN2OPROD( void ) { return n2oprod; }

     inline void setN2OPROD( const double& pn2oprod ) 
     { 
       n2oprod = pn2oprod; 
     }

     // n2prod *************************************************
     
     inline double getN2PROD( void ) { return n2prod; }

     inline void setN2PROD( const double& pn2prod ) 
     { 
       n2prod = pn2prod; 
     }

     // ndecomp ************************************************
     
     inline double getNDECOMP( void ) { return ndecomp; }


     // netnmin ************************************************
     
     inline double getNETNMIN( void ) { return netnmin; }

     inline void setNETNMIN( const double& pnetnmin ) 
     { 
       netnmin = pnetnmin; 
     }


     // nfix ***************************************************
     
     inline double getNFIX( void ) { return nfix; }

     inline void setNFIX( const double& pnfix ) 
     { 
       nfix = pnfix; 
     }


     // nfixpar ************************************************
     
     inline double getNFIXPAR( const int& pcmnt ) 
     { 
       return nfixpar[pcmnt]; 
     }

     inline void setNFIXPAR( const double& pnfixpar, 
                             const int& pcmnt ) 
     { 
       nfixpar[pcmnt] = pnfixpar; 
     }


     // nh4imm *************************************************
     
     inline double getNH4IMM( void ) { return nh4imm; }

     inline void setNH4IMM( const double& pnh4imm ) 
     { 
       nh4imm = pnh4imm; 
     }


     // nh4immcut **********************************************
     
     inline double getNH4IMMCUT( const int& pcmnt ) 
     { 
       return nh4immcut[pcmnt]; 
     }

     inline void setNH4IMMCUT( const double& pnh4immcut, 
                               const int& pcmnt ) 
     { 
       nh4immcut[pcmnt] = pnh4immcut; 
     }


     // nh4imm1a ***********************************************
     
     inline double getNH4IMM1A( const int& pcmnt ) 
     { 
       return nh4imm1a[pcmnt]; 
     }

     inline void setNH4IMM1A( const double& pnh4imm1a, 
                              const int& pcmnt ) 
     { 
       nh4imm1a[pcmnt] = pnh4imm1a; 
     }


     // nh4imm1b ***********************************************
     
     inline double getNH4IMM1B( const int& pcmnt ) 
     { 
       return nh4imm1b[pcmnt]; 
     }

     inline void setNH4IMM1B( const double& pnh4imm1b, 
                              const int& pcmnt ) 
     { 
       nh4imm1b[pcmnt] = pnh4imm1b; 
     }


     // nh4imm2a ***********************************************
     
     inline double getNH4IMM2A( const int& pcmnt ) 
     { 
       return nh4imm2a[pcmnt]; 
     }

     inline void setNH4IMM2A( const double& pnh4imm2a, 
                              const int& pcmnt ) 
     { 
       nh4imm2a[pcmnt] = pnh4imm2a; 
     }


     // nh4imm2b ***********************************************
     
     inline double getNH4IMM2B( const int& pcmnt ) 
     { 
       return nh4imm2b[pcmnt]; 
     }

     inline void setNH4IMM2B( const double& pnh4imm2b, 
                              const int& pcmnt ) 
     { 
       nh4imm2b[pcmnt] = pnh4imm2b; 
     }


     // nimm ***************************************************
     
     inline double getNIMM( void ) { return nimm; }


     // nitrif *************************************************
     
     inline double getNITRIF( void ) { return nitrif; }


     // noprod *************************************************
     
     inline double getNOPROD( void ) { return noprod; }


     // no3prod ************************************************
     
     inline double getNO3PROD( void ) { return no3prod; }


     // ntrfpar ************************************************
     
     inline double getNTRFPAR( void ) { return ntrfpar; }

     inline void setNTRFPAR( const double& pntrfpar ) 
     { 
       ntrfpar = pntrfpar; 
     }


     // ntrfparcut *********************************************
     
     inline double getNTRFPARCUT( const int& pcmnt ) 
     { 
       return ntrfparcut[pcmnt]; 
     }
    
     inline void setNTRFPARCUT( const double& pntrfparcut, 
                                const int& pcmnt ) 
     { 
       ntrfparcut[pcmnt] = pntrfparcut; 
     }


     // ntrfpar1a **********************************************
     
     inline double getNTRFPAR1A( const int& pcmnt ) 
     { 
       return ntrfpar1a[pcmnt]; 
     }

     inline void setNTRFPAR1A( const double& pntrfpar1a, 
                               const int& pcmnt ) 
     { 
       ntrfpar1a[pcmnt] = pntrfpar1a; 
     }


     // ntrfpar1b **********************************************
     
     inline double getNTRFPAR1B( const int& pcmnt ) 
     { 
       return ntrfpar1b[pcmnt]; 
     }

     inline void setNTRFPAR1B( const double& pntrfpar1b, 
                               const int& pcmnt ) 
     { 
       ntrfpar1b[pcmnt] = pntrfpar1b; 
     }


     // ntrfpar2a **********************************************
     
     inline double getNTRFPAR2A( const int& pcmnt ) 
     { 
       return ntrfpar2a[pcmnt]; 
     }
     
     inline void setNTRFPAR2A( const double& pntrfpar2a, 
                               const int& pcmnt ) 
     { 
       ntrfpar2a[pcmnt] = pntrfpar2a; 
     }

     // ntrfpar2b **********************************************

     inline double getNTRFPAR2B( const int& pcmnt ) 
     { 
       return ntrfpar2b[pcmnt]; 
     }

     inline void setNTRFPAR2B( const double& pntrfpar2b, 
                               const int& pcmnt ) 
     { 
       ntrfpar2b[pcmnt] = pntrfpar2b; 
     }


     // propftos ***********************************************
     
     inline double getPROPFTOS( const int& pcmnt ) 
     { 
       return propftos[pcmnt]; 
     }

     inline void setPROPFTOS( const double& ppropftos, 
                              const int& pcmnt ) 
     { 
       propftos[pcmnt] = ppropftos; 
     }


     // qref *****************************************************

     inline double getQREF( const int& pcmnt )
     {
       return qref[pcmnt];
     };

     inline void setQREF( const double& pqref, const int& pcmnt )
     {
       qref[pcmnt] = pqref;
     };


     // rh *****************************************************
     
     inline double getRH( void ) { return rh; }

     inline double getRH_AGR( void ) { return rh_agr; }
     inline double getRH_AGL( void ) { return rh_agl; }
     inline double getRH_BGR( void ) { return rh_bgr; }
     inline double getRH_BGL( void ) { return rh_bgl; }
     inline double getRH_CWD( void ) { return rh_cwd; }
     inline double getRH_SOIL( void ) { return rh_soil; }
     inline double getSOIL_LITTER( void ) { return tosoil_litter; }
     // rhq10 **************************************************
     
//     inline double getRHQ10( const int& pcmnt ) 
//     { 
//       return rhq10[pcmnt]; 
//     }

//     inline void setRHQ10( const double& prhq10, 
//                           const int& pcmnt ) 
//     { 
//       rhq10[pcmnt] = prhq10; 
//     }


     // tgmpar *************************************************
     
     inline double getTGMPAR( const int& pcmnt ) 
     { 
       return tgmpar[pcmnt]; 
     }

     inline void setTGMPAR( const double& ptgmpar, 
                            const int& pcmnt ) 
     { 
       tgmpar[pcmnt] = ptgmpar; 
     }
 
     
     // tref *****************************************************

     inline double getTREF( const int& pcmnt )
     {
       return tref[pcmnt];
     };

     inline void setTREF( const double& ptref, const int& pcmnt )
     {
       tref[pcmnt] = ptref;
     };

     inline double getF0( const int& pcmnt )
     {
       return factor0[pcmnt];
     };

     inline void setF0( const double& factor, const int& pcmnt )
     {
       factor0[pcmnt] = factor;
     };
     inline double getF1( const int& pcmnt )
     {
       return factor1[pcmnt];
     };

     inline void setF1( const double& factor, const int& pcmnt )
     {
       factor1[pcmnt] = factor;
     };
     inline double getF2( const int& pcmnt )
     {
       return factor2[pcmnt];
     };

     inline void setF2( const double& factor, const int& pcmnt )
     {
       factor2[pcmnt] = factor;
     };
     inline double getF3( const int& pcmnt )
     {
       return factor3[pcmnt];
     };

     inline void setF3( const double& factor, const int& pcmnt )
     {
       factor3[pcmnt] = factor;
     };
     inline double getF4( const int& pcmnt )
     {
       return factor4[pcmnt];
     };

     inline void setF4( const double& factor, const int& pcmnt )
     {
       factor4[pcmnt] = factor;
     };

     inline double getF5( const int& pcmnt )
     {
       return factor5[pcmnt];
     };

     inline void setF5( const double& factor, const int& pcmnt )
     {
       factor5[pcmnt] = factor;
     };
     //inline double getrhmoist( void ) { return temp_rhmoist; }



/* **************************************************************
		 Public Variables
************************************************************** */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif

     double relntrf;         // nitrif / netnmin

     // Annual sum of ammnvol
     double yrammnvol;       // (g N / (sq. meter * year))

     // Annual sum of decomp
     double yrdecomp;     // (g C / (sq. meter * year))

     // Annual sum of denitrif
     double yrdenitrf;       // (g N / (sq. meter * year))

     // Annual sum of DOM
     Biomass yrDOMprod;     // (g C or g N / (sq. meter * year))

     // Annual sum of grossnmin
     double yrgmin;            // (g N / (sq. meter * year))

     // Annual sum of immnh4
     double yrimmnh4;        // (g N / (sq. meter * year))

     // Annual sum of immob
     double yrimmb;        // (g N / (sq. meter * year))

     // Annual sum of n2oprod
     double yrn2oprd;        // (g N / (sq. meter * year))

     // Annual sum of n2prod
     double yrn2prd;        // (g N / (sq. meter * year))

     // Annual sum of ndecomp
     double yrndecomp;      // (g N / (sq. meter * year))

     // Annual sum of netnmin
     double yrnmin;         // (g N / (sq. meter * year))

     // Annual sum of nfix
     double yrnfix;         // (g N / (sq. meter * year))

     // Annual sum of nitrif
     double yrnitrif;       // (g N / (sq. meter * year))

     // Annual sum of noprod
     double yrnoprd;         // (g N / (sq. meter * year))

     // Annual sum of no3prod
     double yrno3prd;        // (g N / (sq. meter * year))

     // Annual sum of rh
     double yrrh;       // (g C / (sq. meter * year))


   private:

/* **************************************************************
		 Private Functions
************************************************************** */

     double gminFastxClm( const double& soilorgc,
                          const double& soilorgn,
                          const double& rh );

     double gminSlowxClm( const int& pdcmnt,
                          const double& soilorgc,
                          const double& soilorgn,
                          const double& rh );

     double nimmFastxClm( const int& pdcmnt,
                          const double& nimm,
                          const double& soilh2o,
                          const double& availn,
                          const double& decay,
                          const double& rh,
                          const double& ksoil );

     double nminxclm( const int& pdcmnt,
                      const double& soilh2o,
                      const double& soilorgc,
                      const double& soilorgn,
                      const double& availn,
                      const double& decay,
                      const double& rh,
                      const double& ksoil );

     double ntrfxclm( const int& pdcmnt,
                      const double& soilh2o,
                      const double& nh4,
                      const double& nmin,
                      const double& ksoil );

     double rhxclm( const double& organicc,
                    const double& dq10,
                    const double& moist,
                    const double& factor );

   
/* **************************************************************
		 Private Variables
************************************************************** */
     
      // Ammonia Volatilization (g N / (sq. meter * month))
      double ammnvol; 

      // Decomposition (g C / (sq. meter * month))
      double decomp;
      
      // Denitrification (g N / (sq. meter * month))
      double denitrif; 

      // Dissolved Organic Matter (DOM) production
      Biomass DOMprod; // (g C or g N / (sq. meter * month))

     // Effect of temperature on decomposition
     double dq10;

     // Gross nitrogen mineralization
     double grossnmin;  // (g N / (sq. meter * month))

     // Ammonium uptake or "immobilzation" by microbes
     double immnh4;   // (g N / (sq. meter * month))

     // Total nitrogen uptake or "immobilzation" by microbes
     double immob;  // (g N / (sq. meter * month))

     // N2O production
     double n2oprod;  // (g N / (sq. meter * month))

     // N2 production
     double n2prod;  // (g N / (sq. meter * month))

     // Nitrogen in Decomposition      
     double ndecomp; // (g N / (sq. meter * month))

     // Net nitrogen mineralization
     double netnmin; // (g N / (sq. meter * month))

     // Nitrogen fixation by free-living microbes
     double nfix;    // (g N / (sq. meter * month))

     // Nitrification
     double nitrif;  // (g N / (sq. meter * month))

     // NO production
     double noprod;   // (g N / (sq. meter * month))

      // NO3 production
      double no3prod;  // (g N / (sq. meter * month))

     // Heterotrophic respiration
     double rh;  // (g C / (sq. meter * month))


/* *************************************************************
		 Private Parameters
************************************************************* */

     double alpha[MAXCMNT];

     // Ammonia volatilization

     double ammnpar[MAXCMNT];

     double beta[MAXCMNT];

     double cnsoil[MAXCMNT];

     // Parameter representing the quality of soil organic matter

     double decay;

     double DOCpar[MAXCMNT];

     double DONpar;
     double DONparcut[MAXCMNT];
     double DONpar1a[MAXCMNT];
     double DONpar1b[MAXCMNT];
     double DONpar2a[MAXCMNT];
     double DONpar2b[MAXCMNT];

     double gamma[MAXCMNT];

     // Biome-specific decomposition parameters for function rhxclm

     double kd;
     double kdcut[MAXCMNT];
     double kd1a[MAXCMNT];
     double kd1b[MAXCMNT];
     double kd2a[MAXCMNT];
     double kd2b[MAXCMNT];
     double kdc;

     double kdin[NUMMSAC];     // kd values read in from file
     double kdsave[NUMMSAC];   // kd values saved to a file

     // Biome-specific half saturation parameter for function
     //   nminxclm describing the effect of available nitrogen
     //   on microbial nitrogen uptake

     double kn2[MAXCMNT];


//     double lcclnc[MAXCMNT];

     // Biome-specific parameters describing the influence of 
     //   soil moisture on decomposition (i.e., moist)

     double moistmin[MAXCMNT];
     double moistopt[MAXCMNT];
     double moistmax[MAXCMNT];

     // Biome-specific nitrification parameters

     double ntrfpar;
     double ntrfparcut[MAXCMNT];
     double ntrfpar1a[MAXCMNT];
     double ntrfpar1b[MAXCMNT];
     double ntrfpar2a[MAXCMNT];
     double ntrfpar2b[MAXCMNT];
     
     double allntrf[MAXCMNT];
     double initntrf[MAXCMNT];


    // N fixation parameter

     double nfixpar[MAXCMNT];


     // Biome-specific microbial immobilization parameters for
     //   function nminxclm

     double nh4imm;
     double nh4immcut[MAXCMNT];
     double nh4imm1a[MAXCMNT];
     double nh4imm1b[MAXCMNT];
     double nh4imm2a[MAXCMNT];
     double nh4imm2b[MAXCMNT];

     double nimm;
     
     double propftos[MAXCMNT];

     double qref[MAXCMNT];

     double rhq10[MAXCMNT];

     double tgmpar[MAXCMNT];

     double tref[MAXCMNT];

     double factor0[MAXCMNT];//soil decomposition
     double factor1[MAXCMNT];//AGR decomposition
     double factor2[MAXCMNT];//AGL decomposition
     double factor3[MAXCMNT];//BGR decomposition
     double factor4[MAXCMNT];//BGL decomposition
     double factor5[MAXCMNT];//CWD decomposition
     double rh_agr, rh_agl, rh_bgr, rh_bgl, rh_cwd, rh_soil, tosoil_litter;
     //double temp_rhmoist;
};

#endif

