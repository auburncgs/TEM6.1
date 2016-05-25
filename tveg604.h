/* **************************************************************
*****************************************************************
TVEG604.H -  Vegetation characteristics used in the Terrestrial
             Ecosystem Model (TEM)

Modifications:

20030718 - DWK created by modifying tveg50b1.h to incorporate
           open-N cycle dynamics into TEM
20030718 - DWK added include inorgn51.hpp
20030718 - DWK changed Class Tveg50 to Class Tveg51
20030718 - DWK changed public function gppio3() to gppxio3()
20030718 - DWK added public functions nupnh4xclm() and
           nupno3xclm()
20030718 - DWK changed public double innuptake[CYCLE] to
           InorgN innuptake[CYCLE]
20030718 - DWK changed public double yrinnup to InorgN yrinnup
20030718 - DWK changed public double nuptake[CYCLE] to
           InorgN nuptake[CYCLE]
20030718 - DWK changed public double yrnup to InorgN yrnup
20030718 - DWK added double innfix, nfix[CYCLE], snfix[CYCLE],
           lnfix[CYCLE], yrinnfix, yrnfix, yrsnfix and yrlnfix
20030718 - DWK added double nfixpara[MAXCMNT], nfixparb[MAXCMNT],
           nupnh4, nupnh41a[MAXCMNT], nupnh41b[MAXCMNT],
           nupnh42a[MAXCMNT], nupnh42b[MAXCMNT], nupno3,
           nupno31a[MAXCMNT], nupno31b[MAXCMNT], nupno32a[MAXCMNT],
           nupno32b[MAXCMNT], kvnh4[MAXCMNT], kvno3[MAXCMNT],
           rroot[MAXCMNT]
20030718 - DWK changed variables double& prevdst5, double& dst5,
           and double& nextdst5 to double& prevdst10,
           double& dst10, and double& nextdst10, respectively,
           in function call to setThawPercent()
20030718 - DWK changed include at bottom of file from tveg50b1.cpp
           to tveg51.cpp
20031014 - BSF added const double& eetpet to function call of
           gppxo3()
20031015 - DWK added VolatileOrganicCarbon voc[CYCLE] and
           VolatileOrganicCarbon yrvoc
20031015 - DWK added double rootResp[CYCLE] and double yrrootResp
20031016 - DWK replaced char ecd[MAXFNAME] with string ecd in
           function calls to getecd() and getleafecd()
20031019 - DWK changed include from tprocessXML431a.h to
           tprocessXML51.h
20031019 - DWK changed inheritance from ProcessXML to 
           ProcessXML51
20040707 - DWK changed class Tveg51 to class Tveg60
20040707 - DWK changed include from tprocessXML51.h to
           tprocessXML60.h
20040707 - DWK changed inheritance from ProcessXML51 to 
           ProcessXML60
20040707 - DWK changed include from tveg51.cpp to tveg60.cpp at
           bottom of file
20040716 - DWK changed include from inorgn51.h to inorgn60.h
20040716 - DWK changed public InorgN inuptake[CYCLE] to         
           InorgN60 inuptake[CYCLE]
20040716 - DWK changed public InorgN yrinnup to InorgN60 yrinnup
20040716 - DWK changed public InorgN nuptake[CYCLE] to         
           InorgN60 nuptake[CYCLE]
20040716 - DWK changed public InorgN yrnup to InorgN60 yrnup
20040717 - DWK changed include from voc51.hpp to voc60.hpp
20040717 - DWJ changed public VolatileOrganicCarbon voc[CYCLE]
           to VolatileOrganicCarbon60 voc[CYCLE]
20040717 - DWJ changed public VolatileOrganicCarbon yrvoc
           to VolatileOrganicCarbon60 yrvoc
20040828 - DWK changed public int temveg to int potveg
20040828 - DWK added public int currentveg
20051117 - DWK changed include from tprocessXML.60.h to
           tprocessXML602.h
20051117 - DWK deleted include tveg60.cpp from bottom of file
20051117 - DWK added include temconsts602.hpp"    
20051124 - DWK changed public double abvgrndResp[CYCLE] to 
           double abvgrndResp
20051124 - DWK changed public double findozone[CYCLE] to
           double findozone
20051124 - DWK changed public double fozone[CYCLE] to
           double fozone
20051124 - DWK changed public double fpc[CYCLE] to 
           double fpc
20051124 - DWK changed public double gpp[CYCLE] to
           double gpp
20051124 - DWK changed public double gpr[CYCLE] to
           double gpr
20051124 - DWK changed public double ingpp[CYCLE] to
           double ingpp
20051124 - DWK changed public double innpp[CYCLE] to 
           double innpp
20051124 - DWK changed public InorgN60 inuptake[CYCLE] to 
           InorgN60 inuptake
20051124 - DWK changed public Biomass labile[CYCLE] to
           Biomass labile
20051124 - DWK changed public double lai[CYCLE] to
           double lai
20051124 - DWK changed public double leaf[CYCLE] to 
           double leaf
20051124 - DWK changed public double lnfix[CYCLE] to 
           double lnfix
20051124 - DWK changed public Biomass ltrfal[CYCLE] to
           Biomass ltrfal
20051124 - DWK changed public double luptake[CYCLE] to
           double luptake
20051124 - DWK changed public double nfix[CYCLE] to double nfix
20051125 - DWK changed public double nmobil[CYCLE] to
           double nmobil
20051125 - DWK changed public double npp[CYCLE] to double npp
20051125 - DWK changed public double nresorb[CYCLE] to
           double nresorb
20051125 - DWK changed public InorgN60 nuptake[CYCLE] to
           InorgN60 nuptake
20051125 - DWK changed public Biomass plant[CYCLE] to
           Biomass plant
20051125 - DWK changed public double rg[CYCLE] to double rg
20051125 - DWK changed public double rm[CYCLE] to double rm
20051125 - DWK changed public double rootResp[CYCLE] to 
           double rootResp
20051125 - DWK changed public double snfix[CYCLE] to 
           double snfix
20051125 - DWK changed public Biomass strctrl[CYCLE] to 
           Biomass strctrl
20051125 - DWK changed public double suptake[CYCLE] to 
           double suptake
20051125 - DWK changed public double thawpercent[CYCLE] to
           double thawpercent
20051125 - DWK changed public double unnormleaf[CYCLE] to 
           double unnormleaf
20051125 - DWK changed public VolatileOrganicCarbon60 voc[CYCLE] 
           to VolatileOrganicCarbon60 voc
20051201 - DWK added public function setGV()
20051201 - DWK added public function setTEMP()
20051201 - DWK added private double temp
20051202 - DWK added public function setRESPQ10()
20051202 - DWK added private double respq10
20051202 - DWK added public function resetMonthlyFluxes()
20051202 - DWK added public function resetYrFluxes()
20051207 - DWK added public function resetECD()
20060608 - DWK added double extraNH4 to function call to 
           updateDynamics()
20070210 - DWK added const int& perennial to updateDynamics()                                                                                                                                                                                                                                                                                          
20070829 - DWK added public function updateFoliage()
20070829 - DWK added private double transpiration
20070829 - DWK addd provate parameter proptrans

****************************************************************
************************************************************* */

#ifndef TVEG604_H
#define TVEG604_H

// Tveg60 also uses the global constants CYCLE and NUMVEG
#include "temconsts602.hpp"

// Tveg60 uses Biomass Class
#include "bioms423.hpp"   

// Tveg60 uses InorgN60 Class
#include "inorgn60.hpp"   

// Tveg60 uses VolatileOrganicCarbon60 class
#include "voc60.hpp"      

// Tveg60 inherits ProcessXML60 Class
#include "tprocessXML602.h" 

class Tveg60 : public ProcessXML60
{

  public:

     Tveg60();

/* **************************************************************
		 Public Functions
************************************************************** */

     // Limit parameter topt to a biome-specific range of 
     //   air temperatures
     
     void boundTOPT( const int& pcmnt );

     void getecd( ofstream& rflog1 );

     void getecd( const string& ecd );

     void getleafecd( ofstream& rflog1 );

     void getleafecd( const string& ecd );

     void leafinit( ofstream& rflog1 );

     void resetEcds( const int& pcmnt, const double& psiplusc );

     void resetMonthlyFluxes( void );

     void resetNEWTOPT( const int& pcmnt, 
                        const double& tair, 
                        const double& unnrmleaf );

     void resetYrFluxes( void );

     double setGV( const double& eet,
                   const double& pet,
                   const int& moistlim );
     double getGV(void)
     {
    	 return gv;
     }

     void setNewRESPQ10( const int& pdcmnt, const double& tair );

//     void setRESPQ10( const int& pdcmnt, const double& tair );
     double getNewRESPQ10( void)
     {
    	 return respq10;
     }

     void setTEMP( const int& pdcmnt, const double& tair );
     double getTEMP( void)
     {
    	 return temp;
     }


     void setThawPercent( const double& pprevdst10,
                          const double& pdst10,
                          const double& pnextdst10 );

     void   showecd( const int& pdcmnt );
     
     void   showleaf( const int& pdcmnt );

     void   updateC2N( const int& pdcmnt,
                       const double& yreet,
                       const double& yrpet,
                       const double& currentco2,
                       const double& initco2 );

     void updateDynamics( const double& lat,
    		 	 	 	  const double& lon,
    		 	 	 	  const int& month,
    		 	 	 	  const int& pdcmnt,
                          const double& co2,
                          const double& aot40,
						  const double& temptair,
                          const InorgN60& ndep,
                          const double& extraNH4,
                          const double& par,
                          const double& pet,
                          const double& prevmaxpet,
                          const double& eet,
                          const double& prevmaxeet,
                          const double& vegc,
                          const double& folc,
                          const double& stemc,
                          const double& crootc,
                          const double& frootc,
                          const double& deadwoodc,
                          const double& structn,
                          const double& labilen,
                          const double& deadwoodn,
                          const double& soilh2o,
                          const double& soilnh4,
                          const double& soilno3,
                          const int& moistlim,
                          const int& nfeed,
                          const int& o3flag,
                          const int& agstate,
                          const int& agperennial,
                          const int& fertflag,
                          const double& ksoil,
                          const double& netnmin,
                          const double& ammnvol,
                          const double& nitrif,
                          const double& no3prod,
                          double& agfertn,
                          double& age);

     void updateFoliage( const int& pdcmnt,
    		 	 	 	 const double& leafc,
						 const double& vegc,
						 const double& eet_mature );



     // "Get" and "Set" private variables and parameters
     
     // abvgrndResp ********************************************

     inline double getABVGPR( void ) { return abvgrndResp; };


     // adjc2n *************************************************

     inline double getADJC2N( void ) { return adjc2n; } 

     inline void setADJC2N( const double& padjc2n ) 
     { 
       adjc2n = padjc2n; 
     };


     // aleaf **************************************************
     
     inline double getALEAF( const int& pcmnt ) 
     { 
       return aleaf[pcmnt]; 
     };

     inline void setALEAF( const double& paleaf, 
                           const int& pcmnt ) 
     { 
       aleaf[pcmnt] = paleaf; 
     };

 
     // alpha *****************************************************

     inline double getALPHA( const int& pcmnt )
     {
       return alpha[pcmnt];
     };

     inline void setALPHA( const double& palpha, const int& pcmnt )
     {
       alpha[pcmnt] = palpha;
     };


     // beta *****************************************************

     inline double getBETA( const int& pcmnt )
     {
       return beta[pcmnt];
     };

     inline void setBETA( const double& pbeta, const int& pcmnt )
     {
       beta[pcmnt] = pbeta;
     };


     // bleaf **************************************************
     
     inline double getBLEAF( const int& pcmnt ) 
     { 
       return bleaf[pcmnt]; 
     };

     inline void setBLEAF( const double& pbleaf, 
                           const int& pcmnt ) 
     { 
       bleaf[pcmnt] = pbleaf; 
     };

 
     // c2n ****************************************************
     
     inline double getC2N( void ) { return c2n; };
     
     inline void setC2N( const double& pc2n ) { c2n = pc2n; };


     // c2na ***************************************************
     
     inline double getC2NA( const int& pcmnt ) 
     {
       return c2na[pcmnt]; 
     };
     
     inline void setC2NA( const double& pc2na, 
                          const int& pcmnt ) 
     { 
       c2na[pcmnt] = pc2na; 
     };


     // c2nb ***************************************************
     
     inline double getC2NB( const int& pcmnt ) 
     { 
       return c2nb[pcmnt]; 
     };
     
     inline void setC2NB( const double& pc2nb, 
                          const int& pcmnt ) 
     { 
       c2nb[pcmnt] = pc2nb; 
     };


     // c2nmin *************************************************
     
     inline double getC2NMIN( const int& pcmnt ) 
     { 
       return c2nmin[pcmnt]; 
     };
     
     inline void setC2NMIN( const double& pc2nmin, 
                            const int& pcmnt ) 
     { 
       c2nmin[pcmnt] = pc2nmin; 
     };


     // cfall **************************************************
     
     inline double getCFALL( const int& pcmnt ) 
     { 
       return cfall[pcmnt]; 
     };

     inline void setCFALL( const double& pcfall, 
                           const int& pcmnt ) 
     { 
       cfall[pcmnt] = pcfall; 
     };

     inline double getLCFALL( const int& pcmnt )
     {
       return lcfall[pcmnt];
     };

     inline void setLCFALL( const double& pcfall,
                           const int& pcmnt )
     {
       lcfall[pcmnt] = pcfall;
     };

     inline double getSCFALL( const int& pcmnt )
     {
       return scfall[pcmnt];
     };

     inline void setSCFALL( const double& pcfall,
                           const int& pcmnt )
     {
       scfall[pcmnt] = pcfall;
     };

     inline double getFCFALL( const int& pcmnt )
     {
       return fcfall[pcmnt];
     };

     inline void setFCFALL( const double& pcfall,
                           const int& pcmnt )
     {
       fcfall[pcmnt] = pcfall;
     };

     inline double getCCFALL( const int& pcmnt )
     {
       return ccfall[pcmnt];
     };

     inline void setCCFALL( const double& pcfall,
                           const int& pcmnt )
     {
       ccfall[pcmnt] = pcfall;
     };

     // cleaf **************************************************
     
     inline double getCLEAF( const int& pcmnt ) 
     { 
       return cleaf[pcmnt]; 
     };

     inline void setCLEAF( const double& pcleaf, 
                           const int& pcmnt ) 
     { 
       cleaf[pcmnt] = pcleaf; 
     };


 
     // cmax ***************************************************
     
     inline double getCMAX( void ) { return cmax; };

     inline void setCMAX( const double& pcmax ) 
     { 
       cmax = pcmax; 
     };


     // cmaxcut ************************************************
     
     inline double getCMAXCUT( const int& pcmnt ) 
     { 
       return cmaxcut[pcmnt]; 
     };

     inline void setCMAXCUT( const double& pcmaxcut, 
                             const int& pcmnt ) 
     { 
       cmaxcut[pcmnt] = pcmaxcut; 
     };


     // cmax1a *************************************************
     
     inline double getCMAX1A( const int& pcmnt ) 
     { 
       return cmax1a[pcmnt]; 
     };

     inline void setCMAX1A( const double& pcmax1a, 
                            const int& pcmnt ) 
     { 
       cmax1a[pcmnt] = pcmax1a; 
     };


     // cmax1b *************************************************
     
     inline double getCMAX1B( const int& pcmnt ) 
     { 
       return cmax1b[pcmnt]; 
     };

     inline void setCMAX1B( const double& pcmax1b, 
                            const int& pcmnt ) 
     { 
       cmax1b[pcmnt] = pcmax1b; 
     };


     // cmax2a *************************************************
     
     inline double getCMAX2A( const int& pcmnt ) 
     { 
       return cmax2a[pcmnt]; 
     };

     inline void setCMAX2A( const double& pcmax2a, 
                            const int& pcmnt ) 
     { 
       cmax2a[pcmnt] = pcmax2a; 
     };


     // cmax2b *************************************************
     
     inline double getCMAX2B( const int& pcmnt ) 
     { 
       return cmax2b[pcmnt]; 
     };

     inline void setCMAX2B( const double& pcmax2b, 
                            const int& pcmnt ) 
     { 
       cmax2b[pcmnt] = pcmax2b; 
     };


     // cneven *************************************************
     
     inline double getCNEVEN( void ) { return cneven; };

     inline void setCNEVEN( const double& pcneven ) 
     { 
       cneven = pcneven; 
     };


     // cnmin **************************************************
     
     inline double getCNMIN( const int& pcmnt ) 
     { 
       return cnmin[pcmnt]; 
     };

     inline void setCNMIN( const double& pcnmin, 
                           const int& pcmnt ) 
     { 
       cnmin[pcmnt] = pcnmin; 
     };


     // cov ****************************************************
     
     inline double getCOV( const int& pcmnt ) 
     { 
       return cov[pcmnt]; 
     };

     inline void setCOV( const double& pcov, const int& pcmnt ) 
     { 
       cov[pcmnt] = pcov; 
     };


     // currentveg *********************************************
     
     inline int getCURRENTVEG( void ) { return currentveg; };

     inline void setCURRENTVEG( const int& ptveg ) 
     { 
       currentveg = ptveg; 
     };
     
     
     // dc2n ***************************************************
     
     inline double getDC2N( void ) { return dc2n; };

     inline void setDC2N( const double& pdc2n ) 
     { 
       dc2n = pdc2n; 
     };
     
     
     // findozone **********************************************
     
     inline double getFINDOZONE( void ) { return findozone; };


     // foliage *************************************************
     
     inline double getFOLIAGE( void ) { return foliage; };


     // fozone *************************************************
     
     inline double getFOZONE( void ) { return fozone; };


     // fpc ****************************************************
     
     inline double getFPC( void ) { return fpc; };


     // fpcmax *************************************************
     
     inline double getFPCMAX( const int& pcmnt ) 
     { 
       return fpcmax[pcmnt]; 
     };

     inline void setFPCMAX( const double& pfpcmx, 
                            const int& pcmnt ) 
     { 
       fpcmax[pcmnt] = pfpcmx; 
     };


     // fprevozone *********************************************
     
     inline double getFPREVOZONE( void ) { return fprevozone; };

     inline void setFPREVOZONE( const double& pfprevozone ) 
     { 
       fprevozone = pfprevozone; 
     };


     // gamma *****************************************************

     inline double getGAMMA( const int& pcmnt )
     {
       return gamma[pcmnt];
     };

     inline void setGAMMA( const double& pgamma, const int& pcmnt )
     {
       gamma[pcmnt] = pgamma;
     };


     // gpp ****************************************************
     
     inline double getGPP( void ) { return gpp; };
     inline double getNPP_leaf( void ) { return npp_leaf; };
     inline double getNPP_stem( void ) { return npp_stem; };
     inline double getNPP_croot( void ) { return npp_croot; };
     inline double getNPP_froot( void ) { return npp_froot; };

     // gpr ****************************************************
     
     inline double getGPR( void ) { return gpr; };


     // gva ****************************************************
     
     inline double getGVA( const int& pcmnt ) 
     { 
       return gva[pcmnt]; 
     };

     inline void setGVA( const double& pgva, 
                         const int& pcmnt ) 
     { 
       gva[pcmnt] = pgva; 
     };

     // ingpp **************************************************
     
     inline double getINGPP( void ) { return ingpp; };


     // initcneven *********************************************
     
     inline double getINITCNEVEN( const int& pcmnt ) 
     { 
       return initcneven[pcmnt]; 
     };

     inline void setINITCNEVEN( const double& pincneven, 
                                const int& pcmnt ) 
     { 
       initcneven[pcmnt] = pincneven; 
     };


     // initleafmx *********************************************
     
     inline double getINITLEAFMX( const int& pcmnt ) 
     { 
      return initleafmx[pcmnt]; 
     };

     inline void setINITLEAFMX( const double& pinleafmx, 
                                const int& pcmnt ) 
     { 
       initleafmx[pcmnt] = pinleafmx; 
     };


     // innpp **************************************************
     
     inline double getINNPP( void ) { return innpp; };


     // inuptake.nh4 *******************************************
     
     inline double getINH4UPTAKE( void ) 
     { 
       return inuptake.nh4; 
     };

     
     // inuptake.no3 *******************************************
     
     inline double getINO3UPTAKE( void ) 
     { 
       return inuptake.no3; 
     };

     
     // inuptake.total *****************************************
     
     inline double getINUPTAKE( void ) 
     { 
       return inuptake.total; 
     };


     // kc *****************************************************
     
     inline double getKC( const int& pcmnt ) 
     { 
       return kc[pcmnt]; 
     };

     inline void setKC( const double& pkc, const int& pcmnt ) 
     { 
       kc[pcmnt] = pkc; 
     };


     // ki *****************************************************
     
     inline double getKI( const int& pcmnt ) 
     { 
       return ki[pcmnt]; 
     };

     inline void setKI( const double& pki, const int& pcmnt ) 
     { 
       ki[pcmnt] = pki; 
     };


     // kn1 ****************************************************
     
     inline double getKN1( const int& pcmnt ) 
     { 
       return kn1[pcmnt]; 
     };

     inline void setKN1( const double& pkn1, const int& pcmnt ) 
     { 
       kn1[pcmnt] = pkn1; 
     };


     // kra ****************************************************
     
     inline double getKRA( const int& pcmnt ) 
     { 
       return kra[pcmnt]; 
     };

     inline void setKRA( const double& pkra, const int& pcmnt ) 
     { 
       kra[pcmnt] = pkra; 
     };


     // krb ****************************************************
     
     inline double getKRB( const int& pcmnt ) 
     { 
       return krb[pcmnt]; 
     };

     inline void setKRB( const double& pkrb, const int& pcmnt ) 
     { 
       krb[pcmnt] = pkrb; 
     };


     // kleaf **************************************************
     
     inline double getKLEAFC( const int& pcmnt ) 
     { 
       return kleafc[pcmnt]; 
     };

     inline void setKLEAFC( const double& pkleafc, 
                            const int& pcmnt ) 
     { 
       kleafc[pcmnt] = pkleafc; 
     };


     // kvnh4 ****************************************************
     
     inline double getKVNH4( const int& pcmnt ) 
     { 
       return kvnh4[pcmnt]; 
     };

     inline void setKVNH4( const double& pkvnh4, 
                           const int& pcmnt ) 
     { 
       kvnh4[pcmnt] = pkvnh4; 
     };


     // kvno3 ****************************************************
     
     inline double getKVNO3( const int& pcmnt ) 
     { 
       return kvno3[pcmnt]; 
     };

     inline void setKVNO3( const double& pkvno3, 
                           const int& pcmnt ) 
     { 
       kvno3[pcmnt] = pkvno3; 
     };


     // labile.nitrogen ****************************************
     
     inline double getLABILEN( void ) 
     { 
       return labile.nitrogen; 
     };


     // labncon ************************************************
     
     inline double getLABNCON( const int& pcmnt ) 
     { 
       return labncon[pcmnt]; 
     };

     inline void setLABNCON( const double& plabncon, 
                             const int& pcmnt ) 
     { 
       labncon[pcmnt] = plabncon; 
     };


     // lai ****************************************************
     
     inline double getLAI( void ) { return lai; };


     // lcclnc *************************************************
     
     inline double getLCCLNC( const int& pcmnt ) 
     { 
       return lcclnc[pcmnt]; 
     }

     inline void setLCCLNC( const double& plcclnc, 
                            const int& pcmnt ) 
     { 
       lcclnc[pcmnt] = plcclnc; 
     }


     // leaf ***************************************************
     
     inline double getLEAF( void ) { return leaf; };
     

     // leafmxc ************************************************
     
     inline double getLEAFMXC( const int& pcmnt ) 
     { 
       return leafmxc[pcmnt]; 
     };

     inline void setLEAFMXC( const double& pleafmxc, 
                             const int& pcmnt ) 
     { 
       leafmxc[pcmnt] = pleafmxc; 
     };


     // ltrfal.carbon ******************************************  
     
     inline double getLTRFALC( void ) { return ltrfal.carbon; };

     
     // ltrfal.nitrogen ****************************************
     
     inline double getLTRFALN( void ) 
     { 
       return ltrfal.nitrogen; 
     };


     // luptake ************************************************
     
     inline double getLUPTAKE( void ) { return luptake; };


     // minleaf ************************************************
     
     inline double getMINLEAF( const int& pcmnt ) 
     { 
       return minleaf[pcmnt]; 
     };

     inline void setMINLEAF( const double& pminleaf, 
                             const int& pcmnt ) 
     { 
       minleaf[pcmnt] = pminleaf; 
     };

 
     // newleafmx **********************************************
     
     inline double getNEWLEAFMX( void ) { return newleafmx; };

     inline void setNEWLEAFMX( const double& pnewleafmx ) 
     { 
       newleafmx = pnewleafmx; 
     };


     // newtopt ************************************************
     
     inline double getNEWTOPT( void ) { return newtopt; };

     inline void setNEWTOPT( const double& pnewtopt ) 
     { 
       newtopt = pnewtopt; 
     };


     // nfall **************************************************
     
     inline double getNFALL( const int& pcmnt ) 
     { 
       return nfall[pcmnt]; 
     };

     inline void setNFALL( const double& pnfall, 
                           const int& pcmnt ) 
     { 
       nfall[pcmnt] = pnfall; 
     };


     // nfix ***************************************************
     
     inline double getNFIX( void ) { return nfix; };

     inline void setNFIX( const double& pnfix ) 
     { 
       nfix = pnfix; 
     };


     // nfixpara ***********************************************
     
     inline double getNFIXPARA( const int& pcmnt ) 
     { 
       return nfixpara[pcmnt]; 
     };

     inline void setNFIXPARA( const double& pnfixpara, 
                              const int& pcmnt ) 
     { 
       nfixpara[pcmnt] = pnfixpara; 
     };


     // nfixparb ***********************************************
     
     inline double getNFIXPARB( const int& pcmnt ) 
     { 
       return nfixparb[pcmnt]; 
     };

     inline void setNFIXPARB( const double& pnfixparb, 
                              const int& pcmnt ) 
     { 
       nfixparb[pcmnt] = pnfixparb; 
     };
 
 
     // nmobil *************************************************
     
     inline double getNMOBIL( void ) { return nmobil; };


     // npp ****************************************************
     
     inline double getNPP( void ) { return npp; };


     // nresorb ************************************************
     
     inline double getNRESORB( void ) { return nresorb; };


     // nupnh4 *************************************************
     
     inline double getNUPNH4( void ) { return nupnh4; };

     inline void setNUPNH4( const double& pnupnh4 ) 
     { 
       nupnh4 = pnupnh4; 
     };


     // nupnh4cut **********************************************
     
     inline double getNUPNH4CUT( const int& pcmnt ) 
     { 
       return nupnh4cut[pcmnt]; 
     };

     inline void setNUPNH4CUT( const double& pnupnh4cut, 
                               const int& pcmnt ) 
     { 
       nupnh4cut[pcmnt] = pnupnh4cut; 
     };


     // nupnh41a ***********************************************
     
     inline double getNUPNH41A( const int& pcmnt ) 
     { 
       return nupnh41a[pcmnt]; 
     };

     inline void setNUPNH41A( const double& pnupnh41a, 
                              const int& pcmnt ) 
     { 
       nupnh41a[pcmnt] = pnupnh41a; 
     };


     // nupnh41b ***********************************************
     
     inline double getNUPNH41B( const int& pcmnt ) 
     { 
       return nupnh41b[pcmnt]; 
     };

     inline void setNUPNH41B( const double& pnupnh41b, 
                              const int& pcmnt ) 
     { 
       nupnh41b[pcmnt] = pnupnh41b; 
     };


     // nupnh42a ***********************************************
     
     inline double getNUPNH42A( const int& pcmnt ) 
     { 
       return nupnh42a[pcmnt]; 
     };

     inline void setNUPNH42A( const double& pnupnh42a, 
                              const int& pcmnt ) 
     { 
       nupnh42a[pcmnt] = pnupnh42a; 
     };


     // nupnh42b ***********************************************
     
     inline double getNUPNH42B( const int& pcmnt ) 
     { 
       return nupnh42b[pcmnt]; 
     };

     inline void setNUPNH42B( const double& pnupnh42b, 
                              const int& pcmnt ) 
     { 
       nupnh42b[pcmnt] = pnupnh42b; 
     };


     // nupno3 *************************************************
     
     inline double getNUPNO3( void ) { return nupno3; };

     inline void setNUPNO3( const double& pnupno3 ) 
     { 
       nupno3 = pnupno3; 
     };


     // nupno3cut **********************************************
     
     inline double getNUPNO3CUT( const int& pcmnt ) 
     { 
       return nupno3cut[pcmnt]; 
     };

     inline void setNUPNO3CUT( const double& pnupno3cut, 
                               const int& pcmnt ) 
     { 
       nupno3cut[pcmnt] = pnupno3cut; 
     };


     // nupno31a ***********************************************
     
     inline double getNUPNO31A( const int& pcmnt ) 
     { 
       return nupno31a[pcmnt]; 
     };

     inline void setNUPNO31A( const double& pnupno31a, 
                              const int& pcmnt ) 
     { 
       nupno31a[pcmnt] = pnupno31a; 
     };
  
  
     // nupno31b ***********************************************
     
     inline double getNUPNO31B( const int& pcmnt ) 
     { 
       return nupno31b[pcmnt]; 
     };

     inline void setNUPNO31B( const double& pnupno31b, 
                              const int& pcmnt ) 
     { 
       nupno31b[pcmnt] = pnupno31b; 
     };


     // nupno32a ***********************************************
     
     inline double getNUPNO32A( const int& pcmnt ) 
     { 
       return nupno32a[pcmnt]; 
     };

     inline void setNUPNO32A( const double& pnupno32a, 
                              const int& pcmnt ) 
     { 
       nupno32a[pcmnt] = pnupno32a; 
     };


     // nupno32b ***********************************************
     
     inline double getNUPNO32B( const int& pcmnt ) 
     { 
       return nupno32b[pcmnt]; 
     };

     inline void setNUPNO32B( const double& pnupno32b, 
                              const int& pcmnt ) 
     { 
       nupno32b[pcmnt] = pnupno32b; 
     };


     // nuptake.nh4 ********************************************
     
     inline double getNH4UPTAKE( void ) { return nuptake.nh4; };

     
     // nuptake.no3 ********************************************
     
     inline double getNO3UPTAKE( void ) { return nuptake.no3; };

     inline void setNO3UPTAKE( const double& pno3uptake  ) 
     { 
       nuptake.no3 = pno3uptake; 
     };

     // nuptake.total ******************************************
     
     inline double getNUPTAKE( void ) { return nuptake.total; };


     // o3para *************************************************
     
     inline double getO3PARA( const int& pcmnt ) 
     { 
       return o3para[pcmnt]; 
     };

     inline void setO3PARA( const double& po3para, 
                            const int& pcmnt ) 
     { 
       o3para[pcmnt] = po3para; 
     };


     // o3parb *************************************************
     
     inline double getO3PARB( const int& pcmnt ) 
     { 
       return o3parb[pcmnt]; 
     };

     inline void setO3PARB( const double& po3parb, 
                            const int& pcmnt ) 
     { 
       o3parb[pcmnt] = po3parb; 
     };


     // o3parc *************************************************
     
     inline double getO3PARC( const int& pcmnt ) 
     { 
       return o3parc[pcmnt]; 
     };

     inline void setO3PARC( const double& po3parc, 
                            const int& pcmnt ) 
     { 
       o3parc[pcmnt] = po3parc; 
     };


     // potveg *************************************************
     
     inline int getPOTVEG( void ) { return potveg; };

     inline void setPOTVEG( const int& ptveg ) 
     { 
       potveg = ptveg; 
     };


     // prevunrmleaf *******************************************
     
     inline double getPREVUNRMLEAF( void ) 
     { 
       return prevunrmleaf; 
     };

     inline void setPREVUNRMLEAF( const double& pprevunrmleaf ) 
     { 
       prevunrmleaf = pprevunrmleaf; 
     };


     // proptrans ************************************************
     
     inline double getPROPTRANS( const int& pcmnt ) 
     { 
       return proptrans[pcmnt]; 
     };

     inline void setPROPTRANS( const double& pproptrans, 
                               const int& pcmnt ) 
     { 
       proptrans[pcmnt] = pproptrans; 
     };


     // prvleafmx **********************************************
     
     inline double getPRVLEAFMX( void ) { return prvleafmx; };

     inline void setPRVLEAFMX( const double& pprvleafmx ) 
     { 
       prvleafmx = pprvleafmx; 
     };


     // qref *****************************************************

     inline double getQREF( const int& pcmnt )
     {
       return qref[pcmnt];
     };

     inline void setQREF( const double& pqref, const int& pcmnt )
     {
       qref[pcmnt] = pqref;
     };


     // raq10a0 ************************************************
     
//     inline double getRAQ10A0( const int& pcmnt ) 
//     { 
//       return raq10a0[pcmnt]; 
//     };

//     inline void setRAQ10A0( const double& praq10a0, 
//                             const int& pcmnt ) 
//     { 
//       raq10a0[pcmnt] = praq10a0; 
//     };


     // raq10a1 ************************************************
     
//     inline double getRAQ10A1( const int& pcmnt ) 
//     { 
//       return raq10a1[pcmnt]; 
//     };

//     inline void setRAQ10A1( const double& praq10a1, 
//                             const int& pcmnt ) 
//     { 
//       raq10a1[pcmnt] = praq10a1; 
//     };


     // raq10a2 ************************************************
     
//     inline double getRAQ10A2( const int& pcmnt ) 
//     { 
//       return raq10a2[pcmnt]; 
//     };

//     inline void setRAQ10A2( const double& praq10a2, 
//                             const int& pcmnt ) 
//     { 
//       raq10a2[pcmnt] = praq10a2; 
//     };


     // raq10a3 ************************************************
     
//     inline double getRAQ10A3( const int& pcmnt ) 
//     { 
//       return raq10a3[pcmnt]; 
//     };

//     inline void setRAQ10A3( const double& praq10a3, 
//                             const int& pcmnt ) 
//     { 
//       raq10a3[pcmnt] = praq10a3; 
//     };


     // rg *****************************************************
     
     inline double getRGRWTH( void ) { return rg; };


     // rm *****************************************************
     
     inline double getRMAINT( void ) { return rm; };


 	   // rmmax ****************************************************

	   inline double getRMMAX( const int& pcmnt )
	   {
	     return rmmax[pcmnt];
	   };

	   inline void setRMMAX( const double& prmmax, const int& pcmnt )
	   {
	     rmmax[pcmnt] = prmmax;
       };


     // rootResp ***********************************************
     
     inline double getROOTRESP( void ) { return rootResp; };


     // rroot **************************************************
     
     inline double getRROOT( const int& pcmnt ) 
     { 
       return rroot[pcmnt]; 
     };

     inline void setRROOT( const double& prroot, 
                           const int& pcmnt ) 
     { 
       rroot[pcmnt] = prroot; 
     };


     // sla ****************************************************
     
     inline double getSLA( const int& pcmnt ) 
     { 
       return sla[pcmnt]; 
     };

     inline void setSLA( const double& psla, const int& pcmnt ) 
     { 
       sla[pcmnt] = psla; 
     };


     // strctrl.nitrogen ***************************************
          
     inline double getSTRUCTN( void ) 
     { 
       return strctrl.nitrogen; 
     };


     // subtype ************************************************
          
     inline int getSUBTYPE( void ) { return subtype; };

     inline void setSUBTYPE( const int& psubtype ) 
     { 
       subtype = psubtype; 
     };


     // suptake ************************************************
          
     inline double getSUPTAKE( void ) { return suptake; };


     // thawpercent ********************************************
     
     inline double getTHAWPCT( void ) { return thawpercent; };

     inline void setThawPercent( const double& pthawpct ) 
     { 
       thawpercent = pthawpct; 
     };


     // tmax ***************************************************
     
     inline double getTMAX( const int& pcmnt ) 
     { 
       return tmax[pcmnt]; 
     };

     inline void setTMAX( const double& ptmax, 
                          const int& pcmnt ) 
     { 
       tmax[pcmnt] = ptmax; 
     };


     // tmin ***************************************************
     
     inline double getTMIN( const int& pcmnt ) 
     { 
       return tmin[pcmnt]; 
     };

     inline void setTMIN( const double& ptmin, 
                          const int& pcmnt ) 
     { 
       tmin[pcmnt] = ptmin; 
     };


     // topt ***************************************************
     
     inline double getTOPT( void ) { return topt; };

     inline void setTOPT( const double& ptopt ) 
     { 
       topt = ptopt; 
     };


     // toptmax ************************************************
     
     inline double getTOPTMAX( const int& pcmnt ) 
     { 
       return toptmax[pcmnt]; 
     };

     inline void setTOPTMAX( const double& ptoptmax, 
                             const int& pcmnt ) 
     { 
       toptmax[pcmnt] = ptoptmax; 
     };


     // toptmin ************************************************
     
     inline double getTOPTMIN( const int& pcmnt ) 
     { 
       return toptmin[pcmnt]; 
     };

     inline void setTOPTMIN( const double& ptoptmin, 
                             const int& pcmnt ) 
     { 
       toptmin[pcmnt] = ptoptmin; 
     };


     // transpiration *****************************************
      
     inline double getTRANSPIRATION( void ) 
  	 { 
	     return transpiration; 
	   };

     inline void setTRANSPIRATION( const double& ptrans ) 
     { 
       transpiration = ptrans; 
     };


     // tref *****************************************************

     inline double getTREF( const int& pcmnt )
     {
       return tref[pcmnt];
     };

     inline void setTREF( const double& ptref, const int& pcmnt )
     {
       tref[pcmnt] = ptref;
     };


     // unleaf12 ***********************************************
     
     inline double getUNLEAF12( const int& pcmnt ) 
     { 
       return unleaf12[pcmnt]; 
     };

     inline void setUNLEAF12( const double& punleaf12, 
                              const int& pcmnt ) 
     { 
       unleaf12[pcmnt] = punleaf12; 
     };


     // unnormleaf *********************************************
     
     inline double getUNNORMLEAF( void ) { return unnormleaf; };


     // plant.carbon *******************************************
     
     inline double getVEGC( void ) { return plant.carbon; };


     // plant.nitrogen *****************************************
      
     inline double getVEGN( void ) { return plant.nitrogen; };

     inline void setVEGN( const double& pvegn ) 
     { 
       plant.nitrogen = pvegn; 
     };



/* *************************************************************
		 Public Variables
************************************************************* */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif

     // Index for community type
     int cmnt;


     // Annual sum of abvgrndResp
     double yrabvgrndResp;

     // ratio of yrcarbon to yrnitrogen
     double yrc2n;            

     // Annual sum of plant.carbon
     double yrcarbon;          

     // Sum of monthly FPC
     double yrfpc;        

     // annual sum of monthly GPP
     double yrgpp;             

      // Annual GPR
     double yrgpr;            

     // Annual sum of ingpp
     double yringpp;           

     // Annual sum of innfix
     double yrinnfix;    

     // Annual sum of innpp
     double yrinnpp;           

     // Annual sum of innup
     InorgN60 yrinnup;           

     double yrinpr;

      // Sum of monthly LAI
     double yrlai;            

     // Annual sum of lnfix
     double yrlnfix;           

     // mean annual normalized leaf phenology
     double yrleaf;      

     // Annual sum of ltrfal.carbon
     double yrltrc;      
     
     // Annual sum of ltrfal.nitrogen      
     double yrltrn;            

     // Annual sum of luptake
     double yrlup;             

     // Annual sum of nfix
     double yrnfix;            

     // Annual sum of plant.nitrogen
     double yrnitrogen;        

     // Annual sum of nmobil
     double yrnmobil;          

     // Annual sum of npp
     double yrnpp;             

     // Annual sum of nresorb
     double yrnrsorb;          

     // Annual sum of nuptake
     InorgN60 yrnup;             

     double yrprod;

     // Annual sum of rg
     double yrrgrowth;

     // Annual sum of rm
     double yrrmaint;

     // Annual sum of rootResp
     double yrrootResp;

     // Annual sum of snfix
     double yrsnfix;           

     // Annual sum of labile.nitrogen
     double yrstoren;          
     
      // Annual sum of strctrl.nitrogen
     double yrstructn;       
     
     // Annual sum of suptake
     double yrsup;             

     // Annual sum of thawpercent
     double yrthawpct;

     // Annual sum of unnormleaf
     double yrunleaf;          

     // Annual sum of VOC fluxes
     VolatileOrganicCarbon60 yrvoc;
     double lfraction;

     double rootstemratio; //allocation coefficient between root and stem C

     double yrdeadwoodc;
     double yrdeadwoodn;
     inline double getrsratio( const int& pcmnt )
     {
       return rsratio[pcmnt];
     };

     inline void setrsratio( const double& prsratio,
                               const int& pcmnt )
     {
    	 rsratio[pcmnt] = prsratio;
     };

     inline int getifwoody( const int& pcmnt )
     {
       return ifwoody[pcmnt];
     };

     inline void setifwoody( const int& pifwoody,
                               const int& pcmnt )
     {
    	 ifwoody[pcmnt] = pifwoody;
     };
     inline double getlfratiomax( const int& pcmnt )
     {
       return lfratiomax[pcmnt];
     };

     inline void setlfratiomax( const double& plfratiomax,
                               const int& pcmnt )
     {
    	 lfratiomax[pcmnt] = plfratiomax;
     };
     inline double getwoodfall( const int& pcmnt )
     {
       return woodfall[pcmnt];
     };

     inline void setwoodfall( const double& pwoodfall,
                               const int& pcmnt )
     {
    	 woodfall[pcmnt] = pwoodfall;
     };
     /*
     inline void setdeadwoodltc( const double& pdeadwoodltc)
     {
    	 deadwoodltc = pdeadwoodltc;
     };
     inline double getdeadwoodltc(void)
     {
    	 return deadwoodltc;
     }

     inline void setdeadwoodltn( const double& pdeadwoodltn)
     {
    	 deadwoodltn = pdeadwoodltn;
     };
     inline double getdeadwoodltn(void)
     {
    	 return deadwoodltn;
     }
     */
     inline double getstanddeadpar( const int& pcmnt )
     {
         return standdeadpar[pcmnt];
     };

     inline void setstanddeadpar( const double& pstanddeadpar,
                                   const int& pcmnt )
     {
         standdeadpar[pcmnt] = pstanddeadpar;     };

     inline double getMRCOEF( const int& pcmnt )
     {
        return mrcoef[pcmnt];
     };

     inline void setMRCOEF( const double& mcoef, const int& pcmnt )
     {
        	 mrcoef[pcmnt] = mcoef;
     };
     inline double getAGECOEF( const int& pcmnt )
     {
        return agecoef[pcmnt];
     };

     inline void setAGECOEF( const double& fagecoef,
                               const int& pcmnt )
     {
         agecoef[pcmnt] = fagecoef;
     };

     inline double getMORTCOEFA( const int& pcmnt )
     {
        return mortcoefa[pcmnt];
     };

     inline void setMORTCOEFA( const double& fmortcoefa,
                               const int& pcmnt )
     {
         mortcoefa[pcmnt] = fmortcoefa;
     };
     inline double getMORTCOEFB( const int& pcmnt )
     {
        return mortcoefb[pcmnt];
     };

     inline void setMORTCOEFB( const double& fmortcoefb,
                               const int& pcmnt )
     {
         mortcoefb[pcmnt] = fmortcoefb;
     };
     double deadwoodltc, deadwoodltn;



  private:

/* **************************************************************
		 Private Functions
************************************************************** */

     double deltaleaf( const int& pdcmnt,
                       const double& eet,
                       const double& prveetmx,
                       const double& prvleaf );

     double gppxclm( const int& pdcmnt,
                     const double& co2,
                     const double& par,
                     const double& temp,
                     const double& gv,
                     const double& leaf,
                     const double& foliage,
                     const double& thawpercent );

     double gppxio3( const double& fozone, 
                     const double& eetpet );

     double gppxo3( const int& pdcmnt,
                    const double& gpp,
                    const double& d40,
                    const double& eetpet );

     double nupnh4xclm( const int& pdcmnt,
                        const double& soilh2o,
                        const double& nh4,
                        const double& respq10,
                        const double& ksoil,
                        const double& foliage,
                        const double& fozone );

     double nupno3xclm( const int& pdcmnt,
                        const double& soilh2o,
                        const double& no3,
                        const double& respq10,
                        const double& ksoil,
                        const double& foliage,
                        const double& fozone );

     double nupxclm( const int& pdcmnt,
                     const double& soilh2o,
                     const double& availn,
                     const double& respq10,
                     const double& ksoil,
                     const double& foliage,
		     const double& fozone );

//     double rq10( const int& pdcmnt, const double& tair );

     double rmxclm( const int& pdcmnt,
                    const double& tcarbon,
                    const double& mrcoef,
                    const double& respq10 );

/* **************************************************************
		 Private Variables
************************************************************** */
     
     // Aboveground plant respiration
     double abvgrndResp; 

     double alleaf;

     // Index for current vegetation type
     int currentveg;
     
     // Multiplier of indirect ozone effects on GPP
     double findozone;

     double foliage;

     // Multiplier of direct ozone effects on GPP
     double fozone;

     // Monthly foliar projective cover
     double fpc;        

     // Previous month's value of fozone
     double fprevozone;

     // Monthly gross primary productivity (GPP)
     double gpp;        

     // Monthly gross plant respiration (rm + rg)
     double gpr;        

      // Initial monthly gross primary productivity
     double ingpp;     

     // Initial symbiotic N fixation
     double innfix;            

     // Initial net primary productivity
     double innpp;      

     // Initial C/N of biomass production
     double inprodcn;          

     // Initial N uptake by plants
     InorgN60 inuptake;   

     // Labile plant biomass
     Biomass labile;    

     // Monthly leaf area index
     double lai;        

     // monthly normalized leaf phenology
     double leaf; 

     // Number of annual iterations for determining monthly
     //   phenology
     int leafyrs;
     
     double lnfix;

     // Monthly litterfall
     Biomass ltrfal;    

     // Monthly N uptake by plants for labile N
     double luptake;    

     // Updated maximum leaf for current year 
     double newleafmx;

     // Updated optimum air temperature for current year
     double newtopt;

     // Monthly symbiotic N fixation
     double nfix;       

     // Monthly N mobilization by plants
     double nmobil;     

     // Monthly net primary productivity (NPP)
     double npp;        

     // Monthly N resorption by plants
     double nresorb;    

     // Monthly N uptake by plants
     InorgN60 nuptake;    

     // whole plant biomass (structural + labile)
     Biomass plant;     

     // Index for potential vegetation biome type
     int potveg;

     // Unnormalized relative leaf area of previous month
     double prevunrmleaf;

     // Maximum relative leaf area of previous year
     double prvleafmx;

     // Effect of air temperature on plant respiration
     double respq10;

     // Monthly growth respiration
     double rg;         

     // Monthly maintenance respiration
     double rm;         

     // Monthly root respiration
     double rootResp;   

     // Structural plant biomass
     Biomass strctrl;   

     // Index for vegetation subtype
     int subtype;
     
     // Monthly N uptake by plants for structural N
     double suptake;    

     // Effect of air temperature on GPP (calculated based on tair)
     double temp;

     // Index for vegetation biome type
//     int temveg;

     // Percent of month with thawed ground
     double thawpercent; 

     // Canopy transpiration
	 double transpiration;

     // Monthly unnormalized leaf phenology
     double unnormleaf; 

     // Monthly VOC fluxes
     VolatileOrganicCarbon60 voc;
     double gv;

     double npp_leaf, npp_stem, npp_croot,npp_froot;


/* *************************************************************
		 Private Parameters
************************************************************* */

     // Biome-specific plant respiration parameters

     double alpha[MAXCMNT];
     double beta[MAXCMNT];
     double gamma[MAXCMNT];
     double qref[MAXCMNT];
     double tref[MAXCMNT];

     // Biome-specific vegetation C/N parameters

     double adjc2n;
     double c2n;
     double c2na[MAXCMNT];
     double c2nb[MAXCMNT];
     double c2nmin[MAXCMNT];
     double cnmin[MAXCMNT];
     double dc2n;

     double cneven;
     double initcneven[MAXCMNT];

     double cfall[MAXCMNT];  // proportion of vegetation carbon (stem+coarse root)
     double lcfall[MAXCMNT];  // proportion of leaf carbon

     double scfall[MAXCMNT];  // proportion of stem carbon
     double ccfall[MAXCMNT];  // proportion of coarse root carbon
     double fcfall[MAXCMNT];  // proportion of fine root carbon
     // Biome-specific carbon uptake parameters for function gppxclm

     double cmax;
     double cmaxcut[MAXCMNT];
     double cmax1a[MAXCMNT];
     double cmax1b[MAXCMNT];
     double cmax2a[MAXCMNT];
     double cmax2b[MAXCMNT];

     // Biome-specific respiration parameters for function rmxclm

     double kr;
     double kra[MAXCMNT];
     double krb[MAXCMNT];

     double lcclnc[MAXCMNT];

     // Biome-specific phenology parameters

     double aleaf[MAXCMNT];
     double bleaf[MAXCMNT];
     double cleaf[MAXCMNT];
     double initleafmx[MAXCMNT];
     double minleaf[MAXCMNT];
     double unleaf12[MAXCMNT];


     // Biome-specific foliage projection cover parameters

     double cov[MAXCMNT];
     double fpcmax[MAXCMNT];
     double sla[MAXCMNT];

     // Biome-specific parameter to describe the sensitivity of GPP
     //   to evapotranspiration
     
     double gva[MAXCMNT];

     // Biome-specific half saturation parameter for function 
     //   gppxclm describing the effects of solar atmospheric 
     //   carbon dioxide concentration on GPP

     double kc[MAXCMNT];

     // Biome-specific half saturation parameter for function 
     //   gppxclm describing the effects of photosybtheically 
     //   active radiation on GPP

     double ki[MAXCMNT];

     double kn1[MAXCMNT];
     double kvnh4[MAXCMNT];
     double kvno3[MAXCMNT];

     // Biome-specific allocation parameters

     double kleafc[MAXCMNT];
     double leafmxc[MAXCMNT];

     double labncon[MAXCMNT];

     double nfall[MAXCMNT];  // proportion of vegetation nitrogen

     // Biome-specific symbiotic N fixation parameters

     double nfixpara[MAXCMNT];
     double nfixparb[MAXCMNT];

     // Biome-specific nitrogen uptake parameters for function nupxclm

     double nmax;
     double nmaxcut[MAXCMNT];
     double nmax1a[MAXCMNT];
     double nmax1b[MAXCMNT];
     double nmax2a[MAXCMNT];
     double nmax2b[MAXCMNT];

     double nupnh4;
     double nupnh4cut[MAXCMNT];
     double nupnh41a[MAXCMNT];
     double nupnh41b[MAXCMNT];
     double nupnh42a[MAXCMNT];
     double nupnh42b[MAXCMNT];

     double nupno3;
     double nupno3cut[MAXCMNT];
     double nupno31a[MAXCMNT];
     double nupno31b[MAXCMNT];
     double nupno32a[MAXCMNT];
     double nupno32b[MAXCMNT];

     // Biome-specific ozone parameters

     double o3para[MAXCMNT];
     double o3parb[MAXCMNT];
     double o3parc[MAXCMNT];
	 
	 // Proportion of EET that is canopy conductance
	 double proptrans[MAXCMNT];


     // Biome-specific parameters for function rq10 to describe the
     //   effect of temperature on plant respiration

//     double raq10a0[MAXCMNT];
//     double raq10a1[MAXCMNT];
//     double raq10a2[MAXCMNT];
//     double raq10a3[MAXCMNT];

     // Maximum maintenance respiration rate
     
     double rmmax[MAXCMNT];

     // Biome-specific proportion of gpr lost as root respiration

     double rroot[MAXCMNT];

    // Element-specific optimum temperature for GPP

     double tmax[MAXCMNT];
     double tmin[MAXCMNT];
     double topt;
     double toptmax[MAXCMNT];
     double toptmin[MAXCMNT];
     //new parameters added by cgs in veg table
     double rsratio[MAXCMNT]; //root/stem ratio
     int ifwoody[MAXCMNT]; //woody or nonwoody (0,1)
     double lfratiomax[MAXCMNT]; //max allocation proportion to leaf
     double woodfall[MAXCMNT]; //standing deadwood falling rate

     double standdeadpar[MAXCMNT]; //temporarily used to replace LCLUC standdeadpar, will be deleted if LCLUC provide this data
     double mrcoef[MAXCMNT];//maintanence respiration coefficient
     double agecoef[MAXCMNT];//age effects on plant density
     double mortcoefa[MAXCMNT]; //age effects on mortality rate (slope)
     double mortcoefb[MAXCMNT]; //age effects on mortality rate (exponential decay coefficient)
     double mrcoef_leaf[MAXCMNT];//maintanence respiration coefficient
     double mrcoef_stem[MAXCMNT];//maintanence respiration coefficient
     double mrcoef_croot[MAXCMNT];//maintanence respiration coefficient
     double mrcoef_froot[MAXCMNT];//maintanence respiration coefficient
     double q10_cgs[MAXCMNT]; //q10 for gpp
};

#endif

