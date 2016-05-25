/* *************************************************************
****************************************************************
HUMNACT604.H - describes human disturbances to natural ecosystems

20030513 - DWK created by modifying humnact50a.h
20030513 - DWK changed include from tveg50a.h to tveg51.h
20030513 - DWK changed Humanact50 to Humanact51
20031016 - DWK replaced char ecd[MAXFNAME] with string ecd in
           function calls to getecd()
20031016 - DWK added include "tprocessXML431a.h"
20031016 - DWK added inheritance of ProcessXML to class Humanact51
20031016 - DWK deleted private functions getvalue() and readline()
20031019 - DWK changed include from tprocessXML431a.h to
           tprocessXML51.h
20031019 - DWK changed inheritance of ProcessXML to ProcessXML51
20040707 - DWK changed include from tveg51.h to tveg60.h
20040707 - DWK changed include from tprocessXML51.h to
           tprocessXML60.h
20040707 - DWK changed class Humanact51 to class Humanact60
20040707 - DWK changed inheritance of ProcessXML51 to 
           inheritance of ProcessXML60
20040707 - DWK changed include from humnact51.cpp to 
           humnact60.cpp at bottom of file                    
20040716 - DWK added public int FRF to represent fire return
           frequency
20040716 - DWK changed double slashpar[MAXCMNT] to 
           double slashpar
20040716 - DWK changed double vconvert[MAXCMNT] to 
           double vconvert
20040716 - DWK changed public double prod10par[MAXCMNT] to 
           double prod10par
20040716 - DWK changed public double prod100par[MAXCMNT] to 
           double prod100par
20040716 - DWK changed public double sconvert[MAXCMNT] to 
           double sconvert
20040716 - DWK changed public double vrespar[MAXCMNT] to 
           double vrespar
20051117 - DWK deleted include "humnact60.cpp" from bottom of file
20051117 - DWK changed include from tprocessXML60.h to
           tprocessXML602.h
20051117 - DWK changed include from tveg60.h to tveg602.h
20051117 - DWK added include temconsts602.hpp
20051112 - DWK added public function createWoodProducts()
20051125 - DWK changed public double cflux[CYCLE] to double cflux
20051125 - DWK changed public Biomass convrtflx[CYCLE] to 
           Biomass convrtflx
20051125 - DWK changed public Biomass cropprod[CYCLE] to 
           Biomass cropprod
20051125 - DWK changed public Biomass cropResidue[CYCLE] to 
           Biomass cropResidue
20051125 - DWK changed public Biomass cropResidueFlux[CYCLE] to 
           Biomass cropResidueFlux                             
20051125 - DWK changed public double fertn[CYCLE] to double fertn           
20051125 - DWK changed public Biomass formCropResidue[CYCLE] to 
           Biomass formCropResidue
20051125 - DWK changed public Biomass formPROD10[CYCLE] to 
           Biomass formPROD10
20051125 - DWK changed public Biomass formPROD100[CYCLE] to 
           Biomass formPROD100                      
20051125 - DWK changed public Biomass formTOTPROD[CYCLE] to  
           Biomass formTOTPROD
20051125 - DWK changed public Biomass initPROD1[CYCLE] to 
           Biomass initPROD1 
20051125 - DWK changed public double irrigate[CYCLE] to 
           double irrigate                                
20051125 - DWK changed public double nce[CYCLE] to double nce
20051125 - DWK changed public double nretent[CYCLE] to 
           double nretent
20051125 - DWK changed public double nsretent[CYCLE] to 
           double nsretent
20051125 - DWK changed public double nvretent[CYCLE] to 
           double nvretent
20051125 - DWK changed public Biomass PROD1[CYCLE] to 
           Biomass PROD1
20051125 - DWK changed public Biomass PROD10[CYCLE] to 
           Biomass PROD10
20051125 - DWK changed public Biomass PROD100[CYCLE] to 
           Biomass PROD100
20051125 - DWK changed public Biomass PROD1decay[CYCLE] to 
           Biomass PROD1decay
20051125 - DWK changed public Biomass PROD10decay[CYCLE] to 
           Biomass PROD10decay
20051125 - DWK changed public Biomass PROD100decay[CYCLE] to 
           Biomass PROD100decay
20051125 - DWK changed public Biomass sconvrtflx[CYCLE] to 
           Biomass sconvrtflx                                                       
20051125 - DWK changed public Biomass slash[CYCLE] to 
           Biomass slash
20051125 - DWK changed public Biomass stubble[CYCLE] to 
           Biomass stubble
20051125 - DWK changed public double totec[CYCLE] to 
           double totec                                 
20051125 - DWK changed public double totgc[CYCLE] to 
           double totgc                                 
20051125 - DWK changed public Biomass TOTPROD[CYCLE] to 
           Biomass TOTPROD
20051125 - DWK changed public Biomass TOTPRODdecay[CYCLE] to 
           Biomass TOTPRODdecay
20051125 - DWK changed public Biomass vconvrtflx[CYCLE] to 
           Biomass vconvrtflx 
20051125 - DWK deleted const int& outmon from conversion(),
           decayProducts(), harvest(), updateCropResidue(),
           updateCropResidueFluxes(), updateProducts() and
           updateTotalProductFormation()
20051202 - DWK added public function resetMonthlyFluxes()
20051202 - DWK added public function resetYrFluxes()
20060608 - DWK added public function resetMonthlyDisturbFluxes()
20060608 - DWK added public function setFireNDEP()
20060608 - DWK added private double firendep
20070210 - DWK changed include from tveg602.h to tveg603.h
20070518 - DWK added public function grazing()
20070518 - DWK added public function setNoGrazing()
20070518 - DWK added private Biomass forage
20070518 - DWK added private Biomass manure
20070518 - DWK added private double urine
20070518 - DWK added private double animalresp
20070829 - DWK changed include from tveg603.h to tveg603b.h
20090128 - DWK changed include from tveg603b.h to tveg604.h
										                                   
****************************************************************
************************************************************* */

#ifndef HUMNACT604_H
#define HUMNACT604_H

const double GDDMIN = 5.00000;
const double GDDSEED = 300.0;
const double GDDHARVST = 2000.0;

#include "temconsts602.hpp"
#include "tprocessXML602.h"
#include "bioms423.hpp"  // Humnact60 uses Biomass class
#include "tveg604.h"     // Humanact60 uses Tveg60 class

class Humanact60 : public ProcessXML60
{

  public:

     Humanact60();

/* **************************************************************
		 Public Functions
************************************************************** */

     void conversion( const int& pdcmnt,
                      const double& vegc,
                      const double& leafc,
                      const double& stemc,
                      const double& rootc,
                      const double& deadwoodcarbon,
                      const double& litterc,
                      const double& cwdc,
                      const double& mmdc,
                      const double& deadwoodnitrogen,
                      const double& vstrn,
                      const double& vston, 
                      const double& solc, 
                      const double& soiln,
                      const int& disturbflag,
                      const int& ifwetland,
                      const double& randomnum);

     void createWoodProducts( const int& outyr,
                              const double& vegc,
                              const double& stemc,
                              const double& vstrn,
                              const double& vston );
                      
     void decayProducts( void );

     void frostDamage( const double& vegc, 
                       const double& vstrn,
                       const double& vston );
         
     void getecd( ofstream& rflog1 );
     
     void getecd( const string& ecd );

     void grazing( const double& vegc, const double& strn );

     void harvest( const int& pdm,
                   const double& vegc, 
                   const double& vstrn,
                   const double& vston );
                       
     void resetMonthlyDisturbFluxes( void );

     void resetMonthlyFluxes( void );
         
     void resetPROD( void );

     void resetYrFluxes( void );
     
     void setAgricFlags( ofstream& rflog1 );

     void setFireNDEP( void );

     void setNoCrops( const int& pdm );
     
     void setNoGrazing( void );

     void setNoWoodProducts( const int& pdyr );
     
     void updateCropResidue( void );
     
     void updateCropResidueFluxes( void );
     
     void updateProducts( void );
     
     void updateTotalProductFormation( void );
     double sampleNormal(void);
     double normalDistribution(double x);

     // "Get" and "Set" private variables

     // animalresp *********************************************

     inline double getANIMALRESP( void ) { return animalresp; }

     inline void setANIMALRESP( const double& panimalresp )
     {
       animalresp = panimalresp; }


     // c2n ****************************************************
     
     inline double getC2N( void ) { return c2n; }

     inline void setC2N( const double& pc2n ) { c2n = pc2n; }

     
     // cflux **************************************************
     
     inline double getCFLUX( void ) { return cflux; }
     
     inline void setCFLUX( const double& pcflx ) 
     { 
       cflux = pcflx; 
     }

     
     // convrtflx.carbon ***************************************
     
     inline double getCONVRTFLXC( void ) 
     { 
       return convrtflx.carbon; 
     }

     inline void setCONVRTFLXC( const double& pcnvrtflxc ) 
     { 
       convrtflx.carbon = pcnvrtflxc; 
     }


     // convrtflx.nitrogen *************************************
     
     inline double getCONVRTFLXN( void ) 
     { 
       return convrtflx.nitrogen; 
     }

     inline void setCONVRTFLXN( const double& pcnvrtflxn ) 
     { 
       convrtflx.nitrogen = pcnvrtflxn; 
     }


     // cropprod.carbon ****************************************
     
     inline double getCROPPRODC( void ) 
     { 
       return cropprod.carbon; 
     }

     inline void setCROPPRODC( const double& pcrpprdc ) 
     { 
       cropprod.carbon = pcrpprdc; 
     }


     // cropprod.nitrogen **************************************
     
     inline double getCROPPRODN( void ) 
     { 
       return cropprod.nitrogen; 
     }

     inline void setCROPPRODN( const double& pcrpprdn ) 
     { 
       cropprod.nitrogen = pcrpprdn; 
     }


     // cropprveetmx *******************************************
     
     inline double getCROPPRVEETMX( void ) 
     { 
       return cropprveetmx; 
     }

     inline void setCROPPRVEETMX( const double& pcrpprveetmx ) 
     { 
       cropprveetmx = pcrpprveetmx; 
     }


     // cropprvleafmx ******************************************
     
     inline double getCROPPRVLEAFMX( void ) 
     { 
       return cropprvleafmx; 
     }

     inline void setCROPPRVLEAFMX( const double& pcrpprvleafmx ) 
     { 
       cropprvleafmx = pcrpprvleafmx; 
     }


     // cropprvpetmx *******************************************
     
     inline double getCROPPRVPETMX( void ) 
     { 
       return cropprvpetmx; 
     }

     inline void setCROPPRVPETMX( const double& pcrpprvpetmx ) 
     { 
       cropprvpetmx = pcrpprvpetmx; 
     }


     // cropResidue.carbon *************************************
     
     inline double getCROPRESIDUEC( void ) 
     { 
       return cropResidue.carbon; 
     }

     inline void setCROPRESIDUEC( const double& pcrpresiduec ) 
     { 
       cropResidue.carbon = pcrpresiduec; 
     }


     // cropResidue.nitrogen ***********************************
     
     inline double getCROPRESIDUEN( void ) 
     { 
       return cropResidue.nitrogen; 
     }

     inline void setCROPRESIDUEN( const double& pcrpresiduen ) 
     { 
       cropResidue.nitrogen = pcrpresiduen; 
     }


     // cropResidueFlux.carbon *********************************
     
     inline double getCROPRESIDUEFLXC( void ) 
     { 
       return cropResidueFlux.carbon; 
     }



     // cropResidueFlux.nitrogen *******************************
     
     inline double getCROPRESIDUEFLXN( void ) 
     { 
       return cropResidueFlux.nitrogen; 
     }


     // cropseedC **********************************************
     
     inline double getCROPSEEDC( const int& pcmnt ) 
     { 
       return cropseedC[pcmnt]; 
     }

     inline void setCROPSEEDC( const double& pcropseedc,
                               const int& pcmnt ) 
     { 
       cropseedC[pcmnt] = pcropseedc; 
     }

 
     // cropseedSTON *******************************************
     
     inline double getCROPSEEDSTON( const int& pcmnt ) 
     { 
       return cropseedSTON[pcmnt]; 
     }

     inline void setCROPSEEDSTON( const double& pcropseedston,
                                  const int& pcmnt ) 
     { 
       cropseedSTON[pcmnt] = pcropseedston; 
     }

 
     // cropseedSTRN *******************************************
     
     inline double getCROPSEEDSTRN( const int& pcmnt ) 
     { 
       return cropseedSTRN[pcmnt]; 
     }

     inline void setCROPSEEDSTRN( const double& pcropseedstrn,
                                  const int& pcmnt ) 
     { 
       cropseedSTRN[pcmnt] = pcropseedstrn; 
     }

 
     // croptopt ***********************************************
     
     inline double getCROPTOPT( void ) { return croptopt; }

     inline void setCROPTOPT( const double& pcrptopt ) 
     { 
       croptopt = pcrptopt; 
     }


     // fertn **************************************************
               
//     inline double getFERTN( void ) { return fertn; }

//     inline void setFERTN( const double& pfertn ) 
//     { 
//       fertn = pfertn; 
//     }


     // firendep ***********************************************
     
     inline double getFIRENDEP( void ) 
     { 
       return firendep; 
     }

     inline void setFIRENDEP( const double& pfirendep ) 
     { 
       firendep = pfirendep; 
     }


     // forage.carbon *********************************

     inline double getFORAGEC( void )
     {
       return forage.carbon;
     }

     inline void setFORAGEC( const double& pforagec )
     {
       forage.carbon = pforagec;
     }

     // forage.nitrogen *********************************

     inline double getFORAGEN( void )
     {
       return forage.nitrogen;
     }

     inline void setFORAGEN( const double& pforagen )
     {
       forage.nitrogen = pforagen;
     }


     // formCropResidue.carbon *********************************
     
     inline double getFORMCROPRESIDUEC( void ) 
     { 
       return formCropResidue.carbon; 
     }

     inline void setFORMCROPRESIDUEC( const double& pfrmcrpresiduec ) 
     { 
       formCropResidue.carbon = pfrmcrpresiduec; 
     }


     // formCropResidue.nitrogen *******************************
     
     inline double getFORMCROPRESIDUEN( void ) 
     { 
       return formCropResidue.nitrogen; 
     }

     inline void setFORMCROPRESIDUEN( const double& pfrmcrpresiduen ) 
     { 
       formCropResidue.nitrogen = pfrmcrpresiduen; 
     }


     // formPROD10.carbon **************************************
     
     inline double getFORMPROD10C( void ) 
     { 
       return formPROD10.carbon; 
     }

     inline void setFORMPROD10C( const double& pfrmprod10c ) 
     { 
       formPROD10.carbon = pfrmprod10c; 
     }


     // formPROD10.nitrogen ************************************
     
     inline double getFORMPROD10N( void ) 
     { 
       return formPROD10.nitrogen; 
     }

     inline void setFORMPROD10N( const double& pfrmprod10n ) 
     { 
       formPROD10.nitrogen = pfrmprod10n; 
     }


     // formPROD100.carbon *************************************
     
     inline double getFORMPROD100C( void ) 
     { 
       return formPROD100.carbon; 
     }

     inline void setFORMPROD100C( const double& pfrmprod100c ) 
     { 
       formPROD100.carbon = pfrmprod100c; 
     }


     // formPROD100.nitrogen ***********************************
     
     inline double getFORMPROD100N( void ) 
     { 
       return formPROD100.nitrogen; 
     }

     inline void setFORMPROD100N( const double& pfrmprod100n ) 
     { 
       formPROD100.nitrogen = pfrmprod100n; 
     }


     // formTOTPROD.carbon *************************************
     
     inline double getFORMTOTPRODC( void ) 
     { 
       return formTOTPROD.carbon; 
     }

     inline void setFORMTOTPRODC( const double& pfrmtotprodc ) 
     { 
       formTOTPROD.carbon = pfrmtotprodc; 
     }


     // formTOTPROD.nitrogen ***********************************
     
     inline double getFORMTOTPRODN( void ) 
     { 
       return formTOTPROD.nitrogen; 
     }

     inline void setFORMTOTPRODN( const double& pfrmtotprodn ) 
     { 
       formTOTPROD.nitrogen = pfrmtotprodn; 
     }


     // FRI ****************************************************
     
     inline double getFRI( void ) { return FRI; }

     inline void setFRI( const double& pfri ) { FRI = pfri; }


     // growdd *************************************************
     
     inline double getGROWDD( void ) { return growdd; }

     inline void setGROWDD( const double& pgrowdd ) 
     { 
       growdd = pgrowdd; 
     }


     // harvstC *********************************************
     
     inline double getHARVSTC( const int& pcmnt ) 
     { 
       return harvstC[pcmnt]; 
     }

     inline void setHARVSTC( const double& pharvstC, 
                             const int& pcmnt ) 
     { 
       harvstC[pcmnt] = pharvstC; 
     }


     // harvstN *********************************************
     
     inline double getHARVSTN( const int& pcmnt ) 
     { 
       return harvstN[pcmnt]; 
     }

     inline void setHARVSTN( const double& pharvstN, 
                                const int& pcmnt ) 
     { 
       harvstN[pcmnt] = pharvstN; 
     }


     // initPROD1.carbon ***************************************
     
     inline double getINITPROD1C( const int& pdm ) 
     { 
       return initPROD1[pdm].carbon; 
     }

     inline void setINITPROD1C( const double& pinprod1c,
                                const int& pdm ) 
     { 
       initPROD1[pdm].carbon = pinprod1c; 
     }


     // initPROD1.nitrogen *************************************
     
     inline double getINITPROD1N( const int& pdm ) 
     { 
       return initPROD1[pdm].nitrogen; 
     }

     inline void setINITPROD1N( const double& pinprod1n,
                                const int& pdm ) 
     { 
       initPROD1[pdm].nitrogen = pinprod1n; 
     }


     // initPROD10[].carbon ************************************
     
     inline double getINITPROD10C( const int& i ) 
     { 
       return initPROD10[i].carbon; 
     }

     inline void setINITPROD10C( const double& pinprod10c, 
                                 const int& i ) 
     { 
       initPROD10[i].carbon = pinprod10c; 
     }


     // initPROD10[].nitrogen **********************************
     
     inline double getINITPROD10N( const int& i ) 
     { 
       return initPROD10[i].nitrogen; 
     }

     inline void setINITPROD10N( const double& pinprod10n, 
                                 const int& i ) 
     { 
       initPROD10[i].nitrogen = pinprod10n; 
     }


     // initPROD100[].carbon ***********************************
      
     inline double getINITPROD100C( const int& i ) 
     { 
       return initPROD100[i].carbon; 
     }

     inline void setINITPROD100C( const double& pinprod100c, 
                                  const int& i ) 
     { 
       initPROD100[i].carbon = pinprod100c; 
     }


     // initPROD100[].nitrogen *********************************
     
     inline double getINITPROD100N( const int& i ) 
     { 
       return initPROD100[i].nitrogen; 
     }

     inline void setINITPROD100N( const double& pinprod100n, 
                                  const int& i ) 
     { 
       initPROD100[i].nitrogen = pinprod100n; 
     }


     // irrigate ***********************************************
     
//     inline double getIRRIGATE( void ) { return irrigate; }


     // kd *****************************************************
     
     inline double getKD( void ) { return kd; }

     inline void setKD( const double& pkd ) { kd = pkd; }


     //manure.carbon *********************************

     inline double getMANUREC( void )
     {
       return manure.carbon;
     }

     inline void setMANUREC( const double& pmanurec )
     {
       manure.carbon = pmanurec;
     }

     // manure.nitrogen *********************************

     inline double getMANUREN( void )
     {
       return manure.nitrogen;
     }

     inline void setMANUREN( const double& pmanuren )
     {
       manure.nitrogen = pmanuren;
     }


     // natprveetmx ********************************************
     
     inline double getNATPRVEETMX( void ) 
     { 
       return natprveetmx; 
     }

     inline void setNATPRVEETMX( const double& pnatprveetmx ) 
     { 
       natprveetmx = pnatprveetmx; 
     }


     // natprvleafmx *******************************************
     
     inline double getNATPRVLEAFMX( void ) 
     { 
       return natprvleafmx; 
     }

     inline void setNATPRVLEAFMX( const double& pnatprvleafmx ) 
     { 
       natprvleafmx = pnatprvleafmx; 
     }


     // natprvpetmx ********************************************
     
     inline double getNATPRVPETMX( void ) 
     { 
       return natprvpetmx; 
     }

     inline void setNATPRVPETMX( const double& pnatprvpetmx ) 
     { 
       natprvpetmx = pnatprvpetmx; 
     }


     // natseedC ***********************************************
     
     inline double getNATSEEDC( void ) { return natseedC; }

     inline void setNATSEEDC( const double& pnatseedc ) 
     { 
       natseedC = pnatseedc; 
     }


     // natseedSTRN ********************************************
     
     inline double getNATSEEDSTRN( void ) 
     { 
       return natseedSTRN; 
     }

     inline void setNATSEEDSTRN( const double& pnatseedstrn ) 
     { 
       natseedSTRN = pnatseedstrn; 
     }


     // natseedSTON ********************************************
     
     inline double getNATSEEDSTON( void ) 
     { 
       return natseedSTON; 
     }

     inline void setNATSEEDSTON( const double& pnatseedston ) 
     { 
       natseedSTON = pnatseedston; 
     }


     // natsoil ************************************************
     
     inline double getNATSOIL( void ) { return natsoil; }

     inline void setNATSOIL( const double& pnatsoil ) 
     { 
       natsoil = pnatsoil; 
     }


     // nattopt ************************************************
     
     inline double getNATTOPT( void ) { return nattopt; }

     inline void setNATTOPT( const double& pnattopt ) 
     { 
       nattopt = pnattopt; 
     }


     // natyreet ***********************************************

     inline double getNATYREET( void ) { return natyreet; }

     inline void setNATYREET( const double& pnatyreet ) 
     { 
       natyreet = pnatyreet; 
     }


     // natyrpet ***********************************************

     inline double getNATYRPET( void ) { return natyrpet; }

     inline void setNATYRPET( const double& pnatyrpet ) 
     { 
       natyrpet = pnatyrpet; 
     }

     
     // nretent ************************************************
     
     inline double getNRETENT( void ) { return nretent; }

     inline void setNRETENT( const double& pnretent ) 
     { 
       nretent = pnretent; 
     }


     // nsretconv **********************************************
     
     inline double getNSRETCONV( const int& pcmnt ) 
     { 
       return nsretconv[pcmnt]; 
     }

     inline void setNSRETCONV( const double& pnsretconv, 
                               const int& pcmnt ) 
     { 
       nsretconv[pcmnt] = pnsretconv; 
     }


     // nsretent ***********************************************
     
     inline double getNSRETENT( void ) { return nsretent; }

     inline void setNSRETENT( const double& pnsretent ) 
     { 
       nsretent = pnsretent; 
     }

 
     // nvretconv **********************************************
     
     inline double getNVRETCONV( const int& pcmnt ) 
     { 
       return nvretconv[pcmnt]; 
     }

     inline void setNVRETCONV( const double& pnvretconv, 
                               const int& pcmnt ) 
     { 
       nvretconv[pcmnt] = pnvretconv; 
     }


     // nvretent ***********************************************
     
     inline double getNVRETENT( void ) { return nvretent; }

     inline void setNVRETENT( const double& pnvretent ) 
     { 
       nvretent = pnvretent; 
     }


     // prevCropResidue.carbon *********************************
     
     inline double getPREVCROPRESIDUEC( void ) 
     { 
       return prevCropResidue.carbon; 
     }

     inline void setPREVCROPRESIDUEC( const double& pcrpresc ) 
     { 
       prevCropResidue.carbon = pcrpresc; 
     }


     // prevCropResidue.nitrogen *******************************
     
     inline double getPREVCROPRESIDUEN( void ) 
     { 
       return prevCropResidue.nitrogen; 
     }

     inline void setPREVCROPRESIDUEN( const double& pcrpresn ) 
     { 
       prevCropResidue.nitrogen = pcrpresn; 
     }


     // prevPROD1.carbon ***************************************
     
     inline double getPREVPROD1C( void ) 
     { 
       return prevPROD1.carbon; 
     }

     inline void setPREVPROD1C( const double& pprod1c ) 
     { 
       prevPROD1.carbon = pprod1c; 
     }


     // prev.PROD1.nitrogen ************************************
     
     inline double getPREVPROD1N( void ) 
     { 
       return prevPROD1.nitrogen; 
     }

     inline void setPREVPROD1N( const double& pprod1n ) 
     { 
       prevPROD1.nitrogen = pprod1n; 
     }


     // prevPROD10.carbon **************************************
     
     inline double getPREVPROD10C( void ) 
     { 
       return prevPROD10.carbon; 
     }

     inline void setPREVPROD10C( const double& pprod10c ) 
     { 
       prevPROD10.carbon = pprod10c; 
     }


     // prevPROD10.nitrogen ************************************
     
     inline double getPREVPROD10N( void ) 
     { 
       return prevPROD10.nitrogen; 
     }

     inline void setPREVPROD10N( const double& pprod10n ) 
     { 
       prevPROD10.nitrogen = pprod10n; 
     }


     // prevPROD100.carbon *************************************
     
     inline double getPREVPROD100C( void ) 
     { 
       return prevPROD100.carbon; 
     }

     inline void setPREVPROD100C( const double& pprod100c ) 
     { 
       prevPROD100.carbon = pprod100c; 
     }


     // prevPROD100.nitrogen ***********************************
     
     inline double getPREVPROD100N( void ) 
     { 
       return prevPROD100.nitrogen; 
     }

     inline void setPREVPROD100N( const double& pprod100n ) 
     { 
       prevPROD100.nitrogen = pprod100n; 
     }


     // PROD1.carbon *******************************************
     
     inline double getPROD1C( void ) { return PROD1.carbon; }


     // PROD1.nitrogen *****************************************
     
     inline double getPROD1N( void ) { return PROD1.nitrogen; }


     // PROD10.carbon ******************************************
     
     inline double getPROD10C( void ) { return PROD10.carbon; }


     // PROD10.nitrogen ****************************************
     
     inline double getPROD10N( void ) 
     { 
       return PROD10.nitrogen; 
     }


     // PROD100.carbon *****************************************
     
     inline double getPROD100C( void ) 
     { 
       return PROD100.carbon; 
     }


     // PROD100.nitrogen ***************************************
     
     inline double getPROD100N( void ) 
     { 
       return PROD100.nitrogen; 
     }


     // PROD1decay.carbon **************************************
     
     inline double getPROD1DECAYC( void ) 
     { 
       return PROD1decay.carbon; 
     }


     // PROD1decay.nitrogen ************************************
     
     inline double getPROD1DECAYN( void ) 
     { 
       return PROD1decay.nitrogen; 
     }


     // PROD10decay.carbon *************************************
     
     inline double getPROD10DECAYC( void ) 
     { 
       return PROD10decay.carbon; 
     }


     // PROD10decay.nitrogen ***********************************
     
     inline double getPROD10DECAYN( void ) 
     { 
       return PROD10decay.nitrogen; 
     }


     // PROD100decay.carbon ************************************
     
     inline double getPROD100DECAYC( void ) 
     { 
       return PROD100decay.carbon; 
     }


     // PROD100decay.nitrogen **********************************
     
     inline double getPROD100DECAYN( void ) 
     { 
       return PROD100decay.nitrogen; 
     }

     /*
     // prod10par **********************************************
     
     inline double getPROD10PAR( void ) { return prod10par; }

     inline void setPROD10PAR( const double& pprod10par ) 
     { 
       prod10par = pprod10par; 
     }


     // prod100par *********************************************
     
     inline double getPROD100PAR( void ) { return prod100par; }

     inline void setPROD100PAR( const double& pprod100par ) 
     { 
       prod100par = pprod100par; 
     }

	*/
     // productYear ********************************************
     
     inline int getPRODUCTYEAR( void ) { return productYear; }

     inline void setPRODUCTYEAR( const int& pprodyr ) 
     { 
       productYear = pprodyr; 
     }


     // prvcropnpp *********************************************
     
     inline double getPRVCROPNPP( void ) { return prvcropnpp; }

     inline void setPRVCROPNPP( const double& pprvcropnpp ) 
     { 
       prvcropnpp = pprvcropnpp; 
     }


     // residueC ***********************************************
     
     inline double getRESIDUEC( const int& pcmnt ) 
     { 
       return residueC[pcmnt]; 
     }

     inline void setRESIDUEC( const double& presidueC, 
                              const int& pcmnt ) 
     { 
       residueC[pcmnt] = presidueC; 
     }


     // residueN ***********************************************
     
     inline double getRESIDUEN( const int& pcmnt ) 
     { 
       return residueN[pcmnt]; 
     }

     inline void setRESIDUEN( const double& presidueN, 
                              const int& pcmnt ) 
     { 
       residueN[pcmnt] = presidueN; 
     }


     // sconvert ***********************************************
     
     inline double getSCONVERT( void ) { return sconvert; }

     inline void setSCONVERT( const double& pscnvrt ) 
     { 
       sconvert = pscnvrt; 
     }


     // sconvrtflx.carbon **************************************
     
     inline double getSCONVRTFLXC( void ) 
     { 
       return sconvrtflx.carbon; 
     }

     inline void setSCONVRTFLXC( const double& pscnvrtflxc ) 
     { 
       sconvrtflx.carbon = pscnvrtflxc; 
     }


     // sconvrtflx.nitrogen ************************************
     
     inline double getSCONVRTFLXN( void ) 
     { 
       return sconvrtflx.nitrogen; 
     }

     inline void setSCONVRTFLXN( const double& pscnvrtflxn ) 
     { 
       sconvrtflx.nitrogen = pscnvrtflxn; 
     }


     // slash.carbon *******************************************
     
     inline double getSLASHC( void ) { return slash.carbon; }

     inline void setSLASHC( const double& pslashc ) 
     { 
       slash.carbon = pslashc; 
     }


     // slash.nitrogen *****************************************

     inline double getSLASHN( void ) { return slash.nitrogen; }

     inline void setSLASHN( const double& pslashn ) 
     { 
       slash.nitrogen = pslashn; 
     }


     // slashpar ***********************************************
     
     inline double getSLASHPAR( void ) { return slashpar; }

     inline void setSLASHPAR( const double& pslashpar ) 
     { 
       slashpar = pslashpar; 
     }


     // stubble.carbon *****************************************
     
     inline double getSTUBBLEC( void ) 
     { 
       return stubble.carbon; 
     }

     inline void setSTUBBLEC( const double& pstubc ) 
     { 
       stubble.carbon = pstubc; 
     }


     // stubble.nitrogen ***************************************
     
     inline double getSTUBBLEN( void ) 
     { 
       return stubble.nitrogen; 
     }

     inline void setSTUBBLEN( const double& pstubn ) 
     { 
       stubble.nitrogen = pstubn; 
     }


     // tillfactor *********************************************
     
     inline double getTILLFACTOR( const int& pcmnt ) 
     { 
       return tillfactor[pcmnt]; 
     }

     inline void setTILLFACTOR( const double& ptillfact, 
                                const int& pcmnt ) 
     { 
       tillfactor[pcmnt] = ptillfact; 
     }


     // totec **************************************************
     
     inline double getTOTEC( void ) { return totec; }

     inline void setTOTEC( const double& ptotec ) 
     { 
       totec = ptotec; 
     }


     // TOTPROD.carbon *****************************************
     
     inline double getTOTPRODC( void ) 
     { 
       return TOTPROD.carbon; 
     }


     // TOTPROD.nitrogen ***************************************
     
     inline double getTOTPRODN( void ) 
     { 
       return TOTPROD.nitrogen; 
     }


     // TOTPRODdecay.carbon ************************************
     
     inline double getTOTPRODDECAYC( void ) 
     { 
       return TOTPRODdecay.carbon; 
     }


     // TOTPRODdecay.nitrogen **********************************
     
     inline double getTOTPRODDECAYN( void ) 
     { 
       return TOTPRODdecay.nitrogen; 
     }


     // urine **************************************************

     inline double getURINE( void )
     {
       return urine;
     }

     inline void setURINE( const double& purine )
     {
       urine = purine;
     }


     // vconvert ***********************************************
     
     inline double getVCONVERT( void ) { return vconvert; }

     inline void setVCONVERT( const double& pvcnvrt ) 
     { 
       vconvert = pvcnvrt; 
     }


     // vconvrtflx.carbon **************************************
     
     inline double getVCONVRTFLXC( void ) 
     { 
       return vconvrtflx.carbon; 
     }

     inline void setVCONVRTFLXC( const double& pvcnvrtflxc ) 
     { 
       vconvrtflx.carbon = pvcnvrtflxc; 
     }


     // vconvrtflx.nitrogen ************************************
     
     inline double getVCONVRTFLXN( void ) 
     { 
       return vconvrtflx.nitrogen; 
     }

     inline void setVCONVRTFLXN( const double& pvcnvrtflxn ) 
     { 
       vconvrtflx.nitrogen = pvcnvrtflxn; 
     }


     // vrespar ************************************************
     
     inline double getVRESPAR( void ) { return vrespar; }

     inline void setVRESPAR( const double& pvrespar ) 
     { 
       vrespar = pvrespar; 
     }
     inline double getdeadwoodc( void )
         {
           return deadwood.carbon;
         }

         inline void setdeadwoodc( const double& pdeadwoodc )
         {
           deadwood.carbon = pdeadwoodc;
         }
         inline double getdeadwoodn( void )
         {
           return deadwood.nitrogen;
         }

         inline void setdeadwoodn( const double& pdeadwoodn )
         {
           deadwood.nitrogen = pdeadwoodn;
         }
/*
     inline double getdeadconvert( void ) { return deadconvert; }

     inline void setdeadconvert( const double& pdeadconv )
     {
       deadconvert = pdeadconv;
     }

     inline double getdeadslash( void ) { return deadslash; }

     inline void setdeadslash( const double& pdeadslash )
     {
       deadslash = pdeadslash;
     }*/
     inline double getfflc( void )
     {
       return FFLC;
     }
     inline double getffdc( void )
     {
       return FFDC;
     }
     inline double getswfc( void )
     {
       return SWFC;
     }

     inline double getdwfc( void )
     {
       return DWFC;
     }
	 
	 inline double getdeadconvc( void )
     {
       return convertdeadc;
     }
	 	 inline double getdeadconvn( void )
     {
       return convertdeadn;
     }
     inline double getdeadslashc( void )
     {
       return deadslashc;
     }
     inline double getdeadslashn( void )
     {
       return deadslashn;
     }
     inline double getstem2slash( void )
     {
       return stem2slash;
     }

     inline double getleaf2slash( void )
     {
       return leaf2slash;
     }
     inline double getroot2slash( void )
     {
       return root2slash;
     }
     inline double getage( void )
     {
       return double(forestage);
     }

/* *************************************************************
		 Public Variables
************************************************************* */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif

     // Index for disturbed community type
    int cmnt;

     // Flag to indicate start of optimum fertilization of 
     //   croplands after 1950
//    int fert1950flag;
    
     // Flag to indicate whether or not optimum
     //   fertilization of croplands occurs
    int fertflag;

    // Monthly fertilizer application 
    double fertn;

    // Flag to indicate whether or not optimum irrigation of
    //   croplands occur
    int irrgflag;

    // Monthly irrigation (mm / month)
    double irrigate;

    int lulcyear;

    int massbalflg;

    int mez;

    // Flag to indicate that element is in agriculture during
    //   previous year
    int prvstate;

    // Flag to indicate that element is in agriculture during
    //   current year
    // state = 1 - row-crop agriculture
    // state = 2 - pasture
    // state = 3 - urban area
    int state;

     // Flag to indicate whether or not croplands are tilled
    int tillflag;

     // Flag to indicate whether or not simulation is conducted 
     //   with transient land use data 
    int tlulcflag;

    // Annual sum of cflux
    double yrcflux;

    // Annual sum of convrtflx.carbon
    double yrconvrtC;

    // Annual sum of convrtflx.nitrogen
    double yrconvrtN;

    // Annual sum of PROD1decay.carbon
    double yrdecayPROD1C;

    // Annual sum of PROD10decay.carbon
    double yrdecayPROD10C;

    // Annual sum of PROD100decay.carbon
    double yrdecayPROD100C;

    // Annual sum of PROD1decay.nitrogen
    double yrdecayPROD1N;

    // Annual sum of PROD10decay.nitrogen
    double yrdecayPROD10N;

    // Annual sum of PROD100decay.nitrogen
    double yrdecayPROD100N;

    // Annual sum of TOTPRODdecay.carbon
    double yrdecayTOTPRODC;

    // Annual sum of TOTPRODdecay.nitrogen
    double yrdecayTOTPRODN;

    // Annual sum of fertn
    double yrfertn;

    // Annual sum of cropResidueFlux.carbon
    double yrfluxResidueC;

    // Annual sum of cropResidueFlux.nitrogen
    double yrfluxResidueN;

    // Annual sum of cropprod.carbon
    double yrformPROD1C;

    // Annual sum of formPROD10.carbon
    double yrformPROD10C;

    // Annual sum of formPROD100.carbon
    double yrformPROD100C;

    // Annual sum of cropprod.nitrogen
    double yrformPROD1N;

    // Annual sum of formPROD10.nitrogen
    double yrformPROD10N;

    // Annual sum of formPROD100.nitrogen
    double yrformPROD100N;

    // Annual sum of formCropResidue.carbon
    double yrformResidueC;

    // Annual sum of formCropResidue.nitrogen
    double yrformResidueN;

    // Annual sum of formTOTPROD.carbon
    double yrformTOTPRODC;

    // Annual sum of formTOTPROD.nitrogen
    double yrformTOTPRODN;

    // Annual sum of irrigate
    double yrirrig;

    // Annual sum of nretent
    double yrnrent;

    // Annual sum of nsretent
    double yrnsrent;

    // Annual sum of nvretent
    double yrnvrent;

    // Annual sum of sconvrtflx.carbon
    double yrsconvrtC;

    // Annual sum of sconvrtflx.nitrogen
    double yrsconvrtN;

    // Annual sum of slash.carbon
    double yrslashC;

    // Annual sum of slash.nitrogen
    double yrslashN;

    // Annual sum of stubble.carbon
    double yrstubC;

    // Annual sum of stubble.nitrogen
    double yrstubN;

    // Annual sum of vconvrtflx.carbon
    double yrvconvrtC;

    // Annual sum of vconvrtflx.nitrogen
    double yrvconvrtN;
    double deadconvert;

    double leafmortpar;
    double stemmortpar;
    double rootmortpar;
    double leafslash;
    double leafconv;
    double rootslash;
    double rootconv;
    double stemslash;
    double stemconv;
    double standdead;
    double deadconv;
    double deadslash;
    double sconvert;
    double prod10par;
    double prod100par;
    int forestage;

   private:
   
/* *************************************************************
		 Private Variables
************************************************************* */

     // Monthly carbon flux to the atmosphere as a result of
     // livestock respiration
     double animalresp;

     // Monthly carbon flux to the atmosphere as a result of
     //   NEP plus agricultural conversion fluxes
     double cflux;

     // Monthly loss of carbon and nitrogen to the atmosphere as a
     //   result of agricultural conversion
     Biomass convrtflx;

     // Monthly carbon and nitrogen fluxes associated with crop 
     //   yield
     Biomass cropprod;

     // Maximum monthly estimated evapotranspiration of crops 
     //   during the previous year 
     double cropprveetmx;

     // Maximum monthly relative leaf area of crops during the
     //   previous year 
     double cropprvleafmx;

     // Maximum monthly potential evapotranspiration of crops 
     //   during the previous year 
     double cropprvpetmx;

     // Stock of crop residue left in crop fields after harvest
     Biomass cropResidue;

     // Monthly loss of carbon from burning crop residue
     Biomass cropResidueFlux;

     double croptopt;

     // Monthly reintroduction of N volatilized in fires
     //   back to the soil (as N deposition)
     double firendep;
    
     // Monthly forage rates
     Biomass forage;

     // Monthly creation of crop residue from harvest
     Biomass formCropResidue;

     // Monthly formation of 10-year wood products
     Biomass formPROD10;

     // Monthly formation of 100-year wood products
     Biomass formPROD100;

     // Monthly formation of all "products"
     Biomass formTOTPROD;

     // Fire return interval
     double FRI;

     // Flag to indicate frost damage
     int frostflag;

     // Growing degree days ( degrees C )
     double growdd;
    
     Biomass initPROD1[CYCLE];
    
     Biomass initPROD10[10];
    
     Biomass initPROD100[100];


     // Intrinsic decomposition rate
     double kd;

     // Monthly manure inputs
     Biomass manure;

     // Maximum monthly estimated evapotranspiration of natural 
     //   vegetation during the previous year 
     double natprveetmx;

     // Maximum monthly relative leaf area of natural vegetation 
     //   during the previous year 
     double natprvleafmx;

     // Maximum monthly potential evapotranspiration of natural
     //   vegetation during the previous year 
     double natprvpetmx;

     double natsoil;

     double natseedC;

     double natseedSTRN;

     double natseedSTON;

     double nattopt;

     double natyreet;

     double natyrpet;
   
     // Total nitrogen retained by ecosystem after conversion
     double nretent;

     // Nitrogen from soil organic matter retained after 
     //   conversion
     double nsretent;

     // Nitrogen from vegetation biomass retained after 
     //   conversion
     double nvretent;

     // Carbon and nitrogen stocks in agricultural products
     //   during previous month
     Biomass prevPROD1;
    
     // Carbon and nitrogen stocks in 10-year wood products
     //   during previous month
     Biomass prevPROD10;

     // Carbon and nitrogen stocks in 100-year wood products
     //   during previous month
     Biomass prevPROD100;

     // Carbon and nitrogen stocks in agricultural products
     //   during current month
     Biomass PROD1;

     // Carbon and nitrogen stocks in 10-year wood products
     //   during current month
     Biomass PROD10;

     // Carbon and nitrogen stocks in 100-year wood products
     //   during current month
     Biomass PROD100;

     // Loss of carbon and nitrogen from "decay" of agricultural
     //   products
     Biomass PROD1decay;

     // Loss of carbon and nitrogen from "decay" of 10-year wood
     //   products
     Biomass PROD10decay;

     // Loss of carbon and nitrogen from "decay" of 100-year wood
     //   products
     Biomass PROD100decay;

     int productYear;

     // Crop net primary production during previous month
     double prvcropnpp;

     // Stock of crop residue from previous month
     Biomass prevCropResidue;

     // Soil organic matter lost during conversion to agriculture
     Biomass sconvrtflx;

     // Stock of slash resulting from conversion to agriculture
     Biomass slash;

     //stock of deadwood from disturbance
     Biomass deadwood;

     // Stock of crop residue that is actually added to soil
     //   (i.e. residue that is not burned off after harvest)
     Biomass stubble;

     // Total carbon storage in ecosystem (i.e. VEGC + TOTSOLC,
     //   does not include carbon in agricultural and wood
     //   products
     double totec;

     // Carbon and nitrogen stocks in all agricultural and wood 
     //   products
     Biomass TOTPROD;

     // Loss of carbon and nitrogen from "decay" of all 
     //   agricultural and wood products
     Biomass TOTPRODdecay;

     // Monthly urine inputs
     double urine;

     // Vegetation biomass lost during conversion to agriculture
     Biomass vconvrtflx;
     double omitzero(const double& data1);
     double stem2slash;
     double leaf2slash;
     double root2slash;

/* *************************************************************
                   Private Parameters
************************************************************* */

     double c2n;

     double cropseedC[MAXCMNT];
     double cropseedSTRN[MAXCMNT];
     double cropseedSTON[MAXCMNT];

     double harvstC[MAXCMNT];
     double harvstN[MAXCMNT];

     double nsretconv[MAXCMNT];
     double nvretconv[MAXCMNT];


//     double RAP;

     double residueC[MAXCMNT];
     double residueN[MAXCMNT];


     double slashpar;

     double tillfactor[MAXCMNT];

     double vconvert;

     double vrespar;
     double FFLC; //forest floor litter C
     double FFDC; //forest floor duff C
     double DWDFC;
     double DWFC; //downed wood (CWD) fuel c
     double SWFC; //standing wood fuel c
     double deadslashc;
     double deadslashn;
     double convertdeadc;
     double convertdeadn;

};

#endif

