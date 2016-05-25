/* *************************************************************
****************************************************************
TVEG603.CPP - object describing characteristics of vegetation 
              used in the Terrestrial Ecosystem Model (TEM)

Modifications:

20030718 - DWK created by modifying tveg50b1.cpp to incorporate
           open-N cycle dynamics into TEM
20030718 - DWK changed Tveg50:: to Tveg51::
20030718 - DWK added nfix[], snfix[], lnfix[], yrinnfix, yrnfix,
           yrsnfix, yrlnfix, innuptake, nuptake, yrinnup, yrnup,
           nh4cut[], nupnh41a[], nupnh41b[], nupnh42a[], nupnh42b[],
           no3cut[], nupno31a[], nupno31b[], nupno32a[], nupno32b[],
           kn1nh4[] and kn1no3[] to Tveg51()
20030718 - DWK changed gppio3() to gppxio3()
20030718 - DWK added nupnh4xclm() and nupno3xclm()
20031014 - BSF added const double& eetpet to gppxo3()
20031016 - DWK replaced char ecd[MAXFNAME] with string ecd to
           getecd() and getleafecd()
20031019 - DWK changed inheritance from ProcessXML to ProcessXML51
20040707 - DWK changed Tveg51:: to Tveg60::
20040707 - DWK changed inheritance from ProcessXML51 to 
           ProcessXML60 in Tveg60()
20040714 - DWK incorporated EE changes to the calculation of
           thawpct in setThawPercent()
20051117 - DWK added include tveg602.h
20051130 - DWK added updateDynamics()
20051201 - DWK added setGV()
20051201 - DWK added setTEMP()
20051202 - DWK added setRESPQ10()
20051202 - DWK added resetMonthlyFluxes()
20051207 - DWK added resetECD()
20060608 - DWK added double extraNH4 to updateDynamics()
20070210 - DWK changed include tveg602.h to tveg603.h
20070210 - DWK added const int& perennial to updateDynamics()
20090128 - DWK added setNewRESPQ10()
20090128 - DWK changed algorithm in rmxclm()
20090128 - DWK modified getecd()
20090128 - DWK changed include from tveg603b.h to tveg604.h
           
****************************************************************
************************************************************* */

#include<cstdio>

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

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std::atoi;

#include<cmath>

  using std::exp;
  using std::pow;
    
#include<string>
  
  using std::string;

#include "tveg604.h"


Tveg60::Tveg60() : ProcessXML60()
{
  int dcmnt;

  plant.carbon = MISSING;
  plant.nitrogen = MISSING;

  strctrl.carbon = MISSING;
  strctrl.nitrogen = MISSING;

  labile.carbon = MISSING;
  labile.nitrogen = MISSING;

  ltrfal.carbon = MISSING;
  ltrfal.nitrogen = MISSING;

  nmobil = MISSING;
  nresorb = MISSING;

  inuptake.total = MISSING;
  inuptake.nh4 = MISSING;
  inuptake.no3 = MISSING;

  nuptake.total = MISSING;
  nuptake.nh4 = MISSING;
  nuptake.no3 = MISSING;

  suptake = MISSING;
  luptake = MISSING;

  ingpp = MISSING;
  gpp = MISSING;

  innpp = MISSING;
  npp = MISSING;

  rm = MISSING;
  rg = MISSING;
  gpr = MISSING;
  rootResp = MISSING;
  npp_leaf = MISSING;
  npp_stem = MISSING;
  npp_croot = MISSING;
  npp_froot = MISSING;

  nfix = MISSING;
  lnfix = MISSING;

  unnormleaf = MISSING;
  leaf = MISSING;

  lai = MISSING;
  fpc = MISSING;


  inprodcn = MISSING;

  yrcarbon = MISSING;
  yrnitrogen = MISSING;

  yrstructn = MISSING;
  yrc2n = MISSING;

  yrstoren = MISSING;

  yrltrc = MISSING;
  yrltrn = MISSING;

  yringpp = MISSING;
  yrgpp = MISSING;

  yrinnpp = MISSING;
  yrnpp = MISSING;

  yrgpr = MISSING;
  yrrmaint = MISSING;
  yrrgrowth = MISSING;

  yrinnfix = MISSING;
  yrnfix = MISSING;
  yrsnfix = MISSING;
  yrlnfix = MISSING;

  yrinnup.total = MISSING;
  yrinnup.nh4 = MISSING;
  yrinnup.no3 = MISSING;

  yrnup.total = MISSING;
  yrnup.nh4 = MISSING;
  yrnup.no3 = MISSING;

  yrsup = MISSING;
  yrlup = MISSING;

  yrnmobil = MISSING;
  yrnrsorb = MISSING;


  yrinpr = MISSING;
  yrprod = MISSING;

  yrunleaf = MISSING;
  yrleaf = MISSING;

  yrfpc = MISSING;
  alleaf = MISSING;
  foliage = MISSING;

  leafyrs = 10;
  deadwoodltc= MISSING;
  deadwoodltn=MISSING;
  yrdeadwoodc= MISSING;
  yrdeadwoodn=MISSING;

  for( dcmnt = 0; dcmnt < MAXCMNT; ++dcmnt )
  {
    c2na[dcmnt] = MISSING;
    c2nb[dcmnt] = MISSING;
    c2nmin[dcmnt] = MISSING;
    cnmin[dcmnt] = MISSING;
    initcneven[dcmnt] = MISSING;

    aleaf[dcmnt] = MISSING;
    bleaf[dcmnt] = MISSING;
    cleaf[dcmnt] = MISSING;
    initleafmx[dcmnt] = MISSING;
    minleaf[dcmnt] = MISSING;
    unleaf12[dcmnt] = MISSING;

    cov[dcmnt] = MISSING;
    fpcmax[dcmnt] = MISSING;
    sla[dcmnt] = MISSING;

    kleafc[dcmnt] = MISSING;
    leafmxc[dcmnt] = MISSING;

    cmaxcut[dcmnt] = MISSING;
    cmax1a[dcmnt] = MISSING;
    cmax1b[dcmnt] = MISSING;
    cmax2a[dcmnt] = MISSING;
    cmax2b[dcmnt] = MISSING;

    kc[dcmnt] = MISSING;
    ki[dcmnt] = MISSING;

    tmax[dcmnt] = MISSING;
    tmin[dcmnt] = MISSING;
    toptmin[dcmnt] = MISSING;
    toptmax[dcmnt] = MISSING;

    gva[dcmnt] = MISSING;

//    kra[dcmnt] = MISSING;
//    krb[dcmnt] = MISSING;

//    raq10a0[dcmnt] = MISSING;
//    raq10a1[dcmnt] = MISSING;
//    raq10a2[dcmnt] = MISSING;
//    raq10a3[dcmnt] = MISSING;

    tref[dcmnt] = MISSING;
    qref[dcmnt] = MISSING;
    alpha[dcmnt] = MISSING;
    beta[dcmnt] = MISSING;
    gamma[dcmnt] = MISSING;
    
    nmaxcut[dcmnt] = MISSING;
    nmax1a[dcmnt] = MISSING;
    nmax1b[dcmnt] = MISSING;
    nmax2a[dcmnt] = MISSING;
    nmax2b[dcmnt] = MISSING;

    nupnh4cut[dcmnt] = MISSING;
    nupnh41a[dcmnt] = MISSING;
    nupnh41b[dcmnt] = MISSING;
    nupnh42a[dcmnt] = MISSING;
    nupnh42b[dcmnt] = MISSING;

    nupno3cut[dcmnt] = MISSING;
    nupno31a[dcmnt] = MISSING;
    nupno31b[dcmnt] = MISSING;
    nupno32a[dcmnt] = MISSING;
    nupno32b[dcmnt] = MISSING;

      kn1[dcmnt] = MISSING;
    kvnh4[dcmnt] = MISSING;
    kvno3[dcmnt] = MISSING;

    cfall[dcmnt] = MISSING;
    nfall[dcmnt] = MISSING;

    labncon[dcmnt] = MISSING;
    //added by cgs
    rsratio[dcmnt] = MISSING;
    ifwoody[dcmnt] = 0;
    lfratiomax[dcmnt] = MISSING;
    woodfall[dcmnt] = MISSING;
    standdeadpar[dcmnt] = MISSING;
  }

  cneven = MISSING;
  c2n = MISSING;
  dc2n = MISSING;
  adjc2n = MISSING;

  newleafmx = MISSING;
  prvleafmx = MISSING;

  cmax = MISSING;

  topt = MISSING;
  newtopt = MISSING;

  kr = MISSING;

  nmax = MISSING;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::boundTOPT( const int& pcmnt )
{
  if( topt > toptmax[pcmnt] ) {	topt = toptmax[pcmnt]; }

  if( topt < toptmin[pcmnt] ) { topt = toptmin[pcmnt]; }
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::deltaleaf( const int& pdcmnt,
                          const double& eet,
                          const double& prveetmx,
                          const double& prvleaf )
{

  double normeet;
  double unnormleaf;
  double maxeet;

  if( prveetmx <= ZERO ) { maxeet = 1.0; }
  else { maxeet = prveetmx; }
  
  normeet = eet / maxeet;
  unnormleaf = (aleaf[pdcmnt] * normeet)
               + (bleaf[pdcmnt] * prvleaf)
               + cleaf[pdcmnt];

  if( unnormleaf < (0.5 * minleaf[pdcmnt]) )
  {
    unnormleaf = 0.5 * minleaf[pdcmnt];
  }

  return unnormleaf;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::getecd( ofstream& rflog1 )
{

  string ecd;

#ifdef PMODE
  *fgo >> ecd;

#else
  cout << "Enter name of the file with the vegetation";
  cout << " parameter values (.ECD):";
  cout << endl;

  cin >> ecd;
#endif

  rflog1 << "Enter name of the file with the vegetation";
  rflog1 << " parameter values (.ECD):";
  rflog1 << ecd << endl << endl;

  getecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::getecd( const string& ecd )
{
  int comtype;
  int dcmnt;
  ifstream infile;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for vegetation ECD input" << endl;
    
    exit( -1 );
  }

  getXMLrootNode( infile, "vegECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "vegECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than ";
      cerr << (MAXCMNT-1) << endl;
      cerr << " in vegECD" << endl;
      
      exit( -1 );
    }

    kc[comtype]= getXMLcmntArrayDouble( infile,
                                        "vegECD",
                                        "kc",
                                        comtype );

    ki[comtype]= getXMLcmntArrayDouble( infile,
                                        "vegECD",
                                        "ki",
                                        comtype );

    gva[comtype]= getXMLcmntArrayDouble( infile,
                                         "vegECD",
                                         "gva",
                                         comtype );

    tmin[comtype]= getXMLcmntArrayDouble( infile,
                                          "vegECD",
                                          "tmin",
                                          comtype );

    toptmin[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "toptmin",
                                             comtype );

    toptmax[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "toptmax",
                                             comtype );

    tmax[comtype]= getXMLcmntArrayDouble( infile,
                                          "vegECD",
                                          "tmax",
                                          comtype );

//    raq10a0[comtype]= getXMLcmntArrayDouble( infile,
//                                             "vegECD",
//                                             "raq10a0",
//                                             comtype );

//    raq10a1[comtype]= getXMLcmntArrayDouble( infile,
//                                             "vegECD",
//                                             "raq10a1",
//                                             comtype );

//    raq10a2[comtype]= getXMLcmntArrayDouble( infile,
//                                             "vegECD",
//                                             "raq10a2",
//                                             comtype );

//    raq10a3[comtype]= getXMLcmntArrayDouble( infile,
//                                             "vegECD",
//                                             "raq10a3",
//                                             comtype );

    tref[comtype]= getXMLcmntArrayDouble( infile,
                                          "vegECD",
                                          "tref",
                                          comtype );

    qref[comtype]= getXMLcmntArrayDouble( infile,
                                          "vegECD",
                                          "qref",
                                          comtype );

    alpha[comtype]= getXMLcmntArrayDouble( infile,
                                           "vegECD",
                                           "alpha",
                                           comtype );

    beta[comtype]= getXMLcmntArrayDouble( infile,
                                          "vegECD",
                                          "beta",
                                          comtype );

    gamma[comtype]= getXMLcmntArrayDouble( infile,
                                           "vegECD",
                                           "gamma",
                                           comtype );


    kvnh4[comtype]= getXMLcmntArrayDouble( infile,
                                           "vegECD",
                                           "kvnh4",
                                           comtype );

    kvno3[comtype]= getXMLcmntArrayDouble( infile,
                                           "vegECD",
                                           "kvno3",
                                           comtype );

    labncon[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "labncon",
                                             comtype );

    leafmxc[comtype]= getXMLcmntArrayDouble( infile,
                                             "vegECD",
                                             "leafmxc",
                                             comtype );

    kleafc[comtype]= getXMLcmntArrayDouble( infile,
                                            "vegECD",
                                            "kleafc",
                                            comtype );

    sla[comtype]= getXMLcmntArrayDouble( infile,
                                         "vegECD",
                                         "sla",
                                         comtype );

    cov[comtype]= getXMLcmntArrayDouble( infile,
                                         "vegECD",
                                         "cov",
                                         comtype );

    fpcmax[comtype]= getXMLcmntArrayDouble( infile,
                                            "vegECD",
                                            "fpcmax",
                                            comtype );

    proptrans[comtype]= getXMLcmntArrayDouble( infile,
                                                "vegECD",
                                                "proptrans",
                                                comtype );
    ifwoody[comtype]= getXMLcmntArrayInt( infile,
                                                "vegECD",
                                                "ifwoody",
                                                comtype );
    rsratio[comtype]= getXMLcmntArrayDouble( infile,
                                                "vegECD",
                                                "rsratio",
                                                comtype );
    lfratiomax[comtype]= getXMLcmntArrayDouble( infile,
                                                "vegECD",
                                                "lfratiomax",
                                                comtype );
    woodfall[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "woodfall",
                                              comtype );
    standdeadpar[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "standdeadpar",
                                              comtype );
    mrcoef[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "mrcoef",
                                              comtype );
    mrcoef_leaf[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "mrcoef_leaf",
                                              comtype );
    mrcoef_stem[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "mrcoef_stem",
                                              comtype );
    mrcoef_croot[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "mrcoef_croot",
                                              comtype );
    mrcoef_froot[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "mrcoef_froot",
                                              comtype );

   agecoef[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "agecoef",
                                              comtype );
   mortcoefa[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "mortcoefa",
                                              comtype );
   mortcoefb[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "mortcoefb",
                                              comtype );
   q10_cgs[comtype]= getXMLcmntArrayDouble( infile,
                                              "vegECD",
                                              "q10_cgs",
                                              comtype );
    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in vegECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::getleafecd( ofstream& rflog1 )
{

  string ecd;

#ifdef PMODE
  *fgo >> ecd;

#else
  cout << "Enter name of the file with leaf parameter";
  cout << " values (.ECD):" << endl;

  cin >> ecd;
#endif

  rflog1 << "Enter name of the file with leaf parameter";
  rflog1 << " values (.ECD):";
  rflog1 << ecd << endl << endl;

  getleafecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::getleafecd( const string& ecd )
{
  int comtype;
  int dcmnt;
  ifstream infile;


  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for leaf ECD input" << endl;
    
    exit( -1 );
  }

  getXMLrootNode( infile, "leafECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "leafECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than " << (MAXCMNT-1);
      cerr << " in leafECD" << endl;
      
      exit( -1 );
    }

    minleaf[comtype]= getXMLcmntArrayDouble( infile,
                                             "leafECD",
                                             "minleaf",
                                             comtype );

    aleaf[comtype]= getXMLcmntArrayDouble( infile,
                                           "leafECD",
                                           "aleaf",
                                           comtype );

    bleaf[comtype]= getXMLcmntArrayDouble( infile,
                                           "leafECD",
                                           "bleaf",
                                           comtype );

    cleaf[comtype]= getXMLcmntArrayDouble( infile,
                                           "leafECD",
                                           "cleaf",
                                           comtype );

    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in leafECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::gppxclm( const int& pdcmnt,
                        const double& co2,
                        const double& par,
                        const double& temp,
                        const double& gv,
                        const double& leaf,
                        const double& foliage,
                        const double& thawpercent )
{

  double gpp;


/* *************************************************************
   gpp:    gpp as influenced by carbon dioxide (co2), moisture
	   (gv), phenology (leaf), photosynthetically active
	   radiation (par), air temperature (temp) and freeze-
	   thaw index (thawpercent)

************************************************************* */

  gpp  = co2 * gv;
  gpp *= cmax * foliage / (kc[pdcmnt] + gpp);
  gpp *= leaf * par / (ki[pdcmnt] + par);
  gpp *= temp;
  //gpp *= thawpercent;
  //closed by cgs2014 to prevent model fluctuation due to unequilibrium soil temperature
  //cout <<" co2: " <<co2<<" gv: " <<gv <<" thawpercent: "<<thawpercent<<" temp: "<<temp<<endl;
  return gpp;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::gppxio3( const double& fozone,
                        const double& eetpet )
{
  double findozone;

  if( ZERO == fozone )
  {
    findozone = 1.0;
  }
  else
  {
    findozone = (1.0-1.0/fozone) * eetpet + (1.0/fozone);
  }

//  findozone = ((1.0/fozone)-1) * exp(-5.0*fprevozone) + 1.0;

  return findozone;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::gppxo3( const int& pdcmnt,
                       const double& gpp,
                       const double& d40,
                       const double& eetpet )
{
  double fozone;

  fozone = 1.0 - (o3para[pdcmnt] * pow( 10.0, -6 )
           * (o3parb[pdcmnt] + o3parc[pdcmnt] * gpp) * d40);

  fozone = 1 - eetpet + (eetpet*fozone);

  fozone = fozone + fprevozone - 1.0;

  // Keep fozone between 0.0 and 1.0
  
  if( fozone <= ZERO )
  {
    fozone = ZERO;
  }
  
  if( fozone >= 1.0 )
  { 
    fozone = 1.000000;
  }


  return fozone;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::leafinit( ofstream& rflog1 )
{

  int lfswtch;

#ifdef PMODE
  *fgo >> leafyrs;
  *fgo >> lfswtch;

#else
  cout << "Enter number of years for model run:  " << endl;

  cin >> leafyrs;

  cout << "Do you have a file containing the phenology parameters?:";
  cout << endl;
  cout << "Enter 0 for no:" << endl;
  cout << "Enter 1 for yes:" << endl;

  cin >> lfswtch;
#endif

  rflog1 << "Enter number of years for model run:  ";
  rflog1 << leafyrs << endl;
  rflog1 << "Do you have a file containing the phenology parameters?:";
  rflog1 << endl;
  rflog1 << "Enter 0 for no:" << endl;
  rflog1 << "Enter 1 for yes:" << endl;
  rflog1 << lfswtch << endl;

  if( 0 == lfswtch )
  {
#ifdef PMODE
    *fgo >> aleaf[0];
    *fgo >> bleaf[0];
    *fgo >> cleaf[0];

#else
    cout << "Enter regression coefficient for relative EET";
    cout << " (i.e. 'a'):  ";

    cin >> aleaf[0];

    cout << "Enter regression coefficient for LAI(t-1)";
    cout << " (i.e. 'b'):";

    cin >> bleaf[0];

    cout << "Enter regression intercept (i.e. 'c'):  ";

    cin >> cleaf[0];
#endif
    rflog1 << "Enter regression coefficient for relative EET";
    rflog1 << " (i.e. 'a'): ";
    rflog1 <<  aleaf[0] << endl;
    rflog1 << "Enter regression coefficient for LAI(t-1)";
    rflog1 << " (i.e. 'b'):  ";
    rflog1 << bleaf[0] << endl;
    rflog1 << "Enter regression intercept (i.e. 'c'):  ";
    rflog1 << cleaf[0] << endl;

    minleaf[0] = 2.0;
    
    while( minleaf[0] >= 1.0 )
    {

#ifdef PMODE
      *fgo >> minleaf[0];

#else
      cout << "Enter minimum LAI (must be less than 1.0):";
      cin >> minleaf[0];
#endif
      rflog1 << "Enter minimum LAI (must be less than 1.0): ";
      rflog1 << minleaf[0];
      rflog1 << endl << endl;
    }
  }
  else { getleafecd( rflog1 ); }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::nupnh4xclm( const int& pdcmnt,
                           const double& soilh2o,
                           const double& nh4,
                           const double& respq10,
                           const double& ksoil,
                           const double& foliage,
                           const double& fozone )
{

/* *************************************************************
   nuptake:  uptake of ammonium by plants as influenced by
	     ammonium concentration (nh4), moisture
	     (ksoil), ozone, and air temperature (respq10)
************************************************************* */

  double nuptake;

  nuptake  = (nh4 * ksoil) / soilh2o;
  nuptake *= nupnh4 * foliage / (kvnh4[pdcmnt] + nuptake);
  nuptake *= respq10;
  nuptake *= fozone;  

  return nuptake;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::nupno3xclm( const int& pdcmnt,
                           const double& soilh2o,
                           const double& no3,
                           const double& respq10,
                           const double& ksoil,
                           const double& foliage,
                           const double& fozone )
{

/* *************************************************************
   nuptake:  uptake of nitrate by plants as influenced by
	     soil nitrate concentration (no3), moisture
	     (ksoil), ozone, and air temperature (respq10)
************************************************************* */

  double nuptake;

  nuptake  = (no3 * ksoil) / soilh2o;
  nuptake *= nupno3 * foliage / (kvno3[pdcmnt] + nuptake);
  nuptake *= respq10;
  nuptake *= fozone;

  return nuptake;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::nupxclm( const int& pdcmnt,
                        const double& soilh2o,
                        const double& availn,
                        const double& respq10,
                        const double& ksoil,
                        const double& foliage,
                        const double& fozone )
{

/* *************************************************************
   nuptake:  uptake of nitrogen by plants as influenced by
	     available nitrogen concentration (availn), moisture
	     (ksoil), ozone, and air temperature (respq10)
************************************************************* */

  double nuptake;

  nuptake  = (availn * ksoil) / soilh2o;
  nuptake *= nmax * foliage / (kn1[pdcmnt] + nuptake);
  nuptake *= respq10;
  nuptake *= fozone;

  return nuptake;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void Tveg60::resetEcds( const int& pcmnt,
                        const double& psiplusc )
{
  // Initialize TEM parameters dependent upon a grid cell's
  //   soil texture

  if( psiplusc <= cmaxcut[pcmnt] )
  {
    cmax = (cmax1a[pcmnt] * psiplusc) + cmax1b[pcmnt];
  }
  else
  {
    cmax = (cmax2a[pcmnt] * psiplusc) + cmax2b[pcmnt];
  }

  if( cmax < ZERO ) { cmax = ZERO; }
  
  if( psiplusc <= nupnh4cut[pcmnt] )
  {
    nupnh4 = (nupnh41a[pcmnt] * psiplusc) + nupnh41b[pcmnt];
  }
  else
  {
    nupnh4 = (nupnh42a[pcmnt] * psiplusc) + nupnh42b[pcmnt];
  }

  if( psiplusc <= nupno3cut[pcmnt] )
  {
    nupno3 = (nupno31a[pcmnt] * psiplusc) + nupno31b[pcmnt];
  }
  else
  {
    nupno3 = (nupno32a[pcmnt] * psiplusc) + nupno32b[pcmnt];
  }
  //if (pcmnt ==7) cout <<"pcmnt: "<<pcmnt<<" nupno3: "<<nupno3<<" nupnh4: "<<nupnh4<<" psiplusc: "<<psiplusc<<" nupno31a: "<<nupno31a[pcmnt]<<" nupno32a: "<<nupno32a[pcmnt]<<" nupno31b: "<<nupno31b[pcmnt]<<" nupno32b: "<<nupno32b[pcmnt]<<endl;
  if( nupnh4 < ZERO ) { nupnh4 = ZERO; }
  if( nupno3 < ZERO ) { nupno3 = ZERO; }
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void Tveg60::resetMonthlyFluxes( void )
{
  // Reset monthly fluxes to zero
  
  ingpp = ZERO;
  gpp = ZERO;
  fozone = ZERO;
  findozone = ZERO;

  innpp = ZERO;
  npp = ZERO;
  
  gpr = ZERO;
  rm = ZERO;
  rg = ZERO;
  abvgrndResp = ZERO;
  rootResp = ZERO;
  npp_leaf = ZERO;
  npp_stem = ZERO;
  npp_croot = ZERO;
  npp_froot = ZERO;
  
  deadwoodltc=ZERO;
  deadwoodltn=ZERO;

//  voc.isoprene = ZERO;
//  voc.monoterpene = ZERO;
//  voc.otherReactive = ZERO;
//  voc.other = ZERO;
//  voc.total = ZERO;
  
  ltrfal.carbon = ZERO;

  nfix = ZERO;

  inuptake.total = ZERO;
  inuptake.nh4 = ZERO;
  inuptake.no3 = ZERO;

  nuptake.total = ZERO;
  nuptake.nh4 = ZERO;
  nuptake.no3 = ZERO;

  suptake = ZERO;
  luptake = ZERO;
  nmobil = ZERO;
  nresorb = ZERO;

  ltrfal.nitrogen = ZERO;
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::resetNEWTOPT( const int& pcmnt, 
                           const double& tair, 
                           const double& unnrmleaf )
{
  if( ZERO == aleaf[pcmnt] 
      && ZERO == bleaf[pcmnt]
      && 1.0 == cleaf[pcmnt] )
  {
    if( tair > newtopt ) { newtopt = tair; }
  }
  else
  {
    if( unnrmleaf > newleafmx )
    {
      newleafmx = unnrmleaf;
      newtopt = tair;
    }
  }
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::resetYrFluxes( void )
{
  // Reset annual fluxes and summary variables to zero

  yrcarbon = ZERO;

  yrstructn = ZERO;
  yrstoren = ZERO;

  yrnitrogen = ZERO;
  yrc2n = ZERO;

  yrthawpct = ZERO;

  // Phenology

  yrunleaf = ZERO;
  yrleaf = ZERO;
  yrlai = ZERO;
  yrfpc = ZERO;

  // Carbon fluxes
  
  yringpp = ZERO;
  yrgpp = ZERO;
  yrinnpp = ZERO;
  yrnpp = ZERO;
  yrgpr = ZERO;
  yrrmaint = ZERO;
  yrrgrowth = ZERO;
  yrabvgrndResp = ZERO;
  yrrootResp = ZERO;

  yrvoc.isoprene = ZERO;
  yrvoc.monoterpene = ZERO;
  yrvoc.otherReactive = ZERO;
  yrvoc.other = ZERO;
  yrvoc.total = ZERO;

  yrltrc = ZERO;

  // Nitrogen fluxes
  
  yrnfix = ZERO;

  yrinnup.total = ZERO;
  yrinnup.nh4 = ZERO;
  yrinnup.no3 = ZERO;

  yrnup.total = ZERO;
  yrnup.nh4 = ZERO;
  yrnup.no3 = ZERO;

  yrsup = ZERO;
  yrlup = ZERO;
  yrnmobil = ZERO;
  yrnrsorb = ZERO;

  yrltrn = ZERO;
  yrdeadwoodc= ZERO;
  yrdeadwoodn= ZERO;

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::rmxclm( const int& pdcmnt,
                       const double& tcarbon,
                       const double& tmrcoef,
                       const double& respq10 )
{

/* *************************************************************
   rm: plant maintenance respiration as influenced by plant
       biomass (vegc) and air temperature (respq10)
************************************************************* */

  const double kr25 = 9.648;
  
  double krm;

  double rmaint;

//  kr = exp( (kra[pdcmnt]*vegc) + krb[pdcmnt] );
//  rmaint  = kr * vegc;
//  rmaint *= respq10;

  //added by cgs to calculate RM

  //double gtempair = powf(respq10, (temptair - 25)/10.0);
  rmaint= tmrcoef * respq10 * tcarbon;

  /*
  krm = rmmax[pdcmnt] * initcneven[pdcmnt] / kr25;
  rmaint = vegc * rmmax[pdcmnt] / (krm + vegc);
  rmaint *= respq10;
  */

  if (rmaint <= 0.00001 && rmaint >-0.00001) rmaint = 0.0;
  //cout <<"pdcmnt: " <<pdcmnt<<" vegc: " <<vegc<<" krm: " <<krm <<"rmmax: "<<rmmax[pdcmnt]<<" initcneven[pdcmnt]: " <<initcneven[pdcmnt] << "rmaint: "<<rmaint<<" respq10: " <<respq10<<endl;
  return rmaint;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

//double Tveg60::rq10( const int& pdcmnt, const double& tair )
//{

//  double raq10;

/* *************************************************************
 rq10: effect of temperature on plant respiration
************************************************************* */

//  raq10 = raq10a0[pdcmnt] + (raq10a1[pdcmnt]*tair)
//          + (raq10a2[pdcmnt] * pow( tair,2.0 ))
//          + (raq10a3[pdcmnt] * pow( tair,3.0 ));

//  return pow( raq10,tair/10.0 );

//};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Tveg60::setGV( const double& eet,
                      const double& pet,
                      const int& moistlim )
{
  double tstgv;

/* tstgv: effect of moisture on primary productivity */

  if( 0 == moistlim )
  {
    tstgv = 1.0;
  }
  else
  {
    if( eet / pet <= 0.1 )
    {
      tstgv = (-10.0 * pow( (eet / pet), 2.0 ))
              + (2.9 * (eet / pet));
      
      if( tstgv < ZERO ) { tstgv = ZERO; }
    }
    else
    {
      tstgv = 0.1 + (0.9 * eet / pet);
    }
  }
  if (tstgv >1.0) tstgv =1.0;
  return tstgv;
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::setNewRESPQ10( const int& pdcmnt, const double& tair )
{
/* *************************************************************
 respq10: effect of temperature on plant respiration

 Amthor, J. (in review)
************************************************************* */

  double raq10;

  raq10 = qref[pdcmnt]
          * exp( -1.0 * alpha[pdcmnt] * (tair - tref[pdcmnt]) );

  respq10 = pow( raq10, (tair - tref[pdcmnt])/10.0 )
            / ( 1.0 + exp( beta[pdcmnt] - tair )
                + exp( tair - gamma[pdcmnt] ));
  //if (respq10 <= 0.00001 && respq10> -0.00001) respq10 = 0.0; //added by cgs
  //cout <<"pdcmnt: " <<pdcmnt<<" raq10: "<<raq10<<" qref: "<<qref[pdcmnt]<<" alpha: " <<alpha[pdcmnt]<< "tref: "<<tref[pdcmnt]<<endl;
  //cout <<"pdcmnt: " <<pdcmnt<<" raq10: "<< raq10<< " tair: " <<tair<<" beta: " <<beta[pdcmnt]<<" respq: " <<respq10<<endl;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

//void Tveg60::setRESPQ10( const int& pdcmnt, const double& tair )
//{
//  double raq10;

//  raq10 = raq10a0[pdcmnt] + (raq10a1[pdcmnt] * tair)
//          + (raq10a2[pdcmnt]*pow( tair, 2.0 ))
//          + (raq10a3[pdcmnt]*pow( tair, 3.0 ));

//  respq10 = pow( raq10, tair/10.0 );

//};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::setTEMP( const int& pdcmnt, const double& tair )
{
  /* temp: effect of temperature on primary productivity */
  //below code is slightly modified by cgs2014 to avoid fluctuation of topt

	double topt_mean = (toptmax[pdcmnt] + toptmin[pdcmnt])/2.0;
  //double topt_mean = topt;

  if( tair <= tmin[pdcmnt] || tair >= tmax[pdcmnt] )
  {
    temp = 0.0;
  }
  else
  {
    if( tair >= topt_mean && tair <= toptmax[pdcmnt] )
    {
      temp = 1.0;
    }
    else
    {
      if( tair > tmin[pdcmnt] && tair < topt_mean )
      {
	temp = (tair - tmin[pdcmnt])
               * (tair - tmax[pdcmnt])
               / ((tair - tmin[pdcmnt])
               * (tair - tmax[pdcmnt])
               - pow( (tair - topt_mean), 2.0 ));
      }
      else
      {
	temp = (tair - tmin[pdcmnt])
               * (tair - tmax[pdcmnt])
               / ((tair - tmin[pdcmnt])
               * (tair - tmax[pdcmnt])
               - pow( (tair - toptmax[pdcmnt]), 2.0 ));
      }
    }
  }

  /*
  //cout <<" pdcmnt: "<<pdcmnt<<" tair: "<<tair<<" temp: "<<temp<<" tmin: "<<tmin[pdcmnt]<<endl;
  //below added by cgs2014 to replace the calculation for temperature effects.
  double topt_mean = (toptmax[pdcmnt] + toptmin[pdcmnt])/2.0;
  double f1_stomata = pow(q10_cgs[pdcmnt], ((tair-topt_mean)/10.));
  double tempvar=1. + exp((-2.2e05+710.*(tair+273.16))/(8.314*(tair+273.16)));
  double f2_stomata=0.0;
  if (tempvar ==0.0) f2_stomata = 0.0;
  else f2_stomata = 1/tempvar;
  temp = f1_stomata * f2_stomata;
  if( tair <= tmin[pdcmnt] || tair >= tmax[pdcmnt])
  {
    temp = 0.0;
  }
  */

	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::setThawPercent( const double& pprevdst10,
                             const double& pdst10,
                             const double& pnextdst10 )
{

  double dst10;
  double prevdst10;
  double nextdst10;
  
  // Calculate the thawing-proportion in early month and later fall

  if( ZERO == pdst10 ) { dst10 = 0.001; }
  else { dst10 = pdst10; }
  
  if( ZERO == pprevdst10 ) { prevdst10 = 0.002; }
  else { prevdst10 = pprevdst10; }
  
  if( ZERO == pnextdst10 ) { nextdst10 = 0.003; }
  else { nextdst10 = pnextdst10; }
  
  
  // Conditions of - - -
  if( (prevdst10 < ZERO)
       && (dst10 < ZERO)
       && (nextdst10 < ZERO) )
  {
    thawpercent = ZERO;
  }

  // Conditions of + + +
  if( (prevdst10 > ZERO)
       && (dst10 > ZERO)
       && (nextdst10 > ZERO) )
  {
    thawpercent = 1.0;
  }

  // Conditions of + - +
  if( (prevdst10 > ZERO)
       && (dst10 < ZERO)
       && (nextdst10 > ZERO) )
  {
     thawpercent = ((prevdst10 + dst10 + nextdst10) / 3.0 ) / 3.0;
  }

  // Conditions of + + -
  if( (prevdst10 > ZERO)
       && (dst10 > ZERO)
       && (nextdst10 < ZERO) )
  {
     thawpercent = ((prevdst10 + dst10 + nextdst10) / 3.0 ) / 4.0;    
  }

  // Conditions of - + +
  if( (prevdst10 < ZERO)
       && (dst10 > ZERO)
       && (nextdst10 > ZERO) )
  {
     thawpercent = ((prevdst10 + dst10 + nextdst10) / 3.0 ) / 6.0;    
  }

  // Conditions of - - +
  if( (prevdst10 < ZERO)
       && (dst10 < ZERO)
       && (nextdst10 > ZERO) )
  {
     thawpercent = ((prevdst10 + dst10 + nextdst10) / 3.0 ) / 5.0;    
  }

  // Conditions of + - -
  if( (prevdst10 > ZERO)
       && (dst10 < ZERO)
       && (nextdst10 < ZERO))
  {
     thawpercent = ((prevdst10 + dst10 + nextdst10) / 3.0 ) / 4.0;    
  }

  // Conditions of - + -
  if( (prevdst10 < ZERO)
       && (dst10 > ZERO)
       && (nextdst10 < ZERO) )
  {
     thawpercent = ((prevdst10 + dst10 + nextdst10) / 3.0 ) / 3.0;
  }


  if( thawpercent > 1.0 ) { thawpercent = 1.000000; }
  if( thawpercent < ZERO ) { thawpercent = ZERO; }
  //cout <<" dst10: " <<dst10<< " prevdst10: "<<prevdst10<<" nextdst10: " <<nextdst10<< " thawpercent: "<<thawpercent<<endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::showecd( const int& pdcmnt )
{

  cout << endl << "             VEGETATION PARAMETERS INFLUENCED BY CLIMATE";
  cout << endl << endl;
  printf( "  KI = %6.2lf  KC = %6.2lf  KVNH4 = %6.4lf  KVNO3 = %6.4lf  GVA = %8.4lf\n\n",
          ki[pdcmnt],
          kc[pdcmnt],
          kvnh4[pdcmnt],
          kvno3[pdcmnt],
          gva[pdcmnt] );

  printf( "   TMIN = %5.1lf    TOPTMIN = %5.1lf   TOPTMAX = %5.1lf   TMAX = %5.1lf\n\n",
          tmin[pdcmnt],
          toptmin[pdcmnt],
          toptmax[pdcmnt],
          tmax[pdcmnt] );

//  printf( " RAQ10A0 = %7.5lf  RAQ10A1 = %7.5lf  RAQ10A2 = %7.5lf  RAQ10A3 = %7.5lf\n",
//          raq10a0[pdcmnt],
//          raq10a1[pdcmnt],
//          raq10a2[pdcmnt],
//          raq10a3[pdcmnt] );

  printf( " TREF = %7.5lf  QREF = %7.5lf  ALPHA = %7.5lf  BETA = %7.5lf  GAMMA = %7.5lf\n",
		  tref[pdcmnt],
		  qref[pdcmnt],
		  alpha[pdcmnt],
		  beta[pdcmnt],
		  gamma[pdcmnt] );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::showleaf( const int& pdcmnt )
{

  cout << endl << "         PARAMETERS FOR THE LEAF PHENOLOGY MODEL";
  cout << endl << endl;
  printf( "     ALEAF = %7.5lf     BLEAFC = %8.5lf        CLEAF = %8.5lf\n",
          aleaf[pdcmnt],
          bleaf[pdcmnt],
          cleaf[pdcmnt] );

  printf( "   MINLEAF = %4.2lf       MAXLEAF = %7.4lf      UNLEAF12 = %7.4lf\n",
          minleaf[pdcmnt],
          initleafmx[pdcmnt],
          unleaf12[pdcmnt] );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::updateC2N( const int& pdcmnt,
                        const double& yreet,
                        const double& yrpet,
                        const double& currentco2,
                        const double& initco2 )
{
  if( yrpet != ZERO )
  {
    c2n = c2nb[pdcmnt] + c2na[pdcmnt]*(yreet/yrpet);
  }
  else { c2n = c2nb[pdcmnt]; }

  if( c2n < c2nmin[pdcmnt] ) { c2n = c2nmin[pdcmnt]; }
  
  adjc2n = 1.0 + (dc2n * (currentco2 - initco2));
  c2n *= adjc2n;
  cneven = initcneven[pdcmnt] * adjc2n;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::updateDynamics( const double& lat,
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
                             double& forestage)
{	
  double eetpet;

  // Determine stage of seasonal phenology

  if( 0 == moistlim )
  {
    unnormleaf = deltaleaf( pdcmnt,
                            pet,
                            prevmaxpet,
                            prevunrmleaf );
  }
  else
  {
    unnormleaf = deltaleaf( pdcmnt,
                            eet,
                            prevmaxeet,
                            prevunrmleaf );                           
  }


  if( prvleafmx <= ZERO ) { leaf = ZERO; }
  else { leaf = unnormleaf / prvleafmx; }

  if( leaf < minleaf[pdcmnt] )
  {
    leaf = minleaf[pdcmnt];
  }

  if( leaf > 1.0 ) { leaf = 1.0; }
  //if (pdcmnt == 7) cout <<"comtype: "<<pdcmnt<<" minleaf: "<< minleaf[pdcmnt] <<" unnomleaf :"<<unnormleaf<<" leaf: "<<leaf<<endl;

  // gv: effect of moisture on primary productivity

  gv = setGV( eet, pet, moistlim );

  //if (pdcmnt ==4) cout <<" cmnt: "<<pdcmnt<<" gv: " <<gv<<" eet: "<<eet<<" pet: "<<pet<<" moistlim: "<<moistlim<<endl;
  // Determine Gross Primary Production if nitrogen is not
  //   limiting (veg.ingpp)

  ingpp = gppxclm( pdcmnt,
                   co2,
                   par,
                   temp,
                   gv,
                   leaf,
                   foliage,
                   thawpercent );

  //added by cgs to include forestage effects on carbon recovery after disturbance
  double ageeff=1.0;

  //int nonforest[6] ={1,2,3,7,8,12};
  int nonforest[4] ={1,2,3,7}; //herbaceous
  int index=0;
  for (int i=0;i<4;i++){
  		if (pdcmnt ==nonforest[i]){
  		index=1;
  	    }
  }
  /*
  if (index ==1) ageeff=1.0;
  else ageeff = 1-exp(agecoef[pdcmnt] *forestage);
  if (ageeff<0.5) ageeff=0.5;
  ingpp*=ageeff;
  //if (index==1) cout <<"index: "<<index<<" ingpp: "<<ingpp<<" ageeff: "<<ageeff<<" standage: "<<forestage<<endl;
 */

  if( ingpp < ZERO ) { ingpp = ZERO; }


  /* ozone: effect of surface ozone on gpp */

  if( 0 == o3flag )
  {
    fozone = 1.000000;
    fprevozone = 1.000000;
    findozone = 1.000000;
  }
  else
  {
    eetpet = eet / pet;

    fozone = gppxo3( pdcmnt,
                     ingpp,
                     aot40,
                     eetpet );

    ingpp *= fozone;

//    findozone = gppxio3( fozone, eetpet );
//    ingpp *= findozone;
  }

  if( ingpp < ZERO ) { ingpp = ZERO; }

  //if (cmnt ==5 && lat==29.5 && lon==-98.5) cout <<" month: "<<month<<" ingpp: " <<ingpp<<" co2: " <<co2<<" par: " <<par<<" temptair: " <<temptair<<" temp: " <<temp<<" gv: " <<gv << " leaf: "<<leaf<<" foliage: "<<foliage<<" thawpercent: " <<thawpercent<<" topt: "<<topt<<" q10: "<<respq10<<endl;


  // Determine Maintenance Respiration of vegetation (veg.rm)

  rm = rmxclm( pdcmnt, vegc, mrcoef[pdcmnt], respq10);
  //cout <<"rm: "<<rm<<endl;

  //modified by cgs to avoid
  if( rm <= 0.00001 && rm >= -0.00001) { rm = 0.0; }


  // Determine Net Primary Production if nitrogen is not
  //   limiting (veg.innpp) and Growth Respiration (veg.rg)

  innpp = ingpp - rm;

  rg = ZERO;
  if( innpp > 0.00001 )
  {
    rg  = 0.2 * innpp;
    innpp *= 0.8;
  }

  //if (pdcmnt ==5) cout <<"cmnt: "<<pdcmnt<<" ingpp: " <<ingpp<<" rm: "<<rm<<" innpp: "<<innpp<<endl;


  /*
  ltrfal.carbon = cfall[pdcmnt] * vegc;
  if( ltrfal.carbon < 0.00001 &&  ltrfal.carbon > -0.00001)
  {
    ltrfal.carbon = ZERO;
  }

  ltrfal.nitrogen = nfall[pdcmnt] * structn;
  
  if( ltrfal.nitrogen < 0.00001 &&  ltrfal.nitrogen > -0.00001 )
  {
    ltrfal.nitrogen = ZERO;
  }
 */

  //added by cgs2014 to simulate age effects on tree mortality
  //default mortcoefa=0.5, mortcoefb=0.08
  double agemortc, agemortn;

  if (cfall[pdcmnt]<=0.0) {agemortc=0.0;agemortn=0.0;}
  else
  {
	  agemortc = (1+mortcoefa[pdcmnt] * exp(mortcoefb[pdcmnt]*forestage))* cfall[pdcmnt];
	  agemortn = (1+mortcoefa[pdcmnt] * exp(mortcoefb[pdcmnt]*forestage))* nfall[pdcmnt];
  }

  //agemortc=cfall[pdcmnt];
  //agemortn=nfall[pdcmnt];
  if (ifwoody[pdcmnt] !=1) agemortc=cfall[pdcmnt]; //herbaceous plants mortality is not related to age
  if (ifwoody[pdcmnt] !=1) agemortn=nfall[pdcmnt];
  ltrfal.carbon = agemortc * vegc;
  ltrfal.nitrogen = agemortn * structn;

  if( ltrfal.carbon < 0.00001 &&  ltrfal.carbon > -0.00001)
  {
    ltrfal.carbon = ZERO;
  }
  if( ltrfal.nitrogen < 0.00001 &&  ltrfal.nitrogen > -0.00001 )
  {
    ltrfal.nitrogen = ZERO;
  }


  //standing deadwood decay. moved to ttem604.cpp
  //deadwoodltc = woodfall[pdcmnt] * deadwoodc;
  //deadwoodltn = woodfall[pdcmnt] * deadwoodn;

  // Determine nitrogen uptake by vegetation if carbon is not
  //   limiting (veg.inuptake)

  inuptake.nh4 = nupnh4xclm( pdcmnt,
                             soilh2o,
                             soilnh4,
                             respq10,
                             ksoil,
                             foliage,
                             fozone );

  inuptake.no3 = nupno3xclm( pdcmnt,
                             soilh2o,
                             soilno3,
                             respq10,
                             ksoil,
                             foliage,
                             fozone );

  //if (pdcmnt==7 && currentveg == 50) cout <<"pdcmnt0: "<<pdcmnt<<" soilno3: "<<soilno3<<" soilnh4: "<<soilnh4<<" inup.no3: " <<inuptake.no3<<" inup.nh4: "<<inuptake.nh4<<endl;

  if( inuptake.nh4 > soilnh4 + ndep.nh4 + extraNH4 + netnmin 
                     - ammnvol - nitrif )
  {
    inuptake.nh4 = soilnh4 
                   + ndep.nh4
                   + extraNH4 
                   + netnmin 
                   - ammnvol 
                   - nitrif;
  }

  //if (cmnt ==7 && currentveg == 50) cout <<" cmnt0.5: "<<cmnt<<" inup_nh4: "<<inuptake.nh4<<" inup_no3: "<<inuptake.no3<<" soilh2o: "<<soilh2o<<" soilno3: "<<soilno3<<" soilnh4: "<<soilnh4<<" nupnh4: "<<nupnh4<<" nupno3: "<<nupno3<<" ammnvol: "<<ammnvol<<endl;
  if( inuptake.nh4 < ZERO )
  {
    inuptake.nh4 = ZERO;
  }


  if( inuptake.no3 > soilno3 + ndep.no3 + no3prod )
  {
    inuptake.no3 = soilno3 + ndep.no3 + no3prod;
  }

  if( inuptake.no3 < ZERO ) { inuptake.no3 = ZERO; }


  inuptake.total = inuptake.nh4 + inuptake.no3;


  // Determine how interactions between carbon and nitrogen
  //   availability influence primary production, litterfall
  //   and nitrogen uptake

  gpp = ingpp;
  npp = innpp;

  // Assume CROP NPP is always positive or zero
/*//closed by cgs2014
  if( 1 == agstate && 0 == agperennial )
  {
     if( npp < ZERO ) { npp = ZERO; }
  }
*/
 
  nuptake.total = inuptake.total;

//  if ( nuptake.total < ZERO )
//  {
//    nuptake.total = ZERO;
//  }

  nuptake.no3 = inuptake.no3;

  nuptake.nh4 = nuptake.total - nuptake.no3;

  if( nuptake.nh4 < ZERO )
  {
    nuptake.nh4 = ZERO;
    nuptake.no3 = nuptake.total;
  }

  suptake = nuptake.total;
  luptake = ZERO;

  nmobil = ZERO;
  nresorb = ZERO;
  //if (cmnt ==7 && currentveg == 50) cout <<" strn1: "<<structn<<" suptake: "<<suptake<<" inuptake: "<<inuptake.total<<endl;
  //if (structn != structn) exit (-1);

  // Nitrogen feedback of GPP ( 1 == nfeed)

  //if( 1 == nfeed ) //closed by cgs2014 to change the nitrogen uptake under nfeed==0
  {
    // Determine nitrogen resorption (veg.nresorb)
    
    if( ltrfal.nitrogen
         <= ltrfal.carbon / cneven )
    {
      nresorb = ltrfal.carbon / cneven - ltrfal.nitrogen;
    }
    else
    {
      nresorb = ZERO;
    }

    if( vegc > ZERO )
    {
      nresorb *= (structn / vegc) * c2n;
    }


    // Determine if monthly primary production is carbon-limited  
    //   or nitrogen-limited
    
    if( (nuptake.total + labilen + nfix) < 0.000001 )
    {
      inprodcn = innpp / 0.000001;
    }
    else 
    {
      inprodcn = innpp / (nuptake.total + labilen + nfix);
    }
    //if (cmnt==7) cout <<"nfix: "<<nfix<<endl;

    
    // If primary production is nitrogen-limited, 
    //   (i.e., veg.inprodcn > veg.cneven) revaluate NPP, RG and 
    //   GPP based on nitrogen availability
     
    if( inprodcn > cneven)
    {
      if( 1 == fertflag || nfeed==0) //modified by cgs2014 to add nitrogen when nfeed==0
      //if( 1 == fertflag)
      {
        // Assume nitrogen is never limiting if fertilized
        // Also assume that fertilized crops are based solely 
        //   on a nitrate economy
          
        nuptake.no3 = (npp / cneven) - labilen - nfix;

        if( nuptake.no3 > inuptake.no3 )
        {
          agfertn = nuptake.no3 - inuptake.no3;
        }

        nuptake.total = nuptake.no3;
        nuptake.nh4 = ZERO;
      }
      else
      {
        npp = cneven * (nuptake.total + labilen + nfix);

        //if( npp < 0.00001 ) { npp = ZERO; }
        double tempnpp = npp;
        if (tempnpp<0.0001) tempnpp=0.0;
        rg = 0.25 * tempnpp;

        gpp = npp + rg + rm;

        if( gpp < ZERO ) { gpp = ZERO; }

        nmobil = labilen;
      }
    }

    // If primary production is carbon-limited, 
    //   (i.e., veg.inprodcn < veg.cneven) revaluate nitrogen
    //   uptake by vegetation based on carbon availability

    if( inprodcn <= cneven )
    {
      nuptake.nh4 = nuptake.nh4
                    * (inprodcn - cnmin[pdcmnt])
                    * (inprodcn - 2 * cneven + cnmin[pdcmnt]);

      nuptake.nh4 /= ((inprodcn - cnmin[pdcmnt])
                     * (inprodcn - 2 * cneven + cnmin[pdcmnt]))
                     - pow( inprodcn - cneven, 2.0 );


      if( nuptake.nh4 > soilnh4 + ndep.nh4 + extraNH4 + netnmin 
                        - ammnvol - nitrif )
      {
        nuptake.nh4 = soilnh4 
                      + ndep.nh4
                      + extraNH4 
                      + netnmin 
                      - ammnvol 
                      - nitrif;
      }

      if ( nuptake.nh4 < ZERO )
      {
        nuptake.nh4 = ZERO;
      }

      nuptake.no3 = nuptake.no3
                    * (inprodcn - cnmin[pdcmnt])
                    * (inprodcn - 2 * cneven + cnmin[pdcmnt]);

      nuptake.no3 /= ((inprodcn - cnmin[pdcmnt])
                     * (inprodcn - 2 * cneven + cnmin[pdcmnt]))
                     - pow( inprodcn - cneven, 2.0 );


      if( nuptake.no3 > soilno3 + ndep.no3 + no3prod )
      {
        nuptake.no3 = soilno3 + ndep.no3 + no3prod;
      }

      if ( nuptake.no3 < ZERO )
      {
        nuptake.no3 = ZERO;
      }

      nuptake.total = nuptake.nh4 + nuptake.no3;


      // Drawdown N from labile N pool first to satisfy N
      //   requirement before taking up any nitrogen from the 
      //   soil
      
      if ( labilen >= (npp / cneven - nfix) )
      {
        nmobil = npp / cneven - nfix;
	
        if ( nmobil < ZERO && vegc > 0.01 )
        {
        	nmobil *= (structn / vegc) * c2n;
        }
	
        suptake = ZERO;
      }
      else
      {
        nmobil = labilen;
	
        suptake = (npp / cneven) - nmobil - nfix;

        if( suptake < ZERO )
        {
        	suptake = ZERO;
        }
	
        if( suptake > nuptake.total )
        {
          suptake = nuptake.total;
        }
      }

      // Determine vegetation nitrogen uptake for the labile
      //   N pool
      
      if( (labilen + nuptake.total - suptake + nresorb - nmobil) 
           < (labncon[pdcmnt] * (structn + suptake 
           - ltrfal.nitrogen - nresorb + nmobil + nfix)) )
      {
        luptake = nuptake.total - suptake;
      }
      else
      {
        luptake = (labncon[pdcmnt] 
                  * (structn + suptake - ltrfal.nitrogen
                  - nresorb + nmobil + nfix))  
                  - (labilen + nresorb - nmobil);

        if( luptake < ZERO ) { luptake = ZERO; }

		nuptake.total = suptake + luptake;

        nuptake.nh4 = nuptake.total - nuptake.no3;

        if( nuptake.nh4 < ZERO )
        {
          nuptake.nh4 = ZERO;
          nuptake.no3 = nuptake.total;
        }
      }
    }
  } //end if nfeed
  double xx=soilnh4 + ndep.nh4 + extraNH4 + netnmin
          - ammnvol - nitrif + soilno3 + ndep.no3 + no3prod;
  //if (pdcmnt==7 && currentveg == 50) cout <<"pdcmnt2: "<<pdcmnt<<" nupno3: "<<nupno3<<" soilno3: "<<soilno3<<" soilh2o: " <<soilh2o<<" nuptakeno3: "<<inuptake.no3<<" nupnh4: "<<nupnh4<<" soilnh4: "<<soilnh4<<" soilh2o: " <<soilh2o<<" nuptakenh4: "<<inuptake.nh4<<endl;
  //if (pdcmnt==7 && currentveg == 50 && inuptake.no3!=inuptake.no3) exit(-1);
  //if (pdcmnt==7 && inprodcn >cneven) cout <<" litterc: "<<ltrfal.carbon<<" littn: "<<ltrfal.nitrogen<<" nresorb: "<<nresorb<<" nitrif: "<<nitrif<<" no3prod: "<<no3prod<<" structn: "<<structn<<" vegc :"<<vegc<<" netnmin: "<<netnmin<<endl;
  //if (cmnt ==7 && inprodcn >cneven) cout <<"cmnt: "<<cmnt<<" ingpp: " <<ingpp<<" gpp: "<<gpp<<" inprodcn: "<<inprodcn<<" cneven: " <<cneven<<" nresorb: " <<nresorb<<" totalavn: " <<xx << " nfix: "<<nfix<<" foliage: "<<foliage<<" nuptake: " <<nuptake.total<<endl;
  //if (cmnt ==7 && inprodcn > 300) exit(-1);


 //added by cgs2014 to address gpp/npp/resp for each components

  double alloc_leaf, alloc_stem, alloc_froot, alloc_croot;// NPP fraction allocated to different organs

  double totGPP = gpp;

  if (totGPP==0.0) {alloc_leaf=0.0; alloc_stem=0.0; alloc_croot=0.0;alloc_froot=0.0;}
  else //NPP>0.0 or NPP<0.0
  {
  alloc_leaf = powf(10., (-0.693+1.011*log10(fabs(totGPP)*2.0/10.0)))/ (fabs(totGPP)*2.0/10.0); //reference: Wolf et al. ecological application. 2011.
  alloc_stem = powf(10., (-0.552+1.097*log10(fabs(totGPP)*2.0/10.0)))/ (fabs(totGPP)*2.0/10.0);
  alloc_froot = powf(10., (-0.608+0.987*log10(fabs(totGPP)*2.0/10.0)))/ (fabs(totGPP)*2.0/10.0);
  alloc_croot = 1.0-alloc_leaf-alloc_stem-alloc_froot;
  }
  if (ifwoody[pdcmnt]==0) //grass
  {
	  alloc_leaf=0.3;
	  alloc_stem=0.2;
	  alloc_croot=0.2;
	  alloc_froot=0.3;
  }

  if (ifwoody[pdcmnt]==2) //shrub
  {
	  alloc_leaf=0.3;
	  alloc_stem=0.2;
	  alloc_croot=0.2;
	  alloc_froot=0.3;
  }

  double rm_leaf,rm_stem,rm_croot, rm_froot;

  rm_leaf = rmxclm( pdcmnt, folc, mrcoef_leaf[pdcmnt], respq10);
  rm_stem = rmxclm( pdcmnt, stemc, mrcoef_stem[pdcmnt], respq10);
  rm_croot = rmxclm( pdcmnt, crootc, mrcoef_croot[pdcmnt], respq10);
  rm_froot = rmxclm( pdcmnt, frootc, mrcoef_froot[pdcmnt], respq10);
  double rg_leaf,rg_stem,rg_croot, rg_froot;
  double temp1,temp2,temp3,temp4;
  temp1 = gpp * alloc_leaf - rm_leaf;
  temp2 = gpp * alloc_stem - rm_stem;
  temp2 = gpp * alloc_croot - rm_croot;
  temp2 = gpp * alloc_froot - rm_froot;
  rg_leaf = temp1 *0.25;
  rg_stem = temp2 *0.25;
  rg_croot = temp3 *0.25;
  rg_froot = temp4 *0.25;
  if (rg_leaf<0.0) rg_leaf =0.0;
  if (rg_stem<0.0) rg_stem =0.0;
  if (rg_croot<0.0) rg_croot =0.0;
  if (rg_froot<0.0) rg_froot =0.0;
  npp_leaf = gpp * alloc_leaf - (rm_leaf + rg_leaf);
  npp_stem = gpp * alloc_stem - (rm_stem + rg_stem);
  npp_croot = gpp * alloc_croot - (rm_croot + rg_croot);
  npp_froot = gpp * alloc_froot - (rm_froot + rg_froot);
  npp =npp_leaf + npp_stem+ npp_croot + npp_froot;

  // Determine Gross Plant Respiration (veg.gpr)

   //gpr = rm + rg;
   gpr =rm_leaf+rg_leaf+rm_stem+rg_stem+rm_croot+rg_croot+rm_froot+rg_froot;
   //cout <<"rm: "<<rm<<" rg: " <<rg<<endl;
   //cout <<"gpp: "<<gpp<<" npp: "<<npp<<" gpr: "<<gpr<<endl;

   // Determine Root Respiration (veg.rootResp)

   //rootResp = gpr * rroot[pdcmnt];
   rootResp = rm_croot+rg_croot+rm_froot+rg_froot;
   //if (forestage ==0 && month==8 && gpp==0.0) cout <<" cmnt: "<<pdcmnt<<" ingpp: " <<ingpp<<" co2: " <<co2<<" par: " <<par<<" temptair: " <<temptair<<" temp: " <<temp<<" gv: " <<gv << " leaf: "<<leaf<<" foliage: "<<foliage<<" thawpercent: " <<thawpercent<<" fozone: "<<fozone<<endl;

 //if (lat==29.5&&lon==-98.5&&pdcmnt==5) cout <<"npp_leaf: "<<npp_leaf<<" rg_leaf: "<<rg_leaf<<" mrcoef_leaf: "<<mrcoef_leaf[pdcmnt]<<" rm_leaf: "<<rm_leaf<<" leafc: "<<folc<<" gpp: "<<gpp<<" npp0: "<<npp<<" gpr: "<<gpr<<endl;

   // Determine Aboveground plant respiration (veg.abvgrdResp)

   abvgrndResp = gpr - rootResp;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tveg60::updateFoliage( const int& pdcmnt,
                            const double& leafc,
                            const double& vegc,
							const double& eet_mature )
{

  double trans_mature;
  
  alleaf = leafmxc[pdcmnt]
           /(1.0 + kleafc[pdcmnt] * exp( cov[pdcmnt] * vegc ));

  lai = sla[pdcmnt] * leafc;//modified by cgs2014

  foliage = alleaf / leafmxc[pdcmnt];

  /*
  //below added by cgs to recalculate foliage (normalized fpc)
  fpc = 1.0 - exp( -0.5 * lai );
  double xx=leafc/leafmxc[pdcmnt];
  double laimax = leafmxc[pdcmnt]*sla[pdcmnt];
  double fpcmax = 1.0 - exp( -0.5 * laimax );
  foliage = fpc/fpcmax;
  if (foliage >1.0) foliage =1.0;
  */
  fpc = 1.0;

  trans_mature = eet_mature * proptrans[pdcmnt];
  
  transpiration = trans_mature * foliage;
  
};

