/* *************************************************************
****************************************************************
HUMANACT604.CPP - describes human disturbances to natural 
                  ecosystems

Modifications:

20030513 - DWK created by modifying humnact50a.cpp
20030513 - DWK changed Humanact50 to Humanact51
20031016 - DWK replaced char ecd[MAXFNAME] with string ecd to
           getecd()
20031016 - DWK added inheritance of ProcessXML()
20031016 - DWK deleted private functions getvalue() and readline()
20031019 - DWK changed inheritance of ProcessXML to ProcessXML51
20040707 - DWK changed Humanact51:: to Humanact60::
20040707 - DWK changed inheritance of ProcessXML51 to 
           ProcessXML60 in Humanact60()
20040716 - DWK changed double slashpar[MAXCMNT] to 
           double slashpar in functions
20040716 - DWK changed double vconvert[MAXCMNT] to 
           double vconvert in functions
20040716 - DWK changed public double prod10par[MAXCMNT] to 
           double prod10par in functions
20040716 - DWK changed public double prod100par[MAXCMNT] to 
           double prod100par in functions
20040716 - DWK changed public double sconvert[MAXCMNT] to 
           double sconvert in functions
20040716 - DWK changed public double vrespar[MAXCMNT] to 
           double vrespar in functions
20051117 - DWK added include humnact602.h and standard includes
20051112 - DWK added public function createWoodProducts()
20051202 - DWK added resetMonthlyFluxes()
20051202 - DWK added resetYrFluxes()
20060608 - DWK added resetMonthlyDisturbFluxes() and 
           setFireNDEP()
20070210 - DWK changed include from humnact602.h to humnact603.h
20070518 - DWK added grazing() and setNoGrazing()
20070829 - DWK changed include from humnact603.h to humnact603b.h           
20090128 - DWK changed include from humnact603b.h to humnact604.h
		   
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

#include "humnact604.h"


/* *************************************************************
************************************************************* */

Humanact60::Humanact60() : ProcessXML60()
{

  int i;
  int dm;
  c2n = 54.29;
//  cfall = 0.20;
//  nfall = 0.20;

  state = 0;
  prvstate = 0;

  massbalflg = 1;
  fertflag = 0;
  irrgflag = 0;
  frostflag = 0;
  
  mez = CROPVEG - 1;

  productYear = 0;

  prevPROD1.carbon = ZERO;
  prevPROD1.nitrogen = ZERO;

  prevPROD10.carbon = ZERO;
  prevPROD10.nitrogen = ZERO;

  prevPROD100.carbon = ZERO;
  prevPROD100.nitrogen = ZERO;

  prevCropResidue.carbon = ZERO;
  prevCropResidue.nitrogen = ZERO;

  cropResidue.carbon = ZERO;
  cropResidue.nitrogen = ZERO;

  for( dm = 0; dm < CYCLE; ++dm )
  {  
    initPROD1[dm].carbon = ZERO;
    initPROD1[dm].nitrogen = ZERO;
  }

  for( i = 0; i < 10; ++i )
  {
    initPROD10[i].carbon = ZERO;
    initPROD10[i].nitrogen = ZERO;
  }
  
  for( i = 0; i < 100; ++i )
  {
    initPROD100[i].carbon = ZERO;
    initPROD100[i].nitrogen = ZERO;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Humanact60::omitzero( const double& data1)
{
	double tempx;
	if (data1 ==0.0) tempx = 0.00001;
	else tempx = data1;
	return tempx;

};

void Humanact60::conversion( const int& pdcmnt,
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
                             const double& soilc,
                             const double& soiln,
                             const int& disturbflag,
                             const int& ifwetland,
                             const double& randomnum)
{
	//to recalculate the disturbance parameters.
	//if the parameter values are assigned as -99.0 else omit this calculation
	int nonforest[6] ={1,2,3,7,8,12};
	int index=0;
	for (int i=0;i<6;i++){
		if (pdcmnt ==nonforest[i]){
		index=1;
	    }
	}

	FFLC=0.0;
	FFDC=0.0;
	DWFC=0.0;
	SWFC =0.0;
    double litterc2 = 0.0;
    double maxduffraction=0.0;
    int tempindex=0;
    if (disturbflag == 3) //fire disturbance
    {
    	double tempmdc=mmdc;
		if (mmdc >500) tempmdc =500.0;
	//calculate nonforest floor fuel comsumption
	if (index==1) {
		//modify leaf, stem, and litter combustion completeness

		leafconv = leafconv + (1.0-leafconv) * tempmdc/500.0;
		leafslash = 1.0-leafconv;
		stemconv = stemconv + (1.0-stemconv) * tempmdc/500.0;
		stemslash=1.0-stemconv;
		sconvert = sconvert + (1.0-sconvert) * tempmdc/500.0;
		FFLC = sconvert * litterc;
		if (prod10par==-99.0)
		{
		maxduffraction = 0.125;
		if (ifwetland==1) maxduffraction =0.068;
		}
		litterc2 = maxduffraction * (soilc-litterc);
		if (litterc2 >0.0) FFDC = 1000.0 * exp(-4.3 + 0.710 * log(tempmdc) + 0.671 * log(litterc2/1000.0 * 2.0)) * 0.5;
		if (FFDC >litterc2) FFDC = litterc2;
		DWFC = sconvert * cwdc;
		SWFC = sconvert * deadwoodcarbon;
		prod10par=0.0;
		prod100par=0.0;

	}

    if (index==0 )//woody plants
    {
    	//prod10par==-99.0;
    	if (prod10par==-99.0) //for canada and alaska
    	{
    		tempindex=1; //1: canada and alaska; 0: Cont_us or mexico
    		//leafmortpar=0.99;
    		//stemmortpar=0.99;
    		//rootmortpar=0.99;
    		//leafslash = 0.76;
    		//leafconv = 0.24;
    		//rootslash=1.0;
    		//rootconv=0.0;
    		//stemslash=0.3;
    		//stemconv=0.24;
    		//standdead=0.46;
    		//deadconv=0.6;
    		//deadslash=0.4;
    		//sconvert=0.05;
    		prod10par=0.0;
    		prod100par=0.0;
    		maxduffraction=0.125;//max duff carbon = 15 cm top soil layer. source: van der Werf et al. 2010 and Turetsky et al. 2010
    		if (ifwetland ==1) maxduffraction =0.068;
    		litterc2 = maxduffraction * (soilc-litterc);
    		//cout <<" litterc2: "<<litterc2<<" litterc: "<<litterc<<" soilc: "<<soilc<<endl;
    	}

        //FFFC: forest floor fuel (organic layer) carbon (including litterc and top soil carbon)
		//DWDFC: dead or downed wood fuel carbon (including CWD and standing dead wood)
    	if (litterc>0.0)
    	{
		if (tempindex==1) FFLC=1000.0 * exp(-4.2 + 0.710 * log(tempmdc) + 0.671 * log(litterc/1000.0 * 2.0)) * 0.5; //burnt groud layer fine litter carbon. Source Groot et al. (2009)
		else FFLC=1000.0 * exp(-4.2 + 0.71 * log(tempmdc) + 0.75 * log(litterc/1000.0 * 2.0)) * 0.5;
    	}
    	if (litterc2 >0.0) FFDC=1000.0 * exp(-4.3 + 0.710 * log(tempmdc) + 0.671 * log(litterc2/1000.0 * 2.0)) * 0.5;

		DWFC = 1000.0 * (-0.13 + 0.28 * cwdc/1000.0 * 2.0 + 0.00144 * tempmdc) * 0.5; //downed woody debris. source: Kasischke and Hoy et al. GCB. 2012
		SWFC = 1000.0 *  (-0.131 + 0.27 * deadwoodcarbon/1000.0 * 2.0 + 0.00144 * tempmdc) * 0.5;
		if (FFLC > litterc) FFLC = litterc;
		if (FFLC<0.0) FFLC = 0.0;
		if (FFDC > litterc2) FFDC = litterc2;
		if (FFDC<0.0) FFDC = 0.0;

		if (DWFC > cwdc) DWFC=cwdc;
		if (DWFC<0.0) DWFC=0.0;
		if (SWFC > deadwoodcarbon) SWFC=deadwoodcarbon;
		if (SWFC<0.0) SWFC=0.0;
		if (deadwoodcarbon ==0.0) SWFC=0.0;
		if (litterc2 ==0.0) FFDC=0.0;
		if (cwdc ==0.0) DWFC=0.0;
		if (litterc ==0.0) FFLC=0.0;

    }

    }//end of disturbflag

  /******************************************* below code for uncertainty run
  //modify litter and deadwood parameters
  FFDC = FFDC + FFDC * 0.5 * (randomnum/3.20);
  FFLC = FFLC + FFLC * 0.25 * (randomnum/3.20); //25% uncertainty range
  SWFC = SWFC + SWFC * 0.25 * (randomnum/3.20);
  DWFC = DWFC + DWFC * 0.25 * (randomnum/3.20);
  //cout <<"fflc: "<<FFLC<<" ffdc: "<<FFDC<<" swfc: "<<SWFC<<" DWFC: "<<DWFC<<" random: "<<randomnum<<" litterc2: "<<litterc2<<" tempindex: "<<tempindex<<" soilc: "<<soilc<<" litc: "<<litterc<<endl;
  //double xx = FFLC+SWFC+DWFC;
  //cout <<"burntlitterc: "<<xx<<" litterc: "<<litterc<<" litterc2: "<<litterc2<<endl;
  if (FFLC >litterc) FFLC=litterc;
  if (FFLC<0.0) FFLC = 0.0;
  if (FFDC >litterc2) FFDC=litterc2;
  if (FFDC<0.0) FFDC = 0.0;
  if (SWFC >deadwoodcarbon) SWFC=deadwoodcarbon;
  if (SWFC<0.0) SWFC = 0.0;
  if (DWFC >cwdc) DWFC=cwdc;
  if (DWFC<0.0) DWFC = 0.0;

  stemconv = stemconv + stemconv * 0.5 * (randomnum/3.20);
  stemslash = (1.0-stemconv) * stemslash/omitzero(stemslash+standdead);
  standdead = (1.0-stemconv) * standdead/omitzero(stemslash+standdead);

  leafmortpar= leafmortpar + leafmortpar * 0.25 * (randomnum/3.20);
  stemmortpar= stemmortpar + stemmortpar * 0.25 * (randomnum/3.20);
  rootmortpar= rootmortpar + rootmortpar * 0.25 * (randomnum/3.20);

 */

  //********************************************
  //slash.carbon = slashpar * vegc ;
  //slash.nitrogen = slashpar * (vstrn + vston) ;
  //standing deadwood goes to the slash carbon pool
  deadslashc = deadslash * deadwoodcarbon;
  deadslashn = deadslash * deadwoodnitrogen;
  //allocation of stemc to standing deadwood (only stem is allocated to standing dead, roots and leave are allocated to slash)
  deadwood.carbon = standdead * stemc * stemmortpar;
  deadwood.nitrogen = deadwood.carbon * (vstrn + vston) /omitzero(vegc) ;

  slash.carbon = leafmortpar * leafslash * leafc + stemmortpar * stemslash * stemc + rootmortpar * rootslash * rootc + deadslashc;
  slash.nitrogen = slash.carbon/omitzero(vegc) * (vstrn + vston);
  //sconvrtflx.carbon = deadwoodcarbon * deadconv ;
  //sconvrtflx.nitrogen = ((1.0 - nsretconv[pdcmnt])
                          //* deadwoodnitrogen*deadconv);

  stem2slash = stemmortpar*stemslash*stemc;
  leaf2slash = leafmortpar * leafslash * leafc;
  root2slash= rootmortpar * rootslash * rootc;
  //cout <<" stem2slash: "<<stem2slash<<" leaf2slash: "<<leaf2slash<<" root2slash: "<<root2slash<<endl;

  if (stem2slash<0.0) stem2slash=0.0;
  if (root2slash<0.0) root2slash=0.0;
  if (leaf2slash<0.0) leaf2slash=0.0;


  //cout <<" pdcmnt: " <<pdcmnt <<"slashcarbon: " <<slash.carbon <<"deadwood.carbon: " <<deadwood.carbon<< " deadcarbon2: "<<getdeadwoodc()<<" deadpar " <<deadpar<<endl;
  vconvrtflx.carbon = leafmortpar * leafconv * leafc + stemmortpar * stemconv * stemc + rootmortpar * rootconv * rootc ;
  //sconvrtflx.carbon +=  (sconvert * soilc);
  if (disturbflag == 3) sconvrtflx.carbon= FFLC+FFDC+DWFC;
  else sconvrtflx.carbon=sconvert * litterc;

  if (disturbflag != 3) SWFC = deadwoodcarbon * deadconv;
  convertdeadc=SWFC;
  convrtflx.carbon = vconvrtflx.carbon + sconvrtflx.carbon + SWFC;


  //if (sconvrtflx.carbon >0.0) cout <<"sconvrtflx: "<<sconvrtflx.carbon<<" deadwoodc: "<<deadwood.carbon<<endl;
  vconvrtflx.nitrogen = ((1.0 - nvretconv[pdcmnt]) 
                        * vconvrtflx.carbon/omitzero(vegc)
                        * (vstrn + vston));

  convertdeadn=0.0;
  convertdeadn = deadwoodnitrogen * convertdeadc/omitzero(deadwoodcarbon);
  double convertsoiln = 0.0;
  convertsoiln = soiln* sconvrtflx.carbon/omitzero(soilc) ;

  sconvrtflx.nitrogen = (1.0 - nsretconv[pdcmnt])
                        * convertsoiln;
                                
  convrtflx.nitrogen = vconvrtflx.nitrogen + sconvrtflx.nitrogen + convertdeadn * (1.0 - nsretconv[pdcmnt]);
                               
  nvretent = nvretconv[pdcmnt] * vconvrtflx.carbon/omitzero(vegc) * (vstrn + vston);

  nsretent = nsretconv[pdcmnt] * (convertdeadn+ convertsoiln);
  
  nretent = nvretent + nsretent;
  //if (disturbflag ==3) cout <<"cwdc: "<<cwdc<<" DWFC: "<<DWFC<<" deadwoodcarbon: "<<deadwoodcarbon<<" SWFC: "<<SWFC<<" litterc: "<<litterc<<" FFLC: "<<FFLC<<" litter2: "<<litterc2<<" FFDC: "<<FFDC<<" vconvertc: "<<vconvrtflx.carbon<<" totconvc: "<<vconvrtflx.carbon+sconvrtflx.carbon+SWFC<< endl;

  //if (pdcmnt == 7) cout <<"cmnt: "<<pdcmnt<<" vegc: "<<vegc <<" soilc: "<<soilc<<" stemc: "<<stemc<<" index: "<<index<<" deadwoodcarbon: "<<deadwoodcarbon<<" convertdeadc:  "<<convertdeadc<<" vconvertc: "<<vconvrtflx.carbon<<" sconvertc: "<<sconvrtflx.carbon<<" tot_convertc: "<<convrtflx.carbon<<endl;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::createWoodProducts( const int& pdyr,
                                     const double& vegc,
                                     const double& stemc,
                                     const double& vstrn,
                                     const double& vston )
{
  int i;

  formPROD10.carbon  = prod10par * stemc;
  formPROD10.nitrogen  = prod10par * stemc/omitzero(vegc)* (vstrn + vston);
  formPROD100.carbon = prod100par * stemc;
  formPROD100.nitrogen = prod100par * stemc/omitzero(vegc)* (vstrn + vston);

  if( pdyr < productYear )
  {
    productYear += pdyr;
  }
  else { productYear = pdyr; }

  i = productYear%10;
  initPROD10[i].carbon = formPROD10.carbon;
  initPROD10[i].nitrogen = formPROD10.nitrogen;

  i = productYear%100;
  initPROD100[i].carbon = formPROD100.carbon;
  initPROD100[i].nitrogen = formPROD100.nitrogen;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::decayProducts( void )
{

  int i;
  int dm;
  double yrtrash;

  PROD1decay.carbon  = ZERO;
  PROD1decay.nitrogen  = ZERO;

  for( dm = 0; dm < CYCLE; ++dm )
  {
    PROD1decay.carbon += initPROD1[dm].carbon / (double) CYCLE;

    PROD1decay.nitrogen += initPROD1[dm].nitrogen 
                          / (double) CYCLE;
  }


  PROD10decay.carbon  = ZERO;
  PROD10decay.nitrogen  = ZERO;


  for( i = 0; i < 10; ++i )
  {
    yrtrash = initPROD10[i].carbon * 0.10 / (double) CYCLE;
    PROD10decay.carbon += yrtrash;
    yrtrash = initPROD10[i].nitrogen * 0.10 / (double) CYCLE;
    PROD10decay.nitrogen += yrtrash;
  }


  PROD100decay.carbon = ZERO;
  PROD100decay.nitrogen = ZERO;
  
  for( i = 0; i < 100; ++i )
  {
    yrtrash = initPROD100[i].carbon * 0.01 / (double) CYCLE;
    PROD100decay.carbon += yrtrash;
    yrtrash = initPROD100[i].nitrogen * 0.01 / (double) CYCLE;
    PROD100decay.nitrogen += yrtrash;
  }

  TOTPRODdecay.carbon = PROD1decay.carbon
                        + PROD10decay.carbon
                        + PROD100decay.carbon;

  TOTPRODdecay.nitrogen = PROD1decay.nitrogen
                          + PROD10decay.nitrogen
                          + PROD100decay.nitrogen;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::frostDamage( const double& vegc,
                              const double& vstrn, 
                              const double& vston )
{ 
  frostflag = 1;
  
  // Assume all crop biomass becomes stubble
  
  stubble.carbon = vegc;
  stubble.nitrogen = (vstrn+vston);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::getecd( ofstream& rflog1 )
{
  string ecd;

#ifdef PMODE
  *fgo >> ecd;
#else

  cout << "Enter name of the data file (.ECD) with agricultural ";
  cout << "parameter values: " << endl;
  
  cin >> ecd;
#endif

  rflog1 << "Enter name of the soil data file (.ECD) with agricultural ";
  rflog1 << "parameter values: " << ecd << endl;

  getecd( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::getecd( const string& ecd )
{
  ifstream infile;
  int dcmnt;
  int comtype;

  infile.open( ecd.c_str(), ios::in );

  if( !infile )
  {
    cerr << "\nCannot open " << ecd;
    cerr << " for agriculture ECD input" << endl;
    
    exit( -1 );
  }
  
  getXMLrootNode( infile, "agECD" );

  for( dcmnt = 1; dcmnt < MAXCMNT; ++dcmnt )
  {
    comtype = getXMLcommunityNode( infile, "agECD" );

    if( comtype >= MAXCMNT )
    {
      cerr << endl << "comtype is >= MAXCMNT" << endl;
      cerr << "comtype cannot be greater than " << (MAXCMNT-1) << endl;
      cerr << " in agECD" << endl;
      
      exit( -1 );
    }
    

    nvretconv[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "nvretconv",
                                                comtype );

    nsretconv[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "nsretconv",
                                                comtype );

    tillfactor[comtype] = getXMLcmntArrayDouble( infile,
                                                 "agECD",
                                                 "tillfactor",
                                                 comtype );

    harvstC[comtype] = getXMLcmntArrayDouble( infile,
                                              "agECD",
                                              "harvstC",
                                              comtype );

    harvstN[comtype] = getXMLcmntArrayDouble( infile,
                                              "agECD",
                                              "harvstN",
                                              comtype );

    residueC[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "residueC",
                                               comtype );

    residueN[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "residueN",
                                               comtype );

    cropseedC[comtype] = getXMLcmntArrayDouble( infile,
                                                "agECD",
                                                "cropseedC",
                                                comtype );

    cropseedSTRN[comtype] = getXMLcmntArrayDouble( infile,
                                                   "agECD",
                                                   "cropseedSTRN",
                                                   comtype );

    cropseedSTON[comtype] = getXMLcmntArrayDouble( infile,
                                               "agECD",
                                               "cropseedSTON",
                                               comtype );
    endXMLcommunityNode( infile );
  }

  if( dcmnt < MAXCMNT )
  {
    cerr << endl << " Parameters found for only " << dcmnt;
    cerr << " community types out of a maximum of ";
    cerr << (MAXCMNT-1) << " types in agECD" << endl;
    
    exit( -1 );
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::grazing( const double& vegc,
                          const double& strn )
{
  // Assume 5 percent of vegetation biomass is consumed each
  //   month

  forage.carbon = vegc * 0.05;
  forage.nitrogen = strn * 0.05;

  // Assume 83 percent of forage carbon is respired and 50
  //   percent of nitrogen is mineralized in livestock
  //   metabolism

  animalresp = forage.carbon * 0.83;
  urine = forage.nitrogen * 0.50;

  // Assume remaineder of forage is returned to soil as manure

  manure.carbon = forage.carbon * 0.17;
  manure.nitrogen = forage.nitrogen * 0.50;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::harvest( const int& pdm,
                          const double& vegc,
                          const double& vstrn, 
                          const double& vston )
{
//  int kdm;
//  int idm;

  double initresidueC;
  double initresidueN;

//  double residueCflux;
//  double residueNflux;
//  double residueNretent;

  // Determine crop production

  cropprod.carbon = vegc * harvstC[cmnt];

  if( cropprod.carbon < ZERO )
  {
    cropprod.carbon = ZERO;
  }

  cropprod.nitrogen = (vstrn+vston) * harvstN[cmnt];

  if( cropprod.nitrogen < ZERO )
  {
    cropprod.nitrogen = ZERO;
  }

  initPROD1[pdm].carbon  = cropprod.carbon;
  initPROD1[pdm].nitrogen  = cropprod.nitrogen;


  // Determine amount of carbon and nitrogen left
  // in crop residue

   initresidueC = vegc - cropprod.carbon;

  if( initresidueC < ZERO )
  {
    initresidueC = ZERO;
  }

  initresidueN = (vstrn+vston) - cropprod.nitrogen;

  if( initresidueN < ZERO )
  {
    initresidueN = ZERO;
  }

  // Determine amount of carbon and nitrogen in
  // crop residue that will be lost from ecosystem

  formCropResidue.carbon = initresidueC
                           * residueC[cmnt];

  formCropResidue.nitrogen = initresidueN
                             * residueN[cmnt];


  // Determine amount of stubble carbon and nitrogen
  // added to soils

  stubble.carbon = initresidueC
                   * ( 1.000000 - residueC[cmnt]);

  stubble.nitrogen = initresidueN
                     * ( 1.000000 - residueN[cmnt]);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::resetMonthlyDisturbFluxes( void )
{
  // Initialize disturbance carbon fluxes to zero

  convrtflx.carbon = ZERO;
  vconvrtflx.carbon = ZERO;
  sconvrtflx.carbon = ZERO;
  slash.carbon = ZERO;
  deadwood.carbon = ZERO;

  // Initialize nitrogen fluxes during conversion to zero

  convrtflx.nitrogen = ZERO;
  vconvrtflx.nitrogen = ZERO;
  sconvrtflx.nitrogen = ZERO;
  slash.nitrogen = ZERO;
  deadwood.nitrogen = ZERO;

  nretent = ZERO;
  nvretent = ZERO;
  nsretent = ZERO;

  cropResidueFlux.carbon = ZERO;
  cropResidueFlux.nitrogen = ZERO;

  FFLC=0.0; //forest floor litter C
  FFDC=0.0; //forest floor duff C
  DWDFC=0.0;
  DWFC=0.0; //downed wood (CWD) fuel c
  SWFC=0.0; //standing wood fuel c
  deadslashc=0.0;
  deadslashn=0.0;
  convertdeadc=0.0;
  convertdeadn=0.0;
  stem2slash=0.0;
  leaf2slash=0.0;
  root2slash=0.0;
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::resetMonthlyFluxes( void )
{
  // Reset monthly fluxes to zero
  
  fertn = ZERO;
  
  irrigate = ZERO;

  // Initialize stubble to zero
  
  stubble.carbon = ZERO;
  stubble.nitrogen = ZERO;

  cflux = ZERO;

  
  // Initialize carbon fluxes related to formation and 
  //   decomposition of products and crop residue to zero

  cropprod.carbon = ZERO;
  cropprod.nitrogen = ZERO;

  formCropResidue.carbon = ZERO;
  formCropResidue.nitrogen = ZERO;

  PROD1decay.carbon = ZERO;
  PROD1decay.nitrogen = ZERO;

  formPROD10.carbon = ZERO;
  formPROD10.nitrogen = ZERO;

  PROD10decay.carbon = ZERO;
  PROD10decay.nitrogen = ZERO;

  formPROD100.carbon = ZERO;
  formPROD100.nitrogen = ZERO;

  PROD100decay.carbon = ZERO;
  PROD100decay.nitrogen = ZERO;

  formTOTPROD.carbon = ZERO;
  formTOTPROD.nitrogen = ZERO;

  TOTPRODdecay.carbon = ZERO;
  TOTPRODdecay.nitrogen = ZERO;

	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::resetPROD( void )
{
  int i;
  int dm;

  prevPROD1.carbon = ZERO;
  prevPROD1.nitrogen = ZERO;

  prevPROD10.carbon = ZERO;
  prevPROD10.nitrogen = ZERO;

  prevPROD100.carbon = ZERO;
  prevPROD100.nitrogen = ZERO;

  prevCropResidue.carbon = ZERO;
  prevCropResidue.nitrogen = ZERO;

  cropResidue.carbon = ZERO;
  cropResidue.nitrogen = ZERO;
 
  for( dm = 0; dm < CYCLE; ++dm )
  {
    initPROD1[dm].carbon = ZERO;
    initPROD1[dm].nitrogen = ZERO;
  }

  for( i = 0; i < 10; ++i )
  {
    initPROD10[i].carbon = ZERO;
    initPROD10[i].nitrogen = ZERO;
  }
  for( i = 0; i < 100; ++i )
  {
    initPROD100[i].carbon = ZERO;
    initPROD100[i].nitrogen = ZERO;
  }

  PROD1.carbon = ZERO;
  PROD1.nitrogen = ZERO;

  PROD10.carbon = ZERO;
  PROD10.nitrogen = ZERO;

  PROD100.carbon = ZERO;
  PROD100.nitrogen = ZERO;

  TOTPROD.carbon = ZERO;
  TOTPROD.nitrogen = ZERO;

  PROD1decay.carbon = ZERO;
  PROD1decay.nitrogen = ZERO;

  formPROD10.carbon  = ZERO;
  formPROD10.nitrogen  = ZERO;

  PROD10decay.carbon  = ZERO;
  PROD10decay.nitrogen  = ZERO;

  formPROD100.carbon = ZERO;
  formPROD100.nitrogen = ZERO;

  PROD100decay.carbon = ZERO;
  PROD100decay.nitrogen = ZERO;

  formTOTPROD.carbon = ZERO;
  formTOTPROD.nitrogen = ZERO;

  TOTPRODdecay.carbon = ZERO;
  TOTPRODdecay.nitrogen = ZERO;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::resetYrFluxes( void )
{
  // Reset annual fluxes and summary variables to zero

  yrfertn = ZERO;

  yrirrig = ZERO;

  yrstubC = ZERO;
  yrstubN = ZERO;

  // Annual carbon fluxes from agricultural conversion

  yrconvrtC = ZERO;
  yrvconvrtC = ZERO;
  yrsconvrtC = ZERO;
  yrslashC = ZERO;
  yrcflux = ZERO;

 // Annual nitrogen fluxes from agricultural conversion

  yrconvrtN = ZERO;
  yrvconvrtN = ZERO;
  yrsconvrtN = ZERO;
  yrslashN = ZERO;
  yrnrent = ZERO;
  yrnvrent = ZERO;
  yrnsrent = ZERO;


  // Annual carbon and nitrogen fluxes in the formation of 
  //   agricultural products

  yrformPROD1C   = ZERO;
  yrformPROD1N   = ZERO;


 // Annual carbon and nitrogen flux from crop residue formation

  yrformResidueC = ZERO;
  yrformResidueN = ZERO;

  // Annual carbon and nitrogen fluxes in the decomposition of 
  //   agricultural products

  yrdecayPROD1C   = ZERO;
  yrdecayPROD1N   = ZERO;

 // Annual carbon and nitrogen fluxes from burning crop residue

  yrfluxResidueC = ZERO;
  yrfluxResidueN = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   formation of 10-year wood products

  yrformPROD10C  = ZERO;
  yrformPROD10N  = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   decomposition of 10-year wood products

  yrdecayPROD10C  = ZERO;
  yrdecayPROD10N  = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   formation of 100-year wood products

  yrformPROD100C = ZERO;
  yrformPROD100N = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   decomposition of 100-year wood products

  yrdecayPROD100C = ZERO;
  yrdecayPROD100N = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   formation of all agricultural and wood products

  yrformTOTPRODC = ZERO;
  yrformTOTPRODN = ZERO;

  // Annual carbon and nitrogen fluxes resulting from 
  //   decomposition of all agricultural and wood products

  yrdecayTOTPRODC = ZERO;
  yrdecayTOTPRODN = ZERO;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::setAgricFlags( ofstream& rflog1 )
{

#ifdef PMODE
  *fgo >> tillflag;
#else
  cout << "Are agricultural soils tilled?" << endl;
  cout << "Enter 0 for no:" << endl;
  cout << "Enter 1 for yes:" << endl;
  
  cin >> tillflag;
#endif

  rflog1 << "Are agricultural soils tilled?" << endl;
  rflog1 << "Enter 0 for no:" << endl;
  rflog1 << "Enter 1 for yes:" << endl;
  rflog1 << tillflag << endl << endl;

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Humanact60::setFireNDEP( void )
{
  double newndep;
  double frftemp;

  // Determine firendep from total nitrogen volatilized during 
  //   fire for reintroduction to the soil (Note: 
  //   vconvrtflx.nitrogen and sconvrtflx.nitrogen are monthly
  //   fluxes that occur over a time period of a year
  
  newndep = vconvrtflx.nitrogen + sconvrtflx.nitrogen;

  
  // Determine the potential number of years until next fire 
  //   event (i.e. Fire Return Interval or FRI)
  
  frftemp = FRI;

  if( frftemp > MAXFRI ) { frftemp = MAXFRI; }

  
  // Assume that (1/FRI) of newndep is deposited each month
  
  firendep = newndep / (frftemp-1);
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::setNoCrops( const int& pdm )
{
  cropprod.carbon = ZERO;
  cropprod.nitrogen = ZERO;

  formCropResidue.carbon = ZERO;
  formCropResidue.nitrogen = ZERO;

  if( frostflag != 1 )
  {
    stubble.carbon = ZERO;
    stubble.nitrogen = ZERO;
  }
  
  frostflag = 0;
  
  initPROD1[pdm].carbon = ZERO;
  initPROD1[pdm].nitrogen = ZERO;
  
  	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::setNoGrazing( void )
{
  forage.carbon = ZERO;
  forage.nitrogen = ZERO;

  manure.carbon = ZERO;
  manure.nitrogen = ZERO;

  animalresp = ZERO;

  urine = ZERO;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::setNoWoodProducts( const int& pdyr )
{
  int i;
  	
  formPROD10.carbon = ZERO;
  formPROD10.nitrogen = ZERO;

  formPROD100.carbon = ZERO;
  formPROD100.nitrogen = ZERO;
  
  if( pdyr < productYear )
  {
    productYear += pdyr;
  }
  else { productYear = pdyr; }

  i = productYear%10;
  initPROD10[i].carbon = ZERO;
  initPROD10[i].nitrogen = ZERO;

  i = productYear%100;
  initPROD100[i].carbon = ZERO;
  initPROD100[i].nitrogen = ZERO;
	
};


/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::updateCropResidue( void )
{

  cropResidue.carbon = prevCropResidue.carbon
                       + formCropResidue.carbon
                       - cropResidueFlux.carbon;

  cropResidue.nitrogen = prevCropResidue.nitrogen
                         + formCropResidue.nitrogen
                         - cropResidueFlux.nitrogen;

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::updateCropResidueFluxes( void )
{
  int kdm;

  cropResidueFlux.carbon = ZERO;
  cropResidueFlux.nitrogen = ZERO;

  if ( 1 == state )
  {
    for( kdm = 0; kdm < CYCLE; ++kdm )
    {
      if( harvstC[cmnt] > 0.000001 )
      {
        cropResidueFlux.carbon += initPROD1[kdm].carbon
                                  * (1.000000 - harvstC[cmnt])
                                  / harvstC[cmnt]
                                  * residueC[cmnt]
                                  / (double) CYCLE;
      }
      else { cropResidueFlux.carbon = ZERO; }

      if( harvstN[cmnt] > 0.000001 )
      {
        cropResidueFlux.nitrogen += initPROD1[kdm].nitrogen
                                    * (1.000000 - harvstN[cmnt])
                                    / harvstN[cmnt]
                                    * residueN[cmnt]
                                    / (double) CYCLE;
      }
      else { cropResidueFlux.nitrogen = ZERO; }
    }
  }

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::updateProducts( void )
{
  int i;

  // Carbon in Products

  PROD1.carbon = prevPROD1.carbon 
                 + cropprod.carbon 
                 - PROD1decay.carbon;
  
  if( PROD1.carbon < 0.1 )
  {
    PROD1.carbon = ZERO;

    for(i = 0; i < CYCLE; ++i )
    {
      initPROD1[i].carbon = ZERO;
      initPROD1[i].nitrogen = ZERO;
    }
  }

  PROD10.carbon = prevPROD10.carbon
                  + formPROD10.carbon
                  - PROD10decay.carbon;
  
  if( PROD10.carbon < 0.1 )
  {
    PROD10.carbon = ZERO;

    for(i = 0; i < 10; ++i )
    {
      initPROD10[i].carbon = ZERO;
      initPROD10[i].nitrogen = ZERO;
    }
  }

  PROD100.carbon = prevPROD100.carbon
                   + formPROD100.carbon
                   - PROD100decay.carbon;
  
  if( PROD100.carbon < 0.1 )
  {
    PROD100.carbon = ZERO;
    for(i = 0; i < 100; ++i )
    {
      initPROD100[i].carbon = ZERO;
      initPROD100[i].nitrogen = ZERO;
    }
  }

  TOTPROD.carbon = PROD1.carbon
                   + PROD10.carbon
                   + PROD100.carbon;

  // Nitrogen in Products

  PROD1.nitrogen = prevPROD1.nitrogen
                   + cropprod.nitrogen
                   - PROD1decay.nitrogen;
  
  if( PROD1.nitrogen < 0.000001 )
  {
    PROD1.nitrogen = ZERO;
  }

  PROD10.nitrogen = prevPROD10.nitrogen
                    + formPROD10.nitrogen
                    - PROD10decay.nitrogen;
  
  if( PROD10.nitrogen < 0.000001 )
  {
    PROD10.nitrogen = ZERO;
  }

  PROD100.nitrogen = prevPROD100.nitrogen
                     + formPROD100.nitrogen
                     - PROD100decay.nitrogen;
  
  if( PROD100.nitrogen < 0.000001 )
  {
    PROD100.nitrogen = ZERO;
  }

  TOTPROD.nitrogen = PROD1.nitrogen
                     + PROD10.nitrogen
                     + PROD100.nitrogen;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Humanact60::updateTotalProductFormation( void )
{
  // Carbon in Total Product Formation

  formTOTPROD.carbon = cropprod.carbon
                       + formPROD10.carbon
                       + formPROD100.carbon;

  // Nitrogen in Total Product Formation

  formTOTPROD.nitrogen = cropprod.nitrogen
                         + formPROD10.nitrogen
                         + formPROD100.nitrogen;


};
double Humanact60::sampleNormal(void)
{
	double u=((double)rand()/(RAND_MAX))*2-1;
	double v=((double)rand()/(RAND_MAX))*2-1;
	double r=u*u+v*v;
	if (r==0||r>1) return sampleNormal();
	double c=sqrt(-2*log(r)/r);
	return u*c;
};
double Humanact60::normalDistribution(double x)
{
  if(x<-10.)return 0.;
  if(x>10.)return 1.;
  // number of steps
  int N=2000;
  // range of integration
  double a=0,b=x;
  // local variables
  double s,h,sum=0.;
  // inialise the variables
  h=(b-a)/N;
  // add in the first few terms
  sum = sum + exp(-a*a/2.) + 4.*exp(-(a+h)*(a+h)/2.);
  // and the last one
  sum = sum + exp(-b*b/2.);
  // loop over terms 2 up to N-1
  for(int i=1;i<N/2;i++)
  {
    s = a + 2*i*h;
    sum = sum + 2.*exp(-s*s/2.);
    s = s + h;
    sum = sum + 4.*exp(-s*s/2.);
  }
  // complete the integral
  sum = 0.5 + h*sum/3./sqrt(8.*atan(1.));
  // return result
  return sum;
};

