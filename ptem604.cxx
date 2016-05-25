/* **************************************************************
PTEM604.CXX - Parallel Version of Extrapolation version of the Terrestrial Ecosystem
               Model Version 6.04
*****************************************************************


Modifications:

20040825 - DWK created by modifying xtem50b3.cxx
20040825 - DWK deleted global const MAXFNAME
20040825 - DWK changed include from tclmlist50b3.h to 
           tclmlist60.h
20040825 - DWK changed include from tlulclist431a.h to 
           tlclclist433.h
20040825 - DWK changed include from temlist50b3.h to temlist60.h
20040825 - DWK changed include from telm50b3.h to telm60.h
20040825 - DWK changed ClmList50 gridclm to ClmList60 gridclm
20040825 - DWK changed TEMlist50 gridtem to TEMlist60 gridtem
20040825 - DWK added global const int MAXCHRTS = 5
20040825 - DWK changed char predmap[MAXPRED][9] to 
           string predmap[MAXPRED]
20040825 - DWK added telmnt[0].clm.setO3Flags(), 
           telmnt[0].clm.initO3(), telmnt[0].clm.setNdepFlags() 
           and telmnt[0].clm.initNdep() to initCLM()
20040827 - DWK added string o3name and string ndepname to
           initializeTCLMGridCell()
20040827 - DWK added global FILE* ifo3 and FILE* ifndep
20040827 - DWK added include tndepdat60.h
20040827 - DWK added global Clmdata50 o3dat[MAXRTIME] and
           Ndepdata60 ndepdat[MAXRTIME]
20040828 - DWK changed global const MAXPRED = 124 to 
           MAXPRED = 240
20040828 - DWK renamed initializeTCLMGridCell() to be
           initializeCLMGridCell()
20040828 - DWK added global FILE* fnumchrts
20040828 - DWK added FILE* fnchrts to initializeLCLUCGridCell() 
20040828 - DWK changed include from tlulclist433.h to
           tlulclist60.h   
20040829 - DWK changed Lulcdata43 lulcdat[MAXRTIME]  
           to Lulcdata60 lulcdat[MAXRTIME]
20040830 - DWK added MaxCohortdata60 mxcohrtdat[MAXRTIME];
20040830 - DWK deleted FILE* fpotveg
20040830 - DWK moved initializeLCLUCGridCell() to initLCLUC()
20040830 - DWK changed include from tclmdat50b3.h to
           tclmdat60.h
20040830 - DWK changed global Clmdata50 variables to 
           Clmdat60 variables
20040925 - DWK changed global Lulcdata60 lulcdat[MAXRTIME]  
           to Lulcdata60 lulcdat[MAXRTIME][MAXCHRTS]
20041002 - DWK changed global const int MAXRTIME = 600 to
           MAXRTIME = 900;
20041003 - DWK changed global const int MAXPRED = 240 to
           MAXPRED = 241
20041003 - DWK added global int spinpdyrs 
20050409 - DWK changed include from temlist60.h to temlist601.h
20050409 - DWK changed include from telm60.h to telm601.h
20050622 - DWK changed include from tlulclist60.h to
           tlulclist601.h
20050622 - DWK changed include from telm601.h to telm601a.h
20050722 - DWK changed include from tlulclist601.h to 
           tlulclist601a.h
20051117 - DWK added include temconsts602.hpp
20051117 - DWK changed include from tclmlist60.h to tclmlist602.h
20051117 - DWK changed include from tclmdat60.h to tclmdat602.h
20051117 - DWK changed include from tndepdat60.h to tndepdat602.h
20051117 - DWK changed include from tlulclist601a.h to 
           tlulclist602.h
20051117 - DWK changed include from temlist601.h to temlist602.h
20051117 - DWK changed include from elmnt425.h to elmnt602.h
20051117 - DWK changed include from latdat425.h to latdat602.h
20051117 - DWK changed include from telm601a.h to telm602.h
20051117 - DWK changed Elmnt elmnt to Elmnt60 elmnt
20051119 - DWK deleted Clmlist60 gridclm and associated includes
20051121 - DWK deleted LULCList60 gridlulc and associated 
           includes
20051121 - DWK deleted TEMlist60 gridtem and associated includes 
20051121 - DWK added global vector<string> clmpredmap( NUMATMS )          
20051121 - DWK added const int& posspred and 
           vector<string>& predmap to function call of
           askpred()
20051121 - DWK deleted int spred from askpred()
20051121 - DWK changed global vector<string> predmap( MAXPRED )
           to vector<string> tempredmap( NUMTEM )
20051121 - DWK added co2name to initializeCLMGridCell()
20051121 - DWK added global ifstream ifco2
20051121 - DWK deleted const int& apred from void setTEMPred()
20060607 - DWK added include tco2dat437.h
20060607 - DWK added global CO2data43 co2dat[MAXRTIME]
20060909 - DWK changed include from telm602.h to telm603.h
20070514 - DWK changed include from telm603.h to telm603a.h
20070829 - DWK changed include from telm603a.h to telm603b.h
20070830 - DWK changed include from telm603b.h to telm603c.h
                        
20130805 - FMY add parallel (mpi) to xtem604.cpp: essentially dividing gridcells over processes
             (so, it's not real parallel computing)


*****************************************************************
************************************************************** */
//#define PMODE     // parallel run mode switch (FMY: Aug 2013)

#define ANSI_CPP

//#define BORLAND_CPP

#include<cstdio>

  using std::fopen;
  using std::fclose;
  using std::printf;
  using std::sprintf;
  using std::fscanf;
  using std::FILE;

#include<iostream>

  using std::cout;
  using std::cin;
  using std::ios;
  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setiosflags;
  using std::setw;
  using std::setprecision;

#include<cstdlib>

  using std::exit;
  
#include<cmath>

#include<vector>

  using std::vector;
  
#include<cctype>

  using std::toupper;

#include<cstring>
  
#include<string>
  
  using std::string;
  
#include<sstream>

  using std::ostringstream;

#ifdef ANSI_CPP

  #include<ctime>

  using std::time_t;
  using std::ctime;

#endif

#ifdef BORLAND_CPP

  #include<time>

  using std::time_t;
  using std::ctime;

#endif

#include "temconsts602.hpp"   // Global constants
#include "tclmdat602.h"       // Clmdata60 class
#include "tco2dat437.h"       // CO2data43 class
#include "tndepdat602.h"      // Ndepdata60 class
#include "elmnt602.h"         // Elmnt60 Class
#include "latdat602.h"        // Latdata60 class
#include "telm604.h"         // TEMelmnt60 Class

int askpred( const vector<string>& pvarname, 
             const int& posspred, 
             const vector<string>& predmap );

void initCLM( void );
// void initElmntBatch( void );
void initializeCLMGridCell( void );
void initializeLCLUCGridCell( void );
void initializeTEMGridCell( void );
void initLCLUC( void );
void initTEM( void );
void setClmPred( void );
void setGridParameters( void );
void setInTEMState( void );
void setOutTEMState( const int& nthproc );   // output to different files for paralleling (fmy: Aug 2013)
void setRunMode( const int& nthproc );       // log to different files for paralleling (fmy: Aug 2013)
void setRunTime( void );
void setTEMPred( void );
long skipRecords( void );
void skipGridCells(const int &grdcnt);

void updateTCLMGridCell( const int& pdyr );
void updateTLCLUCGridCell( const int& pdyr );
void updateTTEMGridCell( const int& pdyr);


ofstream flog1;

Elmnt60 elmnt;
TEMelmnt60 telmnt[MAXGRID];

int equil;
int RTIME;

long mxnumgrid;

int spinflag;
int spinoutfg;
int spinoutyrs;
int numspin;
int spintime;
int totsptime;
int transtime;

int temflag;
int istateflag;
int istateyear;
int ostateflag;
int ostateyear;

vector<string> clmpredmap( NUMATMS );
vector<string> tempredmap( NUMTEM );
    
int fatalerr;
int end1;
int icount;
int glob_count;

Latdata60 latdat;

Clmdata60 girrdat[MAXRTIME];
Clmdata60 cldsdat[MAXRTIME];
Clmdata60 nirrdat[MAXRTIME];
Clmdata60 pardat[MAXRTIME];
Clmdata60 tairdat[MAXRTIME];
Clmdata60 precdat[MAXRTIME];

CO2data43 co2dat[MAXRTIME+1];
Clmdata60 o3dat[MAXRTIME];
Ndepdata60 ndepdat[MAXRTIME];

MaxCohortdata43 mxcohrtdat[MAXRTIME];
Lulcdata44 lulcdat[MAXRTIME][MAXCHRTS];  

FILE* flonlat;
  
FILE* ifgirr;
FILE* ifnirr;
FILE* ifpar;
FILE* ifclds;
FILE* iftair;
FILE* ifprec;
ifstream ifco2;
FILE* ifo3;
FILE* ifndep;

FILE* ifnumchrts;
FILE* iflulc;

FILE* fstxt;
FILE* felev;

ifstream ifstate;  // Use TEMstate from a specified year
ofstream ofstate;  // Save TEMstate for a specified year
ofstream fclmpred[NUMATMS];
ofstream ftempred[NUMTEM];
ofstream outfile;
ofstream outfile1;
ofstream outfile2;
ofstream outfile3;
ofstream outfile4;
ofstream outfile5;
ofstream outfile6;
ofstream outfile7;
ofstream outfile8;
ofstream outfile9;
ofstream outfile10;
ofstream outfile11;

//added by cgs to skip some input data
//int skipflag=0;
int ozonestartyr;
int ndepstartyr;
int co2startyr;
int clmstartyr;
int clmendyr;
int ozoneendyr;
int ndependyr;
int co2endyr;
int lulcstartyr;
int lulcendyr;
ifstream fclmfile;
ifstream fndepfile;
ifstream fco2file;
ifstream fozonefile;
int readspinfile;
int readclmyr;
int readndepyr;
int readco2yr;
int readozoneyr;
int skipclm[MAXRTIME];
int skipndep[MAXRTIME];
int skipozone[MAXRTIME];
int skipco2[MAXRTIME];

int assignCO2 = 0;

// **************fmy: parallalization *************************************/
	int noprocs = 1;
	int nthproc = 0;

#ifdef PMODE
	#include <mpi.h>
	string gofile;     // go filename as input (no more stdin)
	ifstream fgo;

	vector<int> pmxnumgrid;    // numbers of grid cells for individual process
	vector<int> pgrdcnt0;      // starting count of grid cells for individual process

// **************fmy: parallalization *************************************/
#endif


/* *************************************************************
**********************START MAIN PROGRAM************************
************************************************************* */

int main(int argc, char *argv[])
{

  time_t stime;
  stime = time( 0 );

#ifdef PMODE
  // **************fmy: parallalization *************************************/
  //initialization
  MPI::Init();
  noprocs = MPI::COMM_WORLD.Get_size();
  nthproc = MPI::COMM_WORLD.Get_rank();

  cout << "TEM runs on PROCESS: "<<nthproc<< "- Simulation started @ "<< ctime( &stime )<<endl;

  // go file as input, not stdin anymore
  gofile = "";
  if (argc==1) {
	  gofile = "xtem604_demo.go";   // this is the default go file name
  } else if (argc==2) {
	  gofile = argv[1];             // user-defined go file name
  }

  // So, we need to make sure ALL 'stdin' from go file are exactly same for all input modules, in which pointer 'fgo'
  // are defined, and here we MUST initiate all of them and refer to the same file input stream
  elmnt.fgo=&fgo;

  telmnt[0].clm.fgo=&fgo;

  telmnt[0].tem.fgo=&fgo;
  telmnt[0].tem.veg.fgo=&fgo;
  telmnt[0].tem.soil.fgo=&fgo;
  telmnt[0].tem.soil.stm.fgo=&fgo;
  telmnt[0].tem.microbe.fgo=&fgo;
  telmnt[0].tem.ag.fgo=&fgo;

  telmnt[0].lcluc.fgo=&fgo;


  // Here, we let every process read the same go file; but,
  // IDEALLY, it should be read-in by 1 process and then pass to the rest processes -TODO
  fgo.open(gofile.c_str(),ios::in );
  bool isOpen = fgo.is_open();
  if ( !isOpen ) {
    	cout << "\nCannot open GO FILE: " << gofile << "! will exit\n" ;
    	exit( -1 );
  }
  // **************fmy: parallalization *************************************/
#endif

  int dyr;
  int i;

//  int itype;
  
  long grdcnt;


  // Run model in steady state mode or dynamic mode? 

  setRunMode(nthproc);


  // Determine time frame of simulation
  
  setRunTime();

  
  // Specify number of grid cells in the study region and determine 
  // if spatially explicit data is referenced by column and row
  // or by longitude and latitude
  
  setGridParameters();

#ifdef PMODE
  // **************fmy: parallalization *************************************/
    // dividing grid cells over processes
  int numskip0 = elmnt.numskip;    // the skipped number of cells in 'go file'
  int pgridnos = floor((mxnumgrid-numskip0)/noprocs);
  int pgridres= (mxnumgrid-numskip0) - pgridnos*noprocs;
  int pgrd0 = numskip0;            // the starting cell no. of grids on each process
  for (int pid=0; pid<noprocs; pid++) {

  	    int grdno = pgridnos;    // the no. of grid cells for each process
		if (pid<pgridres) grdno=pgridnos+1;  // reminder of grids put on the first 'pgridres' processes by 1 more than average
		pmxnumgrid.push_back(grdno);
		if (pid==nthproc)
			cout << "total no. of cells on Process " << pid << ": "<< grdno<<"\n" ;

		if (pid>0) pgrd0=pgrdcnt0.at(pid-1)+pmxnumgrid.at(pid-1);
		pgrdcnt0.push_back(pgrd0);
		if (pid==nthproc)
			cout << "starting no. of cells on Process " << pid << ": "<< pgrd0<<"\n" ;

  }
  // **************fmy: parallalization *************************************/
#endif

  // Initialize climate module

  initCLM();

  // Initialize land cover module

  initLCLUC();
  

  // Initialize TEM 

  initTEM();  


  telmnt[0].col = MISSING;
  telmnt[0].row = MISSING;
  telmnt[0].tem.totyr = -99;

  cout << endl;
  flog1 << endl << endl;
  /*
  elmnt.show( nthproc, flog1,
              telmnt[0].col, 
              telmnt[0].row,
              telmnt[0].tem.totyr, 
              telmnt[0].tem.inittol );
  */

  // Extrapolate TEM across region

#ifdef PMODE
  // **************fmy: parallalization *************************************/

  int grdcnt0 = pgrdcnt0.at(nthproc);                          // starting grid cell for the 'nthproc' process
  int grdcntx = grdcnt0+pmxnumgrid.at(nthproc);                // ending grid cell for the 'nthproc' process

  elmnt.numskip = grdcnt0;  // here, update the orginal input 'numskip' by 'grdcnt0'
  if (grdcnt0>0) elmnt.strtflag = 0;  // whatever previous value, 'strflag' not starting from begining anymore
  skipGridCells(grdcnt0);   // then, skip cells

  grdcnt = grdcnt0;
  while( grdcnt < grdcntx && 0 == fatalerr )   // Grid cell loop for the 'nthproc' process

  // **************fmy: parallalization *************************************/
#else
  //grdcnt = 0;
  grdcnt = elmnt.numskip;
  if (elmnt.strtflag == 0)
  {
	  skipGridCells(grdcnt);   // then, skip cells
  }

  while( grdcnt < mxnumgrid && fatalerr==0 )   // Grid cell loop
#endif
  {
    
// *************************************************************
// Skip to the desired record in the GIS data sets

//  grdcnt = skipRecords();   // fmy: here seems not right - skip only needs once

//*********************************************************** */
#ifdef PMODE
	 //cout <<" nthproc: " <<nthproc<<" grdcnt: "<<grdcnt <<" maxnum: " <<grdcntx<<"mxnumgrid: " <<mxnumgrid<<endl;
#endif
    // Load grid cell climate data into one node of CLM linked list

    dyr = 0;
    updateTCLMGridCell( dyr );


    // Determine number of land use/land cover cohorts in a 
    //   grid cell and load land cover data into cohorts of 
    //   LULC linked list

    updateTLCLUCGridCell( dyr );


    // Initialize TEM to equilibrium conditions for all cohorts
    // using the baseline climate if starting from calibration data 
    // (i.e. istateflag == 0) or read in initial conditions from 
    // temstate file

    initializeTEMGridCell();


    // Begin simulation of transient climate and terrestrial 
    //   ecosystem response
    if( 0 == equil ) //added by cgs
    {
        /*//below code is added for uncertainty analysis. cgs2015
        //if (telmnt[0].col == -98.0 && telmnt[0].row== 55.50) //can forest
    	//if (telmnt[0].col == -104 && telmnt[0].row== 64.25) //can grass
    	//if (telmnt[0].col == -128.5 && telmnt[0].row== 60.25) //ak forest
    	//if (telmnt[0].col == -127.5 && telmnt[0].row== 63.25) //ak grass
    	//if (telmnt[0].col == -83.75 && telmnt[0].row== 35.0)	//us forest
    	//if (telmnt[0].col == -99.5 && telmnt[0].row== 36.25)	//us grass
    	//if (telmnt[0].col == -104.75 && telmnt[0].row== 34 )	//us shrub
    	//if (telmnt[0].col == -102 && telmnt[0].row== 22)	//mex shrub
    	if (telmnt[0].col == -89.75 && telmnt[0].row== 16.5)	//mex forest
        {
        ifstream inrandom;
        inrandom.open("randomfile.txt", ios::in);
        if (!inrandom) {cout <<" can not find random file: "<<endl;exit(-1);}
        double randomvalue;
        outfile1.open("out_gpp.txt", ios::out);
        outfile2.open("out_npp.txt", ios::out);
        outfile3.open("out_ntcb.txt", ios::out);
        outfile4.open("out_convertc.txt", ios::out);
        outfile5.open("out_ffdc.txt", ios::out);
        outfile6.open("out_fflc.txt", ios::out);
        outfile7.open("out_vconvert.txt", ios::out);
        outfile8.open("out_abovelit.txt", ios::out);
        outfile9.open("out_belowlit.txt", ios::out);
        int countrun=0;
        for (int i=0;i<300;i++)
        {
        	countrun++;
        	inrandom>>randomvalue;
        	telmnt[0].tem.randomnum = randomvalue;
        	cout <<"count: "<<countrun<<endl;
            for( dyr = 1; dyr < RTIME; ++dyr )
            {

              // Run climate module or read climate data from file to update
              // climate for grid cell during year "dyr"

              updateTCLMGridCell( dyr );

              // Run land cover module or read in land cover data from file
              // to update land cover characteristics for grid cell during year "dyr"

              updateTLCLUCGridCell( dyr );

      	    //cout <<"model check: " <<" cmnt: " <<telmnt[0].tem.veg.cmnt<<telmnt[0].tem.veg.getC2N()<<" cfall: " <<telmnt[0].tem.veg.getCFALL(telmnt[0].tem.veg.cmnt)<<" cov: " <<telmnt[0].tem.veg.getCOV(telmnt[0].tem.veg.cmnt)<<endl;

              // Run TEM for grid cell during year "dyr"
              //if (dyr <=2) cout <<"col0: "<<telmnt[0].col<<" row0: " <<telmnt[0].row<<"grdcnt: " <<grdcnt<< endl;

              updateTTEMGridCell( dyr);
              //cout <<"count: "<<countrun<<" dyr: "<<dyr<<endl;

             }

        }
        outfile1.close();
        outfile2.close();
        outfile3.close();
        outfile4.close();
        outfile5.close();
        outfile6.close();
        outfile7.close();
        outfile8.close();
        outfile9.close();
      }
      */
      //else
      {
    	  telmnt[0].tem.randomnum = 0.0;
      for( dyr = 1; dyr < RTIME; ++dyr )
      {

        // Run climate module or read climate data from file to update
        // climate for grid cell during year "dyr"

        updateTCLMGridCell( dyr );

        // Run land cover module or read in land cover data from file
        // to update land cover characteristics for grid cell during year "dyr"

        updateTLCLUCGridCell( dyr );

	    //cout <<"model check: " <<" cmnt: " <<telmnt[0].tem.veg.cmnt<<telmnt[0].tem.veg.getC2N()<<" cfall: " <<telmnt[0].tem.veg.getCFALL(telmnt[0].tem.veg.cmnt)<<" cov: " <<telmnt[0].tem.veg.getCOV(telmnt[0].tem.veg.cmnt)<<endl;

        // Run TEM for grid cell during year "dyr"
        //if (dyr <=2) cout <<"col0: "<<telmnt[0].col<<" row0: " <<telmnt[0].row<<"grdcnt: " <<grdcnt<< endl;

        updateTTEMGridCell( dyr);

       }


      }
    }

   
    elmnt.show( nthproc, flog1,
                telmnt[0].col, 
                telmnt[0].row,
                telmnt[0].tem.totyr, 
                telmnt[0].tem.tol );

    ++grdcnt;
  }
  
  //
  if( 0 == fatalerr )
  {
    cout << "Extrapolation successfully completed - Congratulations! - on Process: " <<nthproc << endl;
    flog1 << "Extrapolation successfully completed - Congratulations! - on Process: " <<nthproc << endl;

    time_t etime;
    etime = time(0);
    cout << "Simulation ended @ "<< ctime( &etime );
    cout << "Total period (seconds): " << difftime(etime, stime) << endl;
    flog1 << "Total period (seconds): " << difftime(etime, stime) << endl;

  }
  else
  {
    if( elmnt.grdcnt != -99 && elmnt.count <= elmnt.grdcnt )
    {
      cout << "FATAL ERROR! Program Terminated - on Process: " <<nthproc << endl;
    }
    flog1 << "FATAL ERROR! Program Terminated - on Process: " <<nthproc << endl;
  }

  // Finished processing all elements - close open output files
 
  if( 1 == telmnt[0].clm.predflag )
  {
    for( i = 0; i < telmnt[0].natmspred; ++i ) 
    { 
      fclmpred[i].close(); 
    }
  }
  
  if( 1 == temflag )
  {

    for( i = 0; i < telmnt[0].ntempred; ++i )
    {
      ftempred[i].close();
    }

    if( ostateflag != 0 ) { ofstate.close(); }
  }

  flog1 << "Closed all files!" << endl << endl;
  flog1.close();

#ifdef PMODE
  // **************fmy: parallalization *************************************/
  MPI::Finalize();
  // **************fmy: parallalization *************************************/
#endif

  // close input files
  if( 0 == telmnt[0].lonlatflag ) { fclose( flonlat ); }
  if( 1 == telmnt[0].clm.cldflag ) { fclose( ifclds ); }
  else { fclose( ifnirr ); }

  if( 0 == telmnt[0].clm.sradflag ) { fclose( ifgirr ); }

  if( 1 == telmnt[0].clm.parflag ) { fclose( ifpar ); }

  fclose( iftair );
  fclose( ifprec );
  
  ifco2.close();
  
  fclose( ifo3 );
  fclose( ifndep );
  
  if( 1 == temflag )
  {
    fclose( ifnumchrts );
    fclose( iflulc ); 
    
    fclose( fstxt );
    fclose( felev );
    
    for( i = 0; i < telmnt[0].ntempred; ++i ) 
    
    if( istateflag != 0 ) { ifstate.close(); }
  }

  cout << "Closed all files!" << endl << endl;

  return 0;

};

/* *************************************************************
******************** End of Main ******************************* 
************************************************************* */

/* *************************************************************
************************************************************* */

int askpred( vector<string>& pvarname, 
             const int& posspred,
             vector<string>& predmap ) 
{
  const int MAXCOLUMNS = 7;
  int count = 0;

  int numpred;
  
  int i;
  int j;
  int cnt;
  int length;
  vector<string>::const_iterator dn;

#ifndef PMODE
  cout << endl << endl;
  cout << "           POSSIBLE OUTPUT VARIABLES:";
  cout << endl << endl;

  
  for( dn = pvarname.begin(); dn != pvarname.end(); ++dn )
  {
    cout << std::setw(10) << *dn << " ";
    
    ++count;
  
    if( 0 == count%MAXCOLUMNS ) { cout << endl; }
  }
  cout << endl << endl;
#endif

  flog1 << endl << endl;
  flog1 << "           POSSIBLE OUTPUT VARIABLES:";
  flog1 << endl << endl;

  count = 0;
  
  for( dn = pvarname.begin(); dn != pvarname.end(); ++dn )
  {
    flog1 << std::setw(10) << *dn << " ";
    
    ++count;
    
    if( 0 == count%MAXCOLUMNS ) { flog1 << endl; }
  }

  flog1 << endl << endl;

#ifdef PMODE
	fgo >> numpred;

#else

  cout << endl << endl;
  cout << "How many variables are to be mapped (max ";
  cout << posspred << ") in output files?  ";

  cin >> numpred;

  cout << numpred << endl;
#endif

  flog1 << endl << endl;
  flog1 << "How many variables are to be mapped (max ";
  flog1 << posspred << ") in output files?";
  flog1 << numpred << endl << endl;

#ifndef PMODE
  cout << "Please enter output variable: " << endl;
#endif
  flog1 << "Please enter output variable: " << endl;

  for( i = 0; i < numpred; ++i )
  {
    cnt = i + 1;

#ifdef PMODE
    fgo >> predmap.at( i );
#else
    cout << cnt << " ";
    cin >> predmap.at( i );
#endif

    length = predmap.at( i ).length();
    
    for( j = 0; j < length; ++j ) 
    { 
      predmap.at( i ).at( j ) = toupper( predmap.at( i ).at( j ) ); 
    }
    
    flog1 << cnt << " " << predmap.at( i ) << endl;
  }

  return numpred;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void initCLM( void )
{
  double initco2;
  double co2level;
  
  // Set flags and filenames for spatially explicit cloudiness 
  //   or solar radiation data
  
  telmnt[0].clm.setCldsFlags( flog1, equil );
  telmnt[0].clm.initSolarRad( flog1 ); 


  // Set flags and filenames for spatially explicit air 
  //   temperature data

  telmnt[0].clm.setTairFlags( flog1, equil );
  telmnt[0].clm.initTair( flog1 );


  // Set flags and filenames for spatially explicit 
  //   precipitation data

  telmnt[0].clm.setPrecFlags( flog1, equil );
  telmnt[0].clm.initPrec( flog1 );
  	

  // Determine initial atmospheric CO2 concentration for 
  //   initial equilibrium portion of the simulation
  
#ifdef PMODE
	fgo >> initco2;
	fgo >> co2level;

#else
  cout << endl << endl;
  cout << "Enter the initial concentration of carbon dioxide ";
  cout << "in ppmv: ";

  cin >> initco2;
  
  cout << "Enter the final equilibrium concentration of ";
  cout << "carbon dioxide in ppmv: ";

  cin >> co2level;
#endif

  telmnt[0].clm.setINITCO2( initco2 );

  flog1 << endl << endl;
  flog1 << "Enter the initial concentration of carbon dioxide ";
  flog1 << "in ppmv: ";
  flog1 << telmnt[0].clm.getINITCO2() << endl << endl;

  telmnt[0].clm.setCO2LEVEL( co2level );

  flog1 << "Enter the final equilibrium concentration of ";
  flog1 << "carbon dioxide in ppmv: ";
  flog1 << telmnt[0].clm.getCO2LEVEL() << endl << endl;
  
  
  // Set flags and filenames for global annual or monthly 
  //   spatially explicit atmospheric CO2 concentration 
  //   data for transient portion of the simulation
  
  telmnt[0].clm.setCO2Flags( flog1, equil );
  telmnt[0].clm.initCO2( flog1 );
 

  // Set flags and filenames for spatially explicit ozone data

  telmnt[0].clm.setO3Flags( flog1, equil );
  telmnt[0].clm.initO3( flog1 );

  // Set flags and filenames for spatially explicit nitrogen
  //  deposition data

  telmnt[0].clm.setNdepFlags( flog1, equil );
  telmnt[0].clm.initNdep( flog1 );

  // Identify output variables to be written out from climate 
  //   module
  
  setClmPred();

  // Initialize model for a grid cell
    
  initializeCLMGridCell();

  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */
/*
void initElmntBatch( void )
{
  int i;
  int j;

  for( i = 1; i < MAXGRID; i++ )
  {
    telmnt[i].clm.tcldsflag = telmnt[0].clm.tcldsflag;

    if( 1 == telmnt[0].clm.sradflag )
    {
      telmnt[i].lonlatflag = telmnt[0].lonlatflag;
    }
    if( 1 == temflag )
    {
      telmnt[i].lonlatflag = telmnt[0].lonlatflag;
      telmnt[i].tem.ag.tlulcflag = telmnt[0].tem.ag.tlulcflag;

// Soil texture specific TEM parameters

      telmnt[i].tem.soil.pctpora = telmnt[0].tem.soil.pctpora;
      telmnt[i].tem.soil.pctporb = telmnt[0].tem.soil.pctporb;
      telmnt[i].tem.soil.fldcapa = telmnt[0].tem.soil.fldcapa;
      telmnt[i].tem.soil.fldcapb = telmnt[0].tem.soil.fldcapb;
      telmnt[i].tem.soil.wiltpta = telmnt[0].tem.soil.wiltpta;
      telmnt[i].tem.soil.wiltptb = telmnt[0].tem.soil.wiltptb;

// Initialize CO2 for all grid cells

      telmnt[i].tem.atms.initco2 = telmnt[0].tem.atms.initco2;
      telmnt[i].tem.atms.co2level = telmnt[0].tem.atms.co2level;

      for( j = 0; j < NUMVEG; j++ )
      {
  	telmnt[i].tem.soil.rootza[j] = telmnt[0].tem.soil.rootza[j];
	telmnt[i].tem.soil.rootzb[j] = telmnt[0].tem.soil.rootzb[j];
	telmnt[i].tem.soil.rootzc[j] = telmnt[0].tem.soil.rootzc[j];
	telmnt[i].tem.soil.minrootz[j] = telmnt[0].tem.soil.minrootz[j];
	telmnt[i].tem.veg.minleaf[j] = telmnt[0].tem.veg.minleaf[j];
	telmnt[i].tem.veg.aleaf[j] = telmnt[0].tem.veg.aleaf[j];
	telmnt[i].tem.veg.bleaf[j] = telmnt[0].tem.veg.bleaf[j];
	telmnt[i].tem.veg.cleaf[j] = telmnt[0].tem.veg.cleaf[j];

//Site specific TEM parameters

	telmnt[i].tem.vegca[j] = telmnt[0].tem.vegca[j];
	telmnt[i].tem.vegcb[j] = telmnt[0].tem.vegcb[j];
	telmnt[i].tem.strna[j] = telmnt[0].tem.strna[j];
	telmnt[i].tem.strnb[j] = telmnt[0].tem.strnb[j];
	telmnt[i].tem.solca[j] = telmnt[0].tem.solca[j];
	telmnt[i].tem.solcb[j] = telmnt[0].tem.solcb[j];
	telmnt[i].tem.solna[j] = telmnt[0].tem.solna[j];
	telmnt[i].tem.solnb[j] = telmnt[0].tem.solnb[j];
	telmnt[i].tem.nh4a[j] = telmnt[0].tem.nh4a[j];
	telmnt[i].tem.nh4b[j] = telmnt[0].tem.nh4b[j];
	telmnt[i].tem.no3cut[j] = telmnt[0].tem.no3cut[j];
	telmnt[i].tem.no31a[j] = telmnt[0].tem.no31a[j];
	telmnt[i].tem.no31b[j] = telmnt[0].tem.no31b[j];
	telmnt[i].tem.no32a[j] = telmnt[0].tem.no32a[j];
	telmnt[i].tem.no32b[j] = telmnt[0].tem.no32b[j];
	telmnt[i].tem.stona[j] = telmnt[0].tem.stona[j];
	telmnt[i].tem.stonb[j] = telmnt[0].tem.stonb[j];
	telmnt[i].tem.veg.unleaf12[j] = telmnt[0].tem.veg.unleaf12[j];
	telmnt[i].tem.veg.initleafmx[j] = telmnt[0].tem.veg.initleafmx[j];
	telmnt[i].tem.veg.cmaxcut[j] = telmnt[0].tem.veg.cmaxcut[j];
	telmnt[i].tem.veg.cmax1a[j] =  telmnt[0].tem.veg.cmax1a[j];
	telmnt[i].tem.veg.cmax1b[j] = telmnt[0].tem.veg.cmax1b[j];
	telmnt[i].tem.veg.cmax2a[j] = telmnt[0].tem.veg.cmax2a[j];
	telmnt[i].tem.veg.cmax2b[j] = telmnt[0].tem.veg.cmax2b[j];
	telmnt[i].tem.veg.cfall[j] = telmnt[0].tem.veg.cfall[j];
	telmnt[i].tem.veg.kra[j] = telmnt[0].tem.veg.kra[j];
	telmnt[i].tem.veg.krb[j] = telmnt[0].tem.veg.krb[j];
	telmnt[i].tem.microbe.kdcut[j] = telmnt[0].tem.microbe.kdcut[j];
	telmnt[i].tem.microbe.kd1a[j] = telmnt[0].tem.microbe.kd1a[j];
	telmnt[i].tem.microbe.kd1b[j] = telmnt[0].tem.microbe.kd1b[j];
	telmnt[i].tem.microbe.kd2a[j] = telmnt[0].tem.microbe.kd2a[j];
	telmnt[i].tem.microbe.kd2b[j] = telmnt[0].tem.microbe.kd2b[j];
	
	telmnt[i].tem.microbe.lcclnc[j] = telmnt[0].tem.microbe.lcclnc[j];
	telmnt[i].tem.microbe.propftos[j] = telmnt[0].tem.microbe.propftos[j];

	telmnt[i].tem.veg.nfixpara[j] = telmnt[0].tem.veg.nfixpara[j];
	telmnt[i].tem.veg.nfixparb[j] = telmnt[0].tem.veg.nfixparb[j];

	telmnt[i].tem.veg.nupnh4cut[j] = telmnt[0].tem.veg.nupnh4cut[j];
	telmnt[i].tem.veg.nupnh41a[j] = telmnt[0].tem.veg.nupnh41a[j];
	telmnt[i].tem.veg.nupnh41b[j] = telmnt[0].tem.veg.nupnh41b[j];
	telmnt[i].tem.veg.nupnh42a[j] = telmnt[0].tem.veg.nupnh42a[j];
	telmnt[i].tem.veg.nupnh42b[j] = telmnt[0].tem.veg.nupnh42b[j];

	telmnt[i].tem.veg.nupno3cut[j] = telmnt[0].tem.veg.nupno3cut[j];
	telmnt[i].tem.veg.nupno31a[j] = telmnt[0].tem.veg.nupno31a[j];
	telmnt[i].tem.veg.nupno31b[j] = telmnt[0].tem.veg.nupno31b[j];
	telmnt[i].tem.veg.nupno32a[j] = telmnt[0].tem.veg.nupno32a[j];
	telmnt[i].tem.veg.nupno32b[j] = telmnt[0].tem.veg.nupno32b[j];

	telmnt[i].tem.veg.nfall[j] = telmnt[0].tem.veg.nfall[j];

	telmnt[i].tem.microbe.nh4immcut[j] = telmnt[0].tem.microbe.nh4immcut[j];
	telmnt[i].tem.microbe.nh4imm1a[j] = telmnt[0].tem.microbe.nh4imm1a[j];
	telmnt[i].tem.microbe.nh4imm1b[j] = telmnt[0].tem.microbe.nh4imm1b[j];
	telmnt[i].tem.microbe.nh4imm2a[j] = telmnt[0].tem.microbe.nh4imm2a[j];
	telmnt[i].tem.microbe.nh4imm2b[j] = telmnt[0].tem.microbe.nh4imm2b[j];

	telmnt[i].tem.soil.lchNO3parcut[j] = telmnt[0].tem.soil.lchNO3parcut[j];
	telmnt[i].tem.soil.lchNO3par1a[j] = telmnt[0].tem.soil.lchNO3par1a[j];
	telmnt[i].tem.soil.lchNO3par1b[j] = telmnt[0].tem.soil.lchNO3par1b[j];
	telmnt[i].tem.soil.lchNO3par2a[j] = telmnt[0].tem.soil.lchNO3par2a[j];
	telmnt[i].tem.soil.lchNO3par2b[j] = telmnt[0].tem.soil.lchNO3par2b[j];

	telmnt[i].tem.veg.initcneven[j] = telmnt[0].tem.veg.initcneven[j];
	telmnt[i].tem.veg.cnmin[j] = telmnt[0].tem.veg.cnmin[j];
	telmnt[i].tem.veg.c2na[j] = telmnt[0].tem.veg.c2na[j];
	telmnt[i].tem.veg.c2nb[j] = telmnt[0].tem.veg.c2nb[j];
	telmnt[i].tem.veg.c2nmin[j] = telmnt[0].tem.veg.c2nmin[j];
	telmnt[i].tem.microbe.cnsoil[j] = telmnt[0].tem.microbe.cnsoil[j];

        // Vegetation specific TEM parameters

	telmnt[i].tem.veg.kc[j] = telmnt[0].tem.veg.kc[j];
	telmnt[i].tem.veg.ki[j] = telmnt[0].tem.veg.ki[j];
	telmnt[i].tem.veg.gva[j] = telmnt[0].tem.veg.gva[j];
	telmnt[i].tem.veg.tmin[j] = telmnt[0].tem.veg.tmin[j];
	telmnt[i].tem.veg.toptmin[j] = telmnt[0].tem.veg.toptmin[j];
	telmnt[i].tem.veg.toptmax[j] = telmnt[0].tem.veg.toptmax[j];
	telmnt[i].tem.veg.tmax[j] = telmnt[0].tem.veg.tmax[j];
	telmnt[i].tem.veg.raq10a0[j] = telmnt[0].tem.veg.raq10a0[j];
	telmnt[i].tem.veg.raq10a1[j] = telmnt[0].tem.veg.raq10a1[j];
	telmnt[i].tem.veg.raq10a2[j] = telmnt[0].tem.veg.raq10a2[j];
	telmnt[i].tem.veg.raq10a3[j] = telmnt[0].tem.veg.raq10a3[j];
	telmnt[i].tem.veg.kn1[j] = telmnt[0].tem.veg.kn1[j];
	telmnt[i].tem.veg.labncon[j] = telmnt[0].tem.veg.labncon[j];

        telmnt[i].tem.veg.leafmxc[j] = telmnt[0].tem.veg.leafmxc[j];
	telmnt[i].tem.veg.kleafc[j] = telmnt[0].tem.veg.kleafc[j];
	telmnt[i].tem.veg.sla[j] = telmnt[0].tem.veg.sla[j];
	telmnt[i].tem.veg.cov[j] = telmnt[0].tem.veg.cov[j];
	telmnt[i].tem.veg.fpcmax[j] = telmnt[0].tem.veg.fpcmax[j];

        // Vegetation specific TEM parameters for microbes

	telmnt[i].tem.microbe.rhq10[j] = telmnt[0].tem.microbe.rhq10[j];
	telmnt[i].tem.microbe.kn2[j] = telmnt[0].tem.microbe.kn2[j];
	telmnt[i].tem.microbe.moistmin[j] = telmnt[0].tem.microbe.moistmin[j];
	telmnt[i].tem.microbe.moistopt[j] = telmnt[0].tem.microbe.moistopt[j];
	telmnt[i].tem.microbe.moistmax[j] = telmnt[0].tem.microbe.moistmax[j];

	telmnt[i].tem.ag.nvretconv[j] = telmnt[0].tem.ag.nvretconv[j];
	telmnt[i].tem.ag.nsretconv[j] = telmnt[0].tem.ag.nsretconv[j];
      }
    }
  }

};
*/

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void initializeCLMGridCell()
{   
  ostringstream tempfname;
  string girrname;
  string cldsname;
  string nirrname;
  string parname;
  string tairname;
  string precname;
  string co2name;
  string o3name;
  string ndepname; 

  
  // Open cloudiness or solar radiation file

  if( 0 == telmnt[0].clm.sradflag )
  {
    tempfname.str( "" );

    tempfname << telmnt[0].clm.igirrfname
              << telmnt[0].clm.igirrend;

    girrname = tempfname.str();

    ifgirr = fopen( girrname.c_str(), "r" );

    if( !ifgirr )
    {
      flog1 << endl << "Cannot open " << girrname;
      flog1 << " for GIRR data input" << endl << endl;

      exit( -1 );
    }
  }      

  if( 1 == telmnt[0].clm.cldflag )
  {
    tempfname.str( "" );

    tempfname << telmnt[0].clm.icldsfname
              << telmnt[0].clm.icldsend;

    cldsname = tempfname.str();

    ifclds = fopen( cldsname.c_str(), "r" );

    if( !ifclds )
    {
      flog1 << endl << "Cannot open " << cldsname;
      flog1 << " for CLDS data input" << endl << endl;

      exit( -1 );
    }
  }
  else
  {
    tempfname.str( "" );

    tempfname << telmnt[0].clm.inirrfname
              << telmnt[0].clm.inirrend;

    nirrname = tempfname.str();

    ifnirr = fopen( nirrname.c_str(), "r" );

    if( !ifnirr )
    {
      flog1 << endl << "Cannot open " << nirrname;
      flog1 << " for NIRR data input" << endl << endl;

      exit( -1 );
    }
  }

  if( 1 == telmnt[0].clm.parflag )
  {
    tempfname.str( "" );

    tempfname << telmnt[0].clm.iparfname
              << telmnt[0].clm.iparend;

    parname = tempfname.str();

    ifpar = fopen( parname.c_str(), "r" );

    if( !ifpar )
    {
      flog1 << endl << "Cannot open " << parname;
      flog1 << " for PAR data input" << endl << endl;

      exit( -1 );
    }
  }      


  // Open air temperature file

  tempfname.str( "" );

  tempfname << telmnt[0].clm.itairfname
            << telmnt[0].clm.itairend;

  tairname = tempfname.str();

  iftair = fopen( tairname.c_str(), "r" );

  if( !iftair )
  {
    flog1 << endl << "Cannot open " << tairname;
    flog1 << " for TAIR data input" << endl << endl;

    exit( -1 );
  }

  
  // Open precipitation file

  tempfname.str( "" );

  tempfname << telmnt[0].clm.iprecfname
            << telmnt[0].clm.iprecend;

  precname = tempfname.str();

  ifprec = fopen( precname.c_str(), "r" );

  if( !ifprec )
  {
    flog1 << endl << "Cannot open " << precname;
    flog1 << " for PREC data input" << endl << endl;

    exit( -1 );
  }

  if( 1 == telmnt[0].clm.tco2flag )
  {
    co2name = telmnt[0].clm.ico2fname;
    
    ifco2.open( co2name.c_str(), ios::in );

    if( !ifco2 )
    {
      flog1 << endl << "Cannot open " << co2name;
      flog1 << " for CO2 data input" << endl;

     exit( -1 );
    }
  }
  
  // Open ozone file

  tempfname.str( "" );

  tempfname << telmnt[0].clm.io3fname
            << telmnt[0].clm.io3end;

  o3name = tempfname.str();

  ifo3 = fopen( o3name.c_str(), "r" );
  if( !ifo3 )
  {
    flog1 << endl << "Cannot open " << o3name;
    flog1 << " for O3 data input" << endl << endl;

    exit( -1 );
  }

  // Open Atmospheric N deposition file

  tempfname.str( "" );

  tempfname << telmnt[0].clm.indepfname
            << telmnt[0].clm.independ;

  ndepname = tempfname.str();

  ifndep = fopen( ndepname.c_str(), "r" );

  if( !ifndep )
  {
    flog1 << endl << "Cannot open " << ndepname;
    flog1 << " for NDEP data input" << endl << endl;

    exit( -1 );
  } 
  
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void initializeLCLUCGridCell( void )
{
  string mxcohrtname;
  string lulcname;
  ostringstream tempfname;


  // Open maximum cohort file

  tempfname.str( "" );

  tempfname << telmnt[0].lcluc.imxcohrtfname
            << telmnt[0].lcluc.imxcohrtend;

  mxcohrtname = tempfname.str();

  ifnumchrts = fopen( mxcohrtname.c_str(), "r" );

  if( !ifnumchrts )
  {
    flog1 << endl << "Cannot open " << mxcohrtname;
    flog1 << " for MXCOHRTS data input" << endl << endl;

    exit( -1 );
  }

  // Open land use/land cover cohort file

  tempfname.str( "" );

  tempfname << telmnt[0].lcluc.ilulcfname
            << telmnt[0].lcluc.ilulcend;

  lulcname = tempfname.str();

  iflulc = fopen( lulcname.c_str(), "r" );

  if( !iflulc )
  {
    flog1 << endl << "Cannot open " << lulcname;
    flog1 << " for LULCCHRT data input" << endl << endl;

    exit( -1 );
  }
   
};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void initializeTEMGridCell( void )
{
  int dm;
  const int dyr = 0;
  int ichrt;

// *************************************************************
// Initialize TEM parameters for a "batch" of telmnts when
//   using a parallelized version of TEM

//  initElmntBatch();

// *************************************************************

// *************************************************************
// Skip to the desired record in the GIS data sets

//  grdcnt = skipRecords();

// *************************************************************

  // Set elevation and soil texture for TEM

  telmnt[0].setGIStopography( flog1, fatalerr, fstxt, felev );
  ifstream in_mask;
 	  in_mask.open("newzone_mask.txt",ios::in);
 	  int **maskdata;
 	  maskdata = new int*[296];
 	  for(int i=0;i<296;i++){
 	          maskdata[i] = new int[480];
 	  }
 	  for (int i=0;i<296;i++)
 	  {
 	  for (int j=0;j<480;j++)
 	  {
 	          in_mask>>maskdata[i][j];

 	  }
 	  }
 	  for(int i=0;i<296;i++) delete[] maskdata[i];
 	  delete[] maskdata;
     
/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */
 int xx=telmnt[0].maxcohorts;
 if (telmnt[0].maxcohorts==0) xx=1;
 for( ichrt = 0; ichrt < xx; ++ichrt )
  {

	//if (telmnt[0].col > -109 && telmnt[0].col <-76 && telmnt[0].row == 54) cout <<" ichrt1: " <<ichrt<<" maxcohort: " << telmnt[0].maxcohorts <<" col: " <<telmnt[0].col<< " row: "<<telmnt[0].row<<endl;
	if( istateflag > 0 )
    {
      // Read in initial TEM state determined in a previous 
      //  TEM simulation to telmnt[0].cohort
	  //cout <<"maxcohort: "<<telmnt[0].maxcohorts<<" ichrt: "<<ichrt<<endl;
	  //int rrr = (-1.0 * telmnt[0].row + 83.75)*4.0;
	  //int lll = (170.0 + telmnt[0].col) * 4.0;
      //cout <<"maskdata: "<<maskdata[rrr][lll]<<endl;
      //if (maskdata[rrr][lll] > 0)  telmnt[0].readCohortState( ifstate, ichrt );
      telmnt[0].readCohortState( ifstate, ichrt );
      //telmnt[0].readCohortState0( ifstate, ichrt );
      //if (telmnt[0].col == -123.0 && telmnt[0].row== 44.25 && telmnt[0].tem.veg.cmnt==7) cout <<" vegc: " <<telmnt[0].cohort[ichrt].y[2]<<endl;
      //if (telmnt[0].dumflt1 != telmnt[0].col && telmnt[0].dumflt2 != telmnt[0].row)
      //{
    	 // cout <<"intial data is not the same with climate data at: "<<telmnt[0].col<<", "<<telmnt[0].row <<" dump1: " <<telmnt[0].dumflt1 <<" dump2: "<<telmnt[0].dumflt2<<endl;
         // exit(-1);
       //}
      // Pass telmnt[0].cohort information to TEM
      //double xx = telmnt[0].cohort[ichrt].y[telmnt[0].I_SOLC];
      //cout <<"initial_soilc:" <<xx<<" col: " <<telmnt[0].col<< "row: "<<telmnt[0].row<<" ichrt: " <<ichrt<<endl;
      dm = CYCLE - 1;
      
      telmnt[0].getTEMCohortState( ichrt, dm );
      //cout <<" I_soilc1:" <<xx<<" col: " <<telmnt[0].col<< "row: "<<telmnt[0].row<<endl;

      // Determine soil properties of element based on
      //   soil texture


      telmnt[0].tem.soil.xtext( telmnt[0].tem.veg.cmnt,
                                telmnt[0].tem.soil.getPCTSILT(),
                                telmnt[0].tem.soil.getPCTCLAY() );


      // Set texture dependent parameters for TEM 

      telmnt[0].tem.setELMNTecd( telmnt[0].tem.veg.cmnt, 
                                 telmnt[0].tem.soil.getPSIPLUSC() );
    }
    else
    {     	

     telmnt[0].setTEMequilState( flog1,
                                  equil,
                                  totsptime,
                                  ichrt );


      if( telmnt[0].tem.intflag > 0 )
      {
        if( elmnt.count < elmnt.grdcnt )
        {
          cout << "Integration terminated before attaining ";
          cout << "tolerance level" << endl;
        }

        flog1 << "Integration terminated before attaining ";
        flog1 << "tolerance level" << endl;

        telmnt[0].tem.intflag = 0;
      }
      
      // Write out telmnt[0].cohort to output file for
      //   potential use in a future TEM simulation
      
      if( 1 == ostateflag )
      {
        telmnt[0].writeCohortState( ofstate, ichrt );
      }
      // Write selected TEM variables from telmnt[0].output to 
      //   outfile files
  
      if( 1 == spinoutfg || 3 == spinoutfg || 1 == equil )
      {
        telmnt[0].temwritepred( ftempred, 
                                tempredmap, 
                                dyr, 
                                ichrt,
                                telmnt[0].ntempred );
      }

    } // End of istateflag else
  } // End of cohort loop

//  elmnt.show( flog1,
//              telmnt[0].col,
//              telmnt[0].row,
//              telmnt[0].tem.totyr,
//              telmnt[0].tem.tol );
 
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void initLCLUC( void )
{
  // Hand 0ff startyr from clm
  
  telmnt[0].lcluc.startyr = lulcstartyr;


  // Get community type parameterizations
  
  telmnt[0].lcluc.getvtype( flog1 );
 
  
  // Set flags land cover /land-use change data

  telmnt[0].lcluc.setLCLUCFlags( flog1, equil );


  // Set filename for spatially explicit number of cohorts data
 
  telmnt[0].lcluc.initMaxCohorts( flog1 );

  
  // Set filename for spatially explicit land use/ land cover 
  //   cohort data
  
  telmnt[0].lcluc.initCohorts( flog1 );


  initializeLCLUCGridCell();
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void initTEM( void )
{

  int numcmnt;
  string ifilename;

  double dc2n;

#ifdef PMODE
    fgo >> temflag;
#else

  cout << "Do you want to run the terrestrial ecosystem model ";
  cout << "(TEM)?" << endl;
  cout << "Enter 0 for No" << endl;
  cout << "Enter 1 for Yes: " << endl;
  
  cin >> temflag;
#endif

  flog1 << "Do you want to run the terrestrial ecosystem model ";
  flog1 << "(TEM)?" << endl;
  flog1 << "  Enter 0 for No" << endl;
  flog1 << "  Enter 1 for Yes" << endl;
  flog1 << "temflag = " << temflag << endl << endl;

  telmnt[0].ntempred = 0;

  if( 1 == temflag )
  {
    // Hand off model startyr to TEM
    
    telmnt[0].tem.modstartyr = telmnt[0].clm.modstartyr;

    // Hand off initial CO2 data to TEM
    
    telmnt[0].tem.atms.setINITCO2( telmnt[0].clm.getINITCO2() );

    telmnt[0].tem.atms.setCO2LEVEL( telmnt[0].clm.getCO2LEVEL() );


    // Hand off tlulcflag from LCLUC module to TEM
    telmnt[0].tem.ag.tlulcflag = telmnt[0].lcluc.tlulcflag;

#ifdef PMODE
    fgo >> telmnt[0].tem.soil.stmflg;
#else

    cout << endl;
    cout << "Do you want to run the SOIL THERMAL MODEL (STM) ";
    cout << "model for soil temperatures?" << endl;
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes" << endl;

    cin >> telmnt[0].tem.soil.stmflg;
#endif

    flog1 << endl;
    flog1 << "Do you want to run the SOIL THERMAL MODEL (STM) ";
    flog1 << "model for soil temperatures?"<< endl;
    flog1 << "  Enter 0 for No" << endl;
    flog1 << "  Enter 1 for Yes" << endl;
    flog1 << " telmnt[0].tem.soil.stmflg = ";
    flog1 << telmnt[0].tem.soil.stmflg << endl << endl;


    telmnt[0].tem.initrun( flog1 );
    telmnt[0].tem.askODE( flog1 );


    // Get soil texture dependent parameters

    telmnt[0].tem.soil.getecd( flog1 );

    // Get vegetation type dependent parameters

    telmnt[0].tem.soil.getrootz( flog1 );
    telmnt[0].tem.veg.getecd( flog1 );	
    telmnt[0].tem.veg.getleafecd( flog1 );
    telmnt[0].tem.microbe.getvegecd( flog1 );


    // Get parameters associated with human disturbance 
    //   activities

    //below line is closed by cgs2014. ag exists when run at equilibrium status or non-land use change experiments
    //if( 1 == telmnt[0].tem.ag.tlulcflag )
    {
      telmnt[0].tem.ag.getecd( flog1 );
      telmnt[0].tem.ag.setAgricFlags( flog1 );
    }

    if( 1 == telmnt[0].tem.soil.stmflg )
    {
      // Get snow and soil layer dependent parameters for STM
      telmnt[0].tem.soil.stm.getsnowecd( flog1 );   
      telmnt[0].tem.soil.stm.getsoillecd( flog1 );  
//      telmnt[0].tem.soil.stm.getsoiltecd( flog1 );
    }


   // Get calibration site specific parameters
#ifdef PMODE
    fgo >> numcmnt;
#else
    cout << "Please enter the number of community types with ";
    cout << "calibration data for vegetation:";
    
    cin >> numcmnt;
#endif

    flog1 << endl << endl;
    flog1 << "Please enter the number of community types with ";
    flog1 << "calibration data for vegetation:";
    flog1 << numcmnt << endl;

    telmnt[0].tem.getsitecd( numcmnt, flog1 );

#ifdef PMODE
    fgo >> ifilename;
#else
    cout << "Please enter the name of the file containing the ";
    cout << "soil texture data:";
    cout << endl;
    cout << "               (e.g., TEXTURE.GIS) " << endl;

    cin >> ifilename;
#endif

    flog1 << "Please enter the name of the file containing the ";
    flog1 << "soil texture data:";
    flog1 << endl;
    flog1 << "               (e.g., TEXTURE.GIS) " << endl;
    flog1 << ifilename << endl << endl;
    
    fstxt = fopen( ifilename.c_str(), "r" );

    if( !fstxt )
    {
      flog1 << endl << "Cannot open " << ifilename;
      flog1 << " for data input" << endl;
      
      exit( -1 );
    }

#ifdef PMODE
    fgo >> ifilename;
#else
    cout << "Please enter the name of the file containing the ";
    cout << "elevation data: " << endl;
    cout << "               (e.g., ELEV.GIS) " << endl;

    cin >> ifilename;
#endif
    
    flog1 << "Please enter the name of the file containing the ";
    flog1 << "elevation data: " << endl;
    flog1 << "               (e.g., ELEV.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    felev = fopen( ifilename.c_str(), "r" );

    if( !felev )
    {
      flog1 << "\nCannot open " << ifilename;
      flog1 << " for data input" << endl;
      
      exit( -1 );
    }

#ifdef PMODE
    fgo >> dc2n;
#else
    cout << endl << endl;
    cout << "Enter the factor for changing C:N per ppmv of ";
    cout << "enhanced CO2:" << endl;
    cout << "           (Enter 0.0 for no change): " << endl;

    cin >> dc2n;
#endif

    telmnt[0].tem.veg.setDC2N( dc2n );

    flog1 << endl;
    flog1 << "Enter the factor for changing C:N per ppmv of ";
    flog1 << "enhanced CO2:" << endl;
    flog1 << "           (Enter 0.0 for no change): " << endl;
    flog1 << "telmnt[0].tem.veg.dc2n = ";
    flog1 << telmnt[0].tem.veg.getDC2N() << endl << endl;

    setTEMPred();
 
    // Identify file that has extant TEM State from previous 
    //   simulation

    setInTEMState();

    // Identify files to receive TEM States from particular years 
    //   in current simulation
  
    setOutTEMState(nthproc);
  }
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void setClmPred( void )
{

  int i;
  string clmpredfile;

  telmnt[0].totpred = 0;
  telmnt[0].natmspred = 0;

  telmnt[0].clm.predflag = 0;
#ifdef PMODE
    fgo >> telmnt[0].clm.predflag;
#else

  cout << endl;
  cout << "Do you wish output data from the irradiance model? ";
  cout << endl;
  cout << "  Enter 0 for no" << endl;
  cout << "  Enter 1 for yes: ";

  cin >> telmnt[0].clm.predflag;
#endif

  flog1 << endl;
  flog1 << "Do you wish output data from the irradiance model?";
  flog1 << endl;
  flog1 << "  Enter 0 for no" << endl;
  flog1 << "  Enter 1 for yes: " << endl;
  flog1 << "telmnt[0].clm.predflag = ";
  flog1 << telmnt[0].clm.predflag << endl << endl;

  if( 1 == telmnt[0].clm.predflag )
  {
    telmnt[0].natmspred = askpred( telmnt[0].clm.predstr,
                                   NUMATMS, 
                                   clmpredmap );
                                   
    telmnt[0].totpred = telmnt[0].natmspred;

    for( i = 0; i < telmnt[0].natmspred; ++i )
    {

#ifdef PMODE
      fgo >> clmpredfile;

    // **************fmy: parallalization *************************************/
    // this will output data for simulations on each individual process,
    // with process rank number appeneding to the read-in filename, so needs gluing of all data after runs
    // (fmy: not yet figured out how to output into one single file from all processes! - todo!)
      ostringstream ss;
      ss<<nthproc;
      clmpredfile = clmpredfile+ss.str();
      // **************fmy: parallalization *************************************/
#else
      cout << endl;
      cout << "Enter the name of the OUTPUT file to contain ";
      cout << clmpredmap[i] << ":  ";

      cin >> clmpredfile;
#endif

      flog1 << endl;
      flog1 << "Enter the name of the OUTPUT file to contain ";
      flog1 << clmpredmap[i] << ":  " << clmpredfile << endl;

      fclmpred[i].open( clmpredfile.c_str(), ios::out );
    }
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void setGridParameters( void )
{

  string ifilename;
  
#ifdef PMODE
	fgo >> mxnumgrid;
	fgo >> telmnt[0].lonlatflag;

#else
  
  cout << "How many elements are in the gridded data sets?";
  cin >> mxnumgrid;

  cout << "How do you locate your grid cells?" << endl;
  cout << "Enter 0 for column/row:" << endl;
  cout << "Enter 1 for longitude/latitude: ";
  cin >> telmnt[0].lonlatflag;
#endif

  flog1 << endl << endl;
  flog1 << "How many elements are in the gridded data sets? ";
  flog1 << mxnumgrid << endl << endl;

  flog1 << "How do you locate your grid cells?" << endl;
  flog1 << "Enter 0 for column/row:" << endl;
  flog1 << "Enter 1 for longitude/latitude: ";
  flog1 << telmnt[0].lonlatflag << endl << endl;

  if( 0 == telmnt[0].lonlatflag )
  {

#ifdef PMODE
	fgo >> ifilename;

#else
    cout << "Please enter the name of the file containing ";
    cout << "the latitude data: " << endl;
    cout << "               (e.g., LAT.GIS) " << endl;

    cin >> ifilename;
#endif

    flog1 << "Please enter the name of the file containing ";
    flog1 << "the latitude data: " << endl;
    flog1 << "               (e.g., LAT.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    flonlat = fopen( ifilename.c_str(), "r" );

    if( !flonlat )
    {
      cerr << endl << "Cannot open " << ifilename;
      cerr << " for data input" << endl;
      
      exit( -1 );
    }
  }

 // Use all records in GIS data sets (i.e. desired coverage)?

  elmnt.ask( flog1 );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void setInTEMState( void )
{

  string ifilename;
  
	
  istateflag = 0;
  istateyear = -99;
  
#ifdef PMODE
  fgo >> istateflag;
#else
  cout << endl;
  cout << "Do you want to use spatially explicit data ";
  cout << "for initial conditions? " << endl;
  cout << "  Enter 0 for NO:" << endl;
  cout << "  Enter 1 if spatially explicit data represents the ";
  cout << "state of TEM at equilibrium conditions:" << endl;
  cout << "  Enter 2 if spatially explicit data represents the ";
  cout << "state of TEM at the end of a specific year:" << endl;

  cin >> istateflag;
#endif

  flog1 << endl;
  flog1 << "Do you want to use spatially explicit data ";
  flog1 << "for intial conditions? " << endl;
  flog1 << "  Enter 0 for NO:" << endl;
  flog1 << "  Enter 1 if spatially explicit data represents the ";
  flog1 << "state of TEM at equilibrium conditions:" << endl;
  flog1 << "  Enter 2 if spatially explicit data represents the ";
  flog1 << "state of TEM at the end of a specific year:" << endl;
  flog1 << "istateflag = " << istateflag << endl << endl;

  if( istateflag > 0 )
  {

#ifdef PMODE
    fgo >> ifilename;
#else
    cout << "Please enter the name of the file containing the ";
    cout << "initial TEM state data: " << endl;
    cout << "               (e.g., TEMINIT.GIS) " << endl;

    cin >> ifilename;
#endif

    flog1 << "Please enter the name of the file containing the ";
    flog1 << "initial TEM state data: " << endl;
    flog1 << "               (e.g., TEMINIT.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    ifstate.open( ifilename.c_str(), ios::in );

    if( !ifstate )
    {
      cout << endl << "Cannot open " << ifilename;
      cout << " for data input" << endl;
      flog1 << endl << "Cannot open " << ifilename;
      flog1 << " for data input" << endl;
      
      exit( -1 );
    }

    if( 2 == istateflag )
    {

#ifdef PMODE
      fgo >> istateyear;
#else
      cout << "Please enter the year that you wish to use as ";
      cout << "the initial TEM state: " << endl;

      cin >> istateyear;
#endif
      flog1 << "Please enter the year that you wish to use as ";
      flog1 << "the initial TEM state: " << endl;
      flog1 << istateyear << endl << endl;
    }
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void setOutTEMState( const int &nthproc )
{

  string ifilename;
  
  
  ostateflag = 0;
  ostateyear = -99;

#ifdef PMODE
  fgo >> ostateflag;
#else
  cout << endl;
  cout << "Do you want to save the state of TEM as a spatially ";
  cout << "explicit data set for a specified year? " << endl;
  cout << "  Enter 0 for NO:" << endl;
  cout << "  Enter 1 for saving the state of TEM at equilibrium ";
  cout << "conditions:" << endl;
  cout << "  Enter 2 for saving the state of TEM at the end of ";
  cout << "a specific year:" << endl;

  cin >> ostateflag;
#endif

  flog1 << endl;
  flog1 << "Do you want to save the state of TEM as a spatially ";
  flog1 << "explicit data set for a specified year? " << endl;
  flog1 << "  Enter 0 for NO:" << endl;
  flog1 << "  Enter 1 for saving the state of TEM at equilibrium ";
  flog1 << "conditions:" << endl;
  flog1 << "  Enter 2 for saving the state of TEM at the end of ";
  flog1 << "a specific year:" << endl;
  flog1 << ostateflag << endl << endl;

  if( ostateflag != 0 )
  {

#ifdef PMODE
    fgo >> ifilename;
    // **************fmy: parallalization *************************************/
    // this will output data for simulations on each individual process,
    // with process rank number appeneding to the input filename, so needs gluing of all data
    // (fmy: not yet figured out how to output into one single file from all processes! - todo!)
    ostringstream ssnthproc;
    ssnthproc<<nthproc;
    ifilename = ifilename+ssnthproc.str();
    // **************fmy: parallalization *************************************/

#else

	cout << "Please enter the name of the file to contain ";
    cout << "the 'state' of TEM: " << endl;
    cout << "               (e.g., TEMSTATE.GIS) " << endl;

    cin >> ifilename;
#endif


    flog1 << "Please enter the name of the file to contain ";
    flog1 << "the 'state' of TEM: " << endl;
    flog1 << "               (e.g., TEMSTATE.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    ofstate.open( ifilename.c_str(), ios::out );

    if( 2 == ostateflag )
    {

#ifdef PMODE
      fgo >> ostateyear;
#else
      cout << "Please enter the year that you wish to save ";
      cout << "the TEM state: " << endl;

      cin >> ostateyear;
#endif

      flog1 << "Please enter the year that you wish to save ";
      flog1 << "the TEM state: " << endl;
      flog1 << ostateyear << endl << endl;
    }
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void setRunMode( const int& nthproc )
{

  // Assign a log file for the simulation
  string logfile="tem60.log";
#ifdef PMODE
  ostringstream ss;
  ss<<nthproc;
  logfile="tem60.log"+ss.str();
#endif
  flog1.open( logfile.c_str(), ios::out );

/* *************************************************************
  Run equilibrium simulation or transient simulation ?
************************************************************* */

  equil = 1;

#ifdef PMODE
  fgo >> equil;
#else
  cout << endl;
  cout << "Do you want to run the model only for steady state ";
  cout << "conditions ? " << endl;
  cout << " Enter 0 for transient simulation" << endl;
  cout << " Enter 1 for steady state simulation" << endl;
  cin >> equil;

#endif

  flog1 << endl;
  flog1 << "Do you want to run the model only for steady state ";
  flog1 << "conditions ? " << endl;
  flog1 << " Enter 0 for transient simulation" << endl;
  flog1 << " Enter 1 for steady state simulation" << endl;
  flog1 << "equil = " << equil << endl << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void setRunTime( void )
{
	
  RTIME = 1;
  
  if( 0 == equil )
  {

#ifdef PMODE
	fgo >> telmnt[0].clm.modstartyr;
	fgo >> spinflag;

#else
    cout << "What is the first year of the model simulation period?";
    cout << endl;

    cin >> telmnt[0].clm.modstartyr;

    // Start transient TEM run from equilibrium conditions 
    // (i.e. spinflag == 0) or with a "spin up" period to remove 
    // potential artifacts associated with changing from 
    // equilibrium conditions to transient conditions from model
    //  results

    cout << "Do you want to start the transient with a spin up period? ";
    cout << endl;
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes:" << endl;
 
    cin  >> spinflag;
#endif

    flog1 << "What is the first year of the model simulation period?";
    flog1 << endl;
    flog1 << telmnt[0].clm.modstartyr << endl << endl;

    telmnt[0].tem.modstartyr = telmnt[0].clm.modstartyr;

    flog1 << telmnt[0].clm.modstartyr << endl << endl;

    flog1 << "Do you want to start the transient with a spin up period? ";
    flog1 << endl;
    flog1 << "Enter 0 for no:" << endl;
    flog1 << "Enter 1 for yes: " << endl;
    flog1 << "spinflag = " << spinflag << endl << endl;

    totsptime = 0;
    
    if( 1 == spinflag )
    {
      // Specify conditions for initializing TEM with a transient 
      //   "spin up" period for a grid cell
#ifdef PMODE
      fgo >> numspin;
      fgo >> spintime;

#else

      cout << "How many spins do you want in the spin up period? ";
      cin >> numspin;
      
      cout << "How many years per spin? ";
      cin >> spintime;
#endif
      
      flog1 << "How many spins do you want in the spin up period? ";
      flog1 << endl;
      flog1 << "numspin = " << numspin << endl << endl;

      flog1 << "How many years per spin? " << endl;
      flog1 << "spintime = " << spintime << endl << endl;

      totsptime = spintime * numspin;
      flog1 << "totsptime = " << totsptime << endl << endl;
      RTIME += totsptime;
    }
    
    // Specify conditions for the "non-spin up" part of the 
    //   transient TEM run

#ifdef PMODE
	fgo >> transtime;

#else
    cout << endl;
    cout << "How many years do you run for transient simulations ? ";
    cout << endl;
    
    cin >> transtime;
#endif

    flog1 << endl;
    flog1 << "How many years do you run for transient simulations ? ";
    flog1 << endl;
    flog1 << "transtime = " << transtime << endl << endl;

    RTIME += transtime;
    flog1 << "RTIME = " << RTIME << endl << endl;
  }
  else 
  { 
    totsptime = RTIME; 
    spintime = 1;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void setTEMPred( void )
{

  int i;
  string tempredfile;

  telmnt[0].tem.predflag = 0;
  
#ifdef PMODE
  fgo >> telmnt[0].tem.predflag;
#else
  cout << endl;
  cout << "Do you wish spatially explicit output data from TEM?";
  cout << endl;
  cout << "  Enter 0 for no" << endl;
  cout << "  Enter 1 for yes: ";

  cin >> telmnt[0].tem.predflag;
#endif

  flog1 << endl;
  flog1 << "Do you wish spatially explicit output data from TEM?";
  flog1 << endl;
  flog1 << "  Enter 0 for no" << endl;
  flog1 << "  Enter 1 for yes: " << endl;
  flog1 << "telmnt[0].tem.predflag = " << telmnt[0].tem.predflag;
  flog1 << endl << endl;

  if( 1 == telmnt[0].tem.predflag )
  {
    telmnt[0].ntempred = askpred( telmnt[0].tem.predstr, 
                                  NUMTEM,
                                  tempredmap );
                                  
//    telmnt[0].totpred += telmnt[0].ntempred;

    for( i = 0; i < telmnt[0].ntempred; ++i )
    {

#ifdef PMODE
      fgo >> tempredfile;
      // **************fmy: parallalization *************************************/
      // this will output data for simulations on each individual process,
      // with process rank number appeneding to the read-in filename, so needs gluing of all data after runs
      // (fmy: not yet figured out how to output into one single file from all processes! - todo!)
      ostringstream ss;
      ss<<nthproc;
      tempredfile = tempredfile+ss.str();
      // **************fmy: parallalization *************************************/
#else
      cout << endl;
      cout << "Enter the name of the OUTPUT file to contain ";
      cout << tempredmap[i] << ":  ";

      cin >> tempredfile;
#endif

      
      flog1 << "Enter the name of the OUTPUT file to contain ";
      flog1 << tempredmap[i] << ":  " << tempredfile << endl;
      
      ftempred[i].open( tempredfile.c_str(), ios::out );
    }

    spinoutfg = 0;
    spinoutyrs = 0;
    
#ifdef PMODE
    fgo >> spinoutfg;
#else
    cout << "Do you want to save TEM output from ";
    cout << "the spin-up period?" << endl;
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for all spin-up years: " << endl;
    cout << "Enter 2 for some spinup years: " << endl;
    cout << "Enter 3 to save only equilibrium results:" << endl;
      
    cin >> spinoutfg;
#endif

    flog1 << "Do you want to save TEM output from ";
    flog1 << "the spin-up period?" << endl;
    flog1 << "Enter 0 for no:" << endl;
    flog1 << "Enter 1 for all spin-up years: " << endl;
    flog1 << "Enter 2 for some spinup years: " << endl;
    flog1 << "Enter 3 to save only equilibrium results:" << endl;
    flog1 << "spinoutfg = " << spinoutfg << endl << endl;

    if( 2 == spinoutfg )
    {

#ifdef PMODE
      fgo >> spinoutyrs;
#else
      cout << "How many years of spin-up do you want to ";
      cout << "save in the TEM output?" << endl;
        
      cin >> spinoutyrs;
#endif
      flog1 << "How many years of spin-up do you want to ";
      flog1 << "save in the TEM output?" << endl;     
      flog1 << "spinoutyrs = " << spinoutyrs << endl << endl;
    }
  }
  
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

long skipRecords( void )
{

  long elmntcnt = 0;
  string contnent;

  end1 = 1;
  fatalerr = 0;

  if( 0 == elmnt.strtflag )
  {

    for( elmntcnt = 0; elmntcnt < elmnt.numskip; ++elmntcnt )
    {       
      if( 1 == temflag )
      {
        telmnt[0].setGIStopography( flog1, 
                                    fatalerr, 
                                    fstxt, 
                                    felev );
      }      
    }
  }
  
  return elmntcnt;  // Number of grid cells skipped!

};

/* F.-M.Yuan: skip 'grdcnt' cells from the input 'txt' data files
 * NOTE: the txt files must be in lines, each of line is exactly for ONE grid-cell data
 *  (because here we skip lines in those files)
 *  IF this not works for whatever reason, comment out those 'fgets' lines and switch to '*.getdel(*FILE)' lines
 */
void skipGridCells(const int &grdcnt) {

	//int gisend = 1;
	char dummyline[1024];
    for(long elmntcnt = 0; elmntcnt < grdcnt; ++elmntcnt ) {

    	// latlon data file
    	if( 0 == telmnt[0].lonlatflag )
    	  //gisend = latdat.getdel( flonlat );
    		fgets(dummyline, sizeof(dummyline), flonlat );

		// climate data files
	    if( 1 == telmnt[0].clm.tcldsflag ) {

	      if( 0 == telmnt[0].clm.sradflag ) {
	        for( int dyr = 0; dyr < (transtime+1); ++dyr ){
	        	//gisend = girrdat[dyr].getdel( ifgirr );
		    	  fgets(dummyline, sizeof(dummyline), ifgirr);
	        }
	      }

	      if( 1 == telmnt[0].clm.cldflag ) {
	        for( int dyr = 0; dyr < (transtime+1); ++dyr ) {
	        	//gisend = cldsdat[dyr].getdel( ifclds );
		    	  fgets(dummyline, sizeof(dummyline), ifclds );
	        }
	      } else {
	        for( int dyr = 0; dyr < (transtime+1); ++dyr ){
	    	    //gisend = nirrdat[dyr].getdel( ifnirr );
		    	  fgets(dummyline, sizeof(dummyline), ifnirr );
	        }
	      }

	      if( 1 == telmnt[0].clm.parflag ){
	        for( int dyr = 0; dyr < (transtime+1); ++dyr ){
	        	//gisend = pardat[dyr].getdel( ifpar );
		    	  fgets(dummyline, sizeof(dummyline), ifpar );
        	}
	      }

		} else {  // 0 == telmnt[0].clm.tcldsflag

			if( 0 == telmnt[0].clm.sradflag )
				//gisend = girrdat[0].getdel( ifgirr );
	    	  fgets(dummyline, sizeof(dummyline), ifgirr );


			if( 1 == telmnt[0].clm.cldflag ){
				//gisend = cldsdat[0].getdel( ifclds );
		    	  fgets(dummyline, sizeof(dummyline), ifclds );

			}else{
				//gisend = nirrdat[0].getdel( ifnirr );
		    	  fgets(dummyline, sizeof(dummyline), ifnirr );

			}

			if( 1 == telmnt[0].clm.parflag )
				//gisend = pardat[0].getdel( ifpar );
	    	  fgets(dummyline, sizeof(dummyline), ifpar );

		}

	    // air temperature data file
		if( 1 == telmnt[0].clm.ttairflag ) {
	    	for( int dyr = 0; dyr < (transtime+1); ++dyr )
	    		//gisend = tairdat[dyr].getdel( iftair );
		    	  fgets(dummyline, sizeof(dummyline), iftair );

	    } else {
	    	//gisend = tairdat[0].getdel( iftair );
	    	  fgets(dummyline, sizeof(dummyline), iftair );
	    }

	    // precipitation data file
		if( 1 == telmnt[0].clm.tprecflag ) {
	    	for( int dyr = 0; dyr < (transtime+1); ++dyr )
	    	  //gisend = precdat[dyr].getdel( ifprec );
	    	  fgets(dummyline, sizeof(dummyline), ifprec );
	    } else {
	    	//gisend = precdat[0].getdel( ifprec );
	    	  fgets(dummyline, sizeof(dummyline), ifprec );
	    }

		// O3 data file
	    if( 1 == telmnt[0].clm.to3flag ) {
	    	for( int dyr = 0; dyr < (transtime+1); ++dyr )
	    	  //gisend = o3dat[dyr].getdel( ifo3 );
	    	  fgets(dummyline, sizeof(dummyline), ifo3 );

	    } else {
	    	//gisend = o3dat[0].getdel( ifo3 );
	    	  fgets(dummyline, sizeof(dummyline), ifo3 );
	    }

		// ndep data file
	    if( 1 == telmnt[0].clm.tndepflag ) {
	    	for( int dyr = 0; dyr < (transtime+1); ++dyr )
	    		//gisend = ndepdat[dyr].getdel( ifndep );
	    	  fgets(dummyline, sizeof(dummyline), ifndep );

	    } else {
			//gisend = ndepdat[0].getdel( ifndep );
	    	  fgets(dummyline, sizeof(dummyline), ifndep );
	    }

	    // lulc data file
		if( 1 == telmnt[0].lcluc.tlulcflag ){
	      for( int dyr = 0; dyr < (transtime+1); ++dyr ) {
	    	  mxcohrtdat[dyr].getdel( ifnumchrts );   // has to do this way, because need mxcohrtdat[yr].total
		      for( int ichrt = 0; ichrt < mxcohrtdat[dyr].total; ++ichrt ) {
		    	  fgets(dummyline, sizeof(dummyline), iflulc );

		      }
	      }

	    } else {
	    	  mxcohrtdat[0].getdel( ifnumchrts );   // has to do this way, because need mxcohrtdat[yr].total
		      for( int ichrt = 0; ichrt < mxcohrtdat[0].total; ++ichrt ) {
		    	  fgets(dummyline, sizeof(dummyline), iflulc );

		      }
	    }
		if( istateflag > 0 && 0 == elmnt.strtflag) { //added by cgs2014
			int xx=mxcohrtdat[0].total;
			if (mxcohrtdat[0].total==0) xx=1;
			for( int ichrt = 0; ichrt < xx; ++ichrt ) {

				telmnt[0].readCohortState( ifstate, ichrt );
				//telmnt[0].readCohortState0( ifstate, ichrt );
			}
		}

        // soil/elevation data files
		if ( 1 == temflag && 0 == elmnt.strtflag) {
			telmnt[0].fao.getdel(fstxt);
			telmnt[0].elv.getdel(felev);
			//fgets(dummyline, sizeof(dummyline), fstxt);
			//fgets(dummyline, sizeof(dummyline), felev);
	    }

	}

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void updateTCLMGridCell( const int& pdyr )
{

  const double Watts2cal = 1.0 / 0.4845;
  int dyr;
  int dm;
  int gisend;

  double lat = -999.9;

  string SRADname;
   
  fatalerr = 0; 

  if( 0 == pdyr )
  {
    if( 0 == telmnt[0].lonlatflag )
    {
      gisend = latdat.getdel( flonlat );

      if( -1 == gisend ) 
      {
        cout << "Ran out of LONLAT data" << endl << endl;
        flog1 << "Ran out of LONLAT data" << endl << endl;
            
        exit( -1 );
      }
    }

    if( 1 == telmnt[0].clm.tcldsflag )
    {
      if( 0 == telmnt[0].clm.sradflag )
      {
        for( dyr = 0; dyr < (transtime+1); ++dyr )
        { 
          gisend = girrdat[dyr].getdel( ifgirr );
          
          if( -1 == gisend ) 
          {
            cout << "Ran out of GIRR data" << endl << endl;
            flog1 << "Ran out of GIRR data" << endl << endl;
            
            exit( -1 );
          }
        }
      }
    

      if( 1 == telmnt[0].clm.cldflag )
      {
        for( dyr = 0; dyr < (transtime+1); ++dyr )
        { 
          gisend = cldsdat[dyr].getdel( ifclds );
          
          if( -1 == gisend ) 
          {  
            cout << "Ran out of Cloudiness data" << endl << endl;
            flog1 << "Ran out of Cloudiness data" << endl << endl;
          
            exit( -1 );
          }
        }

        telmnt[0].col = cldsdat[0].col;
        telmnt[0].row = cldsdat[0].row;
        SRADname = "CLDINESS";
      }  
      else
      {
        for( dyr = 0; dyr < (transtime+1); ++dyr )
        { 
          gisend = nirrdat[dyr].getdel( ifnirr );
          
          if( -1 == gisend )
          {  
            cout << "Ran out of NIRR data" << endl << endl;
            flog1 << "Ran out of NIRR data" << endl << endl;
  
            exit( -1 );
          }
        }

        telmnt[0].col = nirrdat[0].col;
        telmnt[0].row = nirrdat[0].row;
        SRADname = "NIRR";
      }

      if( 1 == telmnt[0].clm.parflag )
      {
        for( dyr = 0; dyr < (transtime+1); ++dyr )
        { 
          gisend = pardat[dyr].getdel( ifpar );
  
          if( -1 == gisend )
          { 
            cout << "Ran out of PAR data" << endl << endl;
            flog1 << "Ran out of PAR data" << endl << endl;
  
            exit( -1 );
          }
        }
      }
    }
    else  // 0 == telmnt[0].clm.tcldsflag
    { 
      if( 0 == telmnt[0].clm.sradflag )
      {
        gisend = girrdat[0].getdel( ifgirr );
  
        if( -1 == gisend ) 
        { 
          cout << "Ran out of GIRR data" << endl << endl;
          flog1 << "Ran out of GIRR data" << endl << endl;
  
          exit( -1 );
        }
      }

      if( 1 == telmnt[0].clm.cldflag )
      {
        gisend = cldsdat[0].getdel( ifclds );
  
        if( -1 == gisend ) 
        {  
          cout << "Ran out of Cloudiness data" << endl << endl;
          flog1 << "Ran out of Cloudiness data" << endl << endl;
  
          exit( -1 );
        }

        telmnt[0].col = cldsdat[0].col;
        telmnt[0].row = cldsdat[0].row;
        SRADname = "CLDINESS";
      }
      else
      {
        gisend = nirrdat[0].getdel( ifnirr );
  
        if( -1 == gisend ) 
        {  
          cout << "Ran out of NIRR data" << endl << endl;
          flog1 << "Ran out of NIRR data" << endl << endl;
  
          exit( -1 );
        }

        telmnt[0].col = nirrdat[0].col;
        telmnt[0].row = nirrdat[0].row;
        SRADname = "NIRR";
      }

      if( 1 == telmnt[0].clm.parflag )
      {
        gisend = pardat[0].getdel( ifpar );
  
        if( -1 == gisend ) 
        { 
          cout << "Ran out of PAR data" << endl << endl;
          flog1 << "Ran out of PAR data" << endl << endl;
  
          exit( -1 );
        }
      }
    }
    //cout <<"xxxxcol: " <<telmnt[0].col <<" xxxxro: "<<telmnt[0].row<<endl;
    // Look for spatial co-registration problems between 
    //   cloudiness or net irradiance (NIRR) spatially explicit 
    //   data and gross irradiance (GIRR) spatially explicit 
    //   data 

    if( 0 == telmnt[0].clm.sradflag )
    {
      fatalerr = telmnt[0].coregerr( flog1, 
                                     SRADname,
                                     telmnt[0].col,
                                     telmnt[0].row,    
                                     "GIRR",
                                     girrdat[0].col,
                                     girrdat[0].row );
 
      if( fatalerr != 0 ) { exit( -1 ); }
    }

    // Look for spatial co-registration problems between 
    //   cloudiness or net irradiance (NIRR) spatially explicit 
    //   data and photosynthetically active radiation (PAR) 
    //   spatially explicit data 

    if( 1 == telmnt[0].clm.parflag )
    {
      fatalerr = telmnt[0].coregerr( flog1, 
                                     SRADname, 
                                     telmnt[0].col,
                                     telmnt[0].row,    
                                     "PAR",
                                     pardat[0].col,
                                     pardat[0].row );
 
      if( fatalerr != 0 ) { exit( -1 ); }
    }

    // Read in historical monthly air temperatures for grid cell

    if( 1 == telmnt[0].clm.ttairflag )
    {
      for( dyr = 0; dyr < (transtime+1); ++dyr )
      {
        gisend = tairdat[dyr].getdel( iftair );
        
        if( -1 == gisend ) 
        { 
          cout << "Ran out of Air Temperature data" << endl << endl;
          flog1 << "Ran out of Air Temperature data" << endl << endl;
          
          exit( -1 );
        }
      }
    }
    else
    {
      gisend = tairdat[0].getdel( iftair );
      
      if( -1 == gisend ) 
      { 
        cout << "Ran out of Air Temperature data" << endl << endl;
        flog1 << "Ran out of Air Temperature data" << endl << endl;
      
        exit( -1 );
      }
    }

    
    // Look for spatial co-registration problems between cloudiness and
    // air temperature spatially explicit data sets
    
    fatalerr = telmnt[0].coregerr( flog1, 
                                   SRADname, 
                                   telmnt[0].col,
                                   telmnt[0].row,    
                                   "TAIR",
                                   tairdat[0].col,
                                   tairdat[0].row );
 
    if( fatalerr != 0 ) { exit( -1 ); }

    
    // Read in historical monthly precipitation for grid cell
    
    if( 1 == telmnt[0].clm.tprecflag )
    {
      for( dyr = 0; dyr < (transtime+1); ++dyr )
      { 
        gisend = precdat[dyr].getdel( ifprec );
        
        if( -1 == gisend ) 
        { 
          cout << "Ran out of Precipitation data" << endl << endl;
          flog1 << "Ran out of Precipitation data" << endl << endl;
          
          exit( -1 );
        }
      }
    }
    else
    {
      gisend = precdat[0].getdel( ifprec );
      
      if( -1 == gisend ) 
      { 
        cout << "Ran out of Precipitation data" << endl << endl;
        flog1 << "Ran out of Precipitation data" << endl << endl;
      
        exit( -1 );
      }
    }

    
    // Look for spatial co-registration problems between cloudiness and
    // precipitation spatially explicit data sets

    fatalerr = telmnt[0].coregerr( flog1, 
                                   SRADname, 
                                   telmnt[0].col,
                                   telmnt[0].row,    
                                   "PREC",
                                   precdat[0].col,
                                   precdat[0].row );

    if( fatalerr != 0 ) { exit( -1 ); }

  
    // Read in historical annual atmospheric CO2 data for 
    //   globe
    
    if( 1 == telmnt[0].clm.tco2flag && 0 == assignCO2 )
    {
      for( dyr = 0; dyr < (transtime+2); ++dyr )
      {
        co2dat[dyr].get( ifco2 );
      }

      assignCO2 = 1;
    }

    // Read in historical monthly atmospheric CO2 data for
    //   grid cell
    
    if( 2 == telmnt[0].clm.tco2flag )
    {
      cout << "This feature has not been implemented yet ";
      cout << "in this TEM version!"  << endl;

      flog1 << "This feature has not been implemented yet ";
      flog1 << "in this TEM version!"  << endl;

      exit( -1 );
    }

    
    // Read in historical monthly ozone data for grid cell
    
    if( 1 == telmnt[0].clm.to3flag )
    {
      for( dyr = 0; dyr < (transtime+1); ++dyr )
      { 
        gisend = o3dat[dyr].getdel( ifo3 );
        
        if( -1 == gisend ) 
        { 
          cout << "Ran out of Ozone data" << endl << endl;
          flog1 << "Ran out of Ozone data" << endl << endl;
          
          exit( -1 );
        }
      }
    }
    else
    {
      gisend = o3dat[0].getdel( ifo3 );
      
      if( -1 == gisend ) 
      { 
        cout << "Ran out of Ozone data" << endl << endl;
        flog1 << "Ran out of Ozone data" << endl << endl;
        
        exit( -1 );
      }
    }
    
    // Look for spatial co-registration problems between cloudiness and
    // ozone spatially explicit data sets

    fatalerr = telmnt[0].coregerr( flog1, 
                                   SRADname, 
                                   telmnt[0].col,
                                   telmnt[0].row,    
                                   "AOT40",
                                   o3dat[0].col,
                                   o3dat[0].row );

    if( fatalerr != 0 ) { exit( -1 ); }
  


    // Determine monthly N deposition grid cell

    if( 1 == telmnt[0].clm.tndepflag )
    {
      for( dyr = 0; dyr < (transtime+1); ++dyr )
      {
        gisend = ndepdat[dyr].getdel( ifndep );
        
        if( -1 == gisend )
        {
          cout << "Ran out of N deposition data";
          cout << endl << endl;
          flog1 << "Ran out of N deposition data";
          flog1 << endl << endl;
        
          exit( -1 );
        }

        // Check input for valid data
        
        if( ndepdat[dyr].nhx < ZERO ) 
        { 
          flog1 << "Lon = " << ndepdat[dyr].col << " ";
          flog1 << "Lat = " << ndepdat[dyr].row << " ";
          flog1 << "NH4DEP is less than zero during ";
          flog1 << ndepdat[dyr].year << endl;
          
          ndepdat[dyr].nhx = ZERO; 
        }
        
        if( ndepdat[dyr].noy < ZERO ) 
        { 
          flog1 << "Lon = " << ndepdat[dyr].col << " ";
          flog1 << "Lat = " << ndepdat[dyr].row << " ";
          flog1 << "NO3DEP is less than zero during ";
          flog1 << ndepdat[dyr].year << endl;

          ndepdat[dyr].noy = ZERO; 
        }
        
        if( ndepdat[dyr].totndep 
            != ndepdat[dyr].nhx + ndepdat[dyr].noy )
        {
          ndepdat[dyr].totndep = ndepdat[dyr].nhx 
                                 + ndepdat[dyr].noy;
        }

        // Adjust NDEP units from mg N m-2 mo-1 to 
        //   g N m-2 mo-1 for use in TEM

        ndepdat[dyr].totndep *= 0.001;
        ndepdat[dyr].nhx *= 0.001;
        ndepdat[dyr].noy *= 0.001;
      }
    }
    else
    {
      gisend = ndepdat[0].getdel( ifndep );

      if( -1 == gisend )
      {
        cout << "Ran out of N deposition data";
        cout << endl << endl;
        flog1 << "Ran out of N deposition data";
        flog1 << endl << endl;
        
        exit( -1 );
      }

      // Check input for valid data

      if( ndepdat[0].nhx < ZERO ) 
      { 
        flog1 << "Lon = " << ndepdat[0].col << " ";
        flog1 << "Lat = " << ndepdat[0].row << " ";
        flog1 << "NH4DEP is less than zero during ";
        flog1 << ndepdat[0].year << endl;
          
        ndepdat[0].nhx = ZERO; 
      }
        
      if( ndepdat[0].noy < ZERO ) 
      { 
        flog1 << "Lon = " << ndepdat[0].col << " ";
        flog1 << "Lat = " << ndepdat[0].row << " ";
        flog1 << "NO3DEP is less than zero during ";
        flog1 << ndepdat[0].year << endl;

        ndepdat[0].noy = ZERO; 
      }
        
      if( ndepdat[0].totndep 
          != ndepdat[0].nhx + ndepdat[0].noy )
      {
        ndepdat[0].totndep = ndepdat[0].nhx + ndepdat[0].noy;
      }

      // Adjust NDEP units from mg N m-2 mo-1 to 
      //   g N m-2 mo-1 for use in TEM

      ndepdat[0].totndep *= 0.001;
      ndepdat[0].nhx *= 0.001;
      ndepdat[0].noy *= 0.001;
    }

    // Look for spatial co-registration problems between cloudiness and
    // N deposition spatially explicit data sets

    fatalerr = telmnt[0].coregerr( flog1, 
                                   SRADname, 
                                   telmnt[0].col,
                                   telmnt[0].row,    
                                   "NDEP",
                                   ndepdat[0].col,
                                   ndepdat[0].row );

    if( fatalerr != 0 ) { exit( -1 ); }
  }

  if( 0 == pdyr ) { dyr = 0; }
  else if( istateflag < 2 && pdyr < (totsptime+1) )
  {
    dyr = (pdyr-1)%spintime + 1;
  }
  else
  {
    if( istateflag < 2 ) { dyr = pdyr - totsptime; }
    else{ dyr = pdyr; }
  }


  if( 0 == telmnt[0].clm.tcldsflag )
  {
    if( 0 == telmnt[0].clm.sradflag )
    {
      for( dm = 0; dm < (CYCLE+1); ++dm )
      {
        girrdat[dyr].mon[dm] = girrdat[0].mon[dm];
      }
    }

    if( 1 == telmnt[0].clm.cldflag )
    {
      for( dm = 0; dm < (CYCLE+1); ++dm )
      {
        cldsdat[dyr].mon[dm] = cldsdat[0].mon[dm];
      }
    }
    else
    {
      for( dm = 0; dm < (CYCLE+1); ++dm )
      {
        nirrdat[dyr].mon[dm] = nirrdat[0].mon[dm];
      }
    }

    if( 1 == telmnt[0].clm.parflag )
    {
      for( dm = 0; dm < (CYCLE+1); ++dm )
      {
        pardat[dyr].mon[dm] = pardat[0].mon[dm];
      }
    }
  }

  if( 0 == telmnt[0].clm.ttairflag )
  {
    tairdat[dyr].max = tairdat[0].max;
    tairdat[dyr].ave = tairdat[0].ave;

    for( dm = 0; dm < (CYCLE+1); ++dm )
    {
      tairdat[dyr].mon[dm] = tairdat[0].mon[dm];
    }
  }

  if( 0 == telmnt[0].clm.tprecflag )
  {
    precdat[dyr].total = precdat[0].total;

    for( dm = 0; dm < (CYCLE+1); ++dm )
    {
      precdat[dyr].mon[dm] = precdat[0].mon[dm];
    }
  }

  if( 0 == telmnt[0].clm.to3flag )
  {
    for( dm = 0; dm < (CYCLE+1); ++dm )
    {
      o3dat[dyr].mon[dm] = o3dat[0].mon[dm];
    }
  }

  if( 0 == telmnt[0].clm.tndepflag )
  {
    ndepdat[dyr].totndep = ndepdat[0].totndep;
    ndepdat[dyr].nhx = ndepdat[0].nhx;
    ndepdat[dyr].noy = ndepdat[0].noy;
  }


  // Determine maximum monthly air temperature for year for 
  //   each grid cell
    
  telmnt[0].mxtair = tairdat[dyr].max;
  telmnt[0].avetair = tairdat[dyr].ave;


  // Determine annual precipitation for each grid cell
    
  telmnt[0].yrprec = precdat[dyr].total;


  // Interpolate annual atmospheric CO2 concentrations
  //   to a monthly temporal resolution for grid cell
  //   Value in co2dat assumed to represent July atmospheric 
  //   CO2 concentrations    
  
  if( 1 == telmnt[0].clm.tco2flag && dyr > 0 )
  { 
    if( pdyr < (totsptime+1) )
    {
      for( dm = 0; dm < (CYCLE+1); ++dm ) 
      { 
        telmnt[0].climate[telmnt[0].clm.I_CO2][dm] = co2dat[0].mco2; 
      }
    }
    else
    {
      for( dm = 0; dm < CYCLE; ++dm )
      {
        if( dm < 6 )
        {
          telmnt[0].climate[telmnt[0].clm.I_CO2][dm] = co2dat[dyr-1].mco2
                                                       + ((dm + 6) 
                                                       * (co2dat[dyr].mco2 
                                                       - co2dat[dyr-1].mco2)
                                                       / (double) CYCLE);
        }
        else
        {
          telmnt[0].climate[telmnt[0].clm.I_CO2][dm] = co2dat[dyr].mco2
                                                       + ((dm - 6) 
                                                       * (co2dat[dyr+1].mco2 
                                                       - co2dat[dyr].mco2)
                                                       / (double) CYCLE);
        }
      }
    }
    
    telmnt[0].climate[telmnt[0].clm.I_CO2][CYCLE] = telmnt[0].climate[telmnt[0].clm.I_CO2][CYCLE-1];

//    co2dat[dyr].year = telmnt[0].clm.co2year[pdyr] ;
  }
  else // 0 == telmnt[0].clm.tco2flag || 0 == dyr 
  {
    for( dm = 0; dm < (CYCLE+1); ++dm )
    {
      telmnt[0].climate[telmnt[0].clm.I_CO2][dm] = telmnt[0].clm.getCO2LEVEL();
    }
  }


 // Assign information from files to telmnt.climate

  for( dm = 0; dm < (CYCLE+1); ++dm )
  {
    // Air temperature
    
    telmnt[0].climate[telmnt[0].clm.I_TAIR][dm] = tairdat[dyr].mon[dm];

    // Precipitation
    
    telmnt[0].climate[telmnt[0].clm.I_PREC][dm] = precdat[dyr].mon[dm];
    
    // Rain and Snowfall
    
    telmnt[0].clm.precsplt( telmnt[0].climate[telmnt[0].clm.I_PREC][dm], 
                            telmnt[0].climate[telmnt[0].clm.I_TAIR][dm],
                            telmnt[0].climate[telmnt[0].clm.I_RAIN][dm], 
                            telmnt[0].climate[telmnt[0].clm.I_SNWFAL][dm] );

    // AOT40 ozone index
    
    telmnt[0].climate[telmnt[0].clm.I_AOT40][dm] = o3dat[dyr].mon[dm];

    // Calculate monthly N deposition based on proportion of 
    //   annual precipitation that occurs during a particular 
    //   month

    if( precdat[dyr].total > 0.01 )
    {
      telmnt[0].climate[telmnt[0].clm.I_TNDEP][dm] = ndepdat[dyr].totndep
                                                     * precdat[dyr].mon[dm]
                                                     / precdat[dyr].total;

      telmnt[0].climate[telmnt[0].clm.I_NH4DEP][dm] = ndepdat[dyr].nhx
                                                      * precdat[dyr].mon[dm]
                                                      / precdat[dyr].total;

      telmnt[0].climate[telmnt[0].clm.I_NO3DEP][dm] = ndepdat[dyr].noy
                                                      * precdat[dyr].mon[dm]
                                                      / precdat[dyr].total;

    }
    else
    {
      telmnt[0].climate[telmnt[0].clm.I_TNDEP][dm] = ZERO;
      telmnt[0].climate[telmnt[0].clm.I_NH4DEP][dm] = ZERO;
      telmnt[0].climate[telmnt[0].clm.I_NO3DEP][dm] = ZERO;
    }
  }


//**************************************************************

    
  // Calculate GIRR during first year of simulation 
  //   (Note: use same values throughout simulation) 

  if( 1 == telmnt[0].clm.sradflag )
  {
    if( 0 == pdyr )
    {
      if( 0 == telmnt[0].lonlatflag ) { lat = latdat.lat; } 
      else { lat = telmnt[0].row; }

      telmnt[0].clm.yrsumday = ZERO;
      
      for( dm = 0; dm < CYCLE; ++dm )
      {
        telmnt[0].climate[telmnt[0].clm.I_GIRR][dm] = telmnt[0].clm.xgirr( lat, 
                                                                           dm, 
                                                                           telmnt[0].clm.yrsumday );
      }

      telmnt[0].climate[telmnt[0].clm.I_GIRR][CYCLE] = telmnt[0].climate[telmnt[0].clm.I_GIRR][0];
    }
  }
  else 
  { 
    for( dm = 0; dm < (CYCLE+1); ++dm )
    {  
      telmnt[0].climate[telmnt[0].clm.I_GIRR][dm] = girrdat[dyr].mon[dm] * Watts2cal;
    }
  }

  // Calculate NIRR, CLDINESS and PAR or retrieve from earlier calculations

  for( dm = 0; dm < (CYCLE+1); ++dm )
  {
    if( 1 == telmnt[0].clm.cldflag )
    {
      telmnt[0].climate[telmnt[0].clm.I_CLDS][dm] = cldsdat[dyr].mon[dm];

      telmnt[0].climate[telmnt[0].clm.I_NIRR][dm] = telmnt[0].clm.xnirr( telmnt[0].climate[telmnt[0].clm.I_CLDS][dm], 
	                                                                 telmnt[0].climate[telmnt[0].clm.I_GIRR][dm] );
    }
    else
    {
      telmnt[0].climate[telmnt[0].clm.I_NIRR][dm] = nirrdat[dyr].mon[dm] * Watts2cal;

      telmnt[0].climate[telmnt[0].clm.I_CLDS][dm] = telmnt[0].clm.mkclds( telmnt[0].climate[telmnt[0].clm.I_GIRR][dm], 
                                                                          telmnt[0].climate[telmnt[0].clm.I_NIRR][dm] );
    }
      
    if( 0 == telmnt[0].clm.parflag )
    {
      
      telmnt[0].climate[telmnt[0].clm.I_PAR][dm]  = telmnt[0].clm.xpar( telmnt[0].climate[telmnt[0].clm.I_CLDS][dm], 
                                                                        telmnt[0].climate[telmnt[0].clm.I_NIRR][dm] );
    }
    else
    {
      telmnt[0].climate[telmnt[0].clm.I_PAR][dm] = pardat[dyr].mon[dm] * Watts2cal;
    }
  }

  // Set year 

  telmnt[0].year = telmnt[0].clm.modstartyr
                   - totsptime 
                   - 1 
                   + pdyr;
  //#ifdef PMODE
  //cout <<" processor: " <<nthproc<< " year: " <<telmnt[0].year <<" col: " <<telmnt[0].col<<" row: " <<telmnt[0].row<<" tair: " << telmnt[0].climate[telmnt[0].clm.I_TAIR][1]<<endl;
  //#endif
    
  // Copy TEMclm results to output variables

  if( 1 == telmnt[0].clm.predflag )
  {
    if ( (1 == spinoutfg 
         && telmnt[0].year < clmstartyr)
         || (2 == spinoutfg 
         && telmnt[0].year >= (clmstartyr-spinoutyrs))
         || (telmnt[0].year >= clmstartyr) )
    {
      if( 1 == telmnt[0].clm.cldflag )
      {
        //telmnt[0].carea = cldsdat[0].carea; //closed by cgs2014. instead, the carea is obtained from the maxcohort data.
        //telmnt[0].contnent = cldsdat[0].contnent;
      }
      else
      {
        //telmnt[0].carea = nirrdat[0].carea;
        //telmnt[0].contnent = nirrdat[0].contnent;
      }
    
      telmnt[0].atmswritepred( fclmpred, 
                               clmpredmap, 
                               telmnt[0].natmspred );
    }
  }
	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void updateTLCLUCGridCell( const int& pdyr )
{

  int dyr;
  int tstyr=0;
  int gisend;
  int ichrt;

//  gridlulc.gotoFirstNode();   // start at first LULCListNode

  fatalerr = 0;

  if( 0 == pdyr )
  {
    // Get the total number of cohorts in the grid cell

    if( 1 == telmnt[0].lcluc.tlulcflag )
    {
      for( dyr = 0; dyr < (transtime+1); ++dyr )
      {
        gisend = mxcohrtdat[dyr].getdel( ifnumchrts );


        // Check data for spatial coregistration errors

        fatalerr = telmnt[0].coregerr( flog1,
                                       "Climate",
                                       telmnt[0].col,
                                       telmnt[0].row,
                                       "MAXCOHORTS",
                                       mxcohrtdat[dyr].col,
                                       mxcohrtdat[dyr].row );
        //cout <<" "
        if( fatalerr != 0 ) { exit( -1 ); }
      }

      if( -1 == gisend ) 
      { 
        cout << "Ran out of Number of Cohorts data";
        cout << endl << endl;
        flog1 << "Ran out of Number of Cohorts data";
        flog1 << endl << endl;
        
        exit( -1 );
      } 
    }
    else // 0 == telmnt[0].lcluc.tlulcflag
    {
      gisend = mxcohrtdat[0].getdel( ifnumchrts );
      //cout <<"total: "<<mxcohrtdat[dyr].total<<" carea: "<<mxcohrtdat[dyr].carea<<" cohortnum: "<<mxcohrtdat[0].natchrts<<endl;
      if( -1 == gisend ) 
      { 
        cout << "Ran out of Number of Cohorts data";
        cout << endl << endl;
        flog1 << "Ran out of Number of Cohorts data";
        flog1 << endl << endl;
        
        exit( -1 );
      } 

      // Check data for spatial coregistration errors
      //cout <<"col: " <<mxcohrtdat[0].col <<endl;
      //cout <<"row: " <<mxcohrtdat[0].row <<endl;

      fatalerr = telmnt[0].coregerr( flog1, 
                                     "Climate", 
                                     telmnt[0].col,  
                                     telmnt[0].row, 
                                     "MAXCOHORTS", 
                                     mxcohrtdat[0].col, 
                                     mxcohrtdat[0].row );

      if( fatalerr != 0 ) { exit( -1 ); }


      for( dyr = 1; dyr < (transtime+1); ++dyr )
      {
        mxcohrtdat[dyr].col = mxcohrtdat[0].col;

        mxcohrtdat[dyr].row = mxcohrtdat[0].row;
        mxcohrtdat[dyr].carea = mxcohrtdat[0].carea;
        mxcohrtdat[dyr].year = mxcohrtdat[0].year + dyr;
        mxcohrtdat[dyr].total = mxcohrtdat[0].total;
        mxcohrtdat[dyr].natchrts = mxcohrtdat[0].natchrts;
        mxcohrtdat[dyr].contnent = mxcohrtdat[0].contnent;
      }
    }
      	       

    // Get land use/land cover cohort data for  the grid cell 

    if( 1 == telmnt[0].lcluc.tlulcflag )
    {
      for( dyr = 0; dyr < (transtime+1); ++dyr )
      {
    	if (mxcohrtdat[dyr].total == 0)//added by cgs
    	{
    	    ichrt=0;
    	    lulcdat[dyr][ichrt].col = telmnt[0].col;
    	    lulcdat[dyr][ichrt].row = telmnt[0].col;
    	    lulcdat[dyr][ichrt].varname = "LULCHRT";
    	    lulcdat[dyr][ichrt].year = telmnt[0].year + dyr;
    	    lulcdat[dyr][ichrt].icohort = 1;
    	    lulcdat[dyr][ichrt].isrccohort = 1;
    	    lulcdat[dyr][ichrt].standage = 0;
    	    lulcdat[dyr][ichrt].chrtarea = -9999;
    	    lulcdat[dyr][ichrt].potveg = 1;
    	    lulcdat[dyr][ichrt].currentveg = 1;
    	    lulcdat[dyr][ichrt].subtype = 1;
    	    lulcdat[dyr][ichrt].agstate = 0;
    	    lulcdat[dyr][ichrt].agprevstate = 0;
    	    lulcdat[dyr][ichrt].tillflag = 0;
    	    lulcdat[dyr][ichrt].fertflag = 0;
    	    lulcdat[dyr][ichrt].irrgflag = 0;
    	    lulcdat[dyr][ichrt].disturbflag = 0;
    	    lulcdat[dyr][ichrt].disturbmonth = 0;
    	    lulcdat[dyr][ichrt].FRI = 2000;
    	    lulcdat[dyr][ichrt].leafmortpar =0.0;
    	    lulcdat[dyr][ichrt].stemmortpar=0.0;
    	    lulcdat[dyr][ichrt].rootmortpar=0.0;
    	    lulcdat[dyr][ichrt].leafslash =0.0;
    	    lulcdat[dyr][ichrt].leafconv =0.0;
    	    lulcdat[dyr][ichrt].rootslash =0.0;
    	    lulcdat[dyr][ichrt].rootconv =0.0;
    	    lulcdat[dyr][ichrt].stemslash =0.0;
    	    lulcdat[dyr][ichrt].stemconv =0.0;
    	    lulcdat[dyr][ichrt].standdead =0.0;
    	    lulcdat[dyr][ichrt].deadconv =0.0;
    	    lulcdat[dyr][ichrt].deadslash =0.0;
    	    lulcdat[dyr][ichrt].sconvert =0.0;
    	    lulcdat[dyr][ichrt].prod10par =0.0;
    	    lulcdat[dyr][ichrt].prod100par =0.0;
    	    lulcdat[dyr][ichrt].forestage =0;
    	    //lulcdat[dyr][ichrt].slashpar = 0;
    	    //lulcdat[dyr][ichrt].vconvert = 0;
    	    //lulcdat[dyr][ichrt].prod10par = 0;
    	    //lulcdat[dyr][ichrt].prod100par = 0;
    	    //lulcdat[dyr][ichrt].vrespar = 0;
    	    //lulcdat[dyr][ichrt].sconvert = 0;
    	    //lulcdat[dyr][ichrt].deadconvert = 0;
    	    //lulcdat[dyr][ichrt].deadslash = 0;
    	    lulcdat[dyr][ichrt].region = "other";

    	}

    	for( ichrt = 0; ichrt < mxcohrtdat[dyr].total; ++ichrt )
        {
          gisend = lulcdat[dyr][ichrt].getdel( iflulc );

          if( -1 == gisend ) 
          { 
            flog1 << "Ran out of Land cover/land use data";
            flog1 << endl << endl;
            
            exit( -1 );
          }

          // Check data for spatial coregistration errors
    
          fatalerr = telmnt[0].coregerr( flog1, 
                                         "Climate", 
                                         telmnt[0].col,  
                                         telmnt[0].row, 
                                         "LULC", 
                                         lulcdat[dyr][ichrt].col, 
                                         lulcdat[dyr][ichrt].row );
          //if (nthproc==4)
          //if (telmnt[0].col == -164.25 && telmnt[0].row==66.00) cout <<" col1: "<<telmnt[0].col<<" row1: "<< telmnt[0].row <<" col2: "<< lulcdat[dyr][ichrt].col<<" row2: "<< lulcdat[dyr][ichrt].row<<" year: "<<dyr+1800<<" lulcyear: "<<lulcdat[dyr][ichrt].year<<" clmyear: "<<nirrdat[dyr].year<<" maxlulcyear: "<<mxcohrtdat[dyr].year<<endl;
          if( fatalerr != 0 ) 
          { 
            //cout << "year = " << lulcdat[dyr][ichrt].year;
        	cout << "year = " << dyr+1800;
            cout << "total: "<<mxcohrtdat[dyr].total<< " icohort = " << lulcdat[dyr][ichrt].icohort;
            cout << "FATAL ERROR! Program Terminated - on Process: " <<nthproc << endl;
            cout << endl;

            exit( -1 );
          }
        }
      }
    }
    else // 0 == telmnt[0].lcluc.tlulcflag
    {
        if (mxcohrtdat[dyr].total == 0) //added by cgs2014
        {
      	    ichrt=0;
      	    for( dyr = 0; dyr < (transtime+1); ++dyr )
      	    {
      	  	lulcdat[dyr][ichrt].col = telmnt[0].col;
          	lulcdat[dyr][ichrt].row = telmnt[0].row;
          	lulcdat[dyr][ichrt].varname = "LULCHRT";
          	lulcdat[dyr][ichrt].year = telmnt[0].year + dyr;
          	lulcdat[dyr][ichrt].icohort = 1;
          	lulcdat[dyr][ichrt].isrccohort = 1;
          	lulcdat[dyr][ichrt].standage = 0;
          	lulcdat[dyr][ichrt].chrtarea = -9999;
          	lulcdat[dyr][ichrt].potveg = 1;
          	lulcdat[dyr][ichrt].currentveg = 1;
          	lulcdat[dyr][ichrt].subtype = 1;
          	lulcdat[dyr][ichrt].agstate = 0;
          	lulcdat[dyr][ichrt].agprevstate = 0;
          	lulcdat[dyr][ichrt].tillflag = 0;
          	lulcdat[dyr][ichrt].fertflag = 0;
          	lulcdat[dyr][ichrt].irrgflag = 0;
          	lulcdat[dyr][ichrt].disturbflag = 0;
          	lulcdat[dyr][ichrt].disturbmonth = 0;
          	lulcdat[dyr][ichrt].FRI = 2000;
          	lulcdat[dyr][ichrt].leafmortpar =0.0;
          	lulcdat[dyr][ichrt].stemmortpar=0.0;
          	lulcdat[dyr][ichrt].rootmortpar=0.0;
          	lulcdat[dyr][ichrt].leafslash =0.0;
          	lulcdat[dyr][ichrt].leafconv =0.0;
          	lulcdat[dyr][ichrt].rootslash =0.0;
          	lulcdat[dyr][ichrt].rootconv =0.0;
          	lulcdat[dyr][ichrt].stemslash =0.0;
          	lulcdat[dyr][ichrt].stemconv =0.0;
          	lulcdat[dyr][ichrt].standdead =0.0;
          	lulcdat[dyr][ichrt].deadconv =0.0;
          	lulcdat[dyr][ichrt].deadslash =0.0;
          	lulcdat[dyr][ichrt].sconvert =0.0;
          	lulcdat[dyr][ichrt].prod10par =0.0;
          	lulcdat[dyr][ichrt].prod100par =0.0;
          	lulcdat[dyr][ichrt].forestage =0;
          	lulcdat[dyr][ichrt].region = "other";
      	    }

      }

      for( ichrt = 0; ichrt < mxcohrtdat[0].total; ++ichrt )
      {
        gisend = lulcdat[0][ichrt].getdel( iflulc );
        //cout << " col: "<<lulcdat[0][ichrt].col<<" row: "<<lulcdat[0][ichrt].row<< " total: "<<mxcohrtdat[0].total<< " icohort: " << lulcdat[0][ichrt].icohort<<" currentveg: "<<lulcdat[0][ichrt].currentveg <<" tyear: "<<lulcdat[0][ichrt].year<<" standage: "<<lulcdat[0][ichrt].standage<<endl;
        if( -1 == gisend ) 
        { 
          flog1 << "Ran out of Land cover/land use data";
          flog1 << endl << endl;
          
          exit( -1 );
        }

        fatalerr = telmnt[0].coregerr( flog1, 
                                       "Climate", 
                                       telmnt[0].col,  
                                       telmnt[0].row, 
                                       "LULC", 
                                       lulcdat[0][ichrt].col, 
                                       lulcdat[0][ichrt].row );

        if( fatalerr != 0 )
        {

        	cout << "year = " << dyr;
            cout << " total: "<<mxcohrtdat[dyr].total<< " icohort = " << lulcdat[0][ichrt].icohort<<" currentveg: "<<lulcdat[0][ichrt].currentveg <<" tyear: "<<lulcdat[0][ichrt].year;
            cout << endl;

        	exit( -1 );
        }

        for( dyr = 1; dyr < (transtime+1); ++dyr )
        {
          lulcdat[dyr][ichrt].year = lulcdat[0][ichrt].year + dyr;
          lulcdat[dyr][ichrt].isrccohort = lulcdat[0][ichrt].isrccohort;
          lulcdat[dyr][ichrt].standage = lulcdat[0][ichrt].standage;
          lulcdat[dyr][ichrt].chrtarea = lulcdat[0][ichrt].chrtarea;
          lulcdat[dyr][ichrt].potveg = lulcdat[0][ichrt].potveg;
          lulcdat[dyr][ichrt].currentveg = lulcdat[0][ichrt].currentveg;
          lulcdat[dyr][ichrt].subtype = lulcdat[0][ichrt].subtype;
          lulcdat[dyr][ichrt].agstate = lulcdat[0][ichrt].agstate;
          lulcdat[dyr][ichrt].agprevstate = lulcdat[0][ichrt].agprevstate;
          lulcdat[dyr][ichrt].tillflag = lulcdat[0][ichrt].tillflag;                              
          lulcdat[dyr][ichrt].fertflag = lulcdat[0][ichrt].fertflag;                              
          lulcdat[dyr][ichrt].irrgflag = lulcdat[0][ichrt].irrgflag;
          lulcdat[dyr][ichrt].disturbflag = lulcdat[0][ichrt].disturbflag;
          lulcdat[dyr][ichrt].disturbmonth = lulcdat[0][ichrt].disturbmonth;                              
          lulcdat[dyr][ichrt].FRI = lulcdat[0][ichrt].FRI;
          lulcdat[dyr][ichrt].leafmortpar = lulcdat[0][ichrt].leafmortpar;
          lulcdat[dyr][ichrt].rootmortpar = lulcdat[0][ichrt].rootmortpar;
          lulcdat[dyr][ichrt].stemmortpar = lulcdat[0][ichrt].stemmortpar;
          lulcdat[dyr][ichrt].leafslash= lulcdat[0][ichrt].leafslash;
          lulcdat[dyr][ichrt].leafconv= lulcdat[0][ichrt].leafconv;
          lulcdat[dyr][ichrt].rootslash= lulcdat[0][ichrt].rootslash;
          lulcdat[dyr][ichrt].rootconv= lulcdat[0][ichrt].rootconv;
          lulcdat[dyr][ichrt].stemslash= lulcdat[0][ichrt].stemslash;
          lulcdat[dyr][ichrt].stemconv= lulcdat[0][ichrt].stemconv;
          lulcdat[dyr][ichrt].standdead= lulcdat[0][ichrt].standdead;
          lulcdat[dyr][ichrt].deadconv= lulcdat[0][ichrt].deadconv;
          lulcdat[dyr][ichrt].deadslash= lulcdat[0][ichrt].deadslash;
          lulcdat[dyr][ichrt].sconvert= lulcdat[0][ichrt].sconvert;
          lulcdat[dyr][ichrt].prod10par= lulcdat[0][ichrt].prod10par;
          lulcdat[dyr][ichrt].prod100par= lulcdat[0][ichrt].prod100par;
          lulcdat[dyr][ichrt].forestage= lulcdat[0][ichrt].forestage;
          lulcdat[dyr][ichrt].region = lulcdat[0][ichrt].region;

          if (lulcdat[dyr][ichrt].currentveg==50)
          {
        	  lulcdat[dyr][ichrt].potveg=13;
        	  lulcdat[dyr][ichrt].subtype=13;
        	  lulcdat[dyr][ichrt].currentveg=13;
        	  lulcdat[dyr][ichrt].agstate=0;
        	  lulcdat[dyr][ichrt].agprevstate=0;
        	  lulcdat[dyr][ichrt].tillflag=0;
        	  lulcdat[dyr][ichrt].fertflag=0;
        	  lulcdat[dyr][ichrt].irrgflag=0;

          }

        }
        //below code: modified to change site specific land use
        //if (telmnt[0].col==-98.75 && telmnt[0].row==55.75)  //fire disturbance in 1981 at lat, lon.

      }
    }      
  } // end of 0 == pdyr
    
  if( 0 == pdyr ) { tstyr = 0; }
  else if( istateflag < 2 && pdyr < (totsptime+1) )
  {
    tstyr = 1;
  }
  else
  {
    if( istateflag < 2 )
    {
      tstyr = pdyr - totsptime;
    }
    else { tstyr = pdyr; }
  }

  // Check data for temporal coregistration errors in mxcohrtdat
    
  if( 1 == telmnt[0].lcluc.tlulcflag
      && pdyr >= (totsptime+1)
      && telmnt[0].year != mxcohrtdat[tstyr].year )
  { 
    cout << " Year in CLM data does not match ";
    cout << " Year in MAXCOHORTS data" << endl;
    cout << " at Lon = " << telmnt[0].col;
    cout << "  Lat = " << telmnt[0].row << endl;
    cout << "  CLM year = " << telmnt[0].year;
    cout << "  MAXCOHORTS year = " << mxcohrtdat[tstyr].year;
    cout << endl << endl;
      
    flog1 << " Year in CLM data does not match ";
    flog1 << " Year in MAXCOHORTS data" << endl << endl;
    flog1 << " at Lon = " << telmnt[0].col;
    flog1 << "  Lat = " << telmnt[0].row << endl;
    flog1 << " CLM year = " << telmnt[0].year;
    flog1 << " MAXCOHORTS year = " << mxcohrtdat[tstyr].year;
    flog1 << endl << endl;

    exit( -1 );
  }


  // Pass mxcohortdat information to telmnt[0]
  if (mxcohrtdat[tstyr].total == 0)
  {
	  //mxcohrtdat[tstyr].total=1; //added by cgs2014 to omit the grid with total==0
  }
  if (mxcohrtdat[tstyr].carea < 0)
  {
  	 mxcohrtdat[tstyr].carea= -999; //added by cgs2014
  }
  telmnt[0].maxcohorts = mxcohrtdat[tstyr].total;
  telmnt[0].natcohorts = mxcohrtdat[tstyr].natchrts;
  telmnt[0].carea = mxcohrtdat[tstyr].carea;
  telmnt[0].contnent = mxcohrtdat[tstyr].contnent;
  telmnt[0].region = mxcohrtdat[tstyr].contnent;

  if( 0 == tstyr 
      || (istateflag < 2 && pdyr > 1 && pdyr <= (totsptime+1)) )
  {
    telmnt[0].prvmxcohrts = mxcohrtdat[tstyr].total;
  }
  else
  {
    telmnt[0].prvmxcohrts = mxcohrtdat[tstyr-1].total;
  }

 
  for( ichrt = 0; ichrt < telmnt[0].maxcohorts; ++ichrt )
  {
    // Check data for temporal coregistration errors in lulcdat

    if( 1 == telmnt[0].lcluc.tlulcflag
        && pdyr >= (totsptime+1)
        && telmnt[0].year != lulcdat[tstyr][ichrt].year )
    { 
      cout << " Year in CLM data does not match ";
      cout << " Year in LCLUC data" << endl;
      cout << " at Lon = " << telmnt[0].col;
      cout << "  Lat = " << telmnt[0].row << endl;
      cout << "  CLM year = " << telmnt[0].year;
      cout << "  LCLUC year = " << lulcdat[tstyr][ichrt].year;
      cout << " for cohort " << (ichrt+1);
      cout << endl << endl;
      
      flog1 << " Year in CLM data does not match ";
      flog1 << " Year in LCLUC data" << endl << endl;
      flog1 << " at Lon = " << telmnt[0].col;
      flog1 << "  Lat = " << telmnt[0].row << endl;
      flog1 << " CLM year = " << telmnt[0].year;
      flog1 << " LCLUC year = " << lulcdat[tstyr][0].year;
      cout << " for cohort " << (ichrt+1);
      flog1 << endl << endl;

      exit( -1 );
    }
    //if (telmnt[0].col == -164.25 && telmnt[0].row==66.00) cout <<" col1: "<<telmnt[0].col<<" row1: "<< "telmnt[0].year: "<<telmnt[0].year<<" lulcyear: "<<lulcdat[tstyr][ichrt].year<<" maxlulcyear: "<< mxcohrtdat[tstyr].year<<endl;

    // Pass lulcdat information to telmnt[0].cohort
    
    telmnt[0].cohort[ichrt].srcCohort = lulcdat[tstyr][ichrt].isrccohort;
    telmnt[0].cohort[ichrt].standage = lulcdat[tstyr][ichrt].standage;
    telmnt[0].cohort[ichrt].forestage = lulcdat[tstyr][ichrt].forestage;
    if (pdyr ==0) telmnt[0].cohort[ichrt].standage=1000; //assign forestage as 1000 during equilibrium run
    if (pdyr ==0) telmnt[0].cohort[ichrt].forestage=1000; //assign forestage as 1000 during equilibrium run

    if (telmnt[0].carea <0) telmnt[0].cohort[ichrt].chrtarea =0.0; //added by cgs2014
    else telmnt[0].cohort[ichrt].chrtarea = lulcdat[tstyr][ichrt].chrtarea;
    telmnt[0].cohort[ichrt].potveg = lulcdat[tstyr][ichrt].potveg;
    telmnt[0].cohort[ichrt].currentveg = lulcdat[tstyr][ichrt].currentveg;
    telmnt[0].cohort[ichrt].subtype = lulcdat[tstyr][ichrt].subtype;
    telmnt[0].cohort[ichrt].agstate = lulcdat[tstyr][ichrt].agstate;
    telmnt[0].cohort[ichrt].agprvstate = lulcdat[tstyr][ichrt].agprevstate;   
    telmnt[0].cohort[ichrt].tillflag = lulcdat[tstyr][ichrt].tillflag;                           
    telmnt[0].cohort[ichrt].fertflag = lulcdat[tstyr][ichrt].fertflag;                              
    telmnt[0].cohort[ichrt].irrgflag = lulcdat[tstyr][ichrt].irrgflag;                              
    telmnt[0].cohort[ichrt].disturbflag = lulcdat[tstyr][ichrt].disturbflag;                              
    telmnt[0].cohort[ichrt].disturbmonth = lulcdat[tstyr][ichrt].disturbmonth;

    //cgs2014 used to mask out cropland
    if (telmnt[0].cohort[ichrt].currentveg==50)
    {
    	telmnt[0].cohort[ichrt].potveg = 13;
    	telmnt[0].cohort[ichrt].currentveg = 13;
    	telmnt[0].cohort[ichrt].subtype = 13;
    	telmnt[0].cohort[ichrt].agstate = 0;
    	telmnt[0].cohort[ichrt].agprvstate = 0;
    	telmnt[0].cohort[ichrt].tillflag = 0;
    	telmnt[0].cohort[ichrt].fertflag = 0;
    	telmnt[0].cohort[ichrt].irrgflag = 0;
    }


    // If cohort is in agriculture the during first year of the 
    //   historical study period, convert cohort to agriculture 
    //   during the first year of the transient spinup.  Keep
    //   cohort in agriculture throughout spinup period

    if( 1 == pdyr 
        && istateflag < 2 
        && telmnt[0].cohort[ichrt].agstate > 0
        && telmnt[0].cohort[ichrt].agprvstate > 0 )
    {
      telmnt[0].cohort[ichrt].agprvstate = 0;
      telmnt[0].cohort[ichrt].disturbflag = 1;
      telmnt[0].cohort[ichrt].disturbmonth = 1;
    }
 
    telmnt[0].cohort[ichrt].FRI = lulcdat[tstyr][ichrt].FRI;
    //telmnt[0].cohort[ichrt].slashpar = lulcdat[tstyr][ichrt].slashpar;
    telmnt[0].cohort[ichrt].leafmortpar = lulcdat[tstyr][ichrt].leafmortpar;
    telmnt[0].cohort[ichrt].rootmortpar = lulcdat[tstyr][ichrt].rootmortpar;
    telmnt[0].cohort[ichrt].stemmortpar = lulcdat[tstyr][ichrt].stemmortpar;
    //telmnt[0].cohort[ichrt].standdead = lulcdat[tstyr][ichrt].standdead;
    //telmnt[0].cohort[ichrt].vconvert = lulcdat[tstyr][ichrt].vconvert;
    //telmnt[0].cohort[ichrt].prod10par = lulcdat[tstyr][ichrt].prod10par;
    //telmnt[0].cohort[ichrt].prod100par = lulcdat[tstyr][ichrt].prod100par;
    //telmnt[0].cohort[ichrt].vrespar = lulcdat[tstyr][ichrt].vrespar;
    //telmnt[0].cohort[ichrt].sconvert = lulcdat[tstyr][ichrt].sconvert;
    //telmnt[0].cohort[ichrt].deadconvert = lulcdat[tstyr][ichrt].deadconvert;
    //telmnt[0].cohort[ichrt].deadslash = lulcdat[tstyr][ichrt].deadslash;
    telmnt[0].cohort[ichrt].leafslash = lulcdat[tstyr][ichrt].leafslash;
    telmnt[0].cohort[ichrt].leafconv = lulcdat[tstyr][ichrt].leafconv;
    telmnt[0].cohort[ichrt].rootslash = lulcdat[tstyr][ichrt].rootslash;
    telmnt[0].cohort[ichrt].rootconv = lulcdat[tstyr][ichrt].rootconv;
    telmnt[0].cohort[ichrt].stemslash = lulcdat[tstyr][ichrt].stemslash;
    telmnt[0].cohort[ichrt].stemconv = lulcdat[tstyr][ichrt].stemconv;
    telmnt[0].cohort[ichrt].standdead = lulcdat[tstyr][ichrt].standdead;
    telmnt[0].cohort[ichrt].deadconv = lulcdat[tstyr][ichrt].deadconv;
    telmnt[0].cohort[ichrt].deadslash = lulcdat[tstyr][ichrt].deadslash;
    telmnt[0].cohort[ichrt].sconvert = lulcdat[tstyr][ichrt].sconvert;
    telmnt[0].cohort[ichrt].prod10par = lulcdat[tstyr][ichrt].prod10par;
    telmnt[0].cohort[ichrt].prod100par = lulcdat[tstyr][ichrt].prod100par;



    telmnt[0].region = lulcdat[tstyr][ichrt].region;
    //if (telmnt[0].cohort[ichrt].vconvert >0.0) cout <<"col: "<<telmnt[0].col<<" row: "<<telmnt[0].row<<endl;
    //added by cgs2014. temporarily working code and need to delete later
    /*
    if (telmnt[0].cohort[ichrt].disturbflag ==0 && telmnt[0].cohort[ichrt].vrespar >0.0)
    {
    	telmnt[0].cohort[ichrt].disturbflag=2; //assume to harvest disturbance
    	telmnt[0].cohort[ichrt].disturbmonth=6;
    }
    */
    //if (telmnt[0].col == -98.0 && telmnt[0].row== 55.5 && pdyr < 170) cout <<" tstyr: "<<tstyr<<" pdyr: "<<pdyr<<" telm.year: "<<telmnt[0].year<<endl;
/*
    //added by cgs2014 to evaluate site level results
    if (telmnt[0].cohort[ichrt].currentveg ==4 && telmnt[0].col == -98.0 && telmnt[0].row== 55.5 && tstyr == 51) //can spruce //1850, 1930, 1964, 1981, 1989, 1998, 2003
    //if (telmnt[0].cohort[ichrt].currentveg ==3 && telmnt[0].col == -104 && telmnt[0].row== 64.25 && tstyr == 131) //can grass
    //if (telmnt[0].cohort[ichrt].currentveg ==4 && telmnt[0].col == -128.5 && telmnt[0].row== 60.25 && tstyr == 131) //ak spruce
    //if (telmnt[0].cohort[ichrt].currentveg ==3 && telmnt[0].col == -127.5 && telmnt[0].row== 63.25 && tstyr == 131) //ak grass
    //if (telmnt[0].cohort[ichrt].currentveg ==8 && telmnt[0].col == -83.75 && telmnt[0].row== 35.0 && tstyr == 131) //us forest
    //if (telmnt[0].cohort[ichrt].currentveg ==13 && telmnt[0].col == -99.5 && telmnt[0].row== 36.25 && tstyr == 131) //us grass
    //if (telmnt[0].cohort[ichrt].currentveg ==15 && telmnt[0].col == -104.75 && telmnt[0].row== 34 && tstyr == 131) //us shrub
    //if (telmnt[0].cohort[ichrt].currentveg ==15 && telmnt[0].col == -102 && telmnt[0].row== 22 && tstyr == 131) //mex shrub
    //if (telmnt[0].cohort[ichrt].currentveg ==16 && telmnt[0].col == -89.75 && telmnt[0].row== 16.5 && tstyr == 131) //mex forest
    //if (telmnt[0].cohort[ichrt].currentveg ==15 && telmnt[0].col == -107.25 && telmnt[0].row== 41.75 && tstyr == 166) //us shrub in Wyoming for model evaluation//year 2005: 39, 20, 6, 2
    {

    	//us & mexico grass/forest
    	int index=0;
    	int nonforest[20]= {1,2,3,7,12,13,15,21,22,23,24,26,28,30,35,39,40,49,50,51};
    	for (int i = 0;i<20;i++) {
    	if (telmnt[0].cohort[ichrt].currentveg == nonforest[i]) {
    		index=1;
    	}
    	}

    	telmnt[0].cohort[ichrt].disturbflag = 3;
    	telmnt[0].cohort[ichrt].disturbmonth = 8;
    	telmnt[0].cohort[ichrt].standage = 0;

    	 if (index ==1) //nonforest
    	 {
    		 telmnt[0].cohort[ichrt].leafmortpar = 1.0;
    		 telmnt[0].cohort[ichrt].stemmortpar = 1.0;
    		 telmnt[0].cohort[ichrt].rootmortpar =0.5;  //result in 50% of root death
    		 telmnt[0].cohort[ichrt].leafslash=0.2;
    		 telmnt[0].cohort[ichrt].leafconv=0.8;
    		 telmnt[0].cohort[ichrt].rootslash=1.0;
    		 telmnt[0].cohort[ichrt].rootconv =0.0;
    		 telmnt[0].cohort[ichrt].stemslash=0.2;
    		 telmnt[0].cohort[ichrt].stemconv=0.8;
    		 telmnt[0].cohort[ichrt].standdead=0.0;
    		 telmnt[0].cohort[ichrt].deadconv =1.0;
    		 telmnt[0].cohort[ichrt].deadslash=0.0;
    		 telmnt[0].cohort[ichrt].sconvert = 0.8;
    		 telmnt[0].cohort[ichrt].prod10par = 0.0;
    		 telmnt[0].cohort[ichrt].prod100par = 0.0;

    	 }
    	 else if (telmnt[0].cohort[ichrt].FRI <=35) { //if forest, only litter, leaf, and barks are converted
    		 telmnt[0].cohort[ichrt].leafmortpar = 0.15;
    		 telmnt[0].cohort[ichrt].stemmortpar = 0.15;
    		 telmnt[0].cohort[ichrt].rootmortpar =0.15;
    		 telmnt[0].cohort[ichrt].leafslash=0.0;
    		 telmnt[0].cohort[ichrt].leafconv=1.0;
    		 telmnt[0].cohort[ichrt].rootslash=1.0;
    		 telmnt[0].cohort[ichrt].rootconv =0.0;
    		 telmnt[0].cohort[ichrt].stemslash=0.2;
    		 telmnt[0].cohort[ichrt].stemconv=0.2;
    		 telmnt[0].cohort[ichrt].standdead=0.6;
    		 telmnt[0].cohort[ichrt].deadconv =0.2;
    		 telmnt[0].cohort[ichrt].deadslash=0.0;
    		 telmnt[0].cohort[ichrt].sconvert = 0.5;
    		 telmnt[0].cohort[ichrt].prod10par = 0.0;
    		 telmnt[0].cohort[ichrt].prod100par = 0.0;
    	 }
    	 else if (telmnt[0].cohort[ichrt].FRI >35 && telmnt[0].cohort[ichrt].FRI <=200)
    	 {
    		 telmnt[0].cohort[ichrt].leafmortpar = 0.5;
    		 telmnt[0].cohort[ichrt].stemmortpar = 0.5;
    		 telmnt[0].cohort[ichrt].rootmortpar =0.5;
    		 telmnt[0].cohort[ichrt].leafslash=0.2;
    		 telmnt[0].cohort[ichrt].leafconv=0.8;
    		 telmnt[0].cohort[ichrt].rootslash=1.0;
    		 telmnt[0].cohort[ichrt].rootconv =0.0;
    		 telmnt[0].cohort[ichrt].stemslash=0.3;
    		 telmnt[0].cohort[ichrt].stemconv=0.3;
    		 telmnt[0].cohort[ichrt].standdead=0.4;
    		 telmnt[0].cohort[ichrt].deadconv =0.5;
    		 telmnt[0].cohort[ichrt].deadslash=0.0;
    		 telmnt[0].cohort[ichrt].sconvert = 0.5;
    		 telmnt[0].cohort[ichrt].prod10par = 0.0;
    		 telmnt[0].cohort[ichrt].prod100par = 0.0;
    	 }
    	 else //if FRI >200
    	 {
    		 telmnt[0].cohort[ichrt].leafmortpar =0.75;
    		 telmnt[0].cohort[ichrt].stemmortpar =0.75;
    		 telmnt[0].cohort[ichrt].rootmortpar =0.75;
    		 telmnt[0].cohort[ichrt].leafslash=0.1;
    		 telmnt[0].cohort[ichrt].leafconv=0.9;
    		 telmnt[0].cohort[ichrt].rootslash=1.0;
    		 telmnt[0].cohort[ichrt].rootconv =0.0;
    		 telmnt[0].cohort[ichrt].stemslash=0.3;
    		 telmnt[0].cohort[ichrt].stemconv=0.4;
    		 telmnt[0].cohort[ichrt].standdead=0.3;
    		 telmnt[0].cohort[ichrt].deadconv =1.0;
    		 telmnt[0].cohort[ichrt].deadslash=0.0;
    		 telmnt[0].cohort[ichrt].sconvert = 0.5;
    		 telmnt[0].cohort[ichrt].prod10par = 0.0;
    		 telmnt[0].cohort[ichrt].prod100par = 0.0;
    	 }

    	if (index ==1) //nonforest
    	{
    	// ak and can grass
    	telmnt[0].cohort[ichrt].disturbflag = 3;
    	telmnt[0].cohort[ichrt].disturbmonth = 8;
        telmnt[0].cohort[ichrt].leafmortpar = 1.0;
        telmnt[0].cohort[ichrt].stemmortpar = 1.0;
        telmnt[0].cohort[ichrt].rootmortpar = 0.5;
        telmnt[0].cohort[ichrt].leafslash = 0.2;
        telmnt[0].cohort[ichrt].leafconv = 0.8;
        telmnt[0].cohort[ichrt].rootslash = 1.0;
        telmnt[0].cohort[ichrt].rootconv = 0.0;
        telmnt[0].cohort[ichrt].stemslash = 0.2;
        telmnt[0].cohort[ichrt].stemconv = 0.8;
        telmnt[0].cohort[ichrt].standdead = 0.0;
        telmnt[0].cohort[ichrt].deadconv = 1.0;
        telmnt[0].cohort[ichrt].deadslash = 0.0;
        telmnt[0].cohort[ichrt].sconvert = 0.8;
        telmnt[0].cohort[ichrt].standage = 0;
        telmnt[0].cohort[ichrt].prod10par=-99.0;
        telmnt[0].cohort[ichrt].prod100par=0.0;
    	}
    	else
    	{
		// ak and can forest
    	telmnt[0].cohort[ichrt].disturbflag = 3;
    	telmnt[0].cohort[ichrt].disturbmonth = 8;
        telmnt[0].cohort[ichrt].leafmortpar = 0.99;
        telmnt[0].cohort[ichrt].stemmortpar = 0.99;
        telmnt[0].cohort[ichrt].rootmortpar = 0.99;
        telmnt[0].cohort[ichrt].leafslash = 0.76;
        telmnt[0].cohort[ichrt].leafconv = 0.24;
        telmnt[0].cohort[ichrt].rootslash = 1.0;
        telmnt[0].cohort[ichrt].rootconv = 0.0;
        telmnt[0].cohort[ichrt].stemslash = 0.3;
        telmnt[0].cohort[ichrt].stemconv = 0.24;
        telmnt[0].cohort[ichrt].standdead = 0.46;
        telmnt[0].cohort[ichrt].deadconv = 0.6;
        telmnt[0].cohort[ichrt].deadslash = 0.4;
        telmnt[0].cohort[ichrt].sconvert = 0.05;
        telmnt[0].cohort[ichrt].standage = 0;
        telmnt[0].cohort[ichrt].prod10par=-99.0;
        telmnt[0].cohort[ichrt].prod100par=0.0;
    	}

    }
    else
    {
    	telmnt[0].cohort[ichrt].disturbflag = 0;
    	telmnt[0].cohort[ichrt].disturbmonth = 0;
        telmnt[0].cohort[ichrt].leafmortpar = 0;
        telmnt[0].cohort[ichrt].stemmortpar = 0;
        telmnt[0].cohort[ichrt].rootmortpar = 0;
        telmnt[0].cohort[ichrt].leafslash = 0;
        telmnt[0].cohort[ichrt].leafconv = 0;
        telmnt[0].cohort[ichrt].rootslash = 0;
        telmnt[0].cohort[ichrt].rootconv = 0;
        telmnt[0].cohort[ichrt].stemslash = 0;
        telmnt[0].cohort[ichrt].stemconv = 0;
        telmnt[0].cohort[ichrt].standdead = 0;
        telmnt[0].cohort[ichrt].deadconv = 0;
        telmnt[0].cohort[ichrt].deadslash = 0;
        telmnt[0].cohort[ichrt].sconvert = 0;
        telmnt[0].cohort[ichrt].prod100par=0.0;
        telmnt[0].cohort[ichrt].prod10par=0.0;


    }

    if (telmnt[0].cohort[ichrt].currentveg ==4 && telmnt[0].col == -98.0 && telmnt[0].row== 55.5 && tstyr > 51)
    //if (telmnt[0].cohort[ichrt].currentveg ==3 && telmnt[0].col == -104 && telmnt[0].row== 64.25 && tstyr > 131)
    //if (telmnt[0].cohort[ichrt].currentveg ==4 && telmnt[0].col == -128.5 && telmnt[0].row== 60.25 && tstyr > 131)
    //if (telmnt[0].cohort[ichrt].currentveg ==3 && telmnt[0].col == -127.5 && telmnt[0].row== 63.25 && tstyr > 131)
    //if (telmnt[0].cohort[ichrt].currentveg ==8 && telmnt[0].col == -83.75 && telmnt[0].row== 35.0 && tstyr > 131)
    //if (telmnt[0].cohort[ichrt].currentveg ==13 && telmnt[0].col == -99.5 && telmnt[0].row== 36.25 && tstyr > 131)
    //if (telmnt[0].cohort[ichrt].currentveg ==15 && telmnt[0].col == -104.75 && telmnt[0].row== 34  && tstyr > 131)
    //if (telmnt[0].cohort[ichrt].currentveg ==15 && telmnt[0].col == -102 && telmnt[0].row== 22  && tstyr > 131)
    //if (telmnt[0].cohort[ichrt].currentveg ==16 && telmnt[0].col == -89.75 && telmnt[0].row== 16.5 && tstyr > 131)
    //if (telmnt[0].cohort[ichrt].currentveg ==15 && telmnt[0].col == -107.25 && telmnt[0].row== 41.75 && tstyr > 166)
    {
    	telmnt[0].cohort[ichrt].forestage=tstyr-51;
    }

*/
    //if ( telmnt[0].cohort[ichrt].sconvert>0.01)
    	//cout <<"col: "<<telmnt[0].col<<" row: "<<telmnt[0].row<<" year: "<<dyr<<"ichrt: "<<ichrt<<" sconvert: "<<telmnt[0].cohort[ichrt].sconvert<<" disturbflag: "<<telmnt[0].cohort[ichrt].disturbflag<<endl;

    telmnt[0].cohort[ichrt].cmnt = telmnt[0].lcluc.getCommunityType( lulcdat[tstyr][ichrt].currentveg );//modified by cgs to replace ".subtype" with ".currentveg"
    //telmnt[0].cohort[ichrt].cmnt = telmnt[0].lcluc.getCommunityType( lulcdat[tstyr][ichrt].subtype );
    telmnt[0].cohort[ichrt].agcmnt = telmnt[0].cohort[ichrt].cmnt; //closed by cgs2014

    //if( pdyr > 0 && lulcdat[tstyr][ichrt].agstate > 0 )//revised by cgs2014 to simulate crop during equilibrium run
    if (lulcdat[tstyr][ichrt].agstate > 0 ) //modified by cgs2014
    {
       telmnt[0].cohort[ichrt].agcmnt = telmnt[0].lcluc.getCommunityType( lulcdat[tstyr][ichrt].currentveg );
    }

    telmnt[0].cohort[ichrt].wfpsoff = telmnt[0].lcluc.getWFPSOFF( lulcdat[tstyr][ichrt].subtype );
    //cout <<"agstate :"<< telmnt[0].cohort[ichrt].agstate<<" vegcmt: "<<telmnt[0].cohort[ichrt].cmnt<<" agcmnt: "<<telmnt[0].cohort[ichrt].agcmnt<<endl;

  }

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void updateTTEMGridCell( const int& pdyr)
{
  int dm;

  int ichrt;

  int tchrt;

/* *************************************************************
            INITIALIZE TEM STATE FOR NEW COHORTS
************************************************************* */
  
  if( telmnt[0].maxcohorts > telmnt[0].prvmxcohrts )
  {
    for ( ichrt = telmnt[0].prvmxcohrts; 
          ichrt < telmnt[0].maxcohorts; 
          ++ichrt )
    {     
      tchrt = telmnt[0].cohort[ichrt].srcCohort - 1;
      
      telmnt[0].setCohortTEMState( telmnt[0].cohort[tchrt],
                                   telmnt[0].cohort[ichrt] );
    }
  }

/* *************************************************************
                     UPDATE TEM FOR GRID CELL
************************************************************* */

    
/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */
  int tempcohortnum=telmnt[0].maxcohorts;
  if (telmnt[0].maxcohorts==0) tempcohortnum=1;
  for( ichrt = 0; ichrt < tempcohortnum; ++ichrt)
  {
    // Get vegetation community type of cohort



	//if( pdyr > 0 && lulcdat[pdyr][ichrt].agstate > 0 )
	//{
	//telmnt[0].cohort[ichrt].agcmnt = telmnt[0].lcluc.getCommunityType( lulcdat[pdyr][ichrt].currentveg );
	//}

	telmnt[0].tem.veg.cmnt = telmnt[0].cohort[ichrt].cmnt;
    
    // Determine soil characteristics for cohort
    
    telmnt[0].tem.soil.xtext( telmnt[0].tem.veg.cmnt,
                              telmnt[0].tem.soil.getPCTSILT(),
                              telmnt[0].tem.soil.getPCTCLAY() );
    //recalculate standage
    //if (pdyr >1) telmnt[0].cohort[ichrt].standage+=1;

    double tempxx=0.0;
    double tempyy=0.0;
    double tempzz=0.0;
    double tempmm=0.0;
    double tempnn=0.0;
    double tempoo=0.0;
    double temppp=0.0;
    for( dm = 0; dm < CYCLE; ++dm )
    {     
      // Run TEM

      telmnt[0].updateTEMmonth( flog1,
                                equil, 
                                totsptime, 
                                pdyr, 
                                dm,
                                ichrt );
      tempxx += telmnt[0].output[telmnt[0].tem.I_NPP][dm];
      tempyy += telmnt[0].output[telmnt[0].tem.I_NTCB][dm];
      tempzz += telmnt[0].output[telmnt[0].tem.I_GPP][dm];
      tempmm += telmnt[0].output[telmnt[0].tem.I_NEP][dm];
      tempnn += telmnt[0].output[telmnt[0].tem.I_NCE][dm];
      //tempoo += telmnt[0].output[telmnt[0].tem.I_RH][dm];
      tempoo += telmnt[0].output[telmnt[0].tem.I_LTRC][dm];
      temppp += telmnt[0].output[telmnt[0].tem.I_CNVRTC][dm];
      //if (telmnt[0].row==39.25&& telmnt[0].col==-121.75&& (pdyr==144||pdyr==145) && ichrt ==5) cout <<"pdyr: "<<pdyr<<" mon: "<<dm<<" GPP: "<<telmnt[0].output[telmnt[0].tem.I_GPP][dm]<<" NPP: "<<telmnt[0].output[telmnt[0].tem.I_NPP][dm]<<" forestage: "<<telmnt[0].output[telmnt[0].tem.I_AGE][dm] <<" MDC: "<<telmnt[0].output[telmnt[0].tem.I_MDC][dm]
       //<<" AGL: "<<telmnt[0].output[telmnt[0].tem.I_AGL][dm]<<" AGR: "<<telmnt[0].output[telmnt[0].tem.I_AGR][dm]<<" VEGC: "<<telmnt[0].output[telmnt[0].tem.I_VEGC][dm]<<" SOLC: "<<telmnt[0].output[telmnt[0].tem.I_SOLC][dm]<<" SOC: "<<telmnt[0].output[telmnt[0].tem.I_SOC][dm]<<" DEADWOOD: "<<telmnt[0].output[telmnt[0].tem.I_DEADWOODC][dm]
       //<<" CWD: "<<telmnt[0].output[telmnt[0].tem.I_CWD][dm]<<" convertc: "<<telmnt[0].output[telmnt[0].tem.I_CNVRTC][dm]<<endl;
      //if (telmnt[0].col == -80.75 && telmnt[0].row== 51.25 && telmnt[0].tem.veg.cmnt==4) cout << "year: " <<pdyr<<" mon: "<<dm<<" cmnt: "<<telmnt[0].tem.veg.cmnt<<" vegc: "<<telmnt[0].tem.getY(0)<<" GPP: "<< tempnn <<" NPP: "<<tempoo<<" RH: "<<temppp<<" ntcb: "<<tempxx<<endl;

      int mon1 = dm-1;
      if (mon1<0) mon1=0;
      double calculated_storage1=telmnt[0].output[telmnt[0].tem.I_VEGC][dm]+telmnt[0].output[telmnt[0].tem.I_SOLC][dm]+telmnt[0].output[telmnt[0].tem.I_DEADWOODC][dm];
      double calculated_storage2=(telmnt[0].output[telmnt[0].tem.I_VEGC][mon1]+telmnt[0].output[telmnt[0].tem.I_SOLC][mon1]+telmnt[0].output[telmnt[0].tem.I_DEADWOODC][mon1]);
      double calculated_storage =calculated_storage1-calculated_storage2;
      double nce1=telmnt[0].output[telmnt[0].tem.I_NTCB][dm]+0.5;
      double nce2=telmnt[0].output[telmnt[0].tem.I_NTCB][dm]-0.5;
      //if (dm >0 && telmnt[0].output[telmnt[0].tem.I_RH][dm] >=0.0 && (calculated_storage > nce1 ||calculated_storage<nce2)) cout <<" lat: "<<telmnt[0].row<<" lon: "<<telmnt[0].col<<" year: "<<pdyr<<" mon: "<<dm<<" disturbflag: "<<telmnt[0].tem.disturbflag<<" disturbmon: "<<telmnt[0].tem.disturbmonth<<" currentveg: "<< telmnt[0].cohort[ichrt].currentveg<<" cmnt: "<<telmnt[0].tem.veg.cmnt <<" ichrt: "<<ichrt<<" Cstock2: "<< calculated_storage1<<" Cstock1: "<<calculated_storage2<<" ntcb2: "<<telmnt[0].output[telmnt[0].tem.I_NTCB][dm]
         //<< " diffvegc: "<<telmnt[0].output[telmnt[0].tem.I_VEGC][dm]-telmnt[0].output[telmnt[0].tem.I_VEGC][mon1] <<" diffsoilc: "<< telmnt[0].output[telmnt[0].tem.I_SOLC][dm]-telmnt[0].output[telmnt[0].tem.I_SOLC][mon1]<<" diffwoodc: "<< telmnt[0].output[telmnt[0].tem.I_DEADWOODC][dm]-telmnt[0].output[telmnt[0].tem.I_DEADWOODC][mon1]<<" deadwoodc: "<<telmnt[0].output[telmnt[0].tem.I_DEADWOODC][dm]<<" forestage: "<< telmnt[0].cohort[ichrt].forestage<<endl;
      double monvegc = telmnt[0].output[telmnt[0].tem.I_VEGC][dm];
      double monsoc = telmnt[0].output[telmnt[0].tem.I_SOLC][dm];
      double monwood = telmnt[0].output[telmnt[0].tem.I_DEADWOODC][dm];
      double monntcb = telmnt[0].output[telmnt[0].tem.I_NTCB][dm];
      if ( telmnt[0].row==39.25 && telmnt[0].col == -109.75 && pdyr <=190) cout <<" lat: "<<telmnt[0].row<<" lon: "<<telmnt[0].col<<" year: "<<pdyr<<" mon: "<<dm<<" disturbflag: "<<telmnt[0].tem.disturbflag<<" disturbmon: "<<telmnt[0].tem.disturbmonth<<" currentveg: "<< telmnt[0].cohort[ichrt].currentveg<<" cmnt: "<<telmnt[0].tem.veg.cmnt <<" ichrt: "<<ichrt<<" vegc: "<<monvegc<<" soc: "<<monsoc<<" woodc: "<<monwood<<" ntcb: "<<monntcb<<" carbonpool: "<<monvegc+monsoc+monwood<<" subarea: "<<telmnt[0].cohort[ichrt].chrtarea<<endl;


    } // end of CYCLE loop



    double templitter1 = telmnt[0].output[telmnt[0].tem.I_AGL][11] + telmnt[0].output[telmnt[0].tem.I_AGR][11];
    double templitter2 = telmnt[0].output[telmnt[0].tem.I_BGL][11] + telmnt[0].output[telmnt[0].tem.I_BGR][11];
    double totvegc2= telmnt[0].output[telmnt[0].tem.I_FOLC][11]+telmnt[0].output[telmnt[0].tem.I_STEMC][11]+telmnt[0].output[telmnt[0].tem.I_CROOTC][11]+telmnt[0].output[telmnt[0].tem.I_FROOTC][11];
    double templitterconsum=telmnt[0].output[telmnt[0].tem.I_FFLC][7]+telmnt[0].output[telmnt[0].tem.I_SWFC][7]+telmnt[0].output[telmnt[0].tem.I_DWFC][7];
    //if (telmnt[0].col == -80.75 && telmnt[0].row== 51.25 && telmnt[0].tem.veg.cmnt==4)
    //if (telmnt[0].tem.veg.cmnt==4 && telmnt[0].col == -98.0 && telmnt[0].row== 55.50 && ichrt==0) cout <<" year: "<<pdyr<<" dm: "<<dm<<" cohort: "<<ichrt<<" currentveg: "<< telmnt[0].cohort[ichrt].currentveg<<" cmnt: "<<telmnt[0].tem.veg.cmnt<<" vegc: "<< telmnt[0].output[telmnt[0].tem.I_VEGC][11]<<" soc: "<< telmnt[0].output[telmnt[0].tem.I_SOLC][11]<<" NPP: "<<tempxx<<" GPP: "<<tempzz<<" rh: "<<tempoo<<" ntcb: "<<tempyy<<" nep: "<<tempmm<<" nce: "<<tempnn<<" deadwoodc: "<<telmnt[0].output[telmnt[0].tem.I_DEADWOODC][11]<<" CWD: "<<telmnt[0].output[telmnt[0].tem.I_CWD][11]<<" abovelit: "<<templitter1<<" belowlit: "<<templitter2<<" standage: "<<telmnt[0].cohort[ichrt].standage<<endl;
    //output for sensitivity analysis
    //if (telmnt[0].tem.veg.cmnt==4 && telmnt[0].col == -98.0 && telmnt[0].row== 55.50 && ichrt==0) cout <<" year: "<<pdyr<<" currentveg: "<< telmnt[0].cohort[ichrt].currentveg<<" cmnt: "<<telmnt[0].tem.veg.cmnt <<" NPP: "<<tempxx<<" GPP: "<<tempzz<<" litfalc: "<<tempoo<<" ntcb: "<<tempyy<<" ffdc: "<<telmnt[0].output[telmnt[0].tem.I_FFDC][7]<<" burntlit: "<<templitterconsum<<" convertc: "<<temppp<<" vegc: "<<telmnt[0].output[telmnt[0].tem.I_VEGC][11]<<" vegc2: "<<totvegc2<<" soc: "<<telmnt[0].output[telmnt[0].tem.I_SOLC][11]<<" abovlitterc: "<< templitter1<< " belowlitterc: "<< templitter2<<" cwd: "<< telmnt[0].output[telmnt[0].tem.I_CWD][11] <<" standdead: "<< telmnt[0].output[telmnt[0].tem.I_DEADWOODC][11]<<" forestage: "<< telmnt[0].cohort[ichrt].forestage<<endl;

    //below code added to verify if NTCB = carbon storage change between two months

    //double calculated_storage1=telmnt[0].output[telmnt[0].tem.I_VEGC][7]+telmnt[0].output[telmnt[0].tem.I_SOLC][7]+telmnt[0].output[telmnt[0].tem.I_DEADWOODC][7];
    //double calculated_storage2=(telmnt[0].output[telmnt[0].tem.I_VEGC][6]+telmnt[0].output[telmnt[0].tem.I_SOLC][6]+telmnt[0].output[telmnt[0].tem.I_DEADWOODC][6]);
    //double calculated_storage =calculated_storage1-calculated_storage2;
    //double nce1=telmnt[0].output[telmnt[0].tem.I_NTCB][7]+0.5;
    //double nce2=telmnt[0].output[telmnt[0].tem.I_NTCB][7]-0.5;
    //if (telmnt[0].output[telmnt[0].tem.I_RH][7] >=0.0 && (calculated_storage > nce1 ||calculated_storage<nce2)&&telmnt[0].row==39.25&& telmnt[0].col==-121.75 && ichrt ==5) cout <<" year: "<<pdyr<<" lat: "<<telmnt[0].row<<" lon: "<<telmnt[0].col<<" disturbmon: "<<telmnt[0].tem.disturbmonth<<" currentveg: "<< telmnt[0].cohort[ichrt].currentveg<<" cmnt: "<<telmnt[0].tem.veg.cmnt <<" ichrt: "<<ichrt<<" Cstock7: "<< calculated_storage1<<" Cstock6: "<<calculated_storage2<<" ntcb7: "<<telmnt[0].output[telmnt[0].tem.I_NTCB][7]<<" age6: "<<telmnt[0].output[telmnt[0].tem.I_AGE][6]<<" age7: "<<telmnt[0].output[telmnt[0].tem.I_AGE][7]
     //<< " diffvegc: "<<telmnt[0].output[telmnt[0].tem.I_VEGC][7]-telmnt[0].output[telmnt[0].tem.I_VEGC][6] <<" npp7: "<< telmnt[0].output[telmnt[0].tem.I_NPP][7]<<" NPP: "<<tempxx<<" GPP: "<<tempzz<<" litfalc: "<<tempoo<<" ntcb: "<<tempyy<<" ffdc: "<<telmnt[0].output[telmnt[0].tem.I_FFDC][7]<<" burntlit: "<<templitterconsum<<" convertc: "<<temppp<<" vegc: "<<telmnt[0].output[telmnt[0].tem.I_VEGC][11]<<" vegc2: "<<totvegc2<<" soc: "<<telmnt[0].output[telmnt[0].tem.I_SOLC][11]<<" abovlitterc: "<< templitter1<< " belowlitterc: "<< templitter2<<" cwd: "<< telmnt[0].output[telmnt[0].tem.I_CWD][11] <<" standdead: "<< telmnt[0].output[telmnt[0].tem.I_DEADWOODC][11]<<" forestage: "<< telmnt[0].cohort[ichrt].forestage<<endl;
    //if (telmnt[0].row==39.25&& telmnt[0].col==-121.75&&pdyr>145) exit(-1);


    //if (telmnt[0].col == -127.75 && telmnt[0].row ==55.25) cout <<" year: "<<pdyr<<" dm: "<<dm<<" cohort: "<<ichrt<<" currentveg: "<< telmnt[0].cohort[ichrt].currentveg<<" cmnt: "<<telmnt[0].tem.veg.cmnt<<" vegc: "<< telmnt[0].output[telmnt[0].tem.I_VEGC][11]<<" soc: "<< telmnt[0].output[telmnt[0].tem.I_SOLC][11]<<" NPP: "<<tempxx<<" GPP: "<<tempzz<<" rh: "<<tempoo<<" ntcb: "<<tempyy<<" nep: "<<tempmm<<" nce: "<<tempnn<<" deadwoodc: "<<telmnt[0].output[telmnt[0].tem.I_DEADWOODC][11]<<" CWD: "<<telmnt[0].output[telmnt[0].tem.I_CWD][11]<<" abovelit: "<<templitter1<<" belowlit: "<<templitter2<<" standage: "<<telmnt[0].cohort[ichrt].standage<<endl;

    //telmnt[0].tem.disturbflag >0 &&
    //if (telmnt[0].tem.veg.cmnt==8 && telmnt[0].col == -107.25 && telmnt[0].row== 41.75) cout <<" year: "<<pdyr<<" currentveg: "<< telmnt[0].cohort[ichrt].currentveg<<" cmnt: "<<telmnt[0].tem.veg.cmnt <<" NPP: "<<tempxx<<" GPP: "<<tempzz<<" ntcb: "<<tempyy<<" ffdc: "<<telmnt[0].output[telmnt[0].tem.I_FFDC][7]<<" burntlit: "<<templitterconsum<<" convertc: "<<temppp<<" vegc: "<<telmnt[0].output[telmnt[0].tem.I_VEGC][11]<<" soc: "<<telmnt[0].output[telmnt[0].tem.I_SOLC][11]<<" abovlitterc: "<< templitter1<< " belowlitterc: "<< templitter2<<" cwd: "<< telmnt[0].output[telmnt[0].tem.I_CWD][11] <<" standdead: "<< telmnt[0].output[telmnt[0].tem.I_DEADWOODC][11]<<" standage: "<< telmnt[0].cohort[ichrt].forestage<<endl;


    /*//below code is used for uncertainty analysis
    //if (telmnt[0].tem.veg.cmnt==4 && telmnt[0].col == -98.0 && telmnt[0].row== 55.50 && ichrt==0 && pdyr >=120)
    //if (telmnt[0].tem.veg.cmnt==3 && telmnt[0].col == -104 && telmnt[0].row== 64.25 && ichrt==0 && pdyr >=120)
    //if (telmnt[0].tem.veg.cmnt==4 && telmnt[0].col == -128.5 && telmnt[0].row== 60.25 && ichrt==0 && pdyr >=120)
    //if (telmnt[0].tem.veg.cmnt==3 && telmnt[0].col == -127.5 && telmnt[0].row== 63.25 && ichrt==0 && pdyr >=120)
    //if (telmnt[0].tem.veg.cmnt==6 && telmnt[0].col == -83.75 && telmnt[0].row== 35.0 && ichrt==0 && pdyr >=120)
    //if (telmnt[0].tem.veg.cmnt==7 && telmnt[0].col == -99.5 && telmnt[0].row== 36.25 && ichrt == 1 && pdyr >=120)
    //if (telmnt[0].tem.veg.cmnt==8 && telmnt[0].col == -104.75 && telmnt[0].row== 34 && ichrt == 1 && pdyr >=120)
    //if (telmnt[0].tem.veg.cmnt==8 && telmnt[0].col == -102 && telmnt[0].row== 22 && pdyr >=120) cout <<"ichrt: "<<ichrt<<endl;
    //if (telmnt[0].tem.veg.cmnt==8 && telmnt[0].col == -102 && telmnt[0].row== 22 && ichrt == 2 && pdyr >=120)
    if (telmnt[0].tem.veg.cmnt==9 && telmnt[0].col == -89.75 && telmnt[0].row== 16.5 && ichrt == 0 && pdyr >=120)
    {

    	outfile1<<pdyr <<" " <<tempzz<<endl; //gpp
    	outfile2<<pdyr <<" " <<tempxx<<endl; //npp
    	outfile3<<pdyr <<" " <<tempyy<<endl; //ntcb
    	outfile4<<pdyr <<" " <<temppp<<endl; //convertc
    	//cout <<"pdyr: "<<pdyr<<" convertc: "<<temppp<<endl;
    	//outfile5<<pdyr <<" " <<telmnt[0].output[telmnt[0].tem.I_VEGC][11]<<endl; //vegc
    	//outfile6<<pdyr <<" " <<telmnt[0].output[telmnt[0].tem.I_SOLC][11]<<endl; //soc
    	//outfile7<<pdyr <<" " <<telmnt[0].output[telmnt[0].tem.I_DEADWOODC][11]<<endl; //standing deadwood
    	outfile5<<pdyr <<" " <<telmnt[0].output[telmnt[0].tem.I_FFDC][7]<<endl;
    	outfile6<<pdyr <<" " <<telmnt[0].output[telmnt[0].tem.I_FFLC][7]+telmnt[0].output[telmnt[0].tem.I_SWFC][7]+telmnt[0].output[telmnt[0].tem.I_DWFC][7]<<endl;
    	outfile7<<pdyr <<" " <<telmnt[0].output[telmnt[0].tem.I_VCNVRTC][7]<<endl;
    	outfile8<<pdyr <<" " <<templitter1<<endl; //abovelit
    	outfile9<<pdyr <<" " <<templitter2<<endl; //belowlit

    }
   */
    if( 2 == ostateflag && telmnt[0].tem.totyr == ostateyear )
    {
      telmnt[0].writeCohortState( ofstate, ichrt );
    }
    //if (telmnt[0].output[19][6] > 0.0)cout <<"outputgpp0.5: "<<telmnt[0].output[19][6]<< " getgpp0: "<<telmnt[0].tem.getY(19)<<" tem.i_gpp: "<<telmnt[0].tem.I_GPP<<endl;
    //cout <<"clmyr: "<<clmstartyr<<"totyr: "<<telmnt[0].tem.totyr<<" spinoutyr: "<<spinoutyrs<<endl;
    if ( (1 == spinoutfg && telmnt[0].tem.totyr < clmstartyr)
         || (2 == spinoutfg 
         && telmnt[0].tem.totyr >= (clmstartyr-spinoutyrs))
         || (telmnt[0].tem.totyr >= telmnt[0].tem.outputstartyr && telmnt[0].tem.totyr <= telmnt[0].tem.outputendyr
         && (telmnt[0].wrtyr % telmnt[0].tem.diffyr==0)) )
    {
     //if (telmnt[0].output[19][6]>0.0) cout <<"pdyr: "<<pdyr+1800<<"outputgpp1: "<<telmnt[0].output[19][6]<<endl;
      // Output TEM transient results for specified years to files
         telmnt[0].temwritepred( ftempred,
                              tempredmap, 
                              pdyr, 
                              ichrt,
                              telmnt[0].ntempred );
    }
  } // End of cohort loop
  
};

