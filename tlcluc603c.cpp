/* **************************************************************
*****************************************************************
TLCLUC603C.CPP - determines both potential and "actual" land
                 cover, and identifies land uses within a grid
                 cell

20031019 - DWK created by modifying tlcluc431b.cpp
20030624 - DWK changed TEMlcluc43:: to TEMlcluc51::
20030624 - DWK changed inheritance of Biome43() to Biome51()
20040714 - DWK changed TEMlcluc51:: to TEMlcluc60::
20040714 - DWK changed inheritance of Biome43() to Biome60()
20040716 - DWK added getCohorts(), getVegtype() and initCohorts()
20040716 - DWK added cohorts.col, cohorts.row, cohorts.total, 
           and vegtype.temveg to TEMlcluc60()
20040716 - DWK deleted potveg.col, potveg.row and potveg.temveg 
           from TEMlcluc60()
20040716 - DWK deleted getPotentialVeg() and  initPotentialVeg()
20040828 - DWK renamed public getCohorts() to be 
           getNumberOfCohorts()
20040829 - DWK renamed initCohorts() as initMaxCohorts()
20040829 - DWK renamed initLandUse() as initCohorts()
20040829 - DWK renamed getLandUse() as getCohort()
20070514 - DWK changed include tlcluc602.h to tlcluc603a.h
20070830 - DWK changed include from tlcluc603a.h to tlcluc603c.h          
		    
*****************************************************************
************************************************************** */

#include<cstdio>

  using std::fscanf;
  using std::FILE;

#include<iostream>

  using std::cin;
  using std::cout;
  using std::ios;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setprecision;

#include<string>

  using std::string;


#include "tlcluc603c.h"

/* *************************************************************
************************************************************* */

TEMlcluc60::TEMlcluc60() : Biome60()
{

  cohorts.col = -999.9;
  cohorts.row = -999.9;
  cohorts.total = -99;
  cohorts.natchrts = -99;
  
  lulc.year = -999;
  lulc.agstate = -99;
  lulc.agprevstate = -99;
  lulc.standage = -99;
  lulc.forestage=-99;

  maxtype = -999;

  agcmnt = -999;
  cmnt = -999;

  lastyr = -999;


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMlcluc60::getCohort( FILE* flulc )
{

  int gisend = lulc.getdel( flulc );
  
  if( -1 == gisend ) { return gisend; }

  return gisend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMlcluc60::getNumberOfCohorts( FILE* fnchrts )
{

  int gisend = cohorts.getdel( fnchrts );
  
  if( -1 == gisend ) { return gisend; }

  return gisend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc60::initCohorts( ofstream& rflog1 )
{

#ifdef PMODE
  *fgo >> ilulcfname;
  *fgo >> ilulcend;
#else
  cout << "Please enter the first part of the file name ";
  cout << "containing the cohort data: " << endl;
  cout << "               (e.g., COHORT) " << endl;

  cin >> ilulcfname;

  cout << "Please enter the file extension (include the '.'): ";

  cin >> ilulcend;
#endif

  rflog1 << "Please enter the first part of the file name ";
  rflog1 << "containing the cohort data: " << endl;
  rflog1 << "               (e.g., COHORT) " << endl;
  rflog1 << ilulcfname << endl << endl;
  rflog1 << "Please enter the file extension (include the '.'): ";
  rflog1 << endl;
  rflog1 << ilulcend << endl << endl;

  if( 1 == tlulcflag )
  {

#ifdef PMODE
    *fgo >> lastyr;

#else

	cout << "Please enter the last year for which you have ";
    cout << "cohort data: " << endl;

    cin >> lastyr;
#endif

    rflog1 << "Please enter the last year for which you have ";
    rflog1 << "cohort data: " << endl;
    rflog1 << lastyr << endl << endl;
  }


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc60::initMaxCohorts( ofstream& rflog1 )
{
#ifdef PMODE
  *fgo >> imxcohrtfname;
  *fgo >> imxcohrtend;

#else
 
  cout << "Please enter the first part of the file name ";
  cout << "containing the maximun number of cohort data";
  cout << endl;
  cout << "        (e.g., MXCOHRTS) " << endl;

  cin >> imxcohrtfname;
  
  cout << "Please enter the file extension (include the '.'): ";

  cin >> imxcohrtend;
#endif

  rflog1 << "Please enter the first part of the file name ";
  rflog1 << "containing the maximun number of cohort data";
  rflog1 << endl;
  rflog1 << "        (e.g., MXCOHRTS) " << endl;
  rflog1 << imxcohrtfname << endl << endl;
  rflog1 << "Please enter the file extension (include the '.'): ";
  rflog1 << imxcohrtend << endl << endl;

  return;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TEMlcluc60::setLCLUCFlags( ofstream& rflog1, 
                                const int& requil )
{

  tlulcflag = 0;
  
  if( 0 == requil )
  {

#ifdef PMODE
    *fgo >> tlulcflag;

#else
    cout << "Do you have transient land use data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    
    cin >> tlulcflag;
#endif

    rflog1 << "Do you have transient land use data?:" << endl;
    rflog1 << "Enter 0 for No:" << endl;
    rflog1 << "Enter 1 for Yes: ";
    rflog1 << tlulcflag << endl << endl;
  }

};

