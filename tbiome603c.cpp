/* *************************************************************
****************************************************************
TBIOME603C.CPP - object describing general characteristics of
               vegetation mosaic used in the Terrestrial
	       Ecosystem Model (TEM)

Modifications:

20030119 - DWK created by modifying tbiome431a.cpp
20031019 - DWK changed inheritance from ProcessXML to ProcessXML51
20031019 - DWK changed Biome43:: to Biome51::
20031019 - DWK changed char ecd[MAXFGNAME] to const string& ecd
           in getvtype()
20031021 - DWK added  ostringstream tempString to getvtype()
20031021 - DWK added public function getVegSubarea()
20040714 - DWK changed Biome51:: to Biome60::
20040714 - DWK changed inheritance of ProcessXML51() to
           inheritance of ProcessXML60()
20040716 - DWK added getCommunityType()
20051117 - DWK added include tbiom602.h and standard includes
20070830 - DWK added wfpsoff[] to getvtype()
20070830 - DWK changed include from tbiome602.h to tbiome603c.h
           
****************************************************************
************************************************************* */

#include<iostream>

  using std::cin;
  using std::cout;
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
  
#include<string>
  
  using std::string;

#include<sstream>

  using std::ostringstream;


#include "tbiome603c.h"


/* *************************************************************
************************************************************* */

Biome60::Biome60( void ) : ProcessXML60()
{
  temveg = -99;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome60::getCommunityType( const int& tveg )
{
  int mez;
  int communtype;
  
  mez = tveg - 1;
  
  if( mez < 0 || mez >= NUMVEG )
  {
    communtype = 1;
  }
  else { communtype = subtype[mez][0]; }

  return communtype;	
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome60::getVegMosaic( const int& tveg )
{
  int mez;
  int maxtype;

  mez = tveg - 1;
  
  if( mez < 0 || mez >= NUMVEG )
  {
    maxtype = 1;
  }
  else { maxtype = numtype[mez]; }

  return maxtype;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Biome60::getVegSubarea( const int& tveg,
                               const int& dtype,
                               const int& carea )
{
  int mez;
  double sarea;

  mez = tveg - 1;
  
  if( mez < 0 || mez >= NUMVEG )
  {
    sarea = (double) carea;
  }
  else { sarea = (double) carea * pcttype[mez][dtype] * 0.01; }

  return sarea;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Biome60::getVegSubtype( const int& tveg, const int& dtype )
{
  int mez;
  int vegtype;

  mez = tveg - 1;
  
  if( mez < 0 || mez >= NUMVEG )
  {
    vegtype = 1;
  }
  else { vegtype = subtype[mez][dtype]; }

  return vegtype;
};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Biome60::getvtype( ofstream& rflog1 )
{

  string ecd;

#ifdef PMODE
  *fgo >> ecd;
#else
  cout << endl;
  cout << "Enter name of the file prescribing vegetation mosaics (.ECD):";
  cout << endl;
  
  cin >> ecd;
#endif

  rflog1 << endl;
  rflog1 << "Enter name of the file prescribing vegetation mosaics (.ECD):";
  rflog1 << endl << ecd << endl << endl;

  getvtype( ecd );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Biome60::getvtype( const string& ecd )
{
  ifstream infile;
  int dv;
  int dtype;
  int vegtype;
  int ez;

  ostringstream tempString;

  infile.open( ecd.c_str(), ios::in );

  if ( !infile )
  {
    cerr << endl << "Cannot open " << ecd;
    cerr << " for community ECD input" << endl;

    exit( -1 );
  }

  getXMLrootNode( infile, "communityECD" );


  for ( dv = 0; dv < NUMVEG; ++dv )
  {

    vegtype = getXMLtemvegNode( infile, "communityECD" );

    if ( vegtype > NUMVEG )
    {
      cerr << endl << "TEMVEG type " << vegtype << endl;
      cerr << " cannot be greater than " << NUMVEG;
      cerr << " in communityECD" << endl;
      
      exit( -1 );
    }

    ez = vegtype - 1;

    numtype[ez] = getXMLtvegArrayInt( infile,
                                      "communityECD",
                                      "numtype",
                                      vegtype );

    for ( dtype = 0; dtype < NUMMSAC; ++dtype )
    {
      tempString.str( "" );
      tempString << "subtype" << (dtype+1);
      subtype[ez][dtype] = getXMLtvegArrayInt( infile,
                                               "communityECD",
                                               tempString.str(),
                                               vegtype );

      tempString.str( "" );
      tempString << "pcttype" << (dtype+1);
      pcttype[ez][dtype] = getXMLtvegArrayDouble( infile,
                                                  "communityECD",
                                                  tempString.str(),
                                                  vegtype );
    }

    wfpsoff[ez] = getXMLtvegArrayInt( infile,
                                      "communityECD",
                                      "wfpsoff",
                                      vegtype );

    endXMLtvegNode( infile );
  }


  if ( dv < NUMVEG )
  {
    cerr << endl << " Parameters found for only " << dv;
    cerr << " community types out of a maximum of ";
    cerr << (NUMVEG-1) << " types in communityECD" << endl;

    exit( -1 );
  }

  infile.close();

};


/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Biome60::getWFPSOFF( const int& tveg )
{
  int mez;
  
  mez = tveg - 1;
  
  if( mez < 0 || mez >= NUMVEG ) { return ZERO; }
  else { return wfpsoff[mez]; }

};
