/* *************************************************************
TPROCESSMXL602.CPP - Functions to process XML files for TEM

Programmer: David Kicklighter
Creation Date: 17 June 2003

Modifications:

20031019 - DWK changed ProcessXML:: to ProcessXML51::
20031019 - DWK changed char rootnode[80] to string rootnode
           and char varnode[80] to string varnode in functions
20031019 - DWK modified functions to use standard string
           processing algorithms
20031020 - DWK deleted getvalue() and readline()
20040707 - DWK changed ProcessXML51:: to ProcessXML60::
20051117 - DWK added include tprocessXML602.h

************************************************************* */

#include<iostream>

  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std::atoi;
  using std::atol;
  
#include<string>
  
  using std::string;

#include "tprocessXML602.h"


ProcessXML60::ProcessXML60()
{

};

/* *************************************************************
************************************************************* */

void ProcessXML60::endXMLcommunityNode( ifstream& infile )

{
  string line;

  while ( line.find( "</community>" ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void ProcessXML60::endXMLtvegNode( ifstream& infile )

{
  string line;

  while ( line.find( "</temveg>" ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML60::getXMLcommunityNode( ifstream& infile,
                                       const string& rootnode )

{
  string line;
  string value;

  string temp;

  int startString;

  while ( line.find( ">" ) == string::npos
          && !infile.eof() )
  {
    getline( infile, temp );
    if ( temp.size() > 0 )
    {
      line += temp;
      temp.erase();
    }
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find community type in " << rootnode;
    cerr << " !!!" << endl;
    exit( -1 );
  }

  startString = line.find( "<community type = " );
  temp = line.substr( startString, 30 );
  startString = temp.find( '"' );
  value = temp.substr( (startString+1), 5 );
  
  return atoi( value.c_str() );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double ProcessXML60::getXMLcmntArrayDouble( ifstream& infile,
                                            const string& rootnode,
                                            const string& varnode,
                                            const int& index )

{
  string line;
  string value;

  string endVarnode = "</" + varnode + ">";
  unsigned int startString;

  while ( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for community type ";
    cerr << index << " in " << rootnode << endl;
    exit( -1 );
  }

  startString = line.find( ">" );
  if ( startString == string::npos ) { return MISSING; }
  else
  {
    value = line.substr( (startString+1), 20 );
    return atof( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML60::getXMLcmntArrayInt( ifstream& infile,
                                      const string& rootnode,
                                      const string& varnode,
                                      const int& index )
{
  string line;
  string value;

  string endVarnode = "</" + varnode + ">";
  unsigned int startString;

  while ( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for community type ";
    cerr << index << " in " << rootnode << endl;
    exit( -1 );
  }

  startString = line.find( ">" );
  if ( startString == string::npos ) { return -99; }
  else
  {
    value = line.substr( (startString+1), 20 );
    return atoi( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

long ProcessXML60::getXMLcmntArrayLong( ifstream& infile,
                                        const string& rootnode,
                                        const string& varnode,
                                        const int& index )
{
  string line;
  string value;

  string endVarnode = "</" + varnode + ">";
  unsigned int startString;

  while ( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for community type ";
    cerr << index << " in " << rootnode << endl;
    exit( -1 );
  }

  startString = line.find( ">" );
  if ( startString == string::npos ) { return -99; }
  else
  {
    value = line.substr( (startString+1), 20 );
    return atol( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double ProcessXML60::getXMLdouble( ifstream& infile,
                                   const string& rootnode,
                                   const string& varnode )
{
  string line;
  string value;

  string endVarnode = "</" + varnode + ">";
  unsigned int startString;

  while ( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode;
    cerr << " in " << rootnode << endl;
    exit( -1 );
  }

  startString = line.find( ">" );
  if ( startString == string::npos ) { return MISSING; }
  else
  {
    value = line.substr( (startString+1), 20 );
    return atof( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML60::getXMLint( ifstream& infile,
                             const string& rootnode,
                             const string& varnode )
{
  string line;
  string value;

  string endVarnode = "</" + varnode + ">";
  unsigned int startString;

  while ( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode;
    cerr << " in " << rootnode << endl;
    exit( -1 );
  }

  startString = line.find( ">" );
  if ( startString == string::npos ) { return -99; }
  else
  {
    value = line.substr( (startString+1), 20 );
    return atoi( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

long ProcessXML60::getXMLlong( ifstream& infile,
                               const string& rootnode,
                               const string& varnode )
{
  string line;
  string value;

  string endVarnode = "</" + varnode + ">";
  unsigned int startString;

  while ( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode;
    cerr << " in " << rootnode << endl;
    exit( -1 );
  }

  startString = line.find( ">" );
  if ( startString == string::npos ) { return -99; }
  else
  {
    value = line.substr( (startString+1), 20 );
    return atol( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML60::getXMLrootNode( ifstream& infile,
                                  const string& rootnode )

{
  string line;

  while ( line.find( rootnode ) == string::npos && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << rootnode << " !!!" << endl;
    exit( -1 );
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML60::getXMLsiteCommunityNode( ifstream& infile,
                                           const string& rootnode,
                                           string& description )

{
  string line;
  string value;

  string temp;
  string temp2;

  int startString;
  int endString;

  while ( line.find( ">" ) == string::npos
          && !infile.eof() )
  {
    getline( infile, temp );
    if ( temp.size() > 0 )
    {
      line += temp;
      temp.erase();
    }
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find community type in " << rootnode;
    cerr << " !!!" << endl;
    exit( -1 );
  }

  startString = line.find( "<community type = " );
  temp = line.substr( startString, 30 );
  startString = temp.find( '"' );
  value = temp.substr( (startString+1), 5 );

  temp.erase();
  startString = line.find( "description" );
  temp = line.substr( startString, 50 );
  startString = temp.find( '"' );
  temp2 = temp.substr( (startString+1), 50 );
  endString = temp2.find( '"' );
  description = temp2.substr( 0, endString );

  return atoi( value.c_str() );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML60::getXMLsiteRootNode( ifstream& infile,
                                      const string& rootnode,
                                      string& version,
                                      string& sitename,
                                      string& sitecol,
                                      string& siterow,
                                      string& developer,
                                      string& updated )
{

  string line;
  string temp;
  string temp2;

  int startString;
  int endString;

  while ( line.find( rootnode ) == string::npos && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << rootnode << " !!!" << endl;
    exit( -1 );
  }

  while ( line.find( ">" ) == string::npos && !infile.eof() )
  {
    getline( infile, temp );
    if ( temp.size() > 0 )
    {
      line += temp;
      temp.erase();
    }
  }


  startString = line.find( "version" );
  temp = line.substr( startString, 50 );
  startString = temp.find( '"' );
  temp2 = temp.substr( (startString+1), 50 );
  endString = temp2.find( '"' );
  version = temp2.substr( 0, endString );

  temp.erase();
  startString = line.find( "site" );
  temp = line.substr( startString, 80 );
  startString = temp.find( '"' );
  temp2 = temp.substr( (startString+1), 80 );
  endString = temp2.find( '"' );
  sitename = temp2.substr( 0, endString );

  temp.erase();
  startString = line.find( "longitude" );
  temp = line.substr( startString, 50 );
  startString = temp.find( '"' );
  temp2 = temp.substr( (startString+1), 50 );
  endString = temp2.find( '"' );
  sitecol = temp2.substr( 0, endString );

  temp.erase();
  startString = line.find( "latitude" );
  temp = line.substr( startString, 50 );
  startString = temp.find( '"' );
  temp2 = temp.substr( (startString+1), 50 );
  endString = temp2.find( '"' );
  siterow = temp2.substr( 0, endString );

  temp.erase();
  startString = line.find( "developedBy" );
  temp = line.substr( startString, 80 );
  startString = temp.find( '"' );
  temp2 = temp.substr( (startString+1), 80 );
  endString = temp2.find( '"' );
  developer = temp2.substr( 0, endString );

  temp.erase();
  startString = line.find( "updated" );
  temp = line.substr( startString, 50 );
  startString = temp.find( '"' );
  temp2 = temp.substr( (startString+1), 50 );
  endString = temp2.find( '"' );
  updated = temp2.substr( 0, endString );

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML60::getXMLtemvegNode( ifstream& infile,
                                    const string& rootnode )
{
  string line;
  string value;

  string temp;
  int startString;


  while ( line.find( ">" ) == string::npos
          && !infile.eof() )
  {
    getline( infile, temp );
    if ( temp.size() > 0 )
    {
      line += temp;
      temp.erase();
    }
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find TEMVEG type in " << rootnode;
    cerr << " !!!" << endl;
    exit( -1 );
  }

  startString = line.find( "<temveg type = " );
  temp = line.substr( startString, 30 );
  startString = temp.find( '"' );
  value = temp.substr( (startString+1), 5 );

  return atoi( value.c_str() );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double ProcessXML60::getXMLtvegArrayDouble( ifstream& infile,
                                            const string& rootnode,
                                            const string& varnode,
                                            const int& index )

{
  string line;
  string value;

  string endVarnode = "</" + varnode + ">";
  unsigned int startString;

  while ( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for TEMVEG = ";
    cerr << index << " in " << rootnode << endl;
    exit( -1 );
  }

  startString = line.find( ">" );
  if ( startString == string::npos ) { return MISSING; }
  else
  {
    value = line.substr( (startString+1), 20 );
    return atof( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML60::getXMLtvegArrayInt( ifstream& infile,
                                      const string& rootnode,
                                      const string& varnode,
                                      const int& index )
{
  string line;
  string value;

  string endVarnode = "</" + varnode + ">";
  unsigned int startString;

  while ( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if ( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for TEMVEG = ";
    cerr << index << " in " << rootnode << endl;
    exit( -1 );
  }

  startString = line.find( ">" );
  if ( startString == string::npos ) { return (int) MISSING; }
  else
  {
    value = line.substr( (startString+1), 20 );
    return atoi( value.c_str() );
  }

};


