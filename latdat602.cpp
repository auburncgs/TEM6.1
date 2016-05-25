/* **************************************************************
LATDAT425.CPP - object to read and write the structure of land use
               data from/to files

Modifications:

20051117 - DWK created by modifying latdat425.cpp
20051117 - DWK added include latdat602.h and standard includes
20051117 - DWK changed Latdata:: to Latdata60::
20051117 - DWK changed char varname[9] to string varname
20051117 - DWK changed char contnent9] to string contnent
20060607 - DWK changed float col to double col
20060607 - DWK changed float row to double row

************************************************************** */

#include<cstdio>

  using std::fscanf;
  using std::FILE;

#include<iostream>

  using std::ios;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setprecision;

#include<string>

  using std::string;
  

#include "latdat602.h"

Latdata60::Latdata60( void )
{

  latend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Latdata60::get( ifstream& infile )
{

  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> lat;
  infile >> lon;
  infile >> contnent;

  infile.seekg( 0, ios::cur );
  
  curpos = infile.tellg();

  if(curpos < (lagpos + 10)) { latend = -1;}

  return latend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Latdata60::getdel( FILE* infile )
{
  char tmpvarname[80];
  char tmpcontnent[80];
  
  latend = fscanf( infile,
                   "%lf,%lf, %s ,%lf,%lf, %s",
                   &col,
                   &row,
                   tmpvarname,
                   &lat,
                   &lon,
                   tmpcontnent );

  varname = tmpvarname;
  contnent = tmpcontnent;
  
  return latend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Latdata60::out( ofstream& ofile, 
                     const double& col, 
                     const double& row, 
                     const string& varname,
                     const double& lat, 
                     const double& lon, 
                     const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 7 ) << lat << ' ';
  ofile << lon << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Latdata60::outdel( ofstream& ofile, 
                        const double& col, 
                        const double& row, 
                        const string& varname,
                        const double& lat, 
                        const double& lon, 
                        const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision( 7 ) << lat << ",";
  ofile << lon << ", ";
  ofile << contnent;
  ofile << endl;

};

