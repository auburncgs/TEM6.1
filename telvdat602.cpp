/* *************************************************************
TELVDAT602.CPP - object to read and write the structure of 
                 elevation data from/to files used by the Water 
                 Balance Model

Modifications:

20051117 - DWK created by modifying telvdat433.cpp
20051117 - DWK added include telvdat602.h and standard includes
20051117 - DWK changed Elevdata43:: to Elevdata60::
20060607 - DWK changed float col to double col
20060607 - DWK changed float row to double row

****************************************************************
************************************************************* */

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
  

#include "telvdat602.h"



Elevdata60::Elevdata60( void )
{

  elvend = 1;
  lagpos = -99;
  curpos = 0;

};

/* *************************************************************
************************************************************* */

int Elevdata60::get( ifstream& infile )
{

  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> carea;
  infile >> elev;
  infile >> contnent;

  infile.seekg( 0, ios::cur );
  
  curpos = infile.tellg();

  if( curpos < (lagpos + 10) ) { elvend = -1; }

  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Elevdata60::getdel( FILE* infile )
{
  char tmpvarname[80];
  char tmpcontnent[80];
  
  elvend = fscanf( infile,
                   "%lf,%lf, %s ,%lf,%lf, %s",
                   &col,
                   &row,
                   tmpvarname,
                   &carea,
                   &elev,
                   tmpcontnent );

  varname = tmpvarname;
  contnent = tmpcontnent;

  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Elevdata60::out( ofstream& ofile, 
                      const double& col, 
                      const double& row, 
                      const string& varname,
                      const double& carea, 
                      const double& elev, 
                      const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 2 ) << carea << ' ';
  ofile << setprecision( 1 ) << elev << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Elevdata60::outdel( ofstream& ofile, 
                         const double& col, 
                         const double& row, 
                         const string& varname,
                         const double& carea, 
                         const double& elev, 
                         const string& contnent )
{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision( 2 ) << carea << ",";
  ofile << setprecision( 1 ) << elev << ", ";
  ofile << contnent;
  ofile << endl;

};

