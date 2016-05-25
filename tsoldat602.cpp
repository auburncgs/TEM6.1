/* *************************************************************
TSOLDAT433.CPP - object to read and write the structure of soil
                texture data from/to files used by the  
                the Terrestrial Ecosystem Model (TEM)

Modifications:

20051117 - DWK created by modifying tsoldat433.cpp
20051117 - DWK added include tsoldat602.h and standard includes
20051117 - DWK changed Soildata43:: to Soildata60::

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
  

#include "tsoldat602.h"



Soildata60::Soildata60( void )
{

  soilend = 1;
  lagpos = -99;
  curpos = 0;

};

/* *************************************************************
************************************************************* */

int Soildata60::get( ifstream& infile )
{

  lagpos = infile.tellg();

  infile >> col;
  infile >> row;
  infile >> varname;
  infile >> carea;
  infile >> pctsand;
  infile >> pctsilt;
  infile >> pctclay;
  infile >> wsoil;
  infile  >> source;
  infile >> contnent;

  infile.seekg( 0, ios::cur );
  
  curpos = infile.tellg();

  if ( curpos < (lagpos + 10) ) { soilend = -1; }

  return soilend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soildata60::getdel( FILE* infile )
{
  char tmpvarname[80];
  char tmpsource[80];
  char tmpcontnent[80];
  
  soilend = fscanf( infile,
                    "%lf,%lf, %s ,%lf,%lf,%lf,%lf,%d, %s , %s",
                    &col,
                    &row,
                    tmpvarname,
                    &carea,
                    &pctsand,
                    &pctsilt,
                    &pctclay,
                    &wsoil,
                    tmpsource,
                    tmpcontnent );
  
  varname = tmpvarname;
  source = tmpsource;
  contnent = tmpcontnent;
  
  return soilend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soildata60::out( ofstream& ofile, 
                      const double& col, 
                      const double& row, 
                      const string& varname,
                      const double& carea, 
                      const double& pctsand, 
                      const double& pctsilt, 
                      const double& pctclay,
                      const int& wsoil, 
                      const string& source, 
                      const string& contnent )
{

   ofile.setf( ios::fixed,ios::floatfield );
   ofile.setf( ios::showpoint );
   ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 2 ) << carea << ' ';
  ofile << setprecision( 2 ) << pctsand << ' ';
  ofile << pctsilt << ' ';
  ofile << pctclay << ' ';
  ofile << setprecision( 0 ) << wsoil << ' ';
  ofile << source << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soildata60::outdel( ofstream& ofile, 
                         const double& col, 
                         const double& row, 
                         const string& varname,
                         const double& carea, 
                         const double& pctsand, 
                         const double& pctsilt, 
                         const double& pctclay,
                         const int& wsoil, 
                         const string& source, 
                         const string& contnent )
{

   ofile.setf( ios::fixed,ios::floatfield );
   ofile.setf( ios::showpoint );
   ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision( 2 ) << carea << ",";
  ofile << setprecision( 2 ) << pctsand << ",";
  ofile << pctsilt << ",";
  ofile << pctclay << ",";
  ofile << setprecision( 0 ) << wsoil << ", ";
  ofile << source << " , ";
  ofile << contnent;
  ofile << endl;

};

