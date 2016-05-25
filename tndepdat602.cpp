/* *************************************************************
TNDEPDAT602.CPP - object to read and write the structure of the
                N deposition data from files used by the 
                Terrestrial Ecosystem Model (TEM)

Modifications:
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
  

#include "tndepdat602.h"

Ndepdata60::Ndepdata60( void )
{

  ndepend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

int Ndepdata60::get( ifstream& ifile )
{

  lagpos = ifile.tellg();

  ifile >> col;
  ifile >> row;
  ifile >> varname;
  ifile >> carea;
  ifile >> year;

  ifile >> nhx;
  ifile >> noy;
  ifile >> totndep;

  ifile >> contnent;

  ifile.seekg( 0, ios::cur );
  
  curpos = ifile.tellg();

  if(curpos < (lagpos + 10)) { ndepend = -1; }

  return ndepend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Ndepdata60::getdel( FILE* ifile )
{
  char tmpvarname[80];
  char tmpcontnent[80];
  
  ndepend = fscanf( ifile, 
                    "%lf,%lf, %s ,%lf,%d,%lf,%lf,%lf, %s",
                    &col, 
                    &row, 
                    tmpvarname, 
                    &carea, 
                    &year, 
                    &nhx, 
                    &noy, 
                    &totndep, 
                    tmpcontnent );

  varname = tmpvarname;
  contnent = tmpcontnent;
  
  return ndepend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Ndepdata60::out( ofstream& ofile, 
                      const double& col, 
                      const double& row, 
                      const string& varname, 
                      const double& carea, 
                      const int& year, 
                      const double& nhx, 
                      const double& noy, 
                      const double& totndep, 
                      const string& contnent )

{

  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ' ';
  ofile << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision( 2 ) << carea << ' ';
  ofile << year << ' ';

  ofile << setprecision( 5 ) << nhx << ' ';
  ofile << setprecision( 5 ) << noy << ' ';
  ofile << setprecision( 5 ) << totndep << ' ';

  ofile << contnent << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Ndepdata60::outdel( ofstream& ofile, 
                         const double& col, 
                         const double& row, 
                         const string& varname, 
                         const double& carea, 
                         const int& year, 
                         const double& nhx, 
                         const double& noy, 
                         const double& totndep, 
                         const string& contnent )

{
  ofile.setf( ios::fixed,ios::floatfield );
  ofile.setf( ios::showpoint );
  ofile.precision( 1 );

  ofile << col << ",";
  ofile << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision( 2 ) << carea << ",";
  ofile << setprecision( 0 ) << year << ",";
  ofile << setprecision( 5 ) << nhx << ",";
  ofile << setprecision( 5 ) << noy << ",";
  ofile << setprecision( 5 ) << totndep << ",";

  ofile << " " << contnent;
  ofile << endl;

};



