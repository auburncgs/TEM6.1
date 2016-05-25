/* *************************************************************
TCLMDAT602.H - object to read and write the structure of the
                climate data from files used by the Terrestrial
                Ecosystem Model (TEM)

20040516 - DWK created by modifying tclmdat425.h
20040516 - DWK changed class Clmdata to class Clmdata50
20040516 - DWK changed long year to int year in
           out(), outdel(), pctout() and pctoutdel()
20040516 - DWK changed char varname[9] to string varname in
           out(), outdel(), pctout() and pctoutdel()
20040516 - DWK changed char contnent[9] to string contnent in
           out(), outdel(), pctout() and pctoutdel()
20040516 - DWK changed public variable char varname[9] to
           string varname
20040516 - DWK changed public variable char contnent[9] to
           string contnent
20040516 - DWK changed include from tclmdat425.cpp to
           tclmdat50b3.cpp at bottom of file
20040516 - DWK changed arrays of size CYCLE to CYCLE+1
20040830 - DWK renamed Clmdat50 as Clmdat60
20040830 - DWK DWK changed include from tclmdat50b3.cpp to
           tclmdat60.cpp at bottom of file
20051117 - DWK deleted tclmdat60.cpp at bottom of file
20060606 - DWK changed public float col to public double col
20060606 - DWK changed public float row to public double row
                     
************************************************************* */

#ifndef TCLMDAT602_H
#define TCLMDAT602_H

#include "temconsts602.hpp"

class Clmdata60 
{
  
  public:
    
     Clmdata60( void );

/* *************************************************************
		      Public Functions
************************************************************* */

     // read data structure.
     
     int get( ifstream& ifile );
     int getdel( FILE* ifile );
     
     //write data structure.
     
     void out( ofstream& ofile, 
               const double& col, 
               const double& row, 
               const string& varname, 
               const double& carea, 
               const int& year, 
               double mon[CYCLE+1], 
               const string& contnent );
               
     void outdel( ofstream& ofile, 
                  const double& col, 
                  const double& row, 
                  const string& varname, 
                  const double& carea, 
                  const int& year, 
                  double mon[CYCLE+1], 
                  const string& contnent );
                  
     void pctout( ofstream& ofile, 
                  const double& col, 
                  const double& row, 
                  const string& varname, 
                  const double& carea, 
                  const int& year, 
                  double mon[CYCLE+1], 
                  const string& contnent );
                  
     void poutdel( ofstream& ofile, 
                   const double& col, 
                   const double& row, 
                   const string& varname, 
                   const double& carea, 
                   const int& year, 
                   double mon[CYCLE+1], 
                   const string& contnent );


/* *************************************************************
		     Public Variables
************************************************************* */

     // column or longitude of grid cell (degrees)	  
     double col;          

     // row or latitude of grid cell (degrees)
     double row;          
 
     // climate variable name
     string varname;     

     // area covered by grid cell (sq. km)
     double carea;          

      // date (year) of data
     int year;

     // annual sum of monthly data for grid cell
     double total;       

      // maximum monthly value for grid cell
     double max;        

     // mean annual value for grid cell
     double ave;         

     // minimum monthly value for grid cell
     double min;         

     // monthly values for the grid cell
     double mon[CYCLE+1];  

      // name of continent containing grid cell
     string contnent;   


  private:

/* *************************************************************
		      Private Variables
************************************************************* */

     int clmend;
     long curpos;
     long lagpos;

};

#endif
