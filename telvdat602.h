/* *************************************************************
TELVDAT602.H - object to read and write the structure of 
               elevation data from/to files used by the Water 
               Balance Model

Modifications:

20051117 - DWK created by modifying telvdat433.h
20051117 - DWK changed class Elevdata43 to class Elevdata60
20040601 - DWK deleted telvdat433.cpp from bottom of file
20060607 - DWK changed public float col to double col
20060607 - DWK changed public float row to double row

****************************************************************
************************************************************* */

#ifndef TELVDAT602_H
#define TELVDAT602_H

class Elevdata60 
{
   
  public:
          
     Elevdata60( void );

/* *************************************************************
                      Public Functions
************************************************************* */

     // read data structure.
     int get( ifstream& infile );
     
     int getdel( FILE* infile );

     //write data structure.
     void out( ofstream& ofile, 
               const double& col, 
               const double& row, 
               const string& varname, 
               const double& carea, 
               const double& elev, 
               const string& contnent );
     
     void outdel( ofstream& ofile, 
                  const double& col, 
                  const double& row, 
                  const string& varname, 
                  const double& carea, 
                  const double& elev, 
                  const string& contnent );

          
/* *************************************************************
                     Public Variables
************************************************************* */

     // column or longitude of grid cell (degrees)
     double col;          
     
     // row or latitude of grid cell (degrees)
     double row;          
     
     // "ELEV"
     string varname;
     
     // area covered by grid cell (sq. km)    
     double carea;          
     
     // elevation of grid cell (m)
     double elev;        
     
     // name of continent containing grid cell
     string contnent;   


  private:

/* *************************************************************
                      Private Variables
************************************************************* */

     int elvend;
     long curpos;
     long lagpos;

};

#endif

