/* **************************************************************
LATDAT602.H - object to read and write the structure of latitude
               and longitude data from/to files

Modifications:

20051117 - DWK created by modifying latdat425.h
20051117 - DWK changed class Latdata to class Latdata60
20051117 - DWK changed public char varname[9] to string varname
20051117 - DWK changed public char contnent9] to string contnent
20051117 - DWK deleted include latdat425.cpp from bottom of file
20060607 - DWK changed public float col to double col
20060607 - DWK changed public float row to double row

************************************************************** */

#ifndef LATDAT602_H
#define LATDAT602_H

class Latdata60
{

  public:

     Latdata60( void );

/* **************************************************************
                      Public Functions
************************************************************** */

// read data structure.
     int get( ifstream& infile );
     
     int getdel( FILE* infile );

//write data structure.

     void out( ofstream& ofile, 
               const double& col, 
               const double& row, 
               const string& varname, 
               const double& lat, 
               const double& lon, 
               const string& contnent );

     void outdel( ofstream& ofile, 
                  const double& col, 
                  const double& row, 
                  const string& varname, 
                  const double& lat, 
                  const double& lon, 
                  const string& contnent );


/* **************************************************************
                     Public Variables
************************************************************** */

     double col;        // column of grid cell
     double row;        // row or of grid cell
     string varname;    // "LATITUDE?"
     double lat;        // latitude of grid cell (degrees)
     double lon;        // longitude of grid cell (degrees)
     string contnent;   // name of continent containing grid cell


  private:

/* **************************************************************
                      Private Variables
************************************************************** */

     int latend;
     long curpos;
     long lagpos;

};

#endif

