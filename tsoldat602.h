/* *************************************************************
TSOLDAT602.H - object to read and write the structure of soil
               texture data from/to files used by the  
               Terrestrial Ecosystem Model (TEM)

Modifications:

20051117 - DWK created by modifying tsoldat433.h
20051117 - DWK changed class Soildata43 to class Soildata60
20051117 - DWK deleted tsoldat433.cpp from bottom of file
20060607 - DWK changed public float col to double col
20060607 - DWK changed public float row to double row

************************************************************* */

#ifndef TSOLDAT602_H
#define TSOLDAT602_H

class Soildata60
{
  
  public:
	 
     Soildata60( void );

/* *************************************************************
		    Public Functions
************************************************************* */
     
     int get( ifstream& infile );
     
     int getdel( FILE* infile );
     
     void out( ofstream& ofile, 
               const double& col, 
               const double& row, 
               const string& varname, 
               const double& carea, 
               const double& pctsand, 
               const double& pctsilt, 
               const double& pctclay, 
               const int& wsoil, 
               const string& source, 
               const string& contnent );
     
     void outdel( ofstream& ofile, 
                  const double& col, 
                  const double& row, 
                  const string& varname, 
                  const double& carea, 
                  const double& pctsand, 
                  const double& pctsilt, 
                  const double& pctclay, 
                  const int& wsoil, 
                  const string& source, 
                  const string& contnent );

/* *************************************************************
		     Public Variables
************************************************************* */

     // column or longitude of grid cell (degrees)	 
     double col;           
     
     // row or latitude of grid cell (degrees)
     double row;           
     
     // "TEXTURE"
     string varname;
     
     // area covered by grid cell (sq. km)     
     double carea;
     
     // percent sand of grid cell's soil texture           
     double pctsand;      
     
     // percent silt of grid cell's soil texture
     double pctsilt;      
     
     // percent clay of grid cell's soil texture
     double pctclay;      
     
     // wetland soil type designation (categorical data)
     int wsoil;           
     
     // reference to data source
     string source;
     
     // name of continent containing grid cell      
     string contnent;    


  private:

/* *************************************************************
		      Private Variables
************************************************************* */

     int soilend;
     long curpos;
     long lagpos;

};

#endif

