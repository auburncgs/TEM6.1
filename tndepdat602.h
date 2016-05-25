/* *************************************************************
TNDEPDAT602.H - object to read and write the structure of the
                annual N deposition data from files used by the 
                Terrestrial Ecosystem Model (TEM)

Modifications:

20051117 - DWK deleted tndepdat60.cpp from bottom of file
20060607 - DWK changed public float col to public double col
20060607 - DWK changed public float row to public double row

************************************************************* */

#ifndef TNDEPDAT602_H
#define TNDEPDAT602_H

class Ndepdata60 
{
  
  public:
    
     Ndepdata60( void );

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
               const double& nhx, 
               const double& noy, 
               const double& totndep, 
               const string& contnent );
               
     void outdel( ofstream& ofile, 
                  const double& col, 
                  const double& row, 
                  const string& varname, 
                  const double& carea, 
                  const int& year, 
                  const double& nhx, 
                  const double& noy, 
                  const double& totndep, 
                  const string& contnent );


/* *************************************************************
		     Public Variables
************************************************************* */
  
     // column or longitude of grid cell (degrees)
     double col;          

     // row or latitude of grid cell (degrees)
     double row;          

     // Variable name - i.e. NDEP
     string varname;    
     
     // area covered by grid cell (sq. km)
     double carea;          
     
     // date (year) of data
     int year;          
     
      // annual NHx deposition (mg N m-2 yr-1)
     double nhx;        
     
     // annual NOy deposition (mg N m-2 yr-1)
     double noy;         
     
     // annual total N deposition (mg N m-2 yr-1)
     double totndep;     
     
     // name of continent containing grid cell
     string contnent;   


  private:

/* *************************************************************
		      Private Variables
************************************************************* */

     int ndepend;
     long curpos;
     long lagpos;

};

#endif
