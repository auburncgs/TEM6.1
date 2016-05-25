/* *************************************************************
ELMNT602.H - Contains functions to manage elements used in GIS
             (Note: Elements are usually grid cells)
             
Modifications:

20051117 - DWK created by modifying elmnt425.h
20051117 - DWK deleted include elmnt425.cpp from bottom of file
20051117 - DWK changed char varname[9] to const string& varname 
           in public function call to coregrerr
20051117 - DWK changed public char contnent[9] to 
           public string contnent
20051117 - DWK changed public int carea to public long carea
20051117 - DWK change class Elmnt to class Elmnt60
20060607 - DWK changed public float col to public double col
20060607 - DWK changed public float row to public double row
                      
************************************************************** */

#ifndef ELMNT602_H
#define ELMNT602_H

class Elmnt60 
{
   
   public:

     Elmnt60( void );

/* **************************************************************
                 Public Functions
************************************************************** */

     void ask( ofstream& rflog1 );
     
     int coregerr( ofstream& rflog1, 
                   const string& varname1, 
                   const double& col1, 
                   const double& row1, 
                   const string& varname2, 
                   const double& col2, 
                   const double& row2 );
     
     void show( const int &nthproc, ofstream& rflog1,
                const double& col, 
                const double& row );
     
     void show( const int &nthproc, ofstream& rflog1,
                const double& col, 
                const double& row, 
                const long& totyr, 
                const double& tol );

/* **************************************************************
                 Public Variables
************************************************************** */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif
         
     // Column or longitude of element
     double col;

     // Area of element
     long carea;

     // Continent location of element
     string contnent;

     // Count of elements in a region
     long count;

     // Mean Elevation of element
     double elev;

     int end;

     long grdcnt;

     long numskip;

     // Row or latitude of element
     double row;

     int stopflag;

     int strtflag;

     int  totyr;


   private:

/* **************************************************************
                 Private Variables
************************************************************** */

     int endflag;
     long numgrids;
};

#endif

