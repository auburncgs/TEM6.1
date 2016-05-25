/* *************************************************************
****************************************************************
TBIOME603C.H - object describing general characteristics of
             vegetation mosaic used in the Terrestrial Ecosystem
             Model (TEM)

Modifications:

20031019 - DWK created by modifying tbiome431a.h
20031019 - DWK changed include from tprocessXML431a.h to
           tprocessXML51.h
20031019 - DWK changed inheritance from ProcessXML to ProcessXML51
20031019 - DWK changed class Biome43 to class Biome51
20031019 - DWK changed public char cmnt_name[80] to string
           cmnt_name
20031019 - DWK changed include from tbiome431a.cpp to tbiome51.cpp
           at bottom of file
20031019 - DWK changed char ecd[MAXFGNAME] to const string& ecd
           in function call to getvtype()
20031021 - DWK added public function getVegSubarea()
20031021 - DWK added public double subarea
20040714 - DWK changed include from tprocessXML51.h to 
           tprocessXML60.h
20040714 - DWK changed class Biome51 to class Biome60
20040714 - DWK changed inheritance of ProcessXML to 
           inheritance of ProcessXML60
20040714 - DWK changed include from tbiome51.cpp to tbiome60.cpp
           at bottom of file
20040716 - DWK added public function getCommunityType()
20051117 - DWK added include temconsts602.hpp
20051117 - DWK changed include from tprocessXML60.h to
           tprocessXML602.h
20051117 - DWK deleted tbiom60.cpp from bottom of file
20070830 - DWK added public wfpsoff[NUMVEG]
20070830 - DWK added public getWFPSOFF() 
                                           
****************************************************************
************************************************************* */

#ifndef TBIOME603C_H
#define TBIOME603C_H

#include "temconsts602.hpp"
#include "tprocessXML602.h"


class Biome60 : public ProcessXML60
{

  public:

     Biome60( void );

 /* ************************************************************
		 Public Functions
************************************************************* */

     int   getCommunityType( const int& tveg );
     
     int   getVegMosaic( const int& tveg );
     
     double getVegSubarea( const int& tveg,
                           const int& dtype,
                           const int& carea );
     
      int   getVegSubtype( const int& tveg, const int& dtype );
     
     void   getvtype( ofstream& rflog1 );
     
     void   getvtype( const string& ecd );

     double getWFPSOFF( const int& tveg );

/* *************************************************************
		 Public Variables
************************************************************* */
#ifdef PMODE
     ifstream *fgo;    // the go file input stream
#endif

     // vegetation community type (categorical data)
     int cmnt;

     //Description of a vegetation community type
     string cmnt_name;

     // number of community types in a vegetation mosaic
     int numtype[NUMVEG];

     // community types of a vegetation mosaic
     int subtype[NUMVEG][NUMMSAC];

     // percent coverage of a community type in a vegetation
     //   mosaic
     double pcttype[NUMVEG][NUMMSAC];

     // Area covered by a vegetation community type
     double subarea;

     // biome type or ecozone (categorical data)
     int temveg;

     // Water-filled pore space offset
	 double wfpsoff[NUMVEG];
	 
};

#endif

