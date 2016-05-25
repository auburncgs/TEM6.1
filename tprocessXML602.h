/* *************************************************************
TPROCESSMXL602.H - Functions to process XML files for TEM

Programmer: Shaomin Hu and David Kicklighter
Creation Date: 17 June 2003

Modifications:

20031019 - DWK changed class ProcessXML to class ProcessXML51
20031019 - DWK changed char rootnode[80] to string rootnode
           and char varnode[80] to string varnode in function
          calls
20031019 - DWK changed include from tprocessXML431a.cpp to
           tprocessXML51.cpp
20031020 - DWK deleted private functions getvalue() and
           readline()
20040707 - DWK changed class ProcessXML51 to class ProcessXML60
20040707 - DWK changed include from tprocessXML51.cpp to
           tprocessXML60.cpp at bottom of file
20051117 - DWK deleted include tprocessXML60.cpp from bottom of
           file
20051117 - DWK added include temconst602.hpp
                       
************************************************************** */

#ifndef PROCESSXML602_H
#define PROCESSXML602_H

#include "temconsts602.hpp"

class ProcessXML60
{
  public:

  ProcessXML60();

  void endXMLcommunityNode( ifstream& infile );
  void endXMLtvegNode( ifstream& infile );

  double getXMLcmntArrayDouble( ifstream& infile,
                                const string& rootnode,
                                const string& varnode,
                                const int& index );

  int getXMLcmntArrayInt( ifstream& infile,
                          const string& rootnode,
                          const string& varnode,
                          const int& index );

  long getXMLcmntArrayLong( ifstream& infile,
                            const string& rootnode,
                            const string& varnode,
                            const int& index );

  int getXMLcommunityNode( ifstream& infile,
                           const string& rootnode );

  double getXMLdouble( ifstream& infile,
                       const string& rootnode,
                       const string& varnode );

  int getXMLint( ifstream& infile,
                 const string& rootnode,
                 const string& varnode );

  long getXMLlong( ifstream& infile,
                   const string& rootnode,
                   const string& varnode );

  int getXMLrootNode( ifstream& infile, const string& rootnode );

  int getXMLsiteCommunityNode( ifstream& infile,
                               const string& rootnode,
                               string& description );

  int getXMLsiteRootNode( ifstream& infile,
                          const string& rootnode,
                          string& version,
                          string& sitename,
                          string& sitecol,
                          string& siterow,
                          string& developer,
                          string& updated );

  int getXMLtemvegNode( ifstream& infile, const string& rootnode );

  double getXMLtvegArrayDouble( ifstream& infile,
                                const string& rootnode,
                                const string& varnode,
                                const int& index );

  int getXMLtvegArrayInt( ifstream& infile,
                          const string& rootnode,
                          const string& varnode,
                          const int& index );

};

#endif
