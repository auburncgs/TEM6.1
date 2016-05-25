/* **************************************************************
*****************************************************************
VOC60.HPP - describes components comprising inorganic nitrogen
*****************************************************************
************************************************************** */

#ifndef VOC60_H
#define VOC60_H

class VolatileOrganicCarbon60
{

  public:

/* **************************************************************
		 Public Variables
************************************************************** */

     double total;          // (Units are grams C / square meter)
     double isoprene;       // (Units are grams C / square meter)
     double monoterpene;    // (Units are grams C / square meter)
     double other;          // (Units are grams C / square meter)
     double otherReactive;  // (Units are grams C / square meter)

};

#endif
