
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include "gall/gallotl.h"
#include "gutils/gconst.h"

/// @note jdhuang add for linux cmake
#include <math.h>
namespace gnut {  

//Constructor
t_gallotl::t_gallotl()
{
//   cout << "vytvarim t_gallotl" << endl;
   id_type(  t_gdata::ALLOTL );
}

//Destructor
t_gallotl::~t_gallotl()
{
   _mapotl.clear();
}


int t_gallotl::data(Matrix& otldata, double lon,double lat)
{

	const double dlon_eps = 1E4 / 6378235.0 * R2D, dlat_eps = dlon_eps / cos(lat*D2R) ;

	if (lon > 180) lon -= 360.0;

	for (auto otl_iter = _mapotl.begin(); otl_iter != _mapotl.end(); otl_iter++)
	{
		double dlon, dlat;
		dlon = fabs(otl_iter->second.lon() - lon);
		dlat = fabs(otl_iter->second.lat() - lat);
		if (dlon < dlon_eps && dlat < dlat_eps) {
			otldata = otl_iter->second.data();
			return 1;
		}
	}

	return -1;
}



// Add element to map container
// -----------------------------
void t_gallotl::add(t_gotl& otl)
{

  _gmutex.lock();
   
  string site = otl.site();
  _mapotl[site] = otl;
   
  _gmutex.unlock(); return;
}



} // namespace
