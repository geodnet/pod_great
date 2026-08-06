
#ifndef  GALLOTL_H
#define  GALLOTL_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: container along sites for BLQ data
  Version: $Rev:$

  2013-03-29 /PV: created

  -*/

#include <string>
#include <map>

#include "newmat/newmat.h"
#include "newmat//newmatio.h"

#include "gdata/gdata.h"
#include "gutils/gtime.h"
#include "gmodels/gotl.h"

using namespace std;

namespace gnut
{
	class BASE_LIBRARY_EXPORT t_gallotl : public t_gdata
	{
	public:
		t_gallotl();
		~t_gallotl();

		int data(Matrix& data, double lon, double lat);
		void add(t_gotl& otl);
	private:
		map<string, t_gotl> _mapotl;
	};

} // namespace

#endif
