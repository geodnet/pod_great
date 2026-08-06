
#ifndef GSYSCONV_H
#define GSYSCONV_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: reference system conversion utilities
  Version: $ Rev: $

  2011-04-26 /PV: created
  2012-04-06 /JD: extracted coordinate system conversion utilities only

-*/

#include <string>

#include "gutils/gtriple.h"
#include "gdata/gsatdata.h"

#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "ExportHeaders/ExportBase.h"


using namespace std;

namespace gnut {

	BASE_LIBRARY_EXPORT int ell2xyz(const double* Ell, double* XYZ, bool degrees);           // True/False:  LatLon in Degrees/Radians!
	BASE_LIBRARY_EXPORT int xyz2ell(const double* XYZ, double* Ell, bool degrees);           // True/False:  LatLon in Degrees/Radians!
	BASE_LIBRARY_EXPORT int xyz2neu(const double* XYZ, const double* XYZ_Ref, double* neu);
	BASE_LIBRARY_EXPORT int neu2xyz(const double* Ell, const double* neu, double* xyz);  // ! Radians only: Ell[0], Ell[1]

	BASE_LIBRARY_EXPORT int rao2xyz_rot(const ColumnVector& pos, const ColumnVector& vel, Matrix& R);

	BASE_LIBRARY_EXPORT int rao2xyz(const ColumnVector& pos, const ColumnVector& vel,
		const ColumnVector& rao, ColumnVector& xyz);


	BASE_LIBRARY_EXPORT int  ell2xyz(const t_gtriple& ell, t_gtriple& xyz, bool degrees);
	BASE_LIBRARY_EXPORT int  xyz2ell(const t_gtriple& crd, t_gtriple& ell, bool degrees);    // True/False:  LatLon in Degrees/Radians!
	BASE_LIBRARY_EXPORT void xyz2neu(t_gtriple& ell, t_gtriple& xyz, t_gtriple& neu);        // Radians only: Ell[0], Ell[1]
	BASE_LIBRARY_EXPORT void neu2xyz(t_gtriple& ell, t_gtriple& neu, t_gtriple& xyz);        // Radians only: Ell[0], Ell[1]



} // namespace

#endif // GSYSCONV_H
