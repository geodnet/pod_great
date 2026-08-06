
#ifndef GMATRIXCONV_H
#define GMATRIXCONV_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: matrix conversion utilities
  Version: $ Rev: $

  2011-04-26 /PV: created
  2012-04-06 /JD: extracted matrix conversion utilities from old utils.h

-*/

#include <ostream>   // need for NEWMATIO !
#include <vector>

#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "ExportHeaders/ExportBase.h"

using namespace std;
namespace gnut
{

	// probably will be obsolete using gtriple everywhere !!!!!
	// -----
//	BASE_LIBRARY_EXPORT ColumnVector array2colVec(double*);                     // convert array to columnVector

	BASE_LIBRARY_EXPORT void Matrix_remRC(SymmetricMatrix&, int row, int col);  // remove   r_th row and c_th column in SymMatrix
	BASE_LIBRARY_EXPORT void Matrix_remRC(DiagonalMatrix&, int row);  // remove   r_th row and c_th column in DiagonalMatrix

	BASE_LIBRARY_EXPORT void Matrix_swap(SymmetricMatrix&, int a, int b);      // swap a_th and b_th row and column in SymMatrix

	BASE_LIBRARY_EXPORT Matrix rotX(double Angle);
	BASE_LIBRARY_EXPORT Matrix rotY(double Angle);
	BASE_LIBRARY_EXPORT Matrix rotZ(double Angle);

} // namespace

#endif