
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "gutils/gmatrixconv.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

void Matrix_remRC(SymmetricMatrix& Q, int row, int col)
{   
  SymmetricMatrix Qt(Q.Nrows());
  Qt = Q;
  Q.ReSize(Q.Nrows()-1);
  int rr = 1;
    for( int r = 1; r <= Qt.Nrows(); r++ ){
		   
      if( r == row ) continue;
      int cc = 1;
	       
      for( int c = 1; c <= Qt.Ncols(); c++ ){
	 
        if( c == col ) continue;
        Q(rr, cc) = Qt(r,c);
	cc++;
      }	
    rr++;
  }
}

void Matrix_remRC(DiagonalMatrix& D, int row)
{
  SymmetricMatrix Dt(D.Nrows());
  Dt = D;
  D.ReSize(D.Nrows()-1);
  int rr = 1;
    for( int r=1; r<=Dt.Nrows(); r++ ){
	if (r == row) continue;
	D(rr, rr) = Dt(r,r);
        rr++;  
    }	
}


//// Rotation Matrix
//// -----------------------------------
Matrix rotX(double Angle) {
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix UU(3,3);
    UU[0][0] = 1.0;  UU[0][1] = 0.0;  UU[0][2] = 0.0;
    UU[1][0] = 0.0;  UU[1][1] =  +C;  UU[1][2] =  +S;
    UU[2][0] = 0.0;  UU[2][1] =  -S;  UU[2][2] =  +C;
    return UU;
}
  
Matrix rotY(double Angle) {
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix UU(3,3);
    UU[0][0] =  +C;  UU[0][1] = 0.0;  UU[0][2] =  -S;
    UU[1][0] = 0.0;  UU[1][1] = 1.0;  UU[1][2] = 0.0;
    UU[2][0] =  +S;  UU[2][1] = 0.0;  UU[2][2] =  +C;
    return UU;
}
  
Matrix rotZ(double Angle) {
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix UU(3,3);
    UU[0][0] =  +C;  UU[0][1] =  +S;  UU[0][2] = 0.0;
    UU[1][0] =  -S;  UU[1][1] =  +C;  UU[1][2] = 0.0;
    UU[2][0] = 0.0;  UU[2][1] = 0.0;  UU[2][2] = 1.0;
    return UU;
}

} // namespace
