
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <cmath>
#include <iomanip>

#include "gutils/gpair.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// null constructor
// ----------
t_gpair::t_gpair()
{
  _crd[0] = _crd[1] = 0.0;
}


// constructor
// ----------
t_gpair::t_gpair(double x, double y)
{
  _crd[0] = x;
  _crd[1] = y;
}


// constructor
// ----------
t_gpair::t_gpair(double crd[])
{
  _crd[0] = crd[0];
  _crd[1] = crd[1];
}


// constructor
// ----------
t_gpair::t_gpair(const ColumnVector& crd)
{
  _crd[0] = crd(1);
  _crd[1] = crd(2);
}


// destructor
// ----------
t_gpair::~t_gpair(){}


// get a reference of element
// ----------
double& t_gpair::operator[](const size_t idx)
{
  if( idx > 1 ){
    cerr << "Not valid pair index [used 0]\n";
    return _crd[0];
  }
   
  return _crd[idx];
}


// get a value of element
// ----------
double t_gpair::operator[](const size_t idx) const
{
  if( idx < 2 ) return _crd[idx];
   
  return 0.0;
}

// operator +
// ------------------------
t_gpair t_gpair::operator+(const t_gpair& other) const
{
  t_gpair tmp(*this);
  tmp[0] += other[0];
  tmp[1] += other[1];

  return tmp;   
}


// get single element
// ----------
double t_gpair::crd(int idx) const
{
  if( idx >= 0 && idx < 2 ) return _crd[static_cast<unsigned int>(idx)];
   
  return 0.0;
}





// copy operator
// ----------
t_gpair& t_gpair::operator=(const t_gpair& other)
{
  if( this != &other ){
    _crd[0] = other.crd(0);
    _crd[1] = other.crd(1);
  }
  return *this;
}

// equal operator
// ----------
bool t_gpair::operator==(const t_gpair& tr) const
{
  return ( _crd[0] == tr.crd(0) &&
           _crd[1] == tr.crd(1)     );
}


// operator for sorting
// ----------
bool t_gpair::operator<(const t_gpair& tr) const
{
  return ( ( _crd[0]  < tr.crd(0) ) ||
           ( _crd[0] == tr.crd(0) && _crd[1]  < tr.crd(1) ) );
}



// set array by ColumnVector
// ------------------------------
void t_gpair::set(const ColumnVector& crd)
{
  _crd[0] = crd(1);
  _crd[1] = crd(2);
}



// overloading << operator
// -----------------------------
ostream& operator<<(ostream& os, const t_gpair& x)
{
   os << fixed << setprecision(3) 
      << dbl2str(x[0]) + " " + dbl2str(x[1]);
   return os;
}



} // namespace
