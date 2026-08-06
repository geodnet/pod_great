
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <cmath>
#include <iomanip>

#include "gmodels/gtropo.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"
#include "gutils/gsysconv.h"

using namespace std;

namespace gnut {   

t_gtropo::t_gtropo() // : t_gtropo::t_gtropo("")
{}


t_gtropo::~t_gtropo()
{}


// get ZHD
// ----------
double t_gtropo::getZHD(const t_gtriple& ell, const t_gtime& epo) // ell v RADIANECH !! TREBA SJEDNOTIT
{   
  gtrace("t_gtropo::getZHD");  
  
  double pp =  1013.25 * pow(1.0 - 2.26e-5 * ell[2], 5.225) ;
  double res = (0.002277 * pp) / (1.0 - 0.00266 * cos(2.0*ell[0]) - 0.00000028 * ell[2]);

  cerr << "t_gtropo:getZHD not available data " << epo.str_ymdhms() << ", using default: 2.3m \n";

  return res;
}


// get ZWD
// ----------
double t_gtropo::getZWD(const t_gtriple& ell, const t_gtime& epo) // ell v RADIANECH !! TREBA SJEDNOTIT
{ 
  gtrace("t_gtropo::getZWD");

  return 0.0; // NWM_UNKNOWN;
}

// ----------------------------- SITE MODELS ------------------------------------

// ---------
// Saastamoinen 1972
// ----------
double t_saast::getSTD(const double& ele, const double& height)
{    
  gtrace("t_saast::getSTD"); 
   
  double pp =  1013.25 * pow(1.0 - 2.26e-5 * height, 5.225) ;
  double TT =  18.0 - height * 0.0065 + 273.15;
  double hh =  50.0 * exp(-6.396e-4 * height);
  double ee =  hh / 100.0 * exp(-37.2465 + 0.213166*TT - 0.000256908*TT*TT);

      
  double h_km = height / 1000.0;

  if (h_km < 0.0) h_km = 0.0;
  if (h_km > 5.0) h_km = 5.0;
  int    ii   = int(h_km + 1);
  double href = ii - 1;
	
  double bCor[6];
  bCor[0] = 1.156;
  bCor[1] = 1.006;
  bCor[2] = 0.874;
  bCor[3] = 0.757;
  bCor[4] = 0.654;
  bCor[5] = 0.563;

  double BB = bCor[ii-1] + (bCor[ii]-bCor[ii-1]) * (h_km - href);
	
  double zen   = G_PI/2.0 - ele;
  double delay = (0.002277/cos(zen)) * (pp + ((1255.0/TT)+0.05)*ee - BB*(tan(zen)*tan(zen)));
  return delay;
}

// ----------
double t_saast::getZHD(const t_gtriple& Ell, const t_gtime& epoch)
{   
  gtrace("t_saast::getZHD");
   
  double P, T, N;
   
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P,T,N);
   
  double delay = (0.002277 * P) /
                 (1.0 - 0.00266 * cos(2.0*Ell[0]) - 0.00000028 * Ell[2]);

   
  return delay;
}

// ----------
double t_saast::getZWD(const t_gtriple& Ell, const t_gtime& epoch)
{   
  gtrace("t_saast::getZWD");   
   
  double P, T, N;   
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P, T, N);

  //add need test
  double hh = 0.6;  
  double e = hh * 6.11 * pow(10.0, (7.5 * T / (T + 237.3)));

  T += 273.15;   
   

  double delay = (0.0022768*e*(1255.0/T + 0.05))
               / (1.0 - 0.00266 * cos(2.0*Ell[0]) - 0.00000028 * Ell[2]);
  return delay;
}


} // namespace
