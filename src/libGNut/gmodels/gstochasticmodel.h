
#ifndef STOCHASTIC_H
#define STOCHASTIC_H 
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: stochastic models
  Version: $ Rev: $

  2012-05-02 /PV: created

-*/

#include "gutils/gtime.h"
#include "newmat/newmat.h"
#include "gset/gsetproc.h"
#include "Eigen/Dense"

using namespace std;

namespace gnut {   



class BASE_LIBRARY_EXPORT t_statemode
{
public:
	t_statemode();
	t_statemode(int order,double dt,double noise);
	virtual ~t_statemode();

	int order;
	// x(t+1) = M*x(t)  P
	Matrix M;
	SymmetricMatrix P;

private:
	static const double _coeff[6];
	
	
};

class BASE_LIBRARY_EXPORT t_orb_mode
{
public:

  t_orb_mode(pair<int, double> pv, pair<int, double> srp, pair<int, double> drag, pair<int, double> emp);

  
  Eigen::MatrixXd getR(double dt,double ratio=1.0);

  double get_vel_noise(double dt, double ratio=1.0){
      return _pv_par.second  *  ratio;
  }

	private:

	  pair<int, double> _pv_par;
	  pair<int, double> _srp_par;
	  pair<int, double> _drag_par;
	  pair<int, double> _emp_par;

		static const double _coeff[3];
};

} // namespace

#endif // STOCHASTIC_H
