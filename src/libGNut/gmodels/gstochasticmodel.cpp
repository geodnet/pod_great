
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <iomanip>

#include "gutils/gtypeconv.h"
#include "gmodels/gstochasticmodel.h"

namespace gnut {

const double t_statemode::_coeff[6] = { 1.0,-sqrt(3),2*sqrt(3),sqrt(5),-6*sqrt(5),12*sqrt(5) };

t_statemode::t_statemode() :
	order(0)
{

}

t_statemode::t_statemode(int order,double dt,double noise):
	order(order)
{
	// set the Matrix M P
	M = Matrix(order, order);

	// set Matrix M
	for (int Row = order; Row >= 1; Row--) {
		M(Row, Row) = 1.0;
		for (int Col = Row+1; Col <= order; Col++) {
			M(Row, Col) = M(Row, Col - 1) * dt / ((double)(Col - Row));
		}
	}

	// set Matrix P
	Matrix P_temp(order, order);
	int count = 0;
	for (int Row = order; Row >= 1; Row--) {
		for (int Col = order; Col >= Row; Col--) {
			P_temp(Row, Col) = _coeff[count] * pow(dt, Col - order);
			count++;
		}
	}

	P_temp = P_temp / (noise*sqrt(dt));

	P << P_temp.t() * P_temp;

}

t_statemode::~t_statemode()
{

}

const double t_orb_mode::_coeff[3] = { 1.0,-sqrt(3),2*sqrt(3) };

t_orb_mode::t_orb_mode(pair<int,double> pv, pair<int,double> srp, pair<int,double> drag, pair<int,double> emp):
	_pv_par(pv),_srp_par(srp),_drag_par(drag),_emp_par(emp)
{
}


Eigen::MatrixXd t_orb_mode::getR(double dt,double ratio)
{
	int par_num = _pv_par.first + _srp_par.first + _drag_par.first + _emp_par.first;
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(par_num,par_num); 
  // pos and vel part
  for (int i =0;i<3;i++){
    R(i,i) = _coeff[2]*pow(dt,-1.5) / (_pv_par.second*ratio);
    R(i,i+3) = _coeff[1]*pow(dt,-0.5) / (_pv_par.second*ratio);
    R(i+3,i+3) = _coeff[0]*pow(dt,-0.5) / (_pv_par.second*ratio);
  }

  // srp part
  for (int i = 6; i < 6+_srp_par.first; i++) {
	  R(i, i) = 1 / (_srp_par.second * ratio) * pow(dt, -0.5);
  }

  // drag part
  for (int i = 6 + _srp_par.first; i < 6 + _srp_par.first + _drag_par.first; i++) {
	  R(i, i) = 1 / (_drag_par.second * ratio) * pow(dt, -0.5);
  }

  // emp part
  for (int i = 6 + _srp_par.first + _drag_par.first; i < 6 + _srp_par.first + _drag_par.first + _emp_par.first; i++) {
	  R(i, i) = 1 / (_emp_par.second * ratio) * pow(dt, -0.5);
  }




  return R;
}



} // namespace
