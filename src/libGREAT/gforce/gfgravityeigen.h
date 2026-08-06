/**
 *
 * @verbatim
	 History
	  -1.0 Keke Zhang   2020-12-27 Creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gfgravityeigen.h
 * @brief		The class for EIGEN-6S4 time-variable satellite-only gravity model
 *              Reference: Foerste et al. (2016, http://doi.org/10.5880/icgem.2016.004)
 *
 * @author      Keke Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-12-27
 *
 */

#ifndef GFGRAVITYEIGEN_H
#define GFGRAVITYEIGEN_H

#include "gforce/gfgravityegm.h"
#include "gdata/gegmcoef.h"
#include "gforce/gfarg_init.h"

using namespace std;

namespace great
{
	/**
	* @brief EIGEN-6S4  gravity model,
	* A time-variable global satellite-only gravity field model of maximum d/o 150, containing time variable parameters
	* for all spher. harm. coeff. up to degree 80 (drift parameters per year, annual and semiannual terms).
	* The computation of the values of the time variable spherical coefficients G at a certain epoch t within the
	* respective time span between t0 and t1 must be done by applying the drift parameters trnd as well as the periodic
	* terms acos and asin in the following way. In the given rule the two periodic terms
	* are denoted by asin1/acos1 and asin2/acos2 in accordance with the two periods p1 = 1.0 year and p2 = 0.5 year.
	* G(t)=gfct(t0)+trnd*(t-t0)+asin1*sin(2pi/p1*(t-t0))+acos1*cos(2pi/p1*(t-t0))
						 +asin2*sin(2pi/p2*(t-t0))+acos2*cos(2pi/p2*(t-t0))
	* where t0 <= t <= t1 for the respective gfct, trnd and periodic terms.
	*
	* Similar to EGM2008, solid tide, ocean tide and induced pole tide are also considered in the Cnm and Snm coefficient.
	*
	*/
	class t_gfgravity_eigen : public t_gfgravity_egm08
	{
	public:
		/**
		 * @brief constructor.
		 * @remarks
		 * @param [in]  arg		initial argument
		 * @param [in]  log		log file
		 */
		t_gfgravity_eigen(t_init_arg *arg, t_glog *log);

		/** @brief defult destructor. */
		virtual ~t_gfgravity_eigen(){};

		/**
		 * @brief show the information of model.
		 */
		void show();

		/**
		 * @brief check some important pars in the model.
		 * @return
			 @retval true	data has no problem
			 @retval false  some wrong with data
		*/
		bool check();

		/**
		* @brief calculate the acceleration and partial matrix e.
		* @param [in] argument		Parameters to be used in calculation.
		* @return
			@retval true	calulate sucessfully
			@retval false   calulate unsucessfully
		*/
		bool calculate(t_gforce_calarg *argument);

	protected:
		/**
		* @brief get the time variable spherical harmonic coefficients according to selected mjd
		* @remarks Reference: Foerste et al. (2016, http://doi.org/10.5880/icgem.2016.004)
		* @return
			@retval true	get time variable coefficients sucessfully
			@retval false   get time variable coefficients unsucessfully
		*/
		bool _get_time_coeff(double &mjd);
	};
};
#endif