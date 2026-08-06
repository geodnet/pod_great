/**
 *
 * @verbatim
	 History
	  -1.0 Keke Zhang  2020-12-30 creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gfgravityeigen.cpp
 * @brief		    The class for EIGEN-6S4 time-variable satellite-only gravity model
 *                  Reference: Foerste et al. (2016, http://doi.org/10.5880/icgem.2016.004)
 * @author          Keke Zhang, Wuhan University
 * @version		    1.0.0
 * @date		    2020-12-30
 *
 */
#include "gutils/gconst.h"
#include "gutils/gconst_oi.h"
#include "gutils/ginfolog.h"
#include "gforce/gfgravityeigen.h"
#include "gforce/gfarg_cal.h"

using namespace std;

namespace great
{
	t_gfgravity_eigen::t_gfgravity_eigen(t_init_arg *arg, t_glog *log) : t_gfgravity_egm08(arg, log)
	{
		string function_name = "t_gfgravity_eigen::t_gfgravity_eigen";

		t_init_gravity_egm08 *eigen = nullptr;
		eigen = dynamic_cast<t_init_gravity_egm08 *>(arg);
		if (!eigen || !eigen->gracoef || !eigen->gotl || !eigen->gotdata || !eigen->gplanet)
		{
			write_log(log, 0, "ERROR : " + function_name + " ==> init failed");
			if (!eigen)
			{
				write_log(log, 0, "ERROR : " + function_name + " ==> no init argument");
			}

			// put the EGM data
			if (!eigen->gracoef)
			{
				write_log(log, 0, "ERROR : " + function_name + " ==> no EIGEN6S4 data, check your EIGEN6S4 file");
			}

			// put the oceanload data
			if (!eigen->gotl)
			{
				write_log(log, 0, "ERROR : " + function_name + " ==> no ocean load data, check your oceanload file");
			}

			// put the oceantide data
			if (!eigen->gotdata)
			{
				write_log_info(log, 0, "ERROR : " + function_name + " ==> no ocean ocean tide data, check your oceantide file");
			}

			// put the de nav data
			if (!eigen->gplanet)
			{
				write_log_info(log, 0, "ERROR : " + function_name + " ==> no ocean de nav data, check your de file");
			}

			// put the order of model
			if (eigen->gravity_order == 0)
			{
				write_log_info(log, 0, "ERROR : " + function_name + " ==> no degree order, check your xml file");
			}
		}
		else
		{
			_force_model = eigen->force_model;

			_gracoef = eigen->gracoef;
			_gotl = eigen->gotl;
			_gotdata = eigen->gotdata;
			_gplanets = eigen->gplanet;
			_gravi_degree = eigen->gravity_order;
			_ocean_degree = eigen->ocean_degree;
			_poletide = eigen->ispoletide;
			_otmode = eigen->ot_mode;

			_ocean_tide_Cnm.resize(_max_ocean_degree, _max_ocean_degree + 1);
			_ocean_tide_Snm.resize(_max_ocean_degree, _max_ocean_degree + 1);
			_ocean_tide_Fnm.resize(_max_ocean_degree, _max_ocean_degree + 1);
			_solid_tide_Cnm.resize(4, 4);
			_solid_tide_Snm.resize(4, 4);

			// force model for solide tide calculation(list: sun moom...)
			_planet_list.clear();
			_planet_list.insert("SUN");
			_planet_list.insert("MOON");

			// init the Matrix list
			_Cnm_save.resize(_gravi_degree + 1, _gravi_degree + 1);
			_Snm_save.resize(_gravi_degree + 1, _gravi_degree + 1);

			_Anm.resize(_gravi_degree + 3, _gravi_degree + 3);
			_SUM_F.resize(4, 1);

			// get the init Cnm and Snm
			if (!_isGetCnm())
			{
				write_log_info(_glog, 0, "ERROR : " + function_name + " ==> we set the Cnm zero");
			};

			if (!_isGetSnm())
			{
				write_log_info(_glog, 0, "ERROR : " + function_name + " ==> we set the Snm zero");
			};
		}
	}

	/**
	* @brief get the time variable spherical harmonic coefficients according to selected mjd
	* @remarks Reference: Foerste et al. (2016, http://doi.org/10.5880/icgem.2016.004)
	* @param [in]  dmjd
	*/
	bool t_gfgravity_eigen::_get_time_coeff(double &dmjd)
	{
		const t_time_gravity_coef& tcoef = _gracoef->get_time_coefficient(dmjd);
		double mjd_dif = (dmjd - tcoef.st.dmjd()) / 365.25;
		double p1 = 2 * G_PI / 1.0 * mjd_dif;
		double p2 = 2 * G_PI / 0.5 * mjd_dif;
		int degree = min(80, _gravi_degree);
		for (int n = 1; n <= degree; n++)
		{
			for (int m = 0; m <= n; m++)
			{
				_Cnm(n, m) = tcoef.c[n][m] + tcoef.ctrnd[n][m] * mjd_dif + tcoef.csin1[n][m] * sin(p1) +
					tcoef.ccos1[n][m] * cos(p1) + tcoef.csin2[n][m] * sin(p2) + tcoef.ccos2[n][m] * cos(p2);
				_Snm(n, m) = tcoef.s[n][m] + tcoef.strnd[n][m] * mjd_dif + tcoef.ssin1[n][m] * sin(p1) +
					tcoef.scos1[n][m] * cos(p1) + tcoef.ssin2[n][m] * sin(p2) + tcoef.scos2[n][m] * cos(p2);

			}
		}
		return true;
	}

	/** The acceleration and partial matrix will be calculated here.
	* @param [in] argument Parameters to be used in calculation.
	*/
	bool t_gfgravity_eigen::calculate(t_gforce_calarg *argument)
	{
		string function_name = "t_gfgravity_eigen::calculate";

		t_gfcalarg_gravity_egm *arg = NULL;
		arg = dynamic_cast<t_gfcalarg_gravity_egm *>(argument);

		if (!arg)
		{
			write_log_info(_glog, 0, "ERROR : " + function_name + " ==> no gravity_eigen argument");
			return false;
		}

		if (arg->arg_type != ARG_GRAVI_EGM)
		{
			write_log_info(_glog, 0, "ERROR : " + function_name + " ==> error gravity_eigen argument");
			return false;
		}

		try
		{
			bool isValid;
			// get the moon and sun pos
			_planet_crd_crs = *(arg->allXE);

			//get solitide and oceantide Cnm Snm correct for gravity
			isValid = _calculate_solid_tide(*(arg->time), *(arg->gst));
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_eigen::calculate ==> can not calculate solid tide!");
				return false;
			};

	
			isValid = _calculate_ocean_tide(*(arg->time), *(arg->gst));
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_eigen::calculate ==> can not calculate ocean tide!");
				return false;
			};

			// get the distance and normalize sat position
			_distance = arg->sat_pos_trs->NormFrobenius();
			_sat_pos_norm = *(arg->sat_pos_trs) / _distance;

			// get the coef from solid tide and ocean tide
			double mjd = *(arg->time);


			_Cnm = _Cnm_save;
			_Snm = _Snm_save;

			_get_time_coeff(*(arg->time));

            _Cnm(1,0)=0.0;
            _Cnm(1,1)=0.0;
            _Snm(1,1)=0.0;
			//mean pole according to IERS 2010
			double xpm = 0.0;
			double ypm = 0.0;
			_get_mean_pole(mjd, xpm, ypm);
			double m1 = *(arg->xpole) - xpm;
			double m2 = -(*(arg->ypole) - ypm);

			//Solid earth pole tide
			_Cnm(2, 1) = _Cnm(2, 1) - 1.333e-9 * (m1 + 0.0115 * m2);
			_Snm(2, 1) = _Snm(2, 1) - 1.333e-9 * (m2 - 0.0115 * m1);


			//ocean pole tide
			int max_degree;
			
				_Cnm(2, 1) = _Cnm(2, 1) - 2.1778e-10 * (m1 - 0.01724 * m2);
				_Snm(2, 1) = _Snm(2, 1) - 1.7232e-10 * (m2 - 0.03365 * m1);
			

			//Cnm Snm correction due to solid earth tide and/or ocean tide
			max_degree = max(3, _ocean_degree);
			for (int n = 2; n < max_degree + 1; n++)
			{
				for (int m = 0; m < max_degree + 1; m++)
				{
					if (n > _gravi_degree || m > _gravi_degree)
					{
						break;
					}

					_Cnm(n, m) = _Cnm(n, m) + _ocean_tide_Cnm(n - 1, m);
					_Snm(n, m) = _Snm(n, m) + _ocean_tide_Snm(n - 1, m);


					if (n <= 4 && m <= 3)
					{
						_Cnm(n, m) = _Cnm(n, m) + _solid_tide_Cnm(n - 1, m);
						_Snm(n, m) = _Snm(n, m) + _solid_tide_Snm(n - 1, m);

					}

				}
			}

			isValid = _cal_Anm_coeff();
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_eigen::calculate ==> can not calculate Anm!");
				return false;
			}

			isValid = _cal_acc_mat(*(arg->e2j));
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_eigen::calculate ==> can not calculate Acc!");
				return false;
			}

			isValid = _cal_Amat_coeff(*(arg->e2j));
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_eigen::calculate ==> can not calculate Amat!");
				return false;
			}

			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 0, "t_gfgravity_eigen::calculate", "ERROR : throw exception");
			return false;
		}
	};
}
