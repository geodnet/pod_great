#include "gfarg_init.h"
#include "gfarg_init.h"
/**
 *
 * @verbatim
	 History
	  -1.0	Jiande Huang	  2019-01-11	 creat the file.
	  -1.1  Jie Li	      2019-04-01	 add code in doxygen style.
	  -1.2  Yujie Qin     2019-04-15 change sth about leo
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gfarg_init.cpp
 * @brief		    The base class used to initialize the force model parameters
*
 * @author          Jiande Huang, Wuhan University
 * @version		    1.0.0
 * @date		    2019-01-06
 */

#include "gforce/gfarg_init.h"

namespace great
{
	t_init_gravity_egm08::t_init_gravity_egm08(FORCE_MODEL model, int order, t_gegmcoef* egmcoef, t_gnavde* navde, t_gallotl* otl, t_gocean_tide* ot) : t_init_arg(model)
	{
		gravity_order = order; ///< gravity order.
		gracoef = egmcoef;	   ///< EGM-file data.
		gplanet = navde;	   ///< 
		gotl = otl;			   ///< ocean load data.
		gotdata = ot;		   ///< ocean tide data.
	}
	/**
	* @brief constructor, initialize gravity EGM2008 model.
	* @param[in]   model      force model enum
	* @param[in]   order      EGM model order
	* @param[in]   egmcoef    EGM coefficient data from EGM file
	* @param[in]   navde      jpleph_de405 data
	* @param[in]   otl        oceanload data
	* @param[in]   ot         oceantide data
	*/
	t_init_gravity_egm08::t_init_gravity_egm08(FORCE_MODEL model, int order, t_gegmcoef * egmcoef,
		t_gnavde * navde, t_gallotl * otl, t_gocean_tide * ot, string otmode, modeofmeanpole meanpole) : t_init_arg(model)
	{
		gravity_order = order;		///< gravity order.
		gracoef = egmcoef;	///< EGM-file data.
		gplanet = navde;	///< DE.
		gotl = otl;		    ///< ocean load data.
		gotdata = ot;	    ///< ocean tide data.
		ot_mode = otmode;
		_mean_pole_model = meanpole;
	}


	/**
	 * @brief constructor, base initialize force model.
	 * @param[in]   model      force model enum
	 */
	t_init_arg::t_init_arg(FORCE_MODEL model)
	{
		force_model = model;
	}

	/**
	 * @brief constructor, initialize relativity model.
	 * @param[in]   model      force model enum
	 */
	t_init_relativity::t_init_relativity(FORCE_MODEL model) : t_init_arg(model)
	{
	}



	/**
	 * @brief constructor, initialize point_mass(N-body) model.
	 * @param[in]   model      force model enum
	 * @param[in]   navde	   JPL data
	 */
	t_init_point_mass::t_init_point_mass(FORCE_MODEL model,
		t_gnavde * navde) : t_init_arg(model)
	{
		gplanet = navde;
	}


	/**
	* @brief constructor, initialize customer model.
	* @param[in]   model      force model enum
	* @param[in]   pars_map   parameters list of ecom model
	*/
	t_init_solar_ecom2::t_init_solar_ecom2(FORCE_MODEL model,
		std::map<ICS_PAR, double> pars_map) : t_init_arg(model)
	{
		solar_pars = pars_map;
	}


	t_init_antennathrust::t_init_antennathrust(FORCE_MODEL model_, string ATmodel)
	{
		force_model = model_;
		model = ATmodel;
	}
}



