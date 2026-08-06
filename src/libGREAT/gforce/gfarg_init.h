/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang  2019-01-14 creat the file.
	  -1.1 Jie Li	     2019-03-31 add code in doxygen style
	  -1.2 Jiande Huang  2019-04-01 sdding abnormal throwing mechanism of functions
	  -1.3 Yujie Qin	 2019-04-15 change sth about leo
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfarg_init.h
 * @brief		The argument used in the initialization of force model.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */

#ifndef GFARG_INIT_H
#define GFARG_INIT_H

#include <map>
#include "gforce/gfcommon.h"
#include "gdata/gnavde.h"
#include "gdata/gegmcoef.h"
#include "gdata/gotdata.h"
#include "gdata/gdeadata.h"
#include "gdata/gsolargeomag.h"
#include "gall/gallotl.h"
#include "gutils/gconst_oi.h"
#include "gset/gsetbase.h"
#include <set>

namespace great
{
	/** @brief base calss for force model argument */
	class t_init_arg
	{
	public:
		/** @brief defult constructor */
		t_init_arg() { force_model = FORCE_MODEL::NO_DEF_F_MODEL; };

		/**
		 * @brief constructor, set force_model.
		 *
		 * @param[in]  model     force model type
		 */
		t_init_arg(FORCE_MODEL model);

        /** @brief defult destructor */		
		virtual ~t_init_arg() {};

		FORCE_MODEL force_model;	///< force model type.
	};



	/** @brief class for ecom2 solar radiation model argument */
	class t_init_solar_ecom2 : public t_init_arg
	{
	public:
		/** @brief defult constructor */
		t_init_solar_ecom2() {};

		/**
		* @brief constructor, set force_model and ics parameters value.
		*
		* @param[in]  model		force model used
		* @param[in]  pars_map     ics parameters list and value
		*/
		t_init_solar_ecom2(FORCE_MODEL model, std::map<ICS_PAR, double> pars_map);

		/** @brief defult destructor */
		virtual ~t_init_solar_ecom2() {};

		std::map<ICS_PAR, double> solar_pars;	///< initial value of ecom2 model.
	};

	
	/** @brief class for N-body gravity model argument */
	class t_init_point_mass : public t_init_arg
	{
	public:
		/** @brief defult constructor */
		t_init_point_mass() {};

		/**
		 * @brief constructor, set force_model and DE file.
		 *
		 * @param[in]  model		force model used
		 * @param[in]  navde		DE file data
		 */
		t_init_point_mass(FORCE_MODEL model, t_gnavde* navde);

		/** @brief defult destructor */
		virtual ~t_init_point_mass() {};

		t_gnavde* gplanet = nullptr;	///< JPL file data.
	};

	/** @brief class for EGM2008 model argument */
	class t_init_gravity_egm08 : public t_init_arg
	{
	public:
		/** @brief defult constructor */
		t_init_gravity_egm08() {};

		/**
		 * @brief constructor, set force_model, model order, and file data needed.
		 *
		 * @param[in]  model		force model used
		 * @param[in]  order		gravity model order used
		 * @param[in]  egmcoef		EGM-file data
		 * @param[in]  navde		DE file data
		 * @param[in]  otl			ocean load data
		 * @param[in]  ot			ocean tide data
		 */

		// init for default EGM model
		t_init_gravity_egm08(FORCE_MODEL model, int order, t_gegmcoef* egmcoef,
			t_gnavde* navde, t_gallotl* otl, t_gocean_tide* ot);

		// init for GNSS satellites (without set poletide and set modeofmeanpole)
		t_init_gravity_egm08(FORCE_MODEL model, int order, t_gegmcoef* egmcoef, 
			t_gnavde* navde, t_gallotl* otl, t_gocean_tide* ot, string otmode, modeofmeanpole meanpole);

		
		
		
		/** @brief defult destructor */
		virtual ~t_init_gravity_egm08() {};

		int           gravity_order  = 0;   ///< gravity order.
		int           ocean_degree = 30;
		int           pole_degree = 30;
		bool          ispoletide = false;
		t_gegmcoef*   gracoef   = nullptr;	///< EGM-file data.
		t_gnavde*     gplanet   = nullptr;	///< JPL file data.
		t_gallotl*    gotl      = nullptr;	///< ocean load data.
		t_gocean_tide*gotdata   = nullptr;	///< ocean tide data.
		t_gdeadata*   dedata    = nullptr;  ///< ocean pole data
		string ot_mode = "";				///< ocean tide model name
		modeofmeanpole _mean_pole_model = modeofmeanpole::cubic;
	};

	/** @brief class for reletivity model argument */
	class t_init_relativity : public t_init_arg
	{
	public:
		/** @brief defult constructor */
		t_init_relativity() {};

		/**
		 * @brief constructor, set force_model.
		 *
		 * @param[in]  model		force model used
		 */
		t_init_relativity(FORCE_MODEL model);

		/** @brief defult destructor */
		virtual ~t_init_relativity() {};
	};


	class t_init_antennathrust : public t_init_arg
	{
	public:
		t_init_antennathrust(FORCE_MODEL model_, string ATmodel);
		virtual ~t_init_antennathrust() {};
		string model;
	};
}
#endif // !GFINITARG_H
