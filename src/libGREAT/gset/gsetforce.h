/**
*
* The format :
* @verbatim
	<force_model_nav>
		<spacecraft ID="GNS">
			<integrator order="11" step_in="60" step_out="300" />
			<gravity model="EGM_12" valid="True" />
			<relativity model="Default" valid="True" />
			<point_mass model="Default" valid="True" />
			<solar_radiation model="ECOM_5" valid="True" />
		</spacecraft>
	</force_model_nav>
  @endverbatim
*
* @verbatim
	History
	 -1.0 Jiande Huang	 2019-01-04 creat the file.
	 -1.1 Jie Li		 2019-03-31 add doxygen style cade
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			getforce_nav.h
* @brief		force model set for navigation satellite orbit integral from XML
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#ifndef GSETFORCE_H
#define GSETFORCE_H

#define XMLKEY_FORCE             "force_model"
#define XMLKEY_FORCE_ID          "ID"
#define XMLKEY_FORCE_VALID       "valid"
#define XMLKEY_FORCE_MODEL       "model"
#define XMLKEY_FORCE_SPACE_CRAFT "space_craft"

#include <map>
#include <string>
#include <iostream>

#include "gutils/gnss.h"
#include "gset/gsetbase.h"

#include "ExportHeaders/ExportOI.h"
#include "gutils/gconst_oi.h"
#include "gforce/gfcommon.h"

using namespace std;

namespace great
{
	// string : sat name
	typedef map< string, map< FORCE, bool   > >	t_map_sys_force_model_valid;
	// first string : sat name, second string : force model name
	typedef map< string, map< FORCE, string > >	t_map_sys_force_model_name;

	/**
	 *@brief Class for GNS OI force model setting
	 */
	class OI_LIBRARY_EXPORT t_gsetforce : public virtual t_gsetbase
	{
	public:
		/** @brief default constructor. */
		t_gsetforce();

		/** @brief default destructor. */
		virtual ~t_gsetforce() {};

		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();

		/**
		* @brief get step-in of adams orbit integration(unit: s)
		* @param[in]  sat      satellite name
		* @return    OI input step
		*/
		int		getInpStep(string sat);

		/**
		* @brief get step-out of orbit integration(unit: s)
		* @param[in]  sat      satellite name
		* @return    OI output step
		*/
		int		getOutStep(string sat);

		/**
		* @brief get realtime update interval of orbit integration(unit: s)
		* @param[in]  sat      satellite name
		* @return    OI realtime update interval
		*/
		int		getRtUpdateInt(string sat);
		int		getInitTime(string sat);
		
		double	getPvNoise(string sat);
		double  getSrpNoise(string sat);
		double  getDragNoise(string sat);
		double	getEmpNoise(string sat);
		double getEclipseRatio(string sat);
		/**
		* @brief get the integra degree of adams orbit integration
		* @param[in]  sat      satellite name
		* @return    OI adams order
		*/
		int		getIntegraDegree(string sat);

		/**
		* @brief get gravity(ocean_tide and ocean_pole_tide) set
		* @param[in]  sat                  satellite name
		* @param[out] ot_degree            ocean_tide degree
		* @param[out] _opt_valid           ocean_pole_tide valid
		* @param[out] _opt_degree          ocean_pole_tide degree
		* @retrun
		*/
		void    getGravitySet(string sat, int& ot_degree, string& otmode, bool& opt_valid, int& opt_degree);

		

		/**
		* @brief get used force for GNS satellite(valid)
		* @param[in]  sat   satellite name
		* @return  the forces used
		*/
		set<FORCE>			getForceUsedSet(string sat);	

		/**
		*@brief get used force model for GNS satellite
		* @param[in]  sat  satellite name
		* @return  force model used
		*/
		map<FORCE, string>	getForceModelMap(string sat);

		/**
		* @brief get all force models
		* @param[in]  sat              satellite name
		* @return    map of force model(true or false)
		*/
		map<FORCE, bool>	getForcesValidMap(string sat);

		modeofmeanpole getMeanPoleModel(string sat);

	protected:
		t_map_sys_force_model_valid _map_force_valid;		///< if use the model
		t_map_sys_force_model_name  _map_force_model;		///< the model that used
		map<string, int> _inp_step;							///< adams inp step
		map<string, int> _out_step;							///< adams out step
		map<string, int> _integra_degree;					///< adams order
		map<string, int> _ot_degree;   ///< ocean_tide degree
		map<string, int> _opt_degree;  ///< ocean_pole_tide_degree
		map<string, string> _ot_mode;				  ///< ocean tide model name
		map<string, bool>_opt_valid;   ///< if use ocean_pole_tide
		map <string, modeofmeanpole> _mean_pole_model;

		map<string, int> _rt_update_interval;				///< orbit integration update interval in real-time processing
		map<string, int> _init_time;				///< initial time of orbit integratio in real-time processing

		map<string, double> _pv_noise;
		map<string, double> _srp_noise;
		map<string, double> _drag_noise;
		map<string, double> _emp_noise;

		map<string, double> _eclipse_ratio;
        map<string, double> _maneuver_inflate_pv;
        map<string, double> _maneuver_inflate_srp;

	private:
		// defalut settings, list here for the convenience of changing
		const int default_inp_setp = 60;
		const int default_out_step = 300;
		const int default_integra_degree = 11;
		const int default_rt_update_interval = 86400;
		const int default_init_time = 0;

		const double default_pv_nosie = 1.0E-11;
		const double default_srp_nosie = 2.0E-7;
		const double default_drag_noise = 1.0E-5;
		const double default_emp_noise = 1.0E-6;

		const double default_eclipse_ratio = 4.0;

        const double default_pv_maneuver_inflate = 1.0E-8;
        const double default_srp_maneuver_inflate = 1.0E-8;
	

		const string default_model_GRAVITY    = "EIGEN6S4_12";
		const string default_model_SOLAR_RADI = "ECOM2";
		const string default_model_RELATIVITY = "Default";
		const string default_model_POINT_MASS = "Default";
		const string default_model_ANTENNATHRUST = "Default";
		
		
		int _flag;								///< flag for updata hole force model info

		/** @brief updata info if xml defined. */
		void _updataForceWithXml();

		/**
		* @brief updata all GNSS sats settings
		*
		* reset force models for all GNSS satellite according to the GNS commom settings in config file
		*
		* @param[in]  spacecraft  set node in XML file
		*/
		void _resetAll(xml_node spcecraft);

		/**
		* @brief updata one system sats settings
		*
		* reset force models for a single system satellite according to the system commom settings in config file
		*
		* @param[in]  spacecraft  set node in XML file
		*/
		void _resetSys(GSYS sys, xml_node spacecraft);	

		/** 
		* @brief updata one sat settings
		*
		* reset force models for one satellite according to the sat commom settings in config file
		*
		* @param[in]  spacecraft  set node in XML file
		*/
		void _resetSat(string sat, xml_node spacecraft);
	};
}
#endif // !GSETFORCE_H
