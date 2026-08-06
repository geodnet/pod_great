/**
*
* @verbatim
	History
	 -1.0 Jie Li		 2019-01-10 creat the file.
	 -1.1 Jie Li      2019-04-01 add code in doxygen style.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			getforce_nav.cpp
* @brief		force model set from XML
*
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date			2019-01-10
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetforce.h"
#include "gutils/gsys.h"
using namespace std;
using namespace pugi;


namespace great
{
	/** @brief constructor for navigation satellite. */
	t_gsetforce::t_gsetforce()
		: t_gsetbase()
	{
		_gmutex.lock();

		_set.insert(XMLKEY_FORCE);


		//init GNS setting
		t_map_sats gnss_sats = GNSS_SATS();

		for (auto itGNSS = gnss_sats.begin(); itGNSS != gnss_sats.end(); ++itGNSS)
		{
			GSYS gsys = itGNSS->first;
			set<string> sats = itGNSS->second;

			for (auto itSAT = sats.begin(); itSAT != sats.end(); ++itSAT)
			{
				string sat = *itSAT;
				_map_force_valid[sat][GRAVITY]    = true;
				_map_force_valid[sat][SOLAR_RADI] = true;
				_map_force_valid[sat][RELATIVITY] = true;
				_map_force_valid[sat][POINT_MASS] = true;
				_map_force_valid[sat][ANTENNATHRUST] = false;

				_map_force_model[sat][GRAVITY]    = default_model_GRAVITY   ;
				_map_force_model[sat][SOLAR_RADI] = default_model_SOLAR_RADI;
				_map_force_model[sat][RELATIVITY] = default_model_RELATIVITY;
				_map_force_model[sat][POINT_MASS] = default_model_POINT_MASS;
				_map_force_model[sat][ANTENNATHRUST] = default_model_ANTENNATHRUST;

				_inp_step[sat] = default_inp_setp;
				_out_step[sat] = default_out_step;
				_integra_degree[sat] = default_integra_degree;

				_rt_update_interval[sat] = default_rt_update_interval;
				_init_time[sat] = default_rt_update_interval;
				_pv_noise[sat] = default_pv_nosie;
				_srp_noise[sat] = default_srp_nosie;
				_drag_noise[sat] = default_drag_noise;
				_emp_noise[sat] = default_emp_noise;
				_eclipse_ratio[sat] = default_eclipse_ratio;
                _maneuver_inflate_pv[sat] = default_pv_maneuver_inflate;
                _maneuver_inflate_srp[sat] = default_srp_maneuver_inflate;

				_ot_degree[sat] = 30;
				_opt_valid[sat] = false;
				_opt_degree[sat] = 30;
				_ot_mode[sat] = "";
				_flag = 0;
			}
		}
		

		_gmutex.unlock();
	}

	/**
	* @brief get used force for GNS satellite(valid)
	* @param[in]  sat   satellite name
	* @return  the forces used
	*/
	set<FORCE> t_gsetforce::getForceUsedSet(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}

		set<FORCE> force_used;
		map< FORCE, bool > force_tmp = _map_force_valid[sat];

		map< FORCE, bool >::iterator itUsed = force_tmp.begin();
		for (itUsed; itUsed != force_tmp.end(); ++itUsed)
		{
			if (itUsed->second)
			{
				force_used.insert(itUsed->first);
			}
		}
		return force_used;

		
	}

	/**
	* @brief get used force model for GNS satellite
	* @param[in]  sat  satellite name
	* @return  force model used
	*/
	map<FORCE, string> t_gsetforce::getForceModelMap(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _map_force_model[sat];
		
	}

	/**
	* @brief get all force models
	* @param[in]  sat              satellite name
	* @return    map of force model(true or false)
	*/
	map<FORCE, bool> t_gsetforce::getForcesValidMap(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		if (_map_force_valid.count(sat) == 0)
		{
			map<FORCE, bool> tmp_vaild;
			tmp_vaild[GRAVITY] = true;
			tmp_vaild[RELATIVITY] = true;
			tmp_vaild[CUSTOMER] = true;
			tmp_vaild[POINT_MASS] = true;
			tmp_vaild[SOLAR_RADI] = true;
			tmp_vaild[ATMOSPHERE] = true;
			tmp_vaild[ANTENNATHRUST] = false;
			return tmp_vaild;
		}
		else
		{
			return _map_force_valid[sat];
		}
		
		
	}

	modeofmeanpole t_gsetforce::getMeanPoleModel(string csat)
	{
        lock_guard<t_gmutex> lck(_gmutex);

		return _mean_pole_model[csat];
	}

	/**
	* @brief update all force models for GNS satellite
	*/
	void t_gsetforce::_updataForceWithXml()
	{
		//get force model node
		xml_node force_model = _doc.child(XMLKEY_ROOT).child(XMLKEY_FORCE);

		//reset forcemodel info by xml
		if (!force_model)
		{
			return;
		}

		map<int, set<string> > id_all;
		for (xml_node spacecraft = force_model.first_child(); spacecraft; spacecraft = spacecraft.next_sibling())
		{
			//get force model ID and classificate
			string ID = spacecraft.first_attribute().value();
			if (ID == "GNS"||ID=="LEO")
			{
				id_all[1].insert(ID);
			}
			else if (ID == "GPS" || ID == "GLO" || ID == "GAL" || ID == "BDS" || ID == "QZS")
			{
				id_all[2].insert(ID);
			}
			else
			{
				id_all[3].insert(ID);
			}
		}

		//from low-level sat group to sat sys to sat
		for (auto itID = id_all.begin(); itID != id_all.end(); ++itID)
		{
			set<string> tmp = itID->second;

			if (tmp.empty())
			{
				continue;
			}

			for (auto itVALUE = tmp.begin(); itVALUE != tmp.end(); ++itVALUE)
			{
				xml_node tmp_sc = force_model.find_child_by_attribute("spacecraft", "ID", (*itVALUE).c_str());
				int switch_case = itID->first;
				switch (switch_case)
				{
				case 1:
					//reset GNS forcemodel
					_resetAll(tmp_sc);
					break;
				case 2:
					//reset single sys forcemodel
					_resetSys(t_gsys::str2gsys(*itVALUE), tmp_sc);
					break;
				case 3:
					//reset single sat forcemodel
					_resetSat(*itVALUE, tmp_sc);
					break;
				default:
					break;
				}
			}
		}
	}

	/**
	* @brief get step-in of adams orbit integration(unit: s)
	* @param[in]  sat      satellite name
	* @return    OI input step
	*/
	int t_gsetforce::getInpStep(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _inp_step[sat];
	
	}

	/**
	* @brief get step-out of orbit integration(unit: s)
	* @param[in]  sat      satellite name
	* @return    OI output step
	*/
	int t_gsetforce::getOutStep(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _out_step[sat];
	
	}

	int t_gsetforce::getRtUpdateInt(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _rt_update_interval[sat];
	}

	int t_gsetforce::getInitTime(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _init_time[sat];
	}

	double t_gsetforce::getPvNoise(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _pv_noise[sat];
	}

	double t_gsetforce::getSrpNoise(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _srp_noise[sat];
	}

	double t_gsetforce::getDragNoise(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _drag_noise[sat];
	}

	double t_gsetforce::getEmpNoise(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _emp_noise[sat];
	}

	double t_gsetforce::getEclipseRatio(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _eclipse_ratio[sat];
	}


	/**
	* @brief get order of adams orbit integration
	* @param[in]  sat      satellite name
	* @return    OI adams order
	*/
	int t_gsetforce::getIntegraDegree(string sat)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _integra_degree[sat];
	
	}
	/**
	* @brief get gravity(ocean_tide and ocean_pole_tide) set 
	* @param[in]  sat                  satellite name
	* @param[out] ot_degree            ocean_tide degree
	* @param[out] _opt_valid           ocean_pole_tide valid
	* @param[out] _opt_degree          ocean_pole_tide degree
	* @retrun
	*/
	void t_gsetforce::getGravitySet(string sat, int& ot_degree, string& ot_mode, bool& opt_valid, int& opt_degree)
	{
		lock_guard<t_gmutex> lck(_gmutex);
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		ot_degree = _ot_degree[sat];
		opt_valid = _opt_valid[sat];
		opt_degree = _opt_degree[sat];
		ot_mode = _ot_mode[sat];
	}

	

	/**
	* @brief reset force models for all GNSS satellite according to the GNS commom settings in config file
	* @param[in]  spacecraft  set node in XML file
	*/
	void t_gsetforce::_resetAll(xml_node spacecraft)
	{
		t_map_sats gnss_data = GNSS_SATS();
		for (auto itGNSS = gnss_data.begin(); itGNSS != gnss_data.end(); ++itGNSS)
		{
			
			_resetSys(itGNSS->first, spacecraft);
		}
	}

	/**
	* @brief reset force models for a single system satellite according to the system commom settings in config file
	* @param[in]  spacecraft  set node in XML file
	*/
	void t_gsetforce::_resetSys(GSYS sys, xml_node spacecraft)
	{
		t_map_sats gnss_data = GNSS_SATS();
		for (auto itSAT = gnss_data[sys].begin(); itSAT != gnss_data[sys].end(); ++itSAT)
		{
			_resetSat(*itSAT, spacecraft);
		}
	}

	/**
	* @brief reset force models for one satellite according to the sat commom settings in config file
	* @param[in]  spacecraft  set node in XML file
	*/
	void t_gsetforce::_resetSat(string sat, xml_node spacecraft)
	{
		//GSYS sys = str2gsys(sat.substr(0, 1));
		for (xml_node tmp_set = spacecraft.first_child(); tmp_set; tmp_set = tmp_set.next_sibling())
		{

			if ((strcmp("integrator", tmp_set.name())) == 0)
			{
				_inp_step[sat] = str2int(tmp_set.attribute("step_in").value());
				_out_step[sat] = str2int(tmp_set.attribute("step_out").value());
				_integra_degree[sat] = str2int(tmp_set.attribute("order").value());
				_rt_update_interval[sat] = tmp_set.attribute("rt_update_interval").as_int(_rt_update_interval[sat]);
				_init_time[sat] = tmp_set.attribute("init_time").as_int(_init_time[sat]);
			}

			if ((strcmp("rt_update_noise", tmp_set.name())) == 0)
			{
				_pv_noise[sat] = tmp_set.attribute("pv").as_double(default_pv_nosie);
				_srp_noise[sat] = tmp_set.attribute("srp").as_double(default_srp_nosie);
				_drag_noise[sat] = tmp_set.attribute("drag").as_double(default_drag_noise);
				_emp_noise[sat] = tmp_set.attribute("emp").as_double(default_emp_noise);
			}

			if ((strcmp("eclipse_ratio", tmp_set.name())) == 0)
			{
				_eclipse_ratio[sat] = str2dbl(tmp_set.child_value());
			}

            if ((strcmp("maneuver_inflate", tmp_set.name())) == 0)
            {
                _maneuver_inflate_pv[sat] = tmp_set.attribute("pv").as_double(default_pv_maneuver_inflate);
                _maneuver_inflate_srp[sat] = tmp_set.attribute("srp").as_double(default_srp_maneuver_inflate);
            }

            if ((strcmp("gravity", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][GRAVITY] = str2bool(tmp_set.attribute("valid").value());
				_map_force_model[sat][GRAVITY] = tmp_set.attribute("model").value();
				for (xml_node child_set = tmp_set.first_child(); child_set; child_set = child_set.next_sibling())
				{
					if ((strcmp("ocean_tide", child_set.name())) == 0)
					{
						_ot_degree[sat] = str2int(child_set.attribute("degree").value());
						_ot_mode[sat] = child_set.attribute("model").value();
					}
					if ((strcmp("ocean_pole_tide", child_set.name())) == 0)
					{
						_opt_valid[sat] = str2bool(child_set.attribute("valid").value());
						_opt_degree[sat] = str2int(child_set.attribute("degree").value());
					}
					if ((strcmp("mean_pole", child_set.name())) == 0)
					{
						string mean_pole_model = child_set.attribute("model").value();
						if (mean_pole_model.find("CUBIC") != string::npos || mean_pole_model.find("cubic") != string::npos) {
							_mean_pole_model[sat] = modeofmeanpole::cubic;
						}
						else if (mean_pole_model.find("LINEAR") != string::npos || mean_pole_model.find("linear") != string::npos) {
							_mean_pole_model[sat] = modeofmeanpole::linear;
						}
						else {
							_mean_pole_model[sat] = modeofmeanpole::cubic;
						}
						
					}
				}
			}

			if ((strcmp("relativity", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][RELATIVITY] = str2bool(tmp_set.attribute("valid").value());
			}
			if ((strcmp("point_mass", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][POINT_MASS] = str2bool(tmp_set.attribute("valid").value());

				string word;
				//_map_point_mass_list[sys][sat].clear();
				istringstream is(tmp_set.first_child().value());
				while (is >> word)
				{
					//_map_point_mass_list[sys][sat].insert(word);
				}
			}
			if ((strcmp("antenna_thrust", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][ANTENNATHRUST] = str2bool(tmp_set.attribute("valid").value());
				_map_force_model[sat][ANTENNATHRUST] = tmp_set.attribute("model").value();
			}


		}
	}

	/**
	* @brief settings check
	*/
	void t_gsetforce::check()
	{
		_gmutex.lock();
		_gmutex.unlock();
		return;
	}

	/**
	* @brief settings help
	*/
	void t_gsetforce::help()
	{
		_gmutex.lock();
		
		cerr << "<force_model_nave>                                                       " << endl
			 << "	<spacecraft ID = \"GNS\">                                             " << endl
			 << "		<integrator step_in = \"60\" step_out = \"300\"  order = \"11\" />" << endl
			 << "		<gravity          valid = \"True\"  model = \"EGM_12\" />         " << endl
			 << "		<relativity       valid = \"True\"  model = \"Default\"/>         " << endl
			 << "		<point_mass       valid = \"True\"  model = \"Default\"/>         " << endl
			 << "		<solar_radiation  valid = \"True\"  model = \"ECOM_5\" />         " << endl
			 << "	</spacecraft>                                                         " << endl
			 << endl
			 << "<!--> if your have special setting for one system, you can set as following : <!-->" << endl
			 << "	<spacecraft ID = \"GPS\">                                             " << endl
			 << "		<integrator step_in = \"60\" step_out = \"300\"  order = \"11\" />" << endl
			 << "		<gravity    valid = \"True\" model = \"EGM_12\" />                " << endl
			 << "		<relativity valid = \"True\" model = \"Default\"  />              " << endl
			 << "		<point_mass valid = \"True\" model = \"Default\"  />              " << endl
			 << "		<solar_radiation  valid = \"True\"  model = \"ECOM2\" />         " << endl
			 << "	</spacecraft>                                                         " << endl 
			 << endl
			 << "<!--> if your have special setting for one satellite, you can set as following : <!-->" << endl
			 << "	<spacecraft ID = \"G01\">                                             " << endl
			 << "		<integrator step_in = \"60\" step_out = \"300\"  order = \"11\" />" << endl
			 << "		<gravity    valid = \"True\" model = \"EGM_12\" />                " << endl
			 << "		<relativity valid = \"True\" model = \"Default\"  />              " << endl
			 << "		<point_mass valid = \"True\" model = \"Default\"  />              " << endl
			 << "		<solar_radiation  valid = \"True\"  model = \"ECOM_5\" />         " << endl
			 << "	</spacecraft>  		                                                  " << endl	
			 << endl;
		_gmutex.unlock();
		return;
	}

}