/**
*
* @verbatim
	History
	 -1.0 Jiande Huang   2019-01-04  creat the file
	 -1.1 Jie Li		 2019-04-01  add code in doxygen style.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    goi_nav.cpp
* @brief		orbit integration processing for GNS satellite
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date		    2019-01-04
*
*/



#include "gproc/goi_nav.h"
#include "gforce/gforce.h"
#include "gforce/gfcommon.h"
#include "gset/gsetforce.h"
#include "gutils/ginfolog.h"
#include "newmat/newmatio.h"

namespace great
{

	/** @brief  constructor.
	* @param[in]	mark   satellite prn
	* @param[in]	set    setting information
	*/
	t_goi_nav::t_goi_nav(string mark, t_gsetbase* set)
	:t_goi(mark,set)
	{
		gtrace("t_goi_nav::t_goi_nav");

		_sat = mark;
		_set = set;

	}

	t_goi_nav::~t_goi_nav()
	{

	}

	/**
	* @brief set  ICS data
	* @param[in] gics    ICS data
	*/
	void t_goi_nav::set_data_used(t_gallics_nav* gics)
	{
		_gics = gics;
	}

	/**
	* @brief   init the sat infomation
	* @param[in] sat_name   satellite name
	*/
	bool t_goi_nav::_init_sat_info(string sat_name)
	{
		string function_name = "t_goi_nav::_init_sat_info";
		try
		{
			// 
			t_gtime beg = _extend_beg + NEXT_POINT_NUM * _step_out;
			t_gtime end = _extend_end - NEXT_POINT_NUM * _step_out;
			if (!_gallsat->getSatinfoUsed(sat_name, beg, end, _gsat))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the sat information of "
					+ sat_name + ", check your sat parameter file or the read process of the file.");
				return false;
			}
			_gnsmass = _gsat.mass();
			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 3, function_name, "ERROR : throw while getting the sat information");
			return false;
		}
	}

	/**
	* @brief      init the state matrix
	* @param[in]  neq_size   number of OI funct params
	* @param[out] state      state data of OI processing
	*/
	bool t_goi_nav::_init_state_matrix(int neq_size, vector<double>& state)
	{
		string function_name = "t_goi_nav::_init_state_matrix";
		try
		{
			// x,y,z,vx,vy,vz and solar pars will be read from ICS
			// 
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PXSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PYSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PZSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VXSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VYSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VZSAT));

			// initial partial can be ensured derectly
			// storaging with column vector is easy to integral
	
			for (int i = 1; i < neq_size / 6; ++i)
			{
				for (int j = 1; j < 7; j++)
				{
					if (i == j)
					{
						state.push_back(1.0);
					}
					else
					{
						state.push_back(0.0);
					}
				}
			}


			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 2, function_name, "ERROR : throw while init_state_matrix");
			return false;
		}
	}



	/**
	* @brief   Get settings from XML file and set local variables.
	* @return  get settings correctly or not
	*/
	bool t_goi_nav::_get_settings()
	{
		string function_name = "t_goi_nav::_get_settings";
		try
		{
			// sys and sat name
			_sys = t_gsys::str2gsys(_sat.substr(0, 1));

			// settings for intergrator
			_step_inp = dynamic_cast<t_gsetforce*>(_set)->getInpStep(_sat);
			_step_out = dynamic_cast<t_gsetforce*>(_set)->getOutStep(_sat);
			_integra_degree = dynamic_cast<t_gsetforce*>(_set)->getIntegraDegree(_sat);

			// settings for forces
			// List of models for considered forces
			_forces = dynamic_cast<t_gsetforce*>(_set)->getForceUsedSet(_sat);
			_force_model = dynamic_cast<t_gsetforce*>(_set)->getForceModelMap(_sat);
			_force_valid = dynamic_cast<t_gsetforce*>(_set)->getForcesValidMap(_sat);

			_mean_pole_model = dynamic_cast<t_gsetforce*>(_set)->getMeanPoleModel(_sat);
			
			bool temp1;
			int	 temp2;
			dynamic_cast<t_gsetforce*>(_set)->getGravitySet(_sat, _ocean_degree, _ot_mode, temp1, temp2);
			

			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 0, function_name, "ERROR : throw exception while get settings form xml file");
			return false;
		}

	}

	bool t_goi_nav::_check_data()
	{
		if (!_gplanet || _gplanet->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no de data");
			return false;
		};			///< data for all planets.

		if (!_gerp || _gerp->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no erp data");
			return false;
		};			///< erp data.

		if (!_gotdata || _gotdata->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no ocean tide data");
			return false;
		};			///< ocean tide data.
		
		if (!_gleapsecond || _gleapsecond->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no leapsecond data");
			return false;
		};		///< leap second data.
		
		if (!_gics || _gics->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no ics data");
			return false;
		};				///< all  ics data.
		
		if (!_gegmcoef || _gegmcoef->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no EGM data");
			return false;
		};			///< EGM data.
		
		if (!_gallsat || _gallsat->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no sat parameters data");
			return false;
		};			///< satellite parameters data.
		
		if (!_gotl)
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no ocean load data");
			return false;
		};				///< ocean load data.

		return true;
	}

	
	/**
	* @brief  calculate the acc and partial in epoch t
	* @param[in]  t			current epoch
	* @param[in]  state     state at current epoch
	* @param[out] partial   partial of the state
	* @return  calculate correctly or not
	*/
	bool t_goi_nav::_calculate(const t_gtime& t, const vector<double>& state, vector<double>& partial)
	{
		//change time sys to tdt
		t_gtime tdt_gtime = t;
		tdt_gtime.tsys(t_gtime::TT);
	
		double tdt = tdt_gtime.dmjd();

		//double tdt;
		//change_time_sys(GPS2TDT, t, tdt);
		if (tdt_gtime != _trs2crs_cord->getCurtEpoch()) {
			_trs2crs_cord->calcRotMat(tdt_gtime, true, true, true);
		}

		// Matrix
		ColumnVector satXJ(0.0, 3);		// satellite position in J2000 
		ColumnVector satVJ(0.0, 3);		// satellite VELOCITY in J2000
		ColumnVector satXE(0.0, 3);		// satellite position in TRS   

		map<string, ColumnVector> allXE;	
		map<string, ColumnVector> allXJ;	

		Matrix rot_scf2crs(0.0, 3, 3);	
		Matrix rot_crs2acr(0.0, 3, 3);	
		Matrix rot_trs2crs(0.0, 3, 3);	

		
		// get information of coordinate transformation from TRS to CRS like rotation matrix of TRS to CRS 
		Matrix _Dgmat = _trs2crs_cord->getMatDu();
		Matrix _Dxmat = _trs2crs_cord->getMatDx();
		Matrix _Dymat = _trs2crs_cord->getMatDy();
		rot_trs2crs = _trs2crs_cord->getRotMat();
		double gst = _trs2crs_cord->getGmst();
		double xpole = _trs2crs_cord->getXpole();
		double ypole = _trs2crs_cord->getYpole();


		// for planet
		ColumnVector sun_pos_crs(0.0, 3);
		ColumnVector sun_vel_crs(0.0, 3);
		ColumnVector moon_pos_crs(0.0, 3);
		ColumnVector moon_vel_crs(0.0, 3);

		for (int i = 0; i < 3; i++)
		{
			satXJ(i + 1) = state[i];
			satVJ(i + 1) = state[(int64_t)i + 3];
		}


		ColumnVector x_earth = (satXJ.t() * _Dgmat).t();
		satXE = (satXJ.t() * rot_trs2crs).t();

		// get the sun and moon pos in crs
		_gplanet->get_pos(tdt, "SUN", sun_pos_crs);
		_gplanet->get_pos(tdt, "MOON", moon_pos_crs);
		allXE["SUN"] = (sun_pos_crs.t() * rot_trs2crs).t();
		allXJ["SUN"] = sun_pos_crs;
		allXE["MOON"] = (moon_pos_crs.t() * rot_trs2crs).t();;

		double shadow_factor_ecom = get_shadow_factor(satXJ, sun_pos_crs, moon_pos_crs);

		_scf2crs_cord->calcRotMat(_gsat.blocktype(),_sat,satXJ, satVJ, allXJ["SUN"], rot_scf2crs);

		// get crs2acr
		_crs2acr_cord->calcRotMat(satXJ, satVJ, rot_crs2acr);

		// calculate acc
		t_gfcalarg_gravity_egm* arg_gravity = new t_gfcalarg_gravity_egm();
		arg_gravity->arg_type = ARG_GRAVI_EGM;
		arg_gravity->xpole = &xpole;
		arg_gravity->ypole = &ypole;
		arg_gravity->sat_pos_trs = &satXE;
		arg_gravity->e2j = &rot_trs2crs;
		arg_gravity->time = &tdt;
		arg_gravity->sat_pos = &satXJ;
		arg_gravity->sat_vel = &satVJ;
		arg_gravity->gst = &gst;
		arg_gravity->allXE = &allXE;
		arg_gravity->meanpole = _mean_pole_model;
		if (!_gforce.gravity->calculate(arg_gravity))
		{
			writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate gravity");
			delete arg_gravity; arg_gravity = NULL;
			return false;
		};
		delete arg_gravity; arg_gravity = NULL;


		t_gfcalarg_point_mass* arg_point_mass = new t_gfcalarg_point_mass();
		arg_point_mass->arg_type = ARG_POINT_MASS;
		arg_point_mass->time = &tdt;
		arg_point_mass->sat_pos = &satXJ;
		arg_point_mass->sat_vel = &satVJ;
		if (!_gforce.point_mass->calculate(arg_point_mass))
		{
			writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate point mass");
			delete arg_point_mass; arg_point_mass = NULL;
			return false;
		};
		delete arg_point_mass; arg_point_mass = NULL;


		t_gfcalarg_relativity* arg_relativity = new t_gfcalarg_relativity();
		arg_relativity->arg_type = ARG_RELATIVITY;
		arg_relativity->time = &tdt;
		arg_relativity->sat_pos = &satXJ;
		arg_relativity->sat_vel = &satVJ;
		if (!_gforce.relativity->calculate(arg_relativity))
		{
			writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate relativity");
			delete arg_relativity; arg_relativity = NULL;
			return false;
		};
		delete arg_relativity; arg_relativity = NULL;

		if (_gforce.solar_radi->force_model() == F_SOLAR_ECOM2)
		{
			t_gfcalarg_solar_ecom2* arg_solar_ecom2 = new t_gfcalarg_solar_ecom2();
			arg_solar_ecom2->arg_type = ARG_SOLAR_ECOM2;
			arg_solar_ecom2->sat_type = _gsat.blocktype();
			arg_solar_ecom2->moon_pos = &moon_pos_crs;
			arg_solar_ecom2->sun_pos = &sun_pos_crs;
			arg_solar_ecom2->shadow_factor = &shadow_factor_ecom;
			arg_solar_ecom2->time = &tdt;
			arg_solar_ecom2->sat_pos = &satXJ;
			arg_solar_ecom2->sat_vel = &satVJ;
			if (!_gforce.solar_radi->calculate(arg_solar_ecom2))
			{
				writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate solar_radi");
				delete arg_solar_ecom2; arg_solar_ecom2 = NULL;
				return false;
			};
			delete arg_solar_ecom2; arg_solar_ecom2 = NULL;
		}
		ColumnVector _acc_antennathrust(0.0, 3);
		if (_force_valid[ANTENNATHRUST] == true)
		{
			t_gfcalarg_antennathrust* arg_antenna = new t_gfcalarg_antennathrust();
			arg_antenna->arg_type = ARG_ANTENNATHRUST;
			arg_antenna->sat_mass = _gsat.mass();
			arg_antenna->sat_blocknum = _gsat.blocktype();
			arg_antenna->sat_pos = &satXJ;
			arg_antenna->sat_power = _gsat.power(); 
			if (!_gforce.antenna->calculate(arg_antenna))
			{
				writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate antenna thrust");
				delete arg_antenna; arg_antenna = NULL;
				return false;
			}
			_acc_antennathrust = _gforce.antenna->acc() * 1e-3;
			delete arg_antenna; arg_antenna = NULL;
		}

		//todo &
		ColumnVector _acc_point_mass = _gforce.point_mass->acc();
		ColumnVector _acc_gravity = _gforce.gravity->acc();
		ColumnVector _acc_relativity = _gforce.relativity->acc();
		ColumnVector _acc_solar_radi = _gforce.solar_radi->acc();

		ColumnVector _acc_all =
			(_acc_point_mass
				+ _acc_gravity
				+ _acc_relativity
				+ _acc_solar_radi
                +_acc_antennathrust);

		// calculate partial
		Matrix _partial_point_mass = _gforce.point_mass->partial();
		Matrix _partial_gravity = _gforce.gravity->partial();
		Matrix _partial_solar_radi = _gforce.solar_radi->partial();

		Matrix _partial_pos = _partial_gravity + _partial_point_mass;
		Matrix _partial_force = _partial_solar_radi;
		Matrix _partial_vel(0.0, 3, 3);

		// update partial 
		int index = 0;
		Matrix mat(0.0, 3, 1);
		mat = 0.0;
		vector<double> partial_tmp;
		for (int i = 1; i < 7; i++)
		{
			ColumnVector tmp(0.0, 3);
			index = 6 * i;
			tmp << state[(int64_t)index] << state[(int64_t)index + 1] << state[(int64_t)index + 2];
			mat = _partial_pos * tmp;
			tmp << state[(int64_t)index + 3] << state[(int64_t)index + 4] << state[(int64_t)index + 5];
			mat = mat + _partial_vel * tmp;
			for (int j = 1; j < 4; j++)
			{
				partial_tmp.push_back(mat(j, 1));
			}
		}

		mat = 0.0;
		for (int i = 1; i < _neq_size / 6 - 6; i++)
		{
			ColumnVector tmp(0.0, 3);
			index = 6 * i + 36;
			mat = _partial_force.Row(i).as_column();
			tmp << state[(int64_t)index] << state[(int64_t)index + 1] << state[(int64_t)index + 2];
			mat = mat + _partial_pos * tmp;
			tmp << state[(int64_t)index + 3] << state[(int64_t)index + 4] << state[(int64_t)index + 5];
			mat = mat + _partial_vel * tmp;
			for (int j = 1; j < 4; j++)
			{
				partial_tmp.push_back(mat(j, 1));
			}
		}
		//velocity
		partial[0] = state[3];
		partial[1] = state[4];
		partial[2] = state[5];

		//acc
		for (int i = 1; i <= _acc_all.size(); ++i)
		{
			partial[i + 2] = _acc_all(i);
		}
		//
		index = _acc_all.size() + 3;
		for (int i = 1; i < _neq_size / 6; i++)
		{
			partial[index] = state[(uint64_t)i * 6 + 3];
			partial[(uint64_t)index + 1] = state[(uint64_t)i * 6 + 4];
			partial[(uint64_t)index + 2] = state[(uint64_t)i * 6 + 5];

			int j = (i - 1) * 3;
			partial[(uint64_t)index + 3] = partial_tmp[j];
			partial[(uint64_t)index + 4] = partial_tmp[(uint64_t)j + 1];
			partial[(uint64_t)index + 5] = partial_tmp[(uint64_t)j + 2];

			index += 6;
		}
		return true;
	}

	/**
	* @brief  creat initial information of force model
	* @param[in]  force_ID			force type
	* @param[in]  force_contain     force model initial information
	* @return  creat force model correctly or not
	*/
	bool t_goi_nav::_creat_force_model(FORCE force_ID, t_gforce_nav * force_contain)
	{
		string function_name = "t_goi_nav::_creat_force_model";
		try
		{
			switch (force_ID)
			{
				// for gravity
			case great::GRAVITY:
			{
				t_init_gravity_egm08* gravity = NULL;


				if (_force_model[force_ID].find("EGM08_12") != string::npos)
				{
					gravity = new t_init_gravity_egm08(F_GRAVI_EGM08, 12, _gegmcoef, _gplanet, _gotl, _gotdata,_ot_mode,_mean_pole_model);
					force_contain->creatForceModel(gravity, _glog);
				}
				else if (_force_model[force_ID].find("EGM08_150") != string::npos)
				{
					gravity = new t_init_gravity_egm08(F_GRAVI_EGM08, 150, _gegmcoef, _gplanet, _gotl, _gotdata,_ot_mode,_mean_pole_model);
					force_contain->creatForceModel(gravity, _glog);
				}
				else if (_force_model[force_ID].find("EIGEN6S4") != string::npos)
				{
					int ipt1 = _force_model[force_ID].find("_");
					int ipt2 = _force_model[force_ID].find("+");
					int gravity_degrees;
					if (ipt2 == -1)
					{
						gravity_degrees = atoi(_force_model[force_ID].substr(ipt1 + 1).c_str());
					}
					else
					{
						gravity_degrees = atoi(_force_model[force_ID].substr(ipt1 + 1, ipt2 - ipt1).c_str());
					}
					gravity = new t_init_gravity_egm08(F_GRAVI_EIGEN6S4, gravity_degrees, _gegmcoef, _gplanet, _gotl, _gotdata, _ot_mode, _mean_pole_model);
					force_contain->creatForceModel(gravity, _glog);
				}
				else
				{
					writeLogInfo(_glog, 1, "WARNING : no setting for gravity model, we used EGM08 08th-order as default");
					gravity = new t_init_gravity_egm08(F_GRAVI_DEF, 8, _gegmcoef, _gplanet, _gotl, _gotdata);
					force_contain->creatForceModel(gravity, _glog);
				}
				delete gravity;
				gravity = NULL;
				break;
			}
			// for relavity
			case great::RELATIVITY:
			{
				t_init_relativity* relativity = NULL;
				if (_force_model[force_ID].find("Default") != string::npos)
				{
					relativity = new t_init_relativity(F_RELAV_DEF);
					force_contain->creatForceModel(relativity, _glog);
				}
				else
				{
					writeLogInfo(_glog, 1, "WARNING : no setting for relativity model, we used default model");
					relativity = new t_init_relativity(F_RELAV_DEF);
					force_contain->creatForceModel(relativity, _glog);
				}
				delete relativity;
				relativity = NULL;
				break;
			}
			case great::POINT_MASS:
			{
				t_init_point_mass* point_mass = NULL;
				if (_force_model[force_ID].find("Default") != string::npos)
				{
					point_mass = new t_init_point_mass(F_NBODY_DEF, _gplanet);
					force_contain->creatForceModel(point_mass, _glog);
				}
				else
				{
					writeLogInfo(_glog, 1, "WARNING : no setting for point_mass model, we used default model");
					point_mass = new t_init_point_mass(F_NBODY_DEF, _gplanet);
					force_contain->creatForceModel(point_mass, _glog);
				}
				delete point_mass;
				point_mass = NULL;
				break;
			}
			case great::SOLAR_RADI:
			{
				t_init_solar_ecom2* ecom2 = NULL;
				ecom2 = new t_init_solar_ecom2(F_SOLAR_ECOM2, _gics->get_nav_ics(_sat)->get_ecom2_solar_pars_values());
				force_contain->creatForceModel(ecom2, _glog);
				delete ecom2;
				ecom2 = NULL;
				break;
			}
			case great::ANTENNATHRUST:
			{
				t_init_antennathrust* antenna = NULL;
				if (_force_valid[force_ID] == true)
				{
					if (_force_model[force_ID].find("Default") != string::npos)
					{
						antenna = new t_init_antennathrust(F_ANTENNATHRUST, "default");
						force_contain->creatForceModel(antenna, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for earth radiation model,we use Default model");
						antenna = new t_init_antennathrust(F_ANTENNATHRUST, "default");
						force_contain->creatForceModel(antenna, _glog);
					}
				}
				delete antenna;
				antenna = NULL;
				break;
			}
			case great::CUSTOMER:
				break;
			case great::ATMOSPHERE:
				break;
			default:
				break;
			}

			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 0, function_name, "ERROR : throw exception while creating force model by the setting.");
			return false;
		}
	}


	bool t_goi_nav::prepareOI()
	{
		gtrace("t_gpppflt::prepareOI");
		string function_name = "t_goi_nav::prepareOI";

		_glog->comment(1, function_name, _sat + " prepareOI : started !");
		try {
			if (!_check_data())
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : no compelete data, check your xml file.");
				return false;
			}
			// get all the setting information
			if (!_get_settings())
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the settings form xml, check your xml file.");
				return false;
			}
			// get the ics for the sat
			_sat_ics = _gics->get_nav_ics(_sat);
			if (!_sat_ics)
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the ics data of " + _sat);
				return false;
			}

			// prepare the trs2crs and crs2scf and so on;
			_trs2crs_cord = new t_gtrs2crs(false, _gerp);
			_scf2crs_cord = new t_gscf2crs();
			_crs2acr_cord = new t_gcrs2acr();

			// init all the force models
			// gravity model
			// 
			if (!_creat_force_model(GRAVITY, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(GRAVITY));
				return false;
			};
			// N-body model
			// 
			if (!_creat_force_model(POINT_MASS, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(POINT_MASS));
				return false;
			};
			// solar radiation model
			//
			if (!_creat_force_model(SOLAR_RADI, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(SOLAR_RADI));
				return false;
			};
			//antenna thrust model
			if (!_creat_force_model(ANTENNATHRUST, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(ANTENNATHRUST));
				return false;
			};
			//relativity
			//
			if (!_creat_force_model(RELATIVITY, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(RELATIVITY));
				return false;
			};

			// get the number of equation and init the partial
			//_par_num = _sat_ics->get_par_num();
			if (_sat_ics->get_par_num() == 0)
			{
				writeLogInfo(_glog, 1, function_name, "ERROR : no data in the " + _sat + " ics file, check your ics file.");
				return false;
			}
			else
			{
				_neq_size = 6 * (_sat_ics->get_par_num() + 1);
			}
			return true;
		}
		catch (...)
		{
			if (_trs2crs_cord)
			{
				delete _trs2crs_cord; _trs2crs_cord = NULL;
			}
			if (_scf2crs_cord)
			{
				delete _scf2crs_cord; _scf2crs_cord = NULL;
			}
			if (_crs2acr_cord)
			{
				delete _crs2acr_cord; _crs2acr_cord = NULL;
			}
			writeLogInfo(_glog, 1, function_name, "ERROR : throw exception while doing procseeBatch of" + _sat);
			return false;

		}
	}

	bool t_goi_nav::integrate4srif(const t_gtime & beg, const t_gtime & end, vector<t_gpar::t_type> state_init_name, vector<double> state_init_value, map<int, Matrix>& phi_map, map<int, Matrix>& phi_i_map, map<int, vector<t_gtriple>>& result, map<t_gtime, ColumnVector>& state_map)
	{
		gtrace("t_goi_nav::integrate4srif");
		string function_name = "t_goi_nav::integrate4srif";

		//integrate to get the state of new reference epoch
		// here set the value of sat_icc according to init state matrix
		_sat_ics = _gics->get_nav_ics(_sat);
		_sat_ics->clear();
		_gics->beg_time = beg;
		_gics->end_time = end;
		for (int i = 0; i < state_init_name.size(); i++) {
			t_ics_record ics_record;
			ics_record.enum_name = t_gics_nav::_gpar2icspar(state_init_name[i]);
			ics_record.str_name = ics2str(ics_record.enum_name);
			ics_record.value = state_init_value[i];
			ics_record.sigma = 0.0; 
			ics_record.beg_time = beg;
			ics_record.end_time = end;
			_sat_ics->set_par_record(ics_record.str_name, ics_record);
		}
		if (!_creat_force_model(SOLAR_RADI, &_gforce))
		{
			writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(SOLAR_RADI));
			return false;
		}

		// set time info
		_beg = beg;
		_end = end;

		// init sat info
		_extend_beg = beg - NEXT_POINT_NUM * _step_out;
		_extend_end = end + NEXT_POINT_NUM * _step_out;
		_init_sat_info(_sat);
		if (_gsat.prn().empty())
		{
			writeLogInfo(_glog, 1, function_name, "ERROR : no information for " + _sat + ", check your sat_paramters file.");
			return false;
		}

		// init state value 
		vector<double> state_init_vec;
		if (!_init_state_matrix(_neq_size, state_init_vec))
		{
			writeLogInfo(_glog, 0, function_name, "ERROR : can not init the state matrix.");
			return false;
		};
		
		// get beta (elevation of the sun above the orbit plane) at the begin time

		_beta_beg = _cal_beta(_beg,state_init_vec);

		// forward integration
		_numerical_integration_for_srif(beg, _extend_beg, state_init_vec, phi_map, phi_i_map, result, state_map);
		phi_map[0] = IdentityMatrix(state_init_name.size());
		phi_i_map[0] = IdentityMatrix(state_init_name.size());
		vector<t_gtriple> orb_record;
		// storage orbit state on current time 
		ColumnVector state_curr(0.0, state_init_vec.size() / 6 - 1);
		// zero time initialization
		for (int i = 0; i < state_init_vec.size() / 6; i++) {
			orb_record.push_back(t_gtriple(state_init_vec[i * 6], state_init_vec[i * 6 + 1], state_init_vec[i * 6 + 2]));
			if (i < state_init_vec.size() / 6 - 1) {
				state_curr(i + 1) = state_init_value[i];
			}
		}
		result[0] = orb_record;
		state_map[beg] = state_curr;

		// backward intergration
		_numerical_integration_for_srif(beg, _extend_end, state_init_vec, phi_map, phi_i_map, result, state_map);

		// Uniformly initialize the parameter values of the mechanical model at all times
		for (auto iter = state_map.begin(); iter != state_map.end(); iter++) {
			for (int i = 7; i < state_init_vec.size() / 6; i++) {
				iter->second(i) = state_init_value[i - 1];
			}
		}
		return true;
	}


	bool t_goi_nav::_numerical_integration_for_srif(t_gtime beg, t_gtime end, const vector<double>& state_init, map<int, Matrix>& phi_map, map<int, Matrix>& phi_i_map, map<int, vector<t_gtriple>>& result, map<t_gtime, ColumnVector>& state_map)
	{
		string function_name = "t_goi::_numerical_integration";
		// init rkf settings
		t_gtime beg_time = beg;
		t_gtime end_time = end;
		t_adams* adams = NULL;
		try
		{
			if (_step_inp <= 0 || _step_out <= 0)
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : invalid integrator step_in/step_out settings");
				return false;
			}
			int step_inp = _step_inp;
			if (_step_out < step_inp || _step_out % step_inp != 0)
			{
				// Keep the Adams step a divisor of step_out so phi_map is populated at every output epoch.
				int d = step_inp < _step_out ? step_inp : _step_out;
				while (d > 1 && _step_out % d != 0) --d;
				writeLogInfo(_glog, 1, function_name, "Warning : step_in does not divide step_out; adjusting step_in to " + to_string(d));
				step_inp = d;
			}

			bool forward = (beg >= end ? true : false);
			int  sign = (forward ? 1 : -1);
			double rkf_h = step_inp * sign / 30.0;

			if (rkf_h == 0.0)
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : the _step in was zero, check your process or xml file settings");
				return false;
			}

			// init adams 
			vector<list<double> > admas_state;
			vector<list<double> > admas_partial;

			// get the adams' state and partial by using rkf
			if (!_prepare_adams_with_rkf(_integra_degree, beg_time, rkf_h, state_init, admas_state, admas_partial))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the adams points with the rkf");
				return false;
			}


			int num_of_all = 0;

			// start adams
			int adams_sign = (forward ? -1 : 1);
			double adams_step = step_inp * 1.0 * adams_sign;
			t_gtime adams_ref_time = beg;
			adams = new t_adams(_integra_degree, adams_step);

			// try at first to get the real epoch num;
			int epoch_num = 0;
			bool isValid = true;
			t_gtime epoch_try = beg;
			int n_inter_point = _step_out / step_inp;

			while (epoch_try.dmjd()* adams_sign < end.dmjd() * adams_sign)
			{
				epoch_num++;
				epoch_try = beg + epoch_num * adams_step;

				if (!_processEpoch(epoch_try, adams, admas_state, admas_partial))
				{
					isValid = false;
					break;
				}
				// here means that reach the output point
				if (epoch_num % n_inter_point == 0) {
					int epo = epoch_num * adams_sign / n_inter_point;
					int norb = _neq_size / 6 - 1;
					// phi storage state transition matrix
					Matrix phi(norb, norb); phi = 0.0;
					// orb_record storage the results format same as the orb file
					vector<t_gtriple> orb_record;
					orb_record.push_back(t_gtriple(admas_state[0].back(), admas_state[1].back(), admas_state[2].back()));
					// state storage the orb state X Y Z VX VY VZ
					ColumnVector state(0.0, norb);
					for (int i = 1; i < _neq_size / 6; i++)
					{
						// phi
						phi(1, i) = admas_state[(int64_t)i * 6].back();
						phi(2, i) = admas_state[(int64_t)i * 6 + 1].back();
						phi(3, i) = admas_state[(int64_t)i * 6 + 2].back();
						phi(4, i) = admas_state[(int64_t)i * 6 + 3].back();
						phi(5, i) = admas_state[(int64_t)i * 6 + 4].back();
						phi(6, i) = admas_state[(int64_t)i * 6 + 5].back();

						if (i > 6) {
							phi(i, i) = 1;
						}
						else {
							state(i) = admas_state[i - 1].back();
						}
						orb_record.push_back(t_gtriple(admas_state[(int64_t)i * 6].back(),
							admas_state[(int64_t)i * 6 + 1].back(), admas_state[(int64_t)i * 6 + 2].back()));
					}
					phi_map[epo] = phi;
					phi_i_map[epo] = phi.i();
					result[epo] = orb_record;
					state_map[epoch_try] = state;
				}
			}

			delete adams;
			adams = NULL;
			if (isValid)
			{
				return true;
			}
			else
			{
				writeLogInfo(_glog, 0, function_name,
					"ERROR : There are some problems in process Epoch while calculating with adams");
				return false;
			}

		}
		catch (...)
		{
			if (adams)
			{
				delete adams;
				adams = NULL;
			}
			writeLogInfo(_glog, 3, function_name, "ERROR : throw exception while numerical_integration");
			return false;
		}
		return true;
	}


}