/**
*
* @verbatim
History
-1.0 glfeng  2019-03-05  creat the file.
-1.2 BoWong  2019-04-09  Adding Doxygen Style Code Remarks.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambfixDd.cpp
* @brief	integer estimation with the LAMBDA method.
*
* @author   glfeng BoWong , Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/

#include "gambfix/gambfixDd.h"
#include "gambfix/gallamb.h"
#include "gset/gsetout.h"
#include "gio/gfile.h"

#include "math.h"
#include <string>
#include <stdio.h>
#include <chrono>
using namespace std;

namespace great
{
	t_gambfixDd::t_gambfixDd()
	{

	}

	t_gambfixDd::t_gambfixDd(t_glog* log, t_gsetbase* set) :
		t_gambfix(log, set)
	{
		cerr << "================> constructor of t_gambfixDd" << endl;


		// add data to gall ambs
		this->_gambs.add_set(_gset);
		this->_gambs.add_log(_glog);
	}


	t_gambfixDd::~t_gambfixDd()
	{

	}

	bool t_gambfixDd::ProcessBatch(t_gtime beg, t_gtime end)
	{
		chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
		string func_id = "ProcessBatch";
		
		if (!this->_gpar || !this->_gset || !this->_gbias)
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, this->_class_id, func_id, "check the gpar || gset || gbias");
			return false;
		}

		// clear all ambs data before combine ambs
		this->_gambs.clear_ambs_data();

		// add data to gambs
		
		this->_gambs.add_data(_gpar, _gobs, _gbias, _gobj,beg,end);
		this->_gambs.add_glo_fid(_glofrq_num);

		// get baselines for data or setting
		t_gambfixDd::_screen_rec_and_sat(this->_gambs.IDs(), _gpar, _gset, _gobj);
		this->_gambs.add_WL_bdeci(this->_evaluator_WL);
		this->_gambs.add_NL_bdeci(this->_evaluator_NL);

		chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
		
		DD_MODEL dd_mode = dynamic_cast<t_gsetamb*>(_gset)->dd_mode();
		for (const auto& ID : this->_gambs.IDs())
		{

			beg_t = chrono::high_resolution_clock::now();
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "get Sd ambs");
			this->_gambs.cmb_amb(t_gallambs::AMB_COMB_OPTION::CMB_SAT_SD, ID, this->_gSatBl[ID].baselines(), this->_min_common_time);
			end_t = chrono::high_resolution_clock::now();

			
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "get Dd ambs");
			this->_gambs.cmb_amb(t_gallambs::AMB_COMB_OPTION::CMB_REC_DD, ID, this->_gRecBl[ID].baselines(false), this->_min_common_time);
		}
		
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "check Dd ambs");
		
		this->_write_ambcon(false);

		return true;
	}

	void t_gambfixDd::update_data(t_gallpar* pars, t_gallobs* obs, t_gallbias* bias, t_gallobj* obj, t_gtime beg_time, t_gtime end_time)
	{
		this->_gambs.clear_ambs_data();
		this->_gambs.add_data(pars, obs, bias, obj, beg_time, end_time);
	}



	bool t_gambfixDd::_screen_rec_and_sat(const set<t_gamb::AMB_ID>& IDs, t_gallpar* pars, t_gsetbase* set, t_gallobj* obj)
	{
		if (!pars && !set && !obj) return false;

		_gRecBl.clear();
		_gSatBl.clear();
		_rec_list.clear();
		_sat_list.clear();

		auto xml_rec_list = dynamic_cast<t_gsetgen*>(set)->recs();
		auto xml_sat_list = dynamic_cast<t_gsetgnss*>(set)->sat();
		for (const auto& ID : IDs)
		{
		
			try
			{
				std::set<string> par_rec_list;
				std::set<string> par_sat_list;

				
				_get_rec_and_sat(ID, pars, par_rec_list, par_sat_list);
				set_intersection(par_rec_list.begin(), par_rec_list.end(), xml_rec_list.begin(), xml_rec_list.end(), inserter(_rec_list, _rec_list.begin()));
				set_intersection(par_sat_list.begin(), par_sat_list.end(), xml_sat_list.begin(), xml_sat_list.end(), inserter(_sat_list, _sat_list.begin()));
			}
			catch (...)
			{
				return false;
			}

			const auto npars = pars->parNumber();
			for (unsigned int i = 0; i < npars; i++)
			{
				auto rec = (*pars)[i].site;
				auto sat = (*pars)[i].prn;
				if (_rec_list.find(rec) == _rec_list.end() || _sat_list.find(sat) == _sat_list.end()) continue;
				_gRecBl[ID].add_rec(rec, sat);
				_gSatBl[ID].add_sat(sat);
			}

		
			if (_rec_list.empty() || _sat_list.empty())
			{
				return false;
			}

			_gRecBl[ID].add_cut(_beg, _max_baseline_length);
			_gRecBl[ID].add_gobj(obj);
			_gRecBl[ID].add_gset(set);
			_gRecBl[ID].add_glog(this->_glog);
		}
		return true;

	}

	void t_gambfixDd::_write_ambcon(bool write_file)
	{
	
		t_gambcon _gambcon;
		auto ddambs = this->_gambs.get_Dds();
		
		int num_dds = 0;
		for (auto& sys : ddambs)
		{
			for (auto& bl : sys.second)
			{
				for (auto& amb : bl.second)
				{
					num_dds++;
					_gambcon.add_ddamb(amb);
				}
			}
		}

		_gambcon.set_header(_beg, _intv, num_dds, _sat_list, _rec_list);
		_gambcon.all_NL_num = this->_gambs.all_NL_num;
		_gambcon.all_WL_num = this->_gambs.all_WL_num;
		_gambcon.fix_NL_num = this->_gambs.fix_NL_num;
		_gambcon.fix_WL_num = this->_gambs.fix_WL_num;

		_ambcon_results = _gambcon;
		if (!write_file) {
			return;
		}

	}

	bool t_gambfixDd::_get_rec_and_sat(const t_gamb::AMB_ID& ID, t_gallpar* pars, set<string>& par_rec_list, set<string>& par_sat_list)
	{
		t_gpar::t_type amb_type = t_gpar::t_type::NO_DEF;
		switch (ID)
		{
			case t_gamb::AMB_ID::IF_L1_L2:  {amb_type = t_gpar::t_type::AMB_IF  ; break;}
			default:break;
		}

		if (amb_type == t_gpar::t_type::NO_DEF)
		{
			par_rec_list = pars->amb_recs();
			par_sat_list = pars->amb_prns();
		}
		else
		{
			par_rec_list = pars->amb_recs(amb_type);
			par_sat_list = pars->amb_prns(amb_type);
		}

		return true;
	}
}



