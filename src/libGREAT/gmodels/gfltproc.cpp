#include "gfltproc.h"
/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-04-04  created
	 -1.1 glfeng   2019-09-05  add function _batch_turboedit and variable _turboedit_software
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gfltproc.cpp
* @brief		virtual class for process
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/

#include "gmodels/gfltproc.h"
#include "gset/gsetpar.h"
#include "gset/gsetamb.h"
#include "gset/gsetout.h"
#include "gset/gsetturboedit.h"
#include "gio/gfile.h"
#include "gutils/ginfolog.h"
#include "gutils/gcycleslip.h"
#include "gutils/gturboedit.h"
#include "gutils/gtime.h"
#include "gcoders/erp.h"
#include "gutils/gfileconv.h"
#include "gutils/gtime.h"
#include "gproc/gupdateparIF.h"
#include "gutils/ginfolog.h"
#include "gutils/gstring.h"
#include "gproc/goi_nav.h"
#include "gutils/gtiming.h"
#include "gset/gsetforce.h"
#include <algorithm>
#include <thread>



// Judge whether file exsit
#if  defined _WIN32 || defined _WIN64
#include <io.h>
#include <direct.h> 
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#include <unistd.h>
#include <sys/stat.h>
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif



namespace great
{
	t_gfltproc::t_gfltproc()
	{}

	t_gfltproc::t_gfltproc(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		_glog(log),
		_gset(set),
		_gall_proc(data),
		_recclk_threshold(1e-8)
	{
		// check log
		if (_glog == nullptr)
		{
			_myLog = _create_log(_glog);
		}

		
		// get settings from xml
		_init_xml_settings(set);
		_init_lsq_settings(set, data, log);

	
	}

	t_gfltproc::~t_gfltproc()
	{
		// Init Paramter in Constructor
		if (_lsq)
		{
			delete _lsq;
			_lsq = nullptr;
		}

		if (_myLog)
		{
			_delete_log(_glog);
		}
	}

	void t_gfltproc::add_coder(const vector<t_gcoder*>& coder)
	{
		_gcoder = coder;
	}

	bool t_gfltproc::ProcessBatch(t_gallproc* data, const t_gtime& beg, const t_gtime& end)
	{
		return true;
	}

	bool t_gfltproc::GenerateProduct()
	{
		return false;
	}

	bool t_gfltproc::_init_crs_erp_pars(t_glsq* lsq)
	{
		int ipar_num = _lsq->_x_solve.parNumber();

		// init EOP parmeter
		// init the Eop dx dy vdx vdy
		t_gpar par_xpole("", t_gpar::XPOLE, ++ipar_num, "");
		t_gpar par_ypole("", t_gpar::YPOLE, ++ipar_num, "");
		t_gpar par_dxpole("", t_gpar::DXPOLE, ++ipar_num, "");
		t_gpar par_dypole("", t_gpar::DYPOLE, ++ipar_num, "");
		t_gpar par_ut1("", t_gpar::UT1, ++ipar_num, "");
		t_gpar par_dut1("", t_gpar::DUT1, ++ipar_num, "");
	
		map<string, double> data1, data2;
		t_gtime eop_beg_t,eop_end_t;
		try
		{
			if (!_gdata_erp) _gdata_erp = dynamic_cast<t_gerp*>((*_gall_proc)[t_gdata::ALLERP]);
			if (!_gdata_erp) return false;
			eop_beg_t.from_mjd(_beg_time.mjd(),0,0.0);
			eop_end_t.from_mjd(_beg_time.mjd()+1,0,0.0);
			_gdata_erp->getEopData(eop_beg_t, data1);
			_gdata_erp->getEopData(eop_end_t, data2);

			eop_end_t = eop_end_t - _obs_intv;
		}
		catch (exception e)
		{
			write_log_info(_glog, 0, "ERROR", e.what());
			return false;
		}

		double sig_xpole = dynamic_cast<t_gsetpar*>(_gset)->sigXpole();
		double sig_ypole = dynamic_cast<t_gsetpar*>(_gset)->sigYpole();
		double sig_dxpole = dynamic_cast<t_gsetpar*>(_gset)->sigDxpole();
		double sig_dypole = dynamic_cast<t_gsetpar*>(_gset)->sigDypole();
		double sig_ut1 = dynamic_cast<t_gsetpar*>(_gset)->sigUt1();
		double sig_dut1 = dynamic_cast<t_gsetpar*>(_gset)->sigDut1();

		par_xpole.setTime(eop_beg_t, eop_end_t);
		par_xpole.value(data1["XPOLE"]);
		par_xpole.apriori(sig_xpole);
		par_ypole.setTime(eop_beg_t, eop_end_t);
		par_ypole.value(data1["YPOLE"]);
		par_ypole.apriori(sig_ypole);
		par_dxpole.setTime(eop_beg_t, eop_end_t);
		par_dxpole.value(data2["XPOLE"] - data1["XPOLE"]);
		par_dxpole.apriori(sig_dxpole);
		par_dypole.setTime(eop_beg_t, eop_end_t);
		par_dypole.value(data2["YPOLE"] - data1["YPOLE"]);
		par_dypole.apriori(sig_dypole);
		par_ut1.setTime(eop_beg_t, eop_end_t);
		par_ut1.value(data1["UT1-TAI"]);
		par_ut1.apriori(sig_ut1);
		par_dut1.setTime(eop_beg_t, eop_end_t);
		par_dut1.value(data2["UT1-TAI"] - data1["UT1-TAI"]);
		par_dut1.apriori(sig_dut1);

		// add the Estimator
		lsq->add_parameter(par_xpole);
		lsq->add_parameter(par_ypole);
		lsq->add_parameter(par_dxpole);
		lsq->add_parameter(par_dypole);
		lsq->add_parameter(par_ut1);
		lsq->add_parameter(par_dut1);

		return true;
	}

	bool t_gfltproc::_init_rec_trop_pars(t_glsq* lsq)
	{
		if (_rec_list.empty() || !_lsq || !_gset) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{
			// init the ztd
			string rec = *it;
			t_gpar par_ztd(rec, t_gpar::TRP, ++ipar_num, "");
			shared_ptr<t_gobj> _grec = _gall_obj->obj(rec);
			double sig_ztd = dynamic_cast<t_gsetpar*>(_gset)->sigZtd(rec);
			double sig_TropPd = dynamic_cast<t_gsetpar*>(_gset)->sigTropPd(rec);
			par_ztd.apriori(sig_ztd);
			par_ztd.value(0.0);
			par_ztd.setTime(_beg_time, _beg_time + _pwc_intv * 60.0);
			lsq->add_parameter(par_ztd);
			if (_ztd_model == PWC)
			{
				lsq->add_par_state_equ(par_ztd.parType, 1, _pwc_intv / 60.0, sig_TropPd);
			}

		}
		return true;
	}


	bool t_gfltproc::_init_rec_clk_pars(t_glsq* lsq)
	{

		if (_rec_list.empty() || !lsq || !_gset)
		{
			return false;
		}

		string site_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{

			int ipar_num = lsq->_x_solve.parNumber();

			string site = *it;
			double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigRecCLK(site);
			double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();

			if (site == site_ref)
			{
				if (sig_refclk == 0.0)
				{
					//it++;
					continue;
				}
				else
				{
					sigclk = sig_refclk;
				}
			}

			//intit the Clock
			t_gpar par_clk(site, t_gpar::CLK, ++ipar_num, "");
			par_clk.value(0.0);
			par_clk.setTime(_beg_time, _beg_time);
			par_clk.apriori(sigclk);
			lsq->add_parameter(par_clk);
		}
		return true;
	}


	bool t_gfltproc::_init_sat_clk_pars(t_glsq* lsq)
	{
		if (!lsq || !_gset || _sat_list.empty()) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		auto it_beg = _sat_list.begin();
		auto it_end = _sat_list.end();
		string sat_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();
		
		for (auto iter = it_beg; iter != it_end; iter++)
		{
			string satID = *iter;
			double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigSatCLK(satID);
			if (satID == sat_ref && sig_refclk == 0.0)
			{
				continue;
			}
			else if (satID == sat_ref && sig_refclk != 0.0)
			{
				sigclk = sig_refclk;
			}
			else
			{
				sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigSatCLK(satID);
			}

			t_gpar par_satclk("", t_gpar::CLK_SAT, ++ipar_num, satID);
			par_satclk.setTime(_beg_time, _beg_time);
			par_satclk.apriori(sigclk);
			par_satclk.value(0.0);
			lsq->add_parameter(par_satclk);
		}
		return true;
	}


	bool t_gfltproc::_init_gns_ics_pars(t_glsq* lsq)
	{
		if (!_gall_gns) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		// cycle sat list

		set<string>::const_iterator sat_end = _sat_list.end();
		set<string>::const_iterator sat_beg = _sat_list.begin();
		for (auto sat_iter = sat_beg; sat_iter != sat_end; sat_iter++)
		{
			auto ics = _gall_gns->get_nav_ics(*sat_iter);
			auto par_list = ics->get_par_list();
			if (par_list.empty())
			{
				return false;
			}
			_orb_par_num[*sat_iter] = par_list.size();
			// cycle par list
			auto it_beg = par_list.begin();
			auto it_end = par_list.end();
			for (auto par_iter = it_beg; par_iter != it_end; par_iter++)
			{
				t_gpar::t_type par_type = t_gics_nav::_icspar2gpar(par_iter->first);

				// init 
				t_gpar par_temp("", par_type, ++ipar_num, *sat_iter);
				//TEST change last_time
				par_temp.setTime(par_iter->second.beg_time, LAST_TIME);
				//par_temp.setTime(par_iter->second.beg_time, par_iter->second.end_time);
				par_temp.value(par_iter->second.value);

				auto sig_pos = dynamic_cast<t_gsetpar*>(_gset)->sigSatPos(*sat_iter);
				auto sig_vel = dynamic_cast<t_gsetpar*>(_gset)->sigSatVel(*sat_iter);
				auto sig_force = dynamic_cast<t_gsetpar*>(_gset)->sigSatEcom(*sat_iter);
				string id = "sig" + par_temp.str_type().substr(0, 2);
				if (par_type >= t_gpar::PXSAT && par_type <= t_gpar::PZSAT)
				{
					par_temp.apriori(sig_pos[id] / 1000.0); // change unit to km
				}
				else if (par_type >= t_gpar::VXSAT && par_type <= t_gpar::VZSAT)
				{
					par_temp.apriori(sig_vel[id] / 1000.0); // change unit to km
				}
				else if (par_type >= t_gpar::ECOM2_D0 && par_type <= t_gpar::ECOM2_X1s)
				{
					string ecom_id = "sig" + par_temp.str_type().substr(0, par_temp.str_type().rfind('_'));
					double tmp = sig_force[ecom_id];
					if (double_eq(tmp, 0.0)) tmp = 0.1;
					par_temp.apriori(tmp);
				}
				lsq->add_parameter(par_temp);
			}
		}

		return true;
	}



	// set sys bia pars
	bool t_gfltproc::_init_sys_bia_pars(t_glsq* lsq)
	{
		if (!lsq || !_gset || _rec_list.empty()) return false;
		if (_sys_list.size() < 2) return true;

		// Keep one constellation (GPS when active) as the receiver-clock reference so ISBs stay full-rank
		string ref_sys = _sys_list.count("GPS") ? "GPS" : *_sys_list.begin();
		int ipar_num = lsq->_x_solve.parNumber();

		static const map<string, t_gpar::t_type> isb_map = {
			{"GLO", t_gpar::GLO_ISB},
			{"GAL", t_gpar::GAL_ISB},
			{"BDS", t_gpar::BDS_ISB},
			{"QZS", t_gpar::QZS_ISB},
		};

		for (const auto& rec : _rec_list) {
			for (const auto& sys : _sys_list) {
				if (sys == ref_sys) continue;
				auto it = isb_map.find(sys);
				if (it == isb_map.end()) continue;
				t_gpar par(rec, it->second, ++ipar_num, "");
				_set_isb_par_setting(lsq, par, 0.0, 9000.0);
			}
		}
		return true;
	}


	bool t_gfltproc::_initfltProcPars(t_glsq* lsq)
	{
		return false;
	}

	bool t_gfltproc::_initfltProcData(t_gallproc* data)
	{
		return false;
	}
	
	bool t_gfltproc::_select_rec_obs(const string& rec, const t_gtime& epoch, vector<t_gsatdata>& all_obs, vector<t_gsatdata>& rec_obs)
	{		
		// ===================================================================================================================
		
		if (!_slip12)
		{
			_glog->logError("t_gfltproc", "_select_rec_obs", "_slip 12 is nullptr");
			return false;
		}


		// ===================================================================================================================
		// get all the obs for the rec
		int  index_cout = -1;                       
		map<string, map <string, int> > valid_index;
		for (auto& iter : all_obs)
		{
			index_cout++;                           
			const string& obs_rec = iter.site();
			const string& obs_sat = iter.sat();

			if (_sat_list.find(obs_sat) == _sat_list.end()) {
				continue;
			}

			if (obs_rec != rec) continue;

		
			bool obs_valid = true;
			if (true)
			{
				set<GOBSBAND> band_avail = iter.band_avail();
				bool freq1_ok = false;
				for (auto itband : band_avail)
				{
					FREQ_SEQ freq2 = _freq_index[iter.gsys()][itband];
					switch (freq2)
					{
					case FREQ_SEQ::FREQ_1:
						freq1_ok = true;
						break;
					case FREQ_SEQ::FREQ_2:
						iter.tb12(true);
						break;
					case FREQ_SEQ::FREQ_3:
						iter.tb13(true);
						break;
					case FREQ_SEQ::FREQ_4:
						iter.tb14(true);
						break;
					case FREQ_SEQ::FREQ_5:
						iter.tb15(true);
						break;
					default:
						break;
					}
				}
				if (!freq1_ok)
				{
					iter.tb12(false);
					iter.tb13(false);
					iter.tb14(false);
					iter.tb15(false);
				}
				
			}
			else
			{
				obs_valid = _slip12->use_of_obs(obs_rec, obs_sat, epoch);
			}
			if (obs_valid)
			{
				rec_obs.push_back(iter);
				valid_index[obs_rec][obs_sat] = index_cout;     
			}
		}

		// ===================================================================================================================
		if (rec_obs.empty())
		{
			_glog->logDebug("t_gfltproc", "_select_rec_obs", "This epoch have no useful data : " + rec);
			return false;
		}

		
			return true;
		

		// ===================================================================================================================
		// set tb log in satdata
		for (auto& iter : rec_obs)
		{
			const string& obs_rec = iter.site();
			const string& obs_sat = iter.sat();

			_apply_carrier_range(obs_rec, obs_sat, epoch, iter);
			bool obs_12 = _slip12->use_of_obs(obs_rec, obs_sat, epoch);
			iter.tb12(obs_12);
			all_obs[valid_index[obs_rec][obs_sat]].tb12(obs_12);              // add for PPP-RTK

			if (_frequency < 3) continue;
			bool obs_13 = _slip13->use_of_obs(obs_rec, obs_sat, epoch);
			iter.tb13(obs_13);
			all_obs[valid_index[obs_rec][obs_sat]].tb13(obs_13);              // add for PPP-RTK

			if (_frequency < 5) continue;
			bool obs_14 = _slip14->use_of_obs(obs_rec, obs_sat, epoch);
			iter.tb14(obs_14);
			all_obs[valid_index[obs_rec][obs_sat]].tb14(obs_14);              // add for PPP-RTK

			bool obs_15 = _slip15->use_of_obs(obs_rec, obs_sat, epoch);
			iter.tb15(obs_15);
			all_obs[valid_index[obs_rec][obs_sat]].tb15(obs_15);              // add for PPP-RTK
		}
		return true;
	}

	bool t_gfltproc::_apply_carrier_range(const string& rec, const string& sat, const t_gtime& epoch, t_gsatdata& obs)
	{
		double c1 = 0, c2 = 0, c3 = 0;
		double upd1 = 0, upd2 = 0, upd3 = 0;
		GSYS gsys = obs.gsys();
		GOBS gobsL1 = obs.select_phase(_band_index[gsys][FREQ_1]);
		GOBS gobsL2 = obs.select_phase(_band_index[gsys][FREQ_2]);
		GOBS gobsL3 = obs.select_phase(_band_index[gsys][FREQ_3]);
		double lambda1 = obs.wavelength(_band_index[gsys][FREQ_1]);
		double lambda2 = obs.wavelength(_band_index[gsys][FREQ_2]);

		if (_slip12->is_carrier_range(rec, sat, epoch, c1, c2)) {
			double value1 = obs.getobs(gobsL1) - c1 - upd1;
			double value2 = obs.getobs(gobsL2) - c2 - upd2;
			obs.resetobs(gobsL1, value1);
			obs.resetobs(gobsL2, value2);
			if (_use_carrier_range) {
				obs.is_carrier_range(_band_index[gsys][FREQ_1], true);
				obs.is_carrier_range(_band_index[gsys][FREQ_2], true);
			}
		}

		return true;
	}


	bool t_gfltproc::_processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs)
	{
		
		for (auto& obs_tmp : crt_obs)
		{
			if (t_gsys::sat2gsys(obs_tmp.sat()) == GSYS::GLO)
			{
				if (obs_tmp.channel() == 255)
				{
					if (!_gdata_satinfo)
					{
						_glog->logInfo("t_gfltproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful frequency id for: " + obs_tmp.sat()));
						continue;
					}
					else
					{
						obs_tmp.channel(_gdata_satinfo->getGloFid(obs_tmp.sat(), crt_epoch));
					}
				}
			}
		}

		// prepare Obs [ select observation according to [sat_list; cycle_slip]
		map < string, vector<t_gsatdata> >map_site_obs;
		for (const auto& crt_rec : _rec_list)
		{
			// Use the proc_int to judge
			if (_rec_last_proc_time[crt_rec] + _get_rec_proc_int(crt_rec) > _crt_time) {
				continue;
			}
			_rec_last_proc_time[crt_rec] = _crt_time;

			// get all the obs for the rec
			bool select_valid = _select_rec_obs(crt_rec, crt_epoch, crt_obs, map_site_obs[crt_rec]);
			if (!select_valid)
			{
				_glog->logInfo("t_gfltproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful data : " + crt_rec));
				map_site_obs.erase(crt_rec);
				continue;
			}

		}

		
		vector<t_gsatdata> selected_obs;
		for (const auto& item : map_site_obs) {
			selected_obs.insert(selected_obs.end(), item.second.begin(), item.second.end());
		}

		//map<GSYS, vector<GOBSBAND> > tmp;
		chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
		
		_time_update_orb_parameters(crt_epoch, selected_obs);
		
	

		chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
		_remove_par_msec+= chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count();
	

		// check obs size()
		if (selected_obs.empty())
		{
			_glog->logError("t_gfltproc", "_processOneEpoch", crt_epoch.str_mjdsod("crt obs is empty"));
			return false;
		}


		beg_t = chrono::high_resolution_clock::now();
		// loop while recs


		vector<pair<int,string> > vec_sites;
		vector<t_glsqEquationMatrix>  vec_equ;
		for (auto item : map_site_obs) {
			vec_sites.push_back(make_pair(item.second.size(),item.first));
			vec_equ.push_back(t_glsqEquationMatrix());
		}
		//sort(vec_sites.begin(),vec_sites.end());
		t_gmutex add_mtx;
		// multi-thread  Process REC [Use OpenMP]
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
		for (int site_i = 0; site_i < vec_sites.size(); site_i++)
		{
			t_glsqEquationMatrix& equ_temp = vec_equ[site_i];
			string crt_rec = vec_sites[site_i].second;

			// Process this site and get the equations
			bool proc_site = _processOneRec(crt_epoch, crt_rec, map_site_obs[crt_rec],equ_temp);
			if (!proc_site)
			{
				_glog->logInfo("t_gfltproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful equations : " + crt_rec));
				continue;
			}

			equ_temp.set_newamb(dynamic_cast<t_gturboedit*>(_slip12.get()));

		}

		for (auto& equ_temp: vec_equ){
			
			_equ4one_epoch.add_equ(equ_temp);
			
		}

		if (vec_equ.size() == 0) {
			_glog->logError("t_gfltproc", "_processOneEpoch", crt_epoch.str_mjdsod("All rec process fail!"));
			return false;
		}



		end_t = chrono::high_resolution_clock::now();
		_cmb_equ_msec +=chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count(); 
		return true;

	}


	bool t_gfltproc::_processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs,t_glsqEquationMatrix& equ_result)
	{
		if (crt_obs.empty() || crt_rec.empty())
		{
			_glog->logInfo("t_gfltproc", "_processOneRec", crt_epoch.str_mjdsod("have no useful data : " + crt_rec));
			return false;
		}
		int valid_obs_count = 0;

		// loop for Combine equation
		vector<double> codeOmc;
		bool while_valid = false;
		do
		{
			// process one rec one sat and get equ
			// ========================================================================================================
			equ_result.clear_allequ();

			string crt_sat("NONE");
			for (auto it:crt_obs) 
			{
				crt_sat = it.sat();
				if (_sat_list.empty() || _sat_list.find(crt_sat) == _sat_list.end())
				{
					_glog->logInfo("t_gfltproc", "_processOneRec", crt_epoch.str_mjdsod("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "));
					continue;
				}
				t_gtime epoch = crt_epoch;
				if(!_base_model->cmb_equ(epoch, _lsq->_x_solve, it , equ_result))
				{
					_glog->logInfo("t_gfltproc", "_processOneRec", crt_epoch.str_mjdsod("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "));
					continue;
				}
				valid_obs_count++;
			}

			if (valid_obs_count == 0) {
				_glog->logInfo("t_gfltproc", "_processOneRec", crt_epoch.str_mjdsod("All sat observation process fail!" + crt_rec));
				return false;
			}

			
			if (_is_prefix_res)
			{
				if (this->_prefix_residuals(crt_obs))
				{
					while_valid = true;
					continue;
				}
			}

			codeOmc = equ_result.get_codeomc();
			// =======================================================================================================

			//reset recclk omc, check range
			// =======================================================================================================
			int idx = _lsq->_x_solve.getParam(crt_rec, t_gpar::CLK, "");
			double recclk = 0.0;
			double sigma = 0.0;
			int k = Check_Range(codeOmc, recclk, sigma);
			
			// fix according LX
			if (idx >= 0 && !while_valid && (fabs(recclk / CLIGHT) > _recclk_threshold/*1e-8*/ ||
				(k >= 1 && fabs(_lsq->_x_solve[idx].value() / CLIGHT) < 1e-15)))
			{
				// reset recclk
				recclk += _lsq->_x_solve[idx].value();
				_lsq->_x_solve[idx].value(recclk);
				_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk / CLIGHT);

				codeOmc.clear();
				while_valid = true;
			}
			else
			{
				while_valid = false;
			}

			// =======================================================================================================
		} while (while_valid);



		return true;				
	
	}


	bool t_gfltproc::_init_sat_clk()
	{
		set<string> sat_list = _crt_equ.get_satlist(_crt_rec);
		for (auto sat_iter = sat_list.begin(); sat_iter != sat_list.end(); sat_iter++)
		{
			// get satclk par idx
			int idx_satclk = _lsq->_x_solve.getParam("", t_gpar::CLK_SAT, *sat_iter);
			if (idx_satclk >= 0 && double_eq(_lsq->_x_solve[idx_satclk].value(),0.0))
			{
				try
				{
					// get all phase code observ equations;
					vector<int> equ_list = _crt_equ.find_equ(_crt_rec, *sat_iter);
					double satclk = 0.0;
					// get satclk delay initial value 
					for (int idx : equ_list)
					{
						if (_crt_equ.get_obscombtype(idx).is_code())
						{
							satclk = -(_crt_equ.l[idx] + 1.0);
							break;
						}
					}
					// clk par init
					_lsq->_x_solve[idx_satclk].value(satclk);
					_bias_model->update_obj_clk(*sat_iter, _crt_time, satclk/CLIGHT);

					// omc = omc + satclk;
					for (int idx : equ_list)
					{
						_crt_equ.l[idx] += satclk;
					}
				}
				catch (exception e)
				{
					write_log_info(_glog, 1, "ERROR", e.what());
					return false;
				}
			}
		}
		return true;
	}

	bool t_gfltproc::_init_xml_settings(t_gsetbase* set)
	{
		// for all values with time
		_beg_time = dynamic_cast<t_gsetgen*>(set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(set)->end();
		_obs_intv = dynamic_cast<t_gsetgen*>(set)->sampling();

		// for all values with mode
		_write_equ = dynamic_cast<t_gsetproc*>(set)->write_equ();
		_fix_mode = dynamic_cast<t_gsetamb*>(set)->fix_mode();
		_fix_dd_mode= dynamic_cast<t_gsetamb*>(set)->fix_dd_mode();

		_slip_model = dynamic_cast<t_gsetproc*>(_gset)->slip_model();
		_iono_order = dynamic_cast<t_gsetproc*>(_gset)->iono_order();
		_ztd_model = dynamic_cast<t_gsetproc*>(_gset)->ztd_model(&_pwc_intv);
		_ifb_model = dynamic_cast<t_gsetproc*>(_gset)->ifb_model();
		_frequency = dynamic_cast<t_gsetproc*>(set)->frequency();
		_ambfix_xml = dynamic_cast<t_gsetproc*>(_gset)->ambfix();
		_tropo_grad = dynamic_cast<t_gsetproc*>(_gset)->tropo_grad();   
		

		// for all settings with sat and rec
		_rec_list = dynamic_cast<t_gsetgen*>(set)->rec_all();
		_sat_list = dynamic_cast<t_gsetgnss*>(set)->sat();
		_sys_list = dynamic_cast<t_gsetgen*>(set)->sys();
		_rover_list = dynamic_cast<t_gsetgen*>(set)->list_rover();

		_proc_int = dynamic_cast<t_gsetgen*>(set)->proc_int();
		for (const auto& rec : _rec_list)
		{
			_rec_last_proc_time[rec] = FIRST_TIME;
		}

		//POD mode
		_mode = dynamic_cast<t_gsetgen*>(set)->mode();
		_num_threads = dynamic_cast<t_gsetproc*>(set)->num_threads();
#ifdef USE_OPENMP
		omp_set_num_threads(_num_threads);
#endif




		
		_is_prefix_res = dynamic_cast<t_gsetproc*>(set)->pre_fix_residuals();
	
		_carrier_range = dynamic_cast<t_gsetproc*>(set)->carrier_range();
		_use_carrier_range = dynamic_cast<t_gsetproc*>(set)->apply_carrier_range();
		_ultrasp3 = dynamic_cast<t_gsetproc*>(set)->ultrasp3(); // only used in realtime
		_maxres_norm = dynamic_cast<t_gsetproc*>(set)->max_res_norm();
		sys_support = dynamic_cast<t_gsetgnss*>(set)->SYS_SUPPORTED();
		for (auto it = sys_support.begin(); it != sys_support.end(); it++) {
			_maxres_C[*it] = dynamic_cast<t_gsetgnss*>(set)->maxres_C(*it);
			_maxres_L[*it] = dynamic_cast<t_gsetgnss*>(set)->maxres_L(*it);
		}
		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::QZS);

		
		std::set<string> tmp_sat_list;
		for (const auto& sat : _sat_list)
		{
			auto gsys = t_gsys::sat2gsys(sat);
			auto ssys = t_gsys::gsys2str(gsys);
			if (_sys_list.count(ssys) == 0) continue;
			tmp_sat_list.insert(sat);
		}
		_sat_list = tmp_sat_list;

		if (_sat_list.empty())
		{
			for (auto sys_iter = _sys_list.begin(); sys_iter != _sys_list.end(); sys_iter++)
			{
				auto sats_temp = dynamic_cast<t_gsetgnss*>(set)->sat(t_gsys::str2gsys(*sys_iter));
				_sat_list.insert(sats_temp.begin(), sats_temp.end());
			}
		}

		std::set<string> sat_rm = dynamic_cast<t_gsetgen*>(set)->sat_rm();
		for (auto iter : sat_rm) { _sat_list.erase(iter); }


		// Get site coordinates from XML/rinex O files
		std::set<string> rec_rm;
		for (auto it : _rec_list)
		{
			
			_rec_crds[it] = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(it);
			_rec_stds[it] = dynamic_cast<t_gsetrec*>(_gset)->get_std_xyz(it);
			if (_rec_crds[it].zero())
			{
				rec_rm.insert(it);
			}
		}
		return true;
	}

	bool t_gfltproc::_init_lsq_settings(t_gsetbase* set, t_gallproc* data, t_glog* log)
	{
		// =========================================================================================================
		// cycle slip
		if (_frequency > 6 || _frequency == 0)
		{
			throw logic_error("you should set the frequency in xml as 1/2/3/4/5.");
		}

		if (_frequency >= 2) _slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log));
        if (_frequency >= 3) _slip13 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 3));
		if (_frequency >= 4) _slip14 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 4));
		if (_frequency >= 5) _slip15 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 5));

		// =========================================================================================================
		// Nothing but cout a string here
		// Init Paramter in Constructor, make sure lsq is new class
		if (_lsq) { delete _lsq; _lsq = nullptr;} _lsq = new t_glsq(_gset);
		
		// Nothing but cout a string here
		 _bias_model = shared_ptr<t_gbiasmodel>(new t_gprecisebias(data, log, set));
		
		if (_frequency == 2)
		{
			
				shared_ptr<t_gupdatepar> updateIF(new t_gupdateparIF(_slip12));
				updateIF->set_amb_update_way(true);
				updateIF->set_interval(_obs_intv);
				_lsq->set_update_par(updateIF);

				_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF(set, log, _bias_model, data));
			

		}
		else
		{
			throw logic_error("Support IF[DOUBLE/TRIPLE]/UDUC[5 FREQ]");
		}

		_lsq->get_update_par()->set_freq_band_index(_band_index, _freq_index);

		return true;
	}

	bool t_gfltproc::_init_rec_crd_pars(t_glsq* lsq, bool isPPP)
	{
		int ipar_num = _lsq->_x_solve.parNumber();

		// Here ensure the [rec_ground] sites only include the gournd sites 
		set<string> rec_tmp = _gall_obj->objects_set(t_gdata::REC);
		set<string> rec_ground;
		set_intersection(rec_tmp.begin(), rec_tmp.end(), _rec_list.begin(), _rec_list.end(), insert_iterator<set<string>>(rec_ground, rec_ground.begin()));


		if (rec_ground.empty())
		{
			write_log_info(_glog, 0, "ERROR : The gournd rec list is empty !");
			return false;
		}

		auto it_beg = rec_ground.begin();
		auto it_end = rec_ground.end();
		int missing = 0;
		for (auto it = it_beg; it != it_end; )
		{
			string site = *it;
			shared_ptr<t_gobj> _grec = _gall_obj->obj(site);

			// get the init value

			t_gtriple rec_crd;
			t_gtriple rec_std;

			if (isPPP )
			{
				_grec->get_recent_crd(_beg_time, 0.1, rec_crd, rec_std);
			}
			else
			{
				rec_crd = _rec_crds[site];
				rec_std = _rec_stds[site];
			}

			if (rec_std.zero() || isPPP)
			{
				rec_std = dynamic_cast<t_gsetpar*>(_gset)->sigRecPosTriple(site);
			}

	
			// only need sinex file 
			t_gtriple snx_crd, snx_std;
			if (!_grec->get_recent_crd(_beg_time,0.1,snx_crd,snx_std))
			{
				// no SINEX coordinate (sigma < 0.1 m): previously the iterator was never advanced (endless loop).
				// Fail with the station names instead; the station list and the SINEX input must match.
				write_log_info(_glog, 0, "ERROR : no SINEX coordinate (sigma < 0.1 m) for station " + site);
				cerr << "ERROR: no SINEX coordinate (sigma < 0.1 m) for station " << site << endl;
				missing++;
			}
			else
			{
				_grec->crd(snx_crd, snx_std,FIRST_TIME,LAST_TIME,true);
			}
			it++;
		}
		return missing == 0;
	}



	bool t_gfltproc::_create_log(t_glog* log)
	{
		if (!log)
		{
			this->_glog = new t_glog();
			this->_glog->cache_size(99);
			this->_glog->tsys(t_gtime::GPS);
			this->_glog->time_stamp(true);
			this->_glog->verb(1);
			return true;
		}
		return false;
	}

	bool t_gfltproc::_delete_log(t_glog* log)
	{
		if (log)
		{
			delete log;
			log = nullptr;
			return true;
		}
		return false;
	}


	bool t_gfltproc::_pod_add_Dd_ambcon(t_SRIF * flt, t_gambcon * ambcon)
	{
		if (flt == nullptr || ambcon == nullptr) return false;

		t_glsqEquationMatrix virtual_equ;
		auto ddambs = ambcon->get_Dds();
		if (ddambs.empty()) return false;
		double op[5] = { 0.0, 1.0, -1.0, -1.0, 1.0 };
		// apriori
		double fix_apriori = this->_fix_dd_mode==t_gsetamb::FIX_DD_MODE::IND?1e8:5e5;

		int debug_constriant_num = 0;
		for (const auto& iter : ddambs)
		{
			const auto& Dd_amb = iter;

			if (!(Dd_amb.MW_fixed && Dd_amb.NL_fixed))
			{
				continue;
			}
			//ddamb=0
			if (double_eq(Dd_amb.MW_value, 0.0) || double_eq(Dd_amb.NL_value, 0.0)) continue;

			string rec_1st = Dd_amb.rec_1st;
			string rec_2nd = Dd_amb.rec_2nd;
			string sat_1st = Dd_amb.sat_1st;
			string sat_2nd = Dd_amb.sat_2nd;



			const t_gtime& beg = Dd_amb.beg;
			const t_gtime& end = Dd_amb.end;

			int par_1st = -1;
			int par_2nd = -1;
			int par_3rd = -1;
			int par_4th = -1;

			if (Dd_amb.ID == t_gamb::AMB_ID::IF_L1_L2)
			{
				par_1st = _lsq->_x_solve.getParam(rec_1st, t_gpar::AMB_IF, sat_1st);
				par_2nd = _lsq->_x_solve.getParam(rec_1st, t_gpar::AMB_IF, sat_2nd);
				par_3rd = _lsq->_x_solve.getParam(rec_2nd, t_gpar::AMB_IF, sat_1st);
				par_4th = _lsq->_x_solve.getParam(rec_2nd, t_gpar::AMB_IF, sat_2nd);
			}
			else
			{
				return false;
			}


			if (par_1st < 0 || par_2nd < 0 || par_3rd < 0 || par_4th < 0)
			{
				continue;
			}


			if (double_eq(_lsq->_x_solve.getAmbParam(par_1st).value(), 0.0)) continue;
			if (double_eq(_lsq->_x_solve.getAmbParam(par_2nd).value(), 0.0)) continue;
			if (double_eq(_lsq->_x_solve.getAmbParam(par_3rd).value(), 0.0)) continue;
			if (double_eq(_lsq->_x_solve.getAmbParam(par_4th).value(), 0.0)) continue;

			t_gobsgnss gnss1,gnss2;
			gnss1.sat(sat_1st);
			gnss2.sat(sat_2nd);

			auto gsys = t_gsys::sat2gsys(sat_1st);

			if (gsys == GSYS::GLO) {
				gnss1.channel(_glofrq_num[sat_1st]);
				gnss2.channel(_glofrq_num[sat_2nd]);
			}

			double f11 = gnss1.frequency(_band_index[gsys][FREQ_1]);
			double f12 = gnss1.frequency(_band_index[gsys][FREQ_2]);

			double wlength_1st = CLIGHT / gnss1.frequency(_band_index[gsys][FREQ_1]);
			double wlength_2nd = 0.0;
			if (Dd_amb.ID == t_gamb::AMB_ID::IF_L1_L2) { wlength_2nd = CLIGHT / gnss1.frequency(_band_index[gsys][FREQ_2]); }
			else
			{
				return false;
			}

			double rw_nl1 = (f11 + f12) / CLIGHT;
			double rw_nl2 = rw_nl1;
			
			if (gsys == GSYS::GLO) {
				double f21 = gnss2.frequency(_band_index[gsys][FREQ_1]);
				double f22 = gnss2.frequency(_band_index[gsys][FREQ_2]);

				rw_nl2 = (f21 + f22) / CLIGHT;
			}

			double k = rw_nl2 / rw_nl1;

			double nc = _lsq->_x_solve.getAmbParam(par_1st).value() * op[1] +
				_lsq->_x_solve.getAmbParam(par_2nd).value() * op[2] * k +
				_lsq->_x_solve.getAmbParam(par_3rd).value() * op[3] +
				_lsq->_x_solve.getAmbParam(par_4th).value() * op[4] * k;

			double lambdaw = wlength_1st / (1 - wlength_1st / wlength_2nd);
			double lambdan = wlength_1st / (1 + wlength_1st / wlength_2nd);
			
			double nc_fix = round(Dd_amb.NL_value) + round(Dd_amb.MW_value) * (f12 / (f11 - f12)); // cycle repsentation
			nc_fix = nc_fix / rw_nl1;


			double nl = nc_fix - nc;


			map<int, double> par_index;
			par_index[par_1st] = op[1];
			par_index[par_2nd] = op[2] * k;
			par_index[par_3rd] = op[3];
			par_index[par_4th] = op[4] * k;

			vector< pair<int, double> > B;
			for (const auto& B_iter : par_index)
			{
				B.push_back(make_pair(B_iter.first + 1, B_iter.second));
			}

			t_gobscombtype type;

			// For IGSO and GEO, 1e6 or 1e5 maybe better
			virtual_equ.add_equ(B, fix_apriori, nl, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, type, false);
			_glog->comment(t_glog::LOG_LV::LOG_WARN, "t_gpodfltIF", "_pod_add_Dd_ambcon", _crt_time.str_mjdsod("one ambcon add success"));

		}

		flt->measurement_update(virtual_equ);


		return true;
	}

	bool t_gfltproc::_pod_add_Dd_ambcon(t_SRIF * flt, t_gambcon * ambcon, great::t_gupdateparinfo remove_info)
	{
		if (flt == nullptr || ambcon == nullptr) return false;
	
		auto ddambs = ambcon->get_Dds();
		if (ddambs.empty()) return true;

		// find removed ambiguity parameters; remove id from 1
		vector<int> remove_id;
		remove_info.get(remove_id);
		set<int> amb_table; // pair(rec,sat)
		for (auto id : remove_id) {
			if (_lsq->_x_solve[id - 1].parType == t_gpar::AMB_IF) {
				amb_table.insert(id - 1);
			}
		}

		if (amb_table.empty()) {
			return true;
		}

		double op[5] = { 0.0, 1.0, -1.0, -1.0, 1.0 };
		t_glsqEquationMatrix virtual_equ;
		for (const auto& iter : ddambs)
		{
			const auto& Dd_amb = iter;
			if (!(Dd_amb.MW_fixed && Dd_amb.NL_fixed))
			{
				continue;
			}
			//ddamb=0
			if (double_eq(Dd_amb.MW_value, 0.0) || double_eq(Dd_amb.NL_value, 0.0)) continue;

			string rec_1st = Dd_amb.rec_1st;
			string rec_2nd = Dd_amb.rec_2nd;
			string sat_1st = Dd_amb.sat_1st;
			string sat_2nd = Dd_amb.sat_2nd;

			// add for single difference


			const t_gtime& beg = Dd_amb.beg;
			const t_gtime& end = Dd_amb.end;

			int par_1st = -1;
			int par_2nd = -1;
			int par_3rd = -1;
			int par_4th = -1;

			if (Dd_amb.ID == t_gamb::AMB_ID::IF_L1_L2)
			{
				par_1st = _lsq->_x_solve.getAmbParam(rec_1st, sat_1st, t_gpar::AMB_IF, beg, end);
				par_2nd = _lsq->_x_solve.getAmbParam(rec_1st, sat_2nd, t_gpar::AMB_IF, beg, end);
				par_3rd = _lsq->_x_solve.getAmbParam(rec_2nd, sat_1st, t_gpar::AMB_IF, beg, end);
				par_4th = _lsq->_x_solve.getAmbParam(rec_2nd, sat_2nd, t_gpar::AMB_IF, beg, end);
			}
			else
			{
				return false;
			}


			if (par_1st < 0 || par_2nd < 0 || par_3rd < 0 || par_4th < 0)
			{
				continue;
			}

			if (!(amb_table.count(par_1st) || amb_table.count(par_2nd) || amb_table.count(par_3rd) || amb_table.count(par_4th))) {
				continue;
			}

			// erase Dd where common time smaller than 1 hour && baseline > 1km
			if (Dd_amb.end - Dd_amb.beg < 7200) {
				continue;
			}


			// constraint X(ix1)-X(ix2)-X(ix3)+X(ix4)=nc_fix-nc

			if (double_eq(_lsq->_x_solve.getAmbParam(par_1st).value(), 0.0)) continue;
			if (double_eq(_lsq->_x_solve.getAmbParam(par_2nd).value(), 0.0)) continue;
			if (double_eq(_lsq->_x_solve.getAmbParam(par_3rd).value(), 0.0)) continue;
			if (double_eq(_lsq->_x_solve.getAmbParam(par_4th).value(), 0.0)) continue;

			double nc = _lsq->_x_solve.getAmbParam(par_1st).value() * op[1] +
				_lsq->_x_solve.getAmbParam(par_2nd).value() * op[2] +
				_lsq->_x_solve.getAmbParam(par_3rd).value() * op[3] +
				_lsq->_x_solve.getAmbParam(par_4th).value() * op[4];

			// jdhuang : fix for muti-GNSS
			t_gobsgnss gnss;
			gnss.sat(sat_1st);
			auto gsys = t_gsys::sat2gsys(sat_1st);
			double wlength_1st = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_1]);
			double wlength_2nd = 0.0;
			if (Dd_amb.ID == t_gamb::AMB_ID::IF_L1_L2) { wlength_2nd = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_2]); }
			else
			{
				return false;
			}

			double lambdaw = wlength_1st / (1 - wlength_1st / wlength_2nd);
			double lambdan = wlength_1st / (1 + wlength_1st / wlength_2nd);
			double nc_fix = round(Dd_amb.NL_value) + round(Dd_amb.MW_value) * lambdaw / wlength_2nd;
			nc_fix = nc_fix * lambdan;
			double nl = nc_fix - nc;


			map<int, double> par_index;
			par_index[par_1st] = op[1];
			par_index[par_2nd] = op[2];
			par_index[par_3rd] = op[3];
			par_index[par_4th] = op[4];

			vector< pair<int, double> > B;
			for (const auto& B_iter : par_index)
			{
				B.push_back(make_pair(B_iter.first + 1, B_iter.second));
			}

			t_gobscombtype type;

			// For IGSO and GEO, 1e6 or 1e5 maybe better
			virtual_equ.add_equ(B, 1e5, nl, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, type, false);

			_glog->comment(t_glog::LOG_LV::LOG_WARN, "t_gpodfltIF", "_pod_add_Dd_ambcon", _crt_time.str_mjdsod("one ambcon add success"));

		}
		flt->measurement_update(virtual_equ);

		return true;
	}

	bool t_gfltproc::_pod_add_one_Sd_ambcon(t_glsq* lsq, const t_gamb_Dd* Dd_amb)
	{
		string rec = Dd_amb->rec_1st;
		string sat_1st = Dd_amb->sat_1st;
		string sat_2nd = Dd_amb->sat_2nd;
		double op[3] = { 0.0, 1.0, -1.0 };
		if (double_eq(Dd_amb->MW_value, 0.0) || double_eq(Dd_amb->NL_value, 0.0)) return false;
		const t_gtime& beg = Dd_amb->beg;
		const t_gtime& end = Dd_amb->end;
		int par_1st = lsq->_x_solve.getAmbParam(rec, sat_1st, t_gpar::AMB_IF, beg, end);
		int par_2nd = lsq->_x_solve.getAmbParam(rec, sat_2nd, t_gpar::AMB_IF, beg, end);
		if (par_1st < 0 || par_2nd < 0)
		{
			return false;
		}
		if (double_eq(lsq->_x_solve.getAmbParam(par_1st).value(), 0.0)) return false;
		if (double_eq(lsq->_x_solve.getAmbParam(par_2nd).value(), 0.0)) return false;
		double nc = lsq->_x_solve.getAmbParam(par_1st).value() * op[1] +
			lsq->_x_solve.getAmbParam(par_2nd).value() * op[2];

		t_gobsgnss gnss;
		gnss.sat(sat_1st);
		auto gsys = t_gsys::sat2gsys(sat_1st);
		double wlength_1st = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_1]);
		double wlength_2nd = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_2]);

		double lambdaw = wlength_1st / (1 - wlength_1st / wlength_2nd);
		double lambdan = wlength_1st / (1 + wlength_1st / wlength_2nd);
		double nc_fix = round(Dd_amb->NL_value) + round(Dd_amb->MW_value) * lambdaw / wlength_2nd;
		nc_fix = nc_fix * lambdan;
		nc_fix = nc_fix - Dd_amb->rnl*lambdan;
		double nl = nc_fix - nc;
		map<int, double> par_index;
		par_index[par_1st] = op[1];
		par_index[par_2nd] = op[2];
		vector< pair<int, double> > B;
		for (const auto& B_iter : par_index)
		{
			B.push_back(make_pair(B_iter.first + 1, B_iter.second));
		}
		t_glsqEquationMatrix virtual_equ;
		t_gobscombtype type;
		virtual_equ.add_equ(B, 1e7, nl, rec, sat_1st + "_" + sat_2nd, type, false);
		lsq->add_equation(virtual_equ);
		return true;
	}

	void t_gfltproc::_set_isb_par_setting(t_glsq* lsq, t_gpar& par,  double value, double apriori)
	{
		par.value(value);
		par.apriori(apriori);
		par.setTime(_beg_time, _beg_time);
		lsq->add_parameter(par);
		lsq->add_par_state_equ(par.parType, 1, 24, 3000.0);
		
	}


	bool t_gfltproc::_prefix_residuals(vector<t_gsatdata>& vec_obs)
	{
		set<string> bad_sats;
		string sat, type;
		GSYS gs;
		for (int iequ = 0; iequ < _crt_equ.num_equ(); iequ++)
		{
			sat = _crt_equ.get_satname(iequ);
			type = _crt_equ.get_obscombtype2str(iequ);
			gs = t_gsys::str2gsys(sat.substr(0, 1));

			if (type.empty())
			{
				continue;
			}

			if (type[0] == 'L')
			{
				if (fabs(_crt_equ.l[iequ]) > _maxres_L[gs] * 600)         
				{
					bad_sats.insert(sat);
					_glog->logInfo("t_gfltproc", "_prefix_residuals:", " outlier (v_orig: " + type +
						" " + _crt_equ.get_sitename(iequ) + " " + sat + ":" + dbl2str(fabs(_crt_equ.l[iequ])) + ")");
				}
			}
			else if (type[0] == 'P' || type[0] == 'C') 
			{
				if (fabs(_crt_equ.l[iequ]) > _maxres_C[gs] * 300)
				{
					bad_sats.insert(sat);
					_glog->logInfo("t_gfltproc", "_prefix_residuals:", " outlier (v_orig: " + type +
						" " + _crt_equ.get_sitename(iequ) + " " + sat + ":" + dbl2str(fabs(_crt_equ.l[iequ])) + ")");
				}
			}
		}

		bool eraseSat = false;
		auto data = vec_obs.begin();
		while (data != vec_obs.end())
		{
			string sat = data->sat();
			if (bad_sats.find(sat) != bad_sats.end())
			{
				data = vec_obs.erase(data);
				eraseSat = true;
			}
			else
			{
				data++;
			}
		}

		return eraseSat;
	}

	double t_gfltproc::_get_rec_proc_int(const string& rec)
	{
		string rec_upper=rec;
		transform(rec_upper.begin(), rec_upper.end(), rec_upper.begin(), ::toupper);
		if (_proc_int.count(rec_upper) != 0)
		{
			return _proc_int[rec_upper];
		}
		else {
			return 0.0;
		}
	}

	void t_gfltproc::_time_update_orb_parameters(const t_gtime& crt_epoch,vector<t_gsatdata>& obsdata)
	{
		// get other update par
		shared_ptr<t_gupdatepar> update_par = _lsq->get_update_par();
		t_gupdateparinfo orig_remove_info = update_par->get_all_update_parameters(crt_epoch,_lsq->_x_solve,obsdata);

		t_gupdateparinfo remove_info;

		set<string> sats_have_obs;
		map<string, int> sats_obs_num;
		for (const auto& item : obsdata) {
			sats_have_obs.insert(item.sat());
			if (sats_obs_num.count(item.sat()) == 0) {
				sats_obs_num[item.sat()] = 1;
			}
			else {
				sats_obs_num[item.sat()] += 1;
			}
		}
		for (const auto& item : sats_obs_num) {
		}
		cout << endl;
		_get_orb_update_parameters(remove_info, sats_have_obs);
		remove_info.add(orig_remove_info);


		// update sat valid obs length
		for (auto& item : _sat_valid_obs_length) {
			if (_convergence[item.first]) {
				continue;
			}
			if (_crt_time !=_orb_state_ref[item.first]) {
				continue;
			}
			if (sats_have_obs.count(item.first) != 0 ) {
				item.second += _update_orb_interval[item.first];
			}
			else {
				item.second = 0.0;
			}
		}


        // update init value for clk par
        vector<t_gpar>& update_clkparlist = remove_info.get_newparlist();
        for (auto& par: update_clkparlist){
            if (par.parType == t_gpar::CLK || par.parType == t_gpar::CLK_SAT){
                int idx_clk = _lsq->_x_solve.getParam(par.site,par.parType,par.prn);

                // update the initial value with estimated value
                par.value( par.value() + _flt._x(idx_clk));
            }
        }


        chrono::high_resolution_clock::time_point beg_t,end_t;
		beg_t = chrono::high_resolution_clock::now();



        _updateinfo_crt_time = remove_info;
		//update filter
		_flt.time_update_use_removeinfo(remove_info);

        //update lsq par
        vector<int> remove_id;
        vector<t_gpar> equparlist;
        vector<t_gpar> newparlist;
        t_glsqEquationMatrix tempequ;
        remove_info.get(remove_id,newparlist,equparlist, tempequ);

        sort(remove_id.begin(), remove_id.end(), greater<int>());
		for (auto item : remove_id) {
			_lsq->remove_parameter(item);
		}
		
		for (auto& par : equparlist) {
			_lsq->add_parameter(par);
		}
		

		
		for (auto& par: newparlist){
			_lsq->add_parameter(par);
		}
		


		end_t = chrono::high_resolution_clock::now();
		_timeupdate_msec +=chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count(); 

	}

	void t_gfltproc::_update_srif_with_updateinfo(t_gupdateparinfo& update_info)
	{

		//update filter
		_flt.time_update_use_removeinfo(update_info);

		//update lsq par
		vector<int> remove_id;
		vector<t_gpar> equparlist;
		vector<t_gpar> newparlist;
		t_glsqEquationMatrix tempequ;
		update_info.get(remove_id, newparlist, equparlist, tempequ);

		sort(remove_id.begin(), remove_id.end(), greater<int>());
		for (auto item : remove_id) {
			_lsq->remove_parameter(item);
		}
		
		for (auto& par : equparlist) {
			_lsq->add_parameter(par);
		}
		

		
		for (auto& par : newparlist) {
			_lsq->add_parameter(par);
		}
		

	}

	void t_gfltproc::_get_orb_update_parameters(t_gupdateparinfo& orb_update_info, set<string> sats_have_obs)
	{

		// add constriant equations
		int baseParNum = _lsq->_x_solve.parNumber()+orb_update_info.get_update_newpar_number();

		vector<int> pxsat_index_list;
		vector<string> sat_arr;
		vector<bool> isLeo;
		
		for (auto sat : _sat_list) {
			int idx = _lsq->_x_solve.getParam("", t_gpar::PXSAT, sat);
			if (idx == -1) {
				continue;
			}
			sat_arr.push_back(sat);
			pxsat_index_list.push_back(idx);
			isLeo.push_back(false);
		}
	

		t_gmutex add_mtx;
		// with process noise is zero;
		Eigen::MatrixXd phi_i_zero;
		Eigen::VectorXd l_zero;
		vector<t_glsqEquationMatrix> orb_equ_sats(sat_arr.size());
		vector<vector<pair<int, t_gpar> > > equparlist_sats(sat_arr.size());

		vector<string> oi_sats;
		vector<int> oi_intervals;
		vector<vector<t_gpar::t_type> > state_init_name_sats;
		vector<vector<double> >  state_init_value_sats;



		for (int i=0;i<sat_arr.size();i++) {
			string sat = sat_arr[i];

			// judge update orb_state
			if (_crt_time < _orb_state_ref[sat] + _update_orb_interval[sat]) {
				continue;
			}

			if (!isLeo[i] &&
                sats_have_obs.count(sat) == 0 &&
                _crt_time < _orb_ref[sat] + _update_oi_interval[sat] &&
                _re_integration_sats.count(sat) == 0 ) {
				// don't update until have observations
				continue;
			}  // else: if crt_time > orb_ref + update_oi , but still have no observations for this sat , still will go orbit integration

			// interval_time 
			int orb_update_interval = _crt_time - _orb_state_ref[sat];
			int update_epo = orb_update_interval / _update_orb_interval[sat];


			if (orb_update_interval > 3600 * 3) {
				cout << "orb_update_interval " << orb_update_interval << " sat: " << sat << endl;
			}

			// index from 0 
			int pxsat_index = pxsat_index_list[i];

			double eclipse_ratio = _eclipse_marker[sat] ? dynamic_cast<t_gsetforce*>(_gset)->getEclipseRatio(sat) : 1.0;

			Matrix R_orb = eigen2newmat(_orb_mode[sat]->getR(orb_update_interval, eclipse_ratio));
            double vel_noise = _orb_mode[sat]->get_vel_noise(orb_update_interval, eclipse_ratio);

			// phi(t+1,t) = phi(t+1,t0)phi(t,t0)^-1;
			Matrix phi = _phi_map[sat][_epo[sat] + update_epo] * _phi_i_map[sat][_epo[sat]];
			//cout << phi << endl;
			vector<SimpleB> B;
			SymmetricMatrix P;
			vector<double> l;

            // speical for velcotiy pulse
            int vel_pulse_index = _lsq->_x_solve.getParam("",t_gpar::VX_Pulse,sat);
            Matrix rotmat(0.0,3,3);
            if (vel_pulse_index != - 1){
                orb_update_info.add(vel_pulse_index + 1 + 0);
                orb_update_info.add(vel_pulse_index + 1 + 1);
                orb_update_info.add(vel_pulse_index + 1 + 2);

                t_gcrs2acr crs2acr;
                ColumnVector sat_pos(0.0, 3);
                ColumnVector sat_vel(0.0, 3);
                ColumnVector maneuver_sat_state = _state_map[sat][_orb_state_ref[sat]];

                sat_pos = maneuver_sat_state.rows(1,3);
                sat_vel = maneuver_sat_state.rows(4,6);

                // change unit from kilometer to meter
                sat_pos *= 1000.0;
                sat_vel *= 1000.0;

                crs2acr.calcRotMat(sat_pos,sat_vel, rotmat);
            }


            // Construct virtual observation equations
			for (int j = 0; j < _orb_par_num[sat]; j++) {
				vector<pair<int, double>> singleB;
				t_gpar newOrbPar(_lsq->_x_solve[pxsat_index + j]);
				newOrbPar.value(_state_map[sat][_orb_state_ref[sat] + orb_update_interval](j + 1));
				baseParNum++;

				for (int k = 0; k < _orb_par_num[sat]; k++) {
					singleB.push_back(make_pair(pxsat_index + 1 + k, phi(j + 1, k + 1)));
				}

                // add vel pulse
                if (vel_pulse_index != -1 && j < 6){
                    ColumnVector trans_xyz(0.0,3);
                    trans_xyz(1) = singleB[3].second;
                    trans_xyz(2) = singleB[4].second;
                    trans_xyz(3) = singleB[5].second;

                    ColumnVector trans_acr = rotmat * trans_xyz;

                    singleB.push_back(make_pair(vel_pulse_index + 1 + 0,trans_acr(1)));
                    singleB.push_back(make_pair(vel_pulse_index + 1 + 1,trans_acr(2)));
                    singleB.push_back(make_pair(vel_pulse_index + 1 + 2,trans_acr(3)));
                }

				singleB.push_back(make_pair(baseParNum, -1.0));
				B.push_back(singleB);
				l.push_back(0.0);
				equparlist_sats[i].push_back(make_pair(pxsat_index + j + 1, newOrbPar));
			}
			P << R_orb.t() * R_orb;


			// update record
			_epo[sat]+=update_epo;
			_orb_state_ref[sat] = _crt_time;

			// update partial matrix and update initial value due to the reinit of satellite
			if (_crt_time >= _orb_ref[sat] + _update_oi_interval[sat] || _re_integration_sats.count(sat)!=0  ) {
				// l = phi(t+1,t) * dx
				ColumnVector dx(_orb_par_num[sat]);
				for (int j=0;j<_orb_par_num[sat];j++){
					dx(j+1) = _flt._x(pxsat_index+j);
				}
				//ColumnVector newl = (phi) * eigen2newmat(Eigen::VectorXd(_flt._x.segment(pxsat_index, _orb_par_num[sat])));
				ColumnVector newl = phi * dx;
				//  orbit integration
				_state_map[sat][_orb_state_ref[sat]] = _state_map[sat][_orb_state_ref[sat]] + newl;

				// oi judge interval
				int init_time = isLeo[i] ? 900: 3600 * 12;
				if (!_convergence[sat] && _sat_valid_obs_length[sat] >= init_time) {
					_convergence[sat] = true;
					_update_oi_interval[sat]  = dynamic_cast<t_gsetforce*>(_gset)->getRtUpdateInt(sat);
				}

				oi_sats.push_back(sat);
				int oi_interval = _update_oi_interval[sat];
				oi_intervals.push_back(oi_interval);
				state_init_name_sats.push_back(vector<t_gpar::t_type>());
				state_init_value_sats.push_back(vector<double>());
				for (int j = 0; j < _orb_par_num[sat]; j++) {
					state_init_name_sats[state_init_name_sats.size()-1].push_back(_lsq->_x_solve[pxsat_index + j].parType);
					state_init_value_sats[state_init_value_sats.size()-1].push_back(_state_map[sat][_orb_state_ref[sat]](j + 1));
				}

                if (_re_integration_sats.count(sat) !=0){
                    _re_integration_sats.erase(sat);
                }


				for (int j = 0; j < _orb_par_num[sat]; j++) {
					l[j] = newl(j + 1);
				}

			}

			orb_equ_sats[i].add_multi_equ(B, P, l,"PXSAT",sat);
			orb_update_info.add(equparlist_sats[i], orb_equ_sats[i]);
            orb_update_info.inflate_info[sat]=vel_noise;
			
		}

		
		// oi integartion
		chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
		for (int i =0;i<oi_sats.size();i++){
			cerr << "orbit integration for " << oi_sats[i] << endl;
			_update_orb(oi_sats[i],_orb_state_ref[oi_sats[i]], _orb_state_ref[oi_sats[i]] + oi_intervals[i], state_init_name_sats[i], state_init_value_sats[i], _orb);
		}
		chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
		_oi_msec +=chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count(); 

	}

	void t_gfltproc::_update_orb(string sat, t_gtime beg, t_gtime end, vector<t_gpar::t_type> state_init_name, vector<double> state_init_value, t_gallorb * orb)
	{
		_orb_ref[sat] = beg;
		_epo[sat] = (_orb_state_ref[sat] - _orb_ref[sat])/_update_orb_interval[sat];
		map<int, vector<t_gtriple>> oi_result;
		Matrix phi;
		vector<double> state_newref;
		_phi_map[sat].clear();
		_phi_i_map[sat].clear();
		// intergration results stored in oi_result and _state_map
		_poi[sat]->integrate4srif(beg, end, state_init_name, state_init_value, _phi_map[sat], _phi_i_map[sat], oi_result, _state_map[sat]);
		// Here update the oi_result to the orb (orb will be used in the Equation construction)
		int oi_admout = dynamic_cast<t_gsetforce*>(_gset)->getOutStep(sat);
		orb->update(sat, oi_result, beg + oi_result.begin()->first * oi_admout, beg + oi_result.rbegin()->first* oi_admout, beg, oi_admout, false);
		
		// update estimated state_map
		for (int i = 0; i <= oi_result.rbegin()->first; i++) {
			_estimated_state_map[sat][beg + i * oi_admout] = _state_map[sat][beg + i * oi_admout];
		}

	

		double beta = _poi[sat]->get_beta_beg();
		if ( abs(beta) <= 4){
			_eclipse_marker[sat] = true;
			cout << "eclipse sat: " << sat  << setw(12) << beta << endl;
		}
		else{
			_eclipse_marker[sat] = false;
		}
	}



}