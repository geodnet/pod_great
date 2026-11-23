/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallamb.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambOw.h"
#include "gset/gsetamb.h"
#include "gset/gsetinp.h"
#include "gutils/gtime.h"
#include "gallamb.h"
#include "gall/gallbias.h"
#include "gambfix/gambcommon.h"
#include "gmodels/gprecisemodel.h"
#include "gutils/gstring.h"
#include "gutils/ginfolog.h"

#include "gio/gfile.h"
#include "gcoders/rinexo.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Eigen/Dense"
#include "Eigen/SVD"
#include <thread>
#include <chrono>
namespace great
{
	t_gallambs::t_gallambs()
	{
	}

	t_gallambs::t_gallambs(t_gsetbase* gset, t_glog* log)
	{
		this->_glog = log;
		this->_gset = gset;
		_dd_mode = dynamic_cast<t_gsetamb*>(gset)->dd_mode();
		_fix_dd_mode = dynamic_cast<t_gsetamb*>(gset)->fix_dd_mode();
		_frequency = dynamic_cast<t_gsetproc*>(gset)->frequency();
		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::QZS);

	}

	t_gallambs::~t_gallambs()
	{
		if (_gprecise)
		{
			delete _gprecise;
			_gprecise = nullptr;

		}
	}

	void t_gallambs::add_set(t_gsetbase* gset)
	{
		_gset = gset;
		_dd_mode = dynamic_cast<t_gsetamb*>(gset)->dd_mode();
		_fix_dd_mode = dynamic_cast<t_gsetamb*>(gset)->fix_dd_mode();


		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::QZS);
														  
		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::QZS);
	}

	void t_gallambs::add_log(t_glog* log)
	{
		_glog = log;
	}


	void t_gallambs::add_data(t_gallpar * pars, t_gallobs * obs, t_gallbias * allbias, t_gallobj * obj, t_gtime beg_time, t_gtime end_time)
	{
		string func_id = "t_gallambs::add_data";
		this->_gobs = obs;
		this->_gbias = allbias;
		this->_gpars = pars;
		this->add_site_crd(obj, beg_time);
		// ADD DATA
		t_gallproc* data = new t_gallproc();
		data->Add_Data(this->_gbias);
		_gprecise = new t_gprecisemodel(data, _glog, _gset);

		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_ows_table_from_pars");
		t_gallambs::_get_ows_table_from_pars();

		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_obs_table_from_gobs");
		if (_dd_mode != DD_MODEL::RAW_CB) t_gallambs::_get_obs_table_from_gobs(beg_time,end_time);
	}
	void t_gallambs::add_site_crd(t_gallobj* obj, t_gtime beg_time)
	{
		_gobj = obj;
		// get site receiver
		for (const auto& rec : _rec_list) {
			_site_rec[rec] = "";
			_rec_crd[rec] = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(rec);
			auto gobj = _gobj->obj(rec);
			if (gobj) {
				string receiver = dynamic_cast<t_grec*>(gobj.get())->rec(beg_time);
				_site_rec[rec] = receiver;
				if (!receiver.empty()) {
					auto n = receiver.find_first_of(" ");
					if (n != string::npos) _site_rec[rec] = receiver.substr(0, n);
					else _site_rec[rec] = receiver;
				}
				if (_rec_crd[rec].zero()) {
					t_gtriple crd(0.0, 0.0, 0.0);
					t_gtriple std(0.0, 0.0, 0.0);
					if (gobj) {
						gobj->get_recent_crd(beg_time, 0.1, crd, std);
						_rec_crd[rec] = crd;
					}
				}
			}
		}
	}

	void t_gallambs::update_data(t_gallpar* pars, t_gallobs* obs, t_gallbias* allbias, t_gallobj* obj, t_gtime beg_time, t_gtime end_time)
	{
		string func_id = "t_gallambs::update_data";

		// get_new data
		this->_gobs = obs;
		this->_gbias = allbias;
		this->_gpars = pars;
		this->_gobj = obj;

		// get MW value
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_obs_table_from_gobs");
		if (_dd_mode != DD_MODEL::RAW_CB) t_gallambs::_get_obs_table_from_gobs(beg_time, end_time);
	}

	void t_gallambs::clear_ambs_data()
	{
		this->_IDs.clear();
		this->_nOw_ambs.clear();
		this->_nSd_ambs.clear();
		this->_nSd_ambs_fix_WL.clear();
		this->_nSd_ambs_fix_NO.clear();
		this->_nDd_ambs.clear();
		this->_nDd_ambs_fix_WL.clear();
		this->_nDd_ambs_fix_NO.clear();

		this->fix_WL_num.clear();
		this->all_WL_num.clear();
		this->fix_NL_num.clear();
		this->all_NL_num.clear();
	}



	void t_gallambs::add_WL_bdeci(t_gbdeci* WL)
	{
		_evaluator_MW = WL;
	}

	void t_gallambs::add_NL_bdeci(t_gbdeci* NL)
	{
		_evaluator_NL = NL;
	}


	t_gallambs::t_dd t_gallambs::get_Dds()
	{
		return this->_nDd_ambs;
	}

	bool t_gallambs::cmb_amb(const AMB_COMB_OPTION& opt, const t_gamb::AMB_ID& ID, const vector<pair<string, string>>& baselines, double shortest_arc)
	{
		string func_id = "combine_Ow_ambs";

		if (baselines.empty())
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "the baselines is empty!");
			return false;
		}

		map<string,t_gmutex> add_mtx_table;
		if (opt == AMB_COMB_OPTION::CMB_SAT_SD){
			for (const auto& sys_key: this->_nOw_ambs){
				for (const auto& site: this->_rec_list){
					add_mtx_table[site];
					this->_nSd_ambs[sys_key.first][site];
				}
			}
		}

		if (opt == AMB_COMB_OPTION::CMB_REC_DD){
			for (const auto& sys_key: this->_nOw_ambs){
				for (const auto& bl: baselines){
					this->_nDd_ambs[sys_key.first][bl];
				}
			}
		}

		chrono::high_resolution_clock::time_point beg_t,end_t;
		if (opt== AMB_COMB_OPTION::CMB_REC_DD){
			beg_t = chrono::high_resolution_clock::now();
		}
		// TODO: add parallel
#pragma omp parallel for
		//for (const auto& iter : baselines)
		for (int i =0;i<baselines.size();i++)
		{
			const auto& iter = baselines[i];
			if (_dd_mode == DD_MODEL::IF_CB_WN)
			{
				if (
					ID == t_gamb::AMB_ID::IF_L1_L2 
					)
				{
					if (opt == AMB_COMB_OPTION::CMB_SAT_SD) this->_cmb_intersat_Sd_amb(iter, ID, shortest_arc, add_mtx_table);
					if (opt == AMB_COMB_OPTION::CMB_REC_DD) this->_cmb_interrec_Dd_amb(iter, ID, shortest_arc);
				}
			}	
			else continue;
		}
		
		if (opt == AMB_COMB_OPTION::CMB_REC_DD){
			end_t = chrono::high_resolution_clock::now();
		}

		if (opt == AMB_COMB_OPTION::CMB_REC_DD)
		{
			// get Dd value
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_Dds_value");

			if (_dd_mode == DD_MODEL::IF_CB_WN)
			{
				if (
					ID == t_gamb::AMB_ID::IF_L1_L2 
					)
				{

					chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
					_get_Dds_value(ID);
					chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();

					beg_t = chrono::high_resolution_clock::now();
					if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_select_Dds_value");
					_select_Dds_value(ID);
					if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_dependent_Dd");
					end_t = chrono::high_resolution_clock::now();
					if (_fix_dd_mode == t_gsetamb::FIX_DD_MODE::IND){
						beg_t = chrono::high_resolution_clock::now();
						_get_dependent_IF_Dd(ID);
						end_t = chrono::high_resolution_clock::now();
					}
					beg_t = chrono::high_resolution_clock::now();
					_get_fix_rate(ID);
					end_t = chrono::high_resolution_clock::now();
				}
			}
			else
			{
				return true;
			}
			
		}
		return true;
	}
	

	bool t_gallambs::_cmb_intersat_Sd_amb(const pair<string, string>& sat_baseline, const t_gamb::AMB_ID& ID, double shortest_arc,map<string,t_gmutex>& add_mtx_table)
	{
		string func_id = "_cmb_intersat_Sd_amb";

		if (sat_baseline.first.substr(0, 1) != sat_baseline.second.substr(0, 1))
		{
			if (this->_glog) this->_glog->logDebug(this->_class_id, func_id, "not the same system " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}
		t_dd_sys_index sys_key = make_pair(t_gsys::sat2gsys(sat_baseline.first), ID);
		if (this->_nOw_ambs.count(sys_key) < 1 || this->_nOw_ambs[sys_key].empty())
		{
			if (this->_glog) this->_glog->logDebug(this->_class_id, func_id, "no Ow " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}

		const auto& sys_Ow = this->_nOw_ambs[sys_key];
		for (const auto& iter_rec : this->_rec_list)
		{
			if (sys_Ow.count(iter_rec) < 1)
			{
				if (this->_glog) this->_glog->logDebug(this->_class_id, func_id, "no Ow " +
					iter_rec);
				continue;
			}
			const auto& Ow_Ambs = this->_nOw_ambs[sys_key][iter_rec];

			vector<int> satA;
			vector<int> satB;
			for (int amb_idx = 0; amb_idx < Ow_Ambs.size(); amb_idx++)
			{
				if (Ow_Ambs[amb_idx].sat == sat_baseline.first)  satA.push_back(amb_idx);
				if (Ow_Ambs[amb_idx].sat == sat_baseline.second) satB.push_back(amb_idx);
			}

			if (satA.empty() || satB.empty())
			{
				if (this->_glog) this->_glog->logDebug(this->_class_id, func_id, "not ambs for " +
					sat_baseline.first + " or " + sat_baseline.second);
				continue;
			}

			for (const auto& iter_satA : satA)
			{
				for (const auto& iter_satB : satB)
				{
					if (Ow_Ambs[iter_satA].rec != Ow_Ambs[iter_satB].rec) continue;

					const auto& amb_1st = Ow_Ambs[iter_satA];
					const auto& amb_2nd = Ow_Ambs[iter_satB];

					auto MAX_BEG = MAX_TIME(amb_1st.beg, amb_2nd.beg);
					auto MIN_END = MIN_TIME(amb_1st.end, amb_2nd.end);
					if (MAX_BEG >= MIN_END) continue;

					double cross_arc = MIN_END - MAX_BEG;
					if (round(cross_arc) < shortest_arc) continue;


					t_gamb_Sd tmp(iter_satA, iter_satB, ID);
					tmp.sgl_obj = Ow_Ambs[iter_satA].rec;
					tmp.dbl_obj_1st = Ow_Ambs[iter_satA].sat;
					tmp.dbl_obj_2nd = Ow_Ambs[iter_satB].sat;
					tmp.beg = (MAX_BEG);
					tmp.end = (MIN_END);

					string sat_rec = tmp.sgl_obj + "_" + tmp.dbl_obj_1st + "_" + tmp.dbl_obj_2nd;
					_glog->logDebug("t_gallambs", "_cmb_intersat_Sd_amb", "combine Sd for " + sat_rec);
					add_mtx_table[tmp.sgl_obj].lock();
					this->_nSd_ambs[sys_key][tmp.sgl_obj].push_back(tmp);
					add_mtx_table[tmp.sgl_obj].unlock();
				}
			}
		}
		return true;
	}

	bool t_gallambs::_cmb_intersat_Sd_amb(const pair<string, string>& sat_baseline, const t_gamb::AMB_ID& ID, double shortest_arc)
	{
		string func_id = "_cmb_intersat_Sd_amb";


		if (sat_baseline.first.substr(0, 1) != sat_baseline.second.substr(0, 1))
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "not the same system " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}

		t_dd_sys_index sys_key = make_pair(t_gsys::sat2gsys(sat_baseline.first), ID);
		if (this->_nOw_ambs.count(sys_key) < 1 || this->_nOw_ambs[sys_key].empty())
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "no Ow " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}

		const auto& sys_Ow = this->_nOw_ambs[sys_key];
		for (const auto& iter_rec : this->_rec_list)
		{
			if (sys_Ow.count(iter_rec) < 1)
			{
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "no Ow " +
					iter_rec);
				continue;
			}

			const auto& Ow_Ambs = this->_nOw_ambs[sys_key][iter_rec];

			vector<int> satA;
			vector<int> satB;
			for (int amb_idx = 0; amb_idx < Ow_Ambs.size(); amb_idx++)
			{
				if (Ow_Ambs[amb_idx].sat == sat_baseline.first)  satA.push_back(amb_idx);
				if (Ow_Ambs[amb_idx].sat == sat_baseline.second) satB.push_back(amb_idx);
			}

			if (satA.empty() || satB.empty())
			{
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "not ambs for " +
					sat_baseline.first + " or " + sat_baseline.second);
				continue;
			}

			for (const auto& iter_satA : satA)
			{
				for (const auto& iter_satB : satB)
				{
					if (Ow_Ambs[iter_satA].rec != Ow_Ambs[iter_satB].rec) continue;

					const auto& amb_1st = Ow_Ambs[iter_satA];
					const auto& amb_2nd = Ow_Ambs[iter_satB];

					auto MAX_BEG = MAX_TIME(amb_1st.beg, amb_2nd.beg);
					auto MIN_END = MIN_TIME(amb_1st.end, amb_2nd.end);
					if (MAX_BEG >= MIN_END) continue;

					double cross_arc = MIN_END - MAX_BEG;
					if (round(cross_arc) < shortest_arc) continue;

					t_gamb_Sd tmp(iter_satA, iter_satB, ID);
					tmp.sgl_obj = Ow_Ambs[iter_satA].rec;
					tmp.dbl_obj_1st = Ow_Ambs[iter_satA].sat;
					tmp.dbl_obj_2nd = Ow_Ambs[iter_satB].sat;
					tmp.beg = (MAX_BEG);
					tmp.end = (MIN_END);

					string sat_rec = tmp.sgl_obj + "_" + tmp.dbl_obj_1st + "_" + tmp.dbl_obj_2nd;
					_glog->logDebug("t_gallambs", "_cmb_intersat_Sd_amb", "combine Sd for " + sat_rec);
					this->_nSd_ambs[sys_key][tmp.sgl_obj].push_back(tmp);
				}
			}
		}
		return true;
	}

	bool t_gallambs::_cmb_interrec_Dd_amb(const pair<string, string>& rec_baseline, const t_gamb::AMB_ID& ID, double shortest_arc)
	{
		string func_id = "_cmb_interrec_Dd_amb";

		for (const auto& iter_key : _nSd_ambs)
		{
			const auto& key = iter_key.first;
			const auto& Ambs = iter_key.second;
			if (key.second != ID) continue;

			if ((Ambs.count(rec_baseline.first) < 1) ||
				(Ambs.count(rec_baseline.second) < 1))
				return false;

			const auto& Sd_1st = this->_nSd_ambs[key][rec_baseline.first];
			const auto& Sd_2nd = this->_nSd_ambs[key][rec_baseline.second];

			if (Sd_1st.empty() || Sd_2nd.empty()) return false;

			// GLONASS should be the same type receiver
			if (iter_key.first.first == GSYS::GLO) {
				string rec1 = rec_baseline.first;
				string rec2 = rec_baseline.second;

				bool same_rectp = (!_site_rec.at(rec1).empty() && _site_rec.at(rec1) == _site_rec.at(rec2));

				if (!same_rectp) {
					return false;
				}
			}

			const auto& Ow_1st = this->_nOw_ambs[key][rec_baseline.first];
			const auto& Ow_2nd = this->_nOw_ambs[key][rec_baseline.second];

			for (const auto& iter_1st : Sd_1st)
			{
				const auto& idx1 = (iter_1st.idx2ow_1st.find(1) != iter_1st.idx2ow_1st.end()) ? iter_1st.idx2ow_1st.find(1)->second : -1;
				const auto& idx2 = (iter_1st.idx2ow_1st.find(2) != iter_1st.idx2ow_1st.end()) ? iter_1st.idx2ow_1st.find(2)->second : -1;
				if (idx1 == -1 || idx2 == -1) continue;
				const auto& line_1st = make_pair(Ow_1st[idx1].sat, Ow_1st[idx2].sat);

				for (const auto& iter_2nd : Sd_2nd)
				{
					const auto& idx3 = (iter_2nd.idx2ow_1st.find(1) != iter_2nd.idx2ow_1st.end()) ? iter_2nd.idx2ow_1st.find(1)->second : -1;
					const auto& idx4 = (iter_2nd.idx2ow_1st.find(2) != iter_2nd.idx2ow_1st.end()) ? iter_2nd.idx2ow_1st.find(2)->second : -1;
					if (idx3 == -1 || idx4 == -1) continue;
					const auto& line_2nd = make_pair(Ow_2nd[idx3].sat, Ow_2nd[idx4].sat);

					if (line_1st != line_2nd) continue;

					const auto& amb_1st = iter_1st;
					const auto& amb_2nd = iter_2nd;

					auto MAX_BEG = MAX_TIME(amb_1st.beg, amb_2nd.beg);
					auto MIN_END = MIN_TIME(amb_1st.end, amb_2nd.end);

					if (MAX_BEG >= MIN_END) continue;

					double cross_arc = MIN_END - MAX_BEG;
					if (round(cross_arc) < shortest_arc) continue;
		

					pair<string, string> idx_key(Ow_1st[idx1].rec, Ow_2nd[idx3].rec);

					t_gamb_Dd tmp(ID, idx1, idx2, idx3, idx4);
					tmp.beg = MAX_BEG;
					tmp.end = MIN_END;

					tmp.sat_1st = Ow_1st[idx1].sat;
					tmp.sat_2nd = Ow_1st[idx2].sat;
					tmp.rec_1st = Ow_1st[idx1].rec;
					tmp.rec_2nd = Ow_2nd[idx3].rec;

					string sat_rec = tmp.sat_1st + "_" + tmp.sat_2nd + "_" + tmp.rec_1st + "_" + tmp.rec_2nd;
					_glog->logDebug("t_gallambs", "_cmb_intersat_Sd_amb", "combine Sd for " + sat_rec);
					this->_nDd_ambs[key][idx_key].push_back(tmp);

				}
			}
		}
		return true;

		//
	}

	void t_gallambs::_get_ows_table_from_pars()
	{
		if (this->_gpars == nullptr) return;
		for (unsigned int i = 0; i < _gpars->parNumber(); ++i)
		{
			const auto& par = (*_gpars)[i];
			//bds geo do not paticipate in ambfix
			if (t_gsys::bds_geo(par.prn)) continue;
			if (_rec_list.find(par.site) == _rec_list.end()) _rec_list.insert(par.site);
			if (_sat_list.find(par.prn) == _sat_list.end()) _sat_list.insert(par.prn);

			t_gamb::AMB_ID ID = t_gamb::AMB_ID::UNDEF;
			switch (par.parType)
			{
			case t_gpar::AMB_IF:  {ID = t_gamb::AMB_ID::IF_L1_L2; break;}
			default:  break;
			}
			if (_IDs.find(ID) == _IDs.end()) _IDs.insert(ID);

			t_gamb_Ow amb(ID, par.prn, par.site, par.beg, par.end);
			amb.value = par.value();
			amb.sigma = par.apriori();
			auto gsys = t_gsys::sat2gsys(par.prn);
			auto key = make_pair(gsys, ID);
			this->_nOw_ambs[key][par.site].push_back(amb);
		}


	}

	void t_gallambs::_get_obs_table_from_gobs()
	{
		t_gtime proc_beg = dynamic_cast<t_gsetgen*>(_gset)->beg();
		t_gtime proc_end = dynamic_cast<t_gsetgen*>(_gset)->end();

		_get_obs_table_from_gobs(proc_beg, proc_end);
	}

	void t_gallambs::_get_obs_table_from_gobs(t_gtime beg, t_gtime end)
	{
		if (_IDs.empty()) return;
		if (!this->_gbias) return;
		if (!this->_gobs) return;
		if (!this->_gset) return;

		bool isRawAll = false;

		double  proc_len = 0.0;
		t_gtime proc_beg = beg;
		t_gtime proc_end = end;
		set<string>  rec_list;
		try
		{
			rec_list = dynamic_cast<t_gsetgen*>(_gset)->recs();
			proc_len = dynamic_cast<t_gsetamb*>(_gset)->wl_interval();
		}
		catch (...)
		{
			return;
		}
		if (double_eq(proc_len, 0.0)) return;

		if (rec_list.empty()) return;

		time_table.clear();
		auto tmp_now = proc_beg ;
		auto tmp_end = (proc_end + 2 * proc_len);
		do
		{
			time_table.insert(tmp_now);
			tmp_now = tmp_now + proc_len;
		} while (tmp_now <= tmp_end);

		if (time_table.empty()) return;

		// insert all sys rec;
		// insert all obs file name;
		vector<string> v_rec;
		vector<map<t_gtime,vector<t_gobsgnss*> > > rec_obs(rec_list.size());
		vector< map<t_dd_sys_index,vector<t_gamb_Ow>*> > rec_amb_for_sys(rec_list.size());
		vector< map<t_dd_sys_index,map<string,MW_OF_AMB*> > >  rec_results(rec_list.size());
		int count =0;
		for (const string&rec : rec_list) {
			for (auto& ambs_for_sys : this->_nOw_ambs) {
				//this->_MW_used[ambs_for_sys.first.first][rec];
				for (auto& sat : _sat_list) {
					rec_results[count][ambs_for_sys.first][sat] = &(this->_MW_mean_sigma[ambs_for_sys.first][rec][sat]);
				}
				rec_amb_for_sys[count][ambs_for_sys.first]=&(ambs_for_sys.second[rec]);
			}
			v_rec.push_back(rec);
			for (const auto& time: time_table){
				rec_obs[count][time]=_gobs->obs_lockfree(true,rec,time);
			}
			count++;
		}
		
#pragma omp parallel for
		for (int i =0;i< v_rec.size();i++){
			this->_get_obs_table_from_gobs_one_site(v_rec[i],isRawAll,rec_amb_for_sys[i],time_table,rec_obs[i],rec_results[i]);
		}

	}

	void t_gallambs::_get_obs_table_from_gobs_one_site(string rec,bool isRawAll,const map<t_dd_sys_index,vector<t_gamb_Ow>*>& ambs_for_sys,const set<t_gtime> local_time_table,const map<t_gtime,vector<t_gobsgnss*>>& rec_obs, map<t_dd_sys_index,map<string,MW_OF_AMB*>>& result){

		for (const auto& ambs_for_rec: ambs_for_sys) // cycle sys
		{
			// ==================================================================================================
			
			t_gtime min = LAST_TIME;
			t_gtime max = FIRST_TIME;
			for (const auto& amb : *(ambs_for_rec.second))
			{
				if (min > amb.beg) min = amb.beg;
				if (max < amb.end) max = amb.end;
			}
			auto it_min = local_time_table.lower_bound(min);
			auto it_max = local_time_table.lower_bound(max);
			if (it_min == local_time_table.end()) it_min = local_time_table.begin();
			if (it_max == local_time_table.end()) it_max = local_time_table.end();

			map<t_gtime, set<string>> sat_used;
			for (const auto& amb : *(ambs_for_rec.second)) // cycle sat
			{
				if (amb.beg > amb.end) continue;
				auto beg = local_time_table.lower_bound(amb.beg);
				auto end = local_time_table.lower_bound(amb.end);
				if (beg == local_time_table.end()) continue;
				if (end != local_time_table.end()) end++;

				auto now = beg;
				for (now; now != end; now++) { sat_used[*now].insert(amb.sat); }

				// reinit MW obj if amb is new
				auto mw_iter = result[ambs_for_rec.first].find(amb.sat);
				if (mw_iter==result[ambs_for_rec.first].end()) continue;
				auto& MW_obj = *(mw_iter->second);
				if (MW_obj.beg_time != amb.beg) {
					// new init
					MW_obj = MW_OF_AMB();
					MW_obj.beg_time = amb.beg;
				}
			}

			auto loop_now = it_min;
			auto loop_end = (it_max == local_time_table.end()) ? it_max : (++it_max);

			for (loop_now; loop_now != loop_end; loop_now++) // cycle time
			{
				//vector<t_gobsgnss*> sat_datas = allobs_ptr->obs_lockfree(true, rec, *loop_now);
				const auto& sat_data_iter = rec_obs.find(*loop_now);
				if (sat_data_iter == rec_obs.end()) continue;
				vector<t_gobsgnss*> sat_datas = sat_data_iter->second;
				if (sat_datas.empty()) continue;
				if (sat_used.empty())  continue;

				for (auto& sat_iter : sat_datas) // cycle sat
				{
					auto sat = sat_iter->sat();
					auto sys = t_gsys::sat2gsys(sat);

					if (sat_used[*loop_now].count(sat) <= 0) continue;

					// update MW mean and sigma
					GOBSBAND b1 = BAND;
					GOBSBAND b2 = BAND;
					int     freq = -1;
					auto amb_ID = ambs_for_rec.first.second;
					if ( amb_ID == t_gamb::AMB_ID::IF_L1_L2) {  freq = 2; b1 = _band_index[sys][FREQ_1]; b2 = _band_index[sys][FREQ_2]; }

					
					t_gobs gobsP1 = t_gobs(sat_iter->select_range(b1, isRawAll));
					t_gobs gobsL1 = t_gobs(sat_iter->select_phase(b1, isRawAll));
					t_gobs gobsP2 = t_gobs(sat_iter->select_range(b2, isRawAll));
					t_gobs gobsL2 = t_gobs(sat_iter->select_phase(b2, isRawAll));

					if (sat_iter->obs_empty()) continue;
					sat_iter->apply_bias(*_gbias);
					double mw_obs_12 = 0.0;
					mw_obs_12 = sat_iter->MW_cycle(gobsL1, gobsL2, gobsP1, gobsP2);
					if (double_eq(mw_obs_12, 0.0)) continue;

					// update MW mean and sigma;
					auto mw_iter = result[ambs_for_rec.first].find(sat);
					if (mw_iter== result[ambs_for_rec.first].end()) continue;
					mw_iter->second->add_data(mw_obs_12);
				}// finsh time
			} // finsh time loop
		}
	}

	void t_gallambs::_get_dependent_IF_Dd(const t_gamb::AMB_ID& ID)
	{
		string func_id = "_get_dependent_Dd";
		vector<pair<GSYS, t_gamb::AMB_ID> > sys_list;
		for (auto&sys : this->_nDd_ambs){
			sys_list.push_back(sys.first);
		}
		// get the sys
#pragma omp parallel for
		for (int i =0;i< sys_list.size();i++)
		{
			auto  sys_first = sys_list[i];	
			auto&  sys_second = this->_nDd_ambs[sys_first];
			if (sys_first.second != ID) continue;
			if (sys_second.empty()) continue;
			// get the bl
			// TODO: check the math dependent together with all baselines
			int nAmb_total = 0;
			unordered_map<string,int> site_table;
			int count = 0;
			for (auto& site_amb: this->_nOw_ambs[sys_first]){
				site_table[site_amb.first] = count;
				nAmb_total += site_amb.second.size();
				count+= site_amb.second.size();
			}

			int ndef = 0;
			vector<double> vector_space((nAmb_total+100)*nAmb_total,0);
			vector<double> crt_vector(nAmb_total+100,0);

			// WL && NL Fix
			for (auto& bl : sys_second){
				vector<t_gamb_Dd> result;
				int num_Ow_1st = site_table[bl.first.first];
				int num_Ow_2nd = site_table[bl.first.second];

				math_dependent_IF_Dd(bl.second,num_Ow_1st,num_Ow_2nd,nAmb_total,result,vector_space,crt_vector,ndef);

				bl.second = result;	
			}

			// Only WL Fix
			for (auto& bl : sys_second){
				vector<t_gamb_Dd> result;
				int num_Ow_1st = site_table[bl.first.first];
				int num_Ow_2nd = site_table[bl.first.second];

				auto& bl_ambs_fix_WL = this->_nDd_ambs_fix_WL[sys_first][bl.first];
				math_dependent_IF_Dd(bl_ambs_fix_WL,num_Ow_1st,num_Ow_2nd,nAmb_total,result,vector_space,crt_vector,ndef);

				bl_ambs_fix_WL = result;
			}

			// No Fix
			for (auto& bl : sys_second){
				vector<t_gamb_Dd> result;
				int num_Ow_1st = site_table[bl.first.first];
				int num_Ow_2nd = site_table[bl.first.second];

				auto& bl_ambs_fix_NO = this->_nDd_ambs_fix_NO[sys_first][bl.first];
				math_dependent_IF_Dd(bl_ambs_fix_NO,num_Ow_1st,num_Ow_2nd,nAmb_total,result,vector_space,crt_vector,ndef);

				bl_ambs_fix_NO = result;
			}

		}
	}


	void t_gallambs::_get_fix_rate(const t_gamb::AMB_ID& ID)
	{

		// fix ALL
		for (auto& sys : this->_nDd_ambs)
		{
			const auto& amb_id = sys.first.second;
			const auto& sys_key = sys.first;
			if (amb_id != ID) continue;
			for (auto& bl : sys.second)
			{
				const auto& rec_1st = bl.first.first;
				const auto& rec_2nd = bl.first.second;

				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
					all_WL_num[sys.first] = all_WL_num[sys.first] + 1;

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}

		// fix WL
		for (auto& sys : this->_nDd_ambs_fix_WL)
		{
			const auto& amb_id = sys.first.second;
			const auto& sys_key = sys.first;
			if (amb_id != ID) continue;
			for (auto& bl : sys.second)
			{
				const auto& rec_1st = bl.first.first;
				const auto& rec_2nd = bl.first.second;

				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					if (!double_eq(amb.NL_value, 0.0) && !double_eq(amb.MW_value, 0.0))
					{
						all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
						all_WL_num[sys.first] = all_WL_num[sys.first] + 1;
					}

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}

		// fix NO
		for (auto& sys : this->_nDd_ambs_fix_NO)
		{
			const auto& amb_id = sys.first.second;
			const auto& sys_key = sys.first;
			if (amb_id != ID) continue;
			for (auto& bl : sys.second)
			{
				const auto& rec_1st = bl.first.first;
				const auto& rec_2nd = bl.first.second;

				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					if (!double_eq(amb.NL_value, 0.0) && !double_eq(amb.MW_value, 0.0))
					{
						all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
						all_WL_num[sys.first] = all_WL_num[sys.first] + 1;
					}

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}
	}

	void t_gallambs::math_dependent_IF_Dd(vector<t_gamb_Dd>& input, int num_Ow_1st, int num_Ow_2nd, int num_total, vector<t_gamb_Dd>& result, vector<double>& vector_space, vector<double>& crt_vector,int& ndef)
	{
		string func_id = "math_dependent_IF_Dd";
		
		// check the depend

		for (auto& amb : input)
		{
			int index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second + 1) : -1;
			int index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second + 1) : -1;
			int index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second + 1) : -1;
			int index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second + 1) : -1;

			if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) continue;

			vector<int> index_locs = {index1+num_Ow_1st,index2+num_Ow_1st,index3+num_Ow_2nd,index4+num_Ow_2nd};
			if (_check_amb_depend(vector_space.data(),crt_vector.data(),num_total,4,index_locs,ndef)){
				continue;
			}
			else{
				result.push_back(amb);
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
					int2str(result.size() + 1) + " " +
					amb.rec_1st + " " +
					amb.rec_2nd + " " +
					amb.sat_1st + " " +
					amb.sat_2nd);	
			}

			
		}
	}

	void t_gallambs::math_dependent_IF_Dd(vector<t_gamb_Dd>& input, int num_Ow_1st, int num_Ow_2nd, vector<t_gamb_Dd>& result, Matrix& vector_space)
	{
		string func_id = "math_dependent_IF_Dd";
		// check the depend
		int ncol = num_Ow_1st + num_Ow_2nd;

		for (auto& amb : input)
		{
			int index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second + 1) : -1;
			int index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second + 1) : -1;
			int index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second + 1) : -1;
			int index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second + 1) : -1;

			if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) continue;

			if (vector_space.is_zero())
			{
				vector_space.resize(1, ncol);
				vector_space = 0.0;
				vector_space(1, index1) = 1;
				vector_space(1, index2) = -1;
				vector_space(1, index3 + num_Ow_1st) = -1;
				vector_space(1, index4 + num_Ow_1st) = 1;
				result.push_back(amb);

				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "Base line is : " +
					int2str(result.size() + 1) + " " +
					amb.rec_1st + " " +
					amb.rec_2nd + " " +
					amb.sat_1st + " " +
					amb.sat_2nd);
				continue;
			}

			RowVector new_vector(0.0, ncol);
			new_vector(index1) = 1;
			new_vector(index2) = -1;
			new_vector(index3 + num_Ow_1st) = -1;
			new_vector(index4 + num_Ow_1st) = 1;

			Matrix A = vector_space & new_vector;
			DiagonalMatrix D;
			if (A.nrows() < A.ncols())
			{
				SVD(A.t(), D);
			}
			else
			{
				SVD(A, D);
			}
			bool isIndepend = (fabs(D(D.nrows(), D.ncols())) > 0.01);

			if (isIndepend)
			{
				vector_space = A;
				result.push_back(amb);
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
					int2str(result.size() + 1) + " " +
					amb.rec_1st + " " +
					amb.rec_2nd + " " +
					amb.sat_1st + " " +
					amb.sat_2nd);
			}
		}
	}

	void t_gallambs::_select_Dds_value(const t_gamb::AMB_ID& ID)
	{
		for (auto& sys : this->_nDd_ambs)
		{
			const auto& amb_id = sys.first.second;
			const auto& sys_key = sys.first;
			if (amb_id != ID) continue;
			vector<decltype(sys.second.begin())> vmap;
			vector<vector<t_gamb_Dd>*> local_Dd_ambs_fix_WL;
			vector<vector<t_gamb_Dd>*> local_Dd_ambs_fix_NO;
			for (auto iter=sys.second.begin();iter!=sys.second.end();iter++){
				vmap.push_back(iter);
				auto& bl = *iter;
				local_Dd_ambs_fix_WL.push_back(&(this->_nDd_ambs_fix_WL[sys.first][bl.first]));
				local_Dd_ambs_fix_NO.push_back(&(this->_nDd_ambs_fix_NO[sys.first][bl.first]));
			}

#pragma omp parallel for
			//for (auto& bl : sys.second)
			for (int i=0;i<vmap.size();i++)
			{
				auto& bl = *(vmap[i]);
				const auto& rec_1st = bl.first.first;
				const auto& rec_2nd = bl.first.second;

				if (bl.second.empty()) continue;
				vector<t_gamb_Dd> result;
				vector<t_gamb_Dd> result_fix_WL;
				vector<t_gamb_Dd> result_fix_NO;

				for (auto& amb : bl.second)
				{		
					if (amb.NL_fixed && amb.MW_fixed) result.push_back(amb);

					if (amb.MW_fixed && !amb.NL_fixed && !double_eq(amb.MW_value,0.0) && !double_eq(amb.NL_value,0.0))
					{
						result_fix_WL.push_back(amb);
					}

					if (!amb.MW_fixed && !amb.NL_fixed && !double_eq(amb.MW_value, 0.0) && !double_eq(amb.NL_value, 0.0))
					{
						result_fix_NO.push_back(amb);
					}
				}

				sort(result_fix_WL.begin(), result_fix_WL.end(), sort_Dd);
				sort(result_fix_NO.begin(), result_fix_NO.end(), sort_Dd);
				*(local_Dd_ambs_fix_WL[i]) = result_fix_WL;
				*(local_Dd_ambs_fix_NO[i]) = result_fix_NO;
				if (result.empty())
				{
					vector<t_gamb_Dd>().swap(bl.second);
					continue;
				}
				sort(result.begin(), result.end(), sort_Dd);	
				bl.second = result;
			}
		}
	}



	void t_gallambs::_get_Dds_value(const t_gamb::AMB_ID& ID)
	{
		for (auto& sys : this->_nDd_ambs)
		{
			const auto& sys_id_key = sys.first;
			const auto& sys_key    = sys.first.first;
			const auto& id_key     = sys.first.second;
			if (id_key == t_gamb::AMB_ID::UNDEF) continue;
			if (id_key != ID) continue;
			vector<decltype(sys.second.begin())> vmap;
			for (auto iter=sys.second.begin();iter!=sys.second.end();iter++)
			{
				vmap.push_back(iter);
			}
#pragma omp parallel for
			// cycle site-site baselines
			for (int i =0;i<vmap.size();i++)
			{
				auto& bl_ambs = *vmap[i];
				const auto& rec_1st = bl_ambs.first.first;
				const auto& rec_2nd = bl_ambs.first.second;
				const auto& bl = bl_ambs.first;

				// cycle Dd ambs on this baseline 
				for (auto& amb : bl_ambs.second)
				{
					if (
						ID == t_gamb::AMB_ID::IF_L1_L2 
						)
					{
						if (!_get_IF_Dd_value(sys_id_key, bl, amb)) continue;
						if (!_get_MW_value(sys_id_key, bl, amb)) continue;
						if (!_get_NL_value(sys_id_key, bl, amb)) continue;
					}
					else
					{
						throw logic_error("wrong type of ambs while _get_Dds_value");
					}
				}// loop baseline
				if (this->_fix_dd_mode==t_gsetamb::FIX_DD_MODE::IND){
					// sort ambs according to wl_nl_fix wl_fix and nofix
					vector<t_gamb_Dd> result_fix_WLNL;
					vector<t_gamb_Dd> result_fix_WL;
					vector<t_gamb_Dd> result_fix_NO;
					for (auto& amb :bl_ambs.second)
					{		
						if (amb.NL_fixed && amb.MW_fixed) result_fix_WLNL.push_back(amb);

						if (amb.MW_fixed && !amb.NL_fixed && !double_eq(amb.MW_value,0.0) && !double_eq(amb.NL_value,0.0))
						{
							result_fix_WL.push_back(amb);
						}

						if (!amb.MW_fixed && !amb.NL_fixed && !double_eq(amb.MW_value, 0.0) && !double_eq(amb.NL_value, 0.0))
						{
							result_fix_NO.push_back(amb);
						}
					}
					// check independence		
					int nAmb_total = 0;
					unordered_map<string,int> site_table;
					int count = 0;
					vector<string> bl_sites={rec_1st,rec_2nd};
					for (const auto& site_amb: bl_sites){
						site_table[site_amb] = count;
						nAmb_total += this->_nOw_ambs[sys.first][site_amb].size();
						count+= this->_nOw_ambs[sys.first][site_amb].size();
					}

					int ndef = 0;
					vector<double> vector_space((nAmb_total+100)*nAmb_total,0);
					vector<double> crt_vector(nAmb_total+100,0);

					// WL && NL Fix
					vector<t_gamb_Dd> result;
					int num_Ow_1st = site_table[rec_1st];
					int num_Ow_2nd = site_table[rec_2nd];

					math_dependent_IF_Dd(result_fix_WLNL,num_Ow_1st,num_Ow_2nd,nAmb_total,result,vector_space,crt_vector,ndef);
					result_fix_WLNL= result;	

					// Only WL Fix
					result.clear();
					math_dependent_IF_Dd(result_fix_WL,num_Ow_1st,num_Ow_2nd,nAmb_total,result,vector_space,crt_vector,ndef);

					result_fix_WL = result;

					// No Fix
					result.clear();
					math_dependent_IF_Dd(result_fix_NO,num_Ow_1st,num_Ow_2nd,nAmb_total,result,vector_space,crt_vector,ndef);
					result_fix_NO = result;

					// reassign the nDd_ambs
					bl_ambs.second.clear();
					bl_ambs.second.insert(bl_ambs.second.end(),result_fix_WLNL.begin(),result_fix_WLNL.end());
					bl_ambs.second.insert(bl_ambs.second.end(),result_fix_WL.begin(),result_fix_WL.end());
					bl_ambs.second.insert(bl_ambs.second.end(),result_fix_NO.begin(),result_fix_NO.end());
				}
			}// loop sys
		} // loop all
	}



	bool t_gallambs::_get_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		if (
			index_key.second == t_gamb::AMB_ID::IF_L1_L2 
			)
		{
			if (!_get_IF_NL_value(index_key, bl, amb)) return false;
		}
		else
		{
			return false;
		}

		return true;
	}

	bool t_gallambs::_get_IF_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;

		if (double_eq(amb.value, 0.0)) return false;
		if (double_eq(amb.MW_value, 0.0)) return false;

		t_gobsgnss gnss;
		gnss.sat(amb.sat_1st);
		auto gsys = t_gsys::sat2gsys(amb.sat_1st);

		double freq_1st = 0.0;
		double freq_2nd = 0.0;

		if (index_key.second == t_gamb::AMB_ID::IF_L1_L2)
		{
			freq_1st = gnss.frequency(_band_index[gsys][FREQ_1]);
			freq_2nd = gnss.frequency(_band_index[gsys][FREQ_2]);
		}
		else
		{
			// no support for others
			return false;
		}

		if (index_key.first == GSYS::GLO) {
			// fix ambiguities in GLONASS systems

			// for a given channel number k, the ratio of the L1 and L2 frequencies attains a fixed value of
			string sat1 = amb.sat_1st;
			string sat2 = amb.sat_2nd;
			t_gobsgnss gnss1,gnss2;
			gnss1.sat(sat1);
			gnss1.channel(_glofrq_num.at(sat1));
			gnss2.sat(sat2);
			gnss2.channel(_glofrq_num.at(sat2));

			double freq_11 = gnss1.frequency(_band_index[GSYS::GLO][FREQ_1]);
			double freq_12 = gnss1.frequency(_band_index[GSYS::GLO][FREQ_2]);
			double freq_21 = gnss2.frequency(_band_index[GSYS::GLO][FREQ_1]);
			double freq_22 = gnss2.frequency(_band_index[GSYS::GLO][FREQ_2]);

			int index1 = -1;
			int index2 = -1;
			int index3 = -1;
			int index4 = -1;

			if (index_key.second == t_gamb::AMB_ID::IF_L1_L2 )
			{
				index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second) : -1;
				index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second) : -1;
				index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second) : -1;
				index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second) : -1;
			}
			else
			{
				index1 = (amb.idx2ow_2nd.find(1) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(1)->second) : -1;
				index2 = (amb.idx2ow_2nd.find(2) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(2)->second) : -1;
				index3 = (amb.idx2ow_2nd.find(3) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(3)->second) : -1;
				index4 = (amb.idx2ow_2nd.find(4) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(4)->second) : -1;
			}

			if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) return false;

			const auto& amb1 = this->_nOw_ambs[index_key][rec_1st][index1];
			const auto& amb2 = this->_nOw_ambs[index_key][rec_1st][index2];
			const auto& amb3 = this->_nOw_ambs[index_key][rec_2nd][index3];
			const auto& amb4 = this->_nOw_ambs[index_key][rec_2nd][index4];

			// single-differenced ambiguity between stations for sat1 and sat2
			double sd_amb_s1 = amb1.value - amb3.value;
			double sd_amb_s2 = amb2.value - amb4.value;
			
			amb.NL_value = (freq_11+freq_12) / CLIGHT * sd_amb_s1 - (freq_21 +freq_22) / CLIGHT * sd_amb_s2 -3.5 * round(amb.MW_value);
			amb.NL_sigma = amb.sigma * (freq_11 + freq_12) / CLIGHT;
		}
		else {

			double length_1st = CLIGHT / freq_1st;
			double length_2nd = CLIGHT / freq_2nd;

			double length_3rd = length_1st * length_2nd / (length_2nd - length_1st);
			// length_4th = c / (f1+f2)
			double length_4th = length_1st * length_2nd / (length_2nd + length_1st);

			amb.NL_value = amb.value / (length_4th)-round(amb.MW_value) * length_3rd / length_2nd;
			amb.NL_sigma = amb.sigma / length_4th;

		}


		//if (double_eq(0.0, amb.NL_sigma)) amb.NL_sigma = 0.05;
		amb.NL_sigma = 0.05;
		amb.NL_fixed = (_evaluator_NL) ? _evaluator_NL->bdeci(amb.NL_value, amb.NL_sigma, 1) : false;

		amb.NL_fixed = (amb.MW_fixed && amb.NL_fixed);
		return true;
	}


	bool t_gallambs::_get_MW_value(const pair<GSYS, t_gamb::AMB_ID>& index_key_input, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;
		
		pair<GSYS, t_gamb::AMB_ID> index_key;
		if (
			index_key_input.second == t_gamb::AMB_ID::IF_L1_L2 )
		{
			index_key = index_key_input;
		}
		else throw logic_error("error amb ID");
			
		int index1 = -1;
		int index2 = -1;
		int index3 = -1;
		int index4 = -1;

		if (
			index_key_input.second == t_gamb::AMB_ID::IF_L1_L2 
			)
		{
			index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second) : -1;
			index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second) : -1;
			index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second) : -1;
			index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second) : -1;
		}
		else
		{
			throw logic_error("The amb ID is Wrong, it should be IF_L1_L2, ALL_L1_L2 and ALL_L1_L3.");
		}

		if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) return false;

		const auto& amb1 = this->_nOw_ambs[index_key][rec_1st][index1];
		const auto& amb2 = this->_nOw_ambs[index_key][rec_1st][index2];
		const auto& amb3 = this->_nOw_ambs[index_key][rec_2nd][index3];
		const auto& amb4 = this->_nOw_ambs[index_key][rec_2nd][index4];

		auto& MW_obj1 = this->_MW_mean_sigma[index_key][amb1.rec][amb1.sat];
		auto& MW_obj2 = this->_MW_mean_sigma[index_key][amb2.rec][amb2.sat];
		auto& MW_obj3 = this->_MW_mean_sigma[index_key][amb3.rec][amb3.sat];
		auto& MW_obj4 = this->_MW_mean_sigma[index_key][amb4.rec][amb4.sat];

		double mw_value_1st =0.0; double mw_sigma_1st = 0.0;
		double mw_value_2nd =0.0; double mw_sigma_2nd = 0.0;
		double mw_value_3rd =0.0; double mw_sigma_3rd = 0.0;
		double mw_value_4th =0.0; double mw_sigma_4th = 0.0;

		MW_obj1.clean_mean_sigma(mw_value_1st, mw_sigma_1st);
		MW_obj2.clean_mean_sigma(mw_value_2nd, mw_sigma_2nd);
		MW_obj3.clean_mean_sigma(mw_value_3rd, mw_sigma_3rd);
		MW_obj4.clean_mean_sigma(mw_value_4th, mw_sigma_4th);


		if (double_eq(mw_value_1st, 0.0)) return false;
		if (double_eq(mw_value_2nd, 0.0)) return false;
		if (double_eq(mw_value_3rd, 0.0)) return false;
		if (double_eq(mw_value_4th, 0.0)) return false;

		
		amb.MW_value = mw_value_1st -
			mw_value_2nd -
			mw_value_3rd +
			mw_value_4th;

		amb.MW_sigma = sqrt(mw_sigma_1st * mw_sigma_1st +
			mw_sigma_2nd * mw_sigma_2nd +
			mw_sigma_3rd * mw_sigma_3rd +
			mw_sigma_4th * mw_sigma_4th);

		amb.Mw_value[1] = mw_value_1st;
		amb.Mw_value[2] = mw_value_2nd;
		amb.Mw_value[3] = mw_value_3rd;
		amb.Mw_value[4] = mw_value_4th;
		amb.Mw_sigma[1] = mw_sigma_1st;
		amb.Mw_sigma[2] = mw_sigma_2nd;
		amb.Mw_sigma[3] = mw_sigma_3rd;
		amb.Mw_sigma[4] = mw_sigma_4th;


		if (double_eq(0.0, amb.MW_sigma)) amb.MW_sigma = 999.0;
		amb.MW_fixed = (_evaluator_MW) ? _evaluator_MW->bdeci(amb.MW_value, amb.MW_sigma, 1) : false;

		return true;
	}

	
	bool t_gallambs::_get_IF_Dd_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;

		int index1 = -1;
		int index2 = -1;
		int index3 = -1;
		int index4 = -1;

		if (index_key.second == t_gamb::AMB_ID::IF_L1_L2 )
		{
			index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second) : -1;
			index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second) : -1;
			index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second) : -1;
			index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second) : -1;
		}
		else
		{
			index1 = (amb.idx2ow_2nd.find(1) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(1)->second) : -1;
			index2 = (amb.idx2ow_2nd.find(2) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(2)->second) : -1;
			index3 = (amb.idx2ow_2nd.find(3) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(3)->second) : -1;
			index4 = (amb.idx2ow_2nd.find(4) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(4)->second) : -1;
		}

		if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) return false;

		const auto& amb1 = this->_nOw_ambs[index_key][rec_1st][index1];
		const auto& amb2 = this->_nOw_ambs[index_key][rec_1st][index2];
		const auto& amb3 = this->_nOw_ambs[index_key][rec_2nd][index3];
		const auto& amb4 = this->_nOw_ambs[index_key][rec_2nd][index4];

		if (double_eq(amb1.value, 0.0)) return false;
		if (double_eq(amb2.value, 0.0)) return false;
		if (double_eq(amb3.value, 0.0)) return false;
		if (double_eq(amb4.value, 0.0)) return false;

		auto ID1 = this->_nOw_ambs[index_key][rec_1st][index1].ID;
		auto ID2 = this->_nOw_ambs[index_key][rec_1st][index2].ID;
		auto ID3 = this->_nOw_ambs[index_key][rec_2nd][index3].ID;
		auto ID4 = this->_nOw_ambs[index_key][rec_2nd][index4].ID;

		if (ID1 != index_key.second)  return false;
		if (ID1 != ID2) return false;
		if (ID1 != ID3) return false;
		if (ID1 != ID4) return false;

		amb.Ow_value[1] = amb1.value;
		amb.Ow_value[2] = amb2.value;
		amb.Ow_value[3] = amb3.value;
		amb.Ow_value[4] = amb4.value;

		amb.Ow_time[1] = make_pair(amb1.beg, amb1.end);
		amb.Ow_time[2] = make_pair(amb2.beg, amb2.end);
		amb.Ow_time[3] = make_pair(amb3.beg, amb3.end);
		amb.Ow_time[4] = make_pair(amb4.beg, amb4.end);

		if (
			index_key.second == t_gamb::AMB_ID::IF_L1_L2 
			)
		{
			amb.value = amb1.value - amb2.value - amb3.value + amb4.value;
		
			amb.sigma = sqrt(amb1.sigma*amb1.sigma + 
				amb2.sigma*amb2.sigma+
				amb3.sigma*amb3.sigma+
				amb4.sigma*amb4.sigma );
		}
		else
		{
			amb.value = 0.0;
			amb.sigma = 0.0;
			amb.DD_LX_value[index_key.second] = 0.0;
			amb.DD_LX_sigma[index_key.second] = 0.0;
			amb.DD_LX_fixed[index_key.second] = false;
			return false;
		}


		return true;
	}

	void t_gallambs::_get_MW(const string& rec, const string& sat, const int& freq, const set<t_gtime>& time_table, set<t_gtime>::iterator beg, set<t_gtime>::iterator end,
		double& value, double& sigma)
	{
		vector<pair<double, double> > mwobs;
		auto sys = t_gsys::sat2gsys(sat);
		auto& mwtable = this->_MW_used[sys][rec];
		auto mw_begin_iter = mwtable.lower_bound(make_pair(sat, *beg));
		auto mw_end_iter = mwtable.upper_bound(make_pair(sat, *end));
		while (mw_begin_iter!=mw_end_iter) {


			double mw_obs = mw_begin_iter->second[freq];

			
			if (double_eq(mw_obs, 0.0))
			{
				mw_begin_iter++;
				continue;
			}
	

			mwobs.push_back(make_pair(mw_obs, 1.0));
			mw_begin_iter++;
		}

		if (mwobs.empty())
		{
			value = 0.0;
			sigma = 0.0;
			return;
		};

		double mean_mw = 0.0, _sigma = 0.0, mean_sigma = 0.0;
		getMeanWgtVec(mwobs, mean_mw, _sigma, mean_sigma);

	
		if (double_eq(mean_sigma, 0.0)) mean_sigma = 999.0;

		value = mean_mw;
		sigma = mean_sigma;

	}

	set<pair<GSYS, t_gamb::AMB_ID>> t_gallambs::get_Dds_index_key(const t_gamb::AMB_ID& ID)
	{
		set<pair<GSYS, t_gamb::AMB_ID>> tmp;

		for (const auto& index_key : _nDd_ambs)
		{
			if (index_key.first.second == ID)
			{
				tmp.insert(index_key.first);
			}
		}
		return tmp;
	}

	set<t_gamb::AMB_ID> t_gallambs::IDs()
	{
		return _IDs;
	}


	bool sort_Dd(const t_gamb_Dd& A, const t_gamb_Dd& B)
	{
		int tmp_A = floor((A.end - A.beg) / 300.0);
		int tmp_B = floor((B.end - B.beg) / 300.0);

		if (tmp_A > tmp_B) return true;
		if (tmp_A < tmp_B) return false;

		if (A.sat_1st < B.sat_1st) return true;
		if (A.sat_1st > B.sat_1st) return false;

		if (A.sat_2nd < B.sat_2nd) return true;
		if (A.sat_2nd > B.sat_2nd) return false;

		if (A.rec_1st < B.rec_1st) return true;
		if (A.rec_1st > B.rec_1st) return false;

		if (A.rec_2nd < B.rec_2nd) return true;
		return false;
	}
	
	void MW_OF_AMB::add_data(const double& new_num)
	{
		if (this->arr.size() < MW_OF_AMB::max_arr_size) {
			this->arr.push_back(new_num);
		}
		else {
			this->arr[this->_next_arr_loc] = new_num;
			this->_next_arr_loc = (this->_next_arr_loc + 1) % this->max_arr_size;
		}
		this->recompute_flag = true;

	}

	void MW_OF_AMB::clean_mean_sigma(double& mean, double& sigma)
	{
		if (this->recompute_flag) {
			this->recompute_flag = false;
			this->_recompute_mean_sigma();
		}
		mean = this->_mean;
		sigma = sqrt(this->_variance);
		return;
		
	}

	void MW_OF_AMB::_recompute_mean_sigma()
	{
		vector<pair<double, double> > arr_priori;
		for (auto& num : arr) {
			arr_priori.push_back(make_pair(num, 1.0));
		}
		double mean_mw = 0.0, _sigma = 0.0, mean_sigma = 0.0;
		int erase_num = 0;
		getMeanWgtVec(arr_priori, mean_mw, _sigma, mean_sigma, &erase_num);

		if (double_eq(mean_sigma, 0.0)) mean_sigma = 999.0;

		this->_mean = mean_mw;
		this->_variance = mean_sigma * mean_sigma;
		this->length = arr_priori.size() - erase_num;
	}

	bool _check_amb_depend(double pdE[], double pdC[], int iNamb, int iN_oneway, vector<int> arriIpt2ow, int& iNdef) {
		double dC_dot = 0.0;
		const double dOper[4] = { 1.0, -1.0, -1.0, 1.0 };
		double dEPS = 1e-12;

		// outrange 
		for (int j = 0; j < iN_oneway; ++j) {
			if (arriIpt2ow[j] > iNamb) return true;
		}

		// schmit-gram
		for (int i = 0; i < iNdef; ++i) {
			pdC[i] = 0.0;
			for (int j = 0; j < iN_oneway; ++j) {
				pdC[i] += pdE[i * iNamb + arriIpt2ow[j] - 1] * dOper[j];
			}
			dC_dot = dC_dot + pdC[i] * pdC[i];
		}
		bool isLdepend = false;
		// change 1e-12 to 1e-8, jqwu
		if (abs(dC_dot - iN_oneway) <= 1e-8) return true;

		// add iNdef vector
		for (int j = 0; j < iNamb; ++j) pdE[iNdef * iNamb + j] = 0.0;
		for (int j = 0; j < iN_oneway; ++j) pdE[iNdef * iNamb + arriIpt2ow[j] - 1] = dOper[j];

		if (dC_dot > dEPS) {
			for (int i = 0; i < iNdef; ++i) {
				if (abs(pdC[i]) < dEPS) continue;
				for (int j = 0; j < iNamb; ++j)
					pdE[iNdef * iNamb + j] -= pdC[i] * pdE[i * iNamb + j];
			}
			//Orthogonalization
			dC_dot = 0.0;
			for (int j = 0; j < iNamb; ++j)
				dC_dot += pdE[iNdef * iNamb + j] * pdE[iNdef * iNamb + j];
			dC_dot = sqrt(dC_dot);
			for (int j = 0; j < iNamb; ++j)
				pdE[iNdef * iNamb + j] /= dC_dot;
		}
		else {
			// Here means the new vector is completely independent from previous vector
			dC_dot = sqrt(iN_oneway * 1.0);
			for (int j = 0; j < iN_oneway; ++j)
				pdE[iNdef * iNamb + arriIpt2ow[j] - 1] /= dC_dot;
		}
		pdC[iNdef] = dC_dot;
		iNdef = iNdef + 1;
		return false;
	}
}


