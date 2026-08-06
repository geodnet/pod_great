/**
*
* @verbatim
	History
	 -1.0 Guolong Feng  2020-03-20 created
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gexeturboedit.cpp
* @brief		main control class for turboedit
* @author       Guolong Feng, Wuhan University
* @date		    2020-03-20
*
*/

#include "gproc/gexeturboedit.h"

namespace great
{

	t_gexeturboedit::t_gexeturboedit(t_gsetbase* set, t_glog* glog, t_gallobs* gobs, t_gallnav* gnav):
		_set(set),
		_log(glog),
		_allobs(gobs),
		_allnav(gnav),
		_allobj(nullptr),
		_check_sf(false),
		_check_pc(true),
		_amb_output(true),
		_simulation(false),
		_check_short(true),
		_site_kinematic(false),
		_check_statistics(true),
		_outliers_proc(set, glog)
	{
		// common intialization
		_initCommonProcess();

		
		_initRealTimeProcess();
		// output control infomation
		_outputControlInfo();
	
	}

	t_gexeturboedit::t_gexeturboedit(t_gsetbase* set, t_glog* glog, t_gallproc* gdata) :
		_set(set),
		_log(glog),
		_check_sf(false),
		_check_pc(true),
		_amb_output(true),
		_simulation(false),
		_check_short(true),
		_site_kinematic(false),
		_check_statistics(true),
		_outliers_proc(set, glog)
	{
		// get data 
		_allobs = dynamic_cast<t_gallobs*>((*gdata)[t_gdata::ALLOBS]);
		_allnav = dynamic_cast<t_gallnav*>((*gdata)[t_gdata::GRP_EPHEM]);
		_allobj = dynamic_cast<t_gallobj*>((*gdata)[t_gdata::ALLOBJ]);

		// common intialization
		_initCommonProcess();

		// init post-processing / real-time-processing
		_initRealTimeProcess();
	
		// output control infomation
		_outputControlInfo();

	}

	t_gexeturboedit::~t_gexeturboedit()
	{
		


		if (_debug_turbo) 
		{ 
			if (_debug_turbo->is_open()) _debug_turbo->close();
			delete _debug_turbo; 
			_debug_turbo = nullptr;
		}

	}

	bool t_gexeturboedit::ProcessBatch(string site, const t_gtime& begT, const t_gtime& endT, double sampling, bool isMarkEnd)
	{
		_crt_rec = site;
		_interval = sampling;

		// Check Time Duration
		if (begT > endT)//lvhb begT>endT
		{
			_log->comment(1, _crt_rec, "Error in time settings : beg time >= end time !");
			_log->comment(1, _crt_rec, "The beg time is : " + begT.str_ymdhms());
			_log->comment(1, _crt_rec, "The end time is : " + endT.str_ymdhms());
			return false;
		}

			if(!_realtimeProcess(begT, endT,isMarkEnd)) return false;
		
		return true;
	}

	

	void t_gexeturboedit::set_site_crd(string site, t_gtriple crd)
	{
		map<string, t_gtriple>::iterator it;
		it = _rec_crds.find(site);
		if (it == _rec_crds.end())
		{
			_rec_crds[site] = crd;
		}
	}

	

	void  t_gexeturboedit::_outputControlInfo()
	{
		stringstream str2;
		if (_log) {
			_log->comment(0, "NOTE", "***********************************************");
			_log->comment(0, "NOTE", "+++  Control Parameters from XML file  +++");

			// Realtime
			_log->comment(0, "NOTE", "RealtimeMode : false " );
			_log->comment(0, "NOTE", "LiteMode     :  true" );
		}
		//// MW
		if (_log)_log->comment(0, "NOTE", "TurboEdit(MW):  " + bl2str(_check_mw));
		if (_check_mw)
		{
			str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _mw_limit;
			if (_log)_log->comment(0, "NOTE", "MWLimitation :  " + str2.str() + " cycles");
		}
		// GF
		if (_log)_log->comment(0, "NOTE", "TurboEdit(GF):  " + bl2str(_check_gf));
		if (_check_gf)
		{
			
			str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _gf_limit;
			if (_log)_log->comment(0, "NOTE", "GFLimitation :  " + str2.str() + " meters");
			
		}

		// check gap
		if (_log)_log->comment(0, "NOTE", "Check Gap    :  " + bl2str(_check_gap));
		if (_check_gap)
		{
			str2.str(""); str2 << fixed << left << setw(7) << _gap_limit;
			if (_log)_log->comment(0, "NOTE", "GapLimit     :  " + str2.str() + " epochs ");
		}

	

	}

	void t_gexeturboedit::_initCommonProcess()
	{
		// Global Settings
		_beg_time = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(_set)->end();
		_sys_list = dynamic_cast<t_gsetgen*>(_set)->sys();        // processing systems
		_rec_list = dynamic_cast<t_gsetgen*>(_set)->recs();       // processing receivers
		for (auto it : _sys_list)                                 // processing satellites
		{
			GSYS gsys = t_gsys::str2gsys(it);
			std::set<string> sat_temp = dynamic_cast<t_gsetgnss*>(_set)->sat(gsys);
			_sat_list.insert(sat_temp.begin(), sat_temp.end());  // insert satellites
		}
		// exclude satellites mainly For simulation
		// real data can control processed sats by module <gps>/<bds>/<gal>/..., also useful for real data
		std::set<string> sat_rm = dynamic_cast<t_gsetgen*>(_set)->sat_rm();
		for (auto iter : sat_rm) { _sat_list.erase(iter); }

		// TurboEdit Settings
		_min_elev = dynamic_cast<t_gsetproc*>(_set)->minimum_elev();
		_check_mw = dynamic_cast<t_gsetturboedit*>(_set)->checkMW(_mw_limit);                // cycle for post/real_time
		_check_gf = dynamic_cast<t_gsetturboedit*>(_set)->checkGF(_gf_limit, _gf_rms_limit); // meter for real_time/cycle for post
		_check_gap = dynamic_cast<t_gsetturboedit*>(_set)->checkGap(_gap_limit);
		_ephemeris = dynamic_cast<t_gsetturboedit*>(_set)->isEphemeris();    
		_site_kinematic = dynamic_cast<t_gsetproc*>(_set)->pos_kin();
		if (!_allnav) _ephemeris = false;

		if (!_check_gap)  _gap_limit = 999999;

	}

	void t_gexeturboedit::_initRealTimeProcess()
	{
		
		_smooth_window = dynamic_cast<t_gsetturboedit*>(_set)->smoothWindows();

		double sample = dynamic_cast<t_gsetgen*>(_set)->sampling();
		if (sample <= 1)
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.05;                     ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 3.0;                      ///< cycles
		}
		else if (sample <= 20)
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.1 / 20.0 * sample + 0.05;   ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 2.5 / 20.0 * sample + 2.5;    ///< cycles
		}
		else if (sample <= 60)
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.15;                    ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 4.0;                     ///< cycles
		}
		else if (sample <= 100)
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.25;                    ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 5.0;                     ///< cycles
		}
		else
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.3;                   ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 6.0;                   ///< cycles
		}

		_frequency = dynamic_cast<t_gsetproc*>(_set)->frequency();

		if (_ephemeris)
		{
			_glonass_freq_id = _allnav->glo_freq_num();
		}


		// Get site coordinates from XML/SNX files
		for (const auto& it : _rec_list)
		{
			// XML 
			_rec_crds[it] = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(it);
			// SNX
			if (_rec_crds[it].zero() && _allobj)
			{
				shared_ptr<t_gobj> _grec = _allobj->obj(it);
				t_gtriple crd(0.0, 0.0, 0.0);
				t_gtriple std(0.0, 0.0, 0.0);
				if (_grec)
				{
					_grec->get_recent_crd(_beg_time, 0.1, crd, std);
					_rec_crds[it] = crd;
				}
			}
		}

		if (1) // lvhb shut down the outputting file in 202011
		{
			if (_debug_turbo)
			{
				if (_debug_turbo->is_open()) _debug_turbo->close();
				delete _debug_turbo;
				_debug_turbo = nullptr;
			}
		}
		else {
			_debug_turbo = new t_giof;
			_debug_turbo->tsys(t_gtime::GPS);
			_debug_turbo->mask("debug_turbo.log");
			_debug_turbo->append(dynamic_cast<t_gsetout*>(_set)->append());
		}
	}

	bool t_gexeturboedit::_realtimeProcess(const t_gtime& begT, const t_gtime& endT,bool isMarkEnd)
	{
		_crt_time = begT;
		vector<shared_ptr<t_gobsgnss>>       epoData;
		map<string, map<GOBS, set<GOBS> > > markEndTable;
		while (_crt_time <= endT)
		{
			
			if (_ephemeris && (_site_kinematic || _rec_crds[_crt_rec].zero()))
			{
			}

			if(_inputEpoData.size() == 0) epoData = _allobs->obs_pt(_crt_rec, _crt_time);
			else                          epoData = _inputEpoData;

			if (epoData.size() == 0)
			{
				_crt_time = _crt_time + _interval;
				continue;
			}

			
			for (auto &oneObs : epoData)
			{
				string sat = oneObs->sat();
				GSYS gsys = oneObs->gsys();
				string site = oneObs->site();

				if (site != _crt_rec) continue;  

				// set glonass Frequency ID
				if (sat.substr(0, 1) == "R" && oneObs->channel() >= DEF_CHANNEL)
				{
					oneObs->channel(_glonass_freq_id[sat]);
				}

				if (_crt_time > _beg_time)
				{
					auto satdata_pre = _allobs->obs_prn_pt(_crt_rec, sat, _crt_time - _interval * _gap_limit, _crt_time - _interval);
					if (!satdata_pre.empty()) {
						auto epodata_pre = satdata_pre.rbegin();
						
						_outliers_proc.flagRangeOutliers(*epodata_pre, oneObs, _interval);

					}
				}
				if (_frequency >= 1)
				{
					set<GOBSBAND> bands = oneObs->band_avail(true);
					vector<GOBSBAND> sorted_bands = sort_band(gsys, bands);

					for (unsigned int i = 1; i < sorted_bands.size(); ++i)
					{
						if (sorted_bands[0] == sorted_bands[i]) continue;

						
						_realTimeMWGF_new(oneObs, sorted_bands[0], sorted_bands[i]);
						
					}
				}
				if (isMarkEnd) {
					auto& crtRtData = _rt_data[_crt_rec][sat];
					for (auto& obs1_iter : crtRtData) {
						for (auto& obs2_iter : obs1_iter.second.last_time) {
							if (obs2_iter.second == _crt_time &&
								oneObs->getlli(obs2_iter.first) >= 1) {
								markEndTable[sat][obs1_iter.first].insert(obs2_iter.first);
							}
						}
					}
				}
			}
			// Next
			_crt_time = _crt_time + _interval;

		}
		if (isMarkEnd) {
			epoData = _allobs->obs_pt(_crt_rec, endT);
			for (auto& oneObs : epoData) {
				auto sat = oneObs->sat();
				auto gsys = oneObs->gsys();
				set<GOBSBAND> bands = oneObs->band_avail(true);
				vector<GOBSBAND> sorted_bands = sort_band(gsys, bands);
				for (unsigned int i = 1; i < sorted_bands.size(); ++i)
				{
					if (sorted_bands[0] == sorted_bands[i]) continue;

					GOBS obs_L1 = oneObs->select_phase(sorted_bands[0], true);
					GOBS obs_L2 = oneObs->select_phase(sorted_bands[i], true);
					if (markEndTable[sat].count(obs_L1) &&
						markEndTable[sat][obs_L1].count(obs_L2)) {
						oneObs->addlli(obs_L1, 1);
						oneObs->addlli(obs_L2, 1);
					}
				}
			}
		}
		return true;
	}

	void t_gexeturboedit::_realTimeMWGF_new(shared_ptr<t_gobsgnss> oneObs, GOBSBAND b1, GOBSBAND b2)
	{
		string sat = oneObs->sat();
		//// get the set<GOBS>
		//set<GOBS> obs_vec1 = oneObs->obs_phase(b1);
		//set<GOBS> obs_vec2 = oneObs->obs_phase(b2);
		// add select Phase
		GOBS obs_L1 = oneObs->select_phase(b1, true);
		GOBS obs_L2 = oneObs->select_phase(b2, true);
		set<GOBS> obs_vec1, obs_vec2;
		if(obs_L1 != GOBS::X) obs_vec1.insert(obs_L1);
		if(obs_L2 != GOBS::X) obs_vec2.insert(obs_L2);

		t_gobs gobs_P1 = t_gobs(oneObs->select_range(b1, true));
		t_gobs gobs_P2 = t_gobs(oneObs->select_range(b2, true));

		double distance, sat_clk;
		double sat_elev = _min_elev;
		double scale_factor = 1;
		double thres_MW = 0, thres_GF = 0;
		if (this->_getSatInfo(sat, distance, sat_elev, sat_clk))
		{
			if (sat_elev < _min_elev) sat_elev = _min_elev;
		}

		ostringstream os;
		os.str("");

		for (auto obs1 : obs_vec1)
		{
			t_gobs gobs_L1 = t_gobs(obs1);

			for (auto obs2 : obs_vec2)
			{
				t_gobs gobs_L2 = t_gobs(obs2);

				double crt_MW_val = oneObs->MW_cycle(gobs_L1, gobs_L2, gobs_P1, gobs_P2);     // MW combination [cycles]
				double crt_GF_val = oneObs->L4(gobs_L1, gobs_L2);                             // Geometry-free combination [meters]

				os << fixed << setw(10) << " MW[c] " << oneObs->sat() << "  " << setw(6) << _crt_time.sod()
					<< "  " << gobs2str(obs1) << "  " << gobs2str(obs2)
					<< "  " << gobs2str(gobs_P1.gobs()) << "  " << gobs2str(gobs_P2.gobs())
					<< setw(15) << setprecision(4) << crt_MW_val << endl;

				os << fixed << setw(10) << " GF[m] " << oneObs->sat() << "  " << setw(6) << _crt_time.sod()
					<< "  " << gobs2str(obs1) << "  " << gobs2str(obs2)
					<< setw(15) << setprecision(4) << crt_GF_val << endl;

				
				if (_rt_data[_crt_rec][sat][obs1].epo_num.find(obs2) == _rt_data[_crt_rec][sat][obs1].epo_num.end()
					|| _rt_data[_crt_rec][sat][obs1].mwslip[obs2] || _rt_data[_crt_rec][sat][obs1].gfslip[obs2]
					|| oneObs->getlli(obs_L1) >= 1 || oneObs->getlli(obs_L2) >= 1)
				{
                    bool is_start = false;
                    if (_rt_data[_crt_rec][sat][obs1].epo_num.find(obs2) == _rt_data[_crt_rec][sat][obs1].epo_num.end()){
                        is_start = true;
                    }
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] = 1;
					_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;
					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;
					_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2] = make_pair(crt_MW_val, 0.0);

					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].clear();
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(0.0);
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
                        if (!is_start){
                            oneObs->addlli(obs1, 1);
                            oneObs->addlli(obs2, 1);
                        }
					}
					continue;
				}

				
				double dt = _crt_time.diff(_rt_data[_crt_rec][sat][obs1].last_time[obs2]);
				double dif_MW = crt_MW_val - _rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].back();
				bool badRange = oneObs->getoutliers(gobs_P1.gobs()) >= 1 || oneObs->getoutliers(gobs_P2.gobs()) >= 1;
				
				if (dt > _interval && abs(dif_MW) > 8 && !badRange || _check_gap && dt > _gap_limit * _interval)
				{
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] = 1;
					_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;
					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;
					_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2] = make_pair(crt_MW_val, 0.0);

					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].clear();
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(0.0);
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);


					ostringstream msg;
					msg << fixed << setw(10) << " Data Gap or Max_diff_MW : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
						<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << dif_MW << setw(16) << setprecision(3) << dt
						<< setw(16) << setprecision(3) << _gap_limit * _interval;
					_log->comment(1, _crt_rec, msg.str());

					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
					continue;
				}

				int gapNum = dt / _interval - 1;

				if (_check_gf)
				{
					double pre_GF_val = _rt_data[_crt_rec][sat][obs1].GF[obs2];

					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;

					if (sat_elev < 15)
					{
						thres_GF = -_gf_limit * 1.0 / 15.0 * sat_elev + 2 * _gf_limit;
					}
					else
					{
						thres_GF = _gf_limit;
					}

					if (dt >= 90 && gapNum >= 3)
					{
						thres_GF = min(thres_GF * gapNum, 0.35);
					}

					if (abs(crt_GF_val - pre_GF_val) > thres_GF)
					{
						_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = true;
						ostringstream msg;
						msg << fixed << setw(10) << " GF cycle slip : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
							<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << pre_GF_val << setw(16) << setprecision(3) << crt_GF_val
							<< setw(16) << setprecision(3) << (crt_GF_val - pre_GF_val) << setw(16) << setprecision(3) << thres_GF;
						_log->comment(1, _crt_rec, msg.str());
					}

				}

				if (_check_mw && !badRange)
				{
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] += 1;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;

					int epo_num = _rt_data[_crt_rec][sat][obs1].epo_num[obs2];
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					double pre_MW_val = _rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first;
					double pre_MW_std = _rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second;
					
					if (_smooth_window == 0 || epo_num <= _smooth_window)
					{
						double smooth_MW_val = pre_MW_val * (epo_num - 1) / epo_num + crt_MW_val / epo_num;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first = smooth_MW_val;
						double smooth_MW_std = pre_MW_std * (epo_num - 1) / epo_num + pow(crt_MW_val - pre_MW_val, 2) / epo_num;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second = smooth_MW_std;
					}
					else
					{
						double smooth_MW_val = 0.0;
						for (int i = epo_num - _smooth_window + 1; i <= epo_num; i++)
						{
							smooth_MW_val += _rt_data[_crt_rec][sat][obs1].origin_MWs[obs2][i] / _smooth_window;
						}
						double smooth_MW_std = 0.0;
						for (int i = epo_num - _smooth_window + 1; i <= epo_num; i++)
						{
							smooth_MW_std += pow(_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2][i] - smooth_MW_val, 2) / _smooth_window;
						}
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first = smooth_MW_val;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second = smooth_MW_std;

						auto& MW_vec = _rt_data[_crt_rec][sat][obs1].origin_MWs[obs2];
						MW_vec.erase(MW_vec.begin());
						_rt_data[_crt_rec][sat][obs1].epo_num[obs2] -= 1;
					}

					if (sat_elev < 20) thres_MW = -_mw_limit * 0.1 * sat_elev + 3 * _mw_limit;
					else               thres_MW = _mw_limit;

					if (gapNum <= 2)		scale_factor = 1.0;
					else if (gapNum <= 4)	scale_factor = 1.25;
					else if (gapNum <= 8)	scale_factor = 1.5;
					else					scale_factor = 2.0;

					if (abs(crt_MW_val - pre_MW_val) > thres_MW * scale_factor && !double_eq(thres_MW, 0.0))
					{
						_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = true;
						if (_log)
						{
							ostringstream msg;
							msg << fixed << setw(10) << " MW cycle slip : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
								<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << pre_MW_val << setw(16) << setprecision(3) << crt_MW_val
								<< setw(16) << setprecision(3) << (crt_MW_val - pre_MW_val) << setw(16) << setprecision(3) << thres_MW;
							_log->comment(1, _crt_rec, msg.str());
						}
					}
				}

				
				if (_rt_data[_crt_rec][sat][obs1].gfslip[obs2] || _rt_data[_crt_rec][sat][obs1].mwslip[obs2])
				{
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] = 1;
					_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;
					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;
					_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2] = make_pair(crt_MW_val, 0.0);

					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].clear();
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(0.0);
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					if(_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
				}

			}
		}

		// Print flt results
		if (_debug_turbo) {
			_debug_turbo->write(os.str().c_str(), os.str().size());
			_debug_turbo->flush();
		}
	}


	bool t_gexeturboedit::_getSatInfo(string sat, double& distance, double& elev, double& sat_clk)
	{
		bool valid = false;
		double sat_crd[3];
		t_gtriple sat_crd_triple, site_sat_vector;

		double delay = 0.0;
		t_gtime sat_time;
		shared_ptr<t_geph> geph;
		while (true)
		{
			sat_time = _crt_time - delay;

			if (_allnav->pos(sat, sat_time, sat_crd, nullptr, nullptr, false) < 0) break;
			sat_crd_triple = t_gtriple(sat_crd);

			if (_rec_crds[_crt_rec].zero()) break;

			site_sat_vector = sat_crd_triple - _rec_crds[_crt_rec];
			distance = site_sat_vector.norm();

			if (abs(distance / CLIGHT - delay) < 1E-9)
			{
				valid = true;
				break;
			}

			delay = distance / CLIGHT;
		}

		// Compute elevation & sat clks
		if (valid)
		{
		
			double var_temp, dclk_temp;
			if (_allnav->clk(sat, sat_time, &sat_clk, &var_temp, &dclk_temp, false) >= 0) {
				sat_clk *= CLIGHT;
				// CRS
				elev = (_rec_crds[_crt_rec][0] * site_sat_vector[0] + _rec_crds[_crt_rec][1] * site_sat_vector[1]
					+ _rec_crds[_crt_rec][2] * site_sat_vector[2]) / _rec_crds[_crt_rec].norm() / distance;
				elev = 90.0 - acos(elev) * 180.0 / G_PI;
				
				return true;
			}
		}

		distance = 0.0; sat_clk = 0.0, elev = 0.0;
		_log->comment(2, _crt_rec, sat + " lack of suitable broadcast ephemeris at Time " + _crt_time.str_ymdhms());
		return false;
	}

}


