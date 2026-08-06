/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-12-23  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gupdatepar.cpp
* @brief		define update parameter IF model
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-12-23
*
*/

#include "gproc/gupdateparIF.h"
#include "gutils/gturboedit.h"

namespace great
{
	t_gupdateparIF::t_gupdateparIF(shared_ptr<t_gcycleslip> cycleslip)
	{
		_cycleslip = cycleslip;
	}


	t_gupdateparIF::~t_gupdateparIF()
	{
	}

	void t_gupdateparIF::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info)
	{
		
		if (_lite_update_amb)
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, t_gpar::AMB_IF, update_info);
		}
		else
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, t_gpar::AMB_IF, _cycleslip.get(), update_info);
		}

	}

	void t_gupdateparIF::_update_one_type_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata,t_gpar::t_type ambtype, t_gcycleslip* _slip, t_gupdateparinfo& update_info)
	{

		map<string, set<string> > mapPrn;
		for (auto obs_data : obsdata)
		{
			string rec = obs_data.site();
			string sat = obs_data.sat();
			mapPrn[rec].insert(sat);


			int idx_amb = allpars.getParam(rec, ambtype, sat);

			// judge cycle  and else
			t_gtime crt_epoch = epoch;
			// first add new amb
			// second update old amb if cycleslip
			if (idx_amb < 0 && _slip->num_of_amb_arc(rec, sat, crt_epoch) ||
				idx_amb >= 0 && _slip->cycle_slip(obs_data, crt_epoch))
			{
				int idx_max = allpars[allpars.parNumber() - 1].index;
				t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1, sat);

				int amb_flag = _slip->num_of_amb_arc(rec, sat, epoch);
				_slip->set_amb_flag(rec, sat, amb_flag);
				auto tb_slip = dynamic_cast<t_gturboedit*>(_slip);

				t_gtime end_time = LAST_TIME;
				if (tb_slip) {
					end_time = tb_slip->get_crt_amb_end(rec, sat);
				}
				par_newamb.value(0.0);
				par_newamb.setTime(epoch, end_time);
				par_newamb.apriori(10000.0);

				// add new amb
				if (idx_amb < 0)
				{
					update_info.add(par_newamb);
				
				}
				// reset the new amb
				else
				{
					allpars[idx_amb].end = epoch - _intv;
					update_info.add(idx_amb + 1, par_newamb);
				}
			}
		}

		// delete old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == ambtype)
			{
				
				if (allpars[i].end < epoch)
				{
					if (!update_info.exist(i + 1))
					{ 
						update_info.add(i + 1);
					}
				}
				
				else if (mapPrn[allpars[i].site].count(allpars[i].prn) == 0)
				{
					if (!update_info.exist(i + 1)) 
					{
						allpars[i].end = epoch - _intv; 
						update_info.add(i + 1);
					}
				}
				else
				{
				
					continue;
				}
			}
		}

	}

	void t_gupdateparIF::_update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gpar::t_type ambtype, t_gupdateparinfo& update_info)
	{
		map<t_gpar::t_type, gnut::FREQ_SEQ> table_ambtype_freq=
		{
			{t_gpar::t_type::AMB_IF,gnut::FREQ_2},
			{t_gpar::t_type::AMB13_IF,gnut::FREQ_3},
			{t_gpar::t_type::AMB14_IF,gnut::FREQ_4},
			{t_gpar::t_type::AMB15_IF,gnut::FREQ_5}
		};

		for (auto obs_data : obsdata)
		{
			string rec = obs_data.site();
			string sat = obs_data.sat();
			GSYS   sys = obs_data.gsys();

			int idx_amb = allpars.getParam(rec, ambtype, sat);
			GOBS obsL1, obsL2;
			GOBSBAND b1, b2;
			b1 = _band_index[sys][FREQ_1];
			if (table_ambtype_freq.count(ambtype)){
				b2 = _band_index[sys][table_ambtype_freq[ambtype]];
			}
			else{
				b2 = BAND;
			}

			obsL1 = obs_data.select_phase(b1);
			obsL2 = obs_data.select_phase(b2);
			
			if (idx_amb < 0)
			{
				// have obs but have no amb pars: add new par
				if (obsL1 != GOBS::X && obsL2 != GOBS::X)
				{
					int idx_max = allpars[allpars.parNumber() - 1].index;
					t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1, sat, false);

					par_newamb.value(0.0);
					par_newamb.setTime(epoch, epoch);
					par_newamb.apriori(9000.0);
					update_info.add(par_newamb);
				}
			}
			else
			{
				if (obsL1 != GOBS::X && obsL2 != GOBS::X)
				{
					// have obs, have amb pars, no cycle slip: update the end time of amb pars
					if (obs_data.getlli(obsL1) < 1 && obs_data.getlli(obsL2) < 1)
					{
						allpars[idx_amb].end = epoch;
					}
					// have obs, have amb pars, have cycle slip: update amb pars to the new one 
					else
					{
						allpars[idx_amb].end = epoch - _get_rec_interval(rec);

						int idx_max = allpars[allpars.parNumber() - 1].index;
						t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1, sat, false);

						par_newamb.value(0.0);
						par_newamb.setTime(epoch, epoch);
						par_newamb.apriori(9000.0);

						update_info.add(idx_amb + 1, par_newamb);
					}
				}
				// no this observation type, have amb pars: delete this amb par
				else  // obsL1==X || obsL2==X
				{
					
					allpars[idx_amb].end = epoch - _get_rec_interval(rec);
					update_info.add(idx_amb + 1);
				}
			}
		}

		// delete old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == ambtype)
			{
				string rec = allpars[i].site;
				// no obs, have amb pars: delete this amb par
				if (allpars[i].end <= epoch - _get_rec_interval(rec))
				{
					if (!update_info.exist(i + 1))
					{
						update_info.add(i + 1);
					}
				}
				else if (allpars[i].amb_ini == true && allpars[i].beg + -1 * _get_rec_interval(rec) < epoch)
				{
					t_gpar newamb = allpars[i];
					newamb.beg = epoch;
					if (!update_info.exist(i + 1))
					{
						update_info.add(i + 1, newamb);
					}
				}
			}
		}

	}
}