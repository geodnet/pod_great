/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambfixDd.h
* @brief	for the use of DD fix ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*
*/

#ifndef GAMBFIXDD_H
#define GAMBFIXDD_H 

#include "ExportHeaders/ExportAMB.h"
#include "gambfix/gambfix.h"
#include "gset/gsetamb.h"
#include "gdata/gambcon.h"

using namespace std;

namespace great
{

	/**
	*@ brief class for LAMADA search finding fixing ambiguity.
	*/
	class AMB_LIBRARY_EXPORT t_gambfixDd : public t_gambfix
	{
	public:
		t_gambfixDd();
		t_gambfixDd(t_glog* log, t_gsetbase* set);
		
		virtual ~t_gambfixDd();

		
	
		
		bool ProcessBatch(t_gtime beg, t_gtime end);

		void update_data(t_gallpar* pars, t_gallobs* obs,t_gallbias* bias,t_gallobj* obj,t_gtime beg_time, t_gtime end_time);


	protected:
		map<t_gamb::AMB_ID, t_gamb_baseline_recs>  _gRecBl;
		map<t_gamb::AMB_ID, t_gamb_baseline_sats>  _gSatBl;

		
		bool _screen_rec_and_sat(const set<t_gamb::AMB_ID>& IDs, t_gallpar* pars, t_gsetbase* set, t_gallobj* obj);
		void _write_ambcon(bool write_file = true);

		string _class_id = "t_gambfixDd";

		bool _get_rec_and_sat(const t_gamb::AMB_ID& ID, t_gallpar* pars, set<string>& rec_list, set<string>& sat_list);
	};
}

#endif