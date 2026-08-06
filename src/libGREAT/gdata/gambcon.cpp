/**
*
* @verbatim
History
	-1.0 Jie Li  2019-11-14  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambcon.cpp
* @brief	Storage the ambigulity constraint files' data
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-11-14
*
*/

#include "gdata/gambcon.h"
 
using namespace std;
namespace great {

	t_gambcon::t_gambcon()
	{
		gtrace("t_gambcon::constructor");
		id_type(t_gdata::AMBCON);
	}

	t_gambcon::~t_gambcon()
	{
	}

	void t_gambcon::set_header(t_gtime beg, double intv, int num_amb, set<string>& sats, set<string>& recs)
	{
		_beg = beg;
		_interval = intv;
		_num_fixamb = num_amb;
		_sats = sats;
		_sites = recs;
	}


	void t_gambcon::add_ddamb(t_gamb_Dd& ddamb)
	{
		_ddamb.push_back(ddamb);
		_nddamb[ddamb.beg].push_back(ddamb);
	}

	vector<t_gamb_Dd> t_gambcon::get_Dds()
	{
		return this->_ddamb;
	}


	vector<t_gamb_Dd>& t_gambcon::get_Dds(const t_gtime& time)
	{
		
		return _nddamb[time];
	}
}//namespace

