/**
*
* @verbatim
	History
	 -1.0 GREAT	    2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			  gturboedit.cpp
* @brief		  CRS coordinate transform to ACR coordinate
*
* @author         GREAT, Wuhan University
* @version		  1.0.0
* @date			  2019-01-04
*
*/

#include "gutils/gturboedit.h"
#include "gset/gsetrec.h"
#include "gset/gsetinp.h"
#include "gutils/gfileconv.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <set>
#include "gio/gfile.h"

#ifdef WIN32
	#include <io.h>
#endif // WIN32

using namespace std;
using namespace gnut;

namespace great
{
	great::t_gturboedit::t_gturboedit()
	{
	}

	//great::t_gturboedit::t_gturboedit(t_gsetbase *set, t_glog *log, bool isTriple) : t_gcycleslip(set, log)
	great::t_gturboedit::t_gturboedit(t_gsetbase * gset, t_glog * glog, int index) : t_gcycleslip(gset, glog)
	{
		_gset = gset;
		_glog = glog;

		
		_index = index;
		// read log information from log file
		t_gtime beg_time = dynamic_cast<t_gsetgen*>(gset)->beg();

		//recs to rec_all for processing both LEO and SITE in combined POD by zhangwei
		set<string> rec_list = dynamic_cast<t_gsetgen*>(gset)->rec_all(); 

		// add lite mode
		_apply_carrier_range = dynamic_cast<t_gsetproc*>(gset)->apply_carrier_range();

	}

	great::t_gturboedit::~t_gturboedit()
	{
	}

	bool great::t_gturboedit::use_of_obs(const string& site, const string& prn, const t_gtime& time) {
		 return false;
		
	}

	

	int great::t_gturboedit::num_of_amb_arc(const string& site, const string& prn, const t_gtime& time)
	{
		 return -1;
		
	}

	bool great::t_gturboedit::cycle_slip(t_gsatdata& obsdata, const t_gtime& time)
	{

		// change some codes
		string sat = obsdata.sat();
		string rec = obsdata.site();

		if (_amb_info_file_exist[rec])
		{
			if (_amb_flag.find(rec) != _amb_flag.end() && _amb_flag[rec].find(sat) != _amb_flag[rec].end() && _amb_flag[rec][sat] == 0)
			{
				// todo
		
				throw runtime_error("Can't judge cycle slip unless already in one amb arc");
			}
			int flag_now = num_of_amb_arc(rec, sat, time);
			if (flag_now && flag_now != _amb_flag[rec][sat])
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		return false;
	}

	int great::t_gturboedit::get_active_amb(string site)
	{
		return _active_amb[site];
	}

	bool great::t_gturboedit::is_carrier_range(const string& site, const string& prn, const t_gtime& time)
	{
		 return false;
		
	}

	bool great::t_gturboedit::is_carrier_range(const string& site, const string& prn, const t_gtime& time, double& c1, double& c2)
	{
		 return false;
	
	}

	void great::t_gturboedit::set_active_amb(string site, int active_num)
	{
		_amb_info_file_exist[site] = true;
		_active_amb[site] = active_num;
	}

	void great::t_gturboedit::set_amb_flag(const string& rec, const string& sat, int flag)
	{
		_amb_flag[rec][sat] = flag;
		_new_amb[rec][sat] = true;
	}

	int great::t_gturboedit::get_amb_flag(const string& rec, const string& sat)
	{
		return _amb_flag[rec][sat];
	}

	bool great::t_gturboedit::new_amb(const string& rec, const string& sat)
	{
		if (_new_amb.find(rec) == _new_amb.end())
		{
			return false;
		}
		if (_new_amb[rec].find(sat) == _new_amb[rec].end())
		{
			return false;
		}
		return _new_amb[rec][sat];
	}

	void great::t_gturboedit::set_new_amb(const string& rec, const string& sat, bool isNew)
	{
		_new_amb[rec][sat] = isNew;
	}

	t_gtime great::t_gturboedit::get_crt_amb_end(const string& site, const string& sat)
	{
		if (_amb_flag[site][sat] <= 0)
		{
			string tmp = "sat : " + sat + " rec : " + site + int2str(_index);
			throw runtime_error("Crt amb isn't exist!!! " + tmp);
		}

		
	}



	void great::t_gturboedit::add_ambflag(string site, string sat, string description, t_gtime beg, t_gtime end)
	{
		if (description == "AMB")
		{
			_cycle_flag[site][sat].push_back(make_pair(beg, end));
		}
		else
		{
			_cycle_flag_unused[site][sat].push_back(make_pair(beg, end));
		}

	}
}

