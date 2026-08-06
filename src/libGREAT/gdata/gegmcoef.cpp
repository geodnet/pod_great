/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li	2019-01-06	creat the file.
	  -1.1  Jie Li	2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gegmcoef.cpp
 * @brief		The base class used to storage EGM file information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include "gdata/gegmcoef.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gegmcoef::t_gegmcoef()
	{
		gtrace("t_gegmcoef::constructor");
		id_type(t_gdata::EGM);
		_name = "EGM2008";
	}

	/** @brief destructor. */
	t_gegmcoef:: ~t_gegmcoef(){}


	/**
	* @brief add time GM and radius data.
	* @param [in]  t       time of EGM data
	* @param [in]  GM      GM of the Earth
	* @param [in]  radius  radius of the Earth
	* @return	if add data successfully
	*/
	int t_gegmcoef::add_head(t_gtime t, double GM, double radius)
	{
		gtrace("t_gegmcoef::add_head");
		_gmutex.lock();

		_t = t;
		_GM = GM;
		_radius = radius;

		_gmutex.unlock(); return 0;
	}

	/**
	* @brief add EGM order data.
	* @param [in]  order    EGM order
	* @return	if add data successfully
	*/
	int t_gegmcoef::add_order(int order)
	{
		gtrace("t_gegmcoef::add_order");
		_gmutex.lock();

		_order = order;

		_gmutex.unlock(); return 0;
	}

	/**
	* @brief add EGM coefficient data.
	* @param [in]  n    index
	* @param [in]  m    index
	* @param [in]  c    c value
	* @param [in]  s    s value
	* @param [in]  sigma_c    sigma_c value
	* @param [in]  sigma_s    sigma_s value
	* @return	index n
	*/
	int t_gegmcoef::add_egmcoef(int n, int m, double c, double s, double sigma_c, double sigma_s)
	{
		gtrace("t_gegmcoef::add_egmcoef");
		_gmutex.lock();

		t_gravity_coef tmp_coef;
		tmp_coef.c = c;
		tmp_coef.s = s;
		tmp_coef.sigma_c = sigma_c;
		tmp_coef.sigma_s = sigma_s;

		_map_egmcoef[n][m] = tmp_coef;

		_gmutex.unlock(); return n;
	}

	/** @brief whether the EGM coefficien data is empty.
	* @return  bool
	*	@retval   true    EGM coefficien data is empty
	*   @retval   false   EGM coefficien data is existent
	*/
	bool t_gegmcoef::is_empty()
	{
		if (_map_egmcoef.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	* @brief get corresponding c value.
	* @param [in]  n    index
	* @param [in]  m    index
	* @return  c value
	*/
	double t_gegmcoef::get_c(int n, int m)
	{
		gtrace("t_gegmcoef::get_c");
		_gmutex.lock();
		double c = 0.0;
		if (_map_egmcoef.count(n) == 0 || _map_egmcoef[n].count(m) == 0)
		{
			_log->comment(0, "ERROR : gravity degreee is out of range! max gravity degree in EGM file is 150");
			_gmutex.unlock();
			return 0.0;
		}
		c = _map_egmcoef[n][m].c;
		_gmutex.unlock();
		return c;
	}

	/**
	* @brief get corresponding s value.
	* @param [in]  n    index
	* @param [in]  m    index
	* @return  s value
	*/
	double t_gegmcoef::get_s(int n, int m)
	{
		gtrace("t_gegmcoef::get_s");
		_gmutex.lock();

		double s = 0.0;
		if (_map_egmcoef.count(n) == 0 || _map_egmcoef[n].count(m) == 0)
		{
			_log->comment(0, "ERROR : gravity degreee is out of range! max gravity degree in EGM file is 150");
			_gmutex.unlock();
			return 0.0;
		}
		s = _map_egmcoef[n][m].s;
		_gmutex.unlock();
		return s;
	}

	/**
	* @brief get reference time of EGM file.
	* @return  mjd of reference time.
	*/
	double t_gegmcoef::mjd_ref()
	{
		return _t.dmjd();
	}
	/**
	* @brief add a new time coefficient(80*80) group into map data
	* @param [in]  st    start time
	* @param [in]  et    end time
	* @return	index ipt
	*/
	int t_gegmcoef::create_timecoef_group(t_gtime st, t_gtime et)
	{
		gtrace("t_gegmcoef::create_timecoef");
		_gmutex.lock();
		t_time_gravity_coef time_egmcoef;
		int ipt;
		time_egmcoef.st = st;
		time_egmcoef.et = et;
		ipt = _map_time_egmcoef.size() + 1;
		_map_time_egmcoef[ipt] = time_egmcoef;
		_gmutex.unlock();
		return ipt;
	}

	/**
	* @brief add one time coefficient into the corresponding group
	* @param [in]  ipt  group index of time coefficient
	* @param [in]  tag  type of time coefficient
	* @param [in]  n    degree index
	* @param [in]  m    order index
	* @param [in]  c    c value
	* @param [in]  s    s value
	* @param [in]  sigma_c    sigma_c value
	* @param [in]  sigma_s    sigma_s value
	* @return	index n
	*/
	int t_gegmcoef::add_timecoef(int ipt, string tag, int n, int m, double c, double s, double sigma_c, double sigma_s, double p)
	{
		gtrace("t_gegmcoef::add_timecoef");
		_gmutex.lock();
		if (tag == "gfct")
		{
			_map_time_egmcoef[ipt].c[n][m] = c;
			_map_time_egmcoef[ipt].s[n][m] = s;
			_map_time_egmcoef[ipt].sigma_c[n][m] = sigma_c;
			_map_time_egmcoef[ipt].sigma_s[n][m] = sigma_s;
		}
		else if (tag == "trnd")
		{
			_map_time_egmcoef[ipt].ctrnd[n][m] = c;
			_map_time_egmcoef[ipt].strnd[n][m] = s;
		}
		else if (tag == "acos")
		{
			if (p == 1)
			{
				_map_time_egmcoef[ipt].ccos1[n][m] = c;
				_map_time_egmcoef[ipt].scos1[n][m] = s;
			}
			else if (p == 0.5)
			{
				_map_time_egmcoef[ipt].ccos2[n][m] = c;
				_map_time_egmcoef[ipt].scos2[n][m] = s;
			}
		}
		else if (tag == "asin")
		{
			if (p == 1)
			{
				_map_time_egmcoef[ipt].csin1[n][m] = c;
				_map_time_egmcoef[ipt].ssin1[n][m] = s;
			}
			else if (p == 0.5)
			{
				_map_time_egmcoef[ipt].csin2[n][m] = c;
				_map_time_egmcoef[ipt].ssin2[n][m] = s;
			}
		}
		else
		{
			_log->comment(0, "ERROR : unknow time gravity coefficient:" + tag);
			_gmutex.unlock();
			return 0;
		}
		_gmutex.unlock();
		return n;
	};
	/**
	* @brief find the index of group which the time coefficient record belongs to (judge by the start time)
	* @param [in]  t     start time of record
	* @return	   ipt   index of group
	*/
	int t_gegmcoef::find_group(t_gtime t)
	{
		gtrace("t_gegmcoef::find_group");
		_gmutex.lock();
		int ipt = 0;
		map<int, t_time_gravity_coef>::iterator it = _map_time_egmcoef.begin();
		for (; it != _map_time_egmcoef.end(); it++)
		{
			if (it->second.st == t)
			{
				ipt = it->first;
				_gmutex.unlock();
				return ipt;
			}
		};
		_gmutex.unlock();
		return 0;
	};

	/**
	* @brief get the time variable spherical harmonic coefficients according to selected mjd
	* @param [in]  dmjd  mjd
	* @return	   bool
	*/
	const t_time_gravity_coef& t_gegmcoef::get_time_coefficient(double dmjd)
	{
		gtrace("t_gegmcoef::get_time_coefficient");
		_gmutex.lock();
		double st = 0.0, et = 0.0;
		map<int, t_time_gravity_coef>::iterator it = _map_time_egmcoef.begin();
		for (; it != _map_time_egmcoef.end(); it++)
		{
			st = it->second.st.dmjd();
			et = it->second.et.dmjd();
			if (st <= dmjd && dmjd <= et)
			{
				_gmutex.unlock();
				return it->second;
				break;
			}
		};
		_gmutex.unlock();
		return t_time_gravity_coef();
	}

}