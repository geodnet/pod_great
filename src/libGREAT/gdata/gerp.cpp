/**
 *
 * @verbatim
	 History
	  -1.0	Yujie Qin	2019-01-06	creat the file.
	  -1.1  Jie Li		2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gerp.cpp
 * @brief		The base class used to storage the erp data.
 *
 * @author      Yujie Qin, Wuhan University
 * @version		1.0.0
 * @date		    2019-01-06
 */


#include "gdata/gerp.h"

using namespace std;
namespace great {

	/** @brief constructor. */
	t_gerp::t_gerp()
		:t_gdata()
	{
		id_type(t_gdata::ALLERP);
	}

	/**
	* @brief set begin and end time for erp data.
	* @param[in]   beg		begin time of the erp data
	* @param[in]   end		end time of the erp data
	*/
	void t_gerp::setBegEndTime(int beg, int end)
	{
		_beg_time = beg;
		_end_end = end;
	}

	/**
	* @brief add erp data.
	* @details add one erp data record to _erpdata and set timetype and interval.
	*
	* @param[in]   mjdtime		(mjd)time
	* @param[in]   data			one erp data record(xpole,ypole and so on)
	* @param[in]   mode			type of UT1(such as UT1R)
	* @param[in]   intv			interal(unit:day)
	*/
	void t_gerp::setEopData(t_gtime mjdtime, map<string, double> data, string mode, double intv)
	{
		_erp_data[mjdtime] = data;
		_UT1_mode = mode;
		_intv = intv;
	}

	/**
	* @brief get corresponding erp data by time.
	* @param[in]   mjdtime		(mjd)time
	* @param[out]   data			one erp data record(xpole,ypole and so on)
	*/
	void t_gerp::getEopData(t_gtime mjdtime, map<string, double>& data)
	{
		_gmutex.lock();
		if (_erp_data.find(mjdtime) != _erp_data.end())
		{
			data = _erp_data[mjdtime];
		}
		// interpolation
		else
		{
			auto iter_right = _erp_data.upper_bound(mjdtime);
			if (iter_right == _erp_data.begin() || iter_right == _erp_data.end())
			{
				cout << "interpolation erp wrong!!" << endl;
				_gmutex.unlock();
				throw runtime_error("interpolation erp wrong!!");
			}
			auto iter_left = iter_right--;

			double alpha = (mjdtime - iter_left->first) / (iter_right->first - iter_left->first);
			data.clear();
			for (auto iter = iter_left->second.begin(); iter != iter_left->second.end(); iter++)
			{
				double temp = iter_left->second[iter->first] + alpha * (iter_right->second[iter->first] - iter_left->second[iter->first]);
				data.insert(make_pair(iter->first, temp));
			}
		}
		_gmutex.unlock();
	}


	/** @brief whether the erp data is empty.
	* @return  bool
	*	@retval   true   erp data is empty
	*   @retval   false  erp data is existent
	*/
	bool t_gerp::isEmpty()
	{
		_gmutex.lock();
		if (_erp_data.size() == 0)
		{
			_gmutex.unlock();
			return true;
		}
		else
		{
			_gmutex.unlock();
			return false;
		}
	}
}