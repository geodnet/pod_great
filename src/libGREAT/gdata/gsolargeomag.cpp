/**
*
* @verbatim
    History
     -1.0 Yujie Qin 2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsolargeomag.cpp
* @brief		    Storage the LEO solar_flux and geomag_kp data for one satellite
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/


#include "gsolargeomag.h"
using namespace std;

namespace great {
	/** @brief default constructor. */
	t_gsolargeomag::t_gsolargeomag()
		:t_gdata()
	{
		id_type(t_gdata::ALLSOLAR);
	}

	/**
	* @brief set one record with time and corresponding value(solar_flux/geomag_kp).
	* @param[in]    datatype      data type(solar/geomag)
	* @param[in]    mjdtime       modified juliay day
	* @param[in]    value         value in the requested time
	*/
	void t_gsolargeomag::set_data(string datatype, double mjdtime, double* value)
	{
		if (datatype == "solar_flux" || datatype == "solar_flux_MSISE")
		{
			_solar_data[mjdtime][0] = value[0];
			_solar_data[mjdtime][1] = value[1];
		}
		else if (datatype == "geomag_kp" || datatype == "geomag_ap")
		{
			_geomag_gdata[mjdtime] = value[0];
		}
		else
		{
			if (_log) _log->comment(1, "WARNING : t_gsolargeomag::setdata Wrong filetype");
			else       cerr << "warning: t_gsolargeomag::setdata unknown filetype:" << endl;
		}
	}

	/**
	* @brief set the interval of data(solar_flux/geomag_kp).
	* @param[in]  datatype   data type(solar/geomag)
	* @param[in]  interval   data interval
	*/
	void t_gsolargeomag::set_intv(string datatype, double interval)
	{
		if (datatype == "solar_flux" || datatype == "solar_flux_MSISE")
		{
			_solar_itv_time = interval;
		}
		if (datatype == "geomag_kp"  || datatype == "geomag_ap")
		{
			_geomag_itv_time = interval;
		}
	}

	/**
	* @brief set begin and end time of data(solar_flux/geomag_kp)
	* @param[in]   datatype  data type(solar/geomag)
	* @param[in]   begin     begin time
	* @param[in]   end       end time
	*/
	void t_gsolargeomag::set_time(string datatype, double begin, double end)
	{
		if (datatype == "solar_flux" || datatype == "solar_flux_MSISE")
		{
			_solar_beg_time = begin;
			_solar_end_time = end;
		}
		if (datatype == "geomag_kp" || datatype == "geomag_ap")
		{
			_geomag_beg_time = begin;
			_geomag_end_time = end;
		}
	}


}