/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li		2019-01-06	creat the file.
	  -1.1  Jie Li		2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gorb.cpp
 * @brief		The base class used to storage the orb data.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include "gdata/gorb.h"
#include "gmodels/gpoly.h"

using namespace std;

namespace great {

	/** @brief constructor. */
	t_gorb::t_gorb()
	{
		gtrace("t_orb::constructor");
		id_type(t_gdata::ORB);
	}

	/** @brief destructor. */
	t_gorb:: ~t_gorb()
	{

	}



	/**
	* @brief get orb head data.
	* @param[out]   orb_hd	orb head data
	*/
	void t_gorb::get_head(t_orb_hd& orb_hd)
	{
		// add jdhuang
		set<string> orb_data_sat_list;
		for (const auto& iter : _funct)
		{
			orb_data_sat_list.insert(iter.first);
		}

		_orb_hd.sat_prn = orb_data_sat_list;
		_orb_hd.sat_num = orb_data_sat_list.size();

		orb_hd = _orb_hd;
	}

	/**
	* @brief get orb funct data.
	* @param[out]   funct	whole orb funct data
	*/
	void t_gorb::get_funct(map<string, map<int, vector< t_gtriple> > >& funct)
	{
		funct = _funct;
	}

	/**
	* @brief get orb funct data of a satellite.
	* @param[in]   sat		satellite prn
	* @param[out]  funct	orb funct data of the sat
	*/
	void t_gorb::get_funct(string sat, map<int, t_sat_funct >& funct)
	{
		funct = _funct[sat];
	}


	void t_gorb::get_parlist(map<string,vector<string> >& parlist)
	{
		parlist = _par_list;
	}
}