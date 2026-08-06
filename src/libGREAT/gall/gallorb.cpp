#include "gallorb.h"
/* ----------------------------------------------------------------------
Purpose: save data of orb file
Version: $ Rev: $ 1.0

2019-04-10 /HJ_Zheng: created

-*/
/**
 *
 * @verbatim
	 History
	  -1.0  Hongjie Zheng   2019-04-10 Creat the file.
	  -1.1  Jie Li			2019-07-04 add doxygen style code remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gorb.h
 * @brief		The class for storaging all orb data.
 *
 * @author      Hongjie Zheng, Wuhan University
 * @version		1.0.0
 * @date		2019-04-10
 *
 */

#include "gall/gallorb.h"
#include "gmodels/gpoly.h"
#include "gutils/ginfolog.h"

#include <set>

using namespace std;
using namespace gnut;

namespace great
{

	/** @brief default constructor. */
	t_gallorb::t_gallorb() :
		t_gdata()
	{
		gtrace("t_gallorb::constructor");
		id_type(t_gdata::ALLORB);
	}
	/** @brief default destructor. */
	t_gallorb:: ~t_gallorb()
	{

	}

	void t_gallorb::getAllPartial(map<string, map<int, t_sat_funct>>& funct)
	{
		for (const auto& item : _orbs) {
			funct.insert(make_pair(item.first, item.second.second));
		}
	}

	void t_gallorb::getFirstHead(t_orb_hd & head)
	{
		head = _orbs.begin()->second.first;
	}

	

	/**
	* @brief interpolate a satellite's position(unit: m) in specific epoch.
	* @param[in]   epoch	epoch time
	* @param[in]   sat		satellite prn
	* @param[out]  px_sat	x coordinates
	* @param[out]  py_sat	y coordinates
	* @param[out]  pz_sat	z coordinates
	* @return   if interpolate successfully
	*	@retval   true   interpolate successfully
	*   @retval   false  interpolate unsuccessfully
	*/
	bool t_gallorb::get_crs_sat_crd(const t_gtime& epoch, const string& sat, double& px_sat, double& py_sat, double& pz_sat) 
	{

		if (_orbs.find(sat) == _orbs.end())
		{
			return false;
		}
		// NOTICE::lagrange order = 8
		const t_orb_hd& hd_info = _orbs[sat].first;
		map<int, vector<t_gtriple> >& funct = _orbs[sat].second;

		// NOTICE::lagrange order = 8
		const int poly_order = 8;
		int epo_left, epo_right;
		double epo_now;
		if (!getInterpPoints(epoch, poly_order, sat, epo_left, epo_right, epo_now))
		{
				return false;
		}

		map<int, t_sat_funct> ::iterator itepo = funct.begin();
		if (epo_left < itepo->first)
		{
			write_log_info(_log, 0, "ERROR", "Cannot perform the lagrange interpolation: input epoch is beyond the epoch table of orb file");
			return false;
		}
		map<int, t_sat_funct> ::reverse_iterator ritepo = funct.rbegin();
		if (epo_right > ritepo->first)
		{
			write_log_info(_log, 0, "ERROR", "Cannot perform the lagrange interpolation: input epoch is beyond the epoch table of orb file");
			return false;
		}

		// cycle for xyz 
		double	xyz[3];
		try
		{
			vector<double> Y;
			vector<double> X;
			Y.resize(epo_right - epo_left + 1);
			X.resize(epo_right - epo_left + 1);
			for (int crd = 0; crd < 3; crd++) 
            {
				for (int i = epo_left; i <= epo_right; i++) {
					Y[i - epo_left] = funct.at(i).at(0)[crd];
					X[i - epo_left] = i;
				}
				xyz[crd] = lagrange_interpolate(X, Y, epo_now);
			}
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << endl;
			return false;
		}
		// result
		px_sat = xyz[0] * 1000;
		py_sat = xyz[1] * 1000;
		pz_sat = xyz[2] * 1000;

		return true;

	}

	/**
	* @brief interpolate a satellite's velocity(unit: m/s) in specific epoch.
	* @param[in]   epoch	epoch time
	* @param[in]   sat		satellite prn
	* @param[out]  vx_sat	x velocity
	* @param[out]  vy_sat	y velocity
	* @param[out]  vz_sat	z velocity
	* @return   if interpolate successfully
	*	@retval   true   interpolate successfully
	*   @retval   false  interpolate unsuccessfully
	*/
	bool t_gallorb::get_crs_sat_vel(const t_gtime& epoch, string sat, double& vx_sat, double& vy_sat, double& vz_sat) {

		if (_orbs.find(sat) == _orbs.end())
		{
			return false;
		}
		// NOTICE::lagrange order = 8
		t_orb_hd& orb_hd = _orbs[sat].first;
		auto& orb_data = _orbs[sat].second;
		const int poly_order = 8;
		int epo_left, epo_right;
		double epo_now;
		if (!getInterpPoints(epoch, poly_order, sat, epo_left, epo_right, epo_now))
		{
				return false;
		}
	
		map<int, t_sat_funct> ::iterator itepo = orb_data.begin();
		if (epo_left < itepo->first)
		{
			write_log_info(_log, 0, "ERROR", "Cannot perform the lagrange interpolation: input epoch is beyond the epoch table of orb file");
			return false;
		}
		map<int, t_sat_funct> ::reverse_iterator ritepo = orb_data.rbegin();
		if (epo_right > ritepo->first)
		{
			write_log_info(_log, 0, "ERROR", "Cannot perform the lagrange interpolation: input epoch is beyond the epoch table of orb file");
			return false;
		}


		// cycle for xyz 
		double	vel[3];
		try
		{
			vector<double> Y;
			vector<double> X;
			Y.resize(epo_right - epo_left + 1);
			X.resize(epo_right - epo_left + 1);
			for (int crd = 0; crd < 3; crd++) {
				
				for (int i = epo_left; i <= epo_right; i++) {
					
					Y[i - epo_left] = orb_data[i][0][crd];
					X[i - epo_left] = i;
				}
				vel[crd] = lagrange_interpolate(X, Y, epo_now, true);
				vel[crd] /= orb_hd.admout;
			}
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << endl;
			return false;
		}
		// result
		vx_sat = vel[0] * 1000;
		vy_sat = vel[1] * 1000;
		vz_sat = vel[2] * 1000;

		return true;

	}

	

	/**
	* @brief interpolate a satellite's partial matrix(only dxdx0..dxdvx0...dxdnp) in specific epoch.
	* @param[in]   epo		epoch time
	* @param[in]   sat		satellite prn
	* @param[out]  partial	partial matrix
	* @return   if interpolate successfully
	*	@retval   true   interpolate successfully
	*   @retval   false  interpolate unsuccessfully
	*/
	bool t_gallorb::getGnsPartial(const t_gtime& epoch, string sat, Matrix& partial) 
	{
		if (_orbs.find(sat) == _orbs.end())
		{
			return false;
		}
		// NOTICE::lagrange order = 8
		const t_orb_hd& hd_info = _orbs[sat].first;
		map<int, vector<t_gtriple> >& funct = _orbs[sat].second;

		const int poly_order = 8;
		int epo_left, epo_right;
		double epo_now;
		if (!getInterpPoints(epoch, poly_order, sat, epo_left, epo_right, epo_now))
		{
		    return false;
		}

		if (funct.find(epo_left) == funct.end() || funct.find(epo_right) == funct.end())
		{
			return false;
		}

		// cycle for funct 
		partial.ReSize(funct[0].size() - 1, 3);
		try
		{
			vector<double> Y;
			vector<double> X;
			Y.resize(epo_right - epo_left + 1);
			X.resize(epo_right - epo_left + 1);
			for (int crd = 0; crd < 3; crd++) 
			{
				for (int par_num = 1; par_num < funct[0].size(); par_num++)
				{
					
					for (int i = epo_left; i <= epo_right; i++) 
					{
						
						Y[i - epo_left] = funct[i][par_num][crd];
						X[i - epo_left] = i;
					}
					partial(par_num, crd + 1) = lagrange_interpolate(X, Y, epo_now);
				}
			}
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << endl;
			return false;
		}
		// result
		return true;
	}


	/**
	* @brief get the interpolate point
	* @param[in]   epo			epoch input(int)
	* @param[in]   norder		lagrange interpolate order
	* @param[in]   sat			satellite prn
	* @param[out]  epo_left		first interpolate point
	* @param[out]  epo_right	last interpolate point
	* @param[out]  epo_now		current epoch time(mjd)
	* @return   if get interpolate point successfully
	*	@retval   true   get interpolate poin successfully
	*   @retval   false  get interpolate poin unsuccessfully
	*/
	bool t_gallorb::getInterpPoints(const t_gtime& epoch, const int& norder, const string& sat, int& epo_left, int& epo_right, double& epo_now)
	{
		if (_orbs.find(sat) == _orbs.end())
		{
			return false;
		}
		t_orb_hd& orb_hd = _orbs[sat].first;
		epo_now = (epoch - orb_hd.t_ref) / (double)orb_hd.admout;
		epo_left = round(epo_now - (int)(norder / 2) + 1.e-5 / orb_hd.admout);
		epo_right = epo_left + norder;
		
		const map<int, t_sat_funct >& funct = _orbs[sat].second;
		if (funct.find(epo_left) == funct.end() || funct.find(epo_right) == funct.end())
		{
			return false;
		}
		return true;
	}

	

	void t_gallorb::test_update(t_gtime beg_time, t_gtime end_time, int sampling,map<string, map<t_gtime, ColumnVector>>& newSatPos){
		set<string> sat_list;
		for (const auto& item : newSatPos) {
			sat_list.insert(item.first);
		}
		for (auto iter = _orbs.begin(); iter != _orbs.end();) {
			string sat = iter->first;
			if (newSatPos.count(sat) == 0) {
				iter = _orbs.erase(iter);
				continue;
			}
			// update header
			iter->second.first.beg = beg_time;
			iter->second.first.end = end_time;
			iter->second.first.t_ref = beg_time;
			iter->second.first.admout = sampling;
			iter->second.first.adminp = iter->second.first.admout;
			iter->second.first.sat_num = newSatPos.size();
			iter->second.first.sat_prn = sat_list;

			iter->second.second.clear();
			for (int i =0;i<=(end_time-beg_time)/iter->second.first.admout;i++){
				t_gtime now = beg_time+i*iter->second.first.admout;
				// poly insert
				if (newSatPos[sat].find(now) == newSatPos[sat].end()) {
					continue;
				}
				iter->second.second[i].push_back(t_gtriple(newSatPos[sat][now](1),newSatPos[sat][now](2),newSatPos[sat][now](3)));
				for (int j =0;j<11;j++){
					iter->second.second[i].push_back(t_gtriple());
				}
			}
			iter++;
		}
	}

	

	bool t_gallorb::update(std::string sat, map<int, vector<t_gtriple>> new_orb, t_gtime beg, t_gtime end, t_gtime t_ref, int oi_admout, bool keep_satpos) {
		_gmutex.lock();
		// Set Orb Header
		_orbs[sat].first.beg = beg;
		_orbs[sat].first.end = end;
		_orbs[sat].first.t_ref = t_ref;
		_orbs[sat].first.adminp = oi_admout;
		_orbs[sat].first.admout = oi_admout;
		// Set Orb body(Sat Pos)
		if (keep_satpos) {
			map<int, vector<t_gtriple> >& funct = _orbs[sat].second;
			for (auto iter = new_orb.begin(); iter != new_orb.end(); iter++) {
				iter->second[0][0] = funct[iter->first + 1][0][0];
				iter->second[0][1] = funct[iter->first + 1][0][1];
				iter->second[0][2] = funct[iter->first + 1][0][2];
			}
		}
		_orbs[sat].second = new_orb;
		_gmutex.unlock();
		return true;
	}

}
