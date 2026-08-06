
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <assert.h>

#include "gall/gallpar.h"

using namespace std;

namespace gnut {

	// Add t_gpar to t_gallpar
	// Parameter is stored at the end of the vector
	// -----------------------------------------------
	void t_gallpar::addParam(t_gpar newPar)
	{
		gtrace("t_gallpar::addParam");
		this->_vParam.push_back(newPar);
		this->_point_par.push_back(_max_point++);
		this->_index_par[newPar.get_head()][newPar.get_timearc()] = this->_point_par[this->_point_par.size() - 1];
	}

	void t_gallpar::addParamHead(t_gpar newPar)
	{
		gtrace("t_gallpar::addParam");

		this->_vParam.insert(this->_vParam.begin(),newPar);
		this->_point_par.insert(this->_point_par.begin(),_min_point--);
		this->_index_par[newPar.get_head()][newPar.get_timearc()] = this->_point_par[0];
	}



	// Delete paremeter
	// Parameter is deleted according to index value
	// ----------------------------------------------------
	void t_gallpar::delParam(int i)
	{
		gtrace("t_gallpar::delParam");

		auto& all = this->_index_par.at(_vParam[i].get_head());
		for (auto iter = all.begin(); iter != all.end(); iter++) 
		{
			if (iter->second == _point_par[i]) 
			{
				all.erase(iter);
				break;
			}
		}
		if (this->_index_par[_vParam[i].get_head()].size() == 0) 
		{
			this->_index_par.erase(_vParam[i].get_head());
		}
		_point_par.erase(_point_par.begin() + i);

		_vParam.erase(_vParam.begin() + i);
	}

	// get position of item according: station name, par type, PRN, begin time, end time
	// -----------------------------------------------------
	int t_gallpar::getParam(string mark, t_gpar::t_type type, string prn,
		t_gtime beg, t_gtime end) const
	{
		auto find_iter = this->_index_par.find(t_gparhead(type, mark, prn));
		if (find_iter == this->_index_par.end()) {
			return -1;
		}
		const auto& all = find_iter->second;
		t_gtimearc dst_timearc(beg, end);

		auto all_end = all.end();
		auto par_beg = _point_par.begin();
		auto par_end = _point_par.end();
		for (auto iter = all.begin(); iter != all_end; iter++)
		{
			if (iter->first.inside(dst_timearc)) 
			{
				auto ans = lower_bound(par_beg, par_end, iter->second);
				assert(ans != par_end && *ans == iter->second);
				return ans - par_beg;
			}
		}
		return -1;
	}


	int t_gallpar::getParIndex(int idx)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			return _vParam[idx].index;
		}
		else
		{
			return -1;
		}
	}

	// do not change, we use the function only for add ambcon
	// the beg and end of par should be Contained in the beg and end of input
	int t_gallpar::getAmbParam(string site, string prn, t_gpar::t_type  type, const t_gtime& beg, const t_gtime& end) const
	{
		gtrace("t_gallpar::getAmbParam");

		if (this->_index_par.count(t_gparhead(type, site, prn)) == 0) {
			return -1;
		}
		else {
			auto all = this->_index_par.find(t_gparhead(type, site, prn))->second;
			t_gtimearc dst_timearc(beg, end);
			//auto iter = lower_bound(all.begin(),all.end(), dst_timearc);
			for (auto iter = all.begin(); iter != all.end(); iter++) {
				auto ans = lower_bound(_point_par.begin(), _point_par.end(), iter->second);
				int par_index = ans -  _point_par.begin();
				if (dst_timearc.inside(this->_vParam[par_index].get_timearc())){
					return par_index;
				}

			}
			return -1;
		}
		return -1;

	}
	// used with t_gallpar::getAmbParam
	t_gpar& t_gallpar::getAmbParam(int idx)
	{
		return _vParam[idx];
	}

	double t_gallpar::getParValue(int idx)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			return _vParam[idx].value();
		}
		else
		{
			return 0.0;
		}
	}

	void t_gallpar::setParValue(int idx, double value)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			_vParam[idx].value(value);
		}
	}



	// Get number of parametres
	// ---------------------------------------------
	unsigned int t_gallpar::parNumber() const
	{
		gtrace("t_gallpar::parNumber");

		return _vParam.size();
	}


// 	// Get number of orbit parametres
   // ---------------------------------------------
	 unsigned int t_gallpar::orbParNumber() const
	 {
	 	gtrace("t_gallpar::orbParNumber");

	 	unsigned int orbparnum = 0;
	 	for (auto sat : _vOrbParam)
	 	{
	 		orbparnum += sat.second.size();
	 	}

	 	return orbparnum;
	 }


	// Strore Coordinates parametres to t_triple crd
	// ------------------------------------------------

	int t_gallpar::getCrdParam(string station, t_gtriple& crd, t_gtime Tbeg, t_gtime Tend) const
	{
		gtrace("t_gallpar::getCrdParam");

		int found = 0;

		int idx = this->getParam(station, t_gpar::CRD_X, "", Tbeg, Tend);
		int idy = this->getParam(station, t_gpar::CRD_Y, "", Tbeg, Tend);
		int idz = this->getParam(station, t_gpar::CRD_Z, "", Tbeg, Tend);

		if (idx != -1) {
			crd.set(0, _vParam[idx].value());
			found++;
		}
		if (idy != -1) {
			crd.set(1, _vParam[idy].value());
			found++;
		}
		if (idz != -1) {
			crd.set(2, _vParam[idz].value());
			found++;
		}

		if (found == 3)  return  1;    // all three crd were found
		if (found == 1)  return -1;    // just one crd was found
		if (found == 2)  return -2;    // just two crd was found

		if (found == 0) return  -3;    // crd not found

		return -1;
	}

	vector<int> t_gallpar::getPartialIndex(string site, string sat) 
	{
		_update_partial_index();

		vector<int> ans;
		pair<string, string> type_list[4] =
		{
			make_pair(site,sat),
			make_pair(site,""),
			make_pair("",sat),
			make_pair("","")
		};

		_allpar_mtx.lock();
		for (int i = 0; i < 4; i++)
		{
			ans.insert(ans.end(), _index_for_parital[type_list[i]].begin(), _index_for_parital[type_list[i]].end());
		}
		_allpar_mtx.unlock();

		return ans;
	}


	t_gpar& t_gallpar::getPar(int idx)
	{
		return _vParam[idx];
	};

	// get single t_gpar element from container
	// ----------------------------------------------
	t_gpar& t_gallpar::operator[](const size_t idx)
	{
		return _vParam[idx];
	}

	
	//Doplnil Gabo
	void t_gallpar::delAllParam()
	{
		gtrace("t_gallpar::delAllParam");

		_vParam.clear();
		this->_index_par.clear();
        this->_index_for_parital.clear();
		this->_point_par.clear();
		this->_max_point = 0;
        this->_min_point = -1;
		this->_last_point = make_pair(0, 0);
	}


	
	set<string> t_gallpar::amb_prns(t_gpar::t_type type)
	{
		gtrace("t_gallpar::amb_prns");

		set<string> prns;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{
			// changed by wangbo
			if (iter->parType == type)
			{
				prns.insert(iter->prn);
			}
		}
		return prns;

	}

// 
// 	// Get prn of all AMB params
	// --------------
	set<string> t_gallpar::amb_prns()
	{
		gtrace("t_gallpar::amb_prns");

		set<string> prns;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{

			if(iter->str_type().find("AMB") != string::npos)
			{
				prns.insert(iter->prn);
			}
		}
		return prns;
	}

	set<string> t_gallpar::amb_recs(t_gpar::t_type type)
	{
		gtrace("t_gallpar::amb_recs");

		set<string> recs;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{
			if (iter->parType == type)
			{
				recs.insert(iter->site);
			}
		}
		return recs;
	}

	set<string> t_gallpar::amb_recs()
	{
		gtrace("t_gallpar::amb_recs");

		set<string> recs;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{

			if (iter->str_type().find("AMB") != string::npos)
			{
				recs.insert(iter->site);
			}
		}
		return recs;
	}

	// overloading << operator
	// -----------------------------
	ostream& operator<<(ostream& os, t_gallpar& x)
	{
		for (unsigned int i = 0; i < x.parNumber(); i++)
		{
			if (x[i].parType == t_gpar::AMB_IF ||
				x[i].parType == t_gpar::AMB_L1 ||
				x[i].parType == t_gpar::AMB_L2 ||
				x[i].parType == t_gpar::AMB_L3 ||
				x[i].parType == t_gpar::AMB_L4 ||
				x[i].parType == t_gpar::AMB_L5 ||
				x[i].parType == t_gpar::SION ||
				x[i].parType == t_gpar::VION ||
				x[i].parType == t_gpar::IFB_C3 ||
				x[i].parType == t_gpar::IFB_C4 ||
				x[i].parType == t_gpar::IFB_C5 ||
				x[i].parType == t_gpar::GPS_REC_IFB_C3||
				x[i].parType == t_gpar::GPS_REC_IFB_C4||
				x[i].parType == t_gpar::GPS_REC_IFB_C5||
				x[i].parType == t_gpar::SAT_IFB_C3||
				x[i].parType == t_gpar::SAT_IFB_C4||
				x[i].parType == t_gpar::SAT_IFB_C5||
				x[i].parType == t_gpar::IFCB_F3 ||
				x[i].parType == t_gpar::IFCB_F4 ||
				x[i].parType == t_gpar::IFCB_F5)
			{
				os << x[i].str_type() << "_" << x[i].prn << " ";
			}
			else
			{
				os << x[i].str_type() << " ";
			}

			if (x[i].parType == t_gpar::GRD_N || x[i].parType == t_gpar::GRD_E)
			{
				os << "value: " << x[i].value() * 1000 << " " << "index:" << x[i].index;
			}
			else
			{
				os << "value: " << x[i].value() << " " << "index:" << x[i].index;
			}
			os << endl;
		}
		return os;
	}


	void t_gallpar::_update_partial_index()
	{
		_allpar_mtx.lock();
		auto point_now = make_pair(_point_par[_point_par.size() - 1], int(_point_par.size() - 1));

		if (point_now!= _last_point)
		{
			_last_point = point_now;
			_index_for_parital.clear();
			set<string> site_list, sat_list;

			for (int i = 0; i < _vParam.size(); i++)
			{
				_index_for_parital[make_pair(_vParam[i].site, _vParam[i].prn)].push_back(i);
			}
		}
		_allpar_mtx.unlock();
		return;
	}




} // namespace
