/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBlRec.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambBlRec.h"
#include "gambfix/gambBl.h"
#include "gall/gallobj.h"
#include "gset/gsetamb.h"
#include "gset/gsetrec.h"
#include "gutils/gstring.h"
#include <algorithm>
#include <chrono>
#include "newmat/newmatap.h"

using namespace std;

namespace great
{
	t_gamb_baseline_recs::t_gamb_baseline_recs()
	{
	}





	void t_gamb_baseline_recs::add_cut(const t_gtime& cut_time, double cut_length)
	{
		_cut_length = cut_length;
		_cut_time = cut_time;
	}

	void t_gamb_baseline_recs::add_gobj(t_gallobj* obj)
	{
		this->_gobj = obj;

	}

	void t_gamb_baseline_recs::add_gset(t_gsetbase* set)
	{
		this->_gset = set;
		_beg_time = dynamic_cast<t_gsetgen*>(set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(set)->end();

	}

	void t_gamb_baseline_recs::add_glog(t_glog* log)
	{
		this->_glog = log;
	}

	void t_gamb_baseline_recs::add_rec(const string& rec, const string& sat)
	{
		if (std::find(_recs.begin(),_recs.end(),rec) == std::end(_recs))
		{
			_recs.push_back(rec);
			std::sort(_recs.begin(), _recs.end());
		}
		
		if (_sys[rec].count(sat.substr(0, 1)) < 1) _sys[rec].insert(sat.substr(0, 1));
	}

	vector<pair<string, string>> t_gamb_baseline_recs::baselines(bool check_depend)
	{

		vector<pair<string, string>> tmp;
		
		if (_baselines.empty())
		{
			this->_get_rec_xyz();


			chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();

			this->_create_baseline(check_depend);
			chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
		}

		for (const auto& iter : _baselines)
		{
			tmp.push_back(iter.bl_pair());
		}


		// cout the base line
		for (const auto& iter : tmp)
		{
			_glog->comment(0, "The baseline is " + iter.first + "   " + iter.second);
		}

		return tmp;

	}


	bool t_gamb_baseline_recs::_check_amb_depend()
	{
		// rec map
		map<string, int> _obj_index;
		for (int i = 0; i < _recs.size(); i++)
		{
			_obj_index[_recs[i]] = i + 1;
		}

		// judge
		vector<t_gamb_baseline> Independent_baseline;
		Matrix	vector_space;
		for (const auto& iter : _baselines)
		{
			const auto& line = iter.bl_pair();
			const auto& tmp1 = _obj_index.count(line.first);
			const auto& tmp2 = _obj_index.count(line.second);
			if (tmp1 <= 0 || tmp2 <= 0) continue;

			int ncol = _obj_index.size();
			if (vector_space.is_zero())
			{
				vector_space.resize(1, ncol);
				vector_space = 0.0;
				int index1 = _obj_index[line.first];
				int index2 = _obj_index[line.second];
				vector_space(1, index1) = 1;
				vector_space(1, index2) = -1;
				
				Independent_baseline.push_back(iter);
				continue;
			}

			RowVector new_vector(0.0, ncol);
			int index1 = _obj_index[line.first];
			int index2 = _obj_index[line.second];
			new_vector(index1) = 1;
			new_vector(index2) = -1;

			// 
			Matrix A = vector_space & new_vector;
			DiagonalMatrix D;
			if (A.nrows() < A.ncols())
			{
				SVD(A.t(), D);
			}
			else
			{
				SVD(A, D);
			}
			bool isDepend = (fabs(D(D.nrows(), D.ncols())) > 0.01);

			if (A.nrows() != (vector_space.nrows() + new_vector.nrows())) return false;
			if (isDepend)
			{
				vector_space = A;	
				Independent_baseline.push_back(iter);
			}
		}

		_baselines = Independent_baseline;

		return true;
	}

	bool t_gamb_baseline_recs::_get_rec_xyz()
	{
		if (this->_gobj == nullptr && this->_gset == nullptr) return false;

		for (const auto& crt : this->_recs)
		{
			bool   isfind = false;
			const  string&  rec = str2upper(crt);

			if (this->_gset)
			{
				t_gtriple tmp = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(rec);
				if (!tmp.zero())
				{
					_glog->logInfo(rec + "  find xyz in xml file  ");
					this->_xyz[rec] = tmp;
					isfind = true;
				}
			}

			if (isfind == false && this->_gobj)
			{
				shared_ptr<t_gobj> grec = this->_gobj->obj(rec);
				if (grec)
				{
					t_gtime tmp = _cut_time;
					while (tmp <= _end_time)
					{
						t_gtriple rec_crd;
						t_gtriple rec_std;
						grec->get_recent_crd(_cut_time, 100, rec_crd, rec_std);
						if (!rec_crd.zero())
						{
							_glog->logInfo(rec + "  find xyz in snx file  ");
							this->_xyz[rec] = rec_crd;
							isfind = true;
						}

						if (isfind) break;
						tmp = tmp + 300;
					}
				}
			}
		}
		return true;
	}

	bool t_gamb_baseline_recs::_create_baseline(bool check_depend)
	{
		if (_recs.empty())
		{
			cout << "  rec list is empty, return false  " << endl;
			return false;
		}

		const auto& loop_1st = _recs;
		const auto& loop_2nd = _recs;

		int index = 0;
		vector<string> tmp;
		for (const auto& iter : _recs)
		{
			if (this->_xyz.count(iter) == 0) continue;
			tmp.push_back(iter);
		}
		_recs = tmp;



		//
		for (int i = 0; i < _recs.size(); i++)
		{
			const string& rec_1st = _recs[i];
			//
			if (_sys[rec_1st].size() == 0)
			{
				cout << "  no system for " + rec_1st + " skiped" << endl;
				continue;
			}
			for (int j = i + 1; j < _recs.size(); j++)
			{
				const string& rec_2nd = _recs[j];
				if (_sys[rec_2nd].size() == 0) continue;

			
				t_gtriple first  = this->_xyz[rec_1st];
				t_gtriple second = this->_xyz[rec_2nd];
				t_gtriple D_value = first - second;
				double distance = D_value.norm() / 1000.0;
				if (distance > fabs(_cut_length)) continue;
				int nsys = 0;
				for (const auto& iter : _sys[rec_1st])
				{
					if (_sys[rec_2nd].find(iter) != std::end(_sys[rec_2nd])) nsys++;
				}
				if (nsys == 0) continue;
				t_gamb_baseline bl(rec_1st, rec_2nd, distance, nsys);
				_baselines.push_back(bl);
			}
		}
		if (_baselines.empty()) return false;
		std::sort(_baselines.begin(), _baselines.end(), compare);

		
		if (check_depend){
			_check_amb_depend();
		}


		return true;
	}

	bool compare(const t_gamb_baseline& bl1, const t_gamb_baseline& bl2)
	{
		if (bl1.nsys > bl2.nsys)
		{
			return true;
		}

		if ((bl1.nsys == bl2.nsys) && (bl1.distance < bl2.distance))
		{
			return true;
		}

		return false;
	}

}

