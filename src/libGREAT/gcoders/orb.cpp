/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li			2019-01-06	creat the file.
	  -1.1  Jie Li			2019-04-01	add code in doxygen style.
	  -1.2  Jie Li			2019-04-02  add abnormal throwing mechanism of functions
	  -1.3  Hongjie Zheng	2019-04-10  add allorb decode
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		orb.cpp
 * @brief		The base class used to decode and encode orb file information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include <algorithm>
#include "gcoders/orb.h"
#include "gdata/gorb.h"
#include "gall/gallorb.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_orb_hd::t_orb_hd()
	{

	}

	/** @brief destructor. */
	t_orb_hd::~t_orb_hd()
	{

	}
	
	/**
	* @brief constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_orb::t_orb(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_orb::constructor");
	}

	/** @brief destructor. */
	t_orb:: ~t_orb(){}

	
	int t_orb::encode_head(char* buff, int sz, vector<string>& errmsg)
	{

		_mutex.lock();

		try
		{
			if (_ss_position == 0)
			{
				map<string, t_gdata*>::iterator it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::ORB)
					{
						//get orb head info
						dynamic_cast<t_gorb*>(it->second)->get_head(_orb_hd);
					}
					if (it->second->id_type() == t_gdata::ALLORB)
					{
						dynamic_cast<t_gallorb*>(it->second)->getFirstHead(_orb_hd);
					}
				}
				//fill head data
				_ss << "Sat Group: " << _orb_hd.sat_group << endl;

				set<char>::iterator ithd = _orb_hd.sat_type.begin();
				_ss << "Sat Type: ";
				for (ithd; ithd != _orb_hd.sat_type.end(); ++ithd)
				{
					_ss << *ithd << " ";
				}
				_ss << endl;
				//sat list
				_ss << "Sat List: " << _orb_hd.sat_num << endl;
				set<string>::iterator itsat = _orb_hd.sat_prn.begin();
				int isat = 0;
				for (itsat; itsat != _orb_hd.sat_prn.end(); ++itsat)
				{
					_ss << "  " << *itsat;
					isat++;
					if (isat == 10)
					{
						_ss << endl;
						isat = 0;
					}
				}
				_ss << endl;
				//intege contral
				string tmp_partial;
				if (_orb_hd.partial)
				{
					tmp_partial = "true";
				}
				else
				{
					tmp_partial = "false";
				}
				_ss << "Intege Control: " << tmp_partial << "  " << _orb_hd.admord << "  " << _orb_hd.adminp << "  " << _orb_hd.admout << endl;
				//Nequ/Nepo
				_ss << "Nequ: " << _orb_hd.nequ << endl;
				_ss << "Nepo: " << _orb_hd.nepo << endl;
				//start time
				int mjd = _orb_hd.beg.mjd();
				double sod = _orb_hd.beg.sod() + _orb_hd.beg.dsec();
				_ss << "Beg Time: " << mjd << " " << sod << endl;
				//end time
				mjd = _orb_hd.end.mjd();
				sod = _orb_hd.end.sod() + _orb_hd.end.dsec();
				_ss << "End Time: " << mjd << " " << sod << endl;
				//ref time
				mjd = _orb_hd.t_ref.mjd();
				sod = _orb_hd.t_ref.sod() + _orb_hd.t_ref.dsec();
				_ss << "Ref Time: " << mjd << " " << sod << endl;

				_ss << "END OF HEAD" << endl;
			}

			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_orb::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_orb::encode_head throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief encode data body of orb file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return  size for data body encoding
	*/
	int t_orb::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{

		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				//get data from _data
				//funct data
				map<string, map<int, t_sat_funct > > funct;
				map<string, t_gdata*>::iterator it = _data.begin();
				map<string, vector<string> > parlist;
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::ORB)
					{
						dynamic_cast<t_gorb*>(it->second)->get_funct(funct);
						dynamic_cast<t_gorb*>(it->second)->get_parlist(parlist);
					}
					if (it->second->id_type() == t_gdata::ALLORB) {
						dynamic_cast<t_gallorb*>(it->second)->getAllPartial(funct);
						for (const auto&item : funct) {
							parlist.insert(make_pair(item.first, vector<string>({
								"PXSAT","PYSAT","PZSAT","VXSAT","VYSAT","VZSAT",
								"D0","Y0","X0","Xc","Xs"
								})));
						}
					}
				}

				//encode
				//fill epo funct data
				map<string, map<int, t_sat_funct > >::iterator itFunct = funct.begin();
				//map<int, t_sat_funct> ::iterator itEpo = (funct.begin()->second).begin();
				//circle by epoch
				int mul_step = _orb_hd.admout / _orb_hd.adminp;
				int nEpo_beg = (_orb_hd.beg - _orb_hd.t_ref) / _orb_hd.admout;
				int nEpo_end = (_orb_hd.end - _orb_hd.t_ref) / _orb_hd.admout;
				int nEpo = nEpo_beg;
				for (nEpo; nEpo <= nEpo_end; nEpo++) {

					_ss << "*" << nEpo << endl;

					int nEpo_inp = nEpo * mul_step;

					//circle by sat
					for (itFunct; itFunct != funct.end(); ++itFunct)
					{
						auto sat_par = parlist[itFunct->first];
						//circle by a record
						t_sat_funct rec = (itFunct->second)[nEpo_inp];
						//write sat pos
						if (rec.size() == 0)
						{
							continue;
						}
						_ss << itFunct->first << endl;
						_ss << "   " << setw(20) << left << "SAT_POS";
						for (int j = 0; j < 3; j++)
						{
							_ss << setw(25) << setprecision(14) << scientific << right << rec[0][j];
						}
						_ss << endl;
						if (sat_par.size() != (rec.size() - 1))
						{
							throw;
						}
						//write partial
						for (int i = 0; i < sat_par.size(); i++)
						{
							_ss << "   " << setw(20) << left << "dPOS/d" + sat_par[i];
							for (int j = 0; j < 3; j++)
							{
								_ss << setw(25) << setprecision(14) << scientific << right << rec[i + 1][j];
							}
							_ss << endl;
						}
						_ss << "END OF SAT" << endl;

					}
					itFunct = funct.begin();
				}
				_ss << "END OF FILE" << endl;
			}
			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_orb::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_orb::encode_data throw exception" << endl;
			}
			return -1;
		}
	}

}