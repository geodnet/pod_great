/**
*
* @verbatim
	History
	 -1.0 Keke Zhang	 2020-12-29 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			eigen6s4.cpp
* @brief		Decode EIGEN6S4 data
*
* @author		Keke Zhang, Wuhan University
* @version		1.0.0
* @date			2020-12-29
*
*/
#include "gcoders/eigen6s4.h"
#include "gdata/gegmcoef.h"
using namespace std;

namespace great
{

	/**
	* @brief constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_eigen6s4::t_eigen6s4(t_gsetbase *s, string version, int sz)
		: t_gcoder(s, version, sz)
	{
		gtrace("t_eigen6s4::constructor");
	}

	/** @brief destructor. */
	t_eigen6s4::~t_eigen6s4() {}

	/**
	* @brief decode header of EIGEN6S4 file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_eigen6s4::decode_head(char *buff, int sz, vector<string> &errmsg)
	{
		gtrace("t_eigen6s4::decode_head");

		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		}

		string tmp;
		int consume = 0;
		int tmpsize = 0;
		double gm = 0.0;
		double radius = 0.0;
		int i = 0, year = 0, mon = 0, day = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0)
			{
				stringstream ss(tmp);
				consume += tmpsize;
				if (tmp.find("+coefficient") == -1)
				{
					if (tmp.substr(0, 5) == "EARTH")
					{
						string str;
						gm = 0.0;
						radius = 0.0;
						i = 0;
						year = 0;
						mon = 0;
						day = 0;

						t_gtime mjd_ref;

						ss >> str >> gm >> radius >> i;
						gm = gm * 1e-9;
						radius = radius * 1e-3;
						year = i / 10000;
						mon = (i - year * 10000) / 100;
						day = i - year * 10000 - mon * 100;
						mjd_ref.from_ymd(year, mon, day, 0, 0.0);

						//fill time gm and radius
						map<string, t_gdata *>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::EGM)
							{
								((t_gegmcoef *)it->second)->add_head(mjd_ref, gm, radius);
							}
							it++;
						}
					}

					if (tmp.substr(0, 3) == "SHM")
					{
						string str;
						int order = 0.0;
						ss >> str >> order;

						//fill max order of egm file
						map<string, t_gdata *>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::EGM)
							{
								((t_gegmcoef *)it->second)->add_order(order);
							}
							it++;
						}
					}
					t_gcoder::_consume(tmpsize);
				}
				else
				{
					if (_log)
						_log->comment(2, "t_eigen6s4::decode_head end of head");
					t_gcoder::_consume(tmpsize);
					_mutex.unlock();
					return -1;
				}
			}
			_mutex.unlock();
			return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_eigen6s4::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_eigen6s4::decode_heada throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief decode data body of EIGEN6S4 file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_eigen6s4::decode_data(char *buff, int sz, int &cnt, vector<string> &errmsg)
	{
		gtrace("t_eigen6s4::decode_data");
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		};
		string tmp;
		int consume = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0)
			{
				stringstream ss(tmp);
				consume += tmpsize;
				if (tmp.find("-coefficient") == -1)
				{
					int n = 0, m = 0, ipt = 0;
					double c = 0.0, s = 0.0, sigma_c = 0.0, sigma_s = 0.0;
					double stime = 0.0, etime = 0.0, p = 0.0;
					int yr = 0, mn = 0, dd = 0, hr = 0, mi = 0;
					string tag;
					t_gtime st, et;
					ss >> tag >> n >> m >> c >> s >> sigma_c >> sigma_s >> stime >> etime >> p;
					if (stime != 0.0)
					{
						yr = stime / 10000;
						mn = (stime - yr * 10000) / 100;
						dd = stime - yr * 10000 - mn * 100;
						hr = (stime - floor(stime)) * 100;
						mi = (stime * 100 - floor(stime * 100)) * 100;
						st.from_ymdhms(yr, mn, dd, hr, mi, 0.0, true);

						yr = etime / 10000;
						mn = (etime - yr * 10000) / 100;
						dd = etime - yr * 10000 - mn * 100;
						hr = (etime - floor(etime)) * 100;
						mi = (etime * 100 - floor(etime * 100)) * 100;
						et.from_ymdhms(yr, mn, dd, hr, mi, 0.0, true);
					}

					//fill coeff data
					map<string, t_gdata *>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::EGM)
						{
							if (tag == "gfc")
								((t_gegmcoef *)it->second)->add_egmcoef(n, m, c, s, sigma_c, sigma_s);
							else
							{
								ipt = ((t_gegmcoef *)it->second)->find_group(st);
								if (ipt != 0)
								{
									((t_gegmcoef *)it->second)->add_timecoef(ipt, tag, n, m, c, s, sigma_c, sigma_s, p);
								}
								else
								{
									ipt = ((t_gegmcoef *)it->second)->create_timecoef_group(st, et);
									((t_gegmcoef *)it->second)->add_timecoef(ipt, tag, n, m, c, s, sigma_c, sigma_s, p);
								}
							}
						}
						it++;
					}
				}
				else
				{
					//cout << "end of data" << endl;
					_mutex.unlock();
					return -1;
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock();
			return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_eigen6s4::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_eigen6s4::decode_data throw exception" << endl;
			}
			return -1;
		}
	}
};