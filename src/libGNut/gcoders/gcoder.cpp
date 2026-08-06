
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <vector>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <algorithm>

#include "gcoders/gcoder.h" 
#include "gutils/gtypeconv.h"
#include "gutils/gfileconv.h"
#include "gutils/gtimesync.h"
#include "gset/gsetrec.h"

using namespace std;

namespace gnut
{

	t_gcoder::t_gcoder()
	{}

	// constructor
	// ----------
	t_gcoder::t_gcoder(t_gsetbase* s, string version, int sz, string id)
		: _class(id),
		_initialized(false),
		_gnss(true),
		_initver(version),
		_close_with_warning(true)
	{
		_init();

		// after init!
	  //  _version = version; 
		_buffsz = sz;
		_log = 0;
		_set = 0;
		_irc = 0;
		_hdr = false;

		_ss_position = 0;

		_out_smp = 0;
		_out_len = 0;
		_out_epo = LAST_TIME;

		// use malloc instead of new due to realocate function!
		// ====================================================
		_buffer = (char*)malloc((_buffsz + 1) * sizeof(char));  // due to realocate function!

		if (s) _gset(s); // READ SETTINGS IN CONSTRUCT ONLY!


		_pgm = "Anubis";
	}
	t_gcoder::t_gcoder(t_gtime beg, t_gtime end, t_gsetbase* s, string version, int sz, string id)
		: _class(id),
		_initialized(false),
		_gnss(true),
		_initver(version),
		_close_with_warning(true),
		beg_epoch(beg),
		end_epoch(end)
	{
		_init();

		// after init!
	  //  _version = version; 
		_buffsz = sz;
		_log = 0;
		_set = 0;
		_irc = 0;
		_hdr = false;

		_ss_position = 0;

		_out_smp = 0;
		_out_len = 0;
		_out_epo = LAST_TIME;

		// use malloc instead of new due to realocate function!
		// ====================================================
		_buffer = (char*)malloc((_buffsz + 1) * sizeof(char));  // due to realocate function!

		if (s) _gset(s); // READ SETTINGS IN CONSTRUCT ONLY!


		_pgm = "Anubis";
	}

	/* ----------
	 * destructor
	 */
	t_gcoder::~t_gcoder()
	{
		// use free instead of delete due to realocate function!
		if (_buffer)
		{
			free(_buffer);
			_buffer = NULL; 
		}

		if (_close_with_warning)
		{
			std::sort(_notes.begin(), _notes.end());
			for (auto it = _notes.begin(); it != _notes.end(); ++it) {
				if (_log) { _log->comment(0, "gcoder", it->str() + " .. " + base_name(_fname)); }
				else { cerr << it->status() << " " << it->str() << endl; }
			}
		}
	}


	/* ----------
	 * read settings
	 */
	int t_gcoder::_gset(t_gsetbase* s)
	{
		gtrace("t_gcoder::gset");

		if (!s) return -1;

		_set = s;
		_beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end = dynamic_cast<t_gsetgen*>(_set)->end();

		_hdr = dynamic_cast<t_gsetgen*>(_set)->thin();

		if (_gnss)
			_sys = dynamic_cast<t_gsetgen*>(_set)->sys();

		/*_rec = dynamic_cast<t_gsetgen*>(_set)->rec();*/
		_rec = dynamic_cast<t_gsetgen*>(_set)->rec_all();//to do,yjqin
		// jdhuang
		if(_rec.empty()) _rec = dynamic_cast<t_gsetgen*>(_set)->recs();
		_int = dynamic_cast<t_gsetgen*>(_set)->sampling();
		_scl = dynamic_cast<t_gsetgen*>(_set)->sampling_scalefc(); // scaling 10^decimal-digits


		for (auto itGNSS = _sys.begin(); itGNSS != _sys.end(); ++itGNSS)
		{
			string gs = *itGNSS;
			GSYS gsys = t_gsys::str2gsys(gs);

			_sat[gsys] = dynamic_cast<t_gsetgnss*>(_set)->sat(gsys, false); // empty if not set (to speed up!)
			_obs[gsys] = dynamic_cast<t_gsetgnss*>(_set)->obs(gsys, false); // empty if not set (to speed up!)
			_nav[gsys] = dynamic_cast<t_gsetgnss*>(_set)->nav(gsys, false); // empty if not set (to speed up!)

			// extend gobs list with completely defined signals
			set<string> sgobs = dynamic_cast<t_gsetgnss*>(_set)->gobs(gsys);
			for (auto it = sgobs.begin(); it != sgobs.end(); it++) _obs[gsys].insert(*it);

		}

		return 0;
	}


	// set path
	// ----------
	void t_gcoder::path(string s)
	{
		substitute(s, GFILE_PREFIX, "");

		_fname = s;
	}


	/* ----------
	 * add message
	 */
	void t_gcoder::mesg(t_note n, string s)
	{
		t_gnote m(n, "gcoder", s);

		//if( find(_notes.begin(), _notes.end(), m) == _notes.end() ){
		//     _notes.push_back(m);

		// eliminate repeating messages
		bool exist = false;
		vector<t_gnote>::const_iterator it;
		for (it = _notes.begin(); it != _notes.end(); ++it) {
			if (*it == m) { exist = true; }
		}

		if (!exist) _notes.push_back(m);
	}

	/* ----------
	 * get message
	 */
	vector<t_gnote> t_gcoder::mesg()
	{
		return _notes;
	}


	/* ----------
	 * initialize
	 */
	void t_gcoder::_init()
	{
		gtrace("t_gcoder::_init");

		_version = _initver;
		_endpos = 0;
	}


	/* ----------
	 * clear
	 */
	void t_gcoder::clear()
	{
		gtrace("t_gcoder::clear");

		// use free instead of delete due to realocate function!
		// modified by zhanglong
		if (_buffer)
		{
			free(_buffer);
			_buffer = NULL;
		}


		_init();

		// use malloc instead of new due to realocate function!
		// ====================================================
		_buffer = (char*)malloc((_buffsz + 1) * sizeof(char));  // due to realocate function!
	}


	/* ----------
	 * get single line from the buffer
	 */
	int t_gcoder::_getline(string& str, int from_pos)
	{
		gtrace("t_gcoder::_getline");

		return _decode_buffer.getline(str, from_pos);


	}



	int t_gcoder::_getbuffer(const char *& buff)
	{
		gtrace("t_gcoder::_getbuffer");

		return _decode_buffer.getbuffer(buff);
	}

	/* ----------
	 * cummulate buffer
	 */
	int t_gcoder::_add2buffer(char* buff, int sz)
	{
		gtrace("t_gcoder::_add2buffer");

		return _decode_buffer.add(buff, sz);

	}



	/* ----------
	 * remove from buffer
	 */
	int t_gcoder::_consume(int bytes_to_eat)
	{
		gtrace("t_gcoder::_consume");

		return _decode_buffer.consume(bytes_to_eat);

	}


	// add specific data structure to the coder
	// ----------
	int t_gcoder::add_data(string data_id, t_gdata* data)
	{
		gtrace("t_gcoder::add_data");

		if (data_id.empty() || data == 0) return -1;

		map<string, t_gdata*>::iterator it = _data.find(data_id);
		if (it != _data.end()) {
			if (_log) _log->comment(1, "gcoder", "warning: structure " + data_id + " already exists !");
			return -1;
		}

		_data[data_id] = data;

		// can be used to directly setup data containers in individual gcoders
		_add_data(data_id, data);

		return 1;
	}



	// remove specific data structure from the coder
	// ----------
	int t_gcoder::_fill_buffer(char* buff, int sz)
	{
		gtrace("t_gcoder::_filter_buffer");

		_ss.seekg(0, _ss.end);
		long len = (long)_ss.tellg();
		len = len - _ss_position;
		_ss.seekg(_ss_position, ios_base::beg);

		int size = (len < sz) ? len : sz;

		_ss.clear();
		_ss.read(buff, size);

		string tmp = buff;
		tmp = tmp.substr(0, size);
		size_t ifirst = 0;
		if ((ifirst = tmp.find_last_of(crlf)) != string::npos) {
			tmp = tmp.substr(0, ifirst + 1);
			size = tmp.size();
			_ss.seekg(_ss_position + size, ios_base::beg);
		}

		if (_ss.fail()) {
			cout << "HEAD: any problem ?\n";

		}
		else if (_ss.gcount() == 0) {
			
			_ss_position = 0;
			_ss.str("");
			_ss.clear();
		}
		else {
			
			_ss_position += size; // _ss.gcount();
		}

		return size;
	}


	// sampling filter for epochs (return true if the epoch fits sampling)
	// ----------
	bool t_gcoder::_filter_epoch(const t_gtime& epo)
	{
		gtrace("t_gcoder::_filter_epoch");

		if (time_sync(epo, _int, _scl, _log)) {
			//   cerr << epo.str_ymdhms("fitting: ")+dbl2str(epo.dsec()) << endl;
			return true;
		}
		return false;
	}


	// GNSS/sat filter (return true if the epoch fits gnss & sat)
	// ----------
	bool t_gcoder::_filter_gnss(const string& prn)
	{
		gtrace("t_gcoder::_filter_gnss");

		string gs = t_gsys::char2str(prn[0]); // DONT USE SYS FROM ARG WHICH MAY BE EMPTY
		GSYS gsys = t_gsys::str2gsys(gs);     // DONT USE SYS FROM ARG WHICH MAY BE EMPTY
		if ((_sys.size() == 0 || _sys.find(gs) != _sys.end()) &&
			(_sat[gsys].size() == 0 || _sat[gsys].find(prn) != _sat[gsys].end()))
		{
			return true;
		}
		return false;
	}

	bool t_gcoder::_filter_gnss_addleo(const string& prn) //add leo
	{
		gtrace("t_gcoder::_filter_gnss_addleo");

		string gs;
		GSYS gsys;
		if (prn[0] == '2' || prn[0] == '3' || prn[0] == '4')
		{
			gs = "GPS";
			gsys = GPS;
		}
		else
		{
			gs = t_gsys::char2str(prn[0]); // DONT USE SYS FROM ARG WHICH MAY BE EMPTY
			gsys = t_gsys::str2gsys(gs);     // DONT USE SYS FROM ARG WHICH MAY BE EMPTY
		}

		if ((_sys.size() == 0 || _sys.find(gs) != _sys.end()) &&
			(_sat[gsys].size() == 0 || _sat[gsys].find(prn) != _sat[gsys].end()))
		{
			return true;
		}
		return false;
	}


} // namespace
