#include "gallnav.h"
#include "gallnav.h"

/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.

-*/

#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "gall/gallnav.h" 
#include "gutils/gsys.h"
#include "gutils/gtimesync.h"
#include "gutils/gtypeconv.h"
#include "gutils/gfileconv.h"
using namespace std;

namespace gnut
{

	// constructor
	// ----------
	t_gallnav::t_gallnav()
		: t_gdata(),
		_com(false),
		_offset(0),
		_nepoch(t_gtime::GPS),
		_multimap(false),
		_overwrite(false),
		_chk_health(false), // jdhuang : change true to false
		_chk_navig(true),
		_chk_tot(false)
	{
		gtrace("t_gallnav::constructor");
		id_type(t_gdata::ALLNAV);
		id_group(t_gdata::GRP_EPHEM);
	}


	// destructor
	// ----------
	t_gallnav::~t_gallnav()
	{
		gtrace("t_gallnav::destructor");

		_mapsat.clear();
	}


	// return gnav element
	// ----------
	shared_ptr<t_geph> t_gallnav::find(string sat, const t_gtime& t, bool chk_mask)
	{
		gtrace("t_gallnav::find");

		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, _chk_health && chk_mask);

		_gmutex.unlock(); return tmp;
	};

	

	// return position
	// ----------
	int t_gallnav::pos(string sat, const t_gtime& t, double  xyz[],
		double  var[], double  vel[], bool chk_mask) // [m]
	{
		gtrace("t_gallnav::pos");

		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, _chk_health && chk_mask);
		// added by yqyuan
		if (tmp == _null) 
		{
			tmp = t_gallnav::_find(sat, t, false);
		}
		// end adding
		if (tmp == _null)
		{
			for (int i = 0; i < 3; i++)
			{
				xyz[i] = 0.0;
				if (var) var[i] = 0.0;
				if (vel) vel[i] = 0.0;
			}
			_gmutex.unlock(); return -1;
		}

		int irc = tmp->pos(t, xyz, var, vel, _chk_health && chk_mask);
        _log->logInfo("t_gallnav","pos", "Sat: " + sat + " used nav t: " + tmp->epoch().str_mjdsod()
                        + " for time: " + t.str_mjdsod());
		// added by yqyuan
		if (irc == -1) 
		{
			irc = tmp->pos(t, xyz, var, vel, false);
		}
		// end adding		

		_gmutex.unlock(); return irc;

		//  return find( sat, t )->pos( t, xyz, var, vel, _chk_health && chk_mask );
	}




	// return position aka navigation
	// ----------
	int t_gallnav::nav(string sat, const t_gtime& t, double  xyz[],
		double  var[], double  vel[], bool chk_mask) // [m]
	{
		gtrace("t_gallnav::nav");

		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, _chk_health && chk_mask);

		if (tmp == _null)
		{
			for (int i = 0; i < 3; i++)
			{
				xyz[i] = 0.0;
				if (var) var[i] = 0.0;
				if (vel) vel[i] = 0.0;
			}
			_gmutex.unlock(); return -1;
		}
		int irc = tmp->nav(t, xyz, var, vel, _chk_health && chk_mask);

		_gmutex.unlock();
		return irc;

	}


	// return clock corrections
	// ----------
	int t_gallnav::clk(string sat, const t_gtime& t, double*  clk,
		double*  var, double* dclk, bool chk_mask) // [s]
	{
		gtrace("t_gallnav::clk");

		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, _chk_health && chk_mask);

		if (tmp == _null) {
			*clk = 0.0;
			if (var) *var = 0.0;
			if (dclk) *dclk = 0.0;
			_gmutex.unlock(); return -1;
		}

		int irc = tmp->clk(t, clk, var, dclk, _chk_health && chk_mask);
		_gmutex.unlock(); return irc;

	}

	int t_gallnav::get_pos_clk_correction(const string & sat, const t_gtime & t, int iod, double * xyz, double * vxyz, double & clk, double & dclk)
	{
		return -1;
	}

	int t_gallnav::get_iod(const string & sat, const t_gtime & t)
	{
		gtrace("t_gallnav::get_iod");

		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t);
		if (tmp == _null) { _gmutex.unlock(); return -1; }

		int iod = dynamic_pointer_cast<t_gnav>(tmp)->iod();
		if (iod < 0) { _gmutex.unlock(); return -1; }

		_gmutex.unlock(); return iod;
	}

	// return first position for satellite
	// ----------
	t_gtime t_gallnav::beg_gnav(string prn)
	{
		gtrace("t_gallnav::beg_gnav");

		_gmutex.lock();

		t_gtime tmp = LAST_TIME;
		if (!prn.empty()) {
			if (_mapsat.find(prn) != _mapsat.end()) tmp = _mapsat[prn].begin()->first;
		}
		else {
			for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT) {
				for (auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it) {
					if (_mapsat[itSAT->first].begin()->first < tmp) tmp = _mapsat[itSAT->first].begin()->first;
				}
			}
		}

		_gmutex.unlock(); return tmp;
	}


	// return last position for satellite
	// ----------
	t_gtime t_gallnav::end_gnav(string prn)
	{
		gtrace("t_gallnav::end_gnav");

		_gmutex.lock();

		t_gtime tmp = FIRST_TIME;

		if (!prn.empty()) {
			if (_mapsat[prn].size() > 0 && _mapsat.find(prn) != _mapsat.end()) {
				tmp = _mapsat[prn].rbegin()->first;
			}

		}
		else {
			for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT) {
				for (auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it) {
					if (_mapsat[itSAT->first].rbegin()->first > tmp) {
						tmp = _mapsat[itSAT->first].rbegin()->first;
					}
				}
			}
		}

		_gmutex.unlock(); return tmp;
	}


	// add navigation message
	// ----------
	int t_gallnav::add(shared_ptr<t_gnav> nav)
	{
		gtrace("t_gallnav::add");

		_gmutex.lock();

		t_gtime ep(nav->epoch());
		string sat = nav->sat();

		if (_chk_navig) {
			set<string> msg; nav->chk(msg);
			if (_log) {
				for (auto it = msg.begin(); it != msg.end(); it++)
					_log->comment(2, "gallnav", *it);
			}
		}

		// test navigation type
		bool add = false;
		if (_multimap) { add = true; }    // multimap
		else if (_mapsat[sat].find(ep) == _mapsat[sat].end()) { add = true; }    // non-existent
		else if (nav->id_type() == t_gdata::EPHGAL) {                  // check nav-type
			auto itB = _mapsat[sat].lower_bound(ep);
			auto itE = _mapsat[sat].upper_bound(ep);
			while (itB != itE) {
				if (dynamic_pointer_cast<t_gnav>(itB->second)->gnavtype() == nav->gnavtype()) {
					add = false; break;
				} // exclude the message and skip!
				else { add = true; } // ok
				++itB;
			}
		}
		if (!nav->valid()) add = false;    // validity test  

		// Remove items with too close time. added by zhj
		if (add && !_multimap && !_mapsat[sat].empty()) {
			const double close_elipse = 180;

			// Get Closest time
			t_gtime closest_time;
			auto itE = _mapsat[sat].upper_bound(ep);
			if (itE == _mapsat[sat].end()) {
				closest_time = _mapsat[sat].rbegin()->first;
			}
			else if (itE == _mapsat[sat].begin()) {
				closest_time = _mapsat[sat].begin()->first;
			}
			else {
				auto itE_pre = itE;
				itE_pre--;
				if (abs(itE_pre->first - ep) > abs(itE->first - ep))
				{
					closest_time = itE->first;
				}
				else {
					closest_time = itE_pre->first;
				}
			}

			double time_diff = abs(closest_time - ep);
			if (time_diff < close_elipse)
			{
				if (ep < closest_time)
				{
					add = false;
				}
				else {
					_mapsat[sat].erase(closest_time);
				}
			}
		}

		if (add) {
			if (_log && _log->verb() >= 3) {
				ostringstream lg;
				lg << "add sat [" << nav->str_type() << "]: " << sat << " " << ep.str_ymdhms()
					<< " iod: " << fixed << setw(3) << nav->iod()
					<< " flg: " << fixed << setw(3) << nav->healthy();
				if (nav->gio()) lg << " " << base_name(nav->gio()->path());
				_log->comment(3, "gallnav", lg.str());
			}
			_mapsat[sat].insert(make_pair(ep, nav));
			if (sat.substr(0, 1) == "R") {
				_glo_freq_num[sat] = nav->freq_num();
			}
			

		}
		else if (_log && _log->verb() >= 3) {
			ostringstream lg;
			lg << "skip sat [" << nav->str_type() << "]: " << sat << " " << ep.str_ymdhms()
				<< " iod: " << fixed << setw(3) << nav->iod()
				<< " flg: " << fixed << setw(3) << nav->healthy();
			if (nav->gio()) lg << " " << base_name(nav->gio()->path());
			_log->comment(3, "gallnav", lg.str());
		}

		_gmutex.unlock(); return 0;
	}


	// list of nav messages
	// ----------
	vector<shared_ptr<t_geph> > t_gallnav::vec_nav(string prn)
	{
		gtrace("t_gallnav::vec_nav");

		_gmutex.lock();

		vector<shared_ptr<t_geph> > v_nav;

		auto itFIRST = _mapsat.begin();
		auto itLAST = _mapsat.end();

		if (!prn.empty()) itLAST++ = itFIRST = _mapsat.find(prn);

		for (auto itSAT = itFIRST; itSAT != itLAST; ++itSAT) {
			for (auto itEPO = itSAT->second.begin(); itEPO != itSAT->second.end(); ++itEPO) {
				v_nav.push_back(itEPO->second);
			}
		}

		_gmutex.unlock(); return v_nav;
	}


	// clean function
	// ----------
	void t_gallnav::clean_outer(const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gallnav::clean_outer");

		if (end < beg) return;
		if (beg == FIRST_TIME) return;
		if (end == LAST_TIME) return;

		_gmutex.lock();



		// loop over all satellites
		auto itPRN = _mapsat.begin();
		while (itPRN != _mapsat.end()) {
			string prn = itPRN->first;

			// find and CLEAN all data (epochs) out of the specified period !
			auto itFirst = _mapsat[prn].begin();
			auto itLast = _mapsat[prn].end();
			auto itBeg = _mapsat[prn].lower_bound(beg);  // greater only   old: // greater|equal
			auto itEnd = _mapsat[prn].upper_bound(end);  // greater only!

	
			// remove before BEGIN request
			if (itBeg != itFirst) {
				if (_log && _log->verb() >= 2) {
					_log->comment(2, "gallnav", "Del NAV before: " + itFirst->first.str_ymdhms() + " " + prn);
				}
				auto it = itFirst;

				while (it != itBeg && it != itLast) {
					if ((it->second).use_count() == 1) { _mapsat[prn].erase(it++); }
					else { it++; }
				}
			}

			// remove after END request
			if (itEnd != itLast) {
				if (_log && _log->verb() >= 2) {
					_log->comment(2, "gallnav", "Del NAV after : " + itEnd->first.str_ymdhms() + " " + prn);
				}
				auto it = itEnd;

				while (it != itLast) {
					if ((it->second).use_count() == 1) { _mapsat[prn].erase(it++); }
					else { it++; }
				}
			}
			itPRN++;
		}

		_gmutex.unlock();
	}


	// get number of satellites
	// ----------
	int t_gallnav::nsat(GSYS gs)
	{
		gtrace("t_gallnav::nsat(GSYS)");

		int nsatell = 0;
		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string sat = itSAT->first;
			if (t_gsys::char2gsys(sat[0]) == gs || gs == GNS) nsatell++;
		}
		return nsatell;
	}


	// get interval between messages   
	// ----------
	int t_gallnav::intv(GSYS gs)
	{
		gtrace("t_gallnav::intv(GSYS)");

		int interval = 0;  if (gs == GNS) return interval;

		// collect sampling intervals
		map<int, int> m_intv;
		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string sat = itSAT->first; if (t_gsys::char2gsys(sat[0]) != gs) continue;

			t_gtime lstEPO;
			for (auto itEPO = _mapsat[sat].begin(); itEPO != _mapsat[sat].end(); ++itEPO)
			{
				if (itEPO != _mapsat[sat].begin()) {
					int diff = int(floor(itEPO->first - lstEPO));
					m_intv[diff]++;
				}
				lstEPO = itEPO->first;
			}
		}

		// auto-detect sampling interval
		int count = 0;
		for (auto it = m_intv.begin(); it != m_intv.end(); ++it) {

			if (it->second > count) {
				interval = it->first;
				count = it->second;
			}
		}

		return interval;
	}


	// get existing number of messages
	// ----------
	int t_gallnav::have(GSYS gs, const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gallnav::have(GSYS,beg,end)");

		int existing = 0;
		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string sat = itSAT->first; if (t_gsys::char2gsys(sat[0]) != gs) continue;

			for (auto itEPO = _mapsat[sat].begin(); itEPO != _mapsat[sat].end(); ++itEPO) {
				if (itEPO->first < beg - 900 || itEPO->first > end + 900) continue;
				else existing++;
			}
		}
		return existing;
	}


	// get expected number of messages
	// ----------
	int t_gallnav::expt(GSYS gs, const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gallnav::expt(GSYS,beg,end)");

		if (end < beg) return 0;
		int xint = intv(gs);
		int xsat = nsat(gs);
		int diff = int(floor(end - beg));
		if (xint == 0 || xsat == 0 || diff == 0.0) return 0;
		return int(floor(diff / xint * xsat));
	}


	// get excluded number of messages
	// ----------
	int t_gallnav::excl(GSYS gs, const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gallnav::excl(GSYS,beg,end)");
		int exclude = 0;
		return exclude;
	}



	// get IONO CORR
    // -----------------
	t_iono_corr t_gallnav::get_iono_corr(const IONO_CORR c) const
	{
		if (_brdc_iono_cor.find(c) != _brdc_iono_cor.end()) {
			return _brdc_iono_cor.at(c);
		}
		t_iono_corr io;
		// default settings , according rtklib->rtkcmn.c->ionmodel
		if (c == IO_GPSA)
		{
			io.x0 = 0.1118e-07;
			io.x1 =-0.7451e-08;
			io.x2 =-0.5961e-07;
			io.x3 = 0.1192e-06;
		}
		else if (c == IO_GPSB)
		{
			io.x0 = 0.1167e+06;
			io.x1 =-0.2294e+06;
			io.x2 =-0.1311e+06;
			io.x3 = 0.1049e+07;
		}
		
		return io;
	}

	// set IONO CORR
	// -----------------
	void t_gallnav::add_iono_corr(const IONO_CORR c, t_iono_corr io)
	{
		if (_brdc_iono_cor.find(c) == _brdc_iono_cor.end()) {
			_brdc_iono_cor[c] = io;
		}
	}

	// return list of available satellites
	// ----------
	shared_ptr<t_geph> t_gallnav::_find(string sat, const t_gtime& t, bool chk_mask)
	{
		gtrace("t_gallnav::_find sat/time");

		if (_mapsat.find(sat) == _mapsat.end()) return _null; // make_shared<t_geph>();

		GSYS gnss = t_gsys::str2gsys(sat.substr(0, 1));

		if (_mapsat[sat].size() == 0) {
			if (_log) _log->comment(2, "gallnav", sat + " no gnav elements found");
			return _null;
		}

		auto it = _mapsat[sat].lower_bound(t);  // greater|equal  (can be still end())
		if (it == _mapsat[sat].end()) it--;                   // size() > 0 already checked above

		double maxdiff = t_gnav::nav_validity(gnss)*1.1;

		if (gnss == GAL) {
			// Galileo ephemerides are valid for toc+10min -> toc+20...180min
			for (int bck = 1; bck <= 5; bck++) {     // max eph for going back is 5
				t_gtime toc = it->second->epoch();
				if ((t < toc + 10 * 60 || t > toc + maxdiff) ||
					(_chk_tot && !it->second->chktot(t))) {        // check ToT for using past messages only
					if (_mapsat[sat].size() > 0 && it != _mapsat[sat].begin()) it--; // one more step back
					else break;
				}
				if ((t >= toc + 10 * 60 && t <= toc + maxdiff) &&
					(_chk_tot && it->second->chktot(t))) break;
			}
		}
		else if (gnss == BDS) {
			// BeiDou ephemerides are valid for toc -> toc+60min
			for (int bck = 1; bck <= 5; bck++) {     // max eph for going back is 5
				t_gtime toc = it->second->epoch();
				if ((t < toc || t > toc + maxdiff) ||
					(_chk_tot && !it->second->chktot(t))) {        // check ToT for using past messages only
					if (_mapsat[sat].size() > 0 && it != _mapsat[sat].begin()) it--; // one more step back
					else break;
				}
				if ((t >= toc && t <= toc + maxdiff) &&
					(_chk_tot && it->second->chktot(t))) break;
			}
		}
		else if (gnss == GLO) {
			t_gtime tt = t - maxdiff;  // time span of Glonass ephemerides is 15 min

			auto itt = _mapsat[sat].lower_bound(tt);  // greater|equal
			if (itt == _mapsat[sat].end()) {                         // size() > 0 already checked above
				if (_log) _log->comment(2, "gallnav", sat + " gnav element not found: " + t.str_ymdhms());
				return _null;
			}

			double dt1 = itt->first - t;
			double dt2 = it->first - t;
			if (fabs(dt1) < fabs(dt2)) it = itt;
		}
		else {
			auto itLAST = _mapsat[sat].end(); itLAST--;
			auto itSAVE = it;
			bool found = false;
			while (_chk_tot && !itSAVE->second->chktot(t) && itSAVE->second->epoch().sod() % 3600 != 0 && itSAVE != itLAST) {
				itSAVE++; // one more step forward (special case for non-regular ephemerides (e.g. XX:59:44))
			}
			if (_chk_tot && itSAVE->second->chktot(t)) {   // check found non-regular ephemeris
				found = true;
				it = itSAVE;
			}

			// 2 conditions for all cases when NAV is in fugure and the iterator should be moved back (it--)
			if (fabs(t - it->second->epoch()) > maxdiff ||  // too far navigation message in future!
				(_chk_tot && !it->second->chktot(t) && !found)        // check ToT for using past messages only
				) {
				if (_mapsat[sat].size() > 0 && it != _mapsat[sat].begin()) it--; // one more step back
			}

		}

		// tested found ephemeris
		if (fabs(t - it->second->epoch()) > maxdiff ||
			(_chk_tot && !it->second->chktot(t))) {   // simulation real-time epoch search 
			if (_log) {
				string lg(sat + " gnav element not found: " + t.str_ymdhms());
				_log->comment(2, "gallnav", lg);
			}
			return _null;
		}

		if (_chk_health && chk_mask) {
			if (it->second->healthy()) return it->second;
			else                      return _null;
		}
		return it->second;

	}

    shared_ptr<t_geph> t_gallnav::find_previous(string sat, const t_gtime &t1) {
        gtrace("t_gallnav::find_preivous");
        lock_guard<t_gmutex> lck_lock(_gmutex);

        auto& sat_items = this->_mapsat[sat];

        auto iter1 = sat_items.upper_bound(t1);
        if (!sat_items.empty() && iter1 != sat_items.begin()){
            iter1--;
            return iter1->second;
        }

        return nullptr;
    }

} // namespace
