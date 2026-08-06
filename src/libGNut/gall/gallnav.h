
#ifndef GALLNAV_H
#define GALLNAV_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose:
  Version: $ Rev: $

  2011-02-14 /JD: created
  2018-08-13 /JD: updated

-*/

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>

#include "gdata/gdata.h"
#include "gdata/geph.h"
#include "gdata/gnavgps.h"
#include "gdata/gnavglo.h"
#include "gdata/gnavgal.h"
#include "gdata/gnavbds.h"
#include "gdata/grxnhdr.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtetrad.h"

#define MAX_GPS_PRN 32
#define MAX_GLO_PRN 24
#define MAX_GAL_PRN 30
#define NAV_BUF   1024


using namespace std;

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_gallnav : public t_gdata
	{

	public:
		t_gallnav();
		virtual ~t_gallnav();

		typedef multimap<t_gtime, shared_ptr<t_geph> > t_map_ref; // all data for a single satellite
		typedef map<string, t_map_ref>                 t_map_sat; // all data for all satellites

		virtual int nav(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true); // [m]
		virtual int pos(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true); // [m]
		virtual int clk(string sat, const t_gtime& t, double* clk, double* var = NULL, double* dclk = NULL, bool chk_mask = true); // [s]  
		virtual int get_pos_clk_correction(const string& sat, const t_gtime& t, int iod, double *xyz, double *vxyz, double &clk, double &dclk);
		virtual int get_iod(const string& sat, const t_gtime& t);  
		virtual void clean_outer(const t_gtime& beg = FIRST_TIME,
			const t_gtime& end = LAST_TIME);

		virtual t_gtime beg_gnav(string prn = "");                            // get first t_gnav epoch
		virtual t_gtime end_gnav(string prn = "");                            // get last  t_gnav epoch
		// IMPROVE beg_time/end_time to distinguish GALLNAV/GALLPREC - t_satview !
		virtual t_gtime beg_time(string prn = "") { return beg_gnav(prn); }    // get first t_gnav epoch
		virtual t_gtime end_time(string prn = "") { return end_gnav(prn); }    // get first t_gnav epoch
		// =========================================================


		virtual int add(shared_ptr<t_gnav> nav);                         // add single navigation message

		virtual vector<shared_ptr<t_geph> > vec_nav(string prn = "");     // get list of nav messages



		virtual void multi(bool b) { _multimap = b; }                       // set/get multimap mode


		virtual void overwrite(bool b) { _overwrite = b; }                  // set/get overwrite mode
                 

		virtual void com(bool b) { _com = b; }                              // position/clock reference point
                          

		virtual void offset(int i) { _offset = i; }                       // set/get offset for satdata


		virtual int nsat(GSYS gs);                                         // get number of satellites
		virtual int intv(GSYS gs);                                         // get interval between messages
		virtual int have(GSYS gs, const t_gtime& beg, const t_gtime& end); // get existing number of messages
		virtual int expt(GSYS gs, const t_gtime& beg, const t_gtime& end); // get expected number of messages
		virtual int excl(GSYS gs, const t_gtime& beg, const t_gtime& end); // get excluded number of messages

	

		virtual shared_ptr<t_geph> find(string sat, const t_gtime& t, bool chk_mask = true);   // find appropriate t_geph element (interface only)
	

        virtual shared_ptr<t_geph> find_previous(string sat, const t_gtime& t1);   // find eph time <= t1
		

		map<string, int> glo_freq_num() { return _glo_freq_num; }
		t_iono_corr get_iono_corr(const IONO_CORR c) const;
		void add_iono_corr(const IONO_CORR c, t_iono_corr io);

	protected:
		virtual shared_ptr<t_geph> _find(string sat, const t_gtime& t, bool chk_mask = true);  // find appropriate t_geph element


		bool               _com;         // position/clock reference point (com = true; apc = false);
		int                _offset;      // offset for RTCM corrections
		int                _nepoch;      // maximum number of epochs (0 = keep all)
		t_map_sat          _mapsat;      // map over all satellites (positions,..?)
		bool               _multimap;    // use multimap for redundant records
		bool               _overwrite;   // overwrite mode (for derived classes with MAP)
		bool               _chk_health;  // check satellite health (navigation)
		bool               _chk_navig;   // check navigation messages (internal)
		shared_ptr<t_geph> _null;        // null pointer

		bool _chk_tot;

		map<string, int> _glo_freq_num;
		t_map_iono       _brdc_iono_cor;

	};

} // namespace

#endif
