
#ifndef GALLPREC_H
#define GALLPREC_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: implementation of SP3 ephemeris container
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include "gall/gallnav.h"
#include "gdata/gephprec.h"
#include "gutils/gtriple.h"
#include "gmodels/gpoly.h"


using namespace std;

namespace gnut
{

	typedef map<string, double>       t_map_dat;     // single data for a single satellite
	typedef map<t_gtime, t_map_dat>       t_map_epo;     // all data for a single satellite
	typedef map<string, t_map_epo>       t_map_prn;     // all data for all satellites

	class BASE_LIBRARY_EXPORT t_gallprec : public t_gallnav
	{
	public:
		enum  clk_type {
			AS,
			AR,
			UNDEF
		};


	public:

		t_gallprec();
		virtual ~t_gallprec();

		typedef map<string, shared_ptr<t_gephprec> > t_map_sp3;



		virtual int pos(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true); // [m]
		virtual int clk_cdr(string sat, const t_gtime& t, double* clk, double* var = NULL, double* dclk = NULL, double* ifcb = NULL, bool chk_mask = true);//add ifcb by xiongyun

		int intv(string sat);
        int intv();

		void add_interval(string sat, int intv);
        void add_interval(int intv);
		void add_agency(string producer);

		int addvel(string sat, const t_gtime& ep, double xyzt[4], double var[4]);  // [m], [m]
		int addclk(string sat, const t_gtime& ep, double  clk[3], double var[3]);  // [s]
		int addclk_tri(string sat, const t_gtime& ep, double  clk[3], double var[3]);  // [s] add by xiongyun

	 	int addpos(string sat, const t_gtime& ep, t_gtriple xyz, double t, t_gtriple dxyz, double dt);  // [m], [usec]

		void use_clkrnx(bool b) { _clkrnx = b; }
		void use_clksp3(bool b) { _clksp3 = b; }
		void use_clknav(bool b) { _clknav = b; }
		void use_posnav(bool b) { _posnav = b; }
		void clean_all() { _mapsp3.clear(); _mapprec.clear(); }
		t_gtime get_beg();
		t_gtime get_end();
		string get_agency();
		vector<string> get_sats();
		vector<string> get_sat3();
		string get_data_type();
		string get_sat_type();
		int get_sp3_size();
		void set_beg(t_gtime beg);
		void set_end(t_gtime end);
		void set_sat(vector<string>sat);
		void set_sat3(vector<string> sat);
		void set_data_type(string type);
		void set_agency(string agency);
		void set_sat_type(string type);

		// // =========================================================
		virtual t_gtime beg_clk(string prn = "");                            // get first precise clocks epoch
		virtual t_gtime end_clk(string prn = "");                            // get last precise clocks epoch

        	








	protected:
		virtual shared_ptr<t_geph> _find(string sat, const t_gtime& t);  // find appropriate t_geph element
		virtual int         _get_crddata(string sat, const t_gtime& t);  // fill PT,X,Y,Z vectors
		virtual int         _get_clkdata(string sat, const t_gtime& t);  // fill CT,C vectors


		// // also gallnav member can be used here !!!
		t_map_sat         _mapprec;     // map of sp3 polynomials
		t_map_prn         _mapsp3;      // precise orbits&clocks (SP3) - full discrete data sets
		t_map_prn         _mapclk;      // precise clocks (CLOCK-RINEX) - full discrete data sets

	private:
		t_map_sp3         _prec;        // CACHE: single SP3 precise ephemeris for all satellites
		unsigned int      _degree_sp3;  // polynom degree for satellite sp3 position and clocks
		double            _sec;         // default polynomial units
		t_gtime           _ref;         // selected reference epoch for crd data/polynomials
		t_gtime           _clkref;      // selected reference epoch for clk data/polynomials
		bool              _clkrnx;      // true: use               clk from Rinex Clocks
		bool              _clksp3;      // true: use alternatively clk from sp3 (~15min!)
		bool              _clknav;      // true: use alternatively nav (low-precise clocks)
		bool              _posnav;      // true: use alternatively nav (low-precise orbits)
	
		unsigned int	  _udclkInt;	// SSR clk interval (added by zhShen)
		unsigned int	  _udorbInt;    // SSR orb interval (added by zhShen)
		int               _intv;
        map<string,int>    _intvm;
		string            _agency;
		t_gtime           _tbeg;
		t_gtime           _tend;
		vector<string>    _sats;
		vector<string>    _sat3;
		string            _datatype;
		string            _sattype;
		string            _frame;
	  // CACHE for approximative solutions
		map<string, t_gtime>   _poly_beg;
		map<string, t_gtime>   _poly_end;
		map<string, t_gpoly>   _poly_x;
		map<string, t_gpoly>   _poly_y;
		map<string, t_gpoly>   _poly_z;
		
		// BEGIN OF TEMPORARY (ALTERNATIVE for direct interpolation)
		vector<double>    _PT;          // vector of time-difference (X for polynomials)
		vector<t_gtime>   _T;           // vector of full time       (X for polynomials)
		vector<double>    _X;           // vector of x-coordinate    (Y for polynomials)
		vector<double>    _Y;           // vector of y-coordinate    (Y for polynomials)
		vector<double>    _Z;           // vector of z-coordinate    (Y for polynomials)

		vector<double>    _CT;          // vector of time-difference (X for polynomials)
		vector<double>    _C;           // vector of clk correction  (Y for polynomials)
		vector<double>    _IFCB_F3;     // vector of ifcb_f3 (add by xiongyun)



		vector<double>    _PTCorr;          // vector of time-difference for SSR Correction(X for polynomials)
		vector<t_gtime>   _TCorr;           // vector of full time for SSR Correction      (X for polynomials)
		vector<double>    _XCorr;           // vector of x-coordinate for SSR Correction   (Y for polynomials)
		vector<double>    _YCorr;           // vector of y-coordinate for SSR Correction   (Y for polynomials)
		vector<double>    _ZCorr;           // vector of z-coordinate for SSR Correction   (Y for polynomials)

		vector<double>    _CTCorr;          // vector of time-difference for SSR Correction(X for polynomials)
		vector<double>    _CCorr;           // vector of clk correction for SSR Correction (Y for polynomials)
	  // END OF TEMPORARY (ALTERNATIVE)
		// CLK TYPE
		set<clk_type> _clk_type_list;
	};

} // namespace

#endif
