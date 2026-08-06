
#ifndef  GSATDATA_H
#define  GSATDATA_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: statistical function (1D)
  Version: $ Rev: $

  2014-04-18 /PV: created
  2018-09-28 /JD: revised

-*/

#include <string>
#include <vector>
#include "newmat/newmat.h"
#include "gdata/gdata.h"
#include "gdata/gobsgnss.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gall/gallnav.h"
#include "gset/gsetproc.h"

using namespace std;

namespace gnut
{

	class BASE_LIBRARY_EXPORT  t_gsatdata : public t_gobsgnss
	{

	public:
		t_gsatdata();
		t_gsatdata(string site, string sat, const t_gtime& t);
		t_gsatdata(const t_gobsgnss& obs);
		virtual ~t_gsatdata();

		void addcrd(const t_gtriple& crd);        // add satellite position

		void addpco(const t_gtriple& pco);        // add satellite pco
		void addvel(const t_gtriple& vel);        // add satellite velocity   

		void addclk(double d);                    // add satellite clocks at the transmision time
		void addreldelay(double d);				   // add satellite relative delay
		void addreccrd(const t_gtriple& crd);	   // add receiver postion

		void adddrate(double drate);

		void addsatindex(int idx);				   // add the sat index
		void addrecTime(const t_gtime& recTime);
		void addsatTime(const t_gtime& satTime);

		void	setecl(bool ecl);						// reset the eclipsing time

		int  addprd(t_gallnav* gnav, bool corrTOT = true, bool msk_health = true); // true  to support PPP/SPP (set by INP:chk_health)
		int  addprd_realtime(t_gallnav* gnav, bool corrTOT = true, bool msk_health = true); // added by zhshen
		int  addprd(const t_gtime T_sat, t_gallnav * gnav, bool corrTOT = true, bool msk_health = true); //add by glfeng
;

		void addele(double d);                // add satellite elevation 
		void addele_leo(double d);              // add satellite elevation_
		void addazi(double d);                // add satellite azimuth                             
		void addrho(double d);                // add satellite rho-vector
		void addmfH(const double& d);
		void addmfW(const double& d);
		void addmfG(const double& d);
		void addwind(const double& d);

		t_gtriple satcrd()    ;                  // get satellite position

		t_gtriple satvel()    ;                     // get satellite velocity

		double  clk();			   // get satellite clocks at the transmision time 

		double  drate();

		Matrix  orbfunct();		   // get satellite funct
		int	   satindex();		   // get sat index
		t_gtriple reccrd();		   // get reciver coord

		t_gtriple sat2reccrs();		   // get reciver coord
		Matrix  rotmat();
		Matrix  drdxpole();
		Matrix  drdypole();
		Matrix  drdut1();

		t_gtime recTime();		       // get the site receiver time

		double  ele();			       // get satellite elevation
		double  ele_leo();			   // get satellite elevation
		double  ele_deg();			   // get satellite elevation [deg]
		double  ele_leo_deg();			   // get satellite elevation [deg]
		double  azi();			   // get satellite azimuth                        
		double  rho();			   // get satellite rho-vector

		double  mfW();

		double  wind();

		double  beta();           // Sun elevation relative to orbital plane
		double  orb_angle();      // Orbit angle
		double  yaw() { return _yaw; }          // get yaw angle
		void    yaw(double yaw) { _yaw = yaw; } // set yaw angle

		vector<double> residuals(RESIDTYPE restype, GOBSTYPE type);
		void           clear_res(RESIDTYPE restype);

		void clear();
		bool valid();

		bool tb12() { return _isTb12Valid; };
		bool tb13() { return _isTb13Valid; };
		bool tb14() { return _isTb14Valid; };
		bool tb15() { return _isTb15Valid; };
		bool is_carrier_range(GOBSBAND band) {
			if (_is_carrier_range.find(band) == _is_carrier_range.end()) return false;
			else return _is_carrier_range[band];
		}

		void tb12(bool tb12) { _isTb12Valid = tb12; };
		void tb13(bool tb13) { _isTb13Valid = tb13; };
		void tb14(bool tb14) { _isTb14Valid = tb14; };
		void tb15(bool tb15) { _isTb15Valid = tb15; };
		void is_carrier_range(GOBSBAND band, bool b) { _is_carrier_range[band] = b; };

		t_gtriple		 _conf_crd;// velocity direction vector
		t_gtriple		 _e;	  // direction vector


	private:
		int  _addprd(t_gallnav* gnav, bool corrTOT = true, bool msk_health = true);         // add satellite pos, clk and ecl (corrTOT is correction of transmition time)
		int  _addprd_realtime(t_gallnav* gnav, bool corrTOT = true, bool msk_health = true);         
		int  _addprd(const t_gtime& T_sat, t_gallnav* gnav, bool corrTOT = true, bool msk_health = true);  
		int _correction(double * xyz, double * vxyz, double & clk, double & dclk, double * xyz_corr, double * vxyz_corr, double & clk_corr, double & dclk_corr);
		double  _b();
		double  _orb_angle();

		virtual void _clear();
		virtual bool _valid() const;

		t_gtriple         _satcrd; // satellite position (X,Y,Z)
		t_gtriple         _satcrdcrs; // satellite position in crs(X,Y,Z)
		t_gtriple         _satpco; // satellite pco
		t_gtriple         _satvel; // satellite velocity
		t_gtriple         _satvel_crs;
		Matrix			 _orbfunct; // satellite funct
		int			     _satindex; // satellite index in lsq
		t_gtriple		 _reccrd; // reciver positoin(TRS)
		t_gtriple		 _reccrdcrs; // reciver positoin(TRS)
		t_gtriple		 _sat2reccrs; // reciver positoin(TRS)
		Matrix			 _rotmat;
		Matrix			 _drdxpole;
		Matrix			 _drdypole;
		Matrix		     _drdut1;
		Matrix		     _scf2crs;
		Matrix		     _scf2trs;
		double			 _drate;
		t_gtriple         _dloudx;
		t_gtime			 _TR;	  // site recieve time
		t_gtime			 _TS;	  // satellite sends
		double            _clk      = 0.0;    // satellite clocks (precise, time of transmision)
		double            _dclk     = 0.0;   // satellite clocks drift (precise, time of transmision)
		double			 _reldelay  = 0.0;//satellite releative delay
		double            _ele      = 0.0;    // satellite elevation
		double            _ele_leo  = 0.0;    // satellite elevation
		double            _azi      = 0.0;    // satellite azimuth
		double            _rho      = 0.0;    // satellite-station geometrical distance
		bool              _eclipse  = 0.0; // eclipsed satellite
		double            _mfH      = 0.0;
		double            _mfW      = 0.0;
		double            _mfG      = 0.0;
		double            _wind     = 0.0;
		bool              _low_prec = 0.0; // low precision of sat pos
		bool              _slipf    = 0.0;	 // cycle slip flag   
		bool              _isTb12Valid = false; 
		bool              _isTb13Valid = false; 
		bool              _isTb14Valid = false; // add to mark is avial in tb log 14;
		bool              _isTb15Valid = false; // add to mark is avial in tb log 15;
		map<GOBSBAND, bool> _is_carrier_range; //add to mark carrier-range

		// normalized residuals
		vector<double>    _code_res_norm;
		vector<double>    _phase_res_norm;

		// original residuals
		vector<double>    _code_res_orig;
		vector<double>    _phase_res_orig;

		double _beta_val      = 0.0;
		double _orb_angle_val = 0.0;
		double _yaw           = 0.0;

	};

} // namespace

#endif
