/**
*
* @verbatim
	History
	 -1.0 ZhengHJ      2019-04-04  created
	 -1.1 glfeng       2019-09-05  add function _batch_turboedit and variable _turboedit_software
	 -1.2 Huang Jiande 2020-04-12  merge UC and IF model
	 -1.3 Huang Jiande 2020-06-21  Remove gfltprocIF and gfltprocUC
	 -1.4 Huang Jiande 2020-09-11  Format the functions
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gfltproc.h
* @brief		virtual class for process
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/


#ifndef GPROCMODEL_H
#define GPROCMODEL_H

#include "ExportHeaders/ExportLSQ.h"
#include "gdata/gdata.h"
#include "gall/gallproc.h"
#include <gall/gallprec.h>
#include "gall/gallobj.h"
#include "gall/gallobs.h"
#include "gall/gallics_nav.h"
#include "gall/gallsatinfo.h"
#include "gproc/gflt.h"
#include "gproc/glsq.h"
#include "gmodels/gmodel.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/gcombmodel.h"
#include "gmodels/gstochasticmodel.h"
#include "gproc/glsqmatrix.h"
#include "gdata/gsatdata.h"
#include "gutils/gcycleslip.h"
#include "gset/gsetamb.h"
#include "gdata/gambcon.h"
#include "gproc/goi.h"

namespace great
{
	typedef vector< tuple < string, string, string > > t_lsq_equ_info;

	/**
	* @brief virtual class for process
	*/
	class LSQ_LIBRARY_EXPORT t_gfltproc
	{
	public:
		t_gfltproc();
		t_gfltproc(t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_gfltproc();

		void add_coder(const vector<t_gcoder*>& coder); 

		virtual bool ProcessBatch(t_gallproc* data, const t_gtime&beg, const t_gtime& end);

		virtual bool GenerateProduct();

		/**
		* @brief batch processing .
		* @param[in] useIon   whether use ion file.
		*/
		void setUseIon(bool useIon) { _useIon_cmd = useIon; };
		
		/**
		* @brief batch processing .
		* @param[in] ambfix   whether fix amb.
		*/
		void setUseAmbcon(bool ambfix) { _ambfix_cmd = ambfix; };


		set<string>& sat_list() { return _sat_list; }
		set<string>& rec_list() { return _rec_list; }
		set<string>& sys_list() { return _sys_list; }

		
		double	             get_obs_intv() { return _obs_intv; };
		set<string>          get_rec_list() { return _rec_list; };
		t_glog*              get_glog()     { return _glog; };
		t_gtime              get_beg_time() { return _beg_time; };
		t_gtime              get_end_time() { return _end_time; };
		t_gtime              get_crt_time() { return _crt_time; };
		t_gallobs*           get_gall_obs() { return _gall_obs; };
		t_glsq*              get_lsq()      { return _lsq; };
		t_glsqEquationMatrix get_crt_equ()  { return _crt_equ; };

		void set_crt_time(t_gtime time) { _crt_time = time; };
		void set_end_time(t_gtime time) { _end_time = time; };
		void set_beg_time(t_gtime time) { _beg_time = time; };
		void set_gall_obs(t_gallobs* obs) { _gall_obs = obs; };

	protected:

		t_glog*     _glog = nullptr;	///< log file for process
		t_gsetbase* _gset = nullptr;	///< setting from xml file

		/// The following can get from xml file
		set<string>             _rec_list;				  ///< all sites of proc (include gournd sites and leo sats)
		set<string>             _sat_list;				  ///< all sats of proc
		set<string>             _sys_list;				  ///< sat systems of proc
		vector<string>          _rover_list;              ///< rover sites for PPPRTK
		map<string, t_gtriple>  _rec_crds;                ///< site coordinates
		map<string, t_gtriple>  _rec_stds;                ///< site coordinates std

		int     _frequency    = 0;
		double	_obs_intv     = 0.0;					  ///< sampling of proc
		double  _pwc_intv     = 0.0;                      ///< sampling of pwc
		double  _grd_pwc_intv = 0.0;                      ///< sampling of pwc for gradients, glfeng 
		t_gtime _beg_time;					            ///< beg time of proc
		t_gtime _end_time;				                ///< end time of proc
		t_gtime _crt_time;					            ///< current epoch of proc

		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;
		map<string, int> _glofrq_num;          ///< glofrq number, set GLO frequency number




		t_gsetamb::FIX_MODE  _fix_mode = t_gsetamb::FIX_MODE::NO;
		t_gsetamb::FIX_DD_MODE _fix_dd_mode;					
		t_gsetamb::UPD_MODE  _upd_mode = t_gsetamb::UPD_MODE::OSB;
		SLIPMODEL            _slip_model = SLIPMODEL::SLIPMODEL_UNDEF;
		ZTDMODEL             _ztd_model = ZTDMODEL::DEF_ZTDMODEL;

		IFCB_MODEL           _ifcb_model = IFCB_MODEL::IFCB_MODEL_UNDEF;
		IONO_ORDER           _iono_order = IONO_ORDER::NONE;
		IFB_MODEL            _ifb_model = IFB_MODEL::NONE;
		DD_MODEL             _dd_model = DD_MODEL::NONE;                        ///< dd fix mode
		GRDMODEL             _grd_model = GRDMODEL::DEF_GRDMODEL;                //glfeng
		bool                 _tropo_grad    = false;    //glfeng
	

		/// The following is public for all lsq
		std::string	_crt_rec;					///< current site of proc
		std::string	_crt_sat;					///< current sat  of proc

		t_gallobs*       _gall_obs     = nullptr;			    ///< all observ info
		t_gallobj*       _gall_obj     = nullptr;			    ///< all obj info include pcv and intial crd
		t_gallbias*      _gall_bias    = nullptr;


		t_gallproc*      _gall_proc    = nullptr;			    ///< all     the data used in proc
		t_gallics_nav*   _gall_gns     = nullptr;			    ///< all gns ics data used in proc 
		t_gerp*			 _gall_erp     = nullptr;
		t_gallnav*       _gall_nav     = nullptr;               ///< used in kin ppp
		t_gallprec*      _gall_prec    = nullptr;


		t_gerp*			 _gdata_erp    = nullptr;			    ///< all erp data for init ERP parameters

		t_gambcon*       _gdata_ambcon = nullptr;

		t_gallsatinfo*   _gdata_satinfo = nullptr;              ///< sat information for glonass fid

		t_glsq*	                    _lsq        = nullptr;		///< least square estimator
		shared_ptr<t_gcycleslip>    _slip12     = nullptr;      ///< turboedit class for IF
		shared_ptr<t_gcycleslip>    _slip13     = nullptr;	    ///< turboedit class for 13 if
		shared_ptr<t_gcycleslip>    _slip14     = nullptr;
		shared_ptr<t_gcycleslip>    _slip15     = nullptr;
		shared_ptr<t_gbasemodel> _base_model    = nullptr;
		shared_ptr<t_gbiasmodel> _bias_model	= nullptr;
		t_gprecisemodel*         _bias_pcemodel = nullptr;

		//

		// filter 
		t_SRIF _flt; //main
     


        t_gupdateparinfo _updateinfo_crt_time; // record update info for current time

		// <sats_name (leo_long_name), interval>  the interval of updating orb parameters
		map<string,int> _update_orb_interval; 
		// the interval of orbit integration
		map<string,int> _update_oi_interval;
		// the number of orbit parameters for each satellite
		map<string, shared_ptr<gnut::t_orb_mode>> _orb_mode;
		map<string, int> _orb_par_num;
		// the eclipse marker of each satellite
		map<string, bool> _eclipse_marker;
		// t_goi for orbit integration
		map<string, shared_ptr<t_goi> > _poi;
		// recor the state transition (dOrb(t) = phi(t,0) dOrb(0)) 
		// its sample should be same with update_orb_interval
		map<string, map<int, Matrix>> _phi_map;
		// record the inverse of state transition matrix 
		// its sample should be same with udpate_orb_interval
		map<string, map<int, Matrix>> _phi_i_map;
		// record the orb state results of OI (same as the t_gallorb)
		map<string, map<t_gtime, ColumnVector>> _state_map;
		map<string, map<t_gtime, ColumnVector>> _estimated_state_map;
		map<string, map<t_gtime, t_gtriple>> _estimated_state_map_trs_sp3;
		// reocrd the orbit integration reference time
		map<string, t_gtime> _orb_ref;
		// record the time of orb_state, the same as the t_gpar PXSAT
		map<string, t_gtime>  _orb_state_ref;
		// record the position of orb_state_ref relative to orb_ref
		map<string, int> _epo;

		// record the sat continous time length
		map<string, int> _sat_valid_obs_length;
		map<string, bool> _convergence;
		set<string> _last_bad_sats;
		set<string> _re_integration_sats;
        map<string,int> _normal_count;
        map<string, vector<double> > quantity_list;

		map<string, t_gtime> _rec_last_proc_time;
		map<string, double> _proc_int;


        

		t_gallorb* _orb;
		t_glsqEquationMatrix _equ4one_epoch;

		bool _is_convergence = false;

		double _get_rec_proc_int(const string& rec);

		void _time_update_orb_parameters(const t_gtime& crt_epoch,vector<t_gsatdata>& obsdata);

		void _update_srif_with_updateinfo(t_gupdateparinfo& update_info);
	
		void _get_orb_update_parameters(t_gupdateparinfo& orb_update_info,set<string> sats_have_obs);

		void _update_orb(string sat, t_gtime beg, t_gtime end, vector<t_gpar::t_type> state_init_name, vector<double> state_init_value, t_gallorb* orb);

		
		// record now process data
		t_glsqEquationMatrix _crt_equ;
		map<string, map<t_gtime, t_gtriple> > _kin_xyz;   ///< xyz coord reference at each epoch, for update parameters
		map<string, map<t_gtime, t_gtriple> > _kin_vel;   ///< xyz velocity reference at each epoch, for inital value

		 

		map<string, vector<string>>  _sat_sta;
		set<string> _mode;  ///< POD mode(S/D/K/SD)

		double             _recclk_threshold = 0.0;

		int _num_threads = 8;

		int _cmb_equ_msec;
		int _remove_par_msec;
		int _preprocess_msec;
		int _timeupdate_msec;
		int _measure_msec;	
		int _oi_msec;
		int _ambfix_msec;

	protected:
		double _maxres_norm = 0.0;                   
		t_glog* _clk_log      = nullptr;
		t_giof* _satclkfile   = nullptr;
		t_giof* _rtnetclkfile = nullptr;
		t_giof* _satclk13file = nullptr;
		bool   _carrier_range = false;                ///< use carrier range data [true=1, false=0] (add by xiongyun)
		bool   _use_carrier_range = false;            

		bool   _ultrasp3 = false;                    
		map<GSYS, double>   _maxres_C;              
		map<GSYS, double>   _maxres_L;                
		set<GSYS> sys_support;                       
		map<string, vector<string>> _outlier_preres;  
		vector<t_gcoder*>           _gcoder;         
		bool _isCreatLog = false;  ///< whether the process batch creat a log file, only when the glog is nullptr
		bool _ambfix_cmd = false;  ///< ambfix  or not, set in cmd pars
		bool _ambfix_xml = false;
		bool _useIon_cmd = false;  ///< use ion or not, set in cmd pars
		bool _write_equ  = true;

		typedef tuple<int, string, string, t_gtime, t_gtime> t_amb_index;
		map<tuple<t_amb_index, t_amb_index, t_amb_index, t_amb_index>, double> _Dd_amb_L12_used; ///< L1 ambcon while fix with L1 and L2, used in checking L1 and L3 ambcon.
		map<tuple<t_amb_index, t_amb_index, t_amb_index, t_amb_index>, double> _Dd_amb_L13_used; ///< L1 ambcon while fix with L1 and L3, used in checking L1 and L2 ambcon.

	protected:
		/// The lsq process module for override
		virtual bool _init_crs_erp_pars(t_glsq* lsq);
		virtual bool _init_rec_trop_pars(t_glsq* lsq);

		virtual bool _init_rec_crd_pars(t_glsq* lsq, bool isPPP = true);
		virtual bool _init_rec_clk_pars(t_glsq* lsq);

		virtual bool _init_sat_clk_pars(t_glsq* lsq);

		virtual bool _init_gns_ics_pars(t_glsq* lsq);
		virtual bool _init_sys_bia_pars(t_glsq* lsq);

		void _set_isb_par_setting(t_glsq* lsq, t_gpar& par, double value, double apriori);


		// -- for whole
		virtual bool _initfltProcPars(t_glsq* lsq);
		virtual bool _initfltProcData(t_gallproc* data);




		bool _select_rec_obs(const string& rec, const t_gtime& epoch, vector<t_gsatdata>& all_obs, vector<t_gsatdata>& rec_obs);
		bool _apply_carrier_range(const string& rec, const string& sat, const t_gtime& epoch, t_gsatdata& obs);

		bool _init_sat_clk();

		
		bool _init_xml_settings(t_gsetbase* set);
		bool _init_lsq_settings(t_gsetbase* set, t_gallproc* data, t_glog* log);





		string _recclk_file;	///< outputs recclk file name

		bool _myLog = false;
		bool _create_log(t_glog* log);
		bool _delete_log(t_glog* log);


		/**
		* @brief add single combine dd ambcon for raw all processs.
		*   L_NL(1) - LNL(2) - L_NL(3) + L_NL(4)
		*   L_WL(1) - LWL(2) - L_WL(3) + L_WL(4)
		* @param[in] crt_epoch   curruent epoch.
		* @param[in] lsq         lsq class.
		* @param[in] ambcon      ambcon data class.
		* @param[in] isL12       add for L1 and L2, L1 and L3.
		* @return	  true - successfully;
		*            false - failure
		*/		
		bool _pod_add_Dd_ambcon(t_SRIF* flt, t_gambcon* ambcon);
		bool _pod_add_Dd_ambcon(t_SRIF* flt, t_gambcon* ambcon, great::t_gupdateparinfo remove_info);


		//single difference
		/**
		* @brief add single combine dd ambcon for raw all processs.
		*   L_NL(1) - LNL(2)
		*   L_WL(1) - LWL(2)
		* @param[in] lsq         lsq class.
		* @param[in] Dd_amb      ambcon data class.
		* @return	  true - successfully;
		*            false - failure
		*/
		bool _pod_add_one_Sd_ambcon(t_glsq* lsq, const t_gamb_Dd* Dd_amb);


		virtual bool _processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs);

		virtual bool _processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs,t_glsqEquationMatrix& equ_result);


		//pre-fix residuals
		bool _prefix_residuals(vector<t_gsatdata>& vec_obs);
		bool _is_prefix_res;
	};


}
#endif