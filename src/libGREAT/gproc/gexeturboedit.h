/**
*
* @verbatim
	History
	 -1.0 Guolong Feng  2020-03-20 created
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gexeturboedit.h
* @brief		main control class for turboedit
* @author       Guolong Feng, Wuhan University
* @date		    2020-03-20
*
*/


#ifndef GEXETURBOEDIT_H
#define GEXETURBOEDIT_H

#include "ExportHeaders/ExportBase.h"
#include "gutils/gfileconv.h"
#include "gutils/gturboedit.h"

#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"
#include "gset/gsetrec.h"
#include "gset/gsetbase.h"
#include "gset/gsetinp.h"
#include "gset/gsetturboedit.h"

#include "gio/glog.h"
#include "gio/gfile.h"

#include "gall/gallobs.h"
#include "gall/gallnav.h"
#include "gall/gallproc.h"


#include "gall/gallobj.h"
#include "gset/gsetout.h"
#include "gproc/gqualitycontrol.h"

#include "gcoders/satinfo.h"  
#include <algorithm>
#include <stdint.h>
#include <string>
#if  defined _WIN32 || defined _WIN64
#include <io.h>
#include <direct.h> 
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#include <unistd.h>
#include <sys/stat.h>
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

using namespace gnut;

namespace great
{
	/*
	* Data flag
	* -2147483648    NODATA
      -1073741824    NO4
        536870912    LOWELE
                8    MWjump
               32    gap
        268435488    short
	*/
	// ---------------------
	enum DATAFLAG {

		GOOD = 0,
		GFJUMP = 2,   // ambiguity inserted due to jump in ionosphere-lane
		MWJUMP = 3,   // ambiguity inserted due to jump in wide-lane		
		GAP = 5,   // ambiguity inserted after long gap
		LLI = 8,   // ambiguity inserted accoording to LLI in rinex-file

		PCBAD = 16,  // not used because of range is bad

		NOBRD = 22,  // not used because not broadcast ephemeris is available
		GFBAD = 24,  // not used because ionosphere - lane is bad
		MWBAD = 25,  // not used because wide - lane is bad

		ASHADOW = 27,  // not used because II / IIA after shadow
		SHRT = 28,  // not used because data piece is too short
		LOWELE = 29,  // not used because of lower elevation
		NO4 = 30,  // not used because the data set is not complete
		NODATA = 31,  // No data in rinex-file
	};


	/** @brief Specified sat data for post-processing turboedit. */
	class BASE_LIBRARY_EXPORT t_gturbo_post
	{
	public:
		/** @brief default constructor. */
		t_gturbo_post() {};
		/** @brief default destructor. */
		virtual ~t_gturbo_post() {};
		
		map<int, int>                iflag;         ///< epoch, flag_index
		map<int, double>             elev;          ///< epoch, elevation
		map<int, map<int, double>>   obs;           ///< epoch, index, obs
		int                          last_epo;      ///< last epoch
	};

	/** @brief Smooth information for real-time turboedit. */
	class BASE_LIBRARY_EXPORT t_gturbo_rt
	{
	public:
		/** @brief default constructor. */
		t_gturbo_rt() {};
		/** @brief default destructor. */
		virtual ~t_gturbo_rt() {};
		map<GOBS, bool>       mwslip;
		map<GOBS, bool>       gfslip;
		map<GOBS, t_gtime>    last_time;
		map<GOBS, int>        epo_num;                      ///< Store the number of consecutive epochs
		map<GOBS, double>     GF;                           ///< store GF value
		map<GOBS, pair<double, double> >   smooth_MW;       ///< store MW smooth value, std
		map<GOBS, vector<double>>          origin_MWs;      ///< sort according _gap
	};

	/**
	* @brief main control class for caculate clk
	*/
	class BASE_LIBRARY_EXPORT t_gexeturboedit
	{
	public:
		
		t_gexeturboedit(t_gsetbase* set, t_glog* glog, t_gallobs* gobs, t_gallnav* gnav);
		t_gexeturboedit(t_gsetbase* set, t_glog* glog, t_gallproc* gdata);
		virtual ~t_gexeturboedit();

		bool ProcessBatch(string site, const t_gtime& begT, const t_gtime& endT, double sampling, bool isMarkEndEpoch=false);


		void set_site_crd(string site, t_gtriple crd);  

	protected:		
		t_gsetbase*                    _set;              ///< set
		t_glog*                        _log;              ///< log file
		///< Global Settings	       
		int                            _epo_num;
		int                            _interval;
		t_gtime                        _beg_time;
		t_gtime                        _end_time;
		set<string>                    _sys_list;
		set<string>                    _sat_list;
		set<string>                    _rec_list;
		map<string, t_gtriple>         _rec_crds;        ///< site coordinates
		

		///< input data
		t_gallobs*       _allobs;
		t_gallnav*       _allnav;
        t_gallobj*       _allobj;                         ///< For sppflt

		///< Post-processing Settings
		map<int, map<GSYS, pair<GOBSBAND, GOBSBAND>>>        _post_all_band;    ///< set post-processing bands
		map<GSYS, pair<GOBSBAND, GOBSBAND>>                  _post_crt_band;    ///< set current post-processing bands
		map<GSYS, map<GOBSBAND, FREQ_SEQ> >                  _freq_index;
		map<GSYS, vector<GOBS>>                              _post_phase_obs;   ///< set phase observation type
		map<string, t_gturbo_post>                           _post_data;
		map<int, t_gtime>                                    _post_time;
		map<string, vector<shared_ptr<t_gcycleslip>>>        _cycleslip;
		map<string, double>                                  _scaling_factor;   ///< For check large cycle, post-processing
		map<string, int>                                     _glonass_freq_id;  ///< glonass frequency id
		

		///< Real-time processing mode
		int                                                 _smooth_window;
		int                                                 _frequency;
		map<string, map<string, map<GOBS, t_gturbo_rt>>>    _rt_data;
		t_goutliers_process                                 _outliers_proc;
		OBSCOMBIN                                           _observ;
		bool												_cycleslip_in_proc;
		t_giof*                                             _debug_turbo = nullptr;

		//< And/F9p
		map<string, map<string, map<GOBSBAND, vector<double>>>>    _epo_diff;
		map<string, map<string, map<GOBSBAND, t_gtime>>>           _epo_diff_time;

		// for npp
		vector<shared_ptr<t_gobsgnss>>                      _inputEpoData;
		vector<t_gsatdata>*                                 _satdata_npp = nullptr;
		map<string, map<string, int>>                       _map_nppdata_idx; // site sat idx

		///< currenct processing site
		t_gtime		                         _crt_time;       
		int                                  _crt_epo;
		string                               _crt_rec;                         ///< For Log file


		///< turboedit settings
		bool   _amb_output;         
		bool   _ephemeris;          ///< whether use ephemeris, default true
		bool   _simulation;         ///< whether simulation data, default False;
		bool   _site_kinematic;    
		bool   _check_pc;           ///< check PC, default 
		bool   _check_mw;           ///< check MW conbination (TURBOEDIT), default True
		bool   _check_gf;           ///< check GF conbination (TURBOEDIT), default True
		bool   _check_sf;           ///< check single frequency data, default is not.
		bool   _check_gap;          ///< check data gap, ambiguity is inserted if data missing longer
		bool   _check_short;        ///< check short arc, data piece shorter considered as short piece

        double _pc_limit;           ///< default 250m		
		double _mw_limit;           ///< widelane limit, default 4 cycles
		double _gf_limit;           ///< geometry_free_limit, default 1 cycles       ///////// in real-time/lite mode , unit is meters
		double _gf_rms_limit;       ///< default 2 cycles
		double _sf_limit;           ///< default ????
		double _min_elev;           ///< cutoff elevation, default is 0
		int    _gap_limit;          ///< epoch number, ambiguity is inserted if data missing longer, default 20 epoch
		int    _short_limit;        ///< epoch number, data piece shorter considered as short piece, default 10 epoch
		bool   _check_statistics;   ///< For Network Processing, PPP should Close; Whether consider following threshold, default not=false
		double _min_percent;        ///< The percentage of epochs(>4 sats) in all useful epochs
		int    _min_mean_nprn;      ///< Minimum mean sat numbers (sum of obs_number(>4 sats) / epoch_number(> 4 sats) )
		int    _max_mean_namb;      ///< Ambiguities Number / sat_number 
	
		///< Functions
		/**
		* @brief  output control informations
		*/
		void _outputControlInfo();

		/**
		* @brief  initial function of processing turboedit
		*/
		void _initCommonProcess();

		/**
		* @brief  initial function of post-processing turboedit
		*/
		void _initRealTimeProcess();

		/**
		* @brief  main function of real-time processing turboedit
		*/
		bool _realtimeProcess(const t_gtime& begT, const t_gtime& endT, bool isMarkEnd=false);



		void _realTimeMWGF_new(shared_ptr<t_gobsgnss> oneObs, GOBSBAND b1, GOBSBAND b2);



		/**
		* @brief  compute distance between satellite and station, satellite
		*         elevation, satellite clock bias,  using the specified brd
        *         record, or the nearest record.
		* @param[in] sat        sat name
		* @param[in] distance   distance between satellite and station, meters 
		* @param[in] elev       elevation, degree
		* @param[in] sat_clk    satellite clock bias, meters
		* @return               Ok is true; not is false
		*/
		bool _getSatInfo(string sat, double& distance, double& elev, double& sat_clk);
		
	};
	
}

#endif /*GTURBOEDIT_H*/