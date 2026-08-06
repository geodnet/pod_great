/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
	-1.1 jdhuang  2020-11-14  add detailed notes
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallamb.h
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#ifndef GALLAMB_H
#define GALLAMB_H 

#include <tuple>
#include "ExportHeaders/ExportAMB.h"
#include "gio/glog.h"
#include "gutils/gnss.h"
#include "gmodels/gpar.h"
#include "gambfix/gambOw.h"
#include "gall/gallpar.h"
#include "gall/gallobs.h"
#include "gall/gallbias.h"
#include "gutils/gmutex.h"
#include "gset/gsetamb.h"
#include "gmodels/gprecisemodel.h"

using namespace std;
using namespace gnut;

namespace great
{
	class t_pair_str_gtime_hash
	{
	public:
		size_t operator()(const pair<string, t_gtime>& a) const
		{
			return std::hash<string>()(a.first) ^ std::hash<int>()(a.second.mjd());
		}
	};

	class MW_OF_AMB
	{
	public:


		void add_data(const double& new_num);
		void clean_mean_sigma(double& mean, double& sigma);
		int size() {
			return this->arr.size() ;
		}

		t_gtime beg_time = FIRST_TIME;

	private:
		double _mean = 0.0;
		double _variance = 0.0;
		int length = 0;

		vector<double> arr;
		const static int max_arr_size = 30;
		int _next_arr_loc = 0;

		bool recompute_flag=false;

		void _recompute_mean_sigma();
	};

	class AMB_LIBRARY_EXPORT t_gallambs
	{
		
	public:

		typedef map<string, set<string> >  t_rec_sats;                                 
		typedef pair<gnut::GSYS, t_gamb::AMB_ID> t_dd_sys_index;                        
		typedef pair<string, string>       t_dd_bl_index;                           
		typedef map<t_dd_sys_index, map<string, vector<t_gamb_Sd> > > t_d;             
		typedef map<t_dd_sys_index, map<t_dd_bl_index, vector<t_gamb_Dd> > > t_dd;     
		enum class AMB_COMB_OPTION 
		{ 
			CMB_REC_SD, 
			CMB_SAT_SD, 
			CMB_REC_DD, 
			CMB_SAT_DD, 
			UNDEF
		};

		//< construct
		t_gallambs();

		//< construct
		//== input:
		//== gset : xml settings ptr
		//== log  : app log ptr
		t_gallambs(t_gsetbase* gset, t_glog* log);

		//< desconstruct
		virtual ~t_gallambs();

		//< add or replace xml setting ptr
		//== input:
		//== gset : xml setting ptr 
		void add_set(t_gsetbase* gset);

		//< add or replace app log ptr
		//== input:
		//== log  : app log file ptr
		void add_log(t_glog* log);

	

		//< add rinexo obs, amb pars, dcb file, sinex files, accroding to need time arc
		//== input:
		//== pars   : amb    data
		//== obs    : rinexo data
		//== allbias: dcb or osb data
		//== obj    : sinex data
		void add_data(t_gallpar* pars, t_gallobs* obs, t_gallbias* allbias, t_gallobj* obj, t_gtime beg_time, t_gtime end_time);
		void add_site_crd(t_gallobj* obj, t_gtime beg_time);
		void update_data(t_gallpar* pars, t_gallobs* obs, t_gallbias* allbias, t_gallobj* obj, t_gtime beg_time, t_gtime end_time);


		void clear_ambs_data();


		//< add WL judgment
		//== input:
		//== WL   : WL Judge
		void add_WL_bdeci(t_gbdeci* WL);

		//< add NL judgment
		//== input:
		//== NL   : NL Judge
		void add_NL_bdeci(t_gbdeci* NL);


		void add_glo_fid(const map<string, int>& glo_fid) { _glofrq_num = glo_fid; }

	

		//< get all Dd ambs
		//< return: Dd ambs
		t_dd get_Dds();

		//< get all amb IDs
		//< return: all amb IDs
		set<t_gamb::AMB_ID> IDs();


		//< combine Sd/Dd or other amb
		bool cmb_amb(const AMB_COMB_OPTION& opt, const t_gamb::AMB_ID& ID, const vector<pair<string, string>>& baselines, double shortest_arc);

	
	public:
		map < pair<GSYS, t_gamb::AMB_ID>, int> all_WL_num; ///< all WL ambs number
		map < pair<GSYS, t_gamb::AMB_ID>, int> fix_WL_num; ///< fix WL ambs number
		map < pair<GSYS, t_gamb::AMB_ID>, int> all_NL_num; ///< all NL ambs number
		map < pair<GSYS, t_gamb::AMB_ID>, int> fix_NL_num; ///< fix NL ambs number
	protected:
		//< combine sat baseline Sd ambs
		bool _cmb_intersat_Sd_amb(const pair<string, string>& sat_baseline, const  t_gamb::AMB_ID& ID,  double shortest_arc,map<string,t_gmutex>& add_mtx_table);
		bool _cmb_intersat_Sd_amb(const pair<string, string>& sat_baseline, const  t_gamb::AMB_ID& ID,  double shortest_arc);
		
		//< combine rec baseline Dd, sat baseline Sd, ambs
		bool _cmb_interrec_Dd_amb(const pair<string, string>& rec_baseline, const  t_gamb::AMB_ID& ID,  double shortest_arc);

		//< get one way ambs table from t_gallpar ptr;
		void _get_ows_table_from_pars();

		//< get rinexo obs used in get MW combine value for ow ambs accroding to time of xml
		void _get_obs_table_from_gobs();

		//< get rinexo obs used in get MW combine value for ow ambs
		void _get_obs_table_from_gobs(t_gtime beg,t_gtime end);
		void _get_obs_table_from_gobs_one_site(string rec,bool isRawAll,const map<t_dd_sys_index,vector<t_gamb_Ow>*>& ambs_for_sys,const set<t_gtime> local_time_table,const map<t_gtime,vector<t_gobsgnss*> >& rec_obs,map<t_dd_sys_index,map<string,MW_OF_AMB*>>& result);

		//< get dependent IF Dd amb baseline
		void _get_dependent_IF_Dd(const t_gamb::AMB_ID& ID);


		//< get the fix rate for < ID > amb 
		void _get_fix_rate(const t_gamb::AMB_ID& ID);



		//< get dependent IF Dd amb baseline, math strategy
		void math_dependent_IF_Dd(vector<t_gamb_Dd>& input, int num_Ow_1st, int num_Ow_2nd, int num_total, vector<t_gamb_Dd>& result, vector<double>& vector_space,vector<double>& crt_vector,int& ndef);
		void math_dependent_IF_Dd(vector<t_gamb_Dd>& input, int num_Ow_1st, int num_Ow_2nd,vector<t_gamb_Dd>& result, Matrix& vector_space);

		//
		//< select <ID> Dds value, fixed!
		void _select_Dds_value(const t_gamb::AMB_ID& ID);



		//< get Dd amb value
		void _get_Dds_value(const t_gamb::AMB_ID& ID);


		//< get Nl amb value
		bool _get_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);
		
		//< get IF NL amb vlaue
		bool _get_IF_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);

		
		//< get MW value 
		bool _get_MW_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);


		//< get IF Dd amb value
		bool _get_IF_Dd_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);



		//< get MW value
		void _get_MW(const string& rec, const string& sat,const int& freq, const set<t_gtime>& time_table, set<t_gtime>::iterator beg, set<t_gtime>::iterator end,
			double& value, double& sigma);
		set< pair<GSYS, t_gamb::AMB_ID>> get_Dds_index_key(const t_gamb::AMB_ID& ID);
	private:
		t_glog*         _glog  = nullptr;
		t_gsetbase*		_gset  = nullptr;
		t_gallobs*      _gobs  = nullptr;
		t_gallbias*     _gbias = nullptr;
		t_gallpar*      _gpars = nullptr;
		t_gallobj*      _gobj  = nullptr;
		double*         _pdE   = nullptr;
		double*         _pdC   = nullptr;

		t_gbdeci*       _evaluator_MW = nullptr;
		t_gbdeci*       _evaluator_NL = nullptr;
		t_gprecisemodel* _gprecise = nullptr;


		set<t_gamb::AMB_ID>                  _IDs;
		DD_MODEL                             _dd_mode;
		t_gsetamb::FIX_DD_MODE 			 _fix_dd_mode;
		DCB_MODEL                            _dcb_model;
		//FREQ_MODEL                           _freq_model;
		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;
		map<string, int> _glofrq_num;          ///< glofrq number, set GLO frequency number
		int _frequency;

		map<string, string> _site_rec;   ///< store receiver type for each station
		map<string, t_gtriple> _rec_crd;

		set<string>     _rec_list;
		set<string>     _sat_list;
		string          _class_id = "t_gallambs";	
		t_gmutex        _gmutex;

		set<t_gtime>     time_table;
		map<GSYS, map<string, map<pair<string, t_gtime>, map<int, double>>>>    _MW_used;  // sys, rec, sat, time
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, map<string,MW_OF_AMB> > >   _MW_mean_sigma; // sys, rec,sat 
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, vector<t_gamb_Ow>> >                                _nOw_ambs; // map< pair<sys,obs_type>, map< rec, vector<ow>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, vector<t_gamb_Sd>> >                                _nSd_ambs; // map< pair<sys,obs_type>, map< rec, vector<sd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, vector<t_gamb_Sd>> >                                _nSd_ambs_fix_WL; // map< pair<sys,obs_type>, map< rec, vector<sd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, vector<t_gamb_Sd>> >                                _nSd_ambs_fix_NO; // map< pair<sys,obs_type>, map< rec, vector<sd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<pair<string, string>, vector<t_gamb_Dd>> >                  _nDd_ambs; // map< pair<sys,obs_type>, map< rec1,rec2,sat1,sat2, vector<Dd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<pair<string, string>, vector<t_gamb_Dd>> >                  _nDd_ambs_fix_WL; // map< pair<sys,obs_type>, map< rec1,rec2,sat1,sat2, vector<Dd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<pair<string, string>, vector<t_gamb_Dd>> >                  _nDd_ambs_fix_NO; // map< pair<sys,obs_type>, map< rec1,rec2,sat1,sat2, vector<Dd>> >

	};

	bool sort_Dd(const t_gamb_Dd& A, const t_gamb_Dd& B);


	bool _check_amb_depend(double pdE[], double pdC[], int iNamb, int iN_oneway, vector<int> arriIpt2ow, int& iNdef);
}

#endif