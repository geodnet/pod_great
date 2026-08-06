/**
*
* @verbatim
History
-1.0 zhangwei  2019-03-18  creat
-1.1 Qinyujie  2019-03-20  optim
-1.2 Qinyujie  2019-04-11  add Doxygen style code
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorblsq_leo.cpp
* @brief		lsq processing for leo
*
* @author       zhangwei, Wuhan University
* @version		1.0.0
* @date		    2019-03-18
*
*/
#include "gproc/gpodfltIF.h"
#include "gset/gsetpar.h"
#include "gset/gsetout.h"
#include "gset/gsetinp.h"
#include "gset/gsetforce.h"
#include "gutils/ginfolog.h"
#include "gutils/gturboedit.h"
#include "gutils/gconst.h"
#include "gio/gfile.h"
#include "gutils/gstring.h"
#include <algorithm>
#include "gcoders/rinexo.h"
#include "gcoders/rinexn.h"
#include "gcoders/sp3.h"
#include "gproc/goi_nav.h"
#include "gmodels/gstochasticmodel.h"
#include <chrono>






namespace great
{
	t_gpodfltIF::t_gpodfltIF()
	{

	}

	/** @brief default constructor. */
	t_gpodfltIF::t_gpodfltIF(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		t_gfltproc(set, data, log),
		_result(0)
	{

		string tmp;
		tmp = dynamic_cast<t_gsetout*>(set)->outputs("flt");
		if (!tmp.empty())
		{
			substitute(tmp, "$(rec)", _crt_rec, false);
			_result = new t_giof;
			_result->tsys(t_gtime::GPS);
			_result->mask(tmp);
			_result->append(dynamic_cast<t_gsetout*>(set)->append());
		}

        _def_ECOM2.insert(ICS_PAR::PXSAT);
        _def_ECOM2.insert(ICS_PAR::PYSAT);
        _def_ECOM2.insert(ICS_PAR::PZSAT);
        _def_ECOM2.insert(ICS_PAR::VXSAT);
        _def_ECOM2.insert(ICS_PAR::VYSAT);
        _def_ECOM2.insert(ICS_PAR::VZSAT);
        _def_ECOM2.insert(ICS_PAR::ECOM2_D0);
        _def_ECOM2.insert(ICS_PAR::ECOM2_Y0);
        _def_ECOM2.insert(ICS_PAR::ECOM2_X0);
        _def_ECOM2.insert(ICS_PAR::ECOM2_D2c);
        _def_ECOM2.insert(ICS_PAR::ECOM2_D2s);
        _def_ECOM2.insert(ICS_PAR::ECOM2_D4c);
        _def_ECOM2.insert(ICS_PAR::ECOM2_D4s);
        _def_ECOM2.insert(ICS_PAR::ECOM2_X1c);
        _def_ECOM2.insert(ICS_PAR::ECOM2_X1s);
    }

	/** @brief default destructor. */
	t_gpodfltIF::~t_gpodfltIF()
	{

		if (_result)
		{
			if (_result->is_open()) _result->close();
			delete _result;
			_result = nullptr;
		}

		if (_satclkfile){
			if (_satclkfile->is_open()) _satclkfile->close();
			delete _satclkfile;
			_satclkfile = nullptr;
		}
	}

	bool t_gpodfltIF::ProcessBatch(t_gallproc * data, const t_gtime & beg, const t_gtime & end)
	{
		// =================================================================================================
		_glog->logInfo("t_gpodfltIF", "ProcessBatch", "Begin t_gpodfltIF::ProcessBatch");

		if (!_check_initial_information(data, beg, end))   return false;

		//Init data (include SITE/GNSSdata)
		if (!_initfltProcData(data))                       return false;

		//Init pars (include SITE/GNSS pars)
		if (!_initfltProcPars(_lsq))                       return false;

		_oi_msec =0;
        
        _initFltProc(data);
        

#ifdef USE_OPENMP
		dynamic_cast<t_gprecisebias*>(_bias_model.get())->set_multi_thread(_rec_list);
		omp_set_num_threads(_num_threads);

#endif


		chrono::high_resolution_clock::time_point beg_total = chrono::high_resolution_clock::now();
		_remove_par_msec = 0;
		_cmb_equ_msec = 0;
		_preprocess_msec = 0;
		_measure_msec = 0;
		_timeupdate_msec=0;
		_ambfix_msec=0;

		// =================================================================================================
		
		OFSTREAM_MODE rom = dynamic_cast<t_gsetproc*>(_gset)->read_ofile_mode();
		int epo_count = 0;
		int epo_valid = 0;
		t_gtime last_proc_time = _crt_time- 300;
		while (_crt_time <= _end_time)
		{
			// ===============================================================================
		
			epo_count++;

	
			if (!_is_convergence && _crt_time - _beg_time >43200) {
				_is_convergence = true;
				_remove_par_msec = 0;
				_cmb_equ_msec = 0;
				_preprocess_msec = 0;
				_measure_msec = 0;
				_timeupdate_msec = 0;
				_ambfix_msec = 0;
			}

			// ===============================================================================
			
			chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
			_glog->logInfo("t_gpodfltIF", "ProcessBatch", _crt_time.str_ymdhms("Process Epoch  ==> ") + " start");

        

			// ===============================================================================
			/* get observation */
			if ((rom == OFSTREAM_MODE::REAL_TIME))
			{
				_read_file_realtime();
			}

			chrono::high_resolution_clock::time_point beg_pre = chrono::high_resolution_clock::now();
			_preprocess_obs();
			chrono::high_resolution_clock::time_point end_pre = chrono::high_resolution_clock::now();
			_preprocess_msec +=chrono::duration_cast<chrono::milliseconds>(end_pre - beg_pre).count(); 
			

			vector<t_gsatdata> crt_obs = _gall_obs->obs(_rec_list, _crt_time);
			
			// ===============================================================================
	
			bool epoch_valid = _processOneEpoch(_crt_time,crt_obs);
			if (!epoch_valid)
			{
				write_log_info(_glog, 0, "ERROR" + _crt_time.str_ymdhms("Processing failed in epoch "));
				_crt_time = _crt_time + _obs_intv;
				continue;
			}

			// ===============================================================================

				if (_equ4one_epoch.num_equ() == 0) {
					write_log_info(_glog, 0, "WARNING" + _crt_time.str_ymdhms("No Observation for this epoch "));
					_crt_time = _crt_time + _obs_intv;
					continue;
				}
				chrono::high_resolution_clock::time_point beg_measure = chrono::high_resolution_clock::now();
				_measure_update();
				chrono::high_resolution_clock::time_point end_measure = chrono::high_resolution_clock::now();
				// ambfix
				if ( _ambfix_cmd || _ambfix_xml) 
				{
					_ambfixDd_SRIF();

				}
				else{
					ColumnVector dx = eigen2newmat(_flt._x);
					_print_clk_solved(_lsq->_x_solve,dx,_crt_time);
					
				}
			

			if ((rom == OFSTREAM_MODE::REAL_TIME))
			{
				_deleteOfiles();
			}

			// ===============================================================================
		
			chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
			_glog->logInfo("t_gpodfltIF", "ProcessBatch",  _crt_time.str_ymdhms("Process Epoch  ==> ") + " successed");
			cout << _crt_time.str_mjdsod("PROCESS ==> ") + "  Totally Spends time is " << chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count() / 1000.0 << " sec " << endl;
			// ===============================================================================
		
			_crt_time = _crt_time + _obs_intv;
			epo_valid++;
		}
		
		chrono::high_resolution_clock::time_point end_total = chrono::high_resolution_clock::now();
		writeLogInfo(_glog, 0, "NOTE", "###AMBFIX_PROC"+dbl2str(_ambfix_msec/1000.0) +" sec.");
		writeLogInfo(_glog, 0, "NOTE", "###COMBINE_EQU "+dbl2str(_cmb_equ_msec/1000.0) +" sec.");
		writeLogInfo(_glog, 0, "NOTE", "###PREPROCESS"+dbl2str(_preprocess_msec/1000.0) +" sec.");
		writeLogInfo(_glog, 0, "NOTE", "###MEASURE"+dbl2str(_measure_msec/1000.0) +" sec.");
		writeLogInfo(_glog, 0, "NOTE", "###TIMEUPDATE"+dbl2str(_timeupdate_msec/1000.0) +" sec.");
		writeLogInfo(_glog, 0, "NOTE", "###OI"+dbl2str(_oi_msec/1000.0) +" sec.");
		writeLogInfo(_glog, 0, "NOTE", "###TOTAL"+dbl2str(chrono::duration_cast<chrono::milliseconds>(end_total - beg_total).count() /1000.0) +" sec.");
		_glog->logInfo("t_gpodfltIF", "ProcessBatch", "valid epoch is :" + dbl2str(1.0*epo_valid/epo_count*100) + " % ");

		// ====================================================================================
		
		// update ics file
		
		_orb->test_update(_beg_time,_crt_time-_obs_intv,_obs_intv,_estimated_state_map);

		t_gfile outorbfile;
		outorbfile.path(dynamic_cast<t_gsetout*>(_gset)->outputs("orb"));
		t_orb coder_orb(_gset);
		coder_orb.add_data("ID1", _orb);
		outorbfile.coder(&coder_orb);
		outorbfile.run_write();

		// output_sp3
		t_gfile outsp3file;
		t_gallprec sp3_data;
		vector<string> sats_temp;
		for (auto& sat_item: _estimated_state_map_trs_sp3){
			auto sat = sat_item.first;
			sats_temp.push_back(sat);
			for (auto& time_item: sat_item.second){
				sp3_data.addpos(sat,time_item.first,time_item.second,0.0,t_gtriple(0,0,0),0.9);
				sp3_data.add_interval(sat,_obs_intv);
				sp3_data.add_agency("GREAT");
			}
		}
		sp3_data.set_beg(_beg_time);
		sp3_data.set_end(_end_time);
		sp3_data.set_sat(sats_temp);
		sp3_data.set_sat3(sats_temp);
		sp3_data.set_data_type("P");
		sp3_data.set_agency("GREAT");
		sp3_data.set_sat_type("GNS");
		outsp3file.path(dynamic_cast<t_gsetout*>(_gset)->outputs("sp3"));
		t_sp3 coder_sp3(_gset);
		coder_sp3.add_data("ID1", &sp3_data);
		outsp3file.coder(&coder_sp3);
		outsp3file.run_write();



		return true;
		
	}


	bool t_gpodfltIF::_check_initial_information(t_gallproc* data, const t_gtime& beg, const t_gtime& end)
	{
		
		if (data == nullptr)
		{
			_glog->logError("t_gpodfltIF", "ProcessBatch", "t_gallproc is nullptr");
			return false;
		}

		if (beg > end)
		{
			_glog->logError("t_gpodfltIF", "ProcessBatch", "beg > end");
			_glog->logError("t_gpodfltIF", "ProcessBatch", beg.str_mjdsod("beg is ==> "));
			_glog->logError("t_gpodfltIF", "ProcessBatch", end.str_mjdsod("end is ==> "));
			return false;
		}

	
		// get the time
		_beg_time = (_beg_time > beg) ? beg : _beg_time;
		_end_time = (_end_time < end) ? end : _end_time;
		_crt_time = beg;
		_glog->logInfo("t_gpodfltIF", "ProcessBatch", _beg_time.str_mjdsod("_beg_time is ==> "));
		_glog->logInfo("t_gpodfltIF", "ProcessBatch", _end_time.str_mjdsod("_end_time is ==> "));
		_glog->logInfo("t_gpodfltIF", "ProcessBatch", _crt_time.str_mjdsod("_crt_time is ==> "));

		return true;
	}

	bool t_gpodfltIF::GenerateProduct()
	{
		
		return true;

	}


	bool t_gpodfltIF::_initfltProcPars(t_glsq* lsq)
	{
		if (_rec_list.empty() || _sat_list.empty())
		{
			write_log_info(_glog, 0, "ERROR", "no rec or no sat in POD_LSQ.");
			return false;
		}
		//Common parts
		bool rec_valid = _init_rec_clk_pars(lsq);
		bool sys_valid = _init_sys_bia_pars(lsq);

		bool ifb_valid = true;


		if (!rec_valid)
		{
			write_log_info(_glog, 0, "ERROR", "Init the rec clk pars failed");
			return false;
		}
		if (!ifb_valid)
		{
			write_log_info(_glog, 0, "ERROR", "Init the ifb pars failed");
			return false;
		}
		if (!sys_valid)
		{
			write_log_info(_glog, 0, "ERROR", "Init the sys bias pars failed");
			return false;
		}


		//GNSS parts
		if (_mode.count("S") > 0)
		{
			bool ics_valid = _init_gns_ics_pars(lsq);
			bool erp_valid = _init_crs_erp_pars(lsq);
			bool crd_valid = _init_rec_crd_pars(lsq, false);
			bool trop_valid = _init_rec_trop_pars(lsq);
			bool clk_valid = _init_sat_clk_pars(lsq);
			if (!ics_valid)
			{
				write_log_info(_glog, 0, "ERROR", "Init the GNSS ICS pars failed");
				return false;
			}
			if (!erp_valid)
			{
				write_log_info(_glog, 0, "ERROR", "Init the ERP pars failed");
				return false;
			}
			if (!crd_valid)
			{
				write_log_info(_glog, 0, "ERROR", "Init the CRD pars failed");
				return false;
			}
			if (!trop_valid)
			{
				write_log_info(_glog, 0, "ERROR", "Init the TROP pars failed");
				return false;
			}
			if (!clk_valid)
			{
				write_log_info(_glog, 0, "ERROR", "Init the sat clk pars failed");
				return false;
			}
		}
		write_log_info(_glog, 0, "INFO", "finish t_gpodfltIF::_initLsqProcPars");
		return true;
	}

	bool t_gpodfltIF::_initfltProcData(t_gallproc * data)
	{
		_glog->logInfo("t_gpodfltIF", "ProcessBatch", "Begin t_gpodfltIF::_initfltProcData");

		// prepare obs data
		_gall_bias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gall_obs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		_gall_nav = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		_gall_obj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		_gdata_erp = dynamic_cast<t_gerp*>((*data)[t_gdata::ALLERP]);
		_gdata_ambcon = dynamic_cast<t_gambcon*>((*data)[t_gdata::AMBCON]);
		_gdata_satinfo= dynamic_cast<t_gallsatinfo*>((*data)[t_gdata::SATPARS]); 

		OFSTREAM_MODE rom = dynamic_cast<t_gsetproc*>(_gset)->read_ofile_mode();
		if ((rom == OFSTREAM_MODE::REAL_TIME))
		{
			_read_file_realtime();
			_gall_obj->read_satinfo(_beg_time);                                                                       	
			_gall_obj->sync_pcvs();                                                                           
		}
		
		if (!_gall_bias)
		{
			_glog->logInfo("t_gpodfltIF", "_initfltProcData", "ALLBIAS is empty!");

			cerr << "t_gpodfltIF  _initfltProcData  ALLBIAS is empty!" << endl;
			return false;
		}
		if (!_gall_obs)
		{
			_glog->logInfo("t_gpodfltIF", "_initfltProcData", "ALLOBS is empty!");
			return false;
		}
		if (!_gall_obj)
		{
			_glog->logInfo("t_gpodfltIF", "_initfltProcData", "ALLOBJ is empty!");
			return false;
		}

		if (!_gdata_erp)
		{
			_glog->logInfo("t_gpodpodIF", "_initfltProcData", "ALLERP is empty!");
			return false;
		}
		if (_mode.count("S") > 0) 
		{
			_gall_gns = dynamic_cast<t_gallics_nav*>((*data)[t_gdata::ALLICS]);
			if (!_gall_gns || _gall_gns->isEmpty())
			{
                bool init_valid = false;
                
                init_valid = _init_gns_ics_with_nav();
                

                if (!init_valid){
                    _glog->logInfo("t_gpodfltIF", "_initfltProcData", "ALLICS is empty!");
                    return false;
                }
			}

		}

		
		set<string> sat_list;
		set<string>::const_iterator sat_beg = _sat_list.begin();
		set<string>::const_iterator sat_end = _sat_list.end();
		set<string>::const_iterator ics_beg = _gall_gns->sat_list.begin();
		set<string>::const_iterator ics_end = _gall_gns->sat_list.end();
		set_intersection(sat_beg, sat_end, ics_beg, ics_end,
			inserter(sat_list, sat_list.begin()));
		_sat_list = sat_list;
		if (_sat_list.empty())
		{
			_glog->logInfo("t_gpodfltIF", "_initfltProcData", "_sat_list is empty!");
			return false;
		}
		

		auto tb_slip = dynamic_cast<t_gturboedit*>(_slip12.get());

		write_log_info(_glog, 0, "INIT : process data in lsq is prepared.");
		return true;
	}

	void t_gpodfltIF::_initFltProc(t_gallproc * data){


		
      

		int norbpar = 0;
		for (auto iter = _sat_list.begin(); iter != _sat_list.end(); iter++) {
			norbpar += _orb_par_num[*iter];
		}

		_ambfix_flt = make_shared<t_gambfixDd>(_glog, _gset);


        t_gtime orb_state_ref = _beg_time;
       
        int npar = _lsq->npar_number();
        Eigen::VectorXd vecR(npar);
        for (int i = 0; i < npar; i++) {
            vecR(i) = 1 / _lsq->_x_solve[i].apriori();
        }
        _flt.init(vecR.asDiagonal(), Eigen::VectorXd::Zero(npar));
		
        

		_flt.set_qr_solver(_num_threads);
		_flt.set_max_res_norm(dynamic_cast<t_gsetproc*>(_gset)->max_res_norm());


	
		for (const auto& sat : _sat_list) {
			_orb_state_ref[sat] =  orb_state_ref;

			_update_orb_interval[sat] = dynamic_cast<t_gsetforce*>(_gset)->getOutStep(sat);
			double pv_noise = dynamic_cast<t_gsetforce*>(_gset)->getPvNoise(sat);
			double srp_noise = dynamic_cast<t_gsetforce*>(_gset)->getSrpNoise(sat);
			double drag_noise = dynamic_cast<t_gsetforce*>(_gset)->getDragNoise(sat);
			double emp_noise = dynamic_cast<t_gsetforce*>(_gset)->getEmpNoise(sat);

			_orb_mode[sat] = make_shared<t_orb_mode>(make_pair(6, pv_noise), make_pair(_orb_par_num[sat] - 6, srp_noise), make_pair(0, drag_noise), make_pair(0, emp_noise));
					
			_sat_valid_obs_length[sat] = 0;
			_convergence[sat] = false;
		}
		// set proc int to update_par
		_lsq->get_update_par()->set_rec_interval(_proc_int);

		_init_OI(data);

		// init rt-turboedit settings
		for (auto it_site = _rec_list.begin(); it_site != _rec_list.end();it_site++) {
			_gturbo[*it_site] = make_shared<t_gexeturboedit>(_gset, _glog, data);
			t_gtriple site_crd = _gall_obj->obj(*it_site)->crd(_beg_time);

			if (!site_crd.zero()) {
				_gturbo[*it_site]->set_site_crd(*it_site, site_crd);
			}
		}
		

		// init rt_clk output file
		string satclk_file = "clk_" + _beg_time.str_yyyydoy() + "_epo";
		_satclkfile = new t_giof;
		_satclkfile->tsys(t_gtime::GPS);
		_satclkfile->mask(satclk_file);
		_satclkfile->append(false);



    }

	void t_gpodfltIF::_init_OI(t_gallproc * data) {

		_orb = dynamic_cast<t_gallorb*>((*data)[t_gdata::ALLORB]);

		vector<string> total_proc_sat_list; 
		vector<int> proc_sat_par_idx;

		if (_mode.count("S") != 0) 
		{
			vector<string> vec_sat_list(_sat_list.begin(), _sat_list.end());
			for (int i = 0; i < vec_sat_list.size(); i++)
			{
				string temp_sat = vec_sat_list[i];
				int idx = _lsq->_x_solve.getParam("", t_gpar::PXSAT, temp_sat);
				if (idx == -1) {
					cerr << "Skip " << temp_sat << " for no PXSAT pars in x_solve" << endl;
					continue;
				}
				total_proc_sat_list.push_back(temp_sat);
				proc_sat_par_idx.push_back(idx);
				
				_update_oi_interval[temp_sat] = dynamic_cast<t_gsetforce*>(_gset)->getInitTime(temp_sat);

				
				_orb_ref[temp_sat] = _orb_state_ref[temp_sat];
				_epo[temp_sat] = (_orb_state_ref[temp_sat] - _orb_ref[temp_sat]) / _update_orb_interval[temp_sat];
				_phi_map[temp_sat].clear();
				_phi_i_map[temp_sat].clear();
				_state_map[temp_sat].clear();
				_eclipse_marker[temp_sat] = false;


				shared_ptr<t_goi> poi = make_shared<t_goi_nav>(temp_sat, _gset);
				auto  gegm = dynamic_cast<t_gegmcoef*>((*data)[t_gdata::ID_TYPE::EGM]);
				auto  gde = dynamic_cast<t_gnavde*>((*data)[t_gdata::ID_TYPE::ALLDE]);
				auto  gerp = dynamic_cast<t_gerp*>((*data)[t_gdata::ID_TYPE::ALLERP]);
				auto  got = dynamic_cast<t_gocean_tide*>((*data)[t_gdata::ID_TYPE::OCEANTIDE]);
				auto  gotl = dynamic_cast<t_gallotl*>((*data)[t_gdata::ID_TYPE::ALLOTL]);
				auto  glep = dynamic_cast<t_gleapsecond*>((*data)[t_gdata::ID_TYPE::LEAPSECOND]);
				auto  gsat = dynamic_cast<t_gallsatinfo*>((*data)[t_gdata::ID_TYPE::SATPARS]);
				auto  gorb = new t_gorb();

				if (gegm) poi->set_data_used(gegm);
				if (gde) poi->set_data_used(gde);
				if (gerp) poi->set_data_used(gerp);
				if (got) poi->set_data_used(got);
				if (gotl) poi->set_data_used(gotl);
				if (glep) poi->set_data_used(glep);
				if (gsat) poi->set_data_used(gsat);
				if (gorb) poi->set_data_used(gorb);
				poi->glog(_glog);

				poi->set_data_used(_gall_gns);
				dynamic_pointer_cast<t_goi_nav>(poi)->prepareOI();
				
				_poi[temp_sat] = poi;
			}
		}

		
		chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
		for (int i =0;i<total_proc_sat_list.size();i++)
		{

			string temp_sat = total_proc_sat_list[i];

			// construct initial orbit state for integration 
			int idx = proc_sat_par_idx[i];
			vector<t_gpar::t_type> state_init_name;
			vector<double> state_init_value;
			for (int j = idx; j < idx + _orb_par_num[temp_sat]; j++) {
				state_init_name.push_back(_lsq->_x_solve[j].parType);
				state_init_value.push_back(_lsq->_x_solve[j].value());
			}

            t_gtime orb_state_ref = _orb_state_ref[temp_sat];

			int oi_interval = dynamic_cast<t_gsetforce*>(_gset)->getInitTime(temp_sat);
			_update_orb(temp_sat, orb_state_ref , orb_state_ref + oi_interval, state_init_name, state_init_value, _orb);
		}
		chrono::high_resolution_clock::time_point end_t= chrono::high_resolution_clock::now();
		_oi_msec +=  chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count();

	}

	void t_gpodfltIF::_measure_update()
	{

		Eigen::MatrixXd A;
		Eigen::VectorXd z;
		
		_equ4one_epoch.changeNomralizedEigen(A, z, _lsq->_x_solve.parNumber());

		int amb_par_num =0;
		for (int i =0;i<_lsq->_x_solve.parNumber();i++){
			if (_lsq->_x_solve[i].parType == t_gpar::AMB_IF){
				amb_par_num++;
			}
		}

        /*
        * the coeff matrix B is constructed according to the init time of orbit integration
        * Here should be changed to the reference tiem of orbit state estimated paramters
        */
		for (const auto& sat : _sat_list) {
			if (_epo[sat] == 0) {
				continue;
			}
			int pxsat_index = _lsq->_x_solve.getParam("", t_gpar::PXSAT, sat);
			if (pxsat_index == -1) {
				continue;
			}
			Eigen::MatrixXd phiMap = newmat2eigen(_phi_i_map[sat][_epo[sat]]);
			A.block(0, pxsat_index, A.rows(), _orb_par_num[sat]) *= phiMap;
		}

        chrono::high_resolution_clock::time_point beg_measure = chrono::high_resolution_clock::now();
		_flt.set_amb_num(amb_par_num);
		_flt.set_phase_num(_equ4one_epoch.phase_equ());
		_flt.set_obs_equ(_equ4one_epoch);
		_flt.set_par_list(_lsq->_x_solve);
		set<int> clk_id;
		for (int i =0;i<_lsq->_x_solve.parNumber();i++){
			if (_lsq->_x_solve[i].parType == t_gpar::CLK 
				|| _lsq->_x_solve[i].parType==t_gpar::CLK_SAT){
				clk_id.insert(i);
			}
		}
		_flt.set_clk_id(clk_id);

      


		Eigen::VectorXd z_e = _flt.measurement_update(A,z);

        set<int> bad_obs_id;
        bad_obs_id = _flt.obs_quality_control(A,z,z_e);
        



		chrono::high_resolution_clock::time_point end_measure = chrono::high_resolution_clock::now();
		_measure_msec+=chrono::duration_cast<chrono::milliseconds>(end_measure- beg_measure).count(); 

       

        // process of bias
        _process_posterior_error(A, z, bad_obs_id);


		_equ4one_epoch.clear_allequ();

		ColumnVector dx = eigen2newmat(_flt._x);

		// change crs2trs
		t_gtrs2crs trs2crs_temp(false,_gdata_erp);
		t_gtime tdt = _crt_time;
		tdt.tsys(t_gtime::TT);
		int eop_index = _lsq->_x_solve.getParam("", t_gpar::XPOLE, "");
		if (eop_index != -1) {
			// EOP parameters are estimated
			ColumnVector eop_corr = dx.SubMatrix(eop_index + 1, eop_index + 6, 1, 1).AsColumn();
			trs2crs_temp.calcRotMat_withCorr(tdt, false, false, false, eop_corr);
		}
		else {
			// don't have eop parameters
			trs2crs_temp.calcRotMat(tdt, false, false, false);
		}
		// save the float result in the estimated_state_map
		for (const auto& sat : _sat_list) {
			int pxsat_index = _lsq->_x_solve.getParam("", t_gpar::PXSAT, sat);
			if (pxsat_index == -1) {
				continue;
			}
			
			if (_orb_state_ref[sat] == _crt_time) {
				_estimated_state_map[sat][_orb_state_ref[sat]] = _state_map[sat][_orb_state_ref[sat]] + dx.SubMatrix(pxsat_index + 1, pxsat_index + _orb_par_num[sat], 1, 1).AsColumn();
				_estimated_state_map_trs_sp3[sat][_orb_state_ref[sat]] = t_gtriple(trs2crs_temp.getRotMat().t() * _estimated_state_map[sat][_orb_state_ref[sat]].rows(1, 3));
			}			
			else {
				int temp_interval = _crt_time - _orb_state_ref[sat];
				if ( temp_interval % _update_orb_interval[sat] !=0 ) {
					// case 1: don't output this epoch; ( this epoch is not in the step_out epoch) 	
					continue;
				}
				int temp_update_epo = temp_interval / _update_orb_interval[sat];
				// case 2: use the results of integration;
				// Here only prcess the [case 2] 
				Matrix phi = _phi_map[sat][_epo[sat] + temp_update_epo] * _phi_i_map[sat][_epo[sat]];
				_estimated_state_map[sat][_crt_time] = _state_map[sat][_crt_time] + phi * dx.SubMatrix(pxsat_index + 1, pxsat_index + _orb_par_num[sat], 1, 1).AsColumn();
				_estimated_state_map_trs_sp3[sat][_crt_time] = t_gtriple(trs2crs_temp.getRotMat().t() * _estimated_state_map[sat][_crt_time].rows(1, 3));
			}
		}
	}


	void t_gpodfltIF::_preprocess_obs()
	{
		vector<string> vrec;
		for (const auto& rec : _rec_list) {
			// no observation
			if (_gall_obs->end_obs(rec) < _crt_time) continue;

			vrec.push_back(rec);
		}
		OFSTREAM_MODE rom = dynamic_cast<t_gsetproc*>(_gset)->read_ofile_mode();
		if (_slip_model == SLIPMODEL::TURBO_EDIT)
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
			for (int i = 0; i < vrec.size(); ++i) 
			{
				_gturbo.at(vrec[i])->ProcessBatch(vrec[i], _crt_time, _crt_time, _obs_intv, true);
			}

	}

	void t_gpodfltIF::_ambfixDd_SRIF()
	{
		
		
		// get current amb pars
		t_gallpar ambPars;
		t_gtime par_beg=LAST_TIME;
		vector<int> debug_amb_index;
		for (int i = 0; i < _lsq->_x_solve.parNumber();i++) {
			if (_lsq->_x_solve[i].parType == t_gpar::AMB_IF) {
				// this ambiguity parameter don't have any observation info
				if (_lsq->_x_solve[i].value() == 0.0) {
					continue;
				}
				t_gpar tempPar = _lsq->_x_solve[i];
				tempPar.value(tempPar.value() + _flt._x(i));
				// RT NL ambiguity sigma use 0.05
				tempPar.apriori(0.0);
				tempPar.end = _crt_time;
				ambPars.addParam(tempPar);
				debug_amb_index.push_back(i);
				par_beg = min(par_beg, tempPar.beg);
			}
		}
		chrono::high_resolution_clock::time_point beg_t, end_t;
		int wl_interval = dynamic_cast<t_gsetamb*>(_gset)->wl_interval();
		if (!_is_convergence) {
			beg_t = chrono::high_resolution_clock::now();
			_ambfix_flt->add_glo_fid(_glofrq_num);
			_ambfix_flt->update_data(&ambPars, _gall_obs, _gall_bias, _gall_obj, _crt_time - _obs_intv + wl_interval, _crt_time);
			end_t = chrono::high_resolution_clock::now();
			return;
		}


		beg_t = chrono::high_resolution_clock::now();
		_ambfix_flt->add_data(_gall_obs, _gall_bias, &ambPars , _gall_obj);
		_ambfix_flt->add_glo_fid(_glofrq_num);
		_ambfix_flt->ProcessBatch(_crt_time - _obs_intv+wl_interval,_crt_time);

		end_t = chrono::high_resolution_clock::now();
		_ambfix_msec+=chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count(); 


		// Add virutal equations
		t_SRIF srif_fix = _flt;
		t_gambcon gdata_ambcon = _ambfix_flt->get_ambcon_results();
		beg_t = chrono::high_resolution_clock::now();

		_pod_add_Dd_ambcon(&srif_fix, &gdata_ambcon);

		end_t = chrono::high_resolution_clock::now();
		_ambfix_msec += chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count();

		ColumnVector dx = eigen2newmat(srif_fix._x);

		// change crs2trs
		t_gtrs2crs trs2crs_temp(false, _gdata_erp);
		t_gtime tdt = _crt_time;
		tdt.tsys(t_gtime::TT);
		int eop_index = _lsq->_x_solve.getParam("", t_gpar::XPOLE, "");
		if (eop_index != -1) {
			// EOP parameters are estimated
			ColumnVector eop_corr = dx.SubMatrix(eop_index + 1, eop_index + 6, 1, 1).AsColumn();
			trs2crs_temp.calcRotMat_withCorr(tdt, false, false, false, eop_corr);
		}
		else {
			// don't have eop parameters
			trs2crs_temp.calcRotMat(tdt, false, false, false);
		}
		// save the fixed result in the estimated_state_map
		for (const auto& sat : _sat_list) {
			int pxsat_index = _lsq->_x_solve.getParam("", t_gpar::PXSAT, sat);
			if (pxsat_index == -1) {
				continue;
			}

			if (_orb_state_ref[sat] == _crt_time) {
				_estimated_state_map[sat][_orb_state_ref[sat]] = _state_map[sat][_orb_state_ref[sat]] + dx.SubMatrix(pxsat_index + 1, pxsat_index + _orb_par_num[sat], 1, 1).AsColumn();
				_estimated_state_map_trs_sp3[sat][_orb_state_ref[sat]] = t_gtriple(trs2crs_temp.getRotMat().t() * _estimated_state_map[sat][_orb_state_ref[sat]].rows(1, 3));
			}
			else {
				int temp_interval = _crt_time - _orb_state_ref[sat];
				if (temp_interval % _update_orb_interval[sat] != 0) {
					// case 1: don't output this epoch; ( this epoch is not in the step_out epoch) 	
					continue;
				}
				int temp_update_epo = temp_interval / _update_orb_interval[sat];
				// case 2: use the results of integration;
				// Here only prcess the [case 2] 
				Matrix phi = _phi_map[sat][_epo[sat] + temp_update_epo] * _phi_i_map[sat][_epo[sat]];
				_estimated_state_map[sat][_crt_time] = _state_map[sat][_crt_time] + phi * dx.SubMatrix(pxsat_index + 1, pxsat_index + _orb_par_num[sat], 1, 1).AsColumn();
				_estimated_state_map_trs_sp3[sat][_crt_time] = t_gtriple(trs2crs_temp.getRotMat().t() * _estimated_state_map[sat][_crt_time].rows(1, 3));
			}
		}
		


		_print_clk_solved(_lsq->_x_solve,dx,_crt_time);

	}


	void t_gpodfltIF::_print_clk_solved(t_gallpar& par, ColumnVector& dx, t_gtime& crt_epoch)
	{
		//set<string>::const_iterator rec_end = _rec_list.end();
		set<string>::const_iterator sat_end = _sat_list.end();
		//set<string>::const_iterator it_rec = _rec_list.begin();
		set<string>::const_iterator it_sat = _sat_list.begin();

		ostringstream os;
		t_gtime creat_t(t_gtime::UTC);
		creat_t.add_secs(28800);

		string sat_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		// write sp3 clk file header
		if (!_satclkfile->is_open()) {
			os << "     2.00           C                                       RINEX VERSION / TYPE" << endl;

			os << left << setw(20) << "G-Nut" << setw(20) << "GREAT" << setw(20) << creat_t.str("%Y%m%d %H%M%S") << "PGM / RUN BY / DATE" << endl;

			os << right << setw(9) << fixed << setprecision(2) << _obs_intv << setw(51) << " " << "INTERVAL" << endl;

			// write # / TYPES OF DATA
			os << "    " << right << setw(2) << "1" << "    " << "AS" << setw(48) << "" << "# / TYPES OF DATA" << endl;

			// write # OF CLK REF
			os << "    " << right << setw(2) << "1" << setw(54) << "" << "# OF CLK REF" << endl;

			// write ANALYSIS CLK REF
			os << setw(4) << sat_ref << setw(56) << " " << "ANALYSIS CLK REF" << endl;

			// write END OF HEADER
			os << setw(60) << " " << "END OF HEADER" << endl;

			_satclkfile->write(os.str().c_str(), os.str().size());
			_satclkfile->flush();

			os.str("");

		}

		// write sp3 clk file data
		double satclk;
		for (it_sat; it_sat != sat_end; it_sat++) {
			int idx = par.getParam("", t_gpar::CLK_SAT, *it_sat, crt_epoch, crt_epoch);
			if (idx < 0) continue;
			int idx_ifcb;

			double x = par[idx].value();
			double dclk = dx(idx + 1);

			double clk;
			if (dclk != 0) {
				//continue;
				clk = x + dclk;
				satclk = clk / CLIGHT;

				os << "AS" << " " << setw(3) << *it_sat << "  " <<
					setw(26) << crt_epoch.str("%Y %C %L %K %O %P") << " " << setw(2)
					<< "1" << "   " << setw(19) << fixed << scientific << setprecision(12) << satclk << endl;
			}

		}

		// Print sp3 clk file data
		if (_satclkfile) {
			_satclkfile->write(os.str().c_str(), os.str().size());
			_satclkfile->flush();
			os.str("");
		}
		else {
			_glog->comment(2, "t_gpcelsqIF::_printout", "have no output file!");
		}
	}


	void t_gpodfltIF::_process_posterior_error(const Eigen::MatrixXd &A, const Eigen::VectorXd &z, set<int> bias_idx)
	{


        // statistic bad data
        set<string> bad_sats;
        map<string,set<string> > sats_bias_code; //Count the number of phase bias for each sats
        map<string, set<string> > sats_bias_phase; // Count the ratio of code bias for each sats
        map<string, int> sats_bias_both;  // Count the number of both phase and code bias for each sats (the union of code and phase)
        for (auto item : bias_idx) {
            string site = _equ4one_epoch.get_sitename(item);
            string sat = _equ4one_epoch.get_satname(item);
            bad_sats.insert(sat);

            if (_equ4one_epoch.get_obscombtype(item).is_code()) {
                sats_bias_code[sat].insert(site);
            }
            if (_equ4one_epoch.get_obscombtype(item).is_phase()) {
                sats_bias_phase[sat].insert(site);
            }
        }
        for (auto& sat_item : bad_sats) {
            if (sats_bias_phase.count(sat_item) || sats_bias_code.count(sat_item)) {
                set<string> merge_sites;
                set_union(sats_bias_code[sat_item].begin(), sats_bias_code[sat_item].end(),
                          sats_bias_phase[sat_item].begin(), sats_bias_phase[sat_item].end(),
                          inserter(merge_sites,merge_sites.begin()));
                sats_bias_both[sat_item] = merge_sites.size();
            }
        }



        t_gallpar par_remove_before = _lsq->_x_solve;

        // update new amb for qualtiy control
		t_gupdateparinfo ambupdate_info;
		if (!bias_idx.empty()) {
			for (auto item : bias_idx) {
				string site = _equ4one_epoch.get_sitename(item);
				string sat = _equ4one_epoch.get_satname(item);
				if (!_equ4one_epoch.get_obscombtype(item).is_phase()) {
					continue;
				}
				
				int par_idx = _lsq->_x_solve.getParam(site, t_gpar::AMB_IF, sat);
				ambupdate_info.add(par_idx + 1);
			}
			vector<int> rm_idx;
			ambupdate_info.get(rm_idx);
			_flt.time_update_use_removeinfo(ambupdate_info);

          

			sort(rm_idx.begin(), rm_idx.end(), greater<int>());
			for (auto item : rm_idx) {
				_lsq->remove_parameter(item);
			}
		}


	}


	void t_gpodfltIF::_deleteOfiles()
	{
		map<string,t_gtime> obs_used_beg;
		for (auto site: this->_rec_list){
			obs_used_beg[site]=LAST_TIME;
		}
		for (int i = 0; i < _lsq->_x_solve.parNumber();i++) {
			if (_lsq->_x_solve[i].parType == t_gpar::AMB_IF) {
				t_gpar& amb_par = _lsq->_x_solve[i];
				obs_used_beg[amb_par.site] = min(obs_used_beg[amb_par.site], amb_par.beg);
			}
		}
		for (auto& items:obs_used_beg){
			_gall_obs->clean_outer(items.first,items.second);
		}
	}

	void t_gpodfltIF::_read_file_realtime()
	{
		

		// read ofile
		auto site_o_file_name = [&](const string& site)->string {
			string filename = dynamic_cast<t_gsetproc*>(this->_gset)->read_ofile_dir();
			filename = str_repaceall(filename, "-YYYY-", _crt_time.str_yyyy());
			filename = str_repaceall(filename, "-DDD-", _crt_time.str_doy());
			filename = str_repaceall(filename, "-YY-", _crt_time.str_yr());
			filename = str_repaceall(filename, "-SITE-", str2lower(site));

			return filename;
			};
		int cnt = 0;

		int obs_interval = dynamic_cast<t_gsetgen*>(this->_gset)->sampling();


		// reading ofiles according to site list 
		for (auto site : _rec_list) {

			// if sites obs time is outdate then read the new file
			if (_gall_obs->end_obs(site) < _crt_time) {

				string obs_f = site_o_file_name(site);

				// already read then skip
				if (_site_used_ofile[site] == obs_f) {
					continue;
				}

				string path(obs_f);

				string id("ID" + int2str(++cnt));

				t_rinexo gcoder(_crt_time, LAST_TIME, _gset, "", 20480, obs_interval);
				gcoder.clear();
				gcoder.path(path);
				gcoder.glog(_glog);
				gcoder.add_data(id, _gall_obs);
				gcoder.add_data("OBJ", _gall_obj);

				t_gfile gio;
				gio.glog(_glog);
				gio.path(path);
				gio.coder(&gcoder);

				gio.run_read();

				_site_used_ofile[site] = obs_f;
			}
		}

		// read navfile
		auto navfile_list = dynamic_cast<t_gsetinp*>(this->_gset)->inputs(RINEXN_INP);
		for (auto nav_f : navfile_list) {
			string original_nav_f = nav_f;
			nav_f = str_repaceall(nav_f, "-YYYY-", _crt_time.str_yyyy());
			nav_f = str_repaceall(nav_f, "-YY-", _crt_time.str_yr());
			nav_f = str_repaceall(nav_f, "-DDD-", _crt_time.str_doy());

			if (_used_navfile[original_nav_f] == nav_f) {
				continue;
			}
			_used_navfile[original_nav_f] = nav_f;
			string id("ID" + int2str(++cnt));
			
			t_rinexn gcoder(_gset);
			gcoder.clear();
			gcoder.path(nav_f);
			gcoder.glog(_glog);
			gcoder.add_data(id, _gall_nav);
			gcoder.add_data("OBJ", _gall_obj);
			
			t_gfile gio;
			gio.glog(_glog);
			gio.path(nav_f);
			gio.coder(&gcoder);

			gio.run_read();
		}



		// update glonass_id
		_glofrq_num = _gall_obs->glo_freq_num();
		if (_glofrq_num.empty())
		{
			_glofrq_num = _gall_nav->glo_freq_num();
		}
	}

    bool t_gpodfltIF::_init_gns_ics_with_nav() {
        if (!_gall_nav) return false;

        if (_gall_gns == nullptr){
            _gall_gns = new t_gallics_nav();
            _gall_gns->glog(_glog);
        }

        // init ics data in _gall_gns with broadcast ephemeris
        _gall_gns->beg_time = _beg_time;
        _gall_gns->end_time = _end_time;
        _gall_gns->sat_num = _sat_list.size();
        _gall_gns->sat_list = _sat_list;
        _gall_gns->sys_list = _sys_list;
        _gall_gns->sys_type = "GNS";

        // init ics data
        for (auto sat_iter = _sat_list.begin(); sat_iter != _sat_list.end();)
        {
            string sat = *sat_iter;

            t_gics_nav ics(sat);
            set<ICS_PAR> par_list;
            auto force = dynamic_cast<t_gsetforce*>(_gset)->getForceModelMap(sat);
            string model = force[SOLAR_RADI];
             par_list = _def_ECOM2;

            // get the vel and pos of sat
            double xyz[3] = { 0.0 };
            double var[3] = { 0.0 };
            double vel[3] = { 0.0 };


            _gall_nav->pos(sat, _beg_time, xyz, var, vel, true);

            if (double_eq(xyz[0], 0.0) && double_eq(xyz[1], 0.0) && double_eq(xyz[2], 0.0))
            {
                _glog->comment(t_glog::LOG_LV::LOG_WARN, "t_gpodfltIF", "_init_gns_ics_with_nav", "can not get the xyz from sp3 or rinexn file " + sat);
                _sat_list.erase(sat_iter++);
                _gall_gns->sat_list = _sat_list;
                _gall_gns->sat_num = _sat_list.size();
                continue;
            }

            t_gtime tdt_beg = _beg_time;
            tdt_beg.tsys(t_gtime::TT);

            // get the current trs2crs
            t_gtrs2crs trs2crs_temp(false, _gdata_erp);
            // don't have eop parameters
            trs2crs_temp.calcRotMat(tdt_beg, true, true, true);

            ColumnVector trs_xyz(xyz, 3);
            ColumnVector trs_vel(vel, 3);
            Matrix trs2crs = trs2crs_temp.getRotMat();
            ColumnVector crs_xyz = trs2crs * trs_xyz;

            ColumnVector x_earth = (crs_xyz.t() * trs2crs_temp.getMatDu() / RAD2TSEC).t();
            ColumnVector crs_vel = trs2crs * (trs_vel + x_earth);

            //unit in ics is km
            crs_xyz = crs_xyz / 1E3;
            crs_vel = crs_vel / 1E3;
            for(auto par_iter = par_list.begin(); par_iter != par_list.end(); par_iter++)
            {
                t_ics_record record;
                ICS_PAR par = *par_iter;
                record.beg_time = _beg_time;
                record.end_time = _end_time;
                record.enum_name = par;
                record.str_name = ics2str(par);

                switch (par)
                {
                    case ICS_PAR::PXSAT:
                        record.value = crs_xyz(1);
                        record.sigma = 0.0;
                        break;
                    case ICS_PAR::PYSAT:
                        record.value = crs_xyz(2);
                        record.sigma = 0.0;
                        break;
                    case ICS_PAR::PZSAT:
                        record.value = crs_xyz(3);
                        record.sigma = 0.0;
                        break;
                    case ICS_PAR::VXSAT:
                        record.value = crs_vel(1);
                        record.sigma = 0.0;
                        break;
                    case ICS_PAR::VYSAT:
                        record.value = crs_vel(2);
                        record.sigma = 0.0;
                        break;
                    case ICS_PAR::VZSAT:
                        record.value = crs_vel(3);
                        record.sigma = 0.0;
                        break;
                    case ICS_PAR::ECOM2_D0:
                    case ICS_PAR::ECOM2_Y0:
                    case ICS_PAR::ECOM2_X0:
                    case ICS_PAR::ECOM2_D2c:
                    case ICS_PAR::ECOM2_D2s:
                    case ICS_PAR::ECOM2_D4c:
                    case ICS_PAR::ECOM2_D4s:
                    case ICS_PAR::ECOM2_X1c:
                    case ICS_PAR::ECOM2_X1s:
                        record.value = 0.0;
                        record.sigma = 0.0;
                        break;
                    default:
                        break;
                }

                ics.set_par_record(record.str_name, record);
            }

            _gall_gns->add_nav_ics(sat,ics);

            // next sat
            sat_iter++;
        }

        return true;
    }

 
    void t_gpodfltIF::_add_vel_pulse(Eigen::MatrixXd &A) {

        set<string> pulse_sat_list = {"C24"};

        // add the vel pulse parameters in the [ now ] epoch
        // add parameters on estimator
        Eigen::MatrixXd vel_apriori  = Eigen::MatrixXd::Zero(3*pulse_sat_list.size(),3*pulse_sat_list.size());

        int sat_count= 0;
        int par_index_num = _lsq->_x_solve.parNumber();

        for (auto sat : pulse_sat_list){
            for (int i =0;i<3;i++){
                vel_apriori(sat_count * 3+ i,sat_count  *3 +i) = 1E6;
                t_gpar::t_type pulse_type;
                if (i == 0){
                    pulse_type = t_gpar::VX_Pulse;
               }
                else if (i == 1){
                    pulse_type = t_gpar::VY_Pulse;
                }
                else if (i == 2){
                    pulse_type = t_gpar::VZ_Pulse;
                }

                t_gpar temp_par("",pulse_type,++par_index_num,sat);
                temp_par.setTime(_crt_time,LAST_TIME);
                _lsq->add_parameter(temp_par);
            }
            sat_count ++;
        }
        _flt.add_parameter(vel_apriori,Eigen::VectorXd::Zero( 3 * pulse_sat_list.size()));


        // modify the observation equation
        Eigen::MatrixXd A_pulse = Eigen::MatrixXd::Zero(A.rows(),A.cols() + 3*pulse_sat_list.size());
        A_pulse.block(0,0,A.rows(),A.cols())= A;

        sat_count = 0;
        for (const auto& sat: pulse_sat_list){

            t_gcrs2acr crs2acr;
            Matrix rotmat(0.0,3,3);
            ColumnVector sat_pos(0.0, 3);
            ColumnVector sat_vel(0.0, 3);
            ColumnVector maneuver_sat_state = _state_map[sat][_crt_time];

            sat_pos = maneuver_sat_state.rows(1,3);
            sat_vel = maneuver_sat_state.rows(4,6);

            // change unit from kilometer to meter
            sat_pos *= 1000.0;
            sat_vel *= 1000.0;

            crs2acr.calcRotMat(sat_pos,sat_vel, rotmat);

            int sat_vel_id = _lsq->_x_solve.getParam("",t_gpar::PXSAT,sat);

            for (int i =0;i< _equ4one_epoch.num_equ();i++){
                if (_equ4one_epoch.get_satname(i) != sat){
                    continue;
                }

                ColumnVector partial_xyz(0.0,3);
                partial_xyz(1) = A(i,sat_vel_id + 3);
                partial_xyz(2) = A(i,sat_vel_id + 4);
                partial_xyz(3) = A(i,sat_vel_id + 5);

                ColumnVector partial_acr = rotmat * partial_xyz;

                A_pulse(i,A.cols() + sat_count * 3 + 0 ) = partial_xyz(1);
                A_pulse(i,A.cols() + sat_count * 3 + 1 ) = partial_xyz(2);
                A_pulse(i,A.cols() + sat_count * 3 + 2 ) = partial_xyz(3);
            }

            sat_count++;
        }

        A = A_pulse;

    }

}
