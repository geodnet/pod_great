/**
*
* @verbatim
	History
	 -1.0 zhangwei     2019-03-18  creat
	 -1.1 Qinyujie     2019-03-20  optim
	 -1.2 Qinyujie     2019-04-11  add Doxygen style code
	 -1.3 Huang Jiande 2019-07-12  add UC model
	 -1.4 Huang Jiande 2020-03-14  fix UC model bug
	 -1.5 Huang Jiande 2020-09-11  format the function
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorblsq_leo.h
* @brief		lsq processing for leo
*
* @author       zhangwei, Wuhan University
* @version		1.0.0
* @date		    2019-03-18
*
*/
#ifndef GPODFLTIF_H
#define GPODFLTIF_H

#include "ExportHeaders/ExportLSQ.h"
#include "gambfix/gallamb.h"
#include "gdata/gambcon.h"
#include "gproc/gexeturboedit.h"
#include "gambfix/gambfixDd.h"
#include "gmodels/gfltproc.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for LSQ processing.
	*/
	class LSQ_LIBRARY_EXPORT t_gpodfltIF : public t_gfltproc
	{
	public:
		//
		t_gpodfltIF();
		//
		t_gpodfltIF(t_gsetbase* base, t_gallproc* data, t_glog* log);
		//
		virtual ~t_gpodfltIF();

	
		bool ProcessBatch(t_gallproc* data, const t_gtime&beg, const t_gtime& end) override;

		// 
		bool GenerateProduct() override;

		// 
		bool initfltProcData(t_gallproc* data) { return _initfltProcData(data); };
		// 
		bool initfltProcPars(t_glsq* lsq) { return _initfltProcPars(lsq); };

	protected:

		// 
		bool _initfltProcPars(t_glsq* lsq)      override;
		//
		bool _initfltProcData(t_gallproc* data) override;
		

		// real-tiemm srif processing
		void _initFltProc(t_gallproc* data);
        bool _init_gns_ics_with_nav();
		void _init_OI(t_gallproc * data);
		void _measure_update();
		void _preprocess_obs();
		void _ambfixDd_SRIF();
		void _print_clk_solved(t_gallpar& par, ColumnVector& dx, t_gtime& crt_epoch);
        void _output_srif_variance();
        void _add_vel_pulse(Eigen::MatrixXd &A);




        // for init gnss ics
        set<ICS_PAR> _def_ECOM2;

		void _process_posterior_error(const Eigen::MatrixXd &A, const Eigen::VectorXd &z, set<int> bias_idx);



	

		map<string, shared_ptr<t_gexeturboedit> > _gturbo;
		shared_ptr<t_gambfixDd> _ambfix_flt;

		map<string, string> _site_used_ofile;
		map<string,string> _used_navfile;

		
        bool _is_init_with_existing_info = false;

		void _read_file_realtime();
		void _deleteOfiles();



		t_giof*	           _result = nullptr;		///< outputs every epoch result file
		t_giof*	           _corr_map = nullptr;


		// 
		bool _check_initial_information(t_gallproc* data, const t_gtime& beg, const t_gtime& end);


	};
}

#endif 