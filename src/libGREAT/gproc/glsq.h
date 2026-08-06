/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2018-12-02  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallprocdata.h
* @brief		header files of least squares estimator
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/
#ifndef GLSQ_H
#define GLSQ_H 


#include "ExportHeaders/ExportBase.h"

#include "newmat/newmat.h"
#include "newmat/newmatap.h"

#include "gall/gallpar.h"
#include "gmodels/gstochasticmodel.h"
#include "gmodels/gmodel.h"
#include "gutils/gsysconv.h"
#include "gutils/gcycleslip.h"
#include "gproc/glsqmatrix.h"
#include "gio/giobigf.h"
#include "gproc/gupdatepar.h"
#include "gproc/ginverse.h"
#include <thread>

namespace great
{

	/**
	*@brief	   Class for lsq estimator include upadte,slove,recover
	*/
	class BASE_LIBRARY_EXPORT t_glsq
	{

	public:

		/** @brief	constructor */
		t_glsq(t_gsetbase* set = nullptr);

		/** @brief  copy constructor */
		t_glsq(const t_glsq& Other);

		/** @brief  default destructor */
		virtual ~t_glsq();

		/**
		* @brief  add now parameter
		* @param[in] par new par preapred to add
		* return -1 for fail 1 for success
		*/
		virtual int add_parameter(const t_gpar& par,bool head=false);

        void clear_all();

		/**
		* @brief  add the partype and the state mode
		* @param[in] par_type type of parameter
		* @param[in] order 0:white noise 1:random walk
		* @param[in] dt [unit:hour]
		* @param[in] noise process noise
		*/
		void add_par_state_equ(t_gpar::t_type par_type, int order, double dt, double noise);

		void set_update_par(shared_ptr<t_gupdatepar> updatepar);

		shared_ptr<t_gupdatepar> get_update_par() {return _update_lsqpar;}


		/**
		* @brief add new observ equations[lsqmatrix format]
		* @param[in] equ equations info
		* @param[in] epoch time of equations
		*/
		int add_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch = t_gtime(),bool write_temp=true);

		

		/**
		* @brief write the equtions information to tempfile for recovering later[lsqmatrix format
		* @param[in] equ equations of lsqmatrix format
		* @param[in] epoch time of equations
		* @param[in] info information of equations
		* @return true for success false for fail
		*/
		bool  write_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch);



		/**
		* @brief remove the specified par in the lsq estimator
		* @note idx begin from 1
		*/
		virtual int remove_parameter(const int& idx);

		/**
		* @brief add apriori for specified parameter
		* @note  idx begin from 1
		*/
		void add_apriori_weight(const int& idx);



		/** @brief reset filename of tempfile */
		void reset_tempfile(string filename);



		/** @brief get nubmer of parameter*/
		int npar_number() const;



		t_gallpar  _x_solve;	///< all parameter in lsq estimator


	protected:
		/**
		* @brief  write the removed NEQ W to tempfile
		* @note idx from1
		*/
		int _write_coefficient(int idx);




		t_glog*	_log;						///< log file

		V_ColumnVector    _W;				///< BTPL Matrix
		V_SymmetricMatrix _NEQ;			    ///< BTPB Matrix

		SymmetricMatrix _Qx;					///< storage Qx after solve
		ColumnVector	_dx;					///< correction of all parameter
		ColumnVector    _dx_final;				///< last correction of all parameter
		ColumnVector	_stdx;					///< covarience of correction of all parameter

		double _res_obs;						///< sum of residuals of currecnt epoch
		double _vtpv;							///< sum of residuals of all epoch
		double _sigma0;							///< error in unit weight
		int _obs_total_num;						///< totoal observation number
		int _npar_tot_num;						///< total par number include ll eliminated parameters 


		t_gtime _epo;						///< current time
		t_gtime	_beg;						///< begin time
		t_gtime	_end;						///< end time
		double  _interval;                     ///< intv 
		t_giobigf*	_tempfile;					///< temp file
		int _buffer_size=1024*1000*10; 		///< tempfile buffer size

		
		shared_ptr<t_gupdatepar> _update_lsqpar;
		shared_ptr<t_ginverse<Matrix> > _solve_matrix;

		t_gmutex _lsq_mtx;

		int size_int = sizeof(int);
		int size_dbl = sizeof(double);

	};

	// mainly for recover or exctract par
	BASE_LIBRARY_EXPORT t_gtime str2gtime(string str_time);


} // namespace

#endif //GLSQ_H

