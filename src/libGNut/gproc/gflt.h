
#ifndef FLT_H
#define FLT_H 

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: implements Kalman filter class
  Version: $ Rev: $

  2011-01-10 /PV: created

-*/

#include "newmat/newmat.h"
#include "newmat/newmatap.h"

#include "gall/gallpar.h"
#include "gproc/gupdatepar.h"
#include "gproc/glsqmatrix.h"
#include "ExportHeaders/ExportBase.h"


#include "Eigen/Dense"

#include "gproc/gqr_srif.h"

using namespace std;

namespace gnut
{

	// Base class for all filtering technique
	// -------------------------------------------------------

	class BASE_LIBRARY_EXPORT t_gflt
	{
	public:
		t_gflt();
		virtual ~t_gflt();
	protected:

		ColumnVector _dx;
		SymmetricMatrix _Qx;
		SymmetricMatrix _Qx_tmp;
		SymmetricMatrix _Qx0;
		Matrix _A, _A_virtual;
		SymmetricMatrix _P, _P_virtual;
		ColumnVector _l, _l_virtual;
		t_gallpar _param;
		double _sigma0, _vtpv;
		int _nobs_total, _npar_number;
		bool _amb = false;
	};

	
	// ------------------------------------------------------

	// Square root information filter
	// ------------------------------------------------------
	class BASE_LIBRARY_EXPORT t_SRIF : public t_gflt
	{
	public:
		t_SRIF();
		t_SRIF(const t_SRIF& Others);
		~t_SRIF();

        t_SRIF operator=(const t_SRIF& Others);

		void set_qr_solver(int num_threads);

		Eigen::VectorXd measurement_update(great::t_glsqEquationMatrix measurement_equ);
		Eigen::VectorXd measurement_update(const Eigen::MatrixXd &A, const Eigen::VectorXd &z);
        set<int> obs_quality_control(const Eigen::MatrixXd &A, const Eigen::MatrixXd &z, const Eigen::VectorXd &z_e);
        void obs_quality_control(const Eigen::MatrixXd &A, const Eigen::MatrixXd &z, const Eigen::VectorXd &z_e, set<int> given_idx);

		void time_update_use_removeinfo(great::t_gupdateparinfo remove_info);


		void measurement_update_with_removepar(Eigen::MatrixXd A,Eigen::VectorXd z,vector<int> remove_id);


		void init(const Eigen::MatrixXd& R, const Eigen::VectorXd& x0);

		void add_parameter(const Eigen::MatrixXd& R, const Eigen::VectorXd& x0);
		void rm_param(vector<int> remove_id);

		void set_amb_num(int value){
			qr_solve->amb_num = value;
		}

		void set_phase_num(int value){
			qr_solve->phase_num = value;
		}

		void set_obs_equ(const great::t_glsqEquationMatrix& obs_equ)
		{
			qr_solve->obs_equ = obs_equ;
		}

		void set_clk_id(const set<int>& clk_id){
			qr_solve->clk_id  = clk_id;
		}

		void set_par_list(const t_gallpar& par_list){
			qr_solve->par_list = par_list;
		}

		void set_max_res_norm(double value){
			_max_res_norm = value;
		}

	

		Eigen::MatrixXd _R;
		Eigen::VectorXd _z;
		Eigen::VectorXd _x;


        int  _qr_solver_thread=1;
	private:
		double _max_res_norm = 3.0;
		double _vtpv_pre = 0.0;
		int _freedom_pre = 0;

		unique_ptr<great::t_gqr_srif> qr_solve;
	};

} // namespace

#endif

