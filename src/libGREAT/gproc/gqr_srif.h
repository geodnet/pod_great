/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2023-09-25  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gqr_srif.h
* @brief		header files of QR factorization
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2023-09-25
*
*/

#ifndef GQR_SRIF_H
#define GQR_SRIF_H

#include <vector>
#include "Eigen/Dense"

#include "gproc/glsqmatrix.h"

namespace great
{
	class t_gqr_srif
	{
	public:

		virtual ~t_gqr_srif(){};
		virtual void M_update(Eigen::MatrixXd &A, Eigen::VectorXd &l) = 0;
		virtual void Apply_M_trans(Eigen::MatrixXd &cm) = 0;
		virtual void T_update(Eigen::MatrixXd &A, Eigen::VectorXd &l, const int &new_num) = 0;

		// information used in dataflow_opt
		std::vector<std::pair<int,int> > update_list; // update list element: <parsize , parnum >
		int remove_par_size;
		int remove_par_before;
		t_glsqEquationMatrix obs_equ;
		set<int> clk_id;
		t_gallpar par_list;

		// information used in dataflow_ambopt
		int phase_num;
		int amb_num;
	};

}
#endif // !GQR_SRIF_H
