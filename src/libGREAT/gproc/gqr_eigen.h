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

#ifndef GQR_EIGEN_H
#define GQR_EIGEN_H

#include <memory>
#include "gproc/gqr_srif.h"

namespace great
{
    class t_gqr_eigen : public t_gqr_srif
    {
    public:
        void M_update(Eigen::MatrixXd &A, Eigen::VectorXd &l) override;
        void Apply_M_trans(Eigen::MatrixXd &cm) override;
        void T_update(Eigen::MatrixXd &A, Eigen::VectorXd &l, const int &new_num) override;

    private:
        std::unique_ptr<Eigen::HouseholderQR<Eigen::MatrixXd>> eigen_qr_solver;
    };

}
#endif // !GQR_EIGEN_H
