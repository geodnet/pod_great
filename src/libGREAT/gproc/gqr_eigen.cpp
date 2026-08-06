#include "gproc/gqr_eigen.h"
#include "gutils/gtiming.h"

namespace great
{
    void t_gqr_eigen::M_update(Eigen::MatrixXd &A, Eigen::VectorXd &l)
    {
        eigen_qr_solver = std::make_unique<Eigen::HouseholderQR<Eigen::MatrixXd>>(A);
        A = eigen_qr_solver->matrixQR().triangularView<Eigen::Upper>();
        l = eigen_qr_solver->householderQ().transpose() * l;
    }

    void t_gqr_eigen::Apply_M_trans(Eigen::MatrixXd &cm)
    {
        cm = eigen_qr_solver->householderQ().transpose() * cm;
    }

    void t_gqr_eigen::T_update(Eigen::MatrixXd &A, Eigen::VectorXd &l, const int &new_num)
    {

        int col_num = A.cols();
        Eigen::MatrixXd AL = A.leftCols(new_num);
        Eigen::MatrixXd AB = A.rightCols(col_num - new_num);


        Eigen::HouseholderQR<Eigen::MatrixXd> qr(AL);
        A.leftCols(new_num) = qr.matrixQR().triangularView<Eigen::Upper>();
        A.rightCols(col_num - new_num) = qr.householderQ().transpose() * AB;

        l = qr.householderQ().transpose() * l;
    }

}