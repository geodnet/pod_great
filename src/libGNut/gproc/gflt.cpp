
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  This file is part of the G-Nut C++ library.

-*/

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <queue>

#include "gproc/gflt.h"
#include "gutils/gmatrixconv.h"
#include "gproc/gqr_srif.h"
#include "gproc/gqr_eigen.h"
#include "gutils/gtiming.h"
#include "gutils/gconst.h"
#include <memory>



using namespace std;
using namespace std::chrono;

namespace gnut {


// Constructors
t_gflt::t_gflt()
{
}

 
t_SRIF::t_SRIF()
{
	qr_solve = make_unique<great::t_gqr_eigen>();
}

t_SRIF::t_SRIF(const t_SRIF& Others):
	_R(Others._R),
	_z(Others._z),
	_x(Others._x),
	_max_res_norm(Others._max_res_norm),
	_vtpv_pre(Others._vtpv_pre),
	_freedom_pre(Others._freedom_pre)
{
    this->set_qr_solver(Others._qr_solver_thread);
}

t_SRIF t_SRIF::operator=(const t_SRIF &Others) {
    _R = Others._R;
    _z = Others._z;
    _x = Others._x;
    _max_res_norm = Others._max_res_norm;
    _vtpv_pre = Others._vtpv_pre;
    _freedom_pre = Others._freedom_pre;

    this->set_qr_solver(Others._qr_solver_thread);

    return *this;
}
   
// Destructors
t_gflt::~t_gflt()
{
}

   
t_SRIF::~t_SRIF()
{
}
   
void t_SRIF::set_qr_solver(int num_threads)
{
	qr_solve = make_unique<great::t_gqr_eigen>();
	_qr_solver_thread = num_threads;
}

void t_SRIF::init(const Eigen::MatrixXd& R, const Eigen::VectorXd& x0) {
	_R = R; // P0 = R0.i() * R0.i().t();
	_z = R * x0;
	_x = x0;
	_vtpv_pre = 0.0;
	_freedom_pre = 0;
}

Eigen::VectorXd t_SRIF::measurement_update(great::t_glsqEquationMatrix measurement_equ)
{

	Eigen::MatrixXd A;
	Eigen::VectorXd z;
	measurement_equ.changeNomralizedEigen(A, z, _R.cols());

	return this->measurement_update(A, z);
}


Eigen::VectorXd t_SRIF::measurement_update(const Eigen::MatrixXd &A, const Eigen::VectorXd &z)
{
	int nObs = A.rows();
	int nPar = _x.rows();


	if (nObs == 0) {
		return Eigen::VectorXd::Zero(0);
	}

	Eigen::MatrixXd A_e(nObs + nPar, nPar);
	Eigen::VectorXd z_e(nObs + nPar);

	A_e.topRows(nPar) = this->_R;
	A_e.bottomRows(nObs) = A;


	z_e.topRows(nPar) = this->_z;
	z_e.bottomRows(nObs) = z;




	
	qr_solve->M_update(A_e,z_e);
	;



	this->_R = A_e.block(0,0,nPar,nPar).triangularView<Eigen::Upper>();
	this->_z.resize(nPar);
	this->_z = z_e.topRows(nPar);
	this->_x.noalias() = this->_R.triangularView<Eigen::Upper>().solve(this->_z);



    Eigen::VectorXd e = z_e.bottomRows(z_e.rows()-nPar);
    double std_sigma0 = sqrt( (e.transpose()*e)(0,0) / nObs );
    this->_sigma0 = std_sigma0;


    return z_e;
}



void reorder_matrix(vector<int> move_id,Eigen::MatrixXd& R);




void t_SRIF::time_update_use_removeinfo(great::t_gupdateparinfo remove_info)
{
	vector<int> remove_id;
	vector<t_gpar> newparlist;
	vector<t_gpar> equparlist;
	great::t_glsqEquationMatrix virtual_equ;
	remove_info.get(remove_id, newparlist, equparlist, virtual_equ);

	// First add virutal equ
	if (equparlist.size() > 0) {

		Eigen::MatrixXd A;
		Eigen::VectorXd z;
		virtual_equ.changeNomralizedEigen(A,z,_R.cols()+equparlist.size());

		this->measurement_update_with_removepar(A,z,remove_id);
	}
	else{
		// Here Since the equpar size == 0 therefore the Matrix is totally a triangular matrix for SRIF Measurement update
		this->qr_solve->remove_par_size = 0;
		// Second Remove Pars
		this->rm_param(remove_id);
	}


		// Third add New
		{
			Eigen::VectorXd vecR(newparlist.size());
			for (int i = 0; i < newparlist.size(); i++) {
				vecR(i) = 1 / newparlist[i].apriori();
			}
			this->add_parameter(vecR.asDiagonal(), Eigen::VectorXd::Zero(newparlist.size()));

		}

}

void t_SRIF::add_parameter(const Eigen::MatrixXd& R, const Eigen::VectorXd& x0) {
	//UpperTriangularMatrix R = Cholesky(P0).t().i(); // P0 = R0.i() * R0.i().t();
	int nPar_add = x0.rows();
	if (nPar_add == 0) {
		return;
	}
	Eigen::VectorXd z = R * x0;

	int nPar_orig = _R.rows();
	int nPar = nPar_orig + nPar_add;

	Eigen::MatrixXd newR = Eigen::MatrixXd::Zero(nPar, nPar);
	newR.block(0, 0, nPar_orig, nPar_orig) = this->_R;
	newR.block(nPar_orig, nPar_orig, nPar_add, nPar_add) = R;
	this->_R = newR;


	Eigen::VectorXd newz = Eigen::VectorXd::Zero(nPar);
	newz.topRows(nPar_orig) = this->_z;
	newz.bottomRows(nPar_add) = z;
	this->_z = newz;

	Eigen::VectorXd newx = Eigen::VectorXd::Zero(nPar);
	newx.topRows(nPar_orig) = this->_x;
	newx.bottomRows(nPar_add) = x0;
	this->_x = newx;
}


void t_SRIF::measurement_update_with_removepar(Eigen::MatrixXd A,Eigen::VectorXd z,vector<int> remove_id)
{
	int nPar_old = _x.rows();
	int nPar_new = A.cols() - nPar_old;

	// Construct A and e
	Eigen::MatrixXd A_e = Eigen::MatrixXd::Zero(nPar_old + nPar_new, nPar_old + nPar_new);
	Eigen::VectorXd z_e = Eigen::VectorXd::Zero(nPar_old + nPar_new);

		A_e.block(0,0,nPar_old,nPar_old) = this->_R;
		A_e.block(nPar_old,0,nPar_new,nPar_old+nPar_new)  =  A;
		z_e.topRows(nPar_old) = this->_z;
		z_e.bottomRows(nPar_new) = z;
	int nPar_rm = remove_id.size();
	sort(remove_id.begin(),remove_id.end());

	int nPar_orig = A_e.cols();
	int nRows = A_e.rows();
	int nPar_keep = nPar_orig - nPar_rm;

	Eigen::MatrixXd newR = Eigen::MatrixXd::Zero(nRows,nPar_orig);
	//Eigen::VectorXd newz = Eigen::VectorXd::Zero(z_e.cols);



	int j = 0;
	int k = nPar_rm ;
	for (int i = 0; i <  nPar_rm; i++) {
		int ncol = (remove_id[i]-1) - j;
		if (ncol > 0) {
			// "don't remove parameters"
			newR.block(0,k,nRows,ncol) = A_e.block(0,j,nRows,ncol);
			k = k + ncol;
		}
		// begin with remove parameters
		newR.block(0,i,nRows,1) = A_e.block(0,remove_id[i]-1,nRows,1);
		j = remove_id[i];
	}
	int ncol = nPar_orig - j;
	if (ncol > 0) {
		newR.block(0,k,nRows,ncol) = A_e.block(0,j,nRows,ncol);
	}



	// QR solve
	qr_solve->T_update(newR,z_e,nPar_rm);

	this->_R = newR.block(nPar_rm,nPar_rm,nPar_keep,nPar_keep);
	this->_z = (z_e.topRows(nPar_orig)).bottomRows(nPar_keep);
	this->_x = Eigen::VectorXd::Zero(nPar_keep);


	return;
}

void reorder_matrix(vector<int> move_id,Eigen::MatrixXd& R)
{
	// the paramter of move_id list should be move to the head part
	int nPar_rm = move_id.size();
	if (nPar_rm == 0){
		return;
	}
	sort(move_id.begin(),move_id.end());

	int nPar_orig = R.cols();
	int nPar_new = nPar_orig - nPar_rm;

	Eigen::MatrixXd newR = Eigen::MatrixXd::Zero(nPar_orig,nPar_orig);

	int j = 0;
	int k = nPar_rm;

	for (int i = 0; i <  nPar_rm; i++) {
		int ncol = (move_id[i]-1) - j;
		if (ncol > 0) {
			// "don't remove parameters"
			newR.block(0,k,nPar_orig,ncol) = R.block(0,j,nPar_orig,ncol);
			k = k + ncol;
		}
		newR.block(0,i,nPar_orig,1) = R.block(0,move_id[i]-1,nPar_orig,1);
		j = move_id[i];
	}
	int ncol = nPar_orig - j;
	if (ncol > 0) {
		newR.block(0,k,nPar_orig,ncol) = R.block(0,j,nPar_orig,ncol);
	}
	R = newR;
}

void t_SRIF::rm_param(vector<int> remove_id) {

	int nPar_rm = remove_id.size();
	if (nPar_rm == 0) {
		return;
	}

	sort(remove_id.begin(),remove_id.end());

	int nPar_orig = _R.cols();
	int nPar_new = nPar_orig - nPar_rm;


	Eigen::MatrixXd newR = Eigen::MatrixXd::Zero(nPar_orig,nPar_orig);
	Eigen::VectorXd newx = Eigen::VectorXd::Zero(nPar_new);
	


	int j = 0;
	int k = nPar_rm ;
	for (int i = 0; i <  nPar_rm; i++) {
		int ncol = (remove_id[i]-1) - j;
		if (ncol > 0) {
			newR.block(0,k,nPar_orig,ncol) = this->_R.block(0,j,nPar_orig,ncol);
			newx.segment(k-nPar_rm,ncol) = this->_x.segment(j,ncol);
			k = k + ncol;
		}
		newR.block(0,i,nPar_orig,1) = this->_R.block(0,remove_id[i]-1,nPar_orig,1);
		j = remove_id[i];
	}
	int ncol = nPar_orig - j;
	if (ncol > 0) {
		newR.block(0,k,nPar_orig,ncol) = this->_R.block(0,j,nPar_orig,ncol);
		newx.segment(k-nPar_rm,ncol) = this->_x.segment(j,ncol);
	}



	Eigen::HouseholderQR<Eigen::MatrixXd> qr(newR);
	Eigen::MatrixXd R_e = qr.matrixQR().triangularView<Eigen::Upper>();
	this->_z = qr.householderQ().transpose() * this->_z;



	this->_R = R_e.block(nPar_rm,nPar_rm, nPar_new, nPar_new).triangularView<Eigen::Upper>();
	this->_z = Eigen::VectorXd(this->_z.bottomRows(nPar_new));

	this->_x = newx;
}

    set<int> t_SRIF::obs_quality_control(const Eigen::MatrixXd &A, const Eigen::MatrixXd &z, const Eigen::VectorXd &z_e) {

        //int nObs = A.rows();
        int nPar = _x.rows();
        int nObs = z_e.rows() - nPar;

        // record idx for big biases, from idx 0
        set<int> detect_idx;
        map<int,double> detect_weight;
        Eigen::VectorXd e = z_e.bottomRows(z_e.rows()-nPar);

        Eigen::VectorXd v = A * this->_x - z;

        double vtpv_pre = (v.transpose()*v)(0,0);
        double etpe_pre = (e.transpose()*e)(0,0);

        double std_sigma0 = sqrt( (e.transpose()*e)(0,0) / nObs );
        double max_v = 0.0;

        int max_detect_size =  1;
        while (true) {

            // sort vtpv
            auto compare_func = [](const pair<double, int>& k1, const pair<double, int >& k2) -> bool {
                return abs(k1.first) > abs(k2.first);
            };
            priority_queue<pair<double, int>, vector<pair<double, int> >, decltype(compare_func)> sort_idx(compare_func);
            for (int i =0;i<v.rows();i++){
                sort_idx.push(make_pair(v(i), i));
                if (sort_idx.size() > max_detect_size) {
                    sort_idx.pop();
                }
            }
            double Max_v = abs(sort_idx.top().first);
            double Mean_v = v.mean();

            vector<int> new_bad_arr;

            // find bias
            while (!sort_idx.empty()) {
                if (!detect_idx.count(sort_idx.top().second) || detect_weight[sort_idx.top().second]!=0.0) {
                    detect_weight[sort_idx.top().second] = 0.0;
                    if (!detect_idx.count(sort_idx.top().second)){
                        new_bad_arr.push_back(sort_idx.top().second);
                        detect_idx.insert(sort_idx.top().second);
                    }
                }
                sort_idx.pop();
            }
            if (new_bad_arr.size() == 0){
                cout << "Duplicate gross observations detected" << endl;
                break;
            }

            Eigen::MatrixXd cm = Eigen::MatrixXd::Zero(nPar + nObs,detect_idx.size());
            Eigen::MatrixXd bad_weight = Eigen::MatrixXd::Zero(detect_idx.size(),detect_idx.size());
            int count_detect=0;
            for (auto idx: detect_idx){
                bad_weight(count_detect,count_detect)=detect_weight[idx];
                cm(nPar+idx,count_detect) = 1;
                count_detect++;
            }

            chrono::high_resolution_clock::time_point end_t1, beg_t1;
            beg_t1 = chrono::high_resolution_clock::now();
            Eigen::MatrixXd Sm = Eigen::MatrixXd::Zero(nPar+nObs,detect_idx.size());

            Sm = cm;

            qr_solve->Apply_M_trans(Sm);
          
            Eigen::MatrixXd SmObs = Sm.bottomRows(Sm.rows()-nPar);

  
            Eigen::VectorXd dz = (SmObs.transpose()*SmObs+bad_weight).inverse()*(SmObs.transpose()*(-(z_e.bottomRows(z_e.rows()-nPar))));

    
            Eigen::VectorXd new_z,new_x,new_e;
            new_z = z_e.topRows(nPar)+Sm.topRows(nPar)*dz;
            new_x = this->_R.triangularView<Eigen::Upper>().solve(new_z);

            new_e = z_e.bottomRows(z_e.rows() - nPar) + SmObs * dz;

            v = A * new_x - (z+(cm*dz).bottomRows(nObs));
            end_t1 = chrono::high_resolution_clock::now();

            double new_vtpv = (v.transpose()*v)(0,0);
            double new_etpe = (new_e.transpose()*new_e)(0,0);
            double new_std_sigma0 = sqrt( (new_e.transpose() * new_e)(0,0) / (nObs - detect_idx.size()) );

            if (etpe_pre - new_etpe > _chi2inv[new_bad_arr.size()]) {
                //if (true){
                std_sigma0 = new_std_sigma0;
                vtpv_pre = new_vtpv;
                etpe_pre = new_etpe;

                // update SRIF z and x
                this->_z = new_z;
                this->_x = new_x;
                //break;
            }
            else{
                for (auto& bad_idx : new_bad_arr) {
                    detect_idx.erase(bad_idx);
                }
                break;
            }

        }
        this->_sigma0 = std_sigma0;


        return detect_idx;
    }

    void t_SRIF::obs_quality_control(const Eigen::MatrixXd &A, const Eigen::MatrixXd &z, const Eigen::VectorXd &z_e,
                                      set<int> given_idx) {

        if (given_idx.empty()){
            return;
        }
        //int nObs = A.rows();
        int nPar = _x.rows();
        int nObs = z_e.rows() - nPar;

        // record idx for big biases, from idx 0
        set<int> detect_idx = given_idx;
        Eigen::VectorXd e = z_e.bottomRows(z_e.rows()-nPar);

        Eigen::VectorXd v = A * this->_x - z;

        double vtpv_pre = (v.transpose()*v)(0,0);
        double etpe_pre = (e.transpose()*e)(0,0);

        double std_sigma0 = sqrt( (e.transpose()*e)(0,0) / nObs );
        double max_v = 0.0;

        Eigen::MatrixXd cm = Eigen::MatrixXd::Zero(nPar + nObs,detect_idx.size());
        Eigen::MatrixXd bad_weight = Eigen::MatrixXd::Zero(detect_idx.size(),detect_idx.size());
        int count_detect=0;
        for (auto idx: detect_idx){
            cm(nPar+idx,count_detect) = 1;
            count_detect++;
        }

        chrono::high_resolution_clock::time_point end_t1, beg_t1;
        beg_t1 = chrono::high_resolution_clock::now();
        Eigen::MatrixXd Sm = Eigen::MatrixXd::Zero(nPar+nObs,detect_idx.size());

        Sm = cm;

        qr_solve->Apply_M_trans(Sm);

        Eigen::MatrixXd SmObs = Sm.bottomRows(Sm.rows()-nPar);

        // dz represent the biase parameter estimated value
        Eigen::VectorXd dz = (SmObs.transpose()*SmObs+bad_weight).inverse()*(SmObs.transpose()*(-(z_e.bottomRows(z_e.rows()-nPar))));


        // DEBUG test recompute v
        Eigen::VectorXd new_z,new_x,new_e;
        new_z = z_e.topRows(nPar)+Sm.topRows(nPar)*dz;
        new_x = this->_R.triangularView<Eigen::Upper>().solve(new_z);

        new_e = z_e.bottomRows(z_e.rows() - nPar) + SmObs * dz;

        v = A * new_x - (z+(cm*dz).bottomRows(nObs));
        end_t1 = chrono::high_resolution_clock::now();

        double new_vtpv = (v.transpose()*v)(0,0);
        double new_etpe = (new_e.transpose()*new_e)(0,0);
        double new_std_sigma0 = sqrt( (new_e.transpose() * new_e)(0,0) / (nObs - detect_idx.size()) );

        std_sigma0 = new_std_sigma0;

        // update SRIF z and x
        this->_z = new_z;
        this->_x = new_x;

        this->_sigma0 = std_sigma0;

    }


} // namespace
