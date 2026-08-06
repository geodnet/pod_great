
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>

#include "gutils/gcommon.h"


using namespace std;

namespace gnut {

// COMMON TRACE FUNCTION

	Eigen::MatrixXd  newmat2eigen(Matrix k) {
		Eigen::MatrixXd ans = Eigen::MatrixXd::Zero(k.Nrows(), k.Ncols());
		for (int i = 0; i < k.Nrows(); i++) {
			for (int j = 0; j < k.Ncols(); j++) {
				ans(i, j) = k(i + 1, j + 1);
			}
		}
		return ans;
	}

	Eigen::VectorXd newmat2eigen(ColumnVector k)
	{
		Eigen::VectorXd ans = Eigen::VectorXd::Zero(k.Nrows());
		for (int i = 0; i < k.Nrows(); i++) {
			ans(i) = k(i + 1);
		}
		return ans;
	}

	ColumnVector eigen2newmat(Eigen::VectorXd k)
	{
		ColumnVector ans(k.rows());
		for (int i = 0; i < k.rows(); i++) {
			ans(i + 1) = k(i);
		}
		return ans;
	}

	Matrix eigen2newmat(Eigen::MatrixXd k)
	{
		Matrix ans(k.rows(),k.cols());
		for (int i =0;i<k.rows();i++){
			for (int j =0;j<k.cols();j++){
				ans(i+1,j+1) = k(i,j);
			}
		}
		return ans;
	}

	size_t t_gpair_string_hash::operator()(const pair<string, string>& a) const
	{
		return std::hash<string>()(a.first+a.second);
	}

    bool _check_amb_depend(vector<double> &pdE, vector<double> &pdC, int iNamb, int iN_oneway, vector<int> arriIpt2ow, int& iNdef) {
        double dC_dot = 0.0;
        const double dOper[4] = { 1.0, -1.0, -1.0, 1.0 };
        double dEPS = 1e-12;

        // outrange
        for (int j = 0; j < iN_oneway; ++j) {
            if (arriIpt2ow[j] > iNamb) return true;
        }

        // schmit-gram
        for (int i = 0; i < iNdef; ++i) {
            pdC[i] = 0.0;
            for (int j = 0; j < iN_oneway; ++j) {
                pdC[i] += pdE[i * iNamb + arriIpt2ow[j] - 1] * dOper[j];
            }
            dC_dot = dC_dot + pdC[i] * pdC[i];
        }
        bool isLdepend = false;
        // change 1e-12 to 1e-8, jqwu
        if (abs(dC_dot - iN_oneway) <= 1e-8) return true;

        // add iNdef vector
        for (int j = 0; j < iNamb; ++j) pdE[iNdef * iNamb + j] = 0.0;
        for (int j = 0; j < iN_oneway; ++j) pdE[iNdef * iNamb + arriIpt2ow[j] - 1] = dOper[j];

        if (dC_dot > dEPS) {
            for (int i = 0; i < iNdef; ++i) {
                if (abs(pdC[i]) < dEPS) continue;
                for (int j = 0; j < iNamb; ++j)
                    pdE[iNdef * iNamb + j] -= pdC[i] * pdE[i * iNamb + j];
            }
            //Orthogonalization
            dC_dot = 0.0;
            for (int j = 0; j < iNamb; ++j)
                dC_dot += pdE[iNdef * iNamb + j] * pdE[iNdef * iNamb + j];
            dC_dot = sqrt(dC_dot);
            for (int j = 0; j < iNamb; ++j)
                pdE[iNdef * iNamb + j] /= dC_dot;
        }
        else {
            // Here means the new vector is completely independent from previous vector
            dC_dot = sqrt(iN_oneway * 1.0);
            for (int j = 0; j < iN_oneway; ++j)
                pdE[iNdef * iNamb + arriIpt2ow[j] - 1] /= dC_dot;
        }
        pdC[iNdef] = dC_dot;
        iNdef = iNdef + 1;
        return false;
    }

} // namespace
