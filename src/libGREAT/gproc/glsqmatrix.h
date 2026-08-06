/**
*
* @verbatim
	History
	 -1.0 ZhengHJ       2018-12-02  created
	 -1.1 Jiande Huang  2020-10-11  Optimize program efficiency
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    glsqmatrix.h
* @brief		Matrix for NEQ,W,observ_equations
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/

#ifndef  GLSQMATRIX_H
#define  GLSQMATRIX_H

#include "ExportHeaders/ExportBase.h"
#include <gset/gsetproc.h>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "newmat/newmat.h"
#include "gutils/gtypeconv.h"
#include "gutils/gturboedit.h"
#include "gmodels/gbasemodel.h"

using namespace std;
using namespace gnut;
namespace great
{
	typedef vector<pair<int, double> > SimpleB;
	/**
	* @brief  Matrix for storage observe equations of B,P,L
	*/
	class BASE_LIBRARY_EXPORT t_glsqEquationMatrix :public t_gbaseEquation
	{
	public:

		/** @brief default constructor */
		t_glsqEquationMatrix();

		/** @brief default destructor  */
		~t_glsqEquationMatrix();

		/**
		* @brief add equations
		* @param[in] B_value coeff of observ equation
		* @param[in] P_value weight of observ equation
		* @param[in] l_value res of observ equaion
		*/
		void add_equ(const vector<pair<int, double> >& B_value, const double& P_value, const double& l_value, const string& stie_name, const string& sat_name, const t_gobscombtype& obscombtype, const bool& is_newamb);
		void add_multi_equ(const vector<SimpleB> & B_value, const SymmetricMatrix& P, vector<double> l_value, string stie_name = "", string sat_name = "", t_gobscombtype obscombtype = t_gobscombtype(), bool is_newamb = false);
		void add_equ(const t_glsqEquationMatrix& Other);
		/**
		* @brief add equations by newmat format
		* @param[in] B_value coeff of observ equation newmat format
		* @param[in] P_value weight of observ equation newmat format
		* @param[in] l_value res of observ equaion newmat format
		*/
		void add_equ(const Matrix& B_value, const DiagonalMatrix& P_value, const ColumnVector& l_value, const vector<string>& site_name, const vector<string>&  sat_name, const vector<t_gobscombtype>& obstype);


		void set_newamb(const int& idx, const bool& is_newamb);
		void set_newamb(t_gturboedit* tb_slip);

		void set_zeroApriori(int idx){
			this->P[idx] = 0.0;
		}
		/**
		* @brief remove specified equations
		* @note idx from 1
		* param[in] idx specified loc
		*/
		void remove(const int& idx);


		void remove_last_equation();

		/**
		* @brief change equations to newmat format
		* @param[out] B_value coeff of observ equation newmat format
		* @param[out] P_value weight of observ equation newmat format
		* @param[out] l_value res of observ equaion newmat format
		*/
		void chageNewMat(Matrix& B_value, DiagonalMatrix& P_value, ColumnVector& l_value, const int& par_num);
		void chageNewMat(Matrix& B_value, SymmetricMatrix& P_value, ColumnVector& l_value, const int& par_num);

		void changeSparseEigen(Eigen::SparseMatrix<double>& SparseB,Eigen::VectorXd& l_value);
		void changeNomralizedEigen(Eigen::MatrixXd& B, Eigen::VectorXd& l, const int& par_num);
        void setP(Eigen::MatrixXd inputP);
	
		/** @brief print B P L Matrix*/
		void print();

		/**
		* @brief get numbers of equations
		* @return size of equations
		*/
		int num_equ()const;
		int phase_equ() const;

		/**
		* @brief get resiuals of equations
		* @return LTPL of equations
		*/
		double res_equ() const;

		/**
		* @brief get resiuals of equations
		* @param[in] phase phase only or both
		* @return LTPL of equations
		*/
		double res_equ(bool phase) const; //add by xiongyun

		string get_sitename(int equ_idx) const;
		set<string> get_satlist(string rec) const;
		set<string> get_sitelist(string sat) const;
        map<pair<string,string>,int> get_site_sat_table(bool is_phase) const;

		string get_satname(int equ_idx) const;

		string get_obscombtype2str(int equ_idx) const;

		t_gobscombtype get_obscombtype(int equ_idx) const;

		vector<double> get_codeomc() const;
		vector<double> get_phaseomc() const;

		bool is_newamb(int equ_idx) const;

		friend BASE_LIBRARY_EXPORT void  print_equ_debinfo(const t_glsqEquationMatrix& equ);


		vector<int> find_equ(const string& site);
		vector<int> find_equ(const string& site, const string & sat);
		int find_equ(const string& site, const string& sat, const t_gobscombtype& obscomtype);

		int get_equ_num(const string& sat,bool onlyPhase = true);

		void clear_allequ();
		void  swap_allequ();
		void reorder_equation2LP(bool P_first=false);

		void move_equ_par(int begin_size,int add_size);

        SymmetricMatrix _P;
	protected:

		void _change_symmP();

		vector<pair<string, string> >  _site_sat_pairlist;
		vector<t_gobscombtype> _obstypelist;
		vector<bool> _newamb_list;
	};


	/**
	* @brief  Matrix for storage Matrix NEQ
	*/
	class BASE_LIBRARY_EXPORT t_glsqSymmetricMatrix
	{
	public:


		/**
		* @brief default destructor
		* @return size of NEQ matrix
		*/
		virtual int num() const = 0;

		virtual double& num(int a, int b) = 0;
		//double  num(int a, int b) const ;

		/**
		* @brief resize the NEQ matrix according to specified size
		* @note clean the data before
		* @param[in] size specified size
		*/
		virtual void resize(int size) = 0;

		/**
		* @brief add observ equations
		* @parma[in] equ observ equations
		*/
		virtual void add(const t_glsqEquationMatrix& equ) = 0;

		/** @brief add one zero dimension at last*/
		virtual void addBackZero() = 0;

		/**
		* @brief removed specified dimension in the NEQ matrix
		* @note idx from 1
		* @param[in] idx specified dimension
		*/
		virtual void remove(int idx) = 0;

		/**
		* @brief print NEQ matrix in std out
		*/
		virtual void print() = 0;

		/**
		* @brief get diagonal value in NEQ matrix by specified loc
		* @note idx from 1
		* @param[in] idx specified loc
		*/
		virtual double center_value(int idx) const = 0;

		/**
		* @brief change NEQ matrix to NewMat format
		* @return NEQ matrix in newmat format
		*/
		virtual SymmetricMatrix changeNewMat() const = 0;

	};


	/**
	* @brief class for store BTPL(W) matrix
	*/
	class BASE_LIBRARY_EXPORT t_glsqColumnVector
	{
	public:

		/**
		* @brief get size of column vector
		* @return size
		*/
		virtual int num() const = 0;

		virtual double& num(int a) = 0;
		//double num(int a) const;

		/**
		* @brief add observ equation
		* @param[in] equ observ equation
		*/
		virtual void add(const t_glsqEquationMatrix& equ) = 0;

		/**
		* @brief add zero dimension in the last loc
		*/
		virtual void addBackZero() = 0;

		/**
		* @brief resize the column vector
		* @note resize will clean data before
		* @param[in] size specified size of column vector
		*/
		virtual void resize(int size) = 0;

		/**
		* @brief remove specified dimension of W matrix
		* @param[in] idx specified loc
		*/
		virtual void remove(int idx) = 0;

		/** @brief print W matrix in stdout */
		virtual void print() = 0;

		virtual ColumnVector changeNewMat() = 0;

	};


	class BASE_LIBRARY_EXPORT V_SymmetricMatrix;

	class BASE_LIBRARY_EXPORT V_ColumnVector :public t_glsqColumnVector
	{
	public:

		int num() const;
		double& num(int a);
		void resize(int num);
		void add(const t_glsqEquationMatrix& equ);
		void del(const t_glsqEquationMatrix& equ);
		void add(const t_glsqEquationMatrix& equ, bool phase); 
		void addBackZero();
		void remove(int idx);
		void print();
		ColumnVector changeNewMat();
		Eigen::VectorXd changeEigen();
		void swap(int a, int b);

		friend bool BASE_LIBRARY_EXPORT remove_lsqmatrix(int idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W);
		friend bool BASE_LIBRARY_EXPORT rearrange_lsqmatrix(const vector<int>& remove_idx, t_gallpar& allpar, V_SymmetricMatrix& NEQ, V_ColumnVector& W, Eigen::MatrixXd& N11, Eigen::MatrixXd& N21, Eigen::MatrixXd& N22, Eigen::VectorXd& W1, Eigen::VectorXd& W2, bool idx_from_zero);
		friend int BASE_LIBRARY_EXPORT rearrange_lsqmatrix(vector<int>& remove_idx, V_SymmetricMatrix & NEQ, V_ColumnVector & W, t_gallpar& allpar);
		friend class t_glsq;
	private:
		vector<double> _element;
	};


	class BASE_LIBRARY_EXPORT V_SymmetricMatrix :public t_glsqSymmetricMatrix
	{
	public:

		int num() const;
		inline double& num(int a, int b) { return (a < b) ? _element[b - 1][a - 1] : _element[a - 1][b - 1]; }
		void resize(int num);
		void add(const t_glsqEquationMatrix& equ);
		void del(const t_glsqEquationMatrix& equ); 
		void add(const t_glsqEquationMatrix& equ, bool phase);
		void addBackZero();
		void remove(int idx);
		void print();
		double center_value(int idx) const;
		SymmetricMatrix changeNewMat() const;
		Eigen::MatrixXd changeEigen();
		void swap(int a, int b);

		friend bool BASE_LIBRARY_EXPORT remove_lsqmatrix(int idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W);
		friend bool BASE_LIBRARY_EXPORT rearrange_lsqmatrix(const vector<int>& remove_idx, t_gallpar& allpar, V_SymmetricMatrix& NEQ, V_ColumnVector& W, Eigen::MatrixXd& N11, Eigen::MatrixXd& N21, Eigen::MatrixXd& N22, Eigen::VectorXd& W1, Eigen::VectorXd& W2, bool idx_from_zero);
		friend int BASE_LIBRARY_EXPORT rearrange_lsqmatrix(vector<int>& remove_idx, V_SymmetricMatrix & NEQ, V_ColumnVector & W, t_gallpar& allpar);
		friend class t_glsq;

	private:
		vector< vector<double> > _element;
	};




	bool BASE_LIBRARY_EXPORT remove_lsqmatrix(int idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W);


	// by matrix part
	bool BASE_LIBRARY_EXPORT rearrange_lsqmatrix(const vector<int>& remove_idx, t_gallpar& allpar, V_SymmetricMatrix& NEQ, V_ColumnVector& W, Eigen::MatrixXd& N11, Eigen::MatrixXd& N21, Eigen::MatrixXd& N22, Eigen::VectorXd& W1, Eigen::VectorXd& W2, bool idx_from_zero = true);
	int BASE_LIBRARY_EXPORT rearrange_lsqmatrix(vector<int>& remove_idx, V_SymmetricMatrix & NEQ, V_ColumnVector & W,t_gallpar& allpar);
}


#endif /*  GLSQMAT_H  */