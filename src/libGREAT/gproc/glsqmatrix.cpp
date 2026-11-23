/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2018-12-02  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    glsqmatrix.cpp
* @brief		Matrix for NEQ,W,observ_equations
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/

#include "gproc/glsqmatrix.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <stdexcept>
#include <algorithm>
#if defined(__x86_64__) || defined(_M_X64)
#include <emmintrin.h>
#endif
#include <assert.h>
#include <chrono>



namespace great
{

	t_glsqEquationMatrix::t_glsqEquationMatrix()
	{
		_P.resize(0);	
	}

	t_glsqEquationMatrix::~t_glsqEquationMatrix()
	{

	}



	void t_glsqEquationMatrix::add_equ(const vector<pair<int, double> >& B_value, const double& P_value, const double& l_value, const string& stie_name, const string& sat_name, const t_gobscombtype& obscombtype, const bool& is_newamb)
	{
		if (this->_P.Nrows() > 0 ){
			B.push_back(B_value);
			l.push_back(l_value);
			this->_site_sat_pairlist.push_back(make_pair(stie_name, sat_name));
			this->_obstypelist.push_back(obscombtype);
			this->_newamb_list.push_back(is_newamb);
			_P.resize_keep(_P.Nrows()+1);	
			_P(_P.Nrows(),_P.Nrows()) = P_value;
		}
		else{
			B.push_back(B_value);
			P.push_back(P_value);
			l.push_back(l_value);
			this->_site_sat_pairlist.push_back(make_pair(stie_name, sat_name));
			this->_obstypelist.push_back(obscombtype);
			this->_newamb_list.push_back(is_newamb);
		}
	}

	void t_glsqEquationMatrix::add_multi_equ(const vector<SimpleB> & B_value, const SymmetricMatrix& P, vector<double> l_value, string site_name, string sat_name , t_gobscombtype obscombtype , bool is_newamb)
	{
		this->_change_symmP();
		B.insert(B.end(), B_value.begin(), B_value.end());
		_P.resize_keep(_P.Nrows() + P.Nrows());
		_P.SubMatrix(_P.Nrows() - P.Nrows() + 1, _P.Nrows(), _P.Nrows() - P.Nrows() + 1, _P.Nrows()) = P;
		l.insert(l.end(), l_value.begin(), l_value.end());

        for (int i =0;i<B_value.size();i++){
            _site_sat_pairlist.push_back(make_pair(site_name,sat_name));
        }
	}

	void t_glsqEquationMatrix::add_equ(const t_glsqEquationMatrix & Other)
	{
		if (Other.B.size() == 0 ){
			return;
		}
		if (Other._P.Nrows() > 0) {
			this->_change_symmP();
			B.insert(B.end(), Other.B.begin(), Other.B.end());
			_P.resize_keep(_P.Nrows() + Other._P.Nrows());
			_P.SubMatrix(_P.Nrows() - Other._P.Nrows() + 1, _P.Nrows(), _P.Nrows() - Other._P.Nrows() + 1, _P.Nrows()) = Other._P;
			l.insert(l.end(), Other.l.begin(), Other.l.end());
            _site_sat_pairlist.insert(_site_sat_pairlist.end(),Other._site_sat_pairlist.begin(),Other._site_sat_pairlist.end());
		}
		else {
			for (int i = 0; i < Other.num_equ(); i++)
			{
				this->add_equ(Other.B[i], Other.P[i], Other.l[i], Other._site_sat_pairlist[i].first, Other._site_sat_pairlist[i].second, Other._obstypelist[i],false);
			}
		}
	}

	void t_glsqEquationMatrix::_change_symmP()
	{
		if (this->_P.Nrows()>0)	{
			return;
		}
		if (B.empty()){
			return;
		}
		this->_P.resize(B.size());
		this->_P=0.0;
		for (int i =0;i<_P.Nrows();i++){
			this->_P(i+1,i+1) = this->P[i];	
		}
		return;	
	}


	void t_glsqEquationMatrix::add_equ(const Matrix& B_value, const DiagonalMatrix& P_value, const ColumnVector& l_value, const vector<string>& site_name, const vector<string>&  sat_name, const vector<t_gobscombtype>& obstype)
	{

		if (B_value.Nrows() != P_value.Nrows() || B_value.Nrows() != l_value.Nrows() ||
			B_value.Nrows() != site_name.size() || B_value.Nrows() != sat_name.size() || B_value.Nrows() != obstype.size())
		{
			throw exception();
		}

		int num_equ = B_value.Nrows();
		//for B
		vector<pair<int, double> > temp;
		for (int row = 1; row <= B_value.Nrows(); row++)
		{
			temp.clear();
			for (int col = 1; col <= B_value.Ncols(); col++)
			{
				if (B_value(row, col) == 0.0)
				{
					continue;
				}
				temp.push_back(make_pair(col, B_value(row, col)));
			}

			B.push_back(temp);


			P.push_back(P_value(row));
			l.push_back(l_value(row));

			// Note: row from 1
			_site_sat_pairlist.push_back(make_pair(site_name[row - 1], sat_name[row - 1]));
			_obstypelist.push_back(obstype[row - 1]);
			_newamb_list.push_back(false);
		}
	}

	void t_glsqEquationMatrix::set_newamb(const int& idx, const bool& is_newamb)
	{
		if (idx >= this->num_equ())
		{
			throw exception(logic_error("input idx is more than number of equ"));
		}

		_newamb_list[idx] = is_newamb;
	}

	void t_glsqEquationMatrix::set_newamb(t_gturboedit* tb_slip)
	{
		if (!tb_slip) {
			return;
		}

		set<pair<string, string> > site_sat_list;
		for (int i = 0; i < _site_sat_pairlist.size(); i++) {
			if (!this->get_obscombtype(i).is_phase()) {
				continue;
			}
			string rec = this->get_sitename(i);
			string sat = this->get_satname(i);

			if (tb_slip->new_amb(rec, sat)) {
				_newamb_list[i] = true;
				site_sat_list.insert(make_pair(rec, sat));
			}
		}
		for (auto iter : site_sat_list) {
			tb_slip->set_new_amb(iter.first, iter.second, false);
		}
	}

	void t_glsqEquationMatrix::remove(const int& idx)
	{
		B.erase(B.begin() + idx - 1);
		P.erase(P.begin() + idx - 1);
		l.erase(l.begin() + idx - 1);

		_site_sat_pairlist.erase(_site_sat_pairlist.begin() + idx - 1);
		_obstypelist.erase(_obstypelist.begin() + idx - 1);
		_newamb_list.erase(_newamb_list.begin() + idx - 1);
	}

	void t_glsqEquationMatrix::remove_last_equation()
	{
		B.pop_back();
		P.pop_back();
		l.pop_back();
		_site_sat_pairlist.pop_back();
		_obstypelist.pop_back();
		_newamb_list.pop_back();
	}


	void t_glsqEquationMatrix::chageNewMat(Matrix& B_value, DiagonalMatrix& P_value, ColumnVector& l_value, const int& par_num)
	{
		B_value.ReSize(B.size(), par_num); B_value = 0.0;
		for (int row = 0; row < B.size(); row++)
		{
			for (int col = 0; col < B[row].size(); col++)
			{
				B_value(row + 1, B[row][col].first) = B[row][col].second;
			}
		}
		P_value.ReSize(P.size()); P_value = 0.0;
		for (int row = 0; row < P.size(); row++) {
			P_value(row + 1) = P[row];
		}
		l_value.ReSize(l.size()); l_value = 0.0;
		for (int row = 0; row < l.size(); row++) {
			l_value(row + 1) = l[row];
		}
	}

	void t_glsqEquationMatrix::chageNewMat(Matrix & B_value, SymmetricMatrix & P_value, ColumnVector & l_value, const int& par_num)
	{
		B_value.ReSize(B.size(), par_num); B_value = 0.0;
		for (int row = 0; row < B.size(); row++)
		{
			for (int col = 0; col < B[row].size(); col++)
			{
				B_value(row + 1, B[row][col].first) = B[row][col].second;
			}
		}
		if (_P.Nrows() > 0) {
			P_value = _P;
		}
		else {
			P_value.ReSize(P.size()); P_value = 0.0;
			for (int row = 0; row < P.size(); row++) {
				P_value(row + 1, row + 1) = P[row];
			}
		}
		l_value.ReSize(l.size()); l_value = 0.0;
		for (int row = 0; row < l.size(); row++) {
			l_value(row + 1) = l[row];
		}
	}

	void t_glsqEquationMatrix::changeSparseEigen(Eigen::SparseMatrix<double>& SparseB,Eigen::VectorXd& l_value)
	{
		vector<Eigen::Triplet<double>> coeffB;
		for (int row = 0; row < this->B.size(); row++) {
			double sqrtP = sqrt(this->P[row]);
			for (int col = 0; col < this->B[row].size(); col++) {
				coeffB.push_back(Eigen::Triplet<double>(row, this->B[row][col].first - 1, this->B[row][col].second*sqrtP));
			}

			l_value(row) = l[row]*sqrtP;
		}
		SparseB.setFromTriplets(coeffB.begin(), coeffB.end());	
	}

	void t_glsqEquationMatrix::changeNomralizedEigen(Eigen::MatrixXd & B_value, Eigen::VectorXd & l_value, const int & par_num)
	{
		B_value = Eigen::MatrixXd::Zero(this->B.size(), par_num);
		l_value = Eigen::VectorXd::Zero(this->B.size());

		Eigen::SparseMatrix<double> SparseB(this->B.size(), par_num);
		vector<Eigen::Triplet<double>> coeffB;
		for (int row = 0; row < B.size(); row++) {
			for (int col = 0; col < B[row].size(); col++) {
				coeffB.push_back(Eigen::Triplet<double>(row, B[row][col].first - 1, B[row][col].second));
			}
		}
		SparseB.setFromTriplets(coeffB.begin(), coeffB.end());

		for (int row = 0; row < l.size(); row++) {
			l_value(row) = l[row];
		}

		
		if (_P.Nrows() > 0) {
			Eigen::MatrixXd SymmP(_P.Nrows(), _P.Nrows());
			// Lower Symm Matrix;
			for (int row = 0; row < _P.Nrows(); row++) {
				for (int col = 0; col <= row; col++) {
					SymmP(row, col) = _P(row + 1, col + 1);
				}
			}
			Eigen::MatrixXd RP = SymmP.selfadjointView<Eigen::Lower>().llt().matrixL().transpose();
			B_value = RP * SparseB;
			l_value = RP * l_value;
		}
		else {
			Eigen::VectorXd vecP(P.size());
			for (int row = 0; row < P.size(); row++) {
				vecP(row) = sqrt(P[row]);
			}
			B_value = vecP.asDiagonal()*SparseB;
			l_value = vecP.asDiagonal()*l_value;
		}
	}

	void t_glsqEquationMatrix::print()
	{
		int count = 0;
		std::cout << "B is" << endl;
		for (int i = 0; i < B.size(); i++)
		{
			count = 0;
			std::cout << setw(6) << _site_sat_pairlist[i].first
				<< setw(6) << _site_sat_pairlist[i].second
				<< endl;
			for (int j = 0; j < B[i].size(); j++)
			{
				count++;
				std::cout << scientific << setprecision(15) << setw(30) << B[i][j].second << "   ";
				if (count % 3 == 0)
					std::cout << endl;
			}
			std::cout << endl;

		}
		std::cout << "P is" << endl;
		for (int i = 0; i < B.size(); i++)
		{
			std::cout << setw(6) << _site_sat_pairlist[i].first
				<< setw(6) << _site_sat_pairlist[i].second
				<< endl;
			std::cout << scientific << setprecision(15) << setw(30) << P[i] << endl;
		}
		std::cout << "l is" << endl;
		for (int i = 0; i < B.size(); i++)
		{
			std::cout << setw(6) << _site_sat_pairlist[i].first
				<< setw(6) << _site_sat_pairlist[i].second
				<< endl;
			std::cout << scientific << setprecision(15) << setw(30) << l[i] << endl;
		}
	}

	int t_glsqEquationMatrix::num_equ() const
	{
		return B.size();
	}

	int t_glsqEquationMatrix::phase_equ() const
	{
		int result_phase_num = 0;
		for (int i =0;i< _obstypelist.size();i++){
			if (_obstypelist[i].is_phase()){
				result_phase_num++;
			}
		}
		return result_phase_num;
	}

	double t_glsqEquationMatrix::res_equ() const
	{
		double ans = 0;
		for (int i = 0; i < num_equ(); i++) {
			ans += P[i] * l[i] * l[i];
		}
		return ans;
	}

	double t_glsqEquationMatrix::res_equ(bool phase) const
	{
		if (!phase)
		{
			return this->res_equ();
		}
		double ans = 0;
		for (int i = 0; i < num_equ(); i++) {
			if (get_obscombtype(i).is_code()) {
				continue;
			}
			ans += P[i] * l[i] * l[i];
		}
		return ans;
	}

	string t_glsqEquationMatrix::get_sitename(int equ_idx) const
	{
		return _site_sat_pairlist[equ_idx].first;
	}

	set<string> t_glsqEquationMatrix::get_satlist(string rec) const
	{
		set<string> sat_list;
		int num_equ = this->num_equ();
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i].first == rec) {
				sat_list.insert(_site_sat_pairlist[i].second);
			}
		}

		return sat_list;
	}

	set<string> t_glsqEquationMatrix::get_sitelist(string sat) const
	{
		set<string> site_list;
		int num_equ = this->num_equ();
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i].second == sat) {
				site_list.insert(_site_sat_pairlist[i].first);
			}
		}

		return site_list;
	}

	string t_glsqEquationMatrix::get_satname(int equ_idx) const
	{
		return _site_sat_pairlist[equ_idx].second;
	}

	string t_glsqEquationMatrix::get_obscombtype2str(int equ_idx) const
	{
		return _obstypelist[equ_idx].convert2str();
	}

	t_gobscombtype t_glsqEquationMatrix::get_obscombtype(int equ_idx) const
	{
		return _obstypelist[equ_idx];
	}

	vector<double> t_glsqEquationMatrix::get_codeomc() const
	{
		int equ_num = this->num_equ();
		vector<double> codeomc;
		for (int i = 0; i < equ_num; i++) 
		{
			if (_obstypelist[i].is_code()) 
			{
				codeomc.push_back(l[i]);
			}
		}
		return codeomc;
	}
	vector<double> t_glsqEquationMatrix::get_phaseomc() const
	{
		int equ_num = this->num_equ();
		vector<double> phaseomc;
		for (int i = 0; i < equ_num; i++)
		{
			if (_obstypelist[i].is_phase())
			{
				phaseomc.push_back(l[i]);
			}
		}
		return phaseomc;
	}

	

	bool t_glsqEquationMatrix::is_newamb(int equ_idx) const
	{
		return _newamb_list[equ_idx];
	}
	void print_equ_debinfo(const t_glsqEquationMatrix& equ)
	{
		for (int i = 0; i < equ.num_equ() / 2.0; i++)
		{
			cout << setw(5) << equ._site_sat_pairlist[2 * i].first
				<< setw(5) << equ._site_sat_pairlist[2 * i].second
				<< setiosflags(ios::fixed) << setprecision(4)
				<< setw(13) << equ.l[2 * i]
				<< setw(13) << equ.l[2 * i + 1]
				<< endl;
		}
	}

	vector<int> t_glsqEquationMatrix::find_equ(const string& site)
	{
		vector<int> ans;
		int num_equ = this->num_equ();
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i].first == site) {
				ans.push_back(i);
			}
		}
		return ans;
	}

	int t_glsqEquationMatrix::get_equ_num(const string& sat,bool onlyPhase)
	{
		int ans = 0;
		int num_equ = this->num_equ();
		for (int i=0;i<num_equ;i++){
			if (_site_sat_pairlist[i].second == sat ){
                if (!onlyPhase || _obstypelist[i].is_phase() ){
                    ans++;
                }
			}
		}
		return ans;
	}

	vector<int> t_glsqEquationMatrix::find_equ(const string& site, const string & sat)
	{
		vector<int> ans;
		int num_equ = this->num_equ();
		auto site_sat = make_pair(site, sat);
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i] == site_sat) {
				ans.push_back(i);
			}
		}
		return ans;
	}

	int t_glsqEquationMatrix::find_equ(const string& site, const string& sat, const t_gobscombtype& obscomtype)
	{
		int ans = -1;
		int num_equ = this->num_equ();
		auto site_sat = make_pair(site, sat);
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i] == site_sat && _obstypelist[i] == obscomtype) {
				ans = i;
				break;
			}
		}
		return ans;
	}

	void t_glsqEquationMatrix::clear_allequ()
	{
		B.clear();
		P.clear();
		l.clear();
		_site_sat_pairlist.clear();
		_obstypelist.clear();
		_newamb_list.clear();
	}

	void t_glsqEquationMatrix::reorder_equation2LP(bool P_first){
		if (_P.Nrows() > 0	){
			throw std::runtime_error("P matrix of observation euations should be diagnol");
		}
		vector<int> code_idx;
		vector<int> phase_idx;
		for (int i =0;i<_obstypelist.size();i++){
			if (_obstypelist[i].is_code()){
				code_idx.push_back(i);	
			}
			else{
				phase_idx.push_back(i);
			}
		}

		// add phase observation
			
		vector<pair<string, string> > new_site_sat_pairlist;
		vector<t_gobscombtype> new_obstypelist;
		vector<bool> new_newamb_list;
		vector<vector<pair<int, double>>> new_B;
        vector<double> new_P;
        vector<double> new_l;
		vector<int> first_idx = P_first?code_idx:phase_idx;
		vector<int> second_idx = P_first?phase_idx:code_idx;
		for (auto idx: first_idx){
			new_site_sat_pairlist.push_back(_site_sat_pairlist[idx]);
			new_obstypelist.push_back(_obstypelist[idx]);
			new_newamb_list.push_back(_newamb_list[idx]);
			new_B.push_back(B[idx]);
			new_P.push_back(P[idx]);
			new_l.push_back(l[idx]);
		}
		for (auto idx: second_idx){
			new_site_sat_pairlist.push_back(_site_sat_pairlist[idx]);
			new_obstypelist.push_back(_obstypelist[idx]);
			new_newamb_list.push_back(_newamb_list[idx]);
			new_B.push_back(B[idx]);
			new_P.push_back(P[idx]);
			new_l.push_back(l[idx]);
		}

		_site_sat_pairlist = new_site_sat_pairlist;
		_obstypelist = new_obstypelist;
		_newamb_list = new_newamb_list;
		B = new_B;
		P = new_P;
		l = new_l;
	}

	void t_glsqEquationMatrix::move_equ_par(int begin_size,int add_size)
	{
		int num_equ = this->num_equ();
		for (int i =0;i<num_equ;i++){
			int coeff_size = this->B[i].size();
			for (int j=0;j<coeff_size;j++){
				if (this->B[i][j].first > begin_size){
					this->B[i][j].first+=add_size;
				}
			}
		}
	}

    map<pair<string,string>,int> t_glsqEquationMatrix::get_site_sat_table(bool is_phase) const {
        map<pair<string,string>,int> result;
        int num_equ = this->num_equ();
        for (int i =0;i<num_equ;i++){
            if (is_phase && _obstypelist[i].is_code()){
                continue;
            }
            if (!is_phase && _obstypelist[i].is_phase()){
                continue;
            }
            result[_site_sat_pairlist[i]] = i;
        }
        return result;
    }

    void t_glsqEquationMatrix::setP(Eigen::MatrixXd inputP) {

        _P.resize(inputP.rows());

        for (int i =0;i<inputP.rows();i++){
            for (int j =0;j<=i;j++){
                _P(i+1,j+1) = inputP(i,j);
            }
        }

    }



	int V_SymmetricMatrix::num() const
	{
		return _element.size();
	}



	void V_SymmetricMatrix::resize(int num)
	{
		_element.clear();
		_element.reserve(num);
		for (int i = 0; i < num; i++)
		{
			_element.push_back(vector<double>(i + 1, 0.0));
		}
	}

	void V_SymmetricMatrix::add(const t_glsqEquationMatrix & equ)
	{
		// cycle equ
		for (int num = 0; num < equ.num_equ(); num++)
		{
			// cycle par
			auto B_temp = equ.B[num];
			sort(B_temp.begin(), B_temp.end());
			for (int ipar = 0; ipar < B_temp.size(); ipar++)
			{
				auto row = B_temp[ipar];
				for (int jpar = 0; jpar <= ipar; jpar++)
				{
					auto col = B_temp[jpar];
					_element[row.first - 1][col.first - 1] += row.second * equ.P[num] * col.second;
				}
			}
		}
	}

	void V_SymmetricMatrix::del(const t_glsqEquationMatrix& equ)
	{
		// cycle equ
		for (int num = 0; num < equ.num_equ(); num++)
		{
			// cycle par
			auto B_temp = equ.B[num];
			sort(B_temp.begin(), B_temp.end());
			for (int ipar = 0; ipar < B_temp.size(); ipar++)
			{
				auto row = B_temp[ipar];
				for (int jpar = 0; jpar <= ipar; jpar++)
				{
					auto col = B_temp[jpar];
					_element[row.first - 1][col.first - 1] -= row.second * equ.P[num] * col.second;
				}
			}
		}
	}

	void V_SymmetricMatrix::add(const t_glsqEquationMatrix& equ, bool phase)
	{
		if (!phase)
		{
			this->add(equ);
			return;
		}
		// cycle equ
		for (int num = 0; num < equ.num_equ(); num++)
		{
			if (equ.get_obscombtype(num).is_code()) {
				continue;
			}
			// cycle par 
			auto B_temp = equ.B[num];
			sort(B_temp.begin(), B_temp.end());
			for (int ipar = 0; ipar < B_temp.size(); ipar++)
			{
				auto row = B_temp[ipar];
				for (int jpar = 0; jpar <= ipar; jpar++)
				{
					auto col = B_temp[jpar];
					_element[row.first - 1][col.first - 1] += row.second * equ.P[num] * col.second;
				}
			}
		}
	}

	void V_SymmetricMatrix::addBackZero()
	{
		_element.push_back(vector<double>(_element.size() + 1, 0.0));
	}

	void V_SymmetricMatrix::remove(int idx)
	{
		for (int row = _element.size(); row > idx; row--) 
		{
			_element[row - 1].erase(_element[row - 1].begin() + idx - 1);
		}
		_element.erase(_element.begin() + idx - 1);
	}

	void V_SymmetricMatrix::print()
	{
		cout << setw(20) << setprecision(5);
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			for (auto Col = Row->begin(); Col != Row->end(); Col++)
			{
				cout << setw(20) << *Col;
			}
			cout << endl;
		}
	}


	SymmetricMatrix V_SymmetricMatrix::changeNewMat() const
	{
		SymmetricMatrix temp(_element.size());
		int row = 1, col = 1;
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			col = 1;
			for (auto Col = Row->begin(); Col != Row->end(); Col++)
			{
				temp(row, col) = *Col;
				col++;
			}
			row++;
		}
		return temp;
	}

	Eigen::MatrixXd V_SymmetricMatrix::changeEigen()
	{
		int row = this->num();
		Eigen::MatrixXd result = Eigen::MatrixXd::Zero(row,row);
		for (int i =0;i<row;i++){
			for (int j=0;j<row;j++){
				result(i,j) = this->num(i+1,j+1);
			}
		}
		return result;
	}

	double V_SymmetricMatrix::center_value(int idx) const
	{
		return _element[idx - 1][idx - 1];
	}

	void V_SymmetricMatrix::swap(int a, int b)
	{
		double temp_a = num(a, a);
		double temp_b = num(b, b);
		double temp_ab = num(a, b);
		for (int i = 0; i < num(); i++)
		{
			double temp = num(a, i + 1);
			num(a, i + 1) = num(b, i + 1);
			num(b, i + 1) = temp;
		}
		num(a, a) = temp_b;
		num(b, b) = temp_a;
		num(a, b) = temp_ab;
	}

	int V_ColumnVector::num() const
	{
		return _element.size();
	}

	double & V_ColumnVector::num(int a)
	{
		return _element[a - 1];
	}

	void V_ColumnVector::resize(int num)
	{
		_element.resize(num, 0.0);
	}

	void V_ColumnVector::add(const t_glsqEquationMatrix & equ)
	{
		for (int num = 0; num < equ.num_equ(); num++)
		{
			for (int ipar = 0; ipar < equ.B[num].size(); ipar++)
			{
				_element[equ.B[num][ipar].first - 1] += equ.B[num][ipar].second * equ.P[num] * equ.l[num];
			}
		}
	}

	void V_ColumnVector::del(const t_glsqEquationMatrix& equ)
	{
		for (int num = 0; num < equ.num_equ(); num++)
		{
			for (int ipar = 0; ipar < equ.B[num].size(); ipar++)
			{
				_element[equ.B[num][ipar].first - 1] -= equ.B[num][ipar].second * equ.P[num] * equ.l[num];
			}
		}
	}

	void V_ColumnVector::add(const t_glsqEquationMatrix& equ, bool phase)
	{
		if (!phase)
		{
			this->add(equ);
			return;
		}
		for (int num = 0; num < equ.num_equ(); num++)
		{
			if (equ.get_obscombtype(num).is_code()) {
				continue;
			}
			for (int ipar = 0; ipar < equ.B[num].size(); ipar++)
			{
				_element[equ.B[num][ipar].first - 1] += equ.B[num][ipar].second * equ.P[num] * equ.l[num];
			}
		}
	}

	void V_ColumnVector::addBackZero()
	{
		_element.push_back(0.0);
	}

	void V_ColumnVector::remove(int idx)
	{
		_element.erase(_element.begin() + idx - 1);
	}

	void V_ColumnVector::print()
	{
		cout << setprecision(5);
		for (auto row = _element.begin(); row != _element.end(); row++)
		{
			cout << setw(20) << *row << endl;
		}
	}

	ColumnVector V_ColumnVector::changeNewMat()
	{
		ColumnVector temp(_element.size());
		int row = 1;
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			temp(row) = *Row;
			row++;
		}
		return temp;
	}
	Eigen::VectorXd V_ColumnVector::changeEigen()
	{
		return Eigen::Map<Eigen::VectorXd>(_element.data(),_element.size());
	}

	void V_ColumnVector::swap(int a, int b)
	{
		double temp;
		temp = _element[a - 1];
		_element[a - 1] = _element[b - 1];
		_element[b - 1] = temp;
	}

	bool remove_lsqmatrix(int idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W)
	{
		if (NEQ.num() != W.num())
			return false;
		if (idx <1 || idx > NEQ.num())
			return false;

		double center_value = NEQ.center_value(idx);

		if (!double_eq(center_value, 0.0))
		{
			int NEQ_num = NEQ.num();

			// record need remove element
			vector<double>  remove_element;
			remove_element.reserve(NEQ_num);
			
			double w_remove = 0.0;
			vector<int> remove_rows;
			for (int col = 1; col <= NEQ_num; col++)
			{
				remove_element.push_back(NEQ.num(idx, col));
				if (col == idx || double_eq(NEQ.num(idx, col), 0.0)) 
				{
					continue;
				}
				remove_rows.push_back(col - 1);
			}

			w_remove = W._element[idx - 1];
			
			// remove the par in NEQ and W by for cycle
			int remove_rows_size = remove_rows.size();

#ifdef USE_OPENMP
#pragma	omp parallel for schedule(dynamic)
			for (int i = 0; i < remove_rows_size; ++i)
			{
				int row = remove_rows[i];
				double temp_coeff = -remove_element[row] / center_value;
				for (int col = 0; col <= row; col++)
				{
					NEQ._element[row][col] += remove_element[col] * temp_coeff;
				}

				W._element[row] += w_remove * temp_coeff;
			}
#else
			for (int i = 0; i < remove_rows_size; ++i)
			{
				int row = remove_rows[i];
				double temp_coeff = -remove_element[row] / center_value;
				vector<double> &row_element=NEQ._element[row];
				for (int col = 0; col <= row; col++)
				{
					//NEQ._element[row][col] += remove_element[col] * temp_coeff;
					row_element[col] += remove_element[col] * temp_coeff;
				}

				W._element[row] += w_remove * temp_coeff;
			}
#endif
		}

		NEQ.remove(idx);
		W.remove(idx);

		return true;
	}
	bool BASE_LIBRARY_EXPORT rearrange_lsqmatrix(const vector<int>& remove_idx, t_gallpar& allpar, V_SymmetricMatrix& NEQ, V_ColumnVector& W, Eigen::MatrixXd& N11, Eigen::MatrixXd& N21, Eigen::MatrixXd& N22, Eigen::VectorXd& W1, Eigen::VectorXd& W2, bool idx_from_zero)
	{	

		auto beg_t = chrono::high_resolution_clock::now();

		int matrix_begin = idx_from_zero ? 0 : 1;
		if (remove_idx.size() == 0) {
			return true;
		}
		int rows = NEQ.num();
		int rows_new = rows - remove_idx.size();
		vector<bool> idx_map(rows, false);
		set<int> zero_idx;
		for (int i = 0; i < remove_idx.size(); i++)
		{
			idx_map[remove_idx[i] - 1] = true;
			if (NEQ.center_value(remove_idx[i]) == 0.0){
				zero_idx.insert(remove_idx[i]-1);
			}
		}

		int n22 = remove_idx.size() - zero_idx.size();

		// init size
		N11 = Eigen::MatrixXd::Zero(rows_new, rows_new);
		N21 = Eigen::MatrixXd::Zero(n22, rows_new);
		N22 = Eigen::MatrixXd::Zero(n22, n22);

		W1 = Eigen::VectorXd::Zero(rows_new);
		W2 = Eigen::VectorXd::Zero(n22);

		auto end_t = chrono::high_resolution_clock::now();
		cout << "Init Matrix Spend time is " << chrono::duration_cast<chrono::milliseconds>(end_t-beg_t).count() / 1000.0 << " sec " << endl;

	
		// rearrange
		int i_skip = 0;
		int i_not_skip = 0;
		for (int i = 0; i < rows; i++)
		{
			if (zero_idx.count(i)!=0){
				continue;
			}
			int j_skip = 0;
			int j_not_skip = 0;
			if (idx_map[i] == true)
			{
				i_skip++;
				W2(i_skip - 1 + matrix_begin) = W._element[i];
				for (int j = 0; j <= i; j++)
				{
					if(zero_idx.count(j)!=0){
						continue;
					}
					if (idx_map[j] == true)
					{
						j_skip++;
						N22(i_skip - 1 + matrix_begin, j_skip - 1 + matrix_begin) = NEQ._element[i][j];
						N22(j_skip - 1 + matrix_begin, i_skip - 1 + matrix_begin) = NEQ._element[i][j];
					}
					else
					{
						j_not_skip++;
						N21(i_skip - 1 + matrix_begin, j_not_skip - 1 + matrix_begin) = NEQ._element[i][j];
					}
				}
			}
			else
			{
				i_not_skip++;
				W1(i_not_skip - 1 + matrix_begin) = W._element[i];
				for (int j = 0; j <= i; j++)
				{
					if (zero_idx.count(j)!=0){
						continue;
					}

					if (idx_map[j] == true)
					{
						j_skip++;
						N21(j_skip - 1 + matrix_begin, i_not_skip - 1+ matrix_begin) = NEQ._element[i][j];
					}
					else
					{
						j_not_skip++;
						N11(i_not_skip - 1 + matrix_begin, j_not_skip - 1 + matrix_begin) = NEQ._element[i][j];
						N11(j_not_skip - 1 + matrix_begin, i_not_skip - 1 + matrix_begin) = NEQ._element[i][j];
					}
				}
			}
		}
		return true;
	}
	int BASE_LIBRARY_EXPORT rearrange_lsqmatrix(vector<int>& remove_idx, V_SymmetricMatrix & NEQ, V_ColumnVector & W, t_gallpar& allpar)
	{
		if (remove_idx.size() == 0) 
		{
			return true;
		}

		int rows = NEQ.num();
		int rows_new = rows - remove_idx.size();


		vector<bool> idx_map(rows, false);
		set<int> zero_idx;
		int remove_size = remove_idx.size();
		for (int i = 0; i < remove_size; i++)
		{
			idx_map[remove_idx[i] - 1] = true;
			if (NEQ.center_value(remove_idx[i])==0.0){
				zero_idx.insert(remove_idx[i]-1);
			}
		}

		// rearrange paramter and init map_idx
		t_gallpar allpar_orig = allpar; allpar.delAllParam();
		vector<int> map_idx(rows);//map from old -> new
		vector<int> remove_idx_new; // from new->old (size is remove size)
		// add not remove parameter
		int count=0;
		for (int i = 0; i < rows; i++) 
		{
			if (idx_map[i]) 
			{
				continue;
			}
			allpar.addParam(allpar_orig[i]);
			map_idx[i]=count++;
		}
		// add remove paramter (but center value not zero)
		for (int i = 0; i < remove_size; i++) 
		{
			if (zero_idx.count(remove_idx[i]-1)!=0){
				continue;
			}
			allpar.addParam(allpar_orig[remove_idx[i]-1]);
			map_idx[remove_idx[i]-1]=count++;
			remove_idx_new.push_back(remove_idx[i]);
		}
		// add remove parameter (but center value is zero)
		for (auto i:zero_idx){
			allpar.addParam(allpar_orig[i]);
			map_idx[i] = count++;
			remove_idx_new.push_back(i+1);
		}

		assert(rows == count);

		// rearrange remove_idx
		assert(remove_idx.size() == remove_idx_new.size());
		remove_idx = remove_idx_new;


		V_SymmetricMatrix NEQ_orig = NEQ;
		V_ColumnVector W_orig = W;
		
		// rearrange

#ifdef USE_OPENMP
#pragma	omp parallel for schedule(dynamic)
#endif
		for (int i =0;i<rows;i++){
			int x = map_idx[i];
			for (int j=0;j<=i;j++){
				int y = map_idx[j];
				(x>y?NEQ._element[x][y]:NEQ._element[y][x]) = NEQ_orig._element[i][j];
				//NEQ.num(map_idx[i]+1,map_idx[j]+1) = NEQ_orig._element[i][j];
			}
			W._element[x]=W_orig._element[i];
		}
		return zero_idx.size();
	}
	

}