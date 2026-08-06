/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2018-12-02  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallprocdata.cpp
* @brief		header files of least squares estimator
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/

#include "gproc/glsq.h"
#include "gutils/ginfolog.h"
#include "gio/gfile.h"
#include "gutils/gturboedit.h"
#include <sstream>
#include <stdlib.h>
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "gproc/ginverse_Eigen.h"
#include <chrono>


// test clk
#include <stack>

//#define DEBUG_COMPARE

using namespace std;

namespace great
{

	// Constructors

	t_glsq::t_glsq(t_gsetbase* set) :
		_log(nullptr),
		_tempfile(0),
		_obs_total_num(0),
		_npar_tot_num(0),//
		_sigma0(0.0), _res_obs(0.0)
	{


		_dx.ReSize(0);

		if (set != nullptr)
		{
			_beg = dynamic_cast<t_gsetgen*>(set)->beg();
			_end = dynamic_cast<t_gsetgen*>(set)->end();
			_interval = dynamic_cast<t_gsetgen*>(set)->sampling(); 
			_buffer_size = dynamic_cast<t_gsetproc*>(set)->lsq_buffer_size() * 1024 * 1000;
		}

		// get random_tempfile
		stringstream this_addr;
		this_addr << this;
		reset_tempfile("tempfile_" + this_addr.str());

		// default using Eigen Cholosky method
		_solve_matrix = shared_ptr<t_ginverse<Matrix> >(new t_ginverse_Eigen_Cholosky<Matrix>(false));

	}

	t_glsq::~t_glsq()
	{

		if (_tempfile)
		{
			string name = _tempfile->name();
			if (_tempfile->is_open())
			{
				_tempfile->close();
			}
			delete _tempfile;
			_tempfile = nullptr;
		
			int remove_err = remove(name.c_str());
		}

	}

	t_glsq::t_glsq(const t_glsq& Other)
		:_epo(Other._epo),
		_beg(Other._beg),
		_end(Other._end),
		_tempfile(nullptr),
		_log(Other._log),
		_NEQ(Other._NEQ),

		_W(Other._W),
		_Qx(Other._Qx),
		_update_lsqpar(Other._update_lsqpar),
		_x_solve(Other._x_solve),
		_stdx(Other._stdx),
		_res_obs(Other._res_obs),
		_npar_tot_num(Other._npar_tot_num),// ztd clk change 
		_interval(Other._interval),   
		_vtpv(Other._vtpv), _sigma0(Other._sigma0), _obs_total_num(Other._obs_total_num),
		_solve_matrix(Other._solve_matrix),
		_buffer_size(Other._buffer_size)
	{
		//TODO:TOO BAD solve way NOW..
		stringstream this_addr;
		this_addr << this;
		reset_tempfile("tempfile_" + this_addr.str());
	}

    void t_glsq::clear_all() {

        _x_solve.delAllParam();
        _NEQ.resize(0);
        _W.resize(0);
        _npar_tot_num = 0;
    }

	int t_glsq::add_parameter(const t_gpar& par,bool head)
	{
		// add the par
		if (head){
			_x_solve.addParamHead(par);
		}
		else{
			_x_solve.addParam(par);
		}
		_NEQ.addBackZero();
		_W.addBackZero();
		_npar_tot_num++;
		return 1;
	}

	void t_glsq::add_par_state_equ(t_gpar::t_type par_type, int order, double dt, double noise)
	{
		_update_lsqpar->set_par_state_mode(par_type, order, dt, noise);
	}

	void t_glsq::set_update_par(shared_ptr<t_gupdatepar> updatepar)
	{
		_update_lsqpar = updatepar;
	}

	int t_glsq::add_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch,bool write_temp)
	{
		if (write_temp) 
		{
			this->write_equation(equ, epoch);
		}

		_epo = epoch;

		_NEQ.add(equ);
		_W.add(equ);

		_obs_total_num += equ.num_equ();
		_res_obs       += equ.res_equ();

		return 1;
	}

	bool  t_glsq::write_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch)
	{
		if (_tempfile == nullptr)
		{
			if (_log) _log->comment(1, "t_glsq::write_equation", "the tempfile don't exist!");
			throw exception();
		}

		int len_size = 0;

		// write coefficinet Matrix
		for (int Row = 1; Row <= equ.num_equ(); Row++)
		{

			// length of record
			len_size = 0;

			// write identifier
			_tempfile->write("obs", 3);
			len_size += 3;

			// write parameter number
			int par_num = equ.B[Row - 1].size();
			_tempfile->write((char*)&par_num, SIZE_INT);
			len_size += SIZE_INT;

			// write time
			_tempfile->write((char*)&epoch, sizeof(t_gtime));
			len_size += sizeof(t_gtime);

			// write ambflag
			int is_newamb = equ.is_newamb(Row - 1) ? 1 : 0;
			_tempfile->write((char*)&is_newamb, SIZE_INT);
			len_size += SIZE_INT;

			// write info
			string out_info;
			out_info = equ.get_sitename(Row - 1) + "    "
				+ equ.get_satname(Row - 1) + "     "
				+ equ.get_obscombtype2str(Row - 1) + "     ";
			int len_out_info = out_info.size();
			_tempfile->write((char*)&len_out_info, SIZE_INT);
			len_size += SIZE_INT;
			_tempfile->write(out_info.c_str(),len_out_info);
			len_size += len_out_info;

			// write coefficient
			for (int Col = 1; Col <= par_num; Col++) {
				double coeff = equ.B[Row - 1][Col - 1].second;
				int loc = equ.B[Row - 1][Col - 1].first;
			

				_tempfile->write((char*)&loc, SIZE_INT);
				len_size += SIZE_INT;
				_tempfile->write((char*)&coeff, SIZE_DBL);
				len_size += SIZE_DBL;

			}

			// Write Max loc 
			int Max_loc = _x_solve.parNumber();
			_tempfile->write((char*)&Max_loc, SIZE_INT);
			len_size += SIZE_INT;

			// write P
			double p = equ.P[Row - 1];
			_tempfile->write((char*)&p, SIZE_DBL);
			len_size += SIZE_DBL;

			// write res
			double res = equ.l[Row - 1];
			_tempfile->write((char*)&res, SIZE_DBL);
			len_size += SIZE_DBL;

			_tempfile->write((char*)&len_size, SIZE_INT);
		}

		return true;
	}


	// remove parameter ! idx from 1
	int t_glsq::remove_parameter(const int& idx)
	{
		
		add_apriori_weight(idx);
		double center_value = _NEQ.center_value(idx);
		double w_remove = _W.num(idx);
		if (center_value!=0.0){
			// compute the ltpl
			_res_obs -= w_remove * w_remove / center_value;

		}
		// write the coeff to tempfile 
		_write_coefficient(idx);
		// remove NEQ W Matrix
		if (!remove_lsqmatrix(idx, _NEQ, _W))
		{
			throw "remove lsq matrix error!";
		}

		// remove the par
		_x_solve.delParam(idx - 1);
		// _npar_tot_num need to plus 1, otherwise the lsq->sigma may be wrong.
		_npar_tot_num--;

		return 1;
	}

	void t_glsq::add_apriori_weight(const int& idx)
	{
		// add apriori weight
		double apriori = _x_solve[idx - 1].apriori();
		if (!double_eq(apriori, 0.0)) 
		{
			_NEQ.num(idx, idx) = _NEQ.num(idx, idx) + 1 / (apriori*apriori);
		}
		return;
	}

	void t_glsq::reset_tempfile(string filename)
	{
		if (_tempfile)
		{
			if (_tempfile->is_open())
			{
				_tempfile->close();
			}
			remove(_tempfile->name().c_str());
			delete _tempfile;
		}

		//_tempfile = new t_giof;
		_tempfile = new t_giobigf(filename,_buffer_size);
		_tempfile->tsys(t_gtime::GPS);
		_tempfile->mask(filename);
		_tempfile->append(false);
		_tempfile->open(filename, ios::out | ios::trunc | ios::binary);

		// for begin
		int identify = -1;
		_tempfile->write((char*)&identify, SIZE_INT);
	}


	int t_glsq::_write_coefficient(int idx)
	{
		ostringstream os;
		if (idx < 1 || idx > _NEQ.num()) 
		{
			if (_log)
			{
				_log->comment(1, "t_glsq::_write_coefficient", "the input idx is Wrong!");
			}
			return -1;
		}
		
		if (!_tempfile)
		{
			if (_log)
			{
				_log->comment(1, "t_glsq::_write_coefficient", "don't init the tempfile");
			}
			return -1;

		}
		int len_size = 0;
		// write identifier
		_tempfile->write("par", 3);
		len_size += 3;

		// write parameter number
		vector<pair<int, double> > par_record;

		
		double w_record = _W.num(idx);
		int NEQ_num = _NEQ.num();
		for (int col = 1; col <= NEQ_num; col++)
		{
			double value = _NEQ.num(idx, col);
			if (!double_eq(value, 0.0))
			{
				par_record.push_back(make_pair(col, value));
			}
		}

		// write par num
		int par_num = par_record.size();
		_tempfile->write((char*)&par_num, SIZE_INT);
		len_size += SIZE_INT;

		// write parameter index
		_tempfile->write((char*)&idx, SIZE_INT);
		len_size += SIZE_INT;

		// write par name station satlltie
		
		const t_gpar& par_tmp = _x_solve[idx - 1];
		string name = par_tmp.site + "_" + par_tmp.str_type();
		if (par_tmp.prn == "") name = name + "_";
		_tempfile->write(name.c_str(), 20);
		len_size += 20;


		// write par time
		//  write start and end time 
		_tempfile->write((char*)&par_tmp.beg, sizeof(t_gtime));
		len_size += sizeof(t_gtime);
		_tempfile->write((char*)&par_tmp.end, sizeof(t_gtime));
		len_size += sizeof(t_gtime);


		// write parameter value
		double value = par_tmp.value();
		_tempfile->write((char*)&value, SIZE_DBL);
		len_size += SIZE_DBL;

		// write zhd for ztd retrieval glfeng
		double zhd = par_tmp.zhd;
		_tempfile->write((char*)&zhd, SIZE_DBL);
		len_size += SIZE_DBL;

		// write coefficient BTPB
		for (int ipar = 1; ipar <= par_num; ipar++)
		{
			_tempfile->write((char*)&(par_record[ipar - 1].first), SIZE_INT);
			len_size += SIZE_INT;
			_tempfile->write((char*)&(par_record[ipar - 1].second), SIZE_DBL);
			len_size += SIZE_DBL;
		}

		// write BTPL
		_tempfile->write((char*)&w_record, SIZE_DBL);
		len_size += SIZE_DBL;

		_tempfile->write((char*)&len_size, SIZE_INT);

		return 1;
	}



	int t_glsq::npar_number()const
	{
		return _npar_tot_num;
	}

    t_gtime str2gtime(string str_time)
	{

		replace(str_time.begin(), str_time.end(), '-', ' ');
		replace(str_time.begin(), str_time.end(), ':', ' ');
		stringstream time(str_time);
		int year, month, day, hour, minute;
		double sec;
		time >> year >> month >> day >> hour >> minute >> sec;
		t_gtime time_t; time_t.from_ymdhms(year, month, day, hour, minute, sec);
		return time_t;
	}
}
