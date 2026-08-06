/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-05-15 created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    t_gcombmodel.cpp
* @brief	    base combine biase model
* @mainpage     C++ version
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-05-18
*
*/

#include "gmodels/gcombmodel.h"
#include "gutils/gcommon.h"
#include <math.h>


namespace great
{
	t_gcombmodel::t_gcombmodel(t_gsetbase* setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data) :
		_bias_model(bias_model),
		_log(log)
	{
		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::QZS);
		
		

		_frequency = dynamic_cast<t_gsetproc*>(setting)->frequency();
		_ion_model = dynamic_cast<t_gsetproc*>(setting)->ion_model();
		_ifb_model = dynamic_cast<t_gsetproc*>(setting)->ifb_model();
		_ifcb_model = dynamic_cast<t_gsetproc*>(setting)->ifcb_model();
		
		_gallbias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gprecisemodel = dynamic_cast<t_gprecisemodel*>(_bias_model->precisemodel());

		

	}

	t_gcombmodel::~t_gcombmodel()
	{
	}



	t_gcombIF::t_gcombIF(t_gsetbase* setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data) :
		t_gcombmodel(setting, log, bias_model, data)
	{
		_clk_type_index[make_pair(FREQ_3, GPS)] = t_gpar::t_type::CLK13_G;
		_clk_type_index[make_pair(FREQ_4, GPS)] = t_gpar::t_type::CLK14_G;
		_clk_type_index[make_pair(FREQ_5, GPS)] = t_gpar::t_type::CLK15_G;

		_clk_type_index[make_pair(FREQ_3, GAL)] = t_gpar::t_type::CLK13_E;
		_clk_type_index[make_pair(FREQ_4, GAL)] = t_gpar::t_type::CLK14_E;
		_clk_type_index[make_pair(FREQ_5, GAL)] = t_gpar::t_type::CLK15_E;

		_clk_type_index[make_pair(FREQ_3, BDS)] = t_gpar::t_type::CLK13_C;
		_clk_type_index[make_pair(FREQ_4, BDS)] = t_gpar::t_type::CLK14_C;
		_clk_type_index[make_pair(FREQ_5, BDS)] = t_gpar::t_type::CLK15_C;

		_clk_type_index[make_pair(FREQ_3, QZS)] = t_gpar::t_type::CLK13_J;
		_clk_type_index[make_pair(FREQ_4, QZS)] = t_gpar::t_type::CLK14_J;
		_clk_type_index[make_pair(FREQ_5, QZS)] = t_gpar::t_type::CLK15_J;
	}

	t_gcombIF::~t_gcombIF()
	{
	}

	bool t_gcombIF::cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result)
	{
		// check Obs type
		GOBSBAND b1 = _band_index[obsdata.gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[obsdata.gsys()][FREQ_2];

		auto band1 = GNSS_BAND_PRIORITY.at(obsdata.gsys())[1];
		auto band3 = GNSS_BAND_PRIORITY.at(obsdata.gsys())[3];

		obsdata.apply_bias(*_gallbias);
		
		return this->cmb_equ_IF(epoch, params, obsdata, b1, b2, result);
	}

	bool t_gcombIF::cmb_equ_IF(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, GOBSBAND b1, GOBSBAND b2, t_gbaseEquation& result)
	{
		// ========================================================================================================
		t_gobs gobsP1 = obsdata.select_range(b1);
		t_gobs gobsP2 = obsdata.select_range(b2);
		t_gobs gobsL1 = obsdata.select_phase(b1);
		t_gobs gobsL2 = obsdata.select_phase(b2);

		vector<pair<t_gobs, t_gobs> > type_list;
		type_list.push_back(make_pair(gobsP1, gobsP2));
		type_list.push_back(make_pair(gobsL1, gobsL2));

		auto gsys = obsdata.gsys();
		// ========================================================================================================
		// IF coef
		// 
		double coef1, coef2;
		obsdata.coef_ionofree(b1, coef1, b2, coef2);

		// ========================================================================================================
		//
		t_glsqEquationMatrix equ_IF;
		for (auto item : type_list)
		{
			t_gobs gobs1 = item.first;
			t_gobs gobs2 = item.second;

			//combine f1 and f2
			t_gbaseEquation temp_equ;
			if (!_bias_model->cmb_equ(epoch, params, obsdata, gobs1, temp_equ)) return false;
			if (!_bias_model->cmb_equ(epoch, params, obsdata, gobs2, temp_equ)) return false;

			// ========================================================================================================
			//
			if (temp_equ.B[0].size() != temp_equ.B[1].size()) throw logic_error("coeff size is not equal in f1 and f2");

			vector<pair<int, double> > coef_IF;
			double P_IF = 0.0, l_IF = 0.0;
			for (int i = 0; i < temp_equ.B[0].size(); i++)
			{
				if (temp_equ.B[0][i].first != temp_equ.B[1][i].first) 
				{
					throw logic_error("coeff par is not the same in f1 and f2");
				}
				// combine coeff
				coef_IF.push_back(make_pair(temp_equ.B[0][i].first, coef1 * temp_equ.B[0][i].second + coef2 * temp_equ.B[1][i].second));
			}
			// combine P and l
			P_IF = 1.0 / (pow(coef1, 2) / temp_equ.P[0] + pow(coef2, 2) / temp_equ.P[1]);
			l_IF = coef1 * temp_equ.l[0] + coef2 * temp_equ.l[1];


			// ========================================================================================================
			// CORRECT AMB
			if (gobs1.is_phase() && (!obsdata.is_carrier_range(b1) || !obsdata.is_carrier_range(b2)))
			{
				t_gpar::t_type amb_type = t_gpar::NO_DEF;
				if (_freq_index[gsys][b2] == FREQ_2)  amb_type = t_gpar::AMB_IF;
				if (_freq_index[gsys][b2] == FREQ_3)  amb_type = t_gpar::AMB13_IF;
				if (_freq_index[gsys][b2] == FREQ_4)  amb_type = t_gpar::AMB14_IF;
				if (_freq_index[gsys][b2] == FREQ_5)  amb_type = t_gpar::AMB15_IF;

				int idx = params.getParam(obsdata.site(), amb_type, obsdata.sat()); 
				if (idx < 0) return false;

				// if amb is new  then init value with Obs - Modelobs				
				if (double_eq(params[idx].value(), 0.0) || params[idx].beg == epoch)
				{
					double obs_P3 = obsdata.P3(gobsP1, gobsP2);
					double obs_L3 = obsdata.L3(gobsL1, gobsL2);					
				
					params[idx].value(obs_L3 - obs_P3);
				}
				// update B l
				coef_IF.push_back(make_pair(idx + 1, 1.0));
				l_IF -= params[idx].value();
			}
			t_gobscombtype type(gobs1, b1, b2, _freq_index[gsys][b1], _freq_index[gsys][b2], OBSCOMBIN::IONO_FREE);
            equ_IF.add_equ(coef_IF,P_IF,l_IF,obsdata.site(),obsdata.sat(), type, false);
        }

        if (dynamic_cast<t_glsqEquationMatrix*>(&result))
        {
            t_glsqEquationMatrix *lsq_result = dynamic_cast<t_glsqEquationMatrix *>(&result);
            lsq_result->add_equ(equ_IF);
        }
        else
        {
            result = result + equ_IF;
        }
	
        return true;
    }
}