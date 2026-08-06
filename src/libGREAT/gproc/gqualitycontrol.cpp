/**
* @verbatim
	History
	 -1.0 Guolong Feng	 2020-07-05  creat the file.
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file			qualityControl.cpp
* @author       Guolong Feng, Wuhan University
* @version		1.0.0
* @date			2020-07-05
*
*/

#include <iomanip>
#include <memory>

#include "gproc/gqualitycontrol.h"

namespace great {


	t_goutliers_process::t_goutliers_process(t_gsetbase* settings, t_glog* glog):
		_log(glog)
	{
		if (1) {
			if (_debug_outliers)
			{
				if (_debug_outliers->is_open()) _debug_outliers->close();
				delete _debug_outliers;
				_debug_outliers = nullptr;
			}
		}
		else {
			_debug_outliers = new t_giof;
			_debug_outliers->tsys(t_gtime::GPS);
			_debug_outliers->mask("debug_outliers.log");
			_debug_outliers->append(false);
		}
	}

	t_goutliers_process::~t_goutliers_process()
	{
		if (_debug_outliers)
		{
			if (_debug_outliers->is_open()) _debug_outliers->close();
			delete _debug_outliers;
			_debug_outliers = nullptr;
		}
	}

	void t_goutliers_process::flagRangeOutliers(shared_ptr<t_gobsgnss> ObsPre, shared_ptr<t_gobsgnss> Obs, double sampling)
	{
		ostringstream os;
		os.str("");
		t_gtime crt_time = Obs->epoch();
		string  crt_site = Obs->site();

		vector<GOBS> obs_vec = Obs->obs();

		// Gao Y et al. Modeling and estimation of C1-P1 bias in GPS receivers[J]. Journal of Geodesy, 2001.
		GSYS gsys = Obs->gsys();
		GOBSBAND b1 = GNSS_BAND_PRIORITY.at(gsys)[1];
		GOBS ref_obsP = Obs->select_range(b1);
		double ref_P = Obs->obs_C(ref_obsP);  // m

		for (auto obs_type : obs_vec)
		{
			// skip not code obs
			if (!t_gobs(obs_type).is_code()) continue;
			// skip ref obs(add case: ref_obsP == X)
			if (obs_type == ref_obsP || ref_obsP == X) continue;

			t_gobs gobs(obs_type);
			GOBSBAND bX = gobs.band();

			double obs_P = Obs->obs_C(obs_type);
			double diff_P = abs(obs_P - ref_P);

			double value = (b1 == bX) ? 10 : 30;  // meters

			os << fixed << setw(10) << " PP[m] " << Obs->sat() << "  " << setw(6) << ObsPre->epoch().sod() << setw(6) << crt_time.sod() <<
				"  " << gobs2str(obs_type) << "  " << gobs2str(ref_obsP) << setw(15) << setprecision(4) << diff_P << endl;

			if (diff_P > value)
			{
				Obs->addoutliers(obs_type, 1);
				Obs->addoutliers(ref_obsP, 1);

				if (_log)
				{
					ostringstream msg;
					msg << fixed << setw(10) << " Outliers in Range : " << crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(ref_obsP)
						<< "   " << gobs2str(obs_type) << "  value  " << setw(16) << setprecision(3) << diff_P << "  >  " << setw(16) << setprecision(3) << value;
					_log->comment(1, crt_site, msg.str());
				}
			}
		}

		// Assume GNSS observations are OK

		if (crt_time.diff(ObsPre->epoch()) > sampling) return;

		for (auto obs_doppler : obs_vec)
		{
			// skip not doppler obs
			if (!t_gobs(obs_doppler).is_doppler()) continue;

			double crt_doppler = Obs->getobs(obs_doppler);      // cycle / s
			double pre_doppler = ObsPre->getobs(obs_doppler);   // cycle / s

			// check Obstype Gap
			if (double_eq(pre_doppler, NULL_GOBS) || double_eq(crt_doppler, NULL_GOBS)) continue;

		
			double delta_doppler = crt_doppler - pre_doppler;

			
			if (obs_doppler == D1N)
			{
				if (abs(delta_doppler) > 10)  Obs->addoutliers(obs_doppler, 1);
				continue;
			}

			string strGOBS = gobs2str(obs_doppler);
			strGOBS.replace(0, 1, "C");
			GOBS obs_range = str2gobs(strGOBS);

			if (Obs->getoutliers(obs_range) > 0)  continue;

			GOBSBAND b = t_gobs(obs_doppler).band();
			crt_doppler *= Obs->wavelength(b);     // m/s
			pre_doppler *= ObsPre->wavelength(b);  // m/s

			double crt_range = Obs->obs_C(obs_range);    // m
			double pre_range = ObsPre->obs_C(obs_range); // m

			if (double_eq(crt_range, NULL_GOBS) || double_eq(pre_range, NULL_GOBS)) continue;

			double delta_PD = crt_range - pre_range + 0.5 * (crt_doppler + pre_doppler) * sampling;  // m

			os << fixed << setw(10) << " D[c/s] " << Obs->sat() << "  " << setw(6) << ObsPre->epoch().sod() << setw(6) << crt_time.sod() <<
			    "  " << gobs2str(obs_doppler) << setw(15) << setprecision(4) << delta_doppler << endl;

			os << fixed << setw(10) << " PD[m] " << Obs->sat() << "  " << setw(6) << ObsPre->epoch().sod() << setw(6) << crt_time.sod() <<
				"  " << gobs2str(obs_range) << "  " << gobs2str(obs_doppler) << setw(15) << setprecision(4) << delta_PD << endl;

			
			if (abs(delta_doppler) > 10 && abs(delta_PD) > 2)     // cycle/s, m
			{
				Obs->addoutliers(obs_doppler, 1);       // doppler
				Obs->addoutliers(obs_range, 1);         // range
				if (_log)
				{
					ostringstream msg;
					msg << fixed << setw(10) << " Outliers in Doppler and Range : " << crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(obs_doppler)
						<< "   " << gobs2str(obs_range) << "  delta_doppler[cycle/s] =  " << setw(16) << setprecision(3) << delta_doppler
						<< "  delta_PD[m] =  " << setw(16) << setprecision(3) << delta_PD;
					_log->comment(1, crt_site, msg.str());
				}
			}
					else if (abs(delta_PD) > 5)  // m
			{
				Obs->addoutliers(obs_range, 1);         // range
				if (_log)
				{
					ostringstream msg;
					msg << fixed << setw(10) << " Outliers in Range : " << crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(obs_range)
						<< "  delta_PD[m] =  " << setw(16) << setprecision(3) << delta_PD;
					_log->comment(1, crt_site, msg.str());
				}
			}
		}

		// Print flt results
		if (_debug_outliers) {
			_debug_outliers->write(os.str().c_str(), os.str().size());
			_debug_outliers->flush();
		}
	}

	void t_goutliers_process::excludeBadObs(vector<t_gsatdata>& obsdata)
	{
		vector<t_gsatdata>::iterator it;
		bool valid_SNR, valid_Code, valid_Phase;
		bool lack_Phase, single_freq;
	
		for (it = obsdata.begin(); it != obsdata.end();)
		{
			 valid_SNR = true;
			 valid_Code = true;
			 valid_Phase = true;
			 lack_Phase = true;
			 single_freq = false;
			GSYS gsys = it->gsys();

			vector<GOBS> obs_vec = it->obs();
			for (auto obs_type : obs_vec)
			{   
				GOBSTYPE type = t_gobs(obs_type).type();
				
				// check_SNR
				if ( type == TYPE_S)
				{
					double obs_snr = it->getobs(obs_type);
					if (obs_snr <= 10.0) valid_SNR = false;
					else valid_SNR = true;
				}
				// check Range
				else if (type == TYPE_C || type == TYPE_P)
				{
					double obs_code = it->getobs(obs_type);
					if (obs_code <= 1.9e7 ) valid_Code = false;
					else valid_Code = true;
				}
				// check Phase
				else if (type == TYPE_L)
				{
					double obs_phase = it->obs_L(obs_type);
					if (obs_phase <= 1.9e7) valid_Phase = false;
					else if (lack_Phase)    lack_Phase = false;
					else  valid_Phase = true; 
				}		
			}
			//LX added for delete the single-frequency obs in dual mode
			/*set<GOBSBAND> bands = it->band_avail(true);
			if ( bands.size() < 2)  single_freq = true;*/

			if (!valid_SNR || !valid_Phase || !valid_Code || lack_Phase)
			{
				it = obsdata.erase(it); 
				continue;
			}

			it++;
		}
	}


} // namespace
