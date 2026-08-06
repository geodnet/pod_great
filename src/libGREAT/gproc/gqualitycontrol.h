
/**
* @verbatim
	History
	 -1.0 Guolong Feng	 2020-07-05  creat the file.
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file			qualityControl.h
* @author       Guolong Feng, Wuhan University
* @version		1.0.0
* @date			2020-07-05
*
*/

#ifndef QUALITYCONTROL_H
#define QUALITYCONTROL_H

#include "gall/gallnav.h"
#include "gall/gallobs.h"

#include "gdata/gsatdata.h"

#include "gset/gsetgen.h"
#include "gset/gsetrec.h"
#include "gset/gsetgnss.h"
#include "gset/gsetproc.h"
#include "gset/gsetturboedit.h"

#include "ExportHeaders/ExportBase.h"

namespace great
{
	enum SMOOTH_MODEL
	{
		SMT_DOPPLER,
		SMT_PHASE,
		SMT_NONE
	};



	
	class BASE_LIBRARY_EXPORT t_goutliers_process
	{
	public:
		t_goutliers_process(t_gsetbase* settings, t_glog* glog);

		virtual ~t_goutliers_process();

		void setLog(t_glog* glog) { _log = glog; };

		void excludeBadObs(vector<t_gsatdata>& obsdata); // Low SNR < 10dB

		void flagRangeOutliers(shared_ptr<t_gobsgnss> ObsPre, shared_ptr<t_gobsgnss> Obs, double sampling);

	private:
		t_glog*            _log;
		t_gsetbase*        _set;
		t_giof*            _debug_outliers=nullptr;

	};

	

} // namespace

#endif
