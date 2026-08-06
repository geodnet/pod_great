/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBlSat.h
* @brief	sattelite baseline of SD and DD ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#ifndef GAMBBLSAT_H
#define GAMBBLSAT_H 

#include "ExportHeaders/ExportAMB.h"
#include "gio/glog.h"
#include "gutils/gnss.h"
#include "gambfix/gambBl.h"

using namespace std;
using namespace gnut;

namespace great
{
	
	class AMB_LIBRARY_EXPORT t_gamb_baseline_sats
	{
	public:
		
		t_gamb_baseline_sats();
		virtual ~t_gamb_baseline_sats();

		
		void add_sat(const string& sat);

		vector<pair<string, string>>  baselines();
	protected:
		vector<t_gamb_baseline>  _baselines;
		vector<string>           _sats;
	};
}

#endif