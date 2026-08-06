/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBlSat.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambBlSat.h"
#include "gset/gsetamb.h"
#include "gset/gsetrec.h"

#include "newmat/newmatap.h"
namespace great
{
	great::t_gamb_baseline_sats::t_gamb_baseline_sats()
	{
	}



	great::t_gamb_baseline_sats::~t_gamb_baseline_sats()
	{

	}

	void t_gamb_baseline_sats::add_sat(const string& sat)
	{
		if (std::find(_sats.begin(), _sats.end(), sat) == std::end(_sats))
		{
			_sats.push_back(sat);
			std::sort(_sats.begin(), _sats.end());
		}
	}

	
	vector<pair<string,string>> great::t_gamb_baseline_sats::baselines()
	{
		vector<pair<string, string>> tmp;
		if (_baselines.empty())
		{
			if (_sats.empty()) return tmp;
			for (int i = 0; i < _sats.size(); i++)
			{
				for (int j = i + 1; j < _sats.size(); j++)
				{
					if (_sats[i].substr(0, 1) != _sats[j].substr(0, 1)) continue;
					t_gamb_baseline line(_sats[i], _sats[j], 0, 1);
					_baselines.push_back(line);
				}
			}
		}

		for (const auto& iter : _baselines)
		{
			tmp.push_back(iter.bl_pair());
		}

		return tmp;
	}


}











