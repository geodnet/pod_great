/**
*
* @verbatim
	History
	 -1.0 GREAT     	 2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gcdr_podflt.cpp
* @brief		coder for lsq
*
* @author		GREAT, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#include "gcdr_podflt.h"
#include "gcfg_podflt.h"

namespace great
{
	t_gcdr_podflt::t_gcdr_podflt() :
		t_gallcoder_base(),
		t_gallcoder_orb()
		
	{
	}


	t_gcdr_podflt::~t_gcdr_podflt()
	{
	}

	bool t_gcdr_podflt::support(const IFMT& ifmt)
	{

		auto gset = dynamic_cast<t_gcfg_podflt*>(_gset);
		if (!_gset)
		{
			_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_gcdr_gns_podflt", "support", "no gset");
			return false;
		}
		// use brdm, skip the rinexc
		if (ifmt == IFMT::AMBFLAG12_INP)   return false;
		if (ifmt == IFMT::AMBFLAG13_INP) return false;
		if (gset->isUseBrdm && ifmt == IFMT::SP3_INP)     return false;
		if (!gset->isUseBrdm && ifmt == IFMT::RINEXN_INP) return false;
		if (gset->isUseBrdm && ifmt == IFMT::RINEXC_INP) return false;
		if (!gset->isUseIon && ifmt == IFMT::ION_INP)    return false;
		if (!gset->isAmbfix && ifmt == IFMT::AMBCON_INP) return false;

		if (t_gallcoder_base::support(ifmt)) return true;
		if (t_gallcoder_orb::support(ifmt)) return true;
	
		return false;
	}

	void t_gcdr_podflt::coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj)
	{
		if (t_gallcoder_base::support(ifmt))
		{
			t_gallcoder_base::coder(ifmt, path, data, obj);
			return;
		}

		if (t_gallcoder_orb::support(ifmt))
		{
			t_gallcoder_orb::coder(ifmt, path, data, obj);
			return;
		}


	}

	void t_gcdr_podflt::coder(t_gtime beg, t_gtime end, const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj)
	{
		if (t_gallcoder_base::support(ifmt))
		{
			t_gallcoder_base::coder(beg, end, ifmt, path, data, obj);
			return;
		}

		if (t_gallcoder_orb::support(ifmt))
		{
			t_gallcoder_orb::coder(ifmt, path, data, obj);
			return;
		}
	}

	void t_gcdr_podflt::creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata)
	{
		if (t_gallcoder_base::support(ifmt))
		{
			t_gallcoder_base::creat(ifmt, map_gdata);
			return;
		}
		if (t_gallcoder_orb::support(ifmt))
		{
			t_gallcoder_orb::creat(ifmt, map_gdata);
			return;
		}

	
	}

}
