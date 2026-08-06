/**
*
* @verbatim
	History
	 -1.0 GREAT     	 2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gcdr_podflt.h
* @brief		coder for lsq
*
* @author		GREAT, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#ifndef GCDR_PODFLT_H
#define GCDR_PODFLT_H

#include "gcoders/gallcoder_base.h"
#include "gcoders/gallcoder_orb.h"

using namespace gnut;
namespace great
{
	class t_gcdr_podflt :
		public t_gallcoder_base, 
		
		public t_gallcoder_orb
	
		
	{
	public:
		t_gcdr_podflt();
		virtual ~t_gcdr_podflt();

		bool support(const IFMT& ifmt) override;
		void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) override;
		void coder(t_gtime beg, t_gtime end, const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj);
		void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) override;
	protected:

	};
}

#endif