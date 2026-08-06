/**
*
* @verbatim
	History
	 -1.0 GREAT  	 2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gcfg_podflt.h
* @brief		LSQ control set from XML
*
* @author		GREAT, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/


#ifndef GCFG_PODFLT_H
#define GCFG_PODFLT_H

#include <string>
#include <iostream>
#include <signal.h>

#include "ExportHeaders/ExportLSQ.h"
#include "gall/gallprec.h"
#include "gall/gallpcv.h"
#include "gall/gallobs.h"
#include "gall/gallotl.h"
#include "gall/gallbias.h"
#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexn.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gmodels/gpcv.h"

#include "gproc/gpodfltIF.h"

#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetrec.h"
#include "gset/gsetpar.h"
#include "gset/gsetamb.h"
#include "gset/gsetturboedit.h"
#include "gset/gsetforce.h"
using namespace std;
using namespace pugi;

namespace great {

	class t_gcfg_podflt : public t_gsetgen,
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgnss,
		public t_gsetproc,
		public t_gsetrec,
		
		public t_gsetpar,
		public t_gsetamb,
		public t_gsetturboedit,
		public t_gsetforce
	
	{

	public:
		t_gcfg_podflt();
		~t_gcfg_podflt();

		void check();                                 // settings check
		void help();                                  // settings help

		void arg(int argc, char* argv[], bool add, bool thin) override;
		bool check_inp();

		bool isUseBrdm = false;
		bool isAmbfix = false;
		bool isUseIon = false;

	protected:

	private:

	};

} // namespace

#endif