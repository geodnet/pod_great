#ifndef GFANTENNATHRUST_H
#define GFANTENNATHRUST_H

#include "gutils/gconst.h"
#include "gforce/gforce.h"
#include "gforce/gfarg_init.h"
#include "gio/glog.h"
#include "gutils/gconst_oi.h"

/**
*
* @verbatim
History
-1.0         Liu Chengbo              2020-09-18 Creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		      gfantennathrust.h
* @brief		    The class for antenna thrust.
*
* @author          Liu Chengbo, Wuhan University
* @version		    1.0.0
* @date		    2020-09-18
*
*/

#include "ExportHeaders/ExportOI.h"
#include "newmat/newmat.h"

using namespace std;

namespace great
{
	/** @brief class for earth radiation model. */
	class OI_LIBRARY_EXPORT t_gfantennathrust : public t_gforce
	{
	public:

		/** @brief constructor */
		t_gfantennathrust(t_init_arg* arg, t_glog* log);

		/** @brief destructor */
		virtual ~t_gfantennathrust();

		
		bool calculate(t_gforce_calarg* argument);

		bool getTransmitPower(string block);

	protected:

		double     _mass;       ///< the mass of sat.
		string blockname;       ///< sat block name
		string model;           ///< antenna model.
		double antennapower = 0.0;   ///< satellite transmit power

	};
}

#endif