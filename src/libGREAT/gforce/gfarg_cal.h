/**
 * This file contains all calculating arguments of diffetent force models used in OI.
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang	 2019-01-14 creat the file.
	  -1.1 Jiande Huang  2019-04-01 Adding abnormal throwing mechanism of functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfarg_cal.h
 * @brief		The argument used in the caculation of force model.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */

#ifndef GFARG_CAL_H
#define GFARG_CAL_H

#include <map>
#include "ExportHeaders/ExportOI.h"
#include "gutils/gconst_oi.h"
#include "gutils/gtime.h"
#include "gforce/gfcommon.h"
#include "gset/gsetbase.h"

#include "newmat/newmat.h"

namespace great
{
	/** @brief argument used for init force model*/
	enum ARG_TYPE
	{	
		ARG_POINT_MASS,	   ///< N body attraction.
		ARG_SOLAR_ECOM2,   ///< ECOM2 model for solar radiation.
		ARG_GRAVI_EGM,	   ///< EGM 2008 gravity model.
		ARG_RELATIVITY,	   ///< Default relativity force model.
		ARG_ATMOSPHERE,     ///< Default atmosphere force model.
		ARG_ANTENNATHRUST,  ///< Antenna thrust model.
		ARG_UNDEF
	};

	/** @brief base calss for caculate argument */
	class t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gforce_calarg() { arg_type = ARG_TYPE::ARG_UNDEF; };
		
		/** @brief default destructor */
		virtual ~t_gforce_calarg() {};

		double*			time = nullptr;		///< current epoch, should be TT or TAT.
		ARG_TYPE		arg_type;	        ///< force model type.
		ColumnVector*	sat_pos = nullptr;	///< sat position in J2000.
		ColumnVector*   sat_vel = nullptr;	///< sat velocity in J2000.
	};

	/** @brief point mass argument */
	class t_gfcalarg_point_mass : public t_gforce_calarg
	{
	public:
	    /** @brief default constructor */
		t_gfcalarg_point_mass() { arg_type = ARG_POINT_MASS; };
		
		/** @brief default destructor */
		virtual ~t_gfcalarg_point_mass() {};
	};

	class t_gfcalarg_antennathrust : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_antennathrust() { arg_type = ARG_ANTENNATHRUST; }

		/** @brief default destructor */
		virtual ~t_gfcalarg_antennathrust() {};
		double sat_mass = 0.0;
		string sat_blocknum;
		double sat_power = 0.0;
	};



	/** @brief ECOM2 model argument */
	class t_gfcalarg_solar_ecom2 : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_solar_ecom2() { arg_type = ARG_SOLAR_ECOM2; }

		/** @brief default destructor */
		virtual ~t_gfcalarg_solar_ecom2() {};

		double Beta() {
			// unit vector of orbit plane
			ColumnVector unit_Nop = crossproduct(*(this->sat_pos), *(this->sat_vel));
			unit_Nop = unit_Nop / unit_Nop.norm_Frobenius();

			// unit vector of _sun_pos
			ColumnVector unit_sun = *(this->sun_pos) / this->sun_pos->norm_Frobenius();

			return asin(dotproduct(unit_Nop, unit_sun) / sqrt(dotproduct(unit_sun, unit_sun)));
		};


		ColumnVector* sun_pos = nullptr;	///< sun pos in J2000.
		ColumnVector* moon_pos= nullptr;	///< moon pos in J2000.
		double*  shadow_factor= nullptr;	///< shadow factor, which is relative with the positions of satellite, moon and sun.
		string  sat_type;		///< type of satellite, specially for BLOCK IIR.

	};





	/** @brief gravity argument */
	class t_gfcalarg_gravity_egm : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_gravity_egm() { arg_type = ARG_GRAVI_EGM; };
		
		/** @brief default destructor */
		virtual ~t_gfcalarg_gravity_egm() {};

		
		ColumnVector *sat_pos_trs = nullptr;			///< sat pos in TRS.
		map<string, ColumnVector> *allXE = nullptr;	///< planet positions in TRS.
		double* xpole = nullptr;						///< X coordinates for pole point. 
		double* ypole = nullptr;						///< Y coordinates for pole point. 
		double* gst   = nullptr;						///< Greenwich Time.
		Matrix* e2j = nullptr;						///< Rot matrix from TRS to J2000.
		modeofmeanpole meanpole = modeofmeanpole::cubic;
	};

	
	/** @brief relativity argument */
	class t_gfcalarg_relativity : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_relativity() { arg_type = ARG_RELATIVITY; };
		
		/** @brief default destructor */
		virtual ~t_gfcalarg_relativity() {};

		
		double  gm = EARTH_GM;				///< GM for earth.
	};
}
#endif // !GFARGUMENT_H
