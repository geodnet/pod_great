/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li	     2019-01-06	creat the file.
	  -1.1	Jie Li	     2019-03-30	add doxygen style code remarks
	  -1.2  Jiande Huang 2020-09-12 Standard code style and Optimize code logic
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gegmcoef.h
 * @brief		The class used to save egm file information.
 *              We consider solar radiation, gravity and so on.\n
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 *
 */

#ifndef GEGMCOEF_H
#define GEGMCOEF_H

#include "gdata/gdata.h"
#include "gutils/gtime.h"
#include <vector>

#include "ExportHeaders/ExportOI.h"

using namespace std;
using namespace gnut;
namespace great
{
	/// The base struct for one record of gravity coeffient
	struct t_gravity_coef { double c = 0.0, s = 0.0, sigma_c = 0.0, sigma_s = 0.0; };

	struct t_time_gravity_coef
	{
		double c[81][81] = { {0.0} }, s[81][81] = { {0.0} };
		double sigma_c[81][81] = { {0.0} }, sigma_s[81][81] = { {0.0} };
		double ctrnd[81][81] = { {0.0} }, strnd[81][81] = { {0.0} };
		double csin1[81][81] = { {0.0} }, ccos1[81][81] = { {0.0} };
		double csin2[81][81] = { {0.0} }, ccos2[81][81] = { {0.0} };
		double ssin1[81][81] = { {0.0} }, scos1[81][81] = { {0.0} };
		double ssin2[81][81] = { {0.0} }, scos2[81][81] = { {0.0} };
		t_gtime st, et; /// start and end time of the record
	};

	/**
	*@brief	   Class for save the EGM file data
	*
	* This class save the EGM file information and provide the interface to get data
	*/
	class  OI_LIBRARY_EXPORT t_gegmcoef :public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gegmcoef();

		/** @brief default destructor. */
		virtual ~t_gegmcoef();

		/**
		* @brief add reference time, GM and radius of the EARTH.
		* @param [in]  t       reference time of EGM data
		* @param [in]  GM      GM of the Earth
		* @param [in]  radius  radius of the Earth
		* @return
			@retval =0 add data successfully
		*/
		int add_head(t_gtime t, double GM, double radius);
		
		/**
		* @brief add max EGM order.
		* @param [in]  order    EGM order
		* @return
			@retval =0 add data successfully
		*/
		int add_order(int order);

		/**
		* @brief add one coefficient recoder into map data
		* @param [in]  n    index
		* @param [in]  m    index
		* @param [in]  c    c value
		* @param [in]  s    s value
		* @param [in]  sigma_c    sigma_c value
		* @param [in]  sigma_s    sigma_s value
		* @return	index n
		*/
		int add_egmcoef(int n, int m, double c, double s, double sigma_c, double sigma_s);

		/**
        * @brief add a new time coefficient group into map data
		* @param [in]  st    start time
		* @param [in]  et    end time
		* @return	index ipt
		*/
		int create_timecoef_group(t_gtime st, t_gtime et);

		/**
		* @brief add one time coefficient into the corresponding group
		* @param [in]  ipt  group index of time coefficient 
		* @param [in]  tag  type of time coefficient 
		* @param [in]  n    degree index
		* @param [in]  m    order index
		* @param [in]  c    c value
		* @param [in]  s    s value
		* @param [in]  sigma_c    sigma_c value
		* @param [in]  sigma_s    sigma_s value
		* @return	index n
		*/
		int add_timecoef(int ipt, string tag, int n, int m, double c, double s, double sigma_c, double sigma_s, double p);

		/**
		* @brief find the index of group which the time coefficient record belongs to (judge by the start time)
		* @param [in]  t     start time of record 
		* @return	   ipt   index of group
		*/
		int find_group(t_gtime t);

		/**
		* @brief get the time variable spherical harmonic coefficient group according to selected mjd
		* @param [in]  dmjd  mjd   
		* @return	   t_time_gravity_coef
		*/
		const t_time_gravity_coef& get_time_coefficient(double dmjd);

		/** @brief whether the EGM coefficien data is empty.
		* @return  bool
		*	@retval   true    EGM coefficien data is empty
		*   @retval   false   EGM coefficien data is existent
		*/
		bool is_empty();

		/**
		* @brief get corresponding c value.
		* @param [in]  n    index
		* @param [in]  m    index
		* @return  c value
		*/
		double get_c(int n, int m);
		/**
		* @brief get corresponding s value.
		* @param [in]  n    index
		* @param [in]  m    index
		* @return  s value
		*/
		double get_s(int n, int m);

		/**
		* @brief get reference time(mjd) of EGM file.
		* @return  mjd of reference time.
		*/
		double mjd_ref();

		/**
		* @brief get GM of the Earth.
		* @return  GM of the Earth.
		*/
		double GM() { return _GM; }

		/**
		* @brief get radius of the Earth.
		* @return  radius of the Earth.
		*/
		double radius() { return _radius; }

	protected:
		map<int, map<int, t_gravity_coef> > _map_egmcoef;		///< map storage EGM coefficient 
		map<int, t_time_gravity_coef> _map_time_egmcoef;	///< map storage time EGM coefficient
		string  _name;			///< file name
		int     _order;			///< EGM coefficient order
		double  _GM;			///< GM of Earth(unit: km/s)
		double  _radius;		///< radius of Earth(unit: km)
		t_gtime _t;				///< reference time(mjd)
	};
}
#endif // !GEGMCOEF_H
