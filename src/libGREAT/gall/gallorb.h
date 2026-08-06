/**
 *
 * @verbatim
	 History
	  -1.0  Hongjie Zheng   2019-04-10 Creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gorb.h
 * @brief		The class for storaging all orb data.
 *
 * @author      Hongjie Zheng, Wuhan University
 * @version		1.0.0
 * @date		2019-04-10
 *
 */


#ifndef GALLORB_H
#define GALLORB_H

#include "ExportHeaders/ExportBase.h"

#include "gdata/gdata.h"
#include "gcoders/orb.h"
#include "gdata/gorb.h"

using namespace gnut;

namespace great
{
	/// the class for storaging orb data
	/**
	*@brief	   Class for storaging orb data
	*
	* The class contains Doodson index and coeffcients of
	*		oceantide in different order and provide the
	*		interface to get data.
	*/
	class BASE_LIBRARY_EXPORT t_gallorb :public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gallorb();

		/** @brief default destructor. */
		virtual ~t_gallorb();

		typedef std::vector<t_gtriple>  t_sat_funct;



		void getAllPartial(map<string, map<int, t_sat_funct > >& funct);
		void getFirstHead(t_orb_hd& head);


		//
		/**
		* @brief interpolate a satellite's position(unit: m) in specific epoch.
		* @param[in]   epoch	epoch time
		* @param[in]   sat		satellite prn
		* @param[out]  px_sat	x coordinates
		* @param[out]  py_sat	y coordinates
		* @param[out]  pz_sat	z coordinates
		* @return   if interpolate successfully
		*	@retval   true   interpolate successfully
		*   @retval   false  interpolate unsuccessfully
		*/
		bool get_crs_sat_crd(const t_gtime& epoch, const string& sat, double& px_sat, double& py_sat, double& pz_sat);

		/**
		* @brief interpolate a satellite's velocity(unit: m/s) in specific epoch.
		* @param[in]   epoch	epoch time
		* @param[in]   sat		satellite prn
		* @param[out]  vx_sat	x velocity
		* @param[out]  vy_sat	y velocity
		* @param[out]  vz_sat	z velocity
		* @return   if interpolate successfully
		*	@retval   true   interpolate successfully
		*   @retval   false  interpolate unsuccessfully
		*/
		bool get_crs_sat_vel(const t_gtime& epoch, std::string sat, double& vx_sat, double& vy_sat, double& vz_sat);

		/*bool get_sat_vel_realtime(const t_gtime& epoch, std::string sat, double& vx_sat, double& vy_sat, double& vz_sat);*/
		/// only get partial dxdx0..dxdvx0...dxdnp by interpolating
		/**
		* @brief interpolate a satellite's partial matrix(only dxdx0..dxdvx0...dxdnp) in specific epoch.
		* @param[in]   epo		epoch time
		* @param[in]   sat		satellite prn
		* @param[out]  partial	partial matrix
		* @return   if interpolate successfully
		*	@retval   true   interpolate successfully
		*   @retval   false  interpolate unsuccessfully
		*/
		bool getGnsPartial(const t_gtime& epo, std::string sat, Matrix& partial);

		
		/**
		* @brief get the interpolate point
		* @param[in]   epo			epoch input(int)
		* @param[in]   norder		lagrange interpolate order
		* @param[in]   sat			satellite prn
		* @param[out]  epo_left		first interpolate point
		* @param[out]  epo_right	last interpolate point
		* @param[out]  epo_now		current epoch time(mjd)
		* @return   if get interpolate point successfully
		*	@retval   true   get interpolate poin successfully
		*   @retval   false  get interpolate poin unsuccessfully
		*/
		bool getInterpPoints(const t_gtime& epo, const int& norder, const string& sat, int& epo_left, int& epo_right, double& epo_now);



		void test_update(t_gtime beg_time,t_gtime end_time, int sampling,map<string, map<t_gtime, ColumnVector>>& newSatPos);


		// overried all results with new_orb
		bool update(std::string sat, map<int, vector<t_gtriple>> new_orb, t_gtime beg, t_gtime end, t_gtime t_ref, int oi_admout, bool keep_satpos);
	


	protected:
		map<string, pair<t_orb_hd,map<int, t_sat_funct > > > _orbs; /// allsat  epo funct
	};

}

#endif