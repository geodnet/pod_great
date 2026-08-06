/**
* The EIGEN6S4 file is used for calculating gravity.\n
* Following is an example for the header of EGM file.
*
* @verbatim
	 FIRST GSM-2_0110_2002214-2003242_EIGEN_G---_0001 SHM     GFZ POTSDAM  20040213
	 CMMNT GRACE-only mean gravity field
	 CMMNT Reported standard deviations are formal (not calibrated)
	 EARTH 0.3986004415E+15 0.6378136460E+07  19970101							! GM(*1E-9 km/s) radius(unit: m) and refrence time of EARTH
	 SHM     150  150 1.00 fully normalized exclusive permanent tide			! max order of the coefficient
   @endverbatim
*
* An example for the record of the EIGEN6S4 file
* @verbatim
  gfc    0    0  1.00000000000E+00  0.00000000000E+00 0.0000E+00 0.0000E+00
  gfct   1    0  7.01776593302E-11  0.00000000000E+00 1.3180E-13 0.0000E+00 19500101.0000 20020815.0817
  trnd   1    0  0.00000000000E+00  0.00000000000E+00 0.0000E+00 0.0000E+00 19500101.0000 20020815.0817
  acos   1    0  7.81500675516E-11  0.00000000000E+00 3.5640E-14 0.0000E+00 19500101.0000 20020815.0817 1.0
  asin   1    0 -2.91423806169E-12  0.00000000000E+00 3.6200E-14 0.0000E+00 19500101.0000 20020815.0817 1.0
  acos   1    0 -3.34601799827E-11  0.00000000000E+00 3.4480E-14 0.0000E+00 19500101.0000 20020815.0817 0.5
  asin   1    0 -2.09359348050E-11  0.00000000000E+00 3.3920E-14 0.0000E+00 19500101.0000 20020815.0817 0.5
  @endverbatim
*
* @verbatim
	History
	 -1.0 Keke Zhang	   2020-12-29 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			eigen6s4.h
* @brief		decode the EIGEN6S4 data
*
* @author		Keke Zhang, Wuhan University
* @version		1.0.0
* @date			2020-12-29
*
*/

#ifndef EIGEN6S4_H
#define EIGEN6S4_H

#include "ExportHeaders/ExportOI.h"
#include "gcoders/gcoder.h"

using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for decoding the EIGEN6S4 data
	*
	* The document contains the coefficiet and sigma values in different order of the gravitational field
	* as well as the corresponding time coefficiet up to d/o 80
	* The gcoder EIGEN6S4 corresponding to the gdata egmcoef like EGM2008
	*/
	class OI_LIBRARY_EXPORT t_eigen6s4 : public t_gcoder
	{
	public:
		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_eigen6s4(t_gsetbase *s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_eigen6s4();

		/**
		* @brief decode header of the EIGEN6S4 data file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		*/
		int decode_head(char *buff, int sz, vector<string> &errmsg);

		/**
		* @brief decode the data body of the EIGEN6S4 data file.
		*
		* decode data body of EIGEN6S4 file, all the data will be stored in the gegmcoef
		*
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
			@retval >=0 consume size of data decoding
			@retval <0  finish reading
		*/
		int decode_data(char *buff, int sz, int &cnt, vector<string> &errmsg);

	protected:
	};
};
#endif
