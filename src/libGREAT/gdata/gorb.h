/**
 *
 * @verbatim
	 History
	  -1.0  Jie Li       2019-01-06 Creat the file.
	  -1.1  Jie Li       2019-03-30 Code comments are added in Doxygen format.
	  -1.2  Huang Jiande 2020-09-12 Add function get sat list to fix the format of the orb file 
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gorb.h
 * @brief		The class for storaging orb data.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-03-30
 *
 */


#ifndef GORB_H
#define GORB_H

#include "ExportHeaders/ExportBase.h"
#include "gdata/gdata.h"
#include "gcoders/orb.h"
#include "gall/gallpar.h"

using namespace gnut;
using namespace great;
namespace great
{
	/**
	*@brief	   Class for storaging one orb file data
	*
	* The class contains orb header information and posation and partial matrix of satellites,
	*		and provide the interface to get orb data.
	*/
	class BASE_LIBRARY_EXPORT t_gorb :public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gorb();

		/** @brief default denstructor. */
		virtual ~t_gorb();

		typedef vector<t_gtriple>  t_sat_funct;

	

		/**
		* @brief get orb head data.
		* @param[out]   orb_hd	orb head data
		*/
		void get_head(t_orb_hd& orb_hd);
		
		/**
		* @brief get whole orb funct data.
		* @param[out]   funct	whole orb funct data
		*/
		void get_funct(map<string, map<int, t_sat_funct > >& funct);

		/**
		* @brief get orb funct data of a satellite.
		* @param[in]   sat		satellite prn
		* @param[out]  funct	orb funct data of the sat
		*/
		void get_funct(string sat, map<int, t_sat_funct >& funct);

		

		/**
		* @brief get orb funct data of a satellite.
		* @param[in]   sat		satellite prn
		* @param[out]  parlist	parlist data of the sat
		*/
		void get_parlist(map<string,vector<string> >& parlist);

	protected:
		t_orb_hd	_orb_hd;	///< orb head data	
		map<string, vector<string> > _par_list; ///< partial order of sat(sat,parlist)
		map<string, map<int, t_sat_funct > > _funct;  ///< orb funct (sat , epo  ,data)
	};

}

#endif