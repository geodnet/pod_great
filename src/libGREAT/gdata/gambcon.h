/**
*
* @verbatim
History
	-1.0 Jie Li  2019-11-14  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambcon.h
* @brief	Storage the ambigulity constraint files' data
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-11-14
*
*/

#ifndef GAMBCON_H
#define GAMBCON_H

#include "gdata/gdata.h"
#include "gambfix/gambOw.h"
#include "ExportHeaders/ExportAMB.h"

using namespace gnut;

namespace great
{
	/**
	*@brief	 Class for storaging ambinp file data
	*/
	class AMB_LIBRARY_EXPORT t_gambcon : public t_gdata
	{
	public:

		/** @brief default constructor. */
		t_gambcon();

		/** @brief default destructor. */
		virtual ~t_gambcon();

		void set_header(t_gtime beg, double intv, int num_amb, set<string>& sats, set<string>& recs);
		void add_ddamb(t_gamb_Dd& ddamb);
		vector<t_gamb_Dd>  get_Dds();
		vector<t_gamb_Dd>& get_Dds(const t_gtime& time);

	public:
		map < pair<GSYS, t_gamb::AMB_ID>, int> all_WL_num;
		map < pair<GSYS, t_gamb::AMB_ID>, int> fix_WL_num;
		map < pair<GSYS, t_gamb::AMB_ID>, int> all_NL_num;
		map < pair<GSYS, t_gamb::AMB_ID>, int> fix_NL_num;

	protected:
		vector<t_gamb_Dd>	            _ddamb;
		map<t_gtime, vector<t_gamb_Dd>> _nddamb;
		vector<t_gamb_Sd>	            _sdamb;
		map<t_gtime, vector<t_gamb_Sd>> _nsdamb;
	private:
		t_gtime		_beg;
		double		_interval;
		int			_num_fixamb;
		set<string> _sats;
		set<string> _sites;
	};

}
#endif // !GAMBCON_H