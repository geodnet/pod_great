/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-12-23  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gupdatepar.h
* @brief		define update parameter IF model
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-12-23
*
*/

#ifndef GUPDATEPARIF_H
#define GUPDATEPARIF_H

#include <map>
#include "gproc/gupdatepar.h"

namespace great {
	class BASE_LIBRARY_EXPORT t_gupdateparIF :public t_gupdatepar
	{
	public:
		t_gupdateparIF(shared_ptr<t_gcycleslip> cycleslip);
		virtual ~t_gupdateparIF();

	protected:

		virtual void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info) override;

		void _update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata,t_gpar::t_type amb_type,t_gcycleslip* _cycleslip,t_gupdateparinfo& update_info);
		void _update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gpar::t_type amb_type, t_gupdateparinfo& update_info);
	
	protected:
		map<pair<string, string>, int> _amb_id;

	};

}
#endif // !GUPDATEPARIF_H
