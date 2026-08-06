/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-05-15 created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    t_gcombmodel.h
* @brief	    base combine biase model
* @mainpage     C++ version 
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-05-18
*
*/

#ifndef GCOMBMODEL_H
#define GCOMBMODEL_H

#include "gmodels/gbasemodel.h"
#include "gmodels/gbiasmodel.h"
#include "gmodels/gprecisemodel.h"
#include "gall/gallbias.h"
#include "gall/gallproc.h"
#include "gutils/gcommon.h"
#include <map>
#include <memory>

using namespace gnut;
using namespace std;

namespace great
{
    class BASE_LIBRARY_EXPORT t_gcombmodel:virtual public t_gbasemodel
    {
    public:
        t_gcombmodel(t_gsetbase* setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data);
        virtual ~t_gcombmodel();

		
    protected:


        map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;


        // baise model
        shared_ptr<t_gbiasmodel> _bias_model;

        t_gprecisemodel* _gprecisemodel = nullptr;
        t_gallbias*      _gallbias      = nullptr;
        t_glog*          _log           = nullptr;

        int             _frequency;
        IFB_MODEL       _ifb_model;
		IFCB_MODEL      _ifcb_model;
        IONMODEL        _ion_model;
		
    };

    class BASE_LIBRARY_EXPORT t_gcombIF:virtual public t_gcombmodel
    {
    public:
        t_gcombIF(t_gsetbase *setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data);
        ~t_gcombIF();

        bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gbaseEquation& result) override;
        bool cmb_equ_IF(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,GOBSBAND b1,GOBSBAND b2, t_gbaseEquation& result);

	private:

		map<pair<FREQ_SEQ, GSYS>, t_gpar::t_type> _clk_type_index;

    };


	
}
#endif