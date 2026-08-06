/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-05-15 created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    t_gbiasemodel.h
* @brief	    mainly about how to cacultae B P l single
* @mainpage     C++ version 
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-05-18
*
*/

#ifndef GBIASEMODEL_H
#define GBIASEMODEL_H

#include "gutils/gtime.h"
#include "gall/gallpar.h"
#include "gall/gallobs.h"
#include "gdata/gsatdata.h"
#include "gdata/gobsgnss.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/gbasemodel.h"
#include "gproc/glsqmatrix.h"

using namespace std;
using namespace gnut;

namespace great
{
    class BASE_LIBRARY_EXPORT t_gbiasmodel
    {
    public:
        t_gbiasmodel();
        virtual ~t_gbiasmodel();

        virtual bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gobs& gobs,t_gbaseEquation& result)=0;
        virtual void update_obj_clk(const string& obj, const t_gtime& epo, double clk) = 0;
        virtual double get_rec_clk(const string& obj) = 0; 

        virtual t_gmodel* precisemodel() { return nullptr; }
    };

    class BASE_LIBRARY_EXPORT t_gprecisebias:public t_gbiasmodel
    {
    public:
        t_gprecisebias(t_gallproc *data, t_glog *log, t_gsetbase *setting);
        ~t_gprecisebias();

		void set_multi_thread(const set<string>& sites);

        bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gobs& gobs,t_gbaseEquation& result) override;
        void update_obj_clk(const string& obj, const t_gtime& epo, double clk) override;
        double get_rec_clk(const string& obj) override;

        t_gmodel* precisemodel() { return &gmodel; };
    protected:
        t_gprecisemodel gmodel;
		tuple<string, string, t_gtime> _rec_sat_before;

		bool _multi_thread_flag;
		map<string, shared_ptr<t_gprecisemodel> > _map_site_model;
		map < string, tuple<string, string, t_gtime> > _map_flag;
		t_gmutex _map_mtx;

    private:
        t_glog *_log;
        t_gsetbase* _setting;
    };


} // namespace great

#endif