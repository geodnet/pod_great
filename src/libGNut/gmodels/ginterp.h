
#ifndef GINTERP_H
#define GINTERP_H 
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: implementation various interpolation tecnique
  Version: $ Rev: $

  2013-06-05 /PV: created

-*/

#include <vector>
#include <map>

#include "newmat/newmat.h"

#include "gutils/gtime.h"
#include "gutils/gpair.h"
#include "gdata/gdata.h"

using namespace std;

namespace gnut {   

// ----------
class BASE_LIBRARY_EXPORT t_ginterp : public t_gdata
{
   public:   
            t_ginterp();
   virtual ~t_ginterp();
   

   
   int linear(        map<double, double>& data, double val, double& fval);

   int bilinear(const map<t_gpair, double>& data, const t_gpair& req_pos, double& fval);

};

} // namespace

#endif
