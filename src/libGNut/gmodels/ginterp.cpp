
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "gmodels/ginterp.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_ginterp::t_ginterp()
{}

// destructor
// ----------
t_ginterp::~t_ginterp(){}



//// ----------------------------
//// Linear Interpolation (double)
// ----------------------------
int t_ginterp::linear(map<double, double>& data, double val, double& fval)
{
   map<double, double>::iterator it2, it1;
   it2 = data.lower_bound(val);

   if(it2!=data.end() && double_eq(it2->first, val) ){
      fval = it2->second;
      return 0;
   }
   if( it2 == data.begin() || it2 == data.end() ){
   
      return -1;
   }
   
   it1 = it2;
   it1--;
   
   Matrix x(2,2), y(2,1);
   x << it1->first << 1
     << it2->first << 1;
   y << it1->second 
     << it2->second;
   
   Matrix coef(2,1);
   coef = x.i()*y;
   
   fval = coef(1,1)*val + coef(2,1);

  return 0;
}

// ----------------------------
// bilinear interpolation
// 
// 0 .. 11 (bottom-left)         21 *---------* 22
// 1 .. 12 (bottom-right)           |         |
// 2 .. 21 (top-left)               |         |
// 3 .. 22 (top-right)           11 *---------* 12
// ----------------------------
int t_ginterp::bilinear(const map<t_gpair, double>& data, const t_gpair& req_pos, double& fval)
{
   vector<t_gpair> v_pt;
   map<t_gpair,double>::const_iterator itDAT = data.begin();
   for( ; itDAT != data.end(); ++itDAT){ 
     v_pt.push_back( itDAT->first );      

   }

   // no interpolation (1 point)
   if( v_pt.size() == 1 ){  fval = data.at(v_pt[0]);
      
   // linear interpolation (2 points)
   }else if( v_pt.size() == 2 ){
     double x = 0.0;
     double a = 0.0;
     double b = 0.0;
      
     if( v_pt[1][0] - v_pt[0][0] == 0 )
     {
       x = req_pos[1];  // lat=const
       a = v_pt[0][1];
       b = v_pt[1][1];	
     }else{
       x = req_pos[0];  // lon=const
       a = v_pt[0][0];
       b = v_pt[1][0];
     }
       
     if( a - b == 0 ){
       cout << "# warning: Bilinear interpolation failed (2 points) \n";
       fval = 0.0;
       return -1;
	 
     }else{
       fval = data.at(v_pt[0]) + (data.at(v_pt[1]) - data.at(v_pt[0])) * (x - a) / (a - b);
     }
      
   // bilinear interpolation (4 points)
   }else if( v_pt.size() == 4 ){
         


     double a = data.at(v_pt[0]) * (v_pt[3][0] - req_pos[0]) * (v_pt[3][1] - req_pos[1]);
     double b = data.at(v_pt[1]) * (req_pos[0] - v_pt[2][0]) * (v_pt[2][1] - req_pos[1]);
     double c = data.at(v_pt[2]) * (v_pt[1][0] - req_pos[0]) * (req_pos[1] - v_pt[1][1]);
     double d = data.at(v_pt[3]) * (req_pos[0] - v_pt[0][0]) * (req_pos[1] - v_pt[0][1]);
        
     if( (v_pt[2][0]-v_pt[1][0]) * (v_pt[1][1]-v_pt[0][1]) == 0 )
     {
       fval = 0.0;
       return -1;
     }else{
       fval = ( a+b+c+d ) / ( (v_pt[2][0]-v_pt[1][0]) * (v_pt[1][1]-v_pt[0][1]) );
     }      
   }else{
     fval = 0.0;
     return -1;
   }

   return 0;
}

} // namespace
