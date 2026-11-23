
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "gutils/gtypeconv.h"
#include "gall/gallprec.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallprec::t_gallprec() : t_gallnav(),
   _degree_sp3(9),
   _sec(3600.0*6),
   _ref(t_gtime::GPS),
   _clkref(t_gtime::GPS),
   _clkrnx(true),
   _clksp3(false),
   _clknav(false),
   _posnav(false)
{
  gtrace("t_gallprec::t_gallprec");

  id_type(  t_gdata::ALLPREC );
  id_group( t_gdata::GRP_EPHEM );

  _com = 1;
}


// destructor
// ----------
t_gallprec::~t_gallprec()
{
   gtrace("t_gallprec::~t_gallprec");
   
  _mapprec.clear();
}





// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// ----------
int t_gallprec::pos( string sat, const t_gtime& t, double xyz[], double var[], double vel[], bool chk_mask )
{       
   gtrace("t_gallprec::pos");   


  _gmutex.lock();


  shared_ptr<t_geph> tmp = t_gallprec::_find( sat, t );
  if( tmp == _null ){
    for(int i = 0; i<3; i++){
                xyz[i] = 0.0;
      if( var ) var[i] = 0.0;
      if( vel ) vel[i] = 0.0;
    }     
    _gmutex.unlock();
    if( _posnav && t_gallnav::pos( sat, t, xyz, var, vel, chk_mask ) >= 0 ){ return 1; } // cout << "USING POS NAV:\n"; return 1; }
    return -1;
  }

  int irc = tmp->pos( t, xyz, var, vel, _chk_health && chk_mask );
  _gmutex.unlock(); return irc;
}



int t_gallprec::clk_cdr(string sat, const t_gtime& t, double* clk, double* var, double* dclk, double* ifcb, bool chk_mask)
{
    if (_mapclk.find(sat) == _mapclk.end()) return -1;

    map<t_gtime, t_map_dat>::iterator itReq = _mapclk[sat].lower_bound(t); // 1st equal|greater [than t]   

    if (itReq == _mapclk[sat].end()) return -1;

    if (itReq->first > t) return -1;
    else
    {
        *clk  = itReq->second["C0"];
        *var  = 0.0;
        *dclk = 0.0;
        *ifcb = itReq->second["IFCB_F3"];
    }
    return 0;
}


void t_gallprec::add_interval(int intv)
{
	_intv = intv;
}


void t_gallprec::add_interval(string sat,int intv)
{
	_intvm[sat] = intv;
}


void t_gallprec::add_agency( string agency)
{
	_agency = agency;
}



int t_gallprec::intv()
{
	return _intv;

}

int t_gallprec::intv(string sat)
{
	return _intvm[sat];

}

string t_gallprec::get_agency()
{
	return _agency;

}

t_gtime t_gallprec::get_beg()
{
	return _tbeg;
}

t_gtime t_gallprec::get_end()
{
	return _tend;
}

vector<string> t_gallprec::get_sats()
{
	return _sats;
}

vector<string> t_gallprec::get_sat3()
{
	return _sat3;
}

string t_gallprec::get_data_type()
{
	return _datatype;
}

string t_gallprec::get_sat_type()
{
	return _sattype;
}

int t_gallprec::get_sp3_size()
{
	return _mapsp3.size();
}


void t_gallprec::set_beg(t_gtime beg)
{
	_tbeg = beg;
}

void t_gallprec::set_end(t_gtime end)
{
	_tend = end;
}

void t_gallprec::set_sat(vector<string> sat)
{
	_sats = sat;
}

void t_gallprec::set_sat3(vector<string> sat)
{
	_sat3 = sat;
}

void t_gallprec::set_data_type(string type)
{
	_datatype = type;
}

void t_gallprec::set_sat_type(string type)
{
	_sattype = type;
}

void t_gallprec::set_agency(string agency)
{
	_agency = agency;
}

// ----------
//int t_gallprec::addpos( string sat, const t_gtime& ep, double xyzt[], double dxyz[] )
 int t_gallprec::addpos( string sat, const t_gtime& ep, t_gtriple  xyz, double t,
 			                               t_gtriple dxyz, double dt )
 {
    gtrace("t_gallprec::addpos");   
   

   _gmutex.lock();
  
   if( _overwrite || _mapsp3[sat].find(ep) == _mapsp3[sat].end() ){

     _mapsp3[sat][ep]["X"]  = xyz[0];
     _mapsp3[sat][ep]["Y"]  = xyz[1];
     _mapsp3[sat][ep]["Z"]  = xyz[2];
     _mapsp3[sat][ep]["C"]  = t;
   
     _mapsp3[sat][ep]["SX"] = dxyz[0];
     _mapsp3[sat][ep]["SY"] = dxyz[1];
     _mapsp3[sat][ep]["SZ"] = dxyz[2];
     _mapsp3[sat][ep]["SC"] = dt;
     
     if( _log && _log->verb() >= 4 ){
       ostringstream lg;
       lg << "add sat xyz, t " << sat << fixed << setprecision(6)
          << " " << ep.str("%Y-%m-%d %H:%M:%S[%T] ")
          << " " << setw(16) << _mapsp3[sat][ep]["X"]
 	 << " " << setw(16) << _mapsp3[sat][ep]["Y"]
 	 << " " << setw(16) << _mapsp3[sat][ep]["Z"] << setprecision(9)
 	 << " " << setw(16) << _mapsp3[sat][ep]["C"];
       _log->comment(4,"gallprec",lg.str());
     }
     
   }else{
     if(_log) _log->comment(3,"gallprec",ep.str_ymdhms(sat + " not overwriting position for "));
     _gmutex.unlock(); return -1;
   }  
   _gmutex.unlock(); return 0;
 }


// add velocity
// ----------
int t_gallprec::addvel( string sat, const t_gtime& ep, double xyzt[], double dxyz[] )
{
   gtrace("t_gallprec::addvel");   
   

  _gmutex.lock();
  
  if( _overwrite || _mapsp3[sat].find(ep) == _mapsp3[sat].end() ){

    _mapsp3[sat][ep]["VX"] = xyzt[0];
    _mapsp3[sat][ep]["VY"] = xyzt[1];
    _mapsp3[sat][ep]["VZ"] = xyzt[2];
    _mapsp3[sat][ep]["VC"] = xyzt[3];
   
    _mapsp3[sat][ep]["SVX"] = dxyz[0];
    _mapsp3[sat][ep]["SVY"] = dxyz[1];
    _mapsp3[sat][ep]["SVZ"] = dxyz[2];
    _mapsp3[sat][ep]["SVC"] = dxyz[3];

  }else{
    if( _log && _log->verb() >= 4 ){
      ostringstream lg;
      lg << "add sat dxyz,dt " << sat << fixed << setprecision(6)
         << " " << ep.str("%Y-%m-%d %H:%M:%S[%T] ")
         << " " << setw(16) << _mapsp3[sat][ep]["VX"]
	     << " " << setw(16) << _mapsp3[sat][ep]["VY"]
	     << " " << setw(16) << _mapsp3[sat][ep]["VZ"] << setprecision(9)
	     << " " << setw(16) << _mapsp3[sat][ep]["VT"];
      _log->comment(4,"gallprec",lg.str());
    }
    if( _log && _log->verb() >= 3 ) _log->comment(3,"gallprec",ep.str_ymdhms(sat + " not overwriting velocity for "));
    _gmutex.unlock(); return -1;
  }
  _gmutex.unlock(); return 0;
}

// ----------
int t_gallprec::addclk( string sat, const t_gtime& ep, double clk[], double dxyz[] )
{
   gtrace("t_gallprec::addclk");   
   

  _gmutex.lock();
  
  if( _overwrite || _mapclk[sat].find(ep) == _mapclk[sat].end() ){

	// change by ZhengHJ for slow
	  t_map_dat data = {
		  {"C0",clk[0]},
		  {"C1",clk[1]},
		  {"C2",clk[2]},
		  {"SC0",dxyz[0]},
		  {"SC1",dxyz[1]},
		  {"SC2",dxyz[2]},
		};
	  _mapclk[sat][ep] = data;



	if (sat.length() == 3) {
		_clk_type_list.insert(AS);
	}
	else {
		_clk_type_list.insert(AR);
	}

  }else{
    if( _log && _log->verb() >= 3 ) _log->comment(3,"gallprec",ep.str_ymdhms(sat + " not overwriting clocks for "));
    _gmutex.unlock(); return 1;
  }


  _gmutex.unlock(); return 0;
}

int t_gallprec::addclk_tri(string sat, const t_gtime& ep, double clk[], double dxyz[])
{
    gtrace("t_gallprec::addclk_tri");


    _gmutex.lock();

    if (_overwrite || _mapclk[sat].find(ep) == _mapclk[sat].end()) {

       
        t_map_dat data = {
            {"C0",clk[0]},
            {"C1",clk[1]},
            {"C2",clk[2]},
            {"IFCB_F3",clk[3]},
            {"SC0",dxyz[0]},
            {"SC1",dxyz[1]},
            {"SC2",dxyz[2]},
        };
        _mapclk[sat][ep] = data;
  


        if (sat.length() == 3) {
            _clk_type_list.insert(AS);
        }
        else {
            _clk_type_list.insert(AR);
        }

    }
    else {
        if (_log && _log->verb() >= 3) _log->comment(3, "gallprec", ep.str_ymdhms(sat + " not overwriting clocks for "));
        _gmutex.unlock(); return 1;
    }


    _gmutex.unlock(); return 0;
}




// return first epoch of rinex clocks
// ----------
t_gtime t_gallprec::beg_clk( string prn )
{
   gtrace("t_gallprec::beg_clk");   
   

  _gmutex.lock();

  t_gtime tmp = LAST_TIME;
  if( ! prn.empty() ){
    if( _mapclk.find(prn) != _mapclk.end() ) tmp = _mapclk[prn].begin()->first;
  }else{
    for( auto itSAT = _mapclk.begin(); itSAT != _mapclk.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
	if( _mapclk[itSAT->first].begin()->first < tmp ) tmp = _mapclk[itSAT->first].begin()->first;
      }
    }
  }

  _gmutex.unlock(); return tmp;
}


// return last epoch of rinex clocks
// ----------
t_gtime t_gallprec::end_clk( string prn )
{
   gtrace("t_gallprec::end_clk");   
   

  _gmutex.lock();

  t_gtime tmp = FIRST_TIME;
  if( ! prn.empty() ){
    if( _mapclk.find(prn) != _mapclk.end() ) tmp = _mapclk[prn].rbegin()->first;
  }else{
    for( auto itSAT = _mapclk.begin(); itSAT != _mapclk.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
	if( _mapclk[itSAT->first].rbegin()->first > tmp ) tmp = _mapclk[itSAT->first].rbegin()->first;
      }
    }
  }

  _gmutex.unlock(); return tmp;
}




shared_ptr<t_geph> t_gallprec::_find( string sat, const t_gtime& t )
{
  gtrace("t_gallprec::_find");

  if( _mapsp3[sat].size() == 0 ) return _null; // make_shared<t_geph>();

  // if not exists satellite not in cache
  t_map_sp3::iterator it = _prec.find(sat);
  if( it == _prec.end() ){
    if( _get_crddata( sat, t ) < 0 ) return _null; // make_shared<t_geph>();
  }

  // could not find the data at all --- SHOULDN'T OCCURE, SINCE _get_crddata already return !
  it = _prec.find(sat);
  if( it == _prec.end() ){
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - gephprec element not found for "));
    return _null; //make_shared<t_geph>();
  }   

  double t_minus_ref = t - (it->second)->epoch();


  // standard case: cache - satellite found and cache still valid!
  if( fabs((float)t_minus_ref) < (it->second)->interval()/_degree_sp3 &&  // more frequently updated cache
      (it->second)->valid(t)                                              // internal gephprec validity
  ){

    return it->second;

  // specific case: cache - satellite is not within its standard validity (try to update cache)
  }else{

    t_gtime beg(_mapsp3[sat].begin()->first);
    t_gtime end(_mapsp3[sat].rbegin()->first);
     
    // update cache only if not close to the prec data boundaries
    if( ( fabs(t.diff(beg)) > (it->second)->interval()/2 &&
          fabs(t.diff(end)) > (it->second)->interval()/2 ) ||
       ! (it->second)->valid(t)
    ){

      if( _get_crddata( sat, t ) < 0 ) return _null; //make_shared<t_geph>();
      it = _prec.find(sat);
      if( _log && _log->verb() >= 3 && _log->verb() >= 3 ) _log->comment(4,"gallprec",t.str_ymdhms(sat + " updated cache for "));
       
    }
  }

  return it->second;

}


// fill PT,X,Y,Z vectors
// ----------
int t_gallprec::_get_crddata( string sat, const t_gtime& t )
{
   gtrace("t_gallprec::_get_crddata");   

  _T.clear();
  _PT.clear(); _X.clear(); _Y.clear(); _Z.clear(); 
  _CT.clear(); _C.clear();

  if( _mapsp3.find(sat) == _mapsp3.end() ) return -1;
   
  map<t_gtime,t_map_dat>::iterator itReq = _mapsp3[sat].lower_bound(t); // 1st equal|greater [than t]
   
  if( itReq == _mapsp3[sat].end() ) return -1;

  // add for nearest epoch add by glfeng
  // (upstream decremented begin() when t precedes the first record: undefined behaviour, hang/crash in the SP3 writer)
  if (itReq != _mapsp3[sat].begin())
  {
	  map<t_gtime, t_map_dat>::iterator itReq_tmp = std::prev(itReq);
	  if (abs(t.diff(itReq_tmp->first)) < abs(t.diff(itReq->first)))  itReq = itReq_tmp;
  }
   
  _ref = itReq->first; // get the nearest epoch after t as reference

  map<t_gtime,t_map_dat>::iterator itBeg = _mapsp3[sat].begin();
  map<t_gtime,t_map_dat>::iterator itEnd = _mapsp3[sat].end();
  map<t_gtime,t_map_dat>::iterator it    = itReq;

  if( itReq == itEnd ){

    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - no ephemeris found for "));
    return -1;
  }

  int limit = static_cast<int>(_degree_sp3/2); // round (floor)

  bool beg_limit = false;
  auto temp_iter = itBeg;
  for (int i = 0; i < limit; i++) {
	  if (temp_iter == itReq) {
		  beg_limit = true;
		  break;
	  }
	  temp_iter++;
  }

  bool end_limit = false;
  temp_iter = itEnd;
  for (int i = 0; i < static_cast<int>(_degree_sp3 - limit); i++) {
	  if (temp_iter == itReq) {
		  end_limit = true;
		  break;
	  }
	  temp_iter--;
  }
   
  // too few data
  if(_mapsp3[sat].size() < static_cast<int>(_degree_sp3) ){
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - not enough eph data found for "));
    return -1;

 // =============================  
  // start from the first item
  }else if( beg_limit ){
    it = itBeg;
	
// =============================  m
// start from the last item
  }else if( end_limit ){
    it = itEnd;
	//
    for(int i=0; i <= static_cast<int>(_degree_sp3); i++) it--;
// =============================
  // around requested item (standard case)
  }else{ 
     for(int i = 0; i < limit; i++ ) it--;
  } 

 
  if (it == itEnd){
	  if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no ephemeris found for "));
	  return -1;
  }

  // vector for polynomial
  for(unsigned int i = 0; i <= _degree_sp3; it++, i++ ){
    double tdiff = it->first - _ref;
    
    // check maximum interval allowed between reference and sta/end epochs
    if( fabs(tdiff) > static_cast<double>(_degree_sp3*MAXDIFF_EPH) ) continue;
       
    if( it->second["X"] != UNDEFVAL_POS ){
       
      _PT.push_back( tdiff );
       _T.push_back( it->first );
       _X.push_back( it->second["X"] );
       _Y.push_back( it->second["Y"] );
       _Z.push_back( it->second["Z"] );
      _CT.push_back( tdiff );
       _C.push_back( it->second["C"] );


    }
  }

  if( _X.size() != _degree_sp3 + 1 ){
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - not enough eph data for "));
    return -1;
  }

  if( _prec.find(sat) != _prec.end() ){
    _prec[sat]->degree(_degree_sp3);
    _prec[sat]->add( sat, _T, _X, _Y, _Z, _C );
    if( _log && _log->verb() >= 3 ) _log->comment(3,"gallprec",_prec[sat]->epoch().str_ymdhms(sat + " updating cache for "));
  }else{
    shared_ptr<t_gephprec> tmp(new t_gephprec); if( _log ) tmp->glog(_log);
    tmp->degree(_degree_sp3);
    tmp->add( sat, _T, _X, _Y, _Z, _C );
     _prec[sat] = tmp;
    if( _log && _log->verb() >= 3 ) _log->comment(3,"gallprec",_prec[sat]->epoch().str_ymdhms(sat + " creating cache for "));
  }

  return 1;
}


// fill CT,C vectors
// ----------
int t_gallprec::_get_clkdata( string sat, const t_gtime& t )
{
   gtrace("t_gallprec::_get_clkdata");   
 

   _CT.clear(); _C.clear(); _IFCB_F3.clear();
   
  if( _mapclk.find(sat) == _mapclk.end() ) return -1;
  map<t_gtime,t_map_dat>::iterator itBeg  = _mapclk[sat].begin();
  map<t_gtime,t_map_dat>::iterator itEnd  = _mapclk[sat].end();
  map<t_gtime,t_map_dat>::iterator itReq  = _mapclk[sat].lower_bound(t); // 1st equal|greater [than t]   

  if( itReq == _mapclk[sat].end() ) return -1;    // too old products
  if( t < itBeg->first ) return -1;               // too new products

  _clkref = itReq->first;                         // get the nearest epoch after t as reference

  map<t_gtime,t_map_dat>::iterator it     = itReq;

  if( itReq == itEnd ){
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - no clock found for "));
    return -1;
  }


   
  unsigned int degree_clk = 1; // MOZNA DAT 3 degree polinom kvuli 1Hz datum (ale cekovat ktery 3)
   

   
  // DISTANCE() NEED TO BE A POSITIVE DIFFERENCE !!!
  int limit = static_cast<int>(degree_clk/2); // round (floor)


  bool flag_left = false;
  auto itleft = itReq;
  for (int i = 0; i <= limit; i++) {
	  if (itleft == itBeg) {
		  flag_left = true;
		  break;
	  }
	  --itleft;
  }
  bool flag_right = false;
  auto itright = itReq;
  for (int i = 0; i < static_cast<int>(degree_clk - limit); i++) {
	  ++itright;
	  if (itright == itEnd) {
		  flag_right = true;
		  break;
	  }
  }

  if (_mapclk[sat].size() <= static_cast<size_t>(degree_clk)) {
    // Allow single-record polynomial clocks when C1/C2 derivatives are valid; otherwise need degree_clk + 1 points
    if (it->second["C1"] >= UNDEFVAL_CLK) {
      if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - not enough clk data found for "));
      return -1;
    }
    it = itReq;
  }else if (flag_left) {
    it = itBeg;
  }else if (flag_right) {
    it = std::prev(itEnd, degree_clk + 1);
  // around requested item (standard case)
  }else{
     for(int i = 0; i <= limit; i++ ) --it;
  }

  // calculate
  if( it->second["C1"] < UNDEFVAL_CLK ){
    if( it->second["C2"] < UNDEFVAL_CLK ){
        // xiongyun:add for ifcb
        if (it->second.find("IFCB_F3") == it->second.end())
        {
           
            double tdiff = it->first - t;  
            _CT.push_back(tdiff);
            _C.push_back(it->second["C0"] + it->second["C1"] * tdiff + it->second["C2"] * tdiff * tdiff);
        }
        else
        {
            auto it_next = std::next(it);
            if (it_next == itEnd) return -1;
            double tdiff = it_next->first - t;
            _CT.push_back(tdiff);
            _C.push_back(it_next->second["C0"]);
            _IFCB_F3.push_back(it_next->second["IFCB_F3"]);
        }
      return 1;
    }else{
      double tdiff = it->first - t;  
      _CT.push_back( tdiff );
      _C.push_back( it->second["C0"] + it->second["C1"]*tdiff);
      return 1;
    }

  // interpolate
  }else{	   
   
    // vector for polynomial
    for(unsigned int i = 0; i <= degree_clk; it++, i++ ){
      double tdiff = it->first - _clkref;
    
      // check maximum interval allowed between reference and sta/end epochs
      if( fabs(tdiff) > static_cast<double>(degree_clk*MAXDIFF_CLK) ){
	continue;
      }
       
      if( it->second["C0"] != UNDEFVAL_CLK ){
        _CT.push_back( tdiff );
         _C.push_back( it->second["C0"] );
	 


      }
    }

    if( _C.size() != degree_clk + 1 ){
		
		if (itReq != itBeg) {
			_C.clear();
			_CT.clear();
			// Add First
			auto itFirst = itReq;
			itFirst--;
			_C.push_back(itFirst->second["C0"]);
			_CT.push_back(itFirst->first - _clkref);

			// Add Second
			_C.push_back(itReq->second["C0"]);
			_CT.push_back(0.0);

			return 1;
		}
		
      if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - not enough clk data found for "));
      return -1;
    }
  }

  return 1;
}

} // namespace
