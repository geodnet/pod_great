
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 -*/

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include "gall/gallobs.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallobs::t_gallobs() 
  : t_gdata(),
    _set(0),
    _nepoch(0),
    _overwrite(false)
{
  gtrace("t_gallobs::constructor");
  id_type(  t_gdata::ALLOBS );
  id_group( t_gdata::GRP_OBSERV );
}


// destructor
// ----------
t_gallobs::~t_gallobs()
{
  gtrace("t_gallobs::destructor");
  _mapobj.clear();
  _filter.clear();
}


// settings
// ----------
void t_gallobs::gset(t_gsetbase* gset)
{  
  _set = gset;
  _sys = dynamic_cast<t_gsetgen*>(_set)->sys();
  _smp = dynamic_cast<t_gsetgen*>(_set)->sampling();
  _scl = dynamic_cast<t_gsetgen*>(_set)->sampling_scalefc(); // scaling 10^decimal-digits



  return;
}


 
   
// protected
// ---------
set<GSYS> t_gallobs::_gsys( const string& site, const t_gtime& t)
{
  set<GSYS> gnss;

  vector<t_gsatdata> satdata = _gobs(site, t);
 
  vector<t_gsatdata>::iterator it;
  for (it = satdata.begin(); it != satdata.end(); it++){
    GSYS gsys = it->gsys();
    gnss.insert(gsys);
  }

  return gnss;
}
   

// return list of available stations
// ----------
set<string> t_gallobs::stations()
{
  gtrace("t_gallobs::stations");


  _gmutex.lock();

  set<string> all_sites;
  t_map_oobj::const_iterator itSITE = _mapobj.begin();
   
  while( itSITE != _mapobj.end() ){
    all_sites.insert(itSITE->first);
    itSITE++;
  }
  _gmutex.unlock(); return all_sites;
}

// return list of available satellites
// ----------
set<string> t_gallobs::_sats(const string& site,
                             const t_gtime& t, GSYS gnss)
{
  set<string> all_sats;

  if( _mapobj.find(site)    == _mapobj.end() ||
      _mapobj[site].find(t) == _mapobj[site].end() ){ return all_sats; }

  t_map_osat::const_iterator itSAT = _mapobj[site][t].begin();  
  while( itSAT != _mapobj[site][t].end() ){
    GSYS sys = itSAT->second->gsys();
    if(gnss == GNS) all_sats.insert(itSAT->first);
    else if(gnss == GPS && sys == GPS) all_sats.insert(itSAT->first);
    else if(gnss == GLO && sys == GLO) all_sats.insert(itSAT->first);
    else if(gnss == GAL && sys == GAL) all_sats.insert(itSAT->first);
    else if(gnss == BDS && sys == BDS) all_sats.insert(itSAT->first);
    else if(gnss == QZS && sys == QZS) all_sats.insert(itSAT->first);
    else if(gnss == SBS && sys == SBS) all_sats.insert(itSAT->first);     
    itSAT++;
  }
  return all_sats;
}

vector<t_gobsgnss*> t_gallobs::obs_lockfree(bool isPtr, const string& site, const t_gtime& t)
{

    return _gobs(true, site, t);
}

// return list of available satellites for all sites
// ----------
vector<t_gsatdata> t_gallobs::obs(const set<string>& sites, const t_gtime& t)
{
	gtrace("t_gallobs::obs");


	_gmutex.lock();

	
	vector<t_gsatdata> all_obs;
	// cerr << sites.size() << endl;
	for (auto site_iter = sites.begin(); site_iter != sites.end(); site_iter++)
	{
		vector<t_gsatdata> site_obs = _gobs(*site_iter, t);
		all_obs.insert(all_obs.end(),site_obs.begin(), site_obs.end());
	}


	_gmutex.unlock(); return all_obs;

}


// protected
// ---------
vector<t_gsatdata> t_gallobs::_gobs(const string& site,
                                    const t_gtime& t)
{
  vector<t_gsatdata> all_obs;
  t_gtime tt(t_gtime::GPS);

  if (_find_epo(site, t, tt) < 0) {
     return all_obs;
  }

  t_map_osat::iterator itSAT = _mapobj[site][tt].begin();
  while( itSAT != _mapobj[site][tt].end() ){    
//    all_obs.push_back( itSAT->second ); 
    all_obs.push_back( *itSAT->second );
    itSAT++;
  }

  return all_obs;
}

vector<t_gobsgnss*> t_gallobs::_gobs(bool isPtr, const string& site,
    const t_gtime& t)
{
    vector<t_gobsgnss*> all_obs;
    t_gtime tt(t_gtime::GPS);

    if (_find_epo(site, t, tt) < 0) {
        return all_obs;
    }
    t_map_osat::iterator itSAT = _mapobj[site][tt].begin();
    while (itSAT != _mapobj[site][tt].end()) {
        all_obs.push_back((*itSAT).second.get());
        itSAT++;
    }

    return all_obs;
}


// TEMPORARY !!!
// return list of available satellites (POINTERS!)
// ----------
vector<t_spt_gobs> t_gallobs::obs_pt(const string& site, const t_gtime& t)
{
  gtrace("t_gallobs::obs_pt");


  _gmutex.lock();

  vector<t_spt_gobs> all_obs;
  t_gtime tt(t_gtime::GPS);
   
  if( _find_epo(site, t, tt) < 0 ){
    _gmutex.unlock(); return all_obs;
  }
               
  t_map_osat::iterator itSAT = _mapobj[site][tt].begin();
  while( itSAT != _mapobj[site][tt].end() ){
     
    string sat = itSAT->first;
     
    // TESTING NEW METHOD
    all_obs.push_back( dynamic_pointer_cast<t_gobsgnss>(itSAT->second) );
 
    itSAT++;
  }

  _gmutex.unlock(); return all_obs;
}



// return list of available satellites
// ----------
vector<t_gtime> t_gallobs::epochs(const string& site)
{
  gtrace("t_gallobs::epochs");

  _gmutex.lock();

  vector<t_gtime> all_epochs;

  if( _mapobj.find(site) == _mapobj.end() ){
    // cout << "site not found: " << site << t.str_ymdhms() << endl; cout.flush();
     _gmutex.unlock(); return all_epochs;
  }
     
  t_map_oref::iterator it = _mapobj[site].begin();

  while( it != _mapobj[site].end() ){
    all_epochs.push_back( it->first );
    it++;
  }
  
  _gmutex.unlock(); return all_epochs;
}


// return first position for satellite
// ----------
t_gtime t_gallobs::beg_obs(const string& site, double smpl)
{
  gtrace("t_gallobs::beg_obs");

  _gmutex.lock();

  t_gtime tmp = LAST_TIME;

  if( _mapobj.find(site)    != _mapobj.end() &&
      _mapobj[site].begin() != _mapobj[site].end() ) tmp = _mapobj[site].begin()->first;

  // get first synchronized obs
  if( smpl > 0.0 ){
    auto itEpoB = _mapobj[site].begin();
    auto itEpoE = _mapobj[site].end();
    int sod = static_cast<int>(dround(itEpoB->first.sod()+itEpoB->first.dsec()));

    while( sod % static_cast<int>(smpl) != 0 && ++itEpoB != itEpoE ){
      sod = static_cast<int>(dround(itEpoB->first.sod()+itEpoB->first.dsec()));
      tmp = itEpoB->first;
      tmp.reset_sod();
      tmp.add_secs(sod); 

    }
  }
   
  _gmutex.unlock(); return tmp;
}


// return last position for satellite
// ----------
t_gtime t_gallobs::end_obs(const string& site)
{
  _gmutex.lock();

  t_gtime tmp = FIRST_TIME;
  if( _mapobj.find(site)    != _mapobj.end() &&
      _mapobj[site].begin() != _mapobj[site].end() )
    tmp = _mapobj[site].rbegin()->first;

  _gmutex.unlock(); return tmp;
}


// add observations ( both P and L in meters !!!)
// ----------
int t_gallobs::addobs(t_spt_gobs obs)
{
  gtrace("t_gallobs::addobs");
   
  // **TEST** 20200924
  // repair small out-sync ( < 5 ms )
  double outsync = fmod(obs->epoch().dsec(), _smp) - round(fmod(obs->epoch().dsec(), _smp));
  if (fabs(outsync) < 0.0051 && fabs(outsync) > 1e-6){
	  obs->epo(obs->epoch() - outsync);
	  vector<GOBS> v_obs = obs->obs();
	  vector<GOBS>::iterator itOBS = v_obs.begin();
	  for (; itOBS != v_obs.end(); ++itOBS){
		  GOBSTYPE obstype = str2gobstype(gobs2str(*itOBS));
		  if (obstype == TYPE_P || obstype == TYPE_C)
			  obs->resetobs(*itOBS, obs->getobs(*itOBS) - CLIGHT * outsync);
		  else if (obstype == TYPE_L)
			  obs->resetobs(*itOBS, obs->getobs(*itOBS) - CLIGHT * outsync/obs->wavelength(str2gobsband(gobs2str(*itOBS))));
          else if(obstype == TYPE_S) // add wh
              obs->resetobs(*itOBS, obs->getobs(*itOBS));
	  }
  }

  t_gtime t(obs->epoch()),   tt = t;
  string site = obs->site();
  string sat  = obs->sat();


  if (_map_sites.find(site) == _map_sites.end()) _map_sites.insert(site);

  _gmutex.lock();

  int  epo_found = _find_epo(site, t, tt);
  auto itSAT     = _mapobj[site][tt].find(sat);


  // ===================================
  if( _overwrite || epo_found < 0                     // epoch exists (smart search)
                 || itSAT == _mapobj[site][tt].end()  // satellite exists

  ){

    // delete old if exists
    // ====================
    if( _overwrite && epo_found > 0                     // epoch exists (smart search)
                   && itSAT != _mapobj[site][tt].end()  // satellite exists    

    ){
        _mapobj[site][tt].erase(sat);
        if( _log && _log->verb() >= 2 ) _log->comment(2, "gallobs", site + tt.str_ymdhms(" obs replaced ") + sat);
    }

    // too many epochs (remove old)
    // ============================
    if( _nepoch > 0 && _mapobj[site].size() > _nepoch + 10 ){  // +10 .. reduce number of removals to every tenths

      auto itEPO = _mapobj[site].begin();
      while( itEPO != _mapobj[site].end() ){

        if( _mapobj[site].size() <= _nepoch ){ break; }
	
        t_gtime t = itEPO->first;
        if( _log && _log->verb() >= 2 )
        {
           string satells(" ");
           for(auto itSAT = _mapobj[site][t].begin(); itSAT != _mapobj[site][t].end(); ++itSAT ){
             satells = satells + " " + itSAT->first;
           }
       
           _log->comment(2, "gallobs", site + t.str_ymdhms(" obs removed ") + satells);
        }

        _mapobj[site][t].erase(_mapobj[site][t].begin(),_mapobj[site][t].end());
        _mapobj[site].erase(itEPO++);
      }
    }

    // distinguish GNSS specific message
    // =================================
    if( obs->id_type() == t_gdata::OBSGNSS ){

      _mapobj[site][tt][sat] = obs;

    }else{
       if( _log ){  _log->comment(0,"gallobs","warning: t_gobsgnss record not identified!"); }
       _gmutex.unlock(); return 1;
    }
  }else{
    if( _log && _log->verb() >= 2 ) _log->comment(2,"gallobs", site + tt.str_ymdhms(" skipped ") + " " + sat);
    _gmutex.unlock(); return 0;
  }

  // comments
  if( _log && _log->verb() >= 3 ){

    vector<GOBS> v_obs = obs->obs();
    vector<GOBS>::iterator itOBS = v_obs.begin();
    ostringstream lg; lg << fixed << setprecision(3);
    for( ; itOBS != v_obs.end(); ++itOBS ) 
      lg << " " << gobs2str( *itOBS ) << ":" << obs->getobs( *itOBS );

    _log->comment(3,"gallobs",site + tt.str_ymdhms(" add obs: ") + " " + sat + lg.str() );
  }
  _gmutex.unlock(); return 0;
}



// clean function
// ----------
void t_gallobs::clean_outer( const  string& obj,
                             const t_gtime& beg, 
                             const t_gtime& end )
{
  gtrace("t_gallobs::clean_outer");

  if( end < beg ) return;

  _gmutex.lock();

  // loop over all object [default]
  t_map_oobj::const_iterator itOBJ  = _mapobj.begin();
  t_map_oobj::const_iterator itOBJ2 = _mapobj.end();
    
  if( _log ){
    _log->comment(2,"gallobs","obs clean request: " + obj + beg.str_ymdhms(": ") + end.str_ymdhms(" - "));
  }

  // loop over a single object (if obj defined !)
  if( ! obj.empty() ) itOBJ = itOBJ2 = _mapobj.find(obj);
  if( itOBJ2 != _mapobj.end() ) itOBJ2++;

  for( ; itOBJ != itOBJ2; ++itOBJ ){

    string site = itOBJ->first;

    // find and CLEAN all data (epochs) out of the specified period !
    t_map_oref::iterator it;
    t_map_oref::iterator itFirst = _mapobj[site].begin();
    t_map_oref::iterator itLast  = _mapobj[site].end();
    t_map_oref::iterator itBeg   = _mapobj[site].lower_bound(beg);  // greater|equal
    t_map_oref::iterator itEnd   = _mapobj[site].upper_bound(end);  // greater only!

    t_map_osat::iterator itOBS;

    // delete before BEGIN request
    if (itBeg != itLast&& itBeg != itFirst){
      if( _log && _log->verb() >= 3  ){
        string lg("obs removed before: " + itBeg->first.str_ymdhms() + " " + site);
        _log->comment(3,"gallobs",lg);
      }

      for( it = itFirst; it != itBeg; ++it ){	 
        t_gtime tt(it->first);
	 
        _mapobj[site][tt].erase(_mapobj[site][tt].begin(),_mapobj[site][tt].end());
      }
       
      _mapobj[site].erase(itFirst,itBeg); // erase all but itBeg
    }

    // delete after END request
    if( itEnd != itLast ){ // && ++itEnd != itLast ){
      if( _log && _log->verb() >= 3 ){
        string lg("obs removed after : " + itEnd->first.str_ymdhms() + " " + site);
        _log->comment(3,"gallobs",lg);
      }

      for( it = itEnd; it != itLast; ++it ){
	    
        t_gtime tt(it->first);
	 
        _mapobj[site][tt].erase(_mapobj[site][tt].begin(),_mapobj[site][tt].end());	 
      }                   
      _mapobj[site].erase(itEnd,itLast);
    }
//    itOBJ++; // WHILE LOOP ONLY
  }
  _gmutex.unlock(); return;
}



// return list of available satellites
// ----------
t_gallobs::t_map_osat t_gallobs::find( string site, const t_gtime& t)
{
  gtrace("t_gallobs::find");

  _gmutex.lock();

  t_gtime tt(t_gtime::GPS);
  t_map_osat  tmp;
  if(_find_epo(site, t, tt) < 0){
     _gmutex.unlock(); return tmp;
  }

//  shared_ptr<t_map_osat> tmp = &_mapobj[site][tt];
  _gmutex.unlock(); return _mapobj[site][tt]; // tmp;
}

// add site-specific filtered data/epochs
// ----------
void t_gallobs::xdata(string site, string file, t_xfilter xflt)
{
  gtrace("t_gallobs::xdata");

  _gmutex.lock();
  
  _filter[site][file].xdat = xflt.xdat;
  _filter[site][file].beg  = xflt.beg;
  _filter[site][file].end  = xflt.end;

  _gmutex.unlock();
}


// get all t_gsatdata for prn from epoch interval
// ----------------------------------------------
vector<t_spt_gobs> t_gallobs::obs_prn_pt(const string& site, const string& prn, 
			  	 	 const t_gtime& beg, const t_gtime& end)
{
  gtrace("t_gallobs::obs_prn_pt");


  _gmutex.lock();
 
  vector<t_spt_gobs> all_obs;
  t_gtime tt(t_gtime::GPS);

  if( _mapobj.find(site) == _mapobj.end() ){

    _gmutex.unlock(); return all_obs;
  }       

  t_map_oref::iterator it1 = _mapobj[site].lower_bound(beg);  // greater || equal
  t_map_oref::iterator it2 = _mapobj[site].lower_bound(end);  // greater || equal   

  for ( ;it1 != it2; ++it1)
  {
     tt = it1->first;
     if( _mapobj[site].find(tt) != _mapobj[site].end() )
     {
       auto itSAT = _mapobj[site][tt].find(prn);
       if( itSAT != _mapobj[site][tt].end() ) all_obs.push_back( itSAT->second );
     }
  }
  _gmutex.unlock(); return all_obs;
}


// find epoch from the map w.r.t. DIFF_SEC
// ---------------------------------------
int t_gallobs::_find_epo(const string& site, const t_gtime& epo, t_gtime& tt)
{
  gtrace("t_gallobs::_find_epo");

  if( _mapobj.find(site) == _mapobj.end() ){
     if( _log && _log->verb() >= 2 ){
        _log->comment(3,"gallobs","site not found: " + site + " " + epo.str_ymdhms());
      }
      return -1;
  }
   
  t_map_oref::iterator it1 = _mapobj[site].lower_bound(epo);  // greater || equal
  t_map_oref::iterator it0 = it1;                             // previous value

  if( it0 != _mapobj[site].begin() ) it0--;                               // set previous value  
  if( it1 == _mapobj[site].end() ) it1 = it0;

  if( it1 == _mapobj[site].end() && it0 == _mapobj[site].end() ){
    cout << site << " " << epo.str_ymdhms() << " " << epo.dsec() << " no observations !\n";
    return -1;
  }
     


  if(      fabs(it1->first - epo) <= DIFF_SEC(_smp) ) tt = it1->first;  // it = it1;                 // set closest value
  else if( fabs(it0->first - epo) <= DIFF_SEC(_smp) ) tt = it0->first;  // it = it0;                 // set closest value
  else{

     return -1;                                     // not found !
  }
   
   return 1;
}

   
} // namespace
