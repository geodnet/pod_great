/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

/* paths:
 *           file://dir/file
 *           http://user:pasw@host:port
 *            ftp://user:pasw@host:port
 *            tcp://host:port
 *          ntrip://mount/user:pasw@host:port
-*/

#include <cstring>
#include <iostream>

#include "gio/gio.h"
#include "gutils/gcommon.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// constructor
// ---------
t_gio::t_gio()
 : _log(0),
   _fd(-1),
   _size(BUF_SIZE),
   _path(""),
   _giof(),
   _count(0),
   _verb(0),
   _stop(0),
   _opened(0),
   _running(0) 
{
  gtrace("t_gio::construct");

  _coder = 0;
}


// destructor
// ---------
t_gio::~t_gio()
{
  gtrace("t_gio::destruct");
}


// set path
// ---------
int t_gio::path( string str )
{
  gtrace("t_gio::path(str)");

  if( str == "" ) return -1;
   
  _path = str;
  return 1;
}


// overloading << operator
// ---------
ostream& operator<<(ostream& os, const t_gio& n)
{
  os << n.path();
  return os;
}



// start reading (could be run in a separate thread)
// ---------
void t_gio::run_read()
{ 
  gtrace("t_gio::run_read");
  cerr << "Reading " << this->_path << endl;


  _gmutex.lock();

  // be sure always clean decoder
  if( _coder ){ _coder->clear(); } 
  // For error message
  vector<string> errmsg;
   
  char* loc_buff = new char[_size];

  if( init_read() < 0 ){
    if( _log ) _log->comment(0,"gio"," warning - initialization failed");
    t_gtime::gmsleep(3000); // [ms]
  }
  int nbytes  = 0;
  int decoded = 0;
  _stop = 0;
  _running = 1;


  while( (( nbytes = _gio_read( loc_buff, _size ) ) > 0 )
	  && _stop != 1 ){
     


    // archive the stream
    _locf_write(loc_buff,nbytes);

    if( _coder && nbytes > 0 ){
       decoded = _coder->decode_data( loc_buff, nbytes, _count, errmsg );
	   if (_coder->end_epoch > t_gtime(0, 0))
	   {
		   if (_coder->epoch > _coder->end_epoch)
		   {
			   break;
		   }
	   }
    }
    else {
       if( _log ) _log->comment(2,"gio","0 data decoded");       
       decoded = 0;
    }

    if( _log && _log->verb() >= 9) _log->comment(9,"gio","read: " + int2str(nbytes)
			                           + " decoded: " + int2str(decoded));
  }

  _stop_common(); 
  delete [] loc_buff;
  _gmutex.unlock(); return;
}


// start reading (could be run in a separate thread)
// ---------
void t_gio::run_write()
{
  gtrace("t_gio::run_write");

  _gmutex.lock();

  if( _coder ){ _coder->clear(); } // be sure always clean decoder
	
  vector<string>  errmsg;
  _stop = 0;

  char* loc_buff = new char[_size];
   
  if (_opened != 1 && init_write() < 0) {
	  _stop_common();
	  if (_log) _log->comment(0, "gio", " warning - initialization failed");
	  delete[] loc_buff;
	  _gmutex.unlock(); return;
  }

  int nbytes = 0;
  _running = 1;
  do{
     nbytes = 0;
     

     nbytes = _locf_read(loc_buff,_size);



     // 3. read from encoder
     if( nbytes < 0 && _coder ){
       nbytes = _coder->encode_data( loc_buff, _size, _count, errmsg );
     }

  } while (_stop == 0 && (_gio_write(loc_buff, nbytes) > 0 || nbytes < -1));

  _stop_common();
  if( _log ) _log->comment(2,"gio","end of read [stop or gio_write problem]");
  delete [] loc_buff;
  _gmutex.unlock(); return;
}


// local log file source
// ---------
int t_gio::_locf_read( char* buff, int size )
{   
  gtrace("t_gio::_locf_read");
  
  if( _giof.mask() == "" ) return -1;

  return _giof.read(buff,size);
}


 // local log file archive
// ---------
int t_gio::_locf_write( const char* buff, int size )
{ 
  gtrace("t_gio::_locf_write");

  if( _giof.mask() == "" ) return -1;

  int irc = _giof.write( buff, size );

  return irc;
}


// common function for initialization
// ---------
int t_gio::_init_common()
{ 
  gtrace("t_gio::_init_common");

  return 1; 
}



// common function for socket/file close
// ---------
int t_gio::_stop_common()
{ 
  gtrace("t_gio::_stop_common");

  _opened  = 0;  // 
 _running = 0;
  return _running;
}
   
} // namespace