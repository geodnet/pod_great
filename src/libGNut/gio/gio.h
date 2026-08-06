
#ifndef GIO_H
#define GIO_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: implements gtime class (day and precise time)
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include <stdio.h>
#include <fstream>
#include <string>

#include "gio/glog.h"
#include "gutils/gmutex.h"
#include "gcoders/gcoder.h"

#define BUF_SIZE 1024

using namespace std;

namespace gnut
{

	class t_gcoder;
	class BASE_LIBRARY_EXPORT t_gio
	{

	public:
		t_gio();
		virtual ~t_gio();

		virtual void run_read();
		virtual void run_write();

		virtual void coder(t_gcoder* coder) { _coder = coder; }

		virtual int init_write() { return _opened = _init_common(); }
		virtual int init_read() { return _opened = _init_common(); }


		virtual int    path(string str);                      // set/get path
		virtual string path()const { return _path; }

	

		void   verb(int i) { _verb = i; }                      // set/get verbosity
		int    verb() { return _verb; }

		

		
		void glog(t_glog* l) { _log = l; }                       // set/get glog pointer
		t_glog* glog() { return _log; }


	protected:
		virtual int _gio_write(const char* buff, int size) = 0;
		virtual int _gio_read(char* buff, int size) = 0;

		virtual int _locf_write(const char* buff, int size);
		virtual int _locf_read(char* buff, int size);

		virtual int _init_common();
		virtual int _stop_common();

		t_glog*         _log;              // log pointer
		int             _fd;               // file descriptor
		size_t          _size;             // buffer size
		string          _path;             // URL-like path (e.g. file:///home/honza/file)
		t_giof          _giof;             // local file
		int             _count;            // record counter
		int             _verb;             // verbosity
		int             _stop;             // require a stop at run() loop
		int             _opened;           // 1: opened/connected
		int             _running;          // running
		t_gcoder*       _coder;            // decoder/encoder
		t_gmutex        _gmutex;



	private:

	};

} // namespace

#endif
