
#ifndef GFILE_H
#define GFILE_H

#include <stdio.h>
#include <fstream>
#include <string>
#include "gio/gio.h"
#include "gio/giof.h"
#define FILEBUF_SIZE 20480
#define FILEHDR_SIZE    48
using namespace std;
namespace gnut {

class BASE_LIBRARY_EXPORT t_gfile : public t_gio {
 public:
   t_gfile();
   virtual ~t_gfile();
   virtual int init_write();
   virtual int init_read();
   virtual int irc()const;
   virtual bool eof();
   virtual string md5sum();
   virtual string mask();
   virtual bool compressed(){return false;}
   virtual int path(string str);
   virtual string path();
   virtual string name();
   virtual void reset();
 protected:
   virtual int _gio_write(const char* buff,int size);
   virtual int _gio_read(char* buff,int size);
   virtual int _stop_common();
   virtual void _set_gzip(string name){}
   virtual int _read(char* b,int s);
   virtual int _write(const char* b,int s);
   int _irc;
   t_giof* _file;
};

} // namespace
#endif
