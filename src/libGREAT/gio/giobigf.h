/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-04-25  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    giobigf.h
* @brief		header files of IO for big file
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/

#ifndef GIOBIGF_H
#define GIOBIGF_H


#include "gio/giof.h"

#define MAX_BUFFER_SIZE 204800000

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_giobigf :public t_giof
	{
	public:

		t_giobigf(string mask="",int buffer_size = 1024*1000*10);
		virtual ~t_giobigf();

		virtual int write(const char* buff, int size) override;
		int flush();
		
		const int buffer_size;

	private:
		char* _buffer;
		int _current;
	};




}


#endif // !GIOBIGF_H