/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-04-25  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    giobigf.cpp
* @brief		header files of IO for big file
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/

#include "giobigf.h"
#include <assert.h>

namespace gnut
{
	t_giobigf::t_giobigf(string mask,int buffer_size)
		:t_giof(mask),
		buffer_size(buffer_size),
		_current(0),
		_buffer(nullptr)
	{
		_buffer = new char[buffer_size];
	}

	t_giobigf::~t_giobigf()
	{
		this->flush();
		if (_buffer) {
			delete[] _buffer;
			_buffer = nullptr;
		}
	}

	int t_giobigf::write(const char * buff, int size)
	{
		if (_current + size > buffer_size) {
			//this->flush();
			if (this->flush() < 0) {
				return -1;
			}
			return (t_giof::write(buff, size));
		}



		memcpy(&this->_buffer[_current], buff, size);
		_current += size;
		return size;
	}

	int t_giobigf::flush()
	{
		if (_current == 0) {
			return 0;
		}

		int state = t_giof::write(_buffer, _current);
		if (state < 0) {
			return state;
		}
		_current = 0;
		return state;
	}


}