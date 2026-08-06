
#ifndef MUTEX_H
#define MUTEX_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: mutual exclusion
  Version: $ Rev: $

  2013-08-14 /JD: created

-*/

#if defined _WIN32 || defined _WIN64
#include <windows.h>
#else
#include <pthread.h>
#endif

#include "ExportHeaders/ExportBase.h"
#include <omp.h>
#include <thread>
#include <mutex>

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_gmutex
	{
	public:
		t_gmutex();
		t_gmutex(const t_gmutex& Other);
		~t_gmutex();

		t_gmutex operator=(const t_gmutex& Other);

		void lock();
		void unlock();


	protected:
#ifdef USE_OPENMP
		omp_lock_t _mutex;
#else
		mutex _mutex;

#endif // USE_OPENMP


	};



} // namespace

#endif // MUTEX_H
