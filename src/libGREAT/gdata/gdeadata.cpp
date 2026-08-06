

#include"gdeadata.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gdeadata::t_gdeadata() :t_gdata()
	{
		id_type(t_gdata::DESAISCOPOLECOEF);
	}


	/**
		* @brief add one oceanide record to _oct.
		* @param[in]   m		index m
		* @param[in]   n		index n
		* @param[in]   A B AI BI	coefficient
	*/
	
	void t_gdeadata::set_rec(const int& n, const int& m, double A, double B, double AI, double BI)
	{

		_gmutex.lock();
		
		_Anm[n][m] = A;
		_Bnm[n][m] = B;
		_AInm[n][m] = AI;
		_BInm[n][m] = BI;
		
		_gmutex.unlock();
		return;
	}

	/**
	* @brief add a load coefficient value to vector.
	* @param[in]   load		load coefficient value
	*/
	void t_gdeadata::set_load( double load)
	{

		_gmutex.lock();
		_load.push_back(load);
		_gmutex.unlock();
		return;
	}

	
}

