/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li	2019-01-11	creat the file.
	  -1.1  Jie Li	2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gsatinfo.cpp
 * @brief		The base class used to storage one satellite information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include "gdata/gsatinfo.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gsatinfo::t_gsatinfo()
	{
		gtrace("t_gsatinfo::t_gsatinfo");
		id_type(t_gdata::SATPARS);
		_mass = 0.0;
	}

	/** @brief destructor. */
	t_gsatinfo:: ~t_gsatinfo()
	{

	}

	/**
	* @brief set satellite prn.
	* @param[in]  prn	satellite prn
	*/
	void t_gsatinfo::set_prn(string prn)
	{
		_prn = prn;
	}

	/**
	* @brief set satellite svn.
	* @param[in]  svn	satellite svn
	*/
	void t_gsatinfo::set_svn(string svn)
	{
		_svn = svn;
	}

	/**
	* @brief set satellite launch time.
	* @param[in]  start 	satellite launch time
	*/
	void t_gsatinfo::set_start(t_gtime start)
	{
		_start = start;
	}

	/**
	* @brief set satellite decommissioned time.
	* @param[in]  end	  satellite decommissioned time
	*/
	void t_gsatinfo::set_end(t_gtime end)
	{
		_end = end;
	}

	/**
	* @brief set satellite cospar_id.
	* @param[in]  cosparid	  satellite cospar_id
	*/
	void t_gsatinfo::set_cosparid(string cosparid)
	{
		_cosparid = cosparid;
	}

	/**
	* @brief set satellite mass.
	* @param[in]  mass	  satellite mass
	*/
	void t_gsatinfo::set_mass(double mass)
	{
		_mass = mass;
	}

	/**
	* @brief set satellite max_yaw.
	* @param[in]  maxyaw	  satellite max_yaw
	*/
	void t_gsatinfo::set_maxyaw(double maxyaw)
	{
		_maxyaw = maxyaw;
	}

	/**
    * @brief set satellite LRA COM correction.
	* @param[in]  lrax,lray,lraz
	*/
	void t_gsatinfo::set_lra(double x, double y, double z)
	{
		_lra[0] = x;
		_lra[1] = y;
		_lra[2] = z;
	}

	/**
	* @brief set the type of LRA.
	* @param[in]  itype
	*/
	void t_gsatinfo::set_lra_type(int itype)
	{
		_lratype = itype;
	}

	/**
	* @brief save satellite antenna transmit power(W).
	* @param[in] lpower
	*/
	void t_gsatinfo::set_lpower(double lpower)
	{
		_lpower = lpower;
	}

	/**
	* @brief set satellite frequency id.
	* @param[in]  fid	  satellite frequency id
	*/
	void t_gsatinfo::set_fid(int fid)
	{
		_fid = fid;
	}

	/**
	* @brief set satellite block type.
	* @param[in]  blocktype	  satellite block type
	*/
	void t_gsatinfo::set_blocktype(string blocktype)
	{
		_blocktype = blocktype;
	}

	/**
	* @brief get satellite prn.
	* @return	 satellite prn.
	*/
	string t_gsatinfo::prn()
	{
		return _prn;
	}

	/**
	* @brief get satellite svn.
	* @return	 satellite svn.
	*/
	string t_gsatinfo::svn()
	{
		return _svn;
	}

	/**
	* @brief get satellite launch time.
	* @return	 satellite launch time.
	*/
	t_gtime t_gsatinfo::start()
	{
		return _start;
	}

	/**
	* @brief get satellite decommissioned time.
	* @return	 satellite decommissioned time.
	*/
	t_gtime t_gsatinfo::end()
	{
		return _end;
	}

	/**
	* @brief get satellite cospar_id.
	* @return	 satellite cospar_id.
	*/
	string t_gsatinfo::cosparid()
	{
		return _cosparid;
	}

	/**
	* @brief get satellite mass.
	* @return	 satellite mass.
	*/
	double t_gsatinfo::mass()
	{
		return _mass;
	}

	/**
	* @brief get satellite max_yaw.
	* @return	 satellite max_yaw.
	*/
	double t_gsatinfo::maxyaw()
	{
		return _maxyaw;
	}

	/**
	* @brief get satellite frequency id.
	* @return	 satellite frequency id.
	*/
	int t_gsatinfo::fid()
	{
		return _fid;
	}

	/**
	* @brief get satellite block type.
	* @return	 satellite block type.
	*/
	string t_gsatinfo::blocktype()
	{
		return _blocktype;
	}

	/**
	* @brief get satellite transmit power.
	* @return	 satellite transmit power.
	* by yqyuan
	*/
	const int t_gsatinfo::power()
	{
		return _lpower;
	}

	/**
    * @brief get satellite's LRA COM correction
	* @return	 LRA COM correction
	*/
	const double *t_gsatinfo::lra()
	{
		return _lra;
	}

	/**
	* @brief get satellite's LRA COM correction
	* @return	 LRA COM correction
	*/
	const int t_gsatinfo::lratype()
	{
		return _lratype;
	}
}