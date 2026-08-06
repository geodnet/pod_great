
#ifndef GALLPAR_H
#define GALLPAR_H 

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: parametres container
  Version: $ Rev: $

  2011-04-18 /PV: created

-*/

#include <string>
#include <set>

#include "gdata/gsatdata.h"
#include "gutils/gtriple.h"
#include "gmodels/gpar.h"
#include <unordered_map>
#include <algorithm>

using namespace std;

namespace gnut
{

	class BASE_LIBRARY_EXPORT  t_gallpar
	{
	public:
		void addParam(t_gpar);
		void addParamHead(t_gpar newPar);
		void delParam(int i);
		void delAllParam();

		int getParam(string site, t_gpar::t_type, string prn,
			t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;
		int getParIndex(int idx);
		double getParValue(int idx);
		void setParValue(int idx, double value);

		unsigned int parNumber() const;

		unsigned int orbParNumber() const;

		int getAmbParam(string site, string prn, t_gpar::t_type  type, 
			const t_gtime& beg = FIRST_TIME, const t_gtime& end = LAST_TIME) const;
		t_gpar& getAmbParam(int idx);

		int getCrdParam(string site, t_gtriple& crd,
			t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;


		vector<int> getPartialIndex(string site, string sat);


		t_gpar&   getPar(int idx);
		t_gpar&   operator[](const size_t idx);

		

		set<string> amb_prns(t_gpar::t_type type);
		set<string> amb_prns();
		set<string> amb_recs(t_gpar::t_type type);
		set<string> amb_recs();

		

		friend ostream& operator<<(ostream& os, t_gallpar& x);

	private:
		vector<t_gpar> _vParam;
		map<string, vector<t_gpar>> _vOrbParam;  // add for orb param

		// Fast get param index
		unordered_map<t_gparhead, map<t_gtimearc, long >, t_gparhead_hash > _index_par;  // index -> point
		long _max_point=0;
		long _min_point=-1;
		vector<long> _point_par;  // point -> vParam

		// Fast get parital  
		unordered_map<pair<string, string>, vector<int>, t_gpair_string_hash > _index_for_parital; // index->param
		pair<long, int> _last_point;
		void _update_partial_index();

		t_gmutex _allpar_mtx;
	};

} // namespace

#endif
