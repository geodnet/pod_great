
#ifndef GTETRAD_H
#define GTETRAD_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: implements 4D representation (e.g. coordinates + time)
  Version: $ Rev: $

  2012-05-11 /JD: created

-*/

#include <iostream>
#include <string.h>

#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "ExportHeaders/ExportBase.h"

using namespace std;

namespace gnut {

	class t_gtetrad {

	public:
		t_gtetrad();
		t_gtetrad(double x, double y, double z, double t);
		t_gtetrad(double crd[4]);
		t_gtetrad(const ColumnVector& crd);
		virtual ~t_gtetrad();

		t_gtetrad& operator=(const t_gtetrad& other);     // assignment operator
		t_gtetrad  operator+(const t_gtetrad& other) const;     //
		bool       operator==(const t_gtetrad& tr) const; // equal operator
		bool       operator<(const t_gtetrad& tr) const;  // equal for sorting
		double& operator[](const size_t idx);          // get a reference of element
		double     operator[](const size_t idx) const;    // get value of element
		friend ostream& operator<<(ostream& os, const t_gtetrad& x);

		double        crd(int idx) const;                 // get single element


	protected:

	private:
		double         _crd[4];

	};

} // namespace

#endif