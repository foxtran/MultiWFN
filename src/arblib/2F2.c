#include "acb_hypgeom.h"

#include "2F2.h"

double hyp2F2_prec(double a1, double a2, double b1, double b2, double z,
                   slong prec) {
  const slong p = 2, q = 2;
  const int regularized = 0;
  double hyp2F2;
  arb_t res_r;
  acb_t res_ri, _z;
  acb_ptr a, b;
  arf_t resf;
  a = _acb_vec_init(p);
  b = _acb_vec_init(q);
  arb_init(res_r);
  acb_init(res_ri);
  acb_init(_z);
  arf_init(resf);
  acb_set_d(a + 0, a1);
  acb_set_d(a + 1, a2);
  acb_set_d(b + 0, b1);
  acb_set_d(b + 1, b2);
  acb_set_d(_z, z);
  acb_hypgeom_pfq(res_ri, a, p, b, q, _z, regularized, prec);
  if (acb_is_real(res_ri)) {
    acb_get_real(res_r, res_ri);
  } else {
    return 0.0;
  }
  arb_get_abs_ubound_arf(resf, res_r, prec);
  hyp2F2 = arf_get_d(resf, ARF_RND_DOWN);
}

// calculate {}_2F_2(a1,a2; b1,b2; z) using Arbitrary precision library (https://arblib.org)
// For our purpose, result should be in (0, 1] range
double hyp2F2(double a1, double a2, double b1, double b2, double z) {
  slong prec = 1000;
  double hyp2F2 = 0.0;
  // for our purpose, 0 < hyp2F2 <= 1
  while (hyp2F2 == 0.0 || hyp2F2 > 1.0) {
    hyp2F2 = hyp2F2_prec(a1, a2, b1, b2, z, prec);
    prec *= 10;
  }
  return hyp2F2;
}
