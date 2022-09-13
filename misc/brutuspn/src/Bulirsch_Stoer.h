#include "Star.h"
#include "Cluster.h"
#include <tuple>

#ifndef __Bulirsch_Stoer_h
#define __Bulirsch_Stoer_h

class Bulirsch_Stoer {
  mpreal tolerance;
  int n_max, k_max;

  public:

  Bulirsch_Stoer();
  Bulirsch_Stoer(mpreal tolerance);
  Bulirsch_Stoer(mpreal tolerance, int n_max, int k_max);

  void set_tolerance(mpreal tolerance);
  void set_n_max(int n_max);
  void set_k_max(int k_max);

  mpreal get_tolerance();
  int get_n_max();
  int get_k_max();

  std::tuple<bool, bool> integrate(Cluster &cl, mpreal &dt);
  std::tuple<bool, bool> step(Cluster &cl, mpreal &dt);
  void extrapol(Cluster &cl_exp, vector<mpreal> &dt, vector<Cluster> &c);
  mpreal extrapolate(vector<mpreal> x, vector<mpreal> y);
  bool error_control(Cluster &c1, Cluster &c2);
};

#endif


