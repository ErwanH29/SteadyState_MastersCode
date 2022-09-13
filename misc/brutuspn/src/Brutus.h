#include "Star.h"
#include "Cluster.h"
#include "Bulirsch_Stoer.h"

#ifndef __Brutus_h
#define __Brutus_h

class Brutus {
  mpreal t;
  int N;  
  vector<mpreal> data;
  array<int, 4> PNtermsB;

  mpreal tolerance;
  int numBits;

  mpreal eta, dt;
  mpreal zzetaB;

  Cluster cl;
  Bulirsch_Stoer bs;

  public:
  static bool lightspeed_flagB;    
  
  Brutus();
  Brutus(vector<mpreal> &data);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits);

  // SETTERS
  void set_data(vector<mpreal> &data);
  void set_eta(mpreal &eta);
  void set_tolerance(mpreal &tolerance);
  void set_numBits(int &numBits);
  void set_t_begin(mpreal &t_begin);
  void set_lightspeed(mpreal &zzeta);
  void set_PN_terms(array<int, 4> &PNterms);
  
  // GETTERS
  mpreal get_eta(mpreal tolerance);
  mpreal get_tolerance();
  mpreal get_lightspeed();
  int get_numBits();
  int get_numBits(mpreal tolerance);
  mpreal get_t();
  vector<mpreal> get_data();
  vector<double> get_data_double();
  vector<string> get_data_string();
  mpreal get_energy();
  bool get_lightspeed_flag();

  void setup();

  void evolve(mpreal t_end);
};

#endif


