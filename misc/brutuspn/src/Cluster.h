#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <numeric> 
#include <cstdlib>

#include "Star.h"
#include "Acceleration.h"

#ifndef __Cluster_h
#define __Cluster_h

class Cluster : public Acceleration { 
public:
    mpreal time, dt_last;
    Cluster() : Acceleration() {}
    
    Cluster(vector<double> data);
    Cluster(vector<mpreal> data);

    vector<double> get_data_double();
    vector<mpreal> get_data();  
    
    void setw0();

    void updatePositions(mpreal dt);
    void updateAuxiliary(mpreal dt);
    bool updateVelocities(mpreal dt);
    void updateAuxiliary2(mpreal dt);

    bool step(mpreal &dt);
  
    mpreal e1(const Star &si, int N);
    mpreal e2(int N);
    mpreal energies();
    
    friend ostream & operator << (ostream &so, Cluster &cl) {
        for (vector<Star>::iterator si = cl.s.begin(); si != cl.s.end(); ++si) {
            so << *si;
        }
        return so;
    }
};

#endif


