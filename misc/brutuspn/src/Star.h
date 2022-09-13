#include <iostream>
using namespace std;

#include <array>

#include "mpreal.h"
using namespace mpfr;

#ifndef __Star_h
#define __Star_h

class Star {
public:
    mpreal m;
    mpreal ener;
    array<mpreal, 3> r;
    array<mpreal, 3> v;
    array<mpreal, 3> a;
    array<mpreal, 3> w;
    
    array<mpreal, 3> vp;

    Star();
    Star(mpreal m, array<mpreal, 3> r, array<mpreal, 3> v);
  
    friend bool operator== (const Star &si, const Star &sj) {
        if(si.r[0] == sj.r[0] && si.r[1] == sj.r[1] && si.r[2] == sj.r[2])
            return true;
        else
            return false;
    }
    
    friend ostream & operator << (ostream &so, const Star &si) {
        so << si.m << " " << si.r[0] << " "<< si.r[1] << " "<< si.r[2] << " "
                      << si.v[0] << " "<< si.v[1] << " "<< si.v[2] << endl;
        return so; 
  }
};

#endif


