#include "Acceleration.h"

// Setting some handy parameters for the PN terms. Used also in Cluster::energies()
void Acceleration::set_PN_prod(array<mpreal, 3> &vj, array<mpreal, 3> &vi, array<mpreal, 3> &dv, array<mpreal, 3> &dr) {
    vj2     = "0";
    viInvj  = "0";
    dvij2   = "0";
    vijvj   = "0";
    
    vnj     = "0";
    vni     = "0";
    vijn    = "0";
    
    for(int k=0; k<3; k++) {
        vj2     += vj[k]*vj[k];
        viInvj  += vi[k]*vj[k];
        dvij2   += dv[k]*dv[k];
        vijvj   += dv[k]*vj[k];
    
        vnj     += vj[k]*dr[k]/Rij;
        vni     += vi[k]*dr[k]/Rij;
        vijn    += dv[k]*dr[k]/Rij;
    }
}

// 1PN pair
void Acceleration::Acceleration_PN1_pair(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv) {
    mpreal PairInProduct, Pairn;

    PairInProduct = c2*mj/(Rij*Rij)*("4"*vni - "3"*vnj);
        
    Pairn = c2*mj*apreij*(("4"/Rij*mj + "5"/Rij*mi) - vi2 + 
            "4"*viInvj - "2"*vj2 + "3"*vnj*vnj/"2");
    
    for(int k=0; k<3; k++) {
        da[k] += PairInProduct*dv[k] + Pairn*dr[k];
    }
}

// 1PN cross
void Acceleration::Acceleration_PN1_cross(const Star &si, const Star &sj, array<mpreal, 3> &da, array<mpreal, 3> &dr) {
    array<mpreal, 3> jk;
    array<mpreal, 3> ik;
    array<mpreal, 3> Cross2;
        
    mpreal aprejk, jk2, ik2, Rjk, Rik;
    mpreal mk;
        
    mpreal Cross1i;
    mpreal InProductxx;
    
    for (vector<Star>::iterator sk = s.begin(); sk != s.end(); ++sk) {
        if(*sk == si || *sk == sj) continue;
        
        mk = sk->m;
        // Set InProductxx, jk2 en ik2 to 0 at the beginning of each iteration in sk
        
        jk2         = "0";
        ik2         = "0";
        InProductxx = "0";
        
        for(int k=0; k<3; k++) {
            jk[k] = sj.r[k] - sk->r[k];
            ik[k] = si.r[k] - sk->r[k];
            jk2 += jk[k]*jk[k];
            ik2 += ik[k]*ik[k];
            InProductxx += dr[k]*jk[k];
        }
        
        Rjk = sqrt(jk2 + eps2);
        Rik = sqrt(ik2 + eps2);
        
        aprejk = "1"/(Rjk*Rjk*Rjk);

        Cross1i += mk*(("1"/Rjk + "4"/Rik) - "1"*aprejk*InProductxx/"2"); 
        
        for(int k=0; k<3; k++) {
            Cross2[k] += mk*aprejk*jk[k];
        }
    }
    
    for(int k=0; k<3; k++) {
        da[k] += c2*mj*(apreij*dr[k]*Cross1i - "3.5"/Rij*Cross2[k]);
    }
}

// 2PN
void Acceleration::Acceleration_PN2(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv) {
    mpreal factor2n, factor2v;
    
    factor2n = c2*c2*apreij*mj*("-2"*vj2*vj2 + ("4"*vj2 - "2"*viInvj)*viInvj + ("1.5"*vi2 + "4.5"*vj2 - "6"*viInvj - 
            "1.875"*vnj*vnj)*vnj*vnj - ("14.25"*mi*mi + "9"*mj*mj + "34.5"*mi*mj)/(Rij*Rij) + 
        mi/Rij*("-3.75"*vi2 + "1.25"*vj2 - "2.5"*viInvj + "19.5"*vni*vni - "39"*vni*vnj + "8.5"*vnj*vnj) + 
        mj/Rij*("4"*vj2 - "8"*viInvj + "2"*vni*vni - "4"*vni*vnj - "6"*vnj*vnj));
      
    factor2v = c2*c2*apreij*Rij*mj*(mi/Rij*("-15.75"*vni + "13.75"*vnj) - "2"*mj/Rij*(vni + vnj) + vi2*vnj + "4"*vj2*vni -
            "5"*vj2*vnj - "4"*viInvj*vijn - "6"*vni*vnj*vnj + "4.5"*vnj*vnj*vnj);
    
    for(int k=0; k<3; k++) {
        da[k] += factor2n*dr[k] + factor2v*dv[k];
    }
}

// 2.5PN
void Acceleration::Acceleration_PN2_5(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv) {
    mpreal factor2_5n, factor2_5v;
    
    factor2_5n = "0.8"*c5*mi*mj*apreij*vijn/Rij*("-6"*mi/Rij + "52"*mj/(Rij*"3") + "3"*dvij2);
    
    factor2_5v = "0.8"*c5*mi*mj*apreij*("2"*mi/Rij - "8"*mj/Rij - dvij2);
    
    for(int k=0; k<3; k++) {
        da[k] += factor2_5n*dr[k] + factor2_5v*dv[k];
    }
}

// Here the acceleration used during Bulirsch-Stoer steps
void Acceleration::calcAcceleration() {
    array<mpreal, 3> vi;
    array<mpreal, 3> dr;
    array<mpreal, 3> dv;
    array<mpreal, 3> da;
    mpreal dr2 = "0";

    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        si->a.fill("0");
    }

    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        array<mpreal, 3> vj;
        vi2 = "0";
        
        for(int k=0; k<3; k++) {
            vi[k] = si->v[k];
            vi2  += vi[k]*vi[k];
        }
        mi = si->m;
        
        for (vector<Star>::iterator sj = s.begin(); sj != s.end(); ++sj) {
            if(*sj == *si) continue;
            
            dr2 = "0";
            da.fill("0");
            
            for(int k=0; k<3; k++) {
                vj[k] = sj->v[k];
                dr[k] = si->r[k]-sj->r[k];
                dv[k] = si->v[k]-sj->v[k];
                dr2  += dr[k]*dr[k];
            }
            mj = sj->m;

            Rij = sqrt(dr2 + eps2);
            apreij = "1"/(Rij*Rij*Rij);

            set_PN_prod(vj, vi, dv, dr);

            if (PN1p)   Acceleration_PN1_pair(da, dr, dv);
            if (PN1c)   Acceleration_PN1_cross(*si, *sj, da, dr);
            if (PN2)    Acceleration_PN2(da, dr, dv);
            if (PN2_5)  Acceleration_PN2_5(da, dr, dv);

            for(int k=0; k<3; k++) {
                da[k] -= mj*apreij*dr[k];
                si->a[k] += da[k];
            }
        }
    }
}

// Here, the acceleration calculated at the very beginning of a new time step.
// The initial stepsize is calculated here
void Acceleration::calcAcceleration_dt() {
    array<mpreal, 3> vi;
    array<mpreal, 3> dr;
    array<mpreal, 3> dv;
    array<mpreal, 3> da;
    mpreal da2 = "0";
    mpreal dr2 = "0";

    dt = "1e100";
    mpreal mydt = "0";
    
    int N = s.size();
    
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        si->a.fill("0");
    }

    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        array<mpreal, 3> vj;
        vi2 = "0";
        
        for(int k=0; k<3; k++) {
            vi[k] = si->v[k];
            vi2  += vi[k]*vi[k];
        }
        mi = si->m;
        
        for (vector<Star>::iterator sj = s.begin(); sj != s.end(); ++sj) {
            if(*sj == *si) continue;
            
            dr2 = "0";
            da.fill("0");
            for(int k=0; k<3; k++) {
                vj[k] = sj->v[k];
                dr[k] = si->r[k]-sj->r[k];
                dv[k] = si->v[k]-sj->v[k];
                dr2  += dr[k]*dr[k];
            }
            mj = sj->m;

            Rij = sqrt(dr2 + eps2);
            apreij = "1"/(Rij*Rij*Rij);

            set_PN_prod(vj, vi, dv, dr);

            if (PN1p)   Acceleration_PN1_pair(da, dr, dv);
            if (PN1c)   Acceleration_PN1_cross(*si, *sj, da, dr);
            if (PN2)    Acceleration_PN2(da, dr, dv);
            if (PN2_5)  Acceleration_PN2_5(da, dr, dv);
            
            da2 = "0";
            for(int k=0; k<3; k++) {
                da[k]    -= mj*apreij*dr[k];
                si->a[k] += da[k];
                da2      += da[k]*da[k];
            }
            mydt = dr2 / da2;
            if(mydt < dt) dt = mydt;
        }
    }
    dt = sqrt(sqrt(dt)); //pow(dt, "0.25");
}
