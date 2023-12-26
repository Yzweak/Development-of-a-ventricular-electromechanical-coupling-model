/* This file defines the functions that can be performed with the Variable object defined in
   Variables.h*/

#include "Header.h"

using namespace std;

Variables::Variables(double V_init, double Cai_init, double CaSR_init, double CaSS_init, double Nai_init, double Ki_init)
{
    Volt = V_init;
    Cai = Cai_init;
    CaSR = CaSR_init;
    CaSS = CaSS_init;
    Nai = Nai_init;
    Ki = Ki_init;
    M = 0.037396266233898176000;  
    H = 0.193531112266500660000;  
    J = 0.170193474094604810000;   
    Xr1 = 0.013720983480392804000;
    Xr2 = 0.318614799111897550000; 
    Xs = 0.009711833626984207200;  
    R = 0.000000318575104674730;   
    S = 0.999952334011770310000;   
    D = 0.000177883400051902350;   
    F = 0.988991405078883500000;   
    F2 = 0.995358028564252040000;   
    FCass = 0.999978028985132990000;  
    RR = 0.992241352890124500000;  
    OO = 0.000000129042183903380; 
    Itot = 0;
    //Force states
    LTRPNCa = 4.3383187618764628e-003;
    HTRPNCa = 1.1347937598465910e-001;
    N0 = 9.9997226457259825e-001;
    N1 = 3.5681318281388550e-006;
    P0 = 4.1098247014058200e-006;
    P1 = 3.5644945690450202e-006;
    P2 = 6.6623946400458660e-006;
    P3 = 5.8053275016232687e-006;

    Force = 0;
    printf("Variables initialized\n");
}



void Variables::writebackup(double* t)
{
    static char filename[300];

    sprintf(filename, "%s", "PointBackupData.dat");

    ofstream oo(filename, ios::app);
    if (!oo)
    {
        printf("cannot open file %s\n", filename);
        exit(1);
    }

    oo << floor(*t - 99000) << "\t";  
    oo << Volt << "\t";
    oo << Cai << "\t";
    oo << CaSR << "\t";
    oo << CaSS << "\t";
    oo << Nai << "\t";
    oo << Ki << "\t";
    oo << M << "\t";
    oo << H << "\t";
    oo << J << "\t";
    oo << Xr1 << "\t";
    oo << Xr2 << "\t";
    oo << Xs << "\t";
    oo << S << "\t";
    oo << R << "\t";
    oo << D << "\t";
    oo << F << "\t";
    oo << F2 << "\t";
    oo << FCass << "\t";
    oo << RR << "\t";
    oo << OO << "\t";
    oo << Force << "\t";
    oo << Itot;    
    oo << endl;
    oo.close();

}





