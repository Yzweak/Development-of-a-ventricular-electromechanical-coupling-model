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
    //Zhai change
    M = 0.033655635139547148000;   
    H = 0.211053852569920620000;  
    J = 0.180147754598375870000;   
    Xr1 = 0.024098542643841897000; 
    Xr2 = 0.323941263226489200000;
    Xs = 0.009339544433470185100; 
    R = 0.000000288926510767676;   
    S = 0.999957603371733010000;  
    D = 0.000164510100997893920;   
    F = 0.980144782112325140000;   
    F2 = 0.995728345214580090000;   
    FCass = 0.999969154239247040000;   
    RR = 0.989633139742973960000;  
    OO = 0.000000171795851160072;  
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





