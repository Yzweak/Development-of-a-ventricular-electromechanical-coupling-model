/* This function performs the actual computations of currents, concentration changes
   voltage changes, gate state changes involved in a single time step. Every 100 timesteps
   the values of all currents are written to file. */

#include "Header.h"

extern double knak;
extern double KmNa;
extern double KmK;
extern double knaca;
extern double KmNai;
extern double KmCa;
extern double ksat;
extern double n;


extern double Ko;
extern double Cao;
extern double Nao;


extern double Bufc;
extern double Kbufc;
extern double Bufsr;
extern double Kbufsr;
extern double Bufss;
extern double Kbufss;

extern double Vmaxup;
extern double Kup;
extern double Vrel;
extern double k1_;
extern double k2_;
extern double k3;
extern double k4;
extern double EC;
extern double maxsr;
extern double minsr;
extern double Vleak;
extern double Vxfer;


extern double pKNa;


extern double CAPACITANCE;
extern double RTONF;
extern double F;
extern double R;
extern double T;


extern double Gkr;
extern double Gks;
extern double GK1;
extern double Gto;
extern double GNa;
extern double GbNa;
extern double GCaL;
extern double GbCa;
extern double GpCa;
extern double KpCa;
extern double GpK;


extern double Vc;
extern double Vsr;
extern double Vss;

extern double period;



void Step(Variables* V, double HT, double* tt, int step, double Istim)
{

    double IKr;
    double IKs;
    double IK1;
    double Ito;
    double INa;
    double IbNa;
    double ICaL;
    double INaL;  //Zhai change
    double IbCa;
    double INaCa;
    double IpCa;
    double IpK;
    double INaK;
    double Irel;
    double Ileak;
    double Iup;
    double Ixfer;
    double k1;
    double k2;
    double kCaSR;


    double dNai;
    double dKi;
    double dCai;
    double dCaSR;
    double dCaSS;
    double dRR;


    double Ek;
    double Ena;
    double Eks;
    double Eca;
    double CaCSQN;
    double bjsr;
    double cjsr;
    double CaSSBuf;
    double bcss;
    double ccss;
    double CaBuf;
    double bc;
    double cc;
    double Ak1;
    double Bk1;
    double rec_iK1;
    double rec_ipK;
    double rec_iNaK;
    double AM;
    double BM;
    double AH_1;
    double BH_1;
    double AH_2;
    double BH_2;
    double AJ_1;
    double BJ_1;
    double AJ_2;
    double BJ_2;
    double M_INF;
    double H_INF;
    double J_INF;
    double TAU_M;
    double TAU_H;
    double TAU_J;
    double axr1;
    double bxr1;
    double axr2;
    double bxr2;
    double Xr1_INF;
    double Xr2_INF;
    double TAU_Xr1;
    double TAU_Xr2;
    double Axs;
    double Bxs;
    double Xs_INF;
    double TAU_Xs;
    double R_INF;
    double TAU_R;
    double S_INF;
    double TAU_S;
    double Ad;
    double Bd;
    double Cd;
    double Af;
    double Bf;
    double Cf;
    double Af2;
    double Bf2;
    double Cf2;
    double TAU_D;
    double D_INF;
    double TAU_F;
    double F_INF;
    double TAU_F2;
    double F2_INF;
    double TAU_FCaSS;
    double FCaSS_INF;


    static double inverseVcF2 = 1 / (2 * Vc * F);
    static double inverseVcF = 1. / (Vc * F);
    static double inversevssF2 = 1 / (2 * Vss * F);

    

    // define all variables
#define       sm          (*V).M
#define       sh          (*V).H
#define       sj          (*V).J
#define       sxr1        (*V).Xr1
#define       sxr2        (*V).Xr2
#define       sxs         (*V).Xs
#define       ss          (*V).S
#define       sr          (*V).R
#define       sd          (*V).D
#define       sf          (*V).F
#define       sf2         (*V).F2
#define       sfcass      (*V).FCass
#define       sRR         (*V).RR
#define       sOO         (*V).OO
#define       svolt       (*V).Volt
#define       Cai         (*V).Cai
#define       CaSR        (*V).CaSR
#define       CaSS        (*V).CaSS
#define       Nai         (*V).Nai
#define       Ki          (*V).Ki
#define       sItot       (*V).Itot
#define       force       (*V).Force
#define       P0          (*V).P0
#define       P1          (*V).P1
#define       P2          (*V).P2
#define       P3          (*V).P3
#define       N0          (*V).N0
#define       N1          (*V).N1
#define       LTRPNCa     (*V).LTRPNCa
  
   

    //Needed to compute currents
    Ek = RTONF * (log((Ko / Ki)));
    Ena = RTONF * (log((Nao / Nai)));
    Eks = RTONF * (log((Ko + pKNa * Nao) / (Ki + pKNa * Nai)));
    Eca = 0.5 * RTONF * (log((Cao / Cai)));
    Ak1 = 0.1 / (1. + exp(0.06 * (svolt - Ek - 200)));
    Bk1 = (3. * exp(0.0002 * (svolt - Ek + 100)) +
        exp(0.1 * (svolt - Ek - 10))) / (1. + exp(-0.5 * (svolt - Ek)));
    rec_iK1 = Ak1 / (Ak1 + Bk1);
    rec_iNaK = (1. / (1. + 0.1245 * exp(-0.1 * svolt * F / (R * T)) + 0.0353 * exp(-svolt * F / (R * T))));
    rec_ipK = 1. / (1. + exp((25 - svolt) / 5.98));


    //Compute currents
    INa = 0.38 * GNa * sm * sm * sm * sh * sj * (svolt - Ena);//Zhai change
    ICaL = 0.62 * GCaL * sd * sf * sf2 * sfcass * 4 * (svolt - 15) * (F * F / (R * T)) *//Zhai change
        (0.25 * exp(2 * (svolt - 15) * F / (R * T)) * CaSS - Cao) / (exp(2 * (svolt - 15) * F / (R * T)) - 1.);
    Ito = 0.37 * Gto * sr * ss * (svolt - Ek);//Zhai change
    IKr = 0.30 * Gkr * sqrt(Ko / 5.4) * sxr1 * sxr2 * (svolt - Ek);//Zhai change
    IKs = 0.20 * Gks * sxs * sxs * (svolt - Eks);//Zhai change
    IK1 = GK1 * rec_iK1 * (svolt - Ek);
    INaCa = knaca * (1. / (KmNai * KmNai * KmNai + Nao * Nao * Nao)) * (1. / (KmCa + Cao)) *
        (1. / (1 + ksat * exp((n - 1) * svolt * F / (R * T)))) *
        (exp(n * svolt * F / (R * T)) * Nai * Nai * Nai * Cao -
            exp((n - 1) * svolt * F / (R * T)) * Nao * Nao * Nao * Cai * 2.5);
    INaK = knak * (Ko / (Ko + KmK)) * (Nai / (Nai + KmNa)) * rec_iNaK;
    IpCa = GpCa * Cai / (KpCa + Cai);
    IpK = GpK * rec_ipK * (svolt - Ek);
    IbNa = GbNa * (svolt - Ena);
    IbCa = GbCa * (svolt - Eca);

    // 补充IKATP电流    Zhai change
    //根据：2011_Experiment-model interaction for analysis of epicardial activation during human ventricular fibrillation with global myocardial ischaemia
    double GKATP = 3.9;                                             
    double H = 2.0;                                                  
    double n = 0.24;                                              
    double ATPi = 4.6;             
    double Khalf = 0.25;         
    double IKATP = GKATP * (1 / (1 + pow(ATPi / Khalf, H))) * pow(Ko / 5.4, n) * (svolt - Ek);//liang_change
    //IKATP=0;//liang_change

     //将ORD模型里的INaL加入    //Zhai change
    double mL = 0;
    double hL = 1;
    double mLss = 1.0 / (1.0 + exp((-(svolt + 42.85)) / 5.264));
    double tmL = 1.0 / (6.765 * exp((svolt + 11.64) / 34.77) + 8.552 * exp(-(svolt + 77.42) / 5.955));
    mL = mLss - (mLss - mL) * exp(-HT / tmL);
    double hLss = 1.0 / (1.0 + exp((svolt + 87.61) / 7.488));
    double thL = 200.0;
    hL = hLss - (hLss - hL) * exp(-HT / thL);
    double GNaL = 0.0075;
    INaL = GNaL * (svolt - Ena) * mL * hL;

    //Determine total current
    (sItot) = IKr +
        IKs +
        IK1 +
        Ito +
        INa +
        IbNa +
        ICaL +
        IbCa +
        INaK +
        INaCa +
        IpCa +
        IpK +
        IKATP +
        INaL +                                                            //Zhai change
        Istim;





    //update concentrations
    kCaSR = maxsr - ((maxsr - minsr) / (1 + (EC / CaSR) * (EC / CaSR)));
    k1 = k1_ / kCaSR;
    k2 = k2_ * kCaSR;
    dRR = k4 * (1 - sRR) - k2 * CaSS * sRR;
    sRR += HT * dRR;
    sOO = k1 * CaSS * CaSS * sRR / (k3 + k1 * CaSS * CaSS);


    Irel = Vrel * sOO * (CaSR - CaSS);
    Ileak = Vleak * (CaSR - Cai);
    Iup = Vmaxup / (1. + ((Kup * Kup) / (Cai * Cai)));
    Ixfer = Vxfer * (CaSS - Cai);


    CaCSQN = Bufsr * CaSR / (CaSR + Kbufsr);
    dCaSR = HT * (Iup - Irel - Ileak);
    bjsr = Bufsr - CaCSQN - dCaSR - CaSR + Kbufsr;
    cjsr = Kbufsr * (CaCSQN + dCaSR + CaSR);
    CaSR = (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2;


    CaSSBuf = Bufss * CaSS / (CaSS + Kbufss);
    dCaSS = HT * (-Ixfer * (Vc / Vss) + Irel * (Vsr / Vss) + (-ICaL * inversevssF2 * CAPACITANCE));
    bcss = Bufss - CaSSBuf - dCaSS - CaSS + Kbufss;
    ccss = Kbufss * (CaSSBuf + dCaSS + CaSS);
    CaSS = (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2;


    CaBuf = Bufc * Cai / (Cai + Kbufc);
    dCai = HT * ((-(IbCa + IpCa - 2 * INaCa) * inverseVcF2 * CAPACITANCE) - (Iup - Ileak) * (Vsr / Vc) + Ixfer);
    bc = Bufc - CaBuf - dCai - Cai + Kbufc;
    cc = Kbufc * (CaBuf + dCai + Cai);
    Cai = (sqrt(bc * bc + 4 * cc) - bc) / 2;


    dNai = -(INa + IbNa + 3 * INaK + 3 * INaCa) * inverseVcF * CAPACITANCE;
    Nai += HT * dNai;

    dKi = -(Istim + IK1 + Ito + IKr + IKs - 2 * INaK + IpK) * inverseVcF * CAPACITANCE;
    Ki += HT * dKi;

      //compute steady state values and time constants
    AM = 1. / (1. + exp((-60. - svolt) / 5.));
    BM = 0.1 / (1. + exp((svolt + 35.) / 5.)) + 0.10 / (1. + exp((svolt - 50.) / 200.));
    TAU_M = AM * BM;
    M_INF = 1. / ((1. + exp((-56.86 - svolt) / 9.03)) * (1. + exp((-56.86 - svolt) / 9.03)));
    if (svolt >= -40.)
    {
        AH_1 = 0.;
        BH_1 = (0.77 / (0.13 * (1. + exp(-(svolt + 10.66) / 11.1))));
        TAU_H = 1.0 / (AH_1 + BH_1);
    }
    else
    {
        AH_2 = (0.057 * exp(-(svolt + 80.) / 6.8));
        BH_2 = (2.7 * exp(0.079 * svolt) + (3.1e5) * exp(0.3485 * svolt));
        TAU_H = 1.0 / (AH_2 + BH_2);
    }

    H_INF = 1. / ((1. + exp((svolt + 74.95) / 7.43)) * (1. + exp((svolt + 74.95) / 7.43)));
    if (svolt >= -40.)
    {
        AJ_1 = 0.;
        BJ_1 = (0.6 * exp((0.057) * svolt) / (1. + exp(-0.1 * (svolt + 32.))));
        TAU_J = 1.0 / (AJ_1 + BJ_1);
    }
    else
    {
        AJ_2 = (((-2.5428e4) * exp(0.2444 * svolt) - (6.948e-6) *
            exp(-0.04391 * svolt)) * (svolt + 37.78) /
            (1. + exp(0.311 * (svolt + 79.23))));
        BJ_2 = (0.02424 * exp(-0.01052 * svolt) / (1. + exp(-0.1378 * (svolt + 40.14))));
        TAU_J = 1.0 / (AJ_2 + BJ_2);
    }
    J_INF = H_INF;

    Xr1_INF = 1. / (1. + exp((-26. - svolt) / 7.));
    axr1 = 450. / (1. + exp((-45. - svolt) / 10.));
    bxr1 = 6. / (1. + exp((svolt - (-30.)) / 11.5));
    TAU_Xr1 = axr1 * bxr1;
    Xr2_INF = 1. / (1. + exp((svolt - (-88.)) / 24.));
    axr2 = 3. / (1. + exp((-60. - svolt) / 20.));
    bxr2 = 1.12 / (1. + exp((svolt - 60.) / 20.));
    TAU_Xr2 = axr2 * bxr2;

    Xs_INF = 1. / (1. + exp((-5. - svolt) / 14.));
    Axs = (1400. / (sqrt(1. + exp((5. - svolt) / 6))));
    Bxs = (1. / (1. + exp((svolt - 35.) / 15.)));
    TAU_Xs = Axs * Bxs + 80;

#ifdef EPI
    R_INF = 1. / (1. + exp((20 - svolt) / 6.));
    S_INF = 1. / (1. + exp((svolt + 20) / 5.));
    TAU_R = 9.5 * exp(-(svolt + 40.) * (svolt + 40.) / 1800.) + 0.8;
    TAU_S = 85. * exp(-(svolt + 45.) * (svolt + 45.) / 320.) + 5. / (1. + exp((svolt - 20.) / 5.)) + 3.;
#endif
#ifdef ENDO
    R_INF = 1. / (1. + exp((20 - svolt) / 6.));
    S_INF = 1. / (1. + exp((svolt + 28) / 5.));
    TAU_R = 9.5 * exp(-(svolt + 40.) * (svolt + 40.) / 1800.) + 0.8;
    TAU_S = 1000. * exp(-(svolt + 67) * (svolt + 67) / 1000.) + 8.;
#endif
#ifdef MCELL
    R_INF = 1. / (1. + exp((20 - svolt) / 6.));
    S_INF = 1. / (1. + exp((svolt + 20) / 5.));
    TAU_R = 9.5 * exp(-(svolt + 40.) * (svolt + 40.) / 1800.) + 0.8;
    TAU_S = 85. * exp(-(svolt + 45.) * (svolt + 45.) / 320.) + 5. / (1. + exp((svolt - 20.) / 5.)) + 3.;
#endif


    D_INF = 1. / (1. + exp((-8 - svolt) / 7.5));
    Ad = 1.4 / (1. + exp((-35 - svolt) / 13)) + 0.25;
    Bd = 1.4 / (1. + exp((svolt + 5) / 5));
    Cd = 1. / (1. + exp((50 - svolt) / 20));
    TAU_D = Ad * Bd + Cd;
    F_INF = 1. / (1. + exp((svolt + 20) / 7));
    Af = 1102.5 * exp(-(svolt + 27) * (svolt + 27) / 225);
    Bf = 200. / (1 + exp((13 - svolt) / 10.));
    Cf = (180. / (1 + exp((svolt + 30) / 10))) + 20;
    TAU_F = Af + Bf + Cf;
    F2_INF = 0.67 / (1. + exp((svolt + 35) / 7)) + 0.33;
    Af2 = 600 * exp(-(svolt + 25) * (svolt + 25) / 170);
    Bf2 = 31 / (1. + exp((25 - svolt) / 10));
    Cf2 = 16 / (1. + exp((svolt + 30) / 10));
    TAU_F2 = Af2 + Bf2 + Cf2;
    FCaSS_INF = 0.6 / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 0.4;
    TAU_FCaSS = 80. / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 2.;



    //Update gates
    sm = M_INF - (M_INF - sm) * exp(-HT / TAU_M);
    sh = H_INF - (H_INF - sh) * exp(-HT / TAU_H);
    sj = J_INF - (J_INF - sj) * exp(-HT / TAU_J);
    sxr1 = Xr1_INF - (Xr1_INF - sxr1) * exp(-HT / TAU_Xr1);
    sxr2 = Xr2_INF - (Xr2_INF - sxr2) * exp(-HT / TAU_Xr2);
    sxs = Xs_INF - (Xs_INF - sxs) * exp(-HT / TAU_Xs);
    ss = S_INF - (S_INF - ss) * exp(-HT / TAU_S);
    sr = R_INF - (R_INF - sr) * exp(-HT / TAU_R);
    sd = D_INF - (D_INF - sd) * exp(-HT / TAU_D);
    sf = F_INF - (F_INF - sf) * exp(-HT / TAU_F);
    sf2 = F2_INF - (F2_INF - sf2) * exp(-HT / TAU_F2);
    sfcass = FCaSS_INF - (FCaSS_INF - sfcass) * exp(-HT / TAU_FCaSS);



    //update voltage
    svolt = svolt + HT * (-sItot);

    //calculate Force 
    double sumPath = g01 * g12 * g23 + f01 * g12 * g23 + f01 * f12 * g23 + f01 * f12 * f23;
    double Pmax1 = (f01 * g12 * g23) / sumPath;
    double Pmax2 = (f01 * f12 * g23) / sumPath;
    double Pmax3 = (f01 * f12 * f23) / sumPath;
    force = 0.1 * alpha_SL * (P1 + N1 + 2 * P2 + 3 * P3) / (Pmax1 + 2 * Pmax2 + 3 * Pmax3);


    double knp = kpn * pow(LTRPNCa * inv_LTRPNtot_Ktrop_half, ntrop);
    double FN_Ca = alpha_SL * (N1 + P1 + 2 * P2 + 3 * P3) / (Pmax1 + Pmax2 + Pmax3);
    double mod_factor = 1 + (2.3 - SL) / pow(0.6, 1.6);

    double dN1 = kpn * P1 - (knp + 0.03 * mod_factor) * N1;
    double dP0 = -1 * (kpn + f01) * P0 + knp * N0 + g01 * P1;
    double dP1 = -1 * (kpn + f12 + g01) * P1 + knp * N1 + f01 * P0 + g12 * P2;
    double dP2 = -1 * (f23 + g12) * P2 + f12 * P1 + g23 * P3;
    double dP3 = -1 * g23 * P3 + f23 * P2;

    N1 += dN1 * HT;
    P0 += dP0 * HT;
    P1 += dP1 * HT;
    P2 += dP2 * HT;
    P3 += dP3 * HT;
    N0 += (-1 * (dP0 + dP1 + dP2 + dP3 + dN1)) * HT;
    LTRPNCa += (100 * Cai * (0.04 - LTRPNCa) - 40e-3 * LTRPNCa * (1.0 - (2.0 / 3.0) * FN_Ca)) * HT;


}
