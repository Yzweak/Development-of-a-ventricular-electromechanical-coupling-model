//
// Created by dell on 2023/5/6.
//

#ifndef MIXVENTRICLE_MODEL_CUH
#define MIXVENTRICLE_MODEL_CUH


#define HT 0.02
#define stimduration 1.0
#define stimstrength (-80)   //-52
#define dx 0.168
#define cellNum 2968
#define stimNum 42
#define endoNum 377
#define myoNum 1335
#define epiNum 1256
//MI_LONG
//#define ischemia_1a_num 0
//#define ischemia_1b_num 0
//#define MI_short_num 0
//#define MI_long_num 322
//圆环
//#define ischemia_1a_num 105
//#define ischemia_1b_num 77
//#define MI_short_num 64
//#define MI_long_num 76
//环形
#define ischemia_1a_num 133
#define ischemia_1b_num 153
#define MI_short_num 159
#define MI_long_num 185

#define rows 2968
#define cols 4
#define XD 64
#define YD 64
#define Dlong 0.007  // 0.167  0.084  0.028  0.008

#define Cao 2.0
#define Nao 140.0
//Intracellular volumes
#define Vc 0.016404
#define Vsr 0.001094
#define Vss 0.00005468
//Calcium buffering dynamics
#define Bufc 0.2
#define Kbufc 0.001
#define Bufsr 10.
#define Kbufsr 0.3
#define Bufss 0.4
#define Kbufss 0.00025
//Intracellular calcium flux dynamics
#define Vmaxup 0.006375
#define Kup 0.00025
#define Vrel 0.102
#define k1_ 0.15
#define k2_ 0.045
#define k3 0.060
#define k4 0.005
#define EC 1.5
#define maxsr 2.5
#define minsr 1.
#define Vleak 0.00036
#define Vxfer 0.0038
//Constants
#define R 8314.472
#define F 96485.3415
#define T 310.0
#define RTONF ((R * T) / F)
//Cellular capacitance
#define CAPACITANCE 0.185
//*****************Parameters for currents
//Parameters for IKr
#define Gkr 0.153
//Parameters for Iks
#define pKNa 0.03
//Parameters for Ik1
#define GK1 5.405
//Parameters for INa
#define GNa 14.838
#define GNaL 0.0075
//Parameters for IbNa
#define GbNa 0.00029
//Parameters for INaK
#define KmK 1.0
#define KmNa 40.0
#define knak 2.724
//Parameters for ICaL
#define GCaL 0.00003980
//Parameters for IbCa
#define GbCa 0.000592
//Parameters for INaCa
#define knaca 1000
#define KmNai 87.5
#define KmCa 1.38
#define ksat 0.1
#define n 0.35
//Parameters for IpCa
#define GpCa 0.1238
#define KpCa 0.0005
//Parameters for IpK
#define GpK 0.0146


enum cell_type
{
    PM,
    SAN,
    Nothing
};

enum cell_location
{
    EPI,
    ENDO,
    MCELL
};


#endif //MIXVENTRICLE_MODEL_CUH
