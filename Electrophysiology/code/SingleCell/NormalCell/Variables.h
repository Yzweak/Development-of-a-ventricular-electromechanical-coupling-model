/* This file contains the object description for the Variables object
   that is intended both to be able to send all state variables as a
   package to functions needing it and allowing efficient handling of
   the set of state variables. */

class Variables
{
public:
	//voltage + backup
	double   Volt;
	double   Cai;
	double   CaSR;
	double   CaSS;
	double   Nai;
	double   Ki;


	//states of voltage and time dependent gates
	//INa
	double   M;
	double   H;
	double   J;
	//IKr
	//IKr1
	double   Xr1;
	//IKr2
	double   Xr2;
	//IKs
	double   Xs;
	//Ito1
	double   R;
	double   S;
	//ICa
	double   D;
	double   F;
	double   F2;
	double   FCass;
	//Irel
	double   RR;
	double   OO;

	//total current
	double   Itot;

	double P0;
	double P1;
	double P2;
	double P3;
	double N0;
	double N1;

	double LTRPNCa;
	double HTRPNCa;

	double Force;

public:
	Variables(double V_init, double Cai_init, double CaSR_init, double CaSS_init, double Nai_init, double Ki_init);
	void writebackup(double* t);
};

#pragma once
