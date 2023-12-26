#include <iostream>

using namespace std;

#include "Header.h"


#define DYNRESTPROTOCOL
//#define  S1S2RESTPROTOCOL

//External concentrations
double Ko = 5.4;              //Zhai change
double Cao = 2.0;
double Nao = 140.0;

//Intracellular volumes
double Vc = 0.016404;
double Vsr = 0.001094;
double Vss = 0.00005468;

//Calcium buffering dynamics
double Bufc = 0.2;
double Kbufc = 0.001;
double Bufsr = 10.;
double Kbufsr = 0.3;
double Bufss = 0.4;
double Kbufss = 0.00025;

//Intracellular calcium flux dynamics
double Vmaxup = 0.006375;
double Kup = 0.00025;
double Vrel = 0.102;
double k1_ = 0.15;
double k2_ = 0.045;
double k3 = 0.060;
double k4 = 0.005;
double EC = 1.5;
double maxsr = 2.5;
double minsr = 1.;
double Vleak = 0.00036;
double Vxfer = 0.0038;



//Constants
double R = 8314.472;
double F = 96485.3415;
double T = 310.0;
double RTONF = (R * T) / F;

//Cellular capacitance
double CAPACITANCE = 0.185; 

//Parameters for currents
//Parameters for IKr
double Gkr = 0.153;
//Parameters for Iks
double pKNa = 0.03;
#ifdef EPI
double Gks = 0.392;
#endif
#ifdef ENDO
double Gks = 0.392;
#endif
#ifdef MCELL
double Gks = 0.098;
#endif
//Parameters for Ik1
double GK1 = 5.405;
//Parameters for Ito
#ifdef EPI
double Gto = 0.294;
#endif
#ifdef ENDO
double Gto = 0.073;
#endif
#ifdef MCELL
double Gto = 0.294;
#endif
//Parameters for INa
double GNa = 14.838;
//Parameters for IbNa
double GbNa = 0.00029;
//Parameters for INaK
double KmK = 1.0;
double KmNa = 40.0;
double knak = 2.724;
//Parameters for ICaL
double GCaL = 0.00003980;
//Parameters for IbCa
double GbCa = 0.000592;
//Parameters for INaCa
double knaca = 1000;
double KmNai = 87.5;
double KmCa = 1.38;
double ksat = 0.1;
double n = 0.35;
//Parameters for IpCa
double GpCa = 0.1238;
double KpCa = 0.0005;
//Parameters for IpK;
double GpK = 0.0146;


/*------------------------------------------------------------------------------
				PARAMETER FOR INTEGRATION
------------------------------------------------------------------------------*/
//timestep (ms)
double  HT = 0.02;

/*-----------------------------------------------------------------------------
				PARAMETERS FOR INITIAL CONDITIONS
------------------------------------------------------------------------------*/
//Initial values of state variables
double V_init = -86.2;
double Cai_init = 0.00007;
double CaSR_init = 1.3;
double CaSS_init = 0.00007;
double Nai_init = 7.67;
double Ki_init = 138.3;

/*--------------------------------------- ------------------------------------
			 PARAMETER FOR SIMULATION DURATION
  ---------------------------------------------------------------------------*/
 double STOPTIME = 100000;
/*-----------------------------------------------------------------------------
  PARAMETERS FOR STIMULATION PROTOCOLS
-----------------------------------------------------------------------------*/

#ifdef DYNRESTPROTOCOL
int i_low = 0, i_high = 1;
int j_low = 0, j_high = 1;
double stimduration = 1.0;
double stimstrength = -38;
double period = 1000;
double tbegin = 100;
double tend = tbegin + stimduration;
#endif


#ifdef S1S2RESTPROTOCOL
int i_low = 0, i_high = 1;
int j_low = 0, j_high = 1;
double stimduration = 1.;
double stimstrength = -52;
double tbegin = 100;
double tend = tbegin + stimduration;
int counter = 1;
double dia = 300;            
double basicperiod = 1000.;
double basicapd = 308;        
int repeats = 10;
#endif

/*----------------------------------------------------------------------------
							OTHER PARAMETERS
  ---------------------------------------------------------------------------*/

  //destination path to put in output files

  /*---------------------------------------------------------------------------*/


int main()
{
	static double time = 0;
	double Istim = 0;


	Variables Var(V_init, Cai_init, CaSR_init, CaSS_init, Nai_init, Ki_init);


	for (int step = 0; time <= STOPTIME; step++)
	{
		

#ifdef DYNRESTPROTOCOL

		if (time >= tbegin && time <= tend)
		{
			Istim = stimstrength;
		}
		if (time > tend)
		{
			printf("Stimulus applied time at %lf\n", floor(tbegin));
			Istim = 0.;
			tbegin = tbegin + period;
			tend = tbegin + stimduration;
		}
#endif


#ifdef S1S2RESTPROTOCOL
		if (counter < repeats)
		{

			if (time >= tbegin && time <= tend)
			{
				Istim = stimstrength;
			}
			if (time > tend)
			{
				Istim = 0.;
				tbegin = tbegin + basicperiod;
				tend = tbegin + stimduration;
				counter++;
			}

		}
		else if (counter == repeats)
		{

			if (time >= tbegin && time <= tend)
			{
				Istim = stimstrength;
			}
			if (time > tend)
			{
				Istim = 0.;
				tbegin = tbegin + basicapd + dia;       
				tend = tbegin + stimduration;
				counter++;
			}
		}
		else if (counter == repeats + 1)
		{
			if (time >= tbegin && time <= tend)
			{
				Istim = stimstrength;
			}
			if (time > tend)
			{
				Istim = 0.;
				tbegin = tbegin + basicperiod;
				tend = tbegin + stimduration;
				counter = 0;
			}
		}
#endif


		if(time >=99000 && step % 50 ==0)
			Var.writebackup(&time);

		Step(&Var, HT, &time, step, Istim);

		time += HT;
	}
	return 0;
}

