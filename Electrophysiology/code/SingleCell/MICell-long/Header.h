/* Header file, forms a generalized header file containing headers of libraries that are used,
   function definitions and the cell type specification */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>

#include "Variables.h"

#define f01 0.15
#define f12 0.5
#define f23 0.35
#define g01 0.1
#define g12 0.2
#define g23 0.3
#define kpn 0.04
#define ntrop 5.525
#define SL 2.15
#define alpha_SL 0.03571428571

#define LTRPNtot 70E-3
#define Ktrop_Ca 40E-3 / 40
#define Ktrop_half 1.0 / (1.0 + (Ktrop_Ca / (1.7 / 1000.0 + ((0.9 / 1000.0 - 1.7 / 1000.0) / (2.3 - 1.7)) * (SL - 1.7))))
#define inv_LTRPNtot_Ktrop_half 1.0 / (LTRPNtot * Ktrop_half)
   /*------------------------------------------------------------------------------
   FLAG TO CHOOSE BETWEEN EPICARDIAL ENDOCARDIAL AND MIDMYOCARDIAL CELL TYPES
   ------------------------------------------------------------------------------*/
#define EPI
   //#define ENDO
   //#define MCELL


void Step(Variables* V, double HT, double* t, int step, double Istim);


