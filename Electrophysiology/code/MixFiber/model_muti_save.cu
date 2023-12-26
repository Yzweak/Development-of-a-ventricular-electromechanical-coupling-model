//
// Created by z on 23-4-19.
//
#include "model.cuh"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace std;

struct Atria_2D {
    double(*V)[YD];
    double(*Cai)[YD];
    double(*CaSR)[YD];
    double(*CaSS)[YD];
    double(*Nai)[YD];
    double(*Ki)[YD];
    double(*M)[YD];
    double(*H)[YD];
    double(*J)[YD];
    double(*Xr1)[YD];
    double(*Xr2)[YD];
    double(*Xs)[YD];
    double(*Rr)[YD];
    double(*S)[YD];
    double(*D)[YD];
    double(*Ff)[YD];
    double(*F2)[YD];
    double(*FCass)[YD];
    double(*RR)[YD];
    double(*OO)[YD];
    double(*Itot)[YD];
    double(*du)[YD];
    cell_type(*type)[YD];
    cell_location(*location)[YD];
    int(*state)[YD];

    double(*dCai)[YD];
    double(*IbCa)[YD];
    double(*IpCa)[YD];
    double(*INaCa)[YD];
    double(*Iup)[YD];
    double(*Ileak)[YD];
    double(*Ixfer)[YD];

};


void size_array(Atria_2D *(&A), int len)
{
    A = new Atria_2D[len];
}


int* splitStr(char* s)
{
    int* res = new int[4];
    int j = 0;
    for(int i=0;i<4;i++)
    {
        string ss = "";
        while(s[j] != ',' && s[j] != '\0')
        {
            ss += s[j];
            j++;
        }
        res[i] = atoi(ss.c_str());
        if(s[j] == '\r')
            break;
        j++;
    }
    return res;
}

void readMatrixFile(int** matrix)
{
    ifstream in;
    in.open(R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\ventricle_file\matrix.txt)",ios::in);
    if(!in){
        printf("File read Failed!\n");
        return;
    }
    char s[30];
    int i = 0;
    while(in.getline(s,sizeof(s)))
    {
        int* res = splitStr(s);
        for(int j=0;j<4;j++)
            matrix[i][j] = res[j];
        i++;
    }
    in.close();
}

void readFile(int* array, int mode)
{
    ifstream in;
    if(mode == 1)
        in.open(R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\ventricle_file\order.txt)",ios::in);
    else if(mode == 2)
        in.open(R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\ventricle_file\ischemia-1a.txt)",ios::in);
    else if(mode == 3)
        in.open(R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\ventricle_file\ischemia-1b.txt)",ios::in);
    else if(mode == 4)
        in.open(R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\ventricle_file\MI-short.txt)",ios::in);
    else if(mode == 5)
        in.open(R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\ventricle_file\MI-long.txt)",ios::in);
    if(!in){
        printf("File read Failed!\n");
        return;
    }
    char s[10];
    int i = 0;
    while(in.getline(s,sizeof(s)))
    {
        int j = 0;
        string ss = "";
        while (s[j] != '\0')
        {
            ss += s[j];
            j++;
        }
        array[i] = atoi(ss.c_str());
        i++;
    }
    in.close();
}

void InitCanine(Atria_2D* , const int*, const int*, const int*, const int*, const int*);
void freeCanine(Atria_2D*);
void cudaInitCanine(Atria_2D*);
void memcpy_hostToDev_Canine(Atria_2D* (&), Atria_2D* (&));
__global__ void kernel_RunCell(Atria_2D*, double);
__device__ void get_Itot(Atria_2D*, int, int, double);
void get_du(Atria_2D*,int**);
void memcpy_devToHost_Canine(Atria_2D* (&), Atria_2D* (&));
void cudaFreeCanine(Atria_2D* );
void WriteBackFile(const double* t,Atria_2D* c0);

/******************************************************************/

int main(int argc, char **argv)
{

    cudaError_t err = cudaSuccess;
    clock_t begin, finish;

    cout << endl <<  "Start simulation!" << endl << endl;
    /*****   near_matrix   ************************************************/
    int** matrix = new int*[rows];
    for(int i=0;i<rows;i++)
    {
        int* matrix_line = new int[cols];
        matrix[i] = matrix_line;
    }
    int* order = new int[stimNum];
    int* state_1a = new int[ischemia_1a_num];
    int* state_1b = new int[ischemia_1b_num];
    int* state_short = new int[MI_short_num];
    int* state_long = new int[MI_long_num];
    readMatrixFile(matrix);
    readFile(order,1);
    readFile(state_1a,2);
    readFile(state_1b,3);
    readFile(state_short,4);
    readFile(state_long,5);

    /*****   host   ************************************************/
    Atria_2D *cc;
    cc = new Atria_2D[1];
    InitCanine(cc,order,state_1a,state_1b,state_short,state_long);

    /*****  用于过渡   **********************************************/
    Atria_2D *cc_temp;
    size_array(cc_temp, 1);
    cudaInitCanine(cc_temp);

    /*****  device   ************************************************/
    Atria_2D *cc_dev;
    cudaMalloc((void**)&cc_dev, 1 * sizeof(Atria_2D));


    /***  将host端数据复制到device端  *********************************/
    memcpy_hostToDev_Canine(cc_temp, cc);
    err = cudaMemcpy(cc_dev, cc_temp, 1 * sizeof(Atria_2D), cudaMemcpyHostToDevice);

    /*****  threads and blocks  ***************************************/
    dim3 threadPerBlock(8, 8);
    dim3 blockPerGrid(8, 8);


    /**  pacing     **************************************************/
    double period = 1000;
    double tbegin = 100;
    double tend = tbegin + stimduration;
    double STOPTIME = 100000;   //100000;
    double time = 0;
    double Istim = 0;

    begin = clock();

    for (int step = 0; time <= STOPTIME; step++)
    {
        if (time >= tbegin && time <= tend)
            Istim = stimstrength;
        if (time > tend)
        {
            Istim = 0.;
            tbegin = tbegin + period;
            tend = tbegin + stimduration;
        }
        if(time >= 99000 && step % 50 == 0)
        {
            err = cudaMemcpy(cc->V, cc_temp->V, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            err = cudaMemcpy(cc->Cai, cc_temp->Cai, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            err = cudaMemcpy(cc->dCai, cc_temp->dCai, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            err = cudaMemcpy(cc->IbCa, cc_temp->IbCa, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            err = cudaMemcpy(cc->IpCa, cc_temp->IpCa, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            err = cudaMemcpy(cc->INaCa, cc_temp->INaCa, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            err = cudaMemcpy(cc->Iup, cc_temp->Iup, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            err = cudaMemcpy(cc->Ileak, cc_temp->Ileak, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            err = cudaMemcpy(cc->Ixfer, cc_temp->Ixfer, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
            WriteBackFile(&time,cc);
        }
        kernel_RunCell <<<blockPerGrid, threadPerBlock >>>(cc_dev, Istim);
        cudaDeviceSynchronize();

        err = cudaMemcpy(cc->V, cc_temp->V, XD*YD*sizeof(double), cudaMemcpyDeviceToHost); //GPU to CPU
        get_du(cc,matrix);
        err = cudaMemcpy(cc_temp->V, cc->V, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
        time += HT;
    }


    finish = clock();
    double time1 = (double)(finish - begin) / CLOCKS_PER_SEC;
    std::cout << "all simulation time (min) =" << time1 / 60 << endl << endl;

    /****** FREE-host ********/
    freeCanine(cc);
    delete[] cc;
    delete[] matrix;
    delete[] order;
    /****** FREE-device ********/
    cudaFreeCanine(cc_temp);
    delete[] cc_temp;
    cudaFree(cc_dev);
    return 0;

}



void WriteBackFile(const double* t,Atria_2D* c0)
{
    char name1[300];
    char name2[300];
    char name3[300];
    char name4[300];
    char name5[300];
    char name6[300];
    char name7[300];
    char name8[300];
    char name9[300];

    sprintf(name1,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\calcium.dat)");
    sprintf(name2,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\volt.dat)");
    sprintf(name3,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\dCai.dat)");
    sprintf(name4,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\IbCa.dat)");
    sprintf(name5,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\IpCa.dat)");
    sprintf(name6,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\INaCa.dat)");
    sprintf(name7,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\Iup.dat)");
    sprintf(name8,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\Ileak.dat)");
    sprintf(name9,"%s",R"(C:\Users\dell\Desktop\CLionProjects\MixFiber\output\Ixfer.dat)");
    ofstream oo1(name1,ios::app);
    ofstream oo2(name2,ios::app);
    ofstream oo3(name3,ios::app);
    ofstream oo4(name4,ios::app);
    ofstream oo5(name5,ios::app);
    ofstream oo6(name6,ios::app);
    ofstream oo7(name7,ios::app);
    ofstream oo8(name8,ios::app);
    ofstream oo9(name9,ios::app);
    if(!oo1){
        printf("cannot open file %s\n",name1);
        exit(1);
    }
    if(!oo2){
        printf("cannot open file %s\n",name2);
        exit(1);
    }

    oo1 << floor(*t-99000) << "\t";
    oo2 << floor(*t-99000) << "\t";
    oo3 << floor(*t-99000) << "\t";
    oo4 << floor(*t-99000) << "\t";
    oo5 << floor(*t-99000) << "\t";
    oo6 << floor(*t-99000) << "\t";
    oo7 << floor(*t-99000) << "\t";
    oo8 << floor(*t-99000) << "\t";
    oo9 << floor(*t-99000) << "\t";
    for (int j = 0;j<YD;j++)
    {
        for(int i=0;i<XD;i++)
        {
            if(j * YD + i < cellNum)
            {
                oo1 << c0->Cai[i][j] << "\t";
                oo2 << c0->V[i][j] << "\t";
                oo3 << c0->dCai[i][j] << "\t";
                oo4 << c0->IbCa[i][j] << "\t";
                oo5 << c0->IpCa[i][j] << "\t";
                oo6 << c0->INaCa[i][j] << "\t";
                oo7 << c0->Iup[i][j] << "\t";
                oo8 << c0->Ileak[i][j] << "\t";
                oo9 << c0->Ixfer[i][j] << "\t";
            }
            else
            {
                oo1 <<"\n";
                oo2 <<"\n";
                oo3 <<"\n";
                oo4 <<"\n";
                oo5 <<"\n";
                oo6 <<"\n";
                oo7 <<"\n";
                oo8 <<"\n";
                oo9 <<"\n";

                oo1.close();
                oo2.close();
                oo3.close();
                oo4.close();
                oo5.close();
                oo6.close();
                oo7.close();
                oo8.close();
                oo9.close();
                return;
            }

        }
    }

}

void InitCanine(Atria_2D* c0, const int* order, const int* state_1a, const int* state_1b, const int* state_short, const int* state_long)
{
    c0->V = new double[XD][YD];
    c0->Cai = new double[XD][YD];
    c0->CaSR = new double[XD][YD];
    c0->CaSS = new double[XD][YD];
    c0->Nai = new double[XD][YD];
    c0->Ki = new double[XD][YD];
    c0->M = new double[XD][YD];
    c0->H = new double[XD][YD];
    c0->J = new double[XD][YD];
    c0->Xr1 = new double[XD][YD];
    c0->Xr2 = new double[XD][YD];
    c0->Xs = new double[XD][YD];
    c0->Rr = new double[XD][YD];
    c0->S = new double[XD][YD];
    c0->D = new double[XD][YD];
    c0->Ff = new double[XD][YD];
    c0->F2 = new double[XD][YD];
    c0->FCass = new double[XD][YD];
    c0->RR = new double[XD][YD];
    c0->OO = new double[XD][YD];
    c0->Itot = new double[XD][YD];
    c0->du = new double[XD][YD];
    c0->type = new cell_type[XD][YD];
    c0->location = new cell_location[XD][YD];
    c0->state = new int[XD][YD];

    c0->dCai = new double[XD][YD];
    c0->IbCa = new double[XD][YD];
    c0->IpCa = new double[XD][YD];
    c0->INaCa = new double[XD][YD];
    c0->Iup = new double[XD][YD];
    c0->Ileak = new double[XD][YD];
    c0->Ixfer = new double[XD][YD];

    for (int y = 0; y < YD; y++) {
        for (int x = 0; x < XD; x++) {

            int temp = y * YD + x;

            if(temp >= cellNum)
            {
                c0->type[x][y] = Nothing;
                c0->location[x][y] = EPI;
                c0->state[x][y] = 0;
            }
            else
            {
                c0->state[x][y] = 0;
                c0->location[x][y] = EPI;
                c0->type[x][y] = PM;
                for (int i=0;i<100;i++)
                {
                    if(temp == order[i])
                    {
                        c0->type[x][y] = SAN;
                        break;
                    }
                }
                bool is_ischemia_1a = false;
                bool is_ischemia_1b = false;
                bool is_MI_short = false;

                for (int i=0;i<ischemia_1a_num;i++)
                {
                    if(temp == state_1a[i])
                    {
                        c0->state[x][y] = 1;
                        is_ischemia_1a = true;
                        break;
                    }
                }
                if(is_ischemia_1a) continue;
                for (int i=0;i<ischemia_1b_num;i++)
                {
                    if(temp == state_1b[i])
                    {
                        c0->state[x][y] = 2;
                        is_ischemia_1b = true;
                        break;
                    }
                }
                if(is_ischemia_1b) continue;
                for (int i=0;i<MI_short_num;i++)
                {
                    if(temp == state_short[i])
                    {
                        c0->state[x][y] = 3;
                        is_MI_short = true;
                        break;
                    }
                }
                if(is_MI_short) continue;
                for (int i=0;i<MI_long_num;i++)
                {
                    if(temp == state_long[i])
                    {
                        c0->state[x][y] = 4;
                        break;
                    }
                }
            }
        }
    }

    for (int y = 0; y < YD; y++)
    {
        for (int x = 0; x < XD; x++)
        {
            //assign cell state
            if(c0->state[x][y] == 0 || c0->state[x][y] == 3 || c0->state[x][y] == 4)
            {
                c0->V[x][y] = -86.2;
                c0->Cai[x][y] = 0.00007;
                c0->CaSR[x][y] = 1.3;
                c0->CaSS[x][y] = 0.00007;
                c0->Nai[x][y] = 7.67;
                c0->Ki[x][y] = 138.3;
                c0->M[x][y] = 0.0;
                c0->H[x][y] = 0.75;
                c0->J[x][y] = 0.75;
                c0->Xr1[x][y] = 0.0;
                c0->Xr2[x][y] = 1.0;
                c0->Xs[x][y] = 0.0;
                c0->Rr[x][y] = 0.0;
                c0->S[x][y] = 1.0;
                c0->D[x][y] = 0.0;
                c0->Ff[x][y] = 1.0;
                c0->F2[x][y] = 1.0;
                c0->FCass[x][y] = 1.0;
                c0->RR[x][y] = 1.0;
                c0->OO[x][y] = 0.0;
                c0->Itot[x][y] = 0.0;
                c0->du[x][y] = 0.0;
            }
            else if(c0->state[x][y] == 1)
            {
                c0->V[x][y] = -70.342786;
                c0->Cai[x][y] = 0.000095738714626814048;
                c0->CaSR[x][y] = 2.874115913016119800000;
                c0->CaSS[x][y] = 0.000302877364231207700;
                c0->Nai[x][y] = 6.525357256767184900000;
                c0->Ki[x][y] = 141.870292899352390000000;
                c0->M[x][y] = 0.033655635139547148000;
                c0->H[x][y] = 0.211053852569920620000;
                c0->J[x][y] = 0.180147754598375870000;
                c0->Xr1[x][y] = 0.024098542643841897000;
                c0->Xr2[x][y] = 0.323941263226489200000;
                c0->Xs[x][y] = 0.009339544433470185100;
                c0->Rr[x][y] = 0.000000288926510767676;
                c0->S[x][y] = 0.999957603371733010000;
                c0->D[x][y] = 0.000164510100997893920;
                c0->Ff[x][y] = 0.980144782112325140000;
                c0->F2[x][y] = 0.995728345214580090000;
                c0->FCass[x][y] = 0.999969154239247040000;
                c0->RR[x][y] = 0.989633139742973960000;
                c0->OO[x][y] = 0.000000171795851160072;
                c0->Itot[x][y] = 0.0;
                c0->du[x][y] = 0.0;
            }
            else if(c0->state[x][y] == 2)
            {
                c0->V[x][y] = -69.756526;
                c0->Cai[x][y] = 0.000110191555337951770;
                c0->CaSR[x][y] = 2.899895833813058000000;
                c0->CaSS[x][y] = 0.000261707275030859110;
                c0->Nai[x][y] = 8.443153663718010900000;
                c0->Ki[x][y] = 139.046086118217600000000;
                c0->M[x][y] = 0.037396266233898176000;
                c0->H[x][y] = 0.193531112266500660000;
                c0->J[x][y] = 0.170193474094604810000;
                c0->Xr1[x][y] = 0.013720983480392804000;
                c0->Xr2[x][y] = 0.318614799111897550000;
                c0->Xs[x][y] = 0.009711833626984207200;
                c0->Rr[x][y] = 0.000000318575104674730;
                c0->S[x][y] = 0.999952334011770310000;
                c0->D[x][y] = 0.000177883400051902350;
                c0->Ff[x][y] = 0.988991405078883500000;
                c0->F2[x][y] = 0.995358028564252040000;
                c0->FCass[x][y] = 0.999978028985132990000;
                c0->RR[x][y] = 0.992241352890124500000;
                c0->OO[x][y] = 0.000000129042183903380;
                c0->Itot[x][y] = 0.0;
                c0->du[x][y] = 0.0;
            }
        }
    }
}

void cudaInitCanine(Atria_2D* c0)
{
    cudaError_t err = cudaMalloc((void **)&c0->V, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Cai, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->CaSR, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->CaSS, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Nai, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Ki, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->M, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->H, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->J, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Xr1, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Xr2, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Xs, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Rr, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->S, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->D, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Ff, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->F2, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->FCass, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->RR, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->OO, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Itot, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->du, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->type, sizeof(cell_type)* XD*YD);
    err = cudaMalloc((void **)&c0->location, sizeof(cell_location)* XD*YD);
    err = cudaMalloc((void **)&c0->state, sizeof(int)* XD*YD);

    err = cudaMalloc((void **)&c0->dCai, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->IbCa, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->IpCa, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->INaCa, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Iup, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Ileak, sizeof(double)* XD*YD);
    err = cudaMalloc((void **)&c0->Ixfer, sizeof(double)* XD*YD);

    if (err != cudaSuccess)
    {
        std::cerr << err << endl;
        // fprintf(stderr, "Failed to allocate device in memcpy_hostToDev_Canine\n", cudaGetErrorString(err));
        printf("Failed to allocate device in cudaInitCanine\n");
        exit(EXIT_FAILURE);
    }
}

void memcpy_hostToDev_Canine(Atria_2D* (&c0_dev), Atria_2D* (&c0))
{
    cudaError_t err = cudaMemcpy(c0_dev->V, c0->V, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Cai, c0->Cai, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->CaSR, c0->CaSR, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->CaSS, c0->CaSS, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Nai, c0->Nai, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Ki, c0->Ki, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->M, c0->M, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->H, c0->H, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->J, c0->J, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Xr1, c0->Xr1, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Xr2, c0->Xr2, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Xs, c0->Xs, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Rr, c0->Rr, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->S, c0->S, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->D, c0->D, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Ff, c0->Ff, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->F2, c0->F2, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->FCass, c0->FCass, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->RR, c0->RR, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->OO, c0->OO, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Itot, c0->Itot, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->du, c0->du, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->type, c0->type, XD*YD*sizeof(cell_type), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->location, c0->location, XD*YD*sizeof(cell_location), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->state, c0->state, XD*YD*sizeof(int), cudaMemcpyHostToDevice);


    err = cudaMemcpy(c0_dev->dCai, c0->dCai, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->IbCa, c0->IbCa, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->IpCa, c0->IpCa, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->INaCa, c0->INaCa, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Iup, c0->Iup, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Ileak, c0->Ileak, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    err = cudaMemcpy(c0_dev->Ixfer, c0->Ixfer, XD*YD*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        std::cerr << err << endl;
        // fprintf(stderr, "Failed to allocate device in memcpy_hostToDev_Canine\n", cudaGetErrorString(err));
        printf("Failed to allocate device in memcpy_hostToDev_Canine\n");
        exit(EXIT_FAILURE);
    }
}

void memcpy_devToHost_Canine(Atria_2D* (&c0), Atria_2D* (&c0_dev))
{
    cudaError_t err = cudaMemcpy(c0->V, c0_dev->V, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Cai, c0_dev->Cai, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->CaSR, c0_dev->CaSR, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->CaSS, c0_dev->CaSS, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Nai, c0_dev->Nai, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Ki, c0_dev->Ki, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->M, c0_dev->M, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->H, c0_dev->H, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->J, c0_dev->J, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Xr1, c0_dev->Xr1, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Xr2, c0_dev->Xr2, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Xs, c0_dev->Xs, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Rr, c0_dev->Rr, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->S, c0_dev->S, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->D, c0_dev->D, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Ff, c0_dev->Ff, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->F2, c0_dev->F2, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->FCass, c0_dev->FCass, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->RR, c0_dev->RR, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->OO, c0_dev->OO, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Itot, c0_dev->Itot, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->du, c0_dev->du, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->type, c0_dev->type, XD*YD*sizeof(cell_type), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->location, c0_dev->location, XD*YD*sizeof(cell_location), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->state, c0_dev->state, XD*YD*sizeof(int), cudaMemcpyDeviceToHost);

    err = cudaMemcpy(c0->dCai, c0_dev->dCai, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->IbCa, c0_dev->IbCa, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->IpCa, c0_dev->IpCa, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->INaCa, c0_dev->INaCa, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Iup, c0_dev->Iup, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Ileak, c0_dev->Ileak, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(c0->Ixfer, c0_dev->Ixfer, XD*YD*sizeof(double), cudaMemcpyDeviceToHost);


}

void get_du(Atria_2D* c0,int** matrix)
{
    for (int j=0;j<YD;j++)
    {
        for (int i=0;i<XD;i++)
        {
            int index = j * YD + i;
            if(index < cellNum)
            {
                int index1 = matrix[index][0];
                int index2 = matrix[index][1];
                int index3 = matrix[index][2];
                int index4 = matrix[index][3];
                if(index3 == -1)
                    c0->du[i][j] = Dlong * (c0->V[index1%YD][index1/YD] + c0->V[index2%YD][index2/YD] - 2*c0->V[index%YD][index/YD]) /  (dx*dx);
                else
                {
                    if(index4 == -1)
                        c0->du[i][j] = Dlong * (c0->V[index1%YD][index1/YD] + c0->V[index2%YD][index2/YD] + c0->V[index3%YD][index3/YD] - 3*c0->V[index%YD][index/YD]) /  (dx*dx);
                    else
                        c0->du[i][j] = Dlong * (c0->V[index1%YD][index1/YD] + c0->V[index2%YD][index2/YD] + c0->V[index3%YD][index3/YD] + c0->V[index4%YD][index4/YD] - 4*c0->V[index%YD][index/YD]) / (dx*dx);
                }
                if(c0->state[i][j] == 2) c0->du[i][j] *= 0.5;
                c0->V[i][j] =  c0->V[i][j] + HT * c0->du[i][j];
            }
            else
            {
                c0->du[i][j] = 0;
            }
        }
    }
}

__global__ void kernel_RunCell(Atria_2D* c0, double Istim)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    if (j * YD + i >= cellNum) //越界保护
        return;

    /****** Voltage update ********/
    get_Itot(c0, i, j, Istim);//已在刺激区域加入Istim
    c0->V[i][j] = c0->V[i][j] - HT * c0->Itot[i][j];
}

__device__ void get_Itot(Atria_2D* s0, int i, int j, double Istim)
{
    double sm = s0->M[i][j];
    double sh = s0->H[i][j];
    double sj = s0->J[i][j];
    double sxr1 = s0->Xr1[i][j];
    double sxr2 = s0->Xr2[i][j];
    double sxs = s0->Xs[i][j];
    double ss = s0->S[i][j];
    double sr = s0->Rr[i][j];
    double sd = s0->D[i][j];
    double sf = s0->Ff[i][j];
    double sf2 = s0->F2[i][j];
    double sfcass = s0->FCass[i][j];
    double sRR = s0->RR[i][j];
    double sOO = s0->OO[i][j];
    double svolt = s0->V[i][j];
    double Cai = s0->Cai[i][j];
    double CaSR = s0->CaSR[i][j];
    double CaSS = s0->CaSS[i][j];
    double Nai = s0->Nai[i][j];
    double Ki = s0->Ki[i][j];
    cell_type type = s0->type[i][j];
    cell_location location = s0->location[i][j];
    int state = s0->state[i][j];

    double Ko = 0;
    if (state == 0 || state == 4)
        Ko = 5.4;
    else
        Ko = 8.0;
    //Needed to compute currents
    double Ek = RTONF * (log((Ko / Ki)));
    double Ena = RTONF * (log((Nao / Nai)));
    double Eks = RTONF * (log((Ko + pKNa * Nao) / (Ki + pKNa * Nai)));
    double Eca = 0.5 * RTONF * (log((Cao / Cai)));
    double Ak1 = 0.1 / (1. + exp(0.06 * (svolt - Ek - 200)));
    double Bk1 = (3. * exp(0.0002 * (svolt - Ek + 100)) + exp(0.1 * (svolt - Ek - 10))) / (1. + exp(-0.5 * (svolt - Ek)));
    double rec_iK1 = Ak1 / (Ak1 + Bk1);
    double rec_iNaK = (1. / (1. + 0.1245 * exp(-0.1 * svolt * F / (R * T)) + 0.0353 * exp(-svolt * F / (R * T))));
    double rec_ipK = 1. / (1. + exp((25 - svolt) / 5.98));
    //Compute currents depend on cell states
    double cof_INa=1.0, cof_ICaL=1.0, cof_Ito=1.0, cof_IKr=1.0, cof_IKs=1.0, cof_INaCa=1.0, cof_INaK=1.0, cof_IbCa=1.0;
    switch (state)
    {
        case 1: // Ischemia-1a
            cof_INa = 0.887;
            cof_ICaL = 0.8;
            cof_Ito = 0.5;
            cof_IKs = 0.781;
            break;
        case 2: // Ischemia-1b
            cof_ICaL = 0.5;
            cof_INaCa = 0.6;
            cof_INaK = 0.46;
            cof_IbCa = 1.3;
            break;
        case 3: // MI-short
            cof_INa = 0.38;
            cof_ICaL = 0.62;
            cof_Ito = 0.37;
            cof_IKr = 0.30;
            cof_IKs = 0.20;
            break;
        case 4: // MI-long
            cof_INa = 0.38;
            cof_ICaL = 0.31;
            cof_Ito = 0.5;
            cof_IKr = 0.30;
            cof_IKs = 0.20;
            break;
        default:
            break;
    }
    double Gks,Gto;
    if (location == EPI)
    {
        Gks = 0.392;
        Gto = 0.294;
    }
    else if (location == ENDO)
    {
        Gks = 0.392;
        Gto = 0.073;
    }
    else
    {
        Gks = 0.098;
        Gto = 0.294;
    }
    double INa = cof_INa * GNa * sm * sm * sm * sh * sj * (svolt - Ena);//Zhai change
    double ICaL = cof_ICaL * GCaL * sd * sf * sf2 * sfcass * 4 * (svolt - 15) * (F * F / (R * T)) * (0.25 * exp(2 * (svolt - 15) * F / (R * T)) * CaSS - Cao) / (exp(2 * (svolt - 15) *             F / (R * T)) - 1.);
    double Ito = cof_Ito * Gto * sr * ss * (svolt - Ek);
    double IKr = cof_IKr * Gkr * sqrt(Ko / 5.4) * sxr1 * sxr2 * (svolt - Ek);
    double IKs = cof_IKs * Gks * sxs * sxs * (svolt - Eks);
    double IK1 = GK1 * rec_iK1 * (svolt - Ek);
    double INaCa = cof_INaCa * knaca * (1. / (KmNai * KmNai * KmNai + Nao * Nao * Nao)) * (1. / (KmCa + Cao)) * (1. / (1 + ksat * exp((n - 1) * svolt * F / (R * T)))) *
                   (exp(n * svolt * F / (R * T)) * Nai * Nai * Nai * Cao - exp((n - 1) * svolt * F / (R * T)) * Nao * Nao * Nao * Cai * 2.5);
    double INaK = cof_INaK * knak * (Ko / (Ko + KmK)) * (Nai / (Nai + KmNa)) * rec_iNaK;
    double IpCa = GpCa * Cai / (KpCa + Cai);
    double IpK = GpK * rec_ipK * (svolt - Ek);
    double IbNa = GbNa * (svolt - Ena);
    double IbCa = cof_IbCa * GbCa * (svolt - Eca);

    // 补充IKATP电流 根据：2011_Experiment-model interaction for analysis of epicardial activation during human ventricular fibrillation with global myocardial ischaemia
    double GKATP = 3.9;
    double H = 2.0;
    double nn = 0.24;
    double ATPi = 4.6;
    double Khalf = 0.25;
    double IKATP = GKATP * (1 / (1 + pow(ATPi / Khalf, H))) * pow(Ko / 5.4, nn) * (svolt - Ek);
    //将ORD模型里的INaL加入初始化INaL用到的参数
    double mL = 0;
    double hL = 1;
    double mLss = 1.0 / (1.0 + exp((-(svolt + 42.85)) / 5.264));
    double tmL = 1.0 / (6.765 * exp((svolt + 11.64) / 34.77) + 8.552 * exp(-(svolt + 77.42) / 5.955));
    mL = mLss - (mLss - mL) * exp(-HT / tmL);
    double hLss = 1.0 / (1.0 + exp((svolt + 87.61) / 7.488));
    double thL = 200.0;
    hL = hLss - (hLss - hL) * exp(-HT / thL);
    double cof_NaL = state == 0 ? 1.5 : 1;
    double INaL = cof_NaL * GNaL * (svolt - Ena) * mL * hL;

    switch (state)
    {
        case 0: // Normal-Cell
            IKATP = 0.0;
            INaL = 0.0;
            break;
        case 1: // Ischemia-1a
            INaL *= 1.5;
            break;
        case 2: // Ischemia-1b
            break;
        case 3: // MI-short
            break;
        case 4: // MI-long
            IKATP = 0.0;
            INaL = 0.0;
            break;
        default:
            break;
    }

    if(type == SAN)
        s0->Itot[i][j] = IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + IKATP + INaL + Istim;
    else
        s0->Itot[i][j] = IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + IKATP + INaL;
    //update concentrations
    double kCaSR = maxsr - ((maxsr - minsr) / (1 + (EC / CaSR) * (EC / CaSR)));
    double k1 = k1_ / kCaSR;
    double k2 = k2_ * kCaSR;
    double dRR = k4 * (1 - sRR) - k2 * CaSS * sRR;
    s0->RR[i][j] += HT * dRR;
    s0->OO[i][j] = k1 * CaSS * CaSS * sRR / (k3 + k1 * CaSS * CaSS);


    double Irel = Vrel * sOO * (CaSR - CaSS);
    double Ileak = Vleak * (CaSR - Cai);
    double Iup = Vmaxup / (1. + ((Kup * Kup) / (Cai * Cai)));
    double Ixfer = Vxfer * (CaSS - Cai);


    double CaCSQN = Bufsr * CaSR / (CaSR + Kbufsr);
    double dCaSR = HT * (Iup - Irel - Ileak);
    double bjsr = Bufsr - CaCSQN - dCaSR - CaSR + Kbufsr;
    double cjsr = Kbufsr * (CaCSQN + dCaSR + CaSR);
    s0->CaSR[i][j] = (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2;

    double inverseVcF2 = 1 / (2 * Vc * F);
    double inverseVcF = 1. / (Vc * F);
    double inversevssF2 = 1 / (2 * Vss * F);
    double CaSSBuf = Bufss * CaSS / (CaSS + Kbufss);
    double dCaSS = HT * (-Ixfer * (Vc / Vss) + Irel * (Vsr / Vss) + (-ICaL * inversevssF2 * CAPACITANCE));
    double bcss = Bufss - CaSSBuf - dCaSS - CaSS + Kbufss;
    double ccss = Kbufss * (CaSSBuf + dCaSS + CaSS);
    s0->CaSS[i][j] = (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2;


    double CaBuf = Bufc * Cai / (Cai + Kbufc);
    double dCai = HT * ((-(IbCa + IpCa - 2 * INaCa) * inverseVcF2 * CAPACITANCE) - (Iup - Ileak) * (Vsr / Vc) + Ixfer);
    double bc = Bufc - CaBuf - dCai - Cai + Kbufc;
    double cc = Kbufc * (CaBuf + dCai + Cai);
    s0->Cai[i][j] = (sqrt(bc * bc + 4 * cc) - bc) / 2;

    double dNai = 0.0;
    double dKi = 0.0;
    if (state == 1)
    {
        dNai = -(INa + INaL + IbNa + 3 * INaK + 3 * INaCa) * inverseVcF * CAPACITANCE;
        dKi = -(Istim + IK1 + Ito + IKr + IKs + IKATP - 2 * INaK + IpK) * inverseVcF * CAPACITANCE;
    }
    else
    {
        dNai = -(INa + IbNa + 3 * INaK + 3 * INaCa) * inverseVcF * CAPACITANCE;
        dKi = -(Istim + IK1 + Ito + IKr + IKs - 2 * INaK + IpK) * inverseVcF * CAPACITANCE;
    }

    s0->Nai[i][j] += HT * dNai;
    s0->Ki[i][j] += HT * dKi;
    //compute steady state values and time constants
    double AM = 1. / (1. + exp((-60. - svolt) / 5.));
    double BM = 0.1 / (1. + exp((svolt + 35.) / 5.)) + 0.10 / (1. + exp((svolt - 50.) / 200.));
    double TAU_M = AM * BM;
    double M_INF = 1. / ((1. + exp((-55.5 - svolt) / 9.03)) * (1. + exp((-55.5 - svolt) / 9.03)));
    double AH_1 = 0;
    double AH_2 = 0;
    double AJ_1 = 0;
    double AJ_2 = 0;
    double BH_1 = 0;
    double BH_2 = 0;
    double BJ_1 = 0;
    double BJ_2 = 0;
    double TAU_H = 0;
    double TAU_J = 0;
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
    double H_INF = 1. / ((1. + exp((svolt + 71.55) / 7.43)) * (1. + exp((svolt + 71.55) / 7.43)));
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
    double J_INF = H_INF;
    double Xr1_INF = 1. / (1. + exp((-26. - svolt) / 7.));
    double axr1 = 450. / (1. + exp((-45. - svolt) / 10.));
    double bxr1 = 6. / (1. + exp((svolt - (-30.)) / 11.5));
    double TAU_Xr1 = axr1 * bxr1;
    double Xr2_INF = 1. / (1. + exp((svolt - (-88.)) / 24.));
    double axr2 = 3. / (1. + exp((-60. - svolt) / 20.));
    double bxr2 = 1.12 / (1. + exp((svolt - 60.) / 20.));
    double TAU_Xr2 = axr2 * bxr2;

    double Xs_INF = 1. / (1. + exp((-5. - svolt) / 14.));
    double Axs = (1400. / (sqrt(1. + exp((5. - svolt) / 6))));
    double Bxs = (1. / (1. + exp((svolt - 35.) / 15.)));
    double TAU_Xs = Axs * Bxs + 80;

    double R_INF,S_INF,TAU_R,TAU_S;
    if (location == EPI)
    {
        if (state == 1)
        {
            R_INF = 1. / (1. + exp((27.2 - svolt) / 6.));
            S_INF = 1. / (1. + exp((svolt + 6.3) / 5.));
        }
        else
        {
            R_INF = 1. / (1. + exp((20 - svolt) / 6.));
            S_INF = 1. / (1. + exp((svolt + 20) / 5.));
        }
        TAU_R = 9.5 * exp(-(svolt + 40.) * (svolt + 40.) / 1800.) + 0.8;
        TAU_S = 85. * exp(-(svolt + 45.) * (svolt + 45.) / 320.) + 5. / (1. + exp((svolt - 20.) / 5.)) + 3.;
    }
    else if (location == ENDO)
    {
        R_INF = 1. / (1. + exp((20 - svolt) / 6.));
        S_INF = 1. / (1. + exp((svolt + 28) / 5.));
        TAU_R = 9.5 * exp(-(svolt + 40.) * (svolt + 40.) / 1800.) + 0.8;
        TAU_S = 1000. * exp(-(svolt + 67) * (svolt + 67) / 1000.) + 8.;
    }
    else
    {
        R_INF = 1. / (1. + exp((20 - svolt) / 6.));
        S_INF = 1. / (1. + exp((svolt + 20) / 5.));
        TAU_R = 9.5 * exp(-(svolt + 40.) * (svolt + 40.) / 1800.) + 0.8;
        TAU_S = 85. * exp(-(svolt + 45.) * (svolt + 45.) / 320.) + 5. / (1. + exp((svolt - 20.) / 5.)) + 3.;
    }

    double D_INF = 1. / (1. + exp((-8 - svolt) / 7.5));
    double Ad = 1.4 / (1. + exp((-35 - svolt) / 13)) + 0.25;
    double Bd = 1.4 / (1. + exp((svolt + 5) / 5));
    double Cd = 1. / (1. + exp((50 - svolt) / 20));
    double TAU_D = Ad * Bd + Cd;
    double F_INF = 1. / (1. + exp((svolt + 20) / 7));
    double Af = 1102.5 * exp(-(svolt + 27) * (svolt + 27) / 225);
    double Bf = 200. / (1 + exp((13 - svolt) / 10.));
    double Cf = (180. / (1 + exp((svolt + 30) / 10))) + 20;
    double TAU_F = Af + Bf + Cf;
    double F2_INF = 0.67 / (1. + exp((svolt + 35) / 7)) + 0.33;
    double Af2 = 600 * exp(-(svolt + 25) * (svolt + 25) / 170);
    double Bf2 = 31 / (1. + exp((25 - svolt) / 10));
    double Cf2 = 16 / (1. + exp((svolt + 30) / 10));
    double TAU_F2 = Af2 + Bf2 + Cf2;
    double FCaSS_INF = 0.6 / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 0.4;
    double TAU_FCaSS = 80. / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 2.;
    //Update gates
    s0->M[i][j] = M_INF - (M_INF - sm) * exp(-HT / TAU_M);
    s0->H[i][j] = H_INF - (H_INF - sh) * exp(-HT / TAU_H);
    s0->J[i][j] = J_INF - (J_INF - sj) * exp(-HT / TAU_J);
    s0->Xr1[i][j] = Xr1_INF - (Xr1_INF - sxr1) * exp(-HT / TAU_Xr1);
    s0->Xr2[i][j] = Xr2_INF - (Xr2_INF - sxr2) * exp(-HT / TAU_Xr2);
    s0->Xs[i][j] = Xs_INF - (Xs_INF - sxs) * exp(-HT / TAU_Xs);
    s0->S[i][j] = S_INF - (S_INF - ss) * exp(-HT / TAU_S);
    s0->Rr[i][j] = R_INF - (R_INF - sr) * exp(-HT / TAU_R);
    s0->D[i][j] = D_INF - (D_INF - sd) * exp(-HT / TAU_D);
    s0->Ff[i][j] = F_INF - (F_INF - sf) * exp(-HT / TAU_F);
    s0->F2[i][j] = F2_INF - (F2_INF - sf2) * exp(-HT / TAU_F2);
    s0->FCass[i][j] = FCaSS_INF - (FCaSS_INF - sfcass) * exp(-HT / TAU_FCaSS);

    s0->dCai[i][j] = dCai;
    s0->IbCa[i][j] = IbCa;
    s0->IpCa[i][j] = IpCa;
    s0->INaCa[i][j] = INaCa;
    s0->Iup[i][j] = Iup;
    s0->Ileak[i][j] = Ileak;
    s0->Ixfer[i][j] = Ixfer;

}

void freeCanine(Atria_2D* c0)
{
    delete[] c0->V;
    delete[] c0->Cai;
    delete[] c0->CaSR;
    delete[] c0->CaSS;
    delete[] c0->Nai;
    delete[] c0->Ki;
    delete[] c0->M;
    delete[] c0->H;
    delete[] c0->J;
    delete[] c0->Xr1;
    delete[] c0->Xr2;
    delete[] c0->Xs;
    delete[] c0->Rr;
    delete[] c0->S;
    delete[] c0->D;
    delete[] c0->Ff;
    delete[] c0->F2;
    delete[] c0->FCass;
    delete[] c0->RR;
    delete[] c0->OO;
    delete[] c0->Itot;
    delete[] c0->du;
    delete[] c0->type;
    delete[] c0->location;
    delete[] c0->state;

    delete[] c0->dCai;
    delete[] c0->IbCa;
    delete[] c0->IpCa;
    delete[] c0->INaCa;
    delete[] c0->Iup;
    delete[] c0->Ileak;
    delete[] c0->Ixfer;
}

void cudaFreeCanine(Atria_2D* c0)
{
    cudaFree(c0->V);
    cudaFree(c0->Cai);
    cudaFree(c0->CaSR);
    cudaFree(c0->CaSS);
    cudaFree(c0->Nai);
    cudaFree(c0->Ki);
    cudaFree(c0->M);
    cudaFree(c0->H);
    cudaFree(c0->J);
    cudaFree(c0->Xr1);
    cudaFree(c0->Xr2);
    cudaFree(c0->Xs);
    cudaFree(c0->Rr);
    cudaFree(c0->S);
    cudaFree(c0->D);
    cudaFree(c0->Ff);
    cudaFree(c0->F2);
    cudaFree(c0->FCass);
    cudaFree(c0->RR);
    cudaFree(c0->OO);
    cudaFree(c0->Itot);
    cudaFree(c0->du);
    cudaFree(c0->type);
    cudaFree(c0->location);
    cudaFree(c0->state);


    cudaFree(c0->dCai);
    cudaFree(c0->IbCa);
    cudaFree(c0->IpCa);
    cudaFree(c0->INaCa);
    cudaFree(c0->Iup);
    cudaFree(c0->Ileak);
    cudaFree(c0->Ixfer);

}
