#include <stdio.h>
#include <math.h>
#include <complex>

#define M_PI       3.14159265358979323846

const int Nm = 120;

double dt = 0.1;
double ts = 1.0;
int tperiod = 500;
int nt = 1;
double f0 = 0.1;
double omegaQ = 6.0;
double omega = 6.001; // field frequency qubit
double gammaF_0 = 0.001;
double gammaE_0 = 0.00065;
double gammaO_0 = 0.01;
double Omega = 1.; // field frequency oscillator
double Ams = 0.02;
double mu = 0.0001;
double lambda = 0.02;
double omegaJ = 1.02;

using namespace std;

unsigned long long state = 22140;
double d_rdistr(unsigned long long & state) {
	state *= 302875106592253ULL;
	return double((state & ((1ULL << 59) - 1)) / double(1ULL << 59));
}


void d_RK2(complex<double>*HH, double Am, double f, complex<double> * c, double t, complex<double> * result, double ddt) {
	complex<double> I(0, 1);
	complex<double> AmCosOmegaTT05 = Am*cos(omega*t)*0.5;
	complex<double> fCosOmegaTT = f*cos(Omega*t);

	result[0] = HH[5*(2*0+0)+2] * c[0] + AmCosOmegaTT05  * c[1] + HH[5*(2*0+0)+4] * fCosOmegaTT * c[2];
	result[1] = AmCosOmegaTT05  * c[0] + HH[5*(2*0+1)+2] * c[1] + HH[5*(2*0+1)+4] * fCosOmegaTT * c[3];
	for (int i = 1; i < Nm - 1; i++){
		result[2*i+0] = HH[5*(2*i+0)+0] * fCosOmegaTT * c[(2*i+0)-2] /*+ AmCosOmegaTT05  * c[(2*i+0)-1]*/ + HH[5*(2*i+0)+2] * c[(2*i+0)] +   AmCosOmegaTT05  * c[(2*i+0)+1]   + HH[5*(2*i+0)+4] * fCosOmegaTT * c[(2*i+0)+2];
		result[2*i+1] = HH[5*(2*i+1)+0] * fCosOmegaTT * c[(2*i+1)-2]   + AmCosOmegaTT05  * c[(2*i+1)-1]   + HH[5*(2*i+1)+2] * c[(2*i+1)] /*+ AmCosOmegaTT05  * c[(2*i+1)+1]*/ + HH[5*(2*i+1)+4] * fCosOmegaTT * c[(2*i+1)+2];
	}
	result[2*(Nm-1)+0] = HH[5*(2*(Nm-1)+0)+0] * fCosOmegaTT * c[(2*(Nm-1)+0)-2] + HH[5*(2*(Nm-1)+0)+2] * c[(2*(Nm-1)+0)] + AmCosOmegaTT05  * c[(2*(Nm-1)+0)+1];
	result[2*(Nm-1)+1] = HH[5*(2*(Nm-1)+1)+0] * fCosOmegaTT * c[(2*(Nm-1)+1)-2] + AmCosOmegaTT05  * c[(2*(Nm-1)+1)-1]  + HH[5*(2*(Nm-1)+1)+2] * c[(2*(Nm-1)+1)];

	for (int i = 0; i < 2 * Nm; i++)
	{
		result[i] *= -I * complex<double>(ddt);
	}
}


int main(void){
    FILE *fout1=fopen("nr_m(t)_stat_Nm=120_f0=0.1_wq=6_w=6.01_gf=0.001_ge=0.00065_go=0.01_Omega=1_Ams=0.02_mu=0.0001_lm=0.02_wJ=1_t=500T_nrealiz=100_puls=PI_t1=500T.txt","w");
	FILE *fout2=fopen("nr_m(t)_stat_Nm=120_f0=0.1_wq=6_w=6.01_gf=0.001_ge=0.00065_go=0.01_Omega=1_Ams=0.02_mu=0.0001_lm=0.02_wJ=1_t=500T_nrealiz=100_puls=PI_t1=220T.txt","w");
	FILE *fout3=fopen("nr_m(t)_stat_Nm=120_f0=0.1_wq=6_w=6.01_gf=0.001_ge=0.00065_go=0.01_Omega=1_Ams=0.02_mu=0.0001_lm=0.02_wJ=1_t=500T_nrealiz=100_puls=PI_t1=180T_500T.txt","w");
  //FILE *fq1out=fopen("qp1_Nm=120_f0=0.2_wq=5_w=5.01_gf=0.0001_ge=0.0001_go=0.02_Omega=1_Ams=0.02_mu=0.0001_lm=0.02_wJ=1.03.txt","w");
  //FILE *fq2out=fopen("qp2_Nm=120_f0=0.2_wq=5_w=5.01_gf=0.0001_ge=0.0001_go=0.02_Omega=1_Ams=0.02_mu=0.0001_lm=0.02_wJ=1.03.txt","w");
	FILE *foutnm=fopen("TESTOsc_stat_Nm=120_f0=0.1_wq=6_w=6.001_gf=0.001_ge=0.00065_go=0.01_Omega=1_Ams=0.02_mu=0.0001_lm=0.02_wJ=1.02_t=500T_nrealiz=100_puls=PI_t1=500T.txt","w");
	FILE *foutqb=fopen("TESTQb_stat_Nm=120_f0=0.1_wq=6_w=6.001_gf=0.001_ge=0.00065_go=0.01_Omega=1_Ams=0.02_mu=0.0001_lm=0.02_wJ=1.02_t=500T_nrealiz=100_puls=PI_t1=500T.txt","w");

  int JumpsCount = 0;

  int Nrealiz = 100;

dt = 0.1;
ts = 1.0;
tperiod = 500;
nt = 1;
f0 = 0.1;
omegaQ = 6.0;
omega = 6.001; // field frequency qubit
gammaF_0 = 0.001;
gammaE_0 = 0.00065;
gammaO_0 = 0.01;
Omega = 1.; // field frequency oscillator
Ams = 0.02;
mu = 0.0001;
lambda = 0.02;
omegaJ = 1.02;

double omegaRabi = sqrt(pow(omegaQ - omega, 2) + pow(Ams, 2));
double period = 2*M_PI/omegaRabi;

  int i, j, k;
	complex<double> I(0, 1);
	complex<double> *host_HH = new complex<double>[5 * 2 * Nm];
    complex<double> *C = new complex<double>[2 * Nm];

  double Nstep = (2*M_PI)/dt;
	int tmax = floor(Nstep)*tperiod;
	double ddt = dt/30;

  double *mas = new double[Nrealiz];
  double *nr_m = new double[tmax];
  double *pq= new double[tmax];
  complex<double> b[2*Nm];
	double a[Nm];
	double a1[Nm];
	double a2[Nm];

  complex<double> k1[2 * Nm], k2[2 * Nm], k3[2 * Nm], k4[2 * Nm], tmp[2 * Nm];

  for(int i = 0; i < tmax; i++){
	nr_m[i] = 0;
	pq[i] = 0;
  }
  for(int ni = 1; ni < Nrealiz+1; ni++) {
		
	   int realiz = ni;
       printf("ni = %d\n", realiz);
	   int k0 = 11304;
  
	for (int i = 0; i < 2 * Nm; ++i)
		C[i] = complex<double>(0, 0);

	C[2*(nt-1)+0] = 0; //evea(0, 0);
	C[2*(nt-1)+1] = 1; //evea(0, 1);

	//for( int i = 0; i < 5 * 2 * Nm; i ++ ){
	//	host_HH[i] = complex<double>(0.0, 0.0);
	//}
	//host_HH[5*(2*0+0) + 0] = 0; // H[1,-2]
	//host_HH[5*(2*0+0) + 1] = 0; // H[1,-1]
	//host_HH[5*(2*0+1) + 0] = 0; // H[2,-1]
	//for(int i = 0; i < Nm; i ++ ){
	//	host_HH[5*(2*i+0) + 2] += lambda*i*omegaJ - mu*lambda*0.25*i*i;
	//	host_HH[5*(2*i+1) + 2] += - lambda*i*omegaJ + mu*lambda*0.25*i*i; // 10.04.15

 //   host_HH[5*(2*i+0) + 2] += omegaQ * complex<double>(0.5,0) - I * complex<double>(gammaF_0*0.5,0) - I * complex<double>(gammaE_0*0.5,0);
	//	host_HH[5*(2*i+0) + 3] += 0;
	//	host_HH[5*(2*i+1) + 1] += 0;
	//	host_HH[5*(2*i+1) + 2] += - omegaQ * complex<double>(0.5,0) - I * complex<double>(gammaF_0*0.5,0);

	//	host_HH[5*(2*i+0) + 2] += i*omegaJ - mu*i*i - i*gammaO_0*0.5*I;
	//	host_HH[5*(2*i+1) + 2] += i*omegaJ - mu*i*i - i*gammaO_0*0.5*I;

	//	if( i < Nm - 1 ){
	//		host_HH[5*(2*i+0) + 4] += sqrt(double(i+1));
	//		host_HH[5*(2*i+1) + 4] += sqrt(double(i+1));
	//		host_HH[5*(2*i+2) + 0] += sqrt(double(i+1));
	//		host_HH[5*(2*i+3) + 0] += sqrt(double(i+1));
	//	}
	//}

	//host_HH[5*(2*(Nm-1)+0) + 4] = complex<double>(0.0, 0.0); // H[Nm-1,Nm+1]
	//host_HH[5*(2*(Nm-1)+1) + 3] = complex<double>(0.0, 0.0); // H[Nm,Nm+1]
	//host_HH[5*(2*(Nm-1)+1) + 4] = complex<double>(0.0, 0.0); // H[Nm,Nm+2]

	int k = k0;
	for(int j = 0; j < tmax; j++) {
		double t = j*dt;     // real Am = (t <= period) ? Ams : 0.0;
   // printf("tmax = %d\n", tmax);
		double Am = 0;
		double gammaF = gammaF_0;
		double gammaE = gammaE_0;
		double gammaO = gammaO_0;
		if(t < period){
			Am = Ams;
			gammaF = 0;
			gammaE = 0;
			gammaO = 0;
			//lambda = 0.;	
		}
		double f = f0;
	    if(t < 2*period){
			 f = 0;
			 gammaF = 0;
			 gammaE = 0;
			//gammaO = 0;
		}

		if(j < 15000){
			 gammaF = 0;
			 gammaE = 0;
		}

		for( int i = 0; i < 5 * 2 * Nm; i ++ ){
		host_HH[i] = complex<double>(0.0, 0.0);
	}
	host_HH[5*(2*0+0) + 0] = 0; // H[1,-2]
	host_HH[5*(2*0+0) + 1] = 0; // H[1,-1]
	host_HH[5*(2*0+1) + 0] = 0; // H[2,-1]
	for(int i = 0; i < Nm; i ++ ){
		host_HH[5*(2*i+0) + 2] += lambda*i*omegaJ - mu*lambda*0.25*i*i;
		host_HH[5*(2*i+1) + 2] += - lambda*i*omegaJ + mu*lambda*0.25*i*i; // 10.04.15

    host_HH[5*(2*i+0) + 2] += omegaQ * complex<double>(0.5,0) - I * complex<double>(gammaF*0.5,0) - I * complex<double>(gammaE*0.5,0);
		host_HH[5*(2*i+0) + 3] += 0;
		host_HH[5*(2*i+1) + 1] += 0;
		host_HH[5*(2*i+1) + 2] += - omegaQ * complex<double>(0.5,0) - I * complex<double>(gammaF*0.5,0);

		host_HH[5*(2*i+0) + 2] += i*omegaJ - mu*i*i - i*gammaO*0.5*I;
		host_HH[5*(2*i+1) + 2] += i*omegaJ - mu*i*i - i*gammaO*0.5*I;

		if( i < Nm - 1 ){
			host_HH[5*(2*i+0) + 4] += sqrt(double(i+1));
			host_HH[5*(2*i+1) + 4] += sqrt(double(i+1));
			host_HH[5*(2*i+2) + 0] += sqrt(double(i+1));
			host_HH[5*(2*i+3) + 0] += sqrt(double(i+1));
		}
	}

	host_HH[5*(2*(Nm-1)+0) + 4] = complex<double>(0.0, 0.0); // H[Nm-1,Nm+1]
	host_HH[5*(2*(Nm-1)+1) + 3] = complex<double>(0.0, 0.0); // H[Nm,Nm+1]
	host_HH[5*(2*(Nm-1)+1) + 4] = complex<double>(0.0, 0.0); // H[Nm,Nm+2]


		//double Am = (t <= period) ? Ams : 0.0; //Am = 0.0;     
		//double f = (t <= 2*period) ? 0.0 : f0; //f = f0;

		for( int q = 0; q < 30; q ++ ){
			ts = j * dt + q * dt / 30;
			d_RK2(host_HH, Am, f, C, ts, k1, ddt);
			for (int i = 0; i < 2 * Nm; i++)
				tmp[i] = C[i] + complex<double>(0.5,0)*k1[i];
			d_RK2(host_HH, Am, f, tmp, ts + 0.5*ddt, k2, ddt);
			for (int i = 0; i < 2 * Nm; i++)
				tmp[i] = C[i] + complex<double>(0.5,0)*k2[i];
			d_RK2(host_HH, Am, f, tmp, ts + 0.5 * ddt, k3, ddt);
			for (int i = 0; i < 2 * Nm; i++)
				tmp[i] = C[i] + k3[i];
			d_RK2(host_HH, Am, f, tmp, ts + ddt, k4, ddt);
			for (int i = 0; i < 2 * Nm; i++)
				C[i] += complex<double>(1.0/6.0,0) * (k1[i] + complex<double>(2.0,0)*k2[i] +
					complex<double>(2.0,0) * k3[i] + k4[i]);
			ts += ddt;
		}

		for(int i=0; i<Nm; i++)
			a[i] = abs(C[2*i]) * abs(C[2*i]) + abs(C[2*i+1]) * abs(C[2*i+1]);

		double Posc = 0.0;
		double Pq1 = 0.0;
		double Pq2 = 0.0;
		for(int i=0; i<Nm; i++){
			Posc += i*a[i];
			Pq1 += abs(C[2*i]) * abs(C[2*i]);
			Pq2 += abs(C[2*i+1]) * abs(C[2*i+1]);
		}
		
        
		
		//fprintf(fq1out, "%f\t%f\n", t/(2*M_PI), Pq1);
		//fprintf(fq2out, "%f\t%f\n", t/(2*M_PI), Pq2);

		nr_m[j] = nr_m[j] + Posc;
		pq[j] = pq[j] + Pq1;

		if(j == 30999){
			fprintf(fout1, "%d\t%f\n", ni, Posc);
		}

		if(j == 13815){
			fprintf(fout2, "%d\t%f\n", ni, Posc);
		}
		
		if(j == k){
		   fprintf(fout3, "%d\t%f\t%f\n", ni, j*dt/(2*M_PI), Posc);
		   k = k + 1256;
		}
		

    double Pe = 0.0;
    for( int i = 0; i < Nm; i ++ ){
      Pe += abs(conj(C[2*i]) * C[2*i]);
    }

    double Pf = 1.0;

    double P = Posc * gammaO + Pe * gammaE + Pf * gammaF;
    double dP = dt * P;

    double ran1 = d_rdistr(state);
//	double ran1 = (rand()%100)/100; 
	if( dP < ran1 ){
      double d = 0.0;
		  for(int i=0; i < 2*Nm; i++) {
			  d += abs(C[i]) * abs(C[i]);
		  }
		  d = 1.0 / sqrt(d);
		  for (int i = 0; i < 2 * Nm; i++)
			  C[i] = C[i] * d;
    }
    else{
      double ran2 = d_rdistr(state);
	  //double ran2 = (rand()%100)/100;
      //printf("Posc =%f\n", Posc);
	  //printf("Pe =%f\n", Pe);
	 // printf("P =%f\n", P);
     // printf("P1 =%f\n", (Posc * gammaO) / P);
	//  printf("P2 =%f\n", (Posc * gammaO + Pe * gammaE) / P);

      if( ran2 < (Posc * gammaO) / P ){
				for(int i = 0; i < Nm - 1; i++) {
					C[2*i+0] = (sqrt((double)(i+1)))*C[2*i+2];
					C[2*i+1] = (sqrt((double)(i+1)))*C[2*i+3];
				}
				C[2*Nm-2] = 0;
				C[2*Nm-1] = 0;
      }
      else if( ran2 < (Posc * gammaO + Pe * gammaE) / P ){
	  
	  JumpsCount ++;
      printf("JumpsCountE = %3d, t=%f\n", JumpsCount, t);
				for (int i = 0; i < Nm; i++) {
					C[2*i+1] = C[2*i+0];
					C[2*i+0] = 0;
				}

		  //double d = 0.0;
		  //for(int i=0; i < 2*Nm; i++) {
			 // d += abs(C[i]) * abs(C[i]);
		  //}
		  //d = 1.0 / sqrt(d);
		  //for (int i = 0; i < 2 * Nm; i++)
			 // C[i] = C[i] * d;

      }
      else{
	  	  JumpsCount ++;
      printf("JumpsCountF = %3d, t=%f\n", JumpsCount, t);
					for (int i = 0; i < Nm; i++) {
						C[2*i+1] = - C[2*i+1];
					}
      }

      double d = 0.0;
		  for(int i=0; i < 2*Nm; i++) {
			  d += abs(C[i]) * abs(C[i]);
		  }
		  d = 1.0 / sqrt(d);
		  for (int i = 0; i < 2 * Nm; i++)
			  C[i] = C[i] * d;
    }

  }

  //mas[ni] = nr_m[tmax];
  //printf("mas = %f\n", mas[ni]);
    //fprintf(fout, "%d\t%f\n", ni, nr_m[tmax]);

  //printf("nr = %f\n", nr_m[tmax-1]);
  //fprintf(fout1, "%d\t%f\n", ni, nr_m[tmax-1]/ni);
  //fprintf(fout2, "%d\t%f\n", ni, nr_m[13816]/ni);




  }
	
//    double Nstep = (2*M_PI)/dt;
//	int tmax = floor(Nstep)*tperiod;

	for(int i=0; i < tmax; i++) {
		double t = i*dt;
		fprintf(foutnm, "%f\t%f\n", t/(2*M_PI), nr_m[i]/Nrealiz);
		fprintf(foutqb, "%f\t%f\n", t/(2*M_PI), pq[i]/Nrealiz);
	}

	delete[] C;
    delete[] host_HH;

	delete[] nr_m;
	delete[] pq;

   
    //do ii = 0, steps
    //    t = ii*dt
    //    write(12,*) t/(2.0*pi), aver_nr(ii+1)/(Nrealiz)
    //end do



  fclose(fout1);
  fclose(fout2);
  fclose(fout3);
  fclose(foutnm);
  fclose(foutqb);
}
