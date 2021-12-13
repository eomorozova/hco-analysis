
// Calculates a spectrogram  using the Burg (1967) method for estimation of the power spectral density in each time window.
// The Burg method fits the autoregressive (AR) model of a specified order p in the time series by minimizing the sum of squares of the residuals.
// The fast-Fourier transform (FFT) spectrum is estimated using the previously calculated AR coefficients.
// This method is characterized by higher resolution in the frequency domain than traditional FFT spectral analysis,
// especially for a relatively short time window (Buttkus, 2000). 


#include <stdio.h>
#include <stdlib.h>
#include "mex.h"        
#include "matrix.h"      
#include <string.h> 
#include <math.h>
#include <time.h>
#include <iostream>



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
// *fdata		data array;
// full_n		number of samples in the data array *fdata;
// window_dim	number of samples in a window array;
// *Spm_data	spectrogram [n_sp_iter*n_sp];
// k_step		moving window step array
// n_sp			number of samples in a spectrogram
// n_sp_iter	n_sp_it = (full_n - window_dim)/k_step
// rate			sampling frequency [Hz];
// cut_off_fr	low cut_off frequency [Hz] (min frequency in a spectrogram);
// cut_off_fr1	high cut_off frequency [Hz] (max frequency in a spectrogram);
// n_p			number of AR-coeffcient to calculate SPD [n_p < window_dim] --> for example [n_p = window_dim/4 + 1] < 501]

//------------inputs---------------
double *fdata; // data
fdata = mxGetPr(prhs[0]);
int full_n;
full_n = *(mxGetPr(prhs[1])); // number of samples in fdata
int window_dim;
window_dim = *(mxGetPr(prhs[2])); // time window
int k_step; // moving window step
k_step = *(mxGetPr(prhs[3])); 
int n_sp; // number of samples in SPM (fr)
n_sp = *(mxGetPr(prhs[4]));
int n_sp_iter; // number of samples in a spectrogram (time)
n_sp_iter = *(mxGetPr(prhs[5]));
int rate; // sampling rate
rate = *(mxGetPr(prhs[6]));
double cut_off_fr; // low cut off frequency
cut_off_fr = *(mxGetPr(prhs[7]));
double cut_off_fr1; // high cut off frequency
cut_off_fr1 = *(mxGetPr(prhs[8]));
int n_p; // number of AR-coefficients
n_p = *(mxGetPr(prhs[9]));

//--------outputs------------
double *Spm_data; // Spectrogram
plhs[0] = mxCreateDoubleMatrix(1,n_sp_iter*n_sp,mxREAL);
Spm_data = mxGetPr(plhs[0]);

double *f;
plhs[1] = mxCreateDoubleMatrix(1,n_sp,mxREAL);
f = mxGetPr(plhs[1]);




int i, kk, k_sp, n_kor, n_sp_it;
int N_x = window_dim;
double* fmas;
double* spm;
double* Ak;
//double* f;
double* b;
double avr;

double		freq_cut = cut_off_fr;
double		freq_cut1 = cut_off_fr1;
		n_kor = n_p;
		if (n_kor > window_dim) n_kor = window_dim;
		if (n_kor > 501) n_kor = 501;
//-------------------------------------
		fmas = new double[window_dim];
		spm = new double[n_sp];	
		Ak = new double[n_kor + 1];
		f = new double[N_x];
		b = new double[N_x];
//-------------------------------------
		n_sp_it = (full_n - window_dim)/k_step;
		if(n_sp_it > n_sp_iter) n_sp_it = n_sp_iter;

// Overlapping step
	for(k_sp = 0; k_sp < n_sp_it; k_sp++)
	{
		for(i = 0; i < window_dim; i++) 
		{
			kk = k_sp*k_step + i; 
			if(kk < full_n) fmas[i] = fdata[kk];
			else fmas[i] = 0;
		}
		avr = 0;
		for(i = 0; i < window_dim; i++) avr += fmas[i];
		avr /= window_dim;
		for(i = 0; i < window_dim; i++) fmas[i] = fmas[i] - avr;

// =========== Calculating autoregression coefficients using Bergs method ===============
//			BurgAlgorithmSPD(fmas, window_dim, spm, n_sp, freq_cut, n_kor);
// GET SIZE FROM INPUT VECTORS 
int j, k, k1, n;
double mu, t1, t2;
int N = N_x - 1; 
	//size_t m = coeffs.size(); 
// INITIALIZE Ak 
	for(i=1; i<n_kor+1; i++) Ak[i] = 0.0;
	Ak[0] = 1.0; 

// INITIALIZE f and b ( x ); 

	for(i=0; i<N_x; i++)
		{
			f[i] = fmas[i];
			b[i] = fmas[i];
		}
// INITIALIZE Dk 
double Dk = 0.0; 
	for(j=0; j<=N; j++ ) Dk += 2.0*f[j]*f[j]; 
    Dk -= f[0]*f[0] + b[N]*b[N]; 
// BURG RECURSION 
	for(k=0; k < n_kor; k++) 
	{ 
		k1 = k + 1;
// COMPUTE MU 
		mu = 0.0; 
		for(n=0; n < N-k; n++) mu += f[n + k1]*b[n]; 
		mu *= -2.0/Dk; 
// UPDATE Ak 
		for(n=0; n <= k1/2; n++) 
		{ 
			t1 = Ak[n] + mu*Ak[k1 - n]; 
			t2 = Ak[k1 - n] + mu*Ak[n]; 
			Ak[n] = t1; 
			Ak[k1 - n] = t2; 
		} 
// UPDATE f and b 
		for (n=0; n < N-k; n++) 
		{ 
			t1 = f[n + k1] + mu*b[n]; 
			t2 = b[n] + mu*f[n + k1]; 
			f[n + k1] = t1; 
			b[n] = t2; 
		} 
// UPDATE Dk 
		Dk = (1.0 - mu*mu)*Dk - f[k1]*f[k1] - b[N - k - 1]*b[N - k - 1]; 
	} 
// ASSIGN COEFFICIENTS coeffs.assign( ++Ak.begin(), Ak.end() ); 
//	for(i=0; i<n_p; i++) coeffs[i] = Ak[i+1];

//============ SPD =====================
double sum, tmp, arg, pi2, df;
//	for(i=0; i<n_p; i++) Ak[i+1] *= (double)(2*n_p-i)/(2*n_p);
	pi2=2*freq_cut*(4.0*atan(1.0))/rate/n_sp;

// 	for(j=0;j<n_sp;j++)
// 	{
// 		sum = 1.05;
// 		tmp = 0;
// 		arg = pi2*j;
// 		for(i=0; i < n_kor; i++)
// 		{
// 			k = i + 1;
// 			sum += Ak[k]*cos(arg*k);
// 			tmp += Ak[k]*sin(arg*k);
// 		}
// //		spm[j]=1./sqrt(sum*sum + tmp*tmp);
// 		Spm_data[k_sp*n_sp + j] = 1./sqrt(sum*sum + tmp*tmp);
// 	}

double pi=8.0*atan(1.0)/rate;       
double f1=freq_cut*pi; 
       df=pi*(freq_cut1-freq_cut)/(n_sp-1);
    	for(j=0;j<n_sp;j++)
	{
		sum = 1.05;
		tmp = 0;
		arg = f1+j*df;
		for(i=0; i < n_kor; i++)
		{
			k = i + 1;
			sum += Ak[k]*cos(arg*k);
			tmp += Ak[k]*sin(arg*k);
		}
//		spm[j]=1./sqrt(sum*sum + tmp*tmp);
		Spm_data[k_sp*n_sp + j] = 1./sqrt(sum*sum + tmp*tmp);
	}

//========================================================================================
//	for(i = 0; i < n_sp; i++) Spm_data[k_sp*n_sp + i] = spm[i];
		}

	delete[] fmas;
	delete[] spm;
	delete[] Ak;
	//delete[] f;
	delete[] b;
return;
}