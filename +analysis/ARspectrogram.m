
function [SP,f,time] = ARspectrogram(V,Fs,window)

Fs1=2;
[b,a]=butter(4,Fs1/(Fs/2));
Vfilt = filter(b,a,V); % filter the data to frequency ~2 times highr than desired frequency 
step1 = Fs/Fs1; % to make new sampling frequency 10Hz
datafilt = Vfilt(1:step1:end); % downsample
n_sp=256; % frequency resolution
%window = n_sp/10; % time window in samples window/Fs1=256/10=25.6 s
n_sp_iter = (length(datafilt)-window)/(window/2);
Fcutoff = 0; % low cutoff frequency
Fcutoff1 = 2; % high cutoff frequency
n_p = window/4+1; % number of AR coefficients to fit the data

%mex AR_Burg_SpGrm_Mex.cpp % uncoment when ruiing got the first time to compile MEX file

SPM = AR_Burg_SpGrm_Mex(datafilt',length(datafilt),window,window/2,n_sp,n_sp_iter,Fs1,Fcutoff,Fcutoff1,n_p);


SP=reshape(SPM,n_sp,[]); 

f=linspace(0,Fcutoff1,n_sp); % frequency
time=window/2/Fs1/60:window/2/Fs1/60:size(SP,2)*window/2/Fs1/60; % time