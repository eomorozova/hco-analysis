
function [hcostat]=hco_stat1(V, Fs, V_th)

% Inputs:
% V     membrane potential recording (mV)
% Fs    sampling frequency (Hz)
% V_th  synaotic threhold
%
% Outputs:
% hcostat
%   .T _mean        mean cycle period (based on spikes)
%   .T_std          standard deviation of the cycle period (based on spikes)
%   .T              cycle periods (based on spikes)
%   .T1 _mean       mean cycle period (based on slow-wave)
%   .T1_std         standard deviation of the cycle period (based on slow-wave)
%   .T1             cycle periods (based on slow-wave)
%   .nSpks_mean     mean number of spikes per burst 
%   .nSpks_std      standard deviation of the number of spikes per burst
%   .SpkFreq_mean   mean frequency of spikes within the burst
%   .SpkFreq_median  mean frequency of spikes within the burst
%   .SpksFreq_std   standard deviation of the frequency of spikes within the
%   .dc             duty cycles (based on spikes) 
%   .dc_mean        mean duty cycle (based on spikes)
%   .dc_std         standard deviation of the duty cycle (based on spikes)
%   .DC             duty cycles (based on slow-wave) 
%   .DC_mean        mean duty cycle (based on slow-wave)
%   .DC_std         standard deviation of the duty cycle (based on slow-wave)
%   .DC1            duty cycle (based on the time above synaptic threshold)
%   .DC1_mean       mean duty cycle (based on the time above synaptic threshold)
%   .A              amplitude of the slow-wave
%   .mean_pks       mean values of the slow-wave peaks
%   .mean_dip       mean values of the slow-wave troughs
%   .BurDur         burst durations based on slow-wave
%   .st             spike times

if nargin<2
    Fs=10000; % sampling frequency (Hz)
end
%% calculates circuit output characteristics based on spikes

[st_v,st]=findpeaks(smooth(V,10),'MinPeakProminence',1.5,'MinPeakHeight',-38,'MinPeakDistance',10);

burststat=analysis.detect_bursts(st'/Fs,mean(diff(st))/Fs,mean(diff(st))/Fs+0.3,1);

T=diff(burststat.center);

if length(burststat.center)>2
    
    T_mean=mean(diff(burststat.center)); T_std=std(diff(burststat.center));
    if length(burststat.center)<3
        nSpks_mean=mean(burststat.nSp(1:end)); nSpks_std=std(burststat.nSp(1:end));
        SpkFreq_mean=mean(burststat.SpFreq(1:end)); SpkFreq_std=std(burststat.SpFreq(1:end));
        SpkFreq_median=median(burststat.SpFreq(1:end));
        dc = (burststat.BuDur(1:end-1))./T(1:end)*100;
    else
        nSpks_mean=mean(burststat.nSp(2:end-1)); nSpks_std=std(burststat.nSp(2:end-1));
        SpkFreq_mean=mean(burststat.SpFreq(2:end-1)); SpkFreq_std=std(burststat.SpFreq(2:end-1));
        SpkFreq_median=median(burststat.SpFreq(2:end-1));
        dc = (burststat.BuDur(2:end-1))./T(2:end)*100;
    end
    dc_mean = mean(dc); dc_std = std(dc);
else
    T_mean=NaN; nSpks_mean=NaN; SpkFreq_mean=NaN; SpkFreq_median=NaN; dc_mean=NaN; dc=NaN;
    T_std=NaN; nSpks_std=NaN; SpkFreq_std=NaN; dc_std=NaN;
end

%% calculate circuit output characteristics based on the slow-wave oscillations

[b,a]=butter(4,1/(Fs/2));
Vfilt = filter(b,a,V); % filter the data to frequency ~2 times higher than desired frequency

Vfilt_smooth=smooth(Vfilt(1:end),1000);
deltaV=2.5;

pks = findpeaks(Vfilt_smooth(Fs:end),'MinPeakProminence',3,'MinPeakHeight',mean(Vfilt_smooth(Fs:end))+deltaV);
dips = findpeaks(-Vfilt_smooth(Fs:end),'MinPeakProminence',3,'MinPeakHeight',mean(-Vfilt_smooth(Fs:end))+deltaV);

if numel(pks)>2
    mean_pks=median(pks); mean_dip=median(-dips);
    A=mean_pks - mean_dip; % slow-wave amplitude
    Vth=(mean_pks + mean_dip)/2;
           
    st_filt = analysis.crossings(Vfilt_smooth',Vth); st_filt=st_filt{1}./Fs; % in sec
    [~,st_filtdown] = analysis.crossings(Vfilt_smooth',Vth); st_filt_down=st_filtdown{1}./Fs; % in sec
    
    % crossings with a synaptic threhold
    st_filt1 = analysis.crossings(Vfilt_smooth',V_th); st_filt1=st_filt1{1}./Fs; % in sec
    [~,st_filtdown1] = analysis.crossings(Vfilt_smooth',V_th); st_filt_down1=st_filtdown1{1}./Fs; % in sec    
    
     if numel(st_filt)>2 && numel(st_filt1)>2 && ~isempty(burststat.nSp) 
        T1=diff(st_filt); T_1=diff(st_filt1);
        
        if st_filt_down(1)<st_filt(1)
            if length(st_filt_down(2:end))<length(st_filt)
                BurDur=st_filt_down(2:end)-st_filt(1:end-1); % burst duration based on the slow-wave  
            else
                BurDur=st_filt_down(2:end)-st_filt; 
            end
        else
            if length(st_filt)> length(st_filt_down)
                BurDur=st_filt_down-st_filt(1:end-1);
            else
                BurDur=st_filt_down-st_filt;            
            end
        end
        
        if st_filt_down1(1)<st_filt1(1)
            if length(st_filt_down1(2:end))<length(st_filt1)
                BurDur1=st_filt_down1(2:end)-st_filt1(1:end-1); % burst duration based on the time above the synaptic threhold
            else
                BurDur1=st_filt_down1(2:end)-st_filt1;
            end
        else
            if length(st_filt1)> length(st_filt_down1)
                BurDur1=st_filt_down1-st_filt1(1:end-1);
            else
                BurDur1=st_filt_down1-st_filt1;
            end
        end
        
        if length(BurDur)>length(T1)
            DC=BurDur(2:end-1)./T1(2:end)*100; % duty cycle based on the slow wave  
        else
            DC=BurDur(2:end)./T1(2:end)*100;             
        end
        
        if length(BurDur1)>length(T_1)
            DC1=BurDur1(2:end-1)./T_1(2:end)*100; % duty cycle based on the time above synaptic threhold
        else
            DC1=BurDur1(2:end)./T_1(2:end)*100; 
        end
        
        T1_mean=mean(T1); T1_std=std(T1);
        DC_mean = mean(DC); DC_std = std(DC); DC1_mean = mean(DC1);
    else
        T1=NaN; T1_mean=NaN; T1_std=NaN;
        DC=NaN; DC_mean=NaN; DC_std=NaN;
        BurDur=NaN; DC1=NaN; DC1_mean=NaN;
    end
    
else
    mean_pks=NaN; mean_dip=NaN; A=NaN; T1=NaN; T1_mean=NaN; T1_std=NaN; DC=NaN; DC_mean=NaN; DC_std=NaN; BurDur=NaN; DC1=NaN; DC1_mean=NaN;
end

%%

hcostat.T_mean=T_mean;
hcostat.T_std=T_std;
hcostat.T=T;
hcostat.T1_mean=T1_mean;
hcostat.T1_std=T1_std;
hcostat.T1=T1;
hcostat.nSpks_mean=nSpks_mean;
hcostat.nSpks_std=nSpks_std;
hcostat.SpkFreq_mean=SpkFreq_mean;
hcostat.SpkFreq_median=SpkFreq_median;
hcostat.SpkFreq_std=SpkFreq_std;
hcostat.dc=dc;
hcostat.dc_mean=dc_mean;
hcostat.dc_std=dc_std;
hcostat.DC=DC;
hcostat.DC_mean=DC_mean;
hcostat.DC_std=DC_std;
hcostat.DC1=DC1;
hcostat.DC1_mean=DC1_mean;
hcostat.A=A;
hcostat.mean_pks=mean_pks;
hcostat.mean_dip=mean_dip;
hcostat.BurDur=BurDur;
hcostat.st=st;
hcostat.burststat=burststat;
