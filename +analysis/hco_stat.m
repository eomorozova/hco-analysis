
function [hcostat]=hco_stat(V, Fs)

% Inputs:
% V     membrane potential recording (mV)
% Fs    sampling frequency (Hz)
%
% Outputs:
% hcostat
%   .T _mean    mean cycle periods based on the spikes
%   .T_std standard deviation of HCO cycle periods  based on the spikes
%   .T HCO cycle periods based on spikes
%   .T1 _mean mean HCO cycle periods based on the slow wave
%   .T1_std standard deviation of HCO cycle periods  based on the slow wave
%   .T1 HCO cycle periods based on the slow wave
%   .nSpks_mean mean number of spikes per burst 
%   .nSpks_std standard deviation of the number of spikes per burst based
%   .SpkFreq_mean mean frequency of spikes withn the burst
%   .SpksFreq_std standard deviation of the frequency of spikes within the
%   .dc duty cycles based on the spikes 
%   .dc_mean mean duty cycle based on the spikes
%   .dc_std standard deviation of the duty cycle based on the spikes
%   .DC duty cycles based on the slow wave 
%   .DC_mean mean duty cycle based on the slow wave
%   .DC_std standard deviation of the duty cycle based on the slow wave
%   .A amplitude of the slow wave
%   .mean_pks  mean slow wave peak
%   .mean_dip mean slow wave trough
%   .st spike times

if nargin<2
    Fs=10000; % sampling frequency (Hz)
end
%% calculates burst period of HCO, number of spikes per burst and spike frequency
%i=find(gh1==350); j=find(gsyn1==450);
%V=Vhco4{i,j}(2,1:end);
%V=Vhco{1}{45,75}(1:end);
%V=maps.Vhco{jj}{1}{75,90}(1:end);
%V=Vseg{1}{9};
%V=Vhco{1}{4}{31}(1,5*Fs(1):end);
 
%[st_v,st]=findpeaks(smooth(V,10),'MinPeakProminence',2,'MinPeakHeight',-38);
[st_v,st]=findpeaks(smooth(V,10),'MinPeakProminence',1.5,'MinPeakHeight',-38);

%burststat=detect_bursts(st'/Fs,mean(diff(st))/Fs,mean(diff(st))/Fs+1.2,1); % for most experiments
%burststat=detect_bursts(st'/Fs,mean(diff(st))/Fs,mean(diff(st))/Fs+0.75,1); % for 936_122
%burststat=detect_bursts(st'/Fs,mean(diff(st))/Fs,mean(diff(st))/Fs+0.5,1);
%burststat=detect_bursts(st'/Fs,mean(diff(st))/Fs,mean(diff(st))/Fs+0.4,1);
burststat=analysis.detect_bursts(st'/Fs,mean(diff(st))/Fs,mean(diff(st))/Fs+0.5,1); % 951_060
%burststat=detect_bursts(st'/Fs,mean(diff(st))/Fs,mean(diff(st))/Fs+1.5,1); % for 936_074
%burststat=detect_bursts(st'/Fs,0.9,0.9,1);
%burststat=detect_bursts(st'/Fs,0.55,0.55,1); % 936_102 release
T=diff(burststat.center);
%if isempty(st)==1
%    st=0;
%end

if length(burststat.center)>2
    T_mean=mean(diff(burststat.center)); T_std=std(diff(burststat.center));
    if length(burststat.center<3)
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
    T_mean=0; nSpks_mean=0; SpkFreq_mean=0; SpkFreq_median=0; dc_mean=0; dc=0;
    T_std=0; nSpks_std=0; SpkFreq_std=0; dc_std=0;
end
%uncomment and run this subsection to see how well bursts are identified
%figure(2),clf, plot(1/Fs:1/Fs:length(V)/Fs,V),
%hold on, plot(burststat.center, repmat(-45, 1, length(burststat.center)),'.','markersize',20)
%hold on, plot(st/Fs,st_v,'.','markersize',15,'color','k')

%% filter data, calculate the amplitude of the slow wave and the frequency by threholding (1st neuron)
%V=maps.Vhco{jj}{1}{75,90}(1:end);
%V=Vhco{1}{1}{12}(2,:);
[b,a]=butter(4,1/(Fs/2));
Vfilt = filter(b,a,V); % filter the data to frequency ~2 times higher than desired frequency
Vfilt_smooth=smooth(Vfilt(1:1:end),1000);
%Vfilt_smooth=smooth(Vfilt(1:1:end),100); % if Fs=1000 Hz
%[pks,lpks]=findpeaks(Vfilt_smooth(Fs:end),'MinPeakDistance',Fs/10*2,'MinPeakProminence',3*std(Vfilt_smooth(Fs:end)));
%[dips,ldips]=findpeaks(-Vfilt_smooth(Fs:end),'MinPeakDistance',Fs/10*2,'MinPeakProminence',3*std(Vfilt_smooth(Fs:end)));

%deltaV=3; % for most experiments
%deltaV=1.8; % for 936_122 (release)
deltaV=2.5; % for 936_074
[pks,lpks]=findpeaks(Vfilt_smooth(Fs:end),'MinPeakProminence',3,'MinPeakHeight',mean(Vfilt_smooth(Fs:end))+deltaV);
[dips,ldips]=findpeaks(-Vfilt_smooth(Fs:end),'MinPeakProminence',3,'MinPeakHeight',mean(-Vfilt_smooth(Fs:end))+deltaV);
if length(pks)>=2
%if length(pks)>2 % 936_074
mean_pks=median(pks); mean_dip=median(-dips);
A=mean_pks - mean_dip; % slow wave amplitude
Vth=(mean_pks + mean_dip)/2; 
st_filt = analysis.crossings(smooth(Vfilt(1:1:end),1000)',Vth); st_filt=st_filt{1}./Fs; % in sec
[~,st_filtdown] = analysis.crossings(smooth(Vfilt(1:1:end),1000)',Vth); st_filt_down=st_filtdown{1}./Fs; % in sec

if length(st_filt)>=2; T1=diff(st_filt); else; T1=0; end
T1_mean=mean(T1); T1_std=std(T1);

% based on the peaks of the slow-wave
%if length(lpks)>=2; T1=diff(lpks)./Fs; else; T1=0; end
%T1_mean=mean(T1); T1_std=std(T1);

if length(st_filt)>=2
%if length(st_filt)>2 % 936_074
if st_filt_down(1)<st_filt(1)
    if length(st_filt_down(2:end))<length(st_filt)
        BurDur=st_filt_down(2:end)-st_filt(1:end-1); % burst duration based on slow wave
    else
        BurDur=st_filt_down(2:end)-st_filt; % burst duration based on slow wave
    end
else
    if length(st_filt)> length(st_filt_down)
        BurDur=st_filt_down-st_filt(1:end-1); % burst duration based on slow wave
    else
        BurDur=st_filt_down-st_filt; % burst duration based on slow wave
    end
end
end

if length(st_filt)>=2
%if length(st_filt)>2 % 936_074
    if length(BurDur)>length(T1)
        DC=BurDur(2:end-1)./T1(2:end); % duty cycle based on the slow wave
    else
        DC=BurDur(2:end)./T1(2:end); % duty cycle based on the slow wave
    end
    DC_mean = mean(DC); DC_std = std(DC);
else
    DC=0; DC_mean=0;  DC_std=0;
end

else
mean_pks=0; mean_dip=0; A=0; T1=0; T1_mean=0; T1_std=0; DC=0; DC_mean=0; DC_std=0;
end

%uncomment to check period calculation
% figure(2),clf, plot(V,'linewidth',1)
% hold on, plot(Vfilt_smooth(Fs:end),'linewidth',1)
% hold on, plot(lpks,pks,'.','markersize',20)
% hold on, plot(ldips,-dips,'.','markersize',20)
% hold on, mmyplothorzline(mean(Vfilt_smooth(Fs:end))+deltaV)
% hold on, mmyplothorzline(-mean(-Vfilt_smooth(Fs:end))-deltaV)
% hold on, mmyplothorzline(Vth);
% ylim([-60 -25])
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
hcostat.SpkFreq_std=SpkFreq_std;
hcostat.dc=dc;
hcostat.dc_mean=dc_mean;
hcostat.dc_std=dc_std;
hcostat.DC=DC;
hcostat.DC_mean=DC_mean;
hcostat.DC_std=DC_std;
hcostat.A=A;
hcostat.mean_pks=mean_pks;
hcostat.mean_dip=mean_dip;
hcostat.st=st;
hcostat.burststat=burststat;

%% frequency based on the spectrogram (uncomment)
%i_gsyn = 75; i_gh = 75;
%V=Vhco{2}{i_gsyn,i_gh}(5*Fs:end);
%V1=Vhco{1}{i_gsyn,i_gh}(5*Fs:end);

%if length(V)/Fs<window
%SPM=NaN; f1=0; fmaxS=0; T11_mean=0; T11_std=0;
%else
%[SPM1,f1,fmaxS]=spectrogram(V,Fs,Fs1,window); % calculate the spectrogram
%T11=1./fmaxS; T11_mean=mean(T11); T11_std=std(T11);
%end

% uncomment to plot spectrogram
% clf, subplot(3,1,1)
% plot(V), hold on, plot(V1)
% subplot(3,1,2)
% Fcutoff=0; Fcutoff1=1; n_sp=256;
% f=linspace(Fcutoff,Fcutoff1,n_sp);
% x=window/2/Fs1:window/2/Fs1:size(SPM1,2)*window/2/Fs1;
% pcolor(x,f,SPM1),shading flat, colorbar
% %caxis([min(prctile(SPM1(1:5:end,:),90)) max(prctile(SPM1(1:5:end,:),90))])
% ylabel('Frequency, Hz'); xlabel('Time, sec'); 
% title('Spectrogram of voltage trace'); set(gca, 'Fontsize',14)
% grid minor; grid on; set(gca, 'layer','top')
% ax = gca; ax.GridLineStyle = '-'; ax.GridAlpha = 1 ;
% % plot maximum frequency from the spectrogram
% subplot(3,1,3), plot(fmaxS,'.','markersize',20)
% ylim([0 Fcutoff1])

