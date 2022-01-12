%% load the data for supplemetary figure S2

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_FigS2')

%% calculate ERQ for each neuron in the circuit

for n=1:2
    ERQ(n,:)=(movmean(data.V(n,1:10:end),data.Fs*60)-data.Vth)./movmean(data.V(n,1:10:end),data.Fs*60);
end

%% Reproduce figure 7-figure supplement
clf

x=1/data.Fs:1/data.Fs:1*length(data.V(1,1*data.Fs*60:2*data.Fs*60))/data.Fs;
for n=1:2
    subplot(8,2,n*2-1), plot(x,data.V(n,1*data.Fs*60:2*data.Fs*60),'k','linewidth',1)
    h=display.plothorzline(data.Vth); set(h,'color','r','linewidth',1)
    ylabel('V_M, mV'); if n==1; title('10^oC'); end; if n==2; xlabel('Time, sec'); end
    ylim([-52 -12]); xlim([0 x(end)]);
    if n==1; set(gca,'xtick',[]); end; box off
    set(gca, 'Fontsize',12,'Fontname','Arial')
end

x=1/data.Fs:1/data.Fs:1*length(data.V(1,end-3*data.Fs*60:end-2*data.Fs*60))/data.Fs;
for n=1:2
    subplot(8,2,n*2), plot(x,data.V(n,end-3*data.Fs*60:end-2*data.Fs*60),'k','linewidth',1)
    h=display.plothorzline(data.Vth); set(h,'color','r','linewidth',1)
    ylabel('V_M, mV'); if n==1; title('20^oC'); end; if n==2; xlabel('Time, sec'); end
    ylim([-52 -12]); xlim([0 x(end)]);
    if n==1; set(gca,'xtick',[]); end; box off
    set(gca, 'Fontsize',12,'Fontname','Arial')
end

x=1/data.Fs/60:1/data.Fs/60:1*length(data.V)/data.Fs/60;
for n=1:2
    subplot(8,1,2+n), plot(x,data.V(n,:),'k','linewidth',1)
    h=display.plothorzline(-38); set(h,'color','r','linewidth',2)
    ylabel('V_M, mV');
    set(gca, 'Fontsize',12,'Fontname','Arial')
    ylim([-52 -12]); xlim([0 x(end)]); set(gca,'xtick',[]); box off
end

subplot(8,1,5) % Temperature
plot(x,data.T,'k','linewidth',3)
ylabel('Temp, ^oC');
xlim([0 x(end)]); ylim([9 21]); set(gca,'xtick',[]); box off
set(gca, 'Fontsize',12,'Fontname','Arial')

subplot(8,1,6) % ERQ
for n=1:2
    plot(x(1:10:end),ERQ(n,:),'linewidth',2), hold on
end
ylabel('ERQ');
set(gca, 'Fontsize',12,'Fontname','Arial')
ylim([-0.05 0.1]);xlim([0 x(end)]); set(gca,'xtick',[]); box off

subplot(8,1,7) % ISI
sthco=analysis.crossings(diff(smooth(data.V(1,:),10))',0.15);
ISIhco=diff(sthco{1})/data.Fs;
plot(sthco{1}(1:end-1)/data.Fs/60,ISIhco,'.','markersize',7,'color','k')
xlim([0 sthco{1}(end-1)/data.Fs/60])
set(gca, 'YScale', 'log'); box off
yticks([10^-2,10^-1,1,10]); ylim([10^-1 10]); set(gca,'xtick',[])
ylabel('log10(ISI), sec'); 
set(gca, 'Fontsize',12,'FontName','Arial')

subplot(8,1,8) % spectrogram
n_sp=256; window = n_sp/6; % time window 
[SP,f,time] = analysis.ARspectrogram(data.V(1,:),data.Fs,window);
colormap(parula)
pcolor(time,f,log10(SP)); shading flat
caxis([min(prctile(log10(SP),95)) max(prctile(log10(SP),95))])
ylabel('Frequency, Hz'); xlabel('Time, min')
set(gca, 'Fontsize',12,'FontName','Arial')
ylim([0 0.25]);
