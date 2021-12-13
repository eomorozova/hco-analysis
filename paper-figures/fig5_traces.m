% load the data to reproduce figure 5A,B

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig5AB.mat')

% reproduce figure 5A,B (release and escape traces and spectrograms)

clf
g=0.02; l=0.06;

% release at 10C
x=1/Fs:1/Fs:1*length(data.V{1}(1,1:1+30*Fs))/data.Fs;
for i=1:2
    subtightplot(12,4,1+(i-1)*4,g,l,l), plot(x,smooth(data.V{1}(i,1:1+30*Fs),10),'k')
    xlim([0 x(end)]); ylim([-60 -18]);
    xticks([0:5:30]); yticks([-60:20:-20])
    if i==1; set(gca,'xtick',[]); title(['10' char(176) 'C']); end
    %if i==2; xlabel('Time, sec'); end
    set(gca, 'Fontsize',12,'FontName','Arial'); box off;
    h=display.plothorzline(data.Vth(1)); set(h,'color','r','linewidth',1)
end

% release at 20C
for i=1:2
    subtightplot(12,4,2+(i-1)*4,g,l,l), plot(x,smooth(data.V{1}(i,3*60*data.Fs:3*60*data.Fs+30*Fs),10),'k')
    xlim([0 x(end)]); ylim([-60 -20]);
    xticks([0:5:30]); yticks([-60:20:-20])
    if i==1; set(gca,'xtick',[]); title(['20' char(176) 'C']); end
    %if i==2; xlabel('Time, sec'); end
    set(gca, 'Fontsize',12,'FontName','Arial'); box off;
    h=display.plothorzline(data.Vth(1)); set(h,'color','r','linewidth',1)
end

% escape at 10C
x=1/Fs:1/Fs:1*length(data.V{2}(1,0.2*60*data.Fs:0.2*60*data.Fs+30*Fs))/data.Fs;
for i=1:2
    subtightplot(12,4,3+(i-1)*4,g,l,l), plot(x,smooth(data.V{2}(i,0.2*60*data.Fs:0.2*60*data.Fs+30*Fs),10),'k')
    xlim([0 x(end)]); ylim([-60 -20]);
    xticks([0:5:30]); yticks([-60:20:-20])
    if i==1; set(gca,'xtick',[]); end
    if i==1; set(gca,'xtick',[]); title(['10' char(176) 'C']); end
    %if i==2; xlabel('Time, sec'); end
    set(gca, 'Fontsize',12,'FontName','Arial'); box off;
    h=display.plothorzline(data.Vth(2)); set(h,'color','r','linewidth',1)
end

% escape at 20C
x=1/Fs:1/Fs:1*length(data.V{2}(1,3.8*60*data.Fs:3.8*60*data.Fs+30*Fs))/data.Fs;
for i=1:2
    subtightplot(12,4,4+(i-1)*4,g,l,l), plot(x,smooth(data.V{2}(i,3.8*60*data.Fs:3.8*60*data.Fs+30*Fs),10),'k')
    xlim([0 x(end)]); ylim([-60 -20]);
    xticks([0:5:30]); yticks([-60:20:-20])
    if i==1; set(gca,'xtick',[]); end
    if i==1; set(gca,'xtick',[]); title(['20' char(176) 'C']); end
    %if i==2; xlabel('Time, sec'); end
    set(gca, 'Fontsize',12,'FontName','Arial'); box off;
    h=display.plothorzline(data.Vth(2)); set(h,'color','r','linewidth',1)
end

g=0.05; l=0.07;
for j=1:2
    x=1/data.Fs/60:1/data.Fs/60:1*length(data.V{j})/data.Fs/60;
    % whole trace
    for i =1:2
        subtightplot(6,2,3+(i-1)*2+j-1,g,l,l), plot(x,data.V{j}(i,:),'k')
        ylabel('V_M, mV');
        if i==1 && j==1; title('Release'); end
        if i==1 && j==2; title('Escape'); end
        set(gca, 'Fontsize',14,'FontName','Arial'); box off
        ylim([-60 -18]); xlim([0 x(end)])
        set(gca,'xtick',[])
        h=display.plothorzline(data.Vth(j)); set(h,'color','r','linewidth',1)
    end
    
    % saline temperature
    subtightplot(6,2,7+j-1,g,l,l), plot(x,data.T{j},'k','linewidth',2)
    ylabel('Temp, ^oC');
    set(gca,'xtick',[])
    set(gca, 'Fontsize',14,'FontName','Arial'); box off;
    xlim([0 x(end)]); ylim([9 21])
    
    % ISIs
    sthco=analysis.crossings(diff(smooth(data.V{j}(1,:),10))',0.4);
    ISIhco=diff(sthco{1})/data.Fs;
    
    subtightplot(6,2,9+j-1,g,l,l), plot(sthco{1}(1:end-1)/Fs/60,ISIhco,'.','markersize',5,'color','k')
    xlim([0 sthco{1}(end-1)/Fs/60])
    set(gca, 'YScale', 'log'); box off
    yticks([10^-2,10^-1,1,10]); ylim([10^-2 10]); xticks([0:1:10])
    ylabel('log10(ISI), sec');
    set(gca,'xtick',[])
    set(gca, 'Fontsize',14,'FontName','Arial')
    
    % Spectrogram
    [SP,f,time] = analysis.ARspectrogram(data.V{j}(1,:),data.Fs);
    
    subtightplot(6,2,11+j-1,g,l,l)
    colormap(parula)
    pcolor(time,f,log10(SP)); shading flat
    caxis([min(prctile(log10(SP),90)) max(prctile(log10(SP),90))])
    ylabel('Frequency, Hz'); yticks([0.0:0.2:0.6])
    set(gca, 'Fontsize',14,'FontName','Arial')
    ylim([0 0.65]); xlim([0 x(end)])
    xlabel('Time, min')
    
end

