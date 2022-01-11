%% load data from the temperature experiments with Q10=1

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig5_suppl.mat')

%% calculate Coefficient of variation for different characteristics of circuit output at 10 and 20C

% release

for i = 1:length(Data.r.V)
    for j=1:2
        [hcostat{i,j}] = analysis.hco_stat(Data.r.V{i}{j}(1,:), Data.r.Fs(i));
         ERQ(i,j) = analysis.erq(Data.r.V{i}{j}(1,:),Data.r.Vth(i));
        
        if (i==1 && j==2) || (j==9 && j==1)  % contains incomplete first burst
            
            meanFR(i,j) = mean(hcostat{i,j}.burststat.BuFreq(2:end));
            stdFR(i,j) = std(hcostat{i,j}.burststat.BuFreq(2:end));
            
            meanSpkFR(i,j) = mean(hcostat{i,j}.burststat.SpFreq(2:end));
            stdSpkFR(i,j) = std(hcostat{i,j}.burststat.SpFreq(2:end));
            
        elseif (i==2 && j==1) || (i==6 && j==1) || i==8 || (j==9 && j==2) % contains incomplete last burst
            
            meanFR(i,j) = mean(hcostat{i,j}.burststat.BuFreq(1:end-1));
            stdFR(i,j) = std(hcostat{i,j}.burststat.BuFreq(1:end-1));
            
            meanSpkFR(i,j) = mean(hcostat{i,j}.burststat.SpFreq(1:end-1));
            stdSpkFR(i,j) = std(hcostat{i,j}.burststat.SpFreq(1:end-1));
            
        elseif (i==2 && j==2) || (i==6 && j==2) || (i==7) || i==10 % contains incomplete first&last burst
            
            meanFR(i,j) = mean(hcostat{i,j}.burststat.BuFreq(2:end-1));
            stdFR(i,j) = std(hcostat{i,j}.burststat.BuFreq(2:end-1));
            
            meanSpkFR(i,j) = mean(hcostat{i,j}.burststat.SpFreq(2:end-1));
            stdSpkFR(i,j) = std(hcostat{i,j}.burststat.SpFreq(2:end-1));
            
        else % all the bursts are complete
            
            meanFR(i,j) = mean(hcostat{i,j}.burststat.BuFreq);
            stdFR(i,j) = std(hcostat{i,j}.burststat.BuFreq);
            
            meanSpkFR(i,j) = mean(hcostat{i,j}.burststat.SpFreq);
            stdSpkFR(i,j) = std(hcostat{i,j}.burststat.SpFreq);
                    
        end
        
        CVFR(i,j) = stdFR(i,j)/meanFR(i,j);
        CVSpkFR(i,j) = stdSpkFR(i,j)/meanSpkFR(i,j);
        
    end
end

outr = [{CVFR},{CVSpkFR}];

%% calculate CV of different metrics of circuit output at 10 and 20C

% escape

for i = 1:length(Data.e.V)
    for j=1:2
        [hcostat{i,j}] = analysis.hco_stat(Data.e.V{i}{j}(1,:), Data.e.Fs(i));
        ERQ(i,j) = analysis.erq(Data.e.V{i}{j}(1,:),Data.e.Vth(i));
        
        if (j==5 && j==2)  % contains incomplete first burst
            
            meanFR(i,j) = mean(hcostat{i,j}.burststat.BuFreq(2:end));
            stdFR(i,j) = std(hcostat{i,j}.burststat.BuFreq(2:end));
            
            meanSpkFR(i,j) = mean(hcostat{i,j}.burststat.SpFreq(2:end));
            stdSpkFR(i,j) = std(hcostat{i,j}.burststat.SpFreq(2:end));
            
        else % all the bursts are complete
            
            meanFR(i,j) = mean(hcostat{i,j}.burststat.BuFreq);
            stdFR(i,j) = std(hcostat{i,j}.burststat.BuFreq);
            
            meanSpkFR(i,j) = mean(hcostat{i,j}.burststat.SpFreq);
            stdSpkFR(i,j) = std(hcostat{i,j}.burststat.SpFreq);
            
        end
        
        CVFR(i,j) = stdFR(i,j)/meanFR(i,j);
        CVSpkFR(i,j) = stdSpkFR(i,j)/meanSpkFR(i,j);
        
    end
end

oute = [{CVFR},{CVSpkFR}];

%% ----------- unity line plot (10C vs 20C) ---------------------

clf
g=0.07; l=0.08;

col = colormaps.linspecer(11);
xx=0:0.01:5; yy=xx;
for ii=1:2 % release/escape
    for i=1:2 % number of characteristics
        display.bigsubplot(2,2,ii,i,g,l)
        if ii==1
            for j = 1:length(Data.r.V) % release
                h1=plot(outr{i}(j,1),outr{i}(j,2),'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.8, 0.3, 0],'MarkerEdgeColor','k'); hold on
                hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
            end
            title('Release')
        else
            for j = 1:length(Data.e.V) % escape
                h1=plot(oute{i}(j,1),oute{i}(j,2),'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.36,0.23,0.6],'MarkerEdgeColor','k'); hold on
                hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
            end
            title('Escape')
        end
        if i==1
            xlabel('CV of cycle frequency at low Temp'); ylabel('CV of cycle frequency at high Temp');
        elseif i==2
            xlabel('CV of spike frequency at low Temp'); ylabel('CV of spike frequency at high Temp');
        else
            xlabel('CV of # spikes/burst at low Temp'); ylabel('CV of # spikes/burst at high Temp');
        end
        ylim([0 0.21]); xlim([0 0.21]);
        box off; axis square
        set(gca,'Fontsize',12,'FontName','Arial')
        
    end
end

%% stats
% Wilcoxon signed rank test (non parametric for paired samples)

for i=1:2 % number of characteristics
    meanoutr(i,:) = [mean(outr{i}(:,1)), mean(outr{i}(:,2))];
    stdoutr(i,:) = [std(outr{i}(:,1)), std(outr{i}(:,2))];
    meanoute(i,:) = [mean(oute{i}(:,1)), mean(oute{i}(:,2))];
    stdoute(i,:) = [std(oute{i}(:,1)), std(oute{i}(:,2))];
    
    pr(i) = signrank(outr{i}(:,1),outr{i}(:,2));
    pe(i) = signrank(oute{i}(:,1),oute{i}(:,2));    
    
end
