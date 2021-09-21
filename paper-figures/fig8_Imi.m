% load data

load('C:\Users\moroz\Documents\code\hco_analysis\data\fig8_data_Imi.mat');

%% calculate properties of the circuit output in control and with Imi

ERQ=[]; hcostat=[];
FR=[]; dc=[]; A=[]; SpkFR=[]; Nspks=[]; ERQ1=[]; stdFR1=[]; stdFR=[];
meanFR=[]; meandc=[]; meanA=[]; meanSpkFR=[]; meanNspks=[]; meanERQ=[];

for i=1:numel(data.V) % experiment
    i
    for j=1:2 % gMI: 0, 150 nS
        for jj=1:numel(data.V{i}{j}) % synaptic threhold
            for ii=1:2 % neuron
                if (i==6 && jj==13) || (i==8 && jj==12) % for short traces, take the whole trace for the analysis
                    ERQ{i}{j}(jj,ii) = analysis.erq(data.V{i}{j}{jj}(ii,:),data.Vth{i}{j}(jj));
                    [hcostat{i}{j}{jj}{ii}] = analysis.hco_stat(data.V{i}{j}{jj}(ii,:),data.Fs{i});
                else
                    ERQ{i}{j}(jj,ii) = analysis.erq(data.V{i}{j}{jj}(ii,5*data.Fs{i}:end),data.Vth{i}{j}(jj));
                    [hcostat{i}{j}{jj}{ii}] = analysis.hco_stat(data.V{i}{j}{jj}(ii,5*data.Fs{i}:end),data.Fs{i});
                end
                FR{i,j}(jj,ii) = 1./hcostat{i}{j}{jj}{ii}.T1_mean;
                dc{i,j}(jj,ii) = hcostat{i}{j}{jj}{ii}.dc_mean;
                A{i,j}(jj,ii) = hcostat{i}{j}{jj}{ii}.A;
                SpkFR{i,j}(jj,ii) = hcostat{i}{j}{jj}{ii}.SpkFreq_mean;
                Nspks{i,j}(jj,ii) = hcostat{i}{j}{jj}{ii}.nSpks_mean;
                ERQ1{i,j}(jj,ii) = ERQ{i}{j}(jj,ii);
                
                stdFR1{i,j}(jj,ii) = std(1./hcostat{i}{j}{jj}{ii}.T);
            end
            
            % mean characteristics across two neurons
            meanFR{i,j}(jj) = mean(FR{i,j}(jj,:),2);
            meandc{i,j}(jj) = mean(dc{i,j}(jj,:),2);
            meanA{i,j}(jj) = mean(A{i,j}(jj,:),2);
            meanSpkFR{i,j}(jj) = mean(SpkFR{i,j}(jj,:),2);
            meanNspks{i,j}(jj) = mean(Nspks{i,j}(jj,:),2);
            meanERQ{i,j}(jj) = mean(ERQ1{i,j}(jj,:),2);
            
            stdFR{i,j}(jj) = mean(stdFR1{i,j}(jj,:),2);
        end
    end
end

out=[];
out = [{meanFR}, {meanA}, {meandc}, {meanNspks}, {meanSpkFR}];

%% circuit properties in release and escape only
% indexes for the synaptic threholds corresponding to escape and release
ictrl_esc=[2,3,2,3,2,3,2,3]; ictrl_rel=[8,8,6,7,8,9,6,8];
iimi_esc=[2,3,2,2,2,2,2,2]; iimi_rel=[6,8,6,6,8,8,8,7];

outesc=[]; outrel=[];

for ii=1:numel(out) % circuit cheracteristics
    for i=1:numel(data.V) % experiment
        outesc{ii}(i,1) = out{ii}{i,1}(ictrl_esc(i));
        outesc{ii}(i,2) = out{ii}{i,2}(iimi_esc(i));
        
        outrel{ii}(i,1) = out{ii}{i,1}(ictrl_rel(i));
        outrel{ii}(i,2) = out{ii}{i,2}(iimi_rel(i));
    end
end

%% plot example traces from one experiment and cycle frequency vs synaptic threhold

% example traces
clf
i=6;
% escape (Vth = -50 mV)

xctrl=1/data.Fs{i}:1/data.Fs{i}:length(data.V{i}{1}{3}(1,5*data.Fs{i}:end))/data.Fs{i};
ximi=1/data.Fs{i}:1/data.Fs{i}:length(data.V{i}{2}{2}(1,5*data.Fs{i}:end))/data.Fs{i};

for ii=1:2  
    subtightplot(8,4,(ii-1)*4+1,g,l,l) % ctrl
    plot(xctrl,data.V{i}{1}{3}(ii,5*data.Fs{i}:end),'k'); hold on;
    display.plothorzline(data.Vth{i}{1}(3));
    ylim([-60 -10]); xlim([5 30]);
    if ii==1; title(['Vth=',num2str(data.Vth{i}{1}(3)),'mV']); end
    set(gca,'Fontsize',12); axis off; box off
    
    subtightplot(8,4,(ii-1)*4+9,g,l,l) % +Imi
    plot(ximi,data.V{i}{2}{2}(ii,5*data.Fs{i}:end),'color',[0 0.55 0.55]); hold on;   
    display.plothorzline(data.Vth{i}{2}(2)); 
    ylim([-60 -10]); xlim([5 30]);
    set(gca,'Fontsize',14); axis off; box off   
end

% mixed mechanism (Vth = -42 mV)

xctrl=1/data.Fs{i}:1/data.Fs{i}:length(data.V{i}{1}{7}(1,5*data.Fs{i}:end))/data.Fs{i};
ximi=1/data.Fs{i}:1/data.Fs{i}:length(data.V{i}{2}{6}(1,5*data.Fs{i}:end))/data.Fs{i};

for ii=1:2   
    subtightplot(8,4,ii*4-2,g,l,l) % ctrl
    plot(xctrl,data.V{i}{1}{7}(ii,5*data.Fs{i}:end),'k'); hold on;
    display.plothorzline(data.Vth{i}{1}(7));
    ylim([-60 -10]); xlim([5 30]);
    if ii==1; title(['Vth=',num2str(data.Vth{i}{1}(7)),'mV']); end
    set(gca,'Fontsize',12); axis off; box off
    
    subtightplot(8,4,ii*4+6,g,l,l) % +Imi
    plot(ximi,data.V{i}{2}{6}(ii,5*data.Fs{i}:end),'color',[0 0.55 0.55]); hold on;   
    display.plothorzline(data.Vth{i}{2}(6)); 
    ylim([-60 -10]); xlim([5 30]);
    set(gca,'Fontsize',14); axis off; box off   
end

% release (Vth = -36 mV)

xctrl=1/data.Fs{i}:1/data.Fs{i}:length(data.V{i}{1}{10}(1,5*data.Fs{i}:end))/data.Fs{i};
ximi=1/data.Fs{i}:1/data.Fs{i}:length(data.V{i}{2}{9}(1,5*data.Fs{i}:end))/data.Fs{i};

for ii=1:2   
    subtightplot(8,4,ii*4-1,g,l,l) % ctrl
    plot(xctrl,data.V{i}{1}{10}(ii,5*data.Fs{i}:end),'k'); hold on;
    display.plothorzline(data.Vth{i}{1}(10));
    ylim([-60 -10]); xlim([5 30]);
    if ii==1; title(['Vth=',num2str(data.Vth{i}{1}(10)),'mV']); end
    set(gca,'Fontsize',12); axis off; box off
    
    subtightplot(8,4,ii*4+7,g,l,l) % +Imi
    plot(ximi,data.V{i}{2}{9}(ii,5*data.Fs{i}:end),'color',[0 0.55 0.55]); hold on;   
    display.plothorzline(data.Vth{i}{2}(9)); 
    ylim([-60 -10]); xlim([5 30]);
    set(gca,'Fontsize',14); axis off; box off   
end


% plot cycle frequency as a function of synaptic threshold

subtightplot(2,4,4,g,l,l)
errorbar(data.Vth{i}{1},meanFR{i,1}, stdFR{i,1},'.-','markersize',20,'linewidth',2,'color','k'), hold on
errorbar(data.Vth{i}{2},meanFR{i,2}, stdFR{i,2},'.-','markersize',20,'linewidth',2,'color',[0 139/255 139/255])
legend('control','gMI=150 nS','location','northwest'); box off
ylabel('Cycle frequency, Hz'); xlabel('Synaptic threhold, mV')
set(gca,'Fontsize',14); 
xlim([-53 -33]); xticks([-53:2:-34])
ylim([0 0.55]); yticks([0:0.1:0.5])


% ----------- unity line plot (ctrl vs modulator) -----------------------

col=display.linspecer(8);
xx=-100:1:150; yy=xx;
% cycle frequency
for ii = 1:numel(out)
    subplot(2,5,ii+5)
    h1=plot(outesc{ii}(:,1),outesc{ii}(:,2),'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.36,0.23,0.6],'MarkerEdgeColor','k'); hold on
    h2=plot(outrel{ii}(:,1),outrel{ii}(:,2),'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.8, 0.3, 0],'MarkerEdgeColor','k'); hold on
    hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
    % legend([h1,h2],'Escape','Release','location','northwest')
    set(gca,'Fontsize',14)
    box off; axis square
    switch ii
        case 1
            xlabel('Cycle freq in ctrl, Hz'); ylabel('Cycle freq with I_{MI}, Hz');
            xticks([0:0.1:0.4]); yticks([0:0.1:0.4])
            ylim([0 0.4]); xlim([0 0.4]);
        case 2
            xlabel('Amplitude in ctrl, mV'); ylabel('Amplitude with I_{MI}, mV');
            xticks([0:5:35]); yticks([0:5:35])
            ylim([0 35]); xlim([0 35]);
        case 3
            xlabel('Duty cycle in ctrl, %'); ylabel('Duty cycle with I_{MI}, %');
            xticks([0:10:50]); yticks([0:10:50])
            ylim([0 50]); xlim([0 50]);
        case 4
            xlabel('# spikes/burst in ctrl'); ylabel('# spikes per burst with I_{MI}');
            xticks([0:5:25]); yticks([0:5:25])
            ylim([0 25]); xlim([0 25]);
        case 5
            xlabel('Spike freq in ctrl, Hz'); ylabel('Spike freq with I_{MI}, Hz');
            xticks([0:2:12]); yticks([0:2:12])
            ylim([0 12]); xlim([0 12]);
    end
end

%% paired t-tests
% cycle frequency
[h,pFR_esc] = ttest(FRctrl_esc,FRimi_esc)
[h,pFR_rel] = ttest(FRctrl_rel,FRimi_rel)

[h,pA_esc] = ttest(meanActrl_esc,meanAimi_esc)
[h,pA_rel] = ttest(meanActrl_rel,meanAimi_rel)
[h,pA_esc_rel] = ttest(meanAimi_esc,meanAimi_rel) % comparison between esc and rel with Imi

[h,p_dc_esc] = ttest(meandcctrl_esc,meandcimi_esc) 
[h,p_dc_rel] = ttest(meandcctrl_rel,meandcimi_rel) 

[h,p_nspk_esc] = ttest(meanNspkctrl_esc,meanNspkimi_esc) 
[h,p_nspk_rel] = ttest(meanNspkctrl_rel,meanNspkimi_rel) 

[h,p_spkfr_esc] = ttest(meanspkFRctrl_esc,meanspkFRimi_esc) 
[h,p_spkfr_rel] = ttest(meanspkFRctrl_rel,meanspkFRimi_rel)

%[h,p_mode] = ttest(dpermodeesc_mean,dpermoderel_mean); % mode
%[h,p_mode] = ttest(dpermodeesc_mean,dpermoderel_mean); % mode