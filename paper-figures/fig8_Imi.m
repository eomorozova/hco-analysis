% load data

load('C:\Users\moroz\Documents\code\hco_analysis\data\fig8_data_Imi.mat');

%% calculate properties of the circuit output in control and with Imi

ERQ=[]; hcostat=[];
FR=[]; dc=[]; A=[]; SpkFR=[]; Nspks=[]; ERQ1=[];
meanFR=[]; meandc=[]; meanA=[]; meanSpkFR=[]; meanNspks=[]; meanERQ=[];

for i=1:numel(data.V) % experiment
    i
    for j=1:2 % gMI: 0, 150 nS
        for jj=1:numel(data.V{i}{j}) % synaptic threhold
            for ii=1:2 % neuron
                ERQ{i}{j}(jj,ii) = analysis.erq(data.V{i}{j}{jj}(ii,:),data.Vth{i}{j}(jj));
                [hcostat{i}{j}{jj}{ii}] = analysis.hco_stat(data.V{i}{j}{jj}(ii,:), data.Fs{j});
                
                FR(i,j,jj,ii) = 1./hcostat{i}{j}{jj}{ii}.T1_mean;
                dc(i,j,jj,ii) = hcostat{i}{j}{jj}{ii}.dc_mean;
                A(i,j,jj,ii) = hcostat{i}{j}{jj}{ii}.A;
                SpkFR(i,j,jj,ii) = hcostat{i}{j}{jj}{ii}.SpkFreq_mean;
                Nspks(i,j,jj,ii) = hcostat{i}{j}{jj}{ii}.nSpks_mean;
                ERQ1(i,j,jj,ii) = ERQ{i}{j}(jj,ii);
            end
            
            % mean characteristics across two neurons
            meanFR(i,j,jj) = mean(FR(i,j,jj,:),4);
            meandc(i,j,jj) = mean(dc(i,j,jj,:),4);
            meanA(i,j,jj) = mean(A(i,j,jj,:),4);
            meanSpkFR(i,j,jj) = mean(SpkFR(i,j,jj,:),4);
            meanNspks(i,j,jj) = mean(Nspks(i,j,jj,:),4);
            meanERQ(i,j,jj) = mean(ERQ1(i,j,jj,:),4);
        end
    end
end

out=[];
out = [{meanFR}, {meanA}, {meanSpkFR}, {meanNspks}, {meandc}, {meanERQ}];
%%
% indexes for the synaptic threholds corresponding to escape and release
ictrl_esc=[2,3,2,3,2,3,2,3]; ictrl_rel=[8,8,6,7,8,9,6,8];
iimi_esc=[2,3,2,2,2,2,2,2]; iimi_rel=[6,8,6,6,8,8,8,7];

outesc=[]; outrel=[];

for ii=1:numel(out) % circuit cheracteristics
    for i=1:numel(data.V) % experiment
        i
        for j=1:2 % gMI: 0, 150 nS
            % cycle frequency
            outesc{ii}(i,j) = out{ii}(i,j,ictrl_esc(i));
            outrel{ii}(i,j) = out{ii}(i,j,ictrl_rel(i));
        end
    end
end
    
%% ----------- unity line plot (ctrl vs modulator) -----------------------
clf
col=display.linspecer(8);
xx=-100:1:150; yy=xx;
% cycle frequency
subplot(2,6,1)
h1=plot(FRctrl_esc,FRimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
h2=plot(FRctrl_rel,FRimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
legend([h1,h2],'Escape','Release','location','northwest')
xlabel('Cycle freq in ctrl, Hz'); ylabel('Cycle freq with I_{MI}, Hz'); 
xticks([0:0.1:0.4]); yticks([0:0.1:0.4])
ylim([0 0.4]); xlim([0 0.4]); box off 
set(gca,'Fontsize',14)
box off; axis square
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

subplot(2,6,2)
%h1=plot(Actrl_esc,Aimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
%h2=plot(Actrl_rel,Aimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
h1=plot(meanActrl_esc,meanAimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
h2=plot(meanActrl_rel,meanAimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
%legend([h1,h2],'Escape','Release','location','northwest')
xlabel('Amplitude in ctrl, Hz'); ylabel('Amplitude with I_{MI}, Hz'); 
xticks([0:5:35]); yticks([0:5:35])
ylim([0 35]); xlim([0 35]); box off 
set(gca,'Fontsize',14)
box off; axis square
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

subplot(2,6,3)
%h1=plot(dcctrl_esc,dcimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
%h2=plot(dcctrl_rel,dcimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
h1=plot(meandcctrl_esc,meandcimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
h2=plot(meandcctrl_rel,meandcimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
%legend([h1,h2],'Escape','Release','location','northwest')
xlabel('Duty cycle in ctrl, %'); ylabel('Duty cycle with I_{MI}, %'); 
xticks([0:10:50]); yticks([0:10:50])
ylim([0 50]); xlim([0 50]); box off 
set(gca,'Fontsize',14)
box off; axis square
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

subplot(2,6,4)
%h1=plot(Nspkctrl_esc,Nspkimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
%h2=plot(Nspkctrl_rel,Nspkimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
h1=plot(meanNspkctrl_esc,meanNspkimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
h2=plot(meanNspkctrl_rel,meanNspkimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
%legend([h1,h2],'Escape','Release','location','northwest')
xlabel('# spikes/burst in ctrl'); ylabel('# spikes per burst with I_{MI}'); 
xticks([0:5:25]); yticks([0:5:25])
ylim([0 25]); xlim([0 25]); box off 
set(gca,'Fontsize',14)
box off; axis square
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

subplot(2,6,5)
%h1=plot(spkFRctrl_esc,spkFRimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
%h2=plot(spkFRctrl_rel,spkFRimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
h1=plot(meanspkFRctrl_esc,meanspkFRimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
h2=plot(meanspkFRctrl_rel,meanspkFRimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
%legend([h1,h2],'Escape','Release','location','northwest')
xlabel('Spike freq in ctrl, Hz'); ylabel('Spike freq with I_{MI}, Hz'); 
xticks([0:2:12]); yticks([0:2:12])
ylim([0 12]); xlim([0 12]); box off 
set(gca,'Fontsize',14)
box off; axis square
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

subplot(2,6,6)
h1=plot(meanmodectrl_esc,meanmodeimi_esc,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0.850, 0.4, 0.098],'MarkerEdgeColor','k'); hold on
h2=plot(meanmodectrl_rel,meanmodeimi_rel,'o','markersize',7,'linewidth',1.5,'MarkerFaceColor',[0 0.7 0.3],'MarkerEdgeColor','k'); hold on
hold on, plot(xx,yy,'k','linewidth',1,'linestyle','--')
h=mmyplothorzline(0); set(h,'color','k','linewidth',1,'linestyle','--')
h=mmyplotvertline(0); set(h,'color','k','linewidth',1,'linestyle','--')
%legend([h1,h2],'Escape','Release','location','northwest')
xlabel('ERQ in ctrl'); ylabel('ERQ with I_{MI}'); 
xticks([-0.1:0.05:0.2]); yticks([-0.1:0.05:0.2])
ylim([-0.12 0.18]); xlim([-0.12 0.18]); box off 
set(gca,'Fontsize',14)
box off; axis square
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

%% paired ttests
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