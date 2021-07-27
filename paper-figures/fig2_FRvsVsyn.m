
% load data for figure 2

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig2.mat')

% break the membrane potential recordings into segments based on the synaptic threshold
Vseg=cell(1,16); Vth1=cell(1,16);
for j = 1:numel(hco.V1)
    [Vseg{j}{1},Vth1{j}] = analysis.Vsegments(hco.V1{j},hco.Vth{j});
    [Vseg{j}{2},Vth1{j}] = analysis.Vsegments(hco.V2{j},hco.Vth{j});
end

% cleanup the dataset -> remove two extra repeated thresholds
Vth1{1}(11:12)=[];
Vseg{1}{1}(12)=[]; Vseg{1}{1}(11)=[];
Vseg{1}{2}(12)=[]; Vseg{1}{2}(11)=[];
% remove repeated -54mV threhold in the middle
Vth1{2}(3)=[]; Vseg{2}{1}(3)=[];  Vseg{2}{2}(3)=[];   

 %% calculate HCO characteristcs as a function of the synaptic threhold
 hcostat=[]; ERQ=[];
 for j = 1:length(folder) %experiment number
     j
     for ii = 1:2 % neuron number
         for i = 1:length(Vseg{j}{ii}) % # of Vth            
             if ~isempty(Vseg{j}{ii}{i})==1 && length(Vseg{j}{ii}{i})>10*Fs(j)
                 V = Vseg{j}{ii}{i}(5*Fs(j):end);
                 ERQ{j}{ii}(i) = (mean(V)-Vth1{j}(i))/mean(V);
                 [hcostat{j}{ii}{i}] = hco_stat(V, Fs(j));
             else
                 hcostat{j}{ii}{i} = [];  ERQ{j}{ii}(i) = NaN;
             end
         end
     end
 end

 %% calulate burst frequency
 FR=[]; T=[]; A=[]; SpkFR=[]; Nspks=[]; DC=[]; burststat=[]; st=[];
 for j=1:length(folder) % experiment
     for ii=1:2 % neurons number
         for i = 1:length(Vseg{j}{ii}) % # of Vth
             if ~isempty(hcostat{j}{ii}{i})
             FR{j}{ii}(i) = 1./hcostat{j}{ii}{i}.T1_mean; % cycle frequency
             T{j}{ii}(i) = hcostat{j}{ii}{i}.T1_mean; % cycle period
             A{j}{ii}(i) = hcostat{j}{ii}{i}.A; % amplitude
             SpkFR{j}{ii}(i) = hcostat{j}{ii}{i}.SpkFreq_mean; % spike frequency
             Nspks{j}{ii}(i) = hcostat{j}{ii}{i}.nSpks_mean; % # spks/burst
             DC{j}{ii}(i) = hcostat{j}{ii}{i}.dc_mean; % #duty cycle
             burststat{j}{ii}{i} = hcostat{j}{ii}{i}.burststat; % #duty cycle
             st{j}{ii}{i} = hcostat{j}{ii}{i}.st; % #duty cycle
             else
             FR{j}{ii}(i)=NaN; T{j}{ii}(i)=NaN;  A{j}{ii}(i)=NaN;
             SpkFR{j}{ii}(i)=NaN; Nspks{j}{ii}(i)=NaN; DC{j}{ii}(i)=NaN;
             st{j}{ii}{i}=[]; burststat{j}{ii}{i}=[]; 
             end
         end
     end
 end

%% calculate all the measures across two neurons
FRmean=[]; Amean=[]; SpkFRmean=[]; Nspksmean=[]; DCmean=[]; ERQmean=[];
for j = 1:length(folder) % experiment
    for i = 1:length(Vseg{j}{1}) % # of Vth
        FRmean{j}(i)=nanmean([FR{j}{1}(i),FR{j}{2}(i)]);
        Amean{j}(i)=nanmean([A{j}{1}(i),A{j}{2}(i)]);
        SpkFRmean{j}(i)=nanmean([SpkFR{j}{1}(i),SpkFR{j}{2}(i)]);
        Nspksmean{j}(i)=nanmean([Nspks{j}{1}(i),Nspks{j}{2}(i)]);        
        DCmean{j}(i)=nanmean([DC{j}{1}(i),DC{j}{2}(i)]);     
        ERQmean{j}(i)=nanmean([ERQ{j}{1}(i),ERQ{j}{2}(i)]);
    end
    % if the circuit not a half-center set measures to NaN
    Amean{j}(isinf(FRmean{j}))=NaN;
    Nspksmean{j}(isinf(FRmean{j}))=NaN;
    SpkFRmean{j}(isinf(FRmean{j}))=NaN;
    DCmean{j}(isinf(FRmean{j}))=NaN;
    ERQmean{j}(isinf(FRmean{j}))=NaN;
end
%% plot ERQ vs Vth
clf
subplot(1,2,1) % plot example ERQ and fit sigmoidal function
plot(Vth1{16},ERQmean{16},'.-','markersize',20,'linewidth',1.6,'color','k'), hold on
axis square
ylim([-0.1 0.2]); xlim([-54 -32])
yticks([-0.1:0.05:0.2]); xticks([-54:2:-34])
ylabel('ERQ'); xlabel('Synaptic threshold (V_{th}), mV')
set(gca,'Fontsize',16,'FontName','Arial')
x=Vth1{i}(~isnan(ERQmean{i}));
f = @(F,x) F(1) + F(2)./(1+exp(-(x-F(3))./F(4))); % fit sigmoidal function
[p,R]= nlinfit(Vth1{16}(~isnan(ERQmean{16})),ERQmean{16}(~isnan(ERQmean{16})),f,[-0.08 0.25 -44 3])
fitnlm(Vth1{16}(~isnan(ERQmean{16})),ERQmean{16}(~isnan(ERQmean{16})),f,[-0.08 0.25])

hold on, plot(x,f(p,x),'linewidth',2,'color','c')
title('951-060')

%col=linspecer(length(Vth1));
col = cbrewer('qual','Set1',length(Vth1)); 
subplot(1,2,2) % all tje data
for i = 1:length(folder) %[4,5,10,12,13,15,16,14]
plot(Vth1{i},ERQmean{i},'.-','markersize',20,'linewidth',1.6,'color',col(i,:,:)), hold on
ylim([-0.1 0.2]); xlim([-54 -32])
end
axis square
text(-53,0.17,['N=',num2str(length(folder))],'Fontsize',16,'FontName','Arial')
yticks([-0.1:0.05:0.2]); xticks([-54:2:-34])
ylabel('ERQ'); xlabel('Synaptic threshold (V_{th}), mV')
set(gca,'Fontsize',16,'FontName','Arial')
%% interpolate all the measures so that they are in the same boundaries (-54 to -28)
Vthint=[]; FRint=[]; Aint=[]; DCint=[];  SpkFRint=[]; Nspksint=[]; ERQint=[];
for n=1:2 % neuron number
    for i = 1:length(folder) % experiment
        Vthall =[-54:2:-28];
        Vthint{n}(i,:) = interp1(Vth1{i},Vth1{i},Vthall,'linear','extrap'); % threshold
        FRint{n}(i,:) = interp1(Vth1{i},FR{i}{n},Vthall); % cycle frequency
        Aint{n}(i,:) = interp1(Vth1{i},A{i}{n},Vthall); % amplitude
        DCint{n}(i,:) = interp1(Vth1{i},DC{i}{n},Vthall); % duty cycle
        SpkFRint{n}(i,:) = interp1(Vth1{i},SpkFR{i}{n},Vthall); % spike frequency
        Nspksint{n}(i,:) = interp1(Vth1{i},Nspks{i}{n},Vthall); % # spikes/burst
        ERQint{n}(i,:) = interp1(Vth1{i},ERQ{i}{n},Vthall); % # ERQ
    end
    FRint{n}(FRint{n}==0)=NaN; Aint{n}(Aint{n}==0)=NaN;
    DCint{n}(DCint{n}==0)=NaN; SpkFRint{n}(SpkFRint{n}==0)=NaN;
    Nspksint{n}(Nspksint{n}==0)=NaN;
end

%% calculate interpolated values across two neurons
FRmeanint=[]; Ameanint=[]; SpkFRmeanint=[]; Nspksmeanint=[]; DCmeanint=[]; ERQmeanint=[];
for j = 1:length(folder) % experiment
        FRmeanint(j,:)=nanmean([FRint{1}(j,:);FRint{2}(j,:)]);
        Ameanint(j,:)=nanmean([Aint{1}(j,:);Aint{2}(j,:)]);
        SpkFRmeanint(j,:)=nanmean([SpkFRint{1}(j,:);SpkFRint{2}(j,:)]);
        Nspksmeanint(j,:)=nanmean([Nspksint{1}(j,:);Nspksint{2}(j,:)]);
        DCmeanint(j,:)=nanmean([DCint{1}(j,:);DCint{2}(j,:)]);
        ERQmeanint(j,:)=nanmean([ERQint{1}(j,:);ERQint{2}(j,:)]);
end
%% if the circuit not a half-center set measures to NaN
for j = 1:length(folder)
Ameanint(isnan(FRmeanint) | isinf(FRmeanint))=NaN;
Nspksmeanint(isnan(FRmeanint) | isinf(FRmeanint))=NaN;
SpkFRmeanint(isnan(FRmeanint) | isinf(FRmeanint))=NaN;
DCmeanint(isnan(FRmeanint) | isinf(FRmeanint))=NaN;
ERQmeanint(isnan(FRmeanint) | isinf(FRmeanint))=NaN;
end

% if values in one of the neurons are NaN
FRmeanint(isnan(FRint{1}) | isinf(FRint{2}))=NaN;
Ameanint(isnan(Aint{1}) | isinf(Aint{2}))=NaN;
Nspksmeanint(isnan(Nspksint{1}) | isinf(Nspksint{2}))=NaN;
SpkFRmeanint(isnan(SpkFRint{1}) | isinf(SpkFRint{2}))=NaN;
DCmeanint(isnan(DCint{1}) | isinf(DCint{2}))=NaN;
ERQmeanint(isnan(ERQint{1}) | isinf(ERQint{2}))=NaN;


%% correlation coefficient for the cycle frequency and the amplitude
corrcoef(FRmeanint(~isnan(Ameanint)),Ameanint(~isnan(Ameanint)))
clf ,scatter(FRmeanint(~isnan(Ameanint)),Ameanint(~isnan(Ameanint)),'filled')

%% stats (one way repeated measures ANOVA) for DC and spike frequency
[p,tbl,stats]=anova1(DCmeanint,[])
%[p,tbl,stats]=anova1(SpkFRmeanint,[],'off')

%comp=multcompare(stats,'display','off');
%group=(comp(:,1:2));
%hcomp=comp(:,6)<0.05;
%pvalsig=comp(hcomp,6);
%sigpoints=([group(hcomp==1,1), group(hcomp==1,2)])
%[~,idx] = unique(sort(sigpoints,2),'rows','stable');
%pvalsig_dc=pvalsig(idx);
%sigpoints1_dc = sigpoints(idx,:)
%groups_dc=(num2cell(sigpoints1_dc,2));

%% repeated measures ANOVA spike frequency
t=table([1:14]',DCmeanint(1,:)',DCmeanint(2,:)',DCmeanint(3,:)',DCmeanint(4,:)',...
    DCmeanint(5,:)',DCmeanint(6,:)',DCmeanint(7,:)',DCmeanint(8,:)',DCmeanint(9,:)',...
    DCmeanint(10,:)',DCmeanint(11,:)',DCmeanint(12,:)',DCmeanint(13,:)',DCmeanint(14,:)',...
    'VariableNames',{'dc','Vth1','Vth2','Vth3','Vth4',...
    'Vth5','Vth6','Vth7','Vth8','Vth9','Vth10','Vth11','Vth12','Vth13','Vth14'});
Meas = table([1 2 3 4 5 6 7 8 9 10 11 12 13 14]','VariableNames',{'Measurements'});
rm=fitrm(t,'Vth1-Vth14~dc','WithinDesign',Meas);
ranovatbl = ranova(rm)
%%
t=table([1:50]',DCint1(1,:)',DCint1(2,:)',DCint1(3,:)',DCint1(4,:)',...
    DCint1(5,:)',DCint1(6,:)',DCint1(7,:)',DCint1(8,:)',DCint1(9,:)',...
    DCint1(10,:)',DCint1(11,:)',DCint1(12,:)',DCint1(13,:)',DCint1(14,:)',...
    'VariableNames',{'dc','Vth1','Vth2','Vth3','Vth4',...
    'Vth5','Vth6','Vth7','Vth8','Vth9','Vth10','Vth11','Vth12','Vth13','Vth14'});
Meas = table([1 2 3 4 5 6 7 8 9 10 11 12 13 14]','VariableNames',{'Measurements'});
rm=fitrm(t,'Vth1-Vth14~1','WithinDesign',Meas);
ranovatbl = ranova(rm)
%% repeated measures ANOVA spike frequency
t=table([1:14]',SpkFRmeanint(1,:)',SpkFRmeanint(2,:)',SpkFRmeanint(3,:)',SpkFRmeanint(4,:)',...
    SpkFRmeanint(5,:)',SpkFRmeanint(6,:)',SpkFRmeanint(7,:)',SpkFRmeanint(8,:)',SpkFRmeanint(9,:)',...
    SpkFRmeanint(10,:)',SpkFRmeanint(11,:)',SpkFRmeanint(12,:)',SpkFRmeanint(13,:)',SpkFRmeanint(14,:)',...
    'VariableNames',{'SpkFR','Vth1','Vth2','Vth3','Vth4',...
    'Vth5','Vth6','Vth7','Vth8','Vth9','Vth10','Vth11','Vth12','Vth13','Vth14'});
Meas = table([1 2 3 4 5 6 7 8 9 10 11 12 13 14]','VariableNames',{'Measurements'});
rm=fitrm(t,'Vth1-Vth14~1','WithinDesign',Meas);
ranovatbl = ranova(rm)
%%
t=table([1:50]',SpkFRint1(1,:)',SpkFRint1(2,:)',SpkFRint1(3,:)',SpkFRint1(4,:)',...
    SpkFRint1(5,:)',SpkFRint1(6,:)',SpkFRint1(7,:)',SpkFRint1(8,:)',SpkFRint1(9,:)',...
    SpkFRint1(10,:)',SpkFRint1(11,:)',SpkFRint1(12,:)',SpkFRint1(13,:)',SpkFRint1(14,:)',...
    'VariableNames',{'SpkFR','Vth1','Vth2','Vth3','Vth4',...
    'Vth5','Vth6','Vth7','Vth8','Vth9','Vth10','Vth11','Vth12','Vth13','Vth14'});
Meas = table([1 2 3 4 5 6 7 8 9 10 11 12 13 14]','VariableNames',{'Measurements'});
rm=fitrm(t,'Vth1-Vth14~1','WithinDesign',Meas);
ranovatbl = ranova(rm)

%% interpolate characteristics for ERQ values
FRint1=[]; ERQint1=[]; Aint1=[]; DCint1=[]; SpkFRint1=[]; Nspksint1=[]; Vthint1=[];
ERQint1 = linspace(min(ERQmeanint(:)),max(ERQmeanint(:)),50);
for i = 1:length(folder) 
    x=ERQmeanint(i,find(~isnan(ERQmeanint(i,:))));
   % Vthint1(i,:) = interp1(x,Vth1{i}(i,find(~isnan(ERQmeanint(i,:)))),ERQint1,'linear');    
    FRint1(i,:) = interp1(x,FRmeanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1,'linear');
    %FRint1(i,:) = interp1(modeint{n}(i,:),FRint{n}(i,:),modeint1);
    Aint1(i,:) = interp1(x,Ameanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1);
    DCint1(i,:) = interp1(x,DCmeanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1);
    SpkFRint1(i,:) = interp1(x,SpkFRmeanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1);
    Nspksint1(i,:) = interp1(x,Nspksmeanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1);
end

FRint1(FRint1==0)=NaN; Aint1(Aint1==0)=NaN; dcint1(dcint1==0)=NaN; 
Spkfrint1(Spkfrint1==0)=NaN; nspksint1(nspksint1==0)=NaN;

% % if sample size is to small for averaging
% for i=1:100
% if length(FRint(~isnan(FRint(:,i)),i))<5
%     FRint(:,i)=NaN; Aint(:,i)=NaN; dcint(:,i)=NaN; Spkfrint(:,i)=NaN; nspksint(:,i)=NaN;
% end
% end
%%
%[rho,pval]=corr(nanmean(SpkFRint1)',ERQint1','Type','Spearman') 
[rho,pval]=corr(nanmean(DCint1)',ERQint1','Type','Spearman') 
%% plot all the characteristics of half-center for all the experiments
%col = cbrewer('qual','Set1',length(Vth1)); 
clf,
subplot(1,5,1) % cycle frequency
plot(ERQmeanint',FRmeanint','linewidth',1,'color','k'), hold on
plot(ERQint1,nanmean(FRint1),'r','linewidth',2)
ylim([0.05 0.65]); xlim([-0.1 0.2])
xticks([-0.1:0.05:0.2]); yticks([0:0.1:0.7])
text(-0.08,0.6,['N=',num2str(length(folder))],'Fontsize',16,'FontName','Arial')
xlabel('ERQ'); ylabel('Cycle freq, Hz')
set(gca,'Fontsize',16,'FontName','Arial')
axis square

subplot(1,5,2) % slow-wave amplitude
plot(ERQmeanint',Ameanint','linewidth',1,'color','k'), hold on
plot(ERQint1,nanmean(Aint1),'r','linewidth',2)
ylim([5 30]); xlim([-0.1 0.2])
xticks([-0.1:0.05:0.2]); yticks([0:5:30])
xlabel('ERQ'); ylabel('Amplitude, mV')
set(gca,'Fontsize',16,'FontName','Arial')
axis square

subplot(1,5,3) % duty cycle
plot(ERQmeanint',DCmeanint','linewidth',1,'color','k'), hold on
plot(ERQint1,nanmean(DCint1),'r','linewidth',2)
ylim([5 50]); xlim([-0.1 0.2])
xticks([-0.1:0.05:0.2]); yticks([0:5:60])
xlabel('ERQ'); ylabel('Duty cycle,%')
set(gca,'Fontsize',16,'FontName','Arial')
axis square

subplot(1,5,4) % number of spikes per burst
plot(ERQmeanint',Nspksmeanint','linewidth',1,'color','k'), hold on
plot(ERQint1,nanmean(Nspksint1),'r','linewidth',2)
xlim([-0.1 0.2]); ylim([0 20])
xticks([-0.1:0.05:0.2]); yticks([0:5:25])
xlabel('ERQ'); ylabel('# spikes/burst')
set(gca,'Fontsize',16,'FontName','Arial')
axis square

subplot(1,5,5) % spike frequency
plot(ERQmeanint',SpkFRmeanint','linewidth',1,'color','k'), hold on
plot(ERQint1,nanmean(SpkFRint1),'r','linewidth',2)
xlim([-0.1 0.2]); ylim([3 10])
xticks([-0.1:0.05:0.2]); yticks([0:1:10])
xlabel('ERQ'); ylabel('Spike freq, Hz')
set(gca,'Fontsize',16,'FontName','Arial')
axis square
%%
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'fig2_characteristics.pdf')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'map_transition','-dpdf','-r0')
%% plot example traces
clf
i=15;
ylim1=-65; ylim2=-5; xlim1=0; xlim2=40;
ii=[3,5,7,9,10]; k=1;
for jj=[1:5]
    j=1; n=1; x=1/Fs(i):1/Fs(i):length(Vseg{i}{n}{ii(k)}(Fs(i)*5:end))/Fs(i);
    subtightplot(2,5,1+(k-1),g,l,l), plot(x,Vseg{i}{n}{ii(k)}(Fs(i)*5:end),'color','k','linewidth',1)
    title([ 'Control, Vth = ',num2str(Vth1{i}(ii(k))), 'mV']);
    ylim([ylim1 ylim2]); xlim([xlim1 xlim2]);
    set(gca, 'Fontsize',16,'FontName','Arial'); box off
    h=mmyplothorzline(Vth1{i}(ii(k))); set(h,'color','k','linewidth',1), hold on
    
    j=1; n=2; x=1/Fs(i):1/Fs(i):length(Vseg{i}{n}{ii(k)}(Fs(i)*5:end))/Fs(i);
    subtightplot(2,5,6+(k-1),g,l,l), plot(x,Vseg{i}{n}{ii(k)}(Fs(i)*5:end),'color','k','linewidth',1)
    title([ 'Control, Vth = ',num2str(Vth1{i}(ii(k))), 'mV']);
    ylim([ylim1 ylim2]); xlim([xlim1 xlim2]); xlabel('Time, sec')
    set(gca,'Fontsize',16,'FontName','Arial'); box off
    h=mmyplothorzline(Vth1{i}(ii(k))); set(h,'color','k','linewidth',1), hold on
    k=k+1;
end