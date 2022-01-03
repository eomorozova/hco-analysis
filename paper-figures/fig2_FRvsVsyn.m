
% load the data for figure 2

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig2.mat')

%% divide the membrane potential recordings into segments based on the synaptic threshold

N = numel(data.V1); % number of experiments

Vseg=cell(1,N); Vth1=cell(1,N);

for j = 1:N
    [Vseg{j}{1},Vth1{j}] = analysis.Vsegments(data.V1{j},data.Vth{j});
    [Vseg{j}{2},Vth1{j}] = analysis.Vsegments(data.V2{j},data.Vth{j});
end

% cleanup the data -> remove two extra repeated thresholds
Vth1{1}(11:12)=[];
Vseg{1}{1}(12)=[]; Vseg{1}{1}(11)=[];
Vseg{1}{2}(12)=[]; Vseg{1}{2}(11)=[];
% remove repeated -54mV threhold in the middle
Vth1{2}(3)=[]; Vseg{2}{1}(3)=[];  Vseg{2}{2}(3)=[];   

%% calculate half-center output characteristcs as a function of the synaptic threhold
 
hcostat=cell(1,N); ERQ=cell(1,N);

for j = 1:N % experiment
    j
    for ii = 1:2 % neuron number
        for i = 1:length(Vseg{j}{ii}) % # of Vth
            if ~isempty(Vseg{j}{ii}{i})==1 && length(Vseg{j}{ii}{i})>10*data.Fs(j)
                V = Vseg{j}{ii}{i}(5*data.Fs(j):end);
                ERQ{j}{ii}(i) = analysis.erq(V, Vth1{j}(i));
                [hcostat{j}{ii}{i}] = analysis.hco_stat(V, data.Fs(j));
            else
                hcostat{j}{ii}{i} = [];  ERQ{j}{ii}(i) = NaN;
            end
        end
    end
end

FR = cell(1,N); A = cell(1,N); SpkFR = cell(1,N); Nspks = cell(1,N); DC = cell(1,N);
SpkFR1 = cell(1,N);

for j = 1:N % experiment
    for ii = 1:2 % neuron number
        for i = 1:length(Vseg{j}{ii}) % # of Vth
            if ~isempty(hcostat{j}{ii}{i})
                FR{j}{ii}(i) = 1./hcostat{j}{ii}{i}.T1_mean; % cycle frequency
                A{j}{ii}(i) = hcostat{j}{ii}{i}.A; % slow-wave amplitude
                SpkFR{j}{ii}(i) = hcostat{j}{ii}{i}.SpkFreq_mean; % spike frequency
                SpkFR1{j}{ii}(i) = hcostat{j}{ii}{i}.SpkFreq_median; % spike frequency
                Nspks{j}{ii}(i) = hcostat{j}{ii}{i}.nSpks_mean; % # spikes/burst
                DC{j}{ii}(i) = hcostat{j}{ii}{i}.dc_mean; % #duty cycle
            else
                FR{j}{ii}(i)=NaN; A{j}{ii}(i)=NaN;
                SpkFR{j}{ii}(i)=NaN; Nspks{j}{ii}(i)=NaN; DC{j}{ii}(i)=NaN;
            end
        end
    end
end

%% calculate ERQ across two neursons in a circuit

ERQmean = cell(1,N);

for j = 1:N % experiment
    for i = 1:length(Vseg{j}{1}) % # of Vth
        ERQmean{j}(i) = nanmean([ERQ{j}{1}(i),ERQ{j}{2}(i)]);
    end
    % if the circuit is not a half-center set ERQ to NaN
    ERQmean{j}(isnan(FR{j}{1}) | isnan(FR{j}{2})) = NaN;
end

%% plot example traces (Figure 2A)

clf
i=15; ii=[3,5,7,9,10]; k=1;
for jj=[1:5] % 5 traces
    for n=1:2 % neuron
    x=1/data.Fs(j):1/data.Fs(j):length(Vseg{i}{n}{ii(k)}(data.Fs(j)*5:end))/data.Fs(j);
    subplot(4,5,1+5*(n-1)+(k-1))
    plot(x,Vseg{i}{n}{ii(k)}(data.Fs(j)*5:end),'color','k','linewidth',1)
    Vmean = mean(Vseg{i}{n}{ii(k)}(data.Fs(j)*5:end));
    hold on, h = display.plothorzline(Vmean); set(h,'linestyle','-')
    if n==1; title(['Vth = ',num2str(Vth1{i}(ii(k))), 'mV']); end
    ylim([-60 -10]); xlim([0 26]);
    set(gca, 'Fontsize',14,'FontName','Arial'); box off
    h=display.plothorzline(Vth1{i}(ii(k))); hold on
    end
    k=k+1;
end

%% plot ERQ vs Vth and fit a sigmoid function (Figure 2C)

clf; i=15;
subplot(1,2,1) % plot example ERQ and fit sigmoidal function
plot(Vth1{i},ERQmean{i},'.','markersize',20,'linewidth',1.6,'color','k'), hold on
axis square; ylim([-0.1 0.2]); xlim([-54 -34])
yticks([-0.1:0.05:0.2]); xticks([-54:2:-34])
ylabel('ERQ'); xlabel('Synaptic threshold (V_{th}), mV')
set(gca,'Fontsize',16,'FontName','Arial')

%x=Vth1{i}(~isnan(ERQmean{i}));
x = [-52:0.2:-36];
ERQmeanint = interp1(Vth1{i},ERQmean{i},x); % interpolate
f = @(F,x) F(1) + F(2)./(1+exp(-(x-F(3))./F(4))); % fit sigmoid function
[p,R]= nlinfit(Vth1{i}(~isnan(ERQmean{i})),ERQmean{i}(~isnan(ERQmean{i})),f,[-0.08 0.25 -44 3]);

hold on, plot(x,f(p,x),'linewidth',2,'color','c')
hold on, plot(x(1:end-2),1000*diff(diff(f(p,x))),'linewidth',2,'color','r')
[~,the] = (max(1000*diff(diff(f(p,x)))));
[~,thr] = (min(1000*diff(diff(f(p,x)))));
hold on, display.plotvertline(x(the));
hold on, display.plotvertline(x(thr));
hold on, display.plothorzline(ERQmeanint(the));
hold on, display.plothorzline(ERQmeanint(thr));
title(data.folder(i,:))

% fit sigmoid function to ERQ vs Vth plot for each experiment
x = cell(1,N); ERQmeanint = cell(1,N);

f = @(F,x) F(1) + F(2)./(1+exp(-(x-F(3))./F(4))); % fit sigmoid function
for i=1:N
    x1 = Vth1{i}(~isnan(ERQmean{i}));
    x{i} = [x1(1):0.2:x1(end)];
    [p(i,:),R]= nlinfit(Vth1{i}(~isnan(ERQmean{i})),ERQmean{i}(~isnan(ERQmean{i})),f,[-0.08 0.25 -44 3]);
    ERQmeanint{i} = interp1(Vth1{i},ERQmean{i},x{i}); % interpolate
    [~,the(i)] = (max(1000*diff(diff(f(p(i,:),x{i})))));
    [~,thr(i)] = (min(1000*diff(diff(f(p(i,:),x{i})))));
    Vthe(i)=x{i}(the(i));
    Vthr(i)=x{i}(thr(i));
    ERQthe(i)=ERQmeanint{i}(the(i));
    ERQthr(i)=ERQmeanint{i}(thr(i));
end

meanERQthe=mean(ERQthe); stdERQthe=std(ERQthe);
meanERQthr=mean(ERQthr); stdERQthr=std(ERQthr);

col = colormaps.cbrewer('qual','Set1',length(Vth1)); 
subplot(1,2,2) % all the data
for i = 1:N % experiment 
plot(Vth1{i},ERQmean{i},'.','markersize',20,'linewidth',1.6,'color',col(i,:,:)), hold on
plot(x{i},f(p(i,:),x{i}),'linewidth',2,'color',col(i,:,:))
ylim([-0.1 0.2]); xlim([-54 -33.8])
end
axis square
text(-53,0.17,['N=',num2str(N)],'Fontsize',16,'FontName','Arial')
yticks([-0.1:0.05:0.2]); xticks([-54:2:-34])
ylabel('ERQ'); xlabel('Synaptic threshold (V_{th}), mV')
set(gca,'Fontsize',16,'FontName','Arial')

%% interpolate all the measures so that they are within the same boundaries (-54 to -28)

Vthint=cell(1,2); FRint=cell(1,2); Aint=cell(1,2);
DCint=cell(1,2);  SpkFRint=cell(1,2); Nspksint=cell(1,2); ERQint=cell(1,2);

for n=1:2 % neuron number
    for i = 1:N % experiment
        Vthall =[-54:2:-28];
        Vthint{n}(i,:) = interp1(Vth1{i},Vth1{i},Vthall,'linear','extrap'); % threshold
        FRint{n}(i,:) = interp1(Vth1{i},FR{i}{n},Vthall); % cycle frequency
        Aint{n}(i,:) = interp1(Vth1{i},A{i}{n},Vthall); % amplitude
        DCint{n}(i,:) = interp1(Vth1{i},DC{i}{n},Vthall); % duty cycle
        SpkFRint{n}(i,:) = interp1(Vth1{i},SpkFR{i}{n},Vthall); % spike frequency
        Nspksint{n}(i,:) = interp1(Vth1{i},Nspks{i}{n},Vthall); % # spikes/burst
        ERQint{n}(i,:) = interp1(Vth1{i},ERQ{i}{n},Vthall); % # ERQ
    end
end

%% calculate interpolated values across two neurons

ERQmeanint=[];

for j = 1:N
        FRmeanint(j,:)=nanmean([FRint{1}(j,:);FRint{2}(j,:)]);
        Ameanint(j,:)=nanmean([Aint{1}(j,:);Aint{2}(j,:)]);
        SpkFRmeanint(j,:)=nanmean([SpkFRint{1}(j,:);SpkFRint{2}(j,:)]);
        Nspksmeanint(j,:)=nanmean([Nspksint{1}(j,:);Nspksint{2}(j,:)]);
        DCmeanint(j,:)=nanmean([DCint{1}(j,:);DCint{2}(j,:)]);
        ERQmeanint(j,:)=nanmean([ERQint{1}(j,:);ERQint{2}(j,:)]);
end

% if the circuit not a half-center set measures to NaN
Ameanint(isnan(FRmeanint)) = NaN;
Nspksmeanint(isnan(FRmeanint)) = NaN;
SpkFRmeanint(isnan(FRmeanint)) = NaN;
DCmeanint(isnan(FRmeanint)) = NaN;
ERQmeanint(isnan(FRmeanint)) = NaN;
    
% if values for one of the neurons are NaN
FRmeanint(isnan(FRint{1}) | isnan(FRint{2}))=NaN;
Ameanint(isnan(Aint{1}) | isnan(Aint{2}))=NaN;
Nspksmeanint(isnan(Nspksint{1}) | isnan(Nspksint{2}))=NaN;
SpkFRmeanint(isnan(SpkFRint{1}) | isnan(SpkFRint{2}))=NaN;
DCmeanint(isnan(DCint{1}) | isnan(DCint{2}))=NaN;
ERQmeanint(isnan(ERQint{1}) | isnan(ERQint{2}))=NaN;

%% interpolate characteristics for ERQ values

%FRint1=[]; ERQint1=[]; Aint1=[]; DCint1=[]; SpkFRint1=[]; Nspksint1=[]; Vthint1=[];

ERQint1 = linspace(min(ERQmeanint(:)),max(ERQmeanint(:)),50);

for i = 1:N 
    x=ERQmeanint(i,find(~isnan(ERQmeanint(i,:))));    
    FRint1(i,:) = interp1(x,FRmeanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1,'linear');
    Aint1(i,:) = interp1(x,Ameanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1);
    DCint1(i,:) = interp1(x,DCmeanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1);
    SpkFRint1(i,:) = interp1(x,SpkFRmeanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1);
    Nspksint1(i,:) = interp1(x,Nspksmeanint(i,find(~isnan(ERQmeanint(i,:)))),ERQint1);
end

%% plot all the characteristics of half-center for all the experiments
 
clf,
subplot(1,5,1) % cycle frequency
plot(ERQmeanint',FRmeanint','linewidth',1,'color','k'), hold on
plot(ERQint1,nanmean(FRint1),'r','linewidth',2)
ylim([0.05 0.65]); xlim([-0.1 0.2])
xticks([-0.1:0.05:0.2]); yticks([0:0.1:0.7])
text(-0.08,0.6,['N=',num2str(N)],'Fontsize',16,'FontName','Arial')
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

%% Pearson correlation coefficient between the cycle frequency and the amplitude

r = corrcoef(FRmeanint(~isnan(Ameanint)),Ameanint(~isnan(Ameanint)))
%clf ,scatter(FRmeanint(~isnan(Ameanint)),Ameanint(~isnan(Ameanint)),'filled')

%% Spearman rank correlation test (duty cycle vs ERQ and spike frequency vs ERQ)

[rhoDC,pval]=corr(nanmean(DCint1)',ERQint1','Type','Spearman') 
[rhoSpkFR,pval]=corr(nanmean(SpkFRint1)',ERQint1','Type','Spearman') 
