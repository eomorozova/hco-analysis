%% load data from temperature experiments

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_temp_pooled.mat')

%% calculate the charackteristics of the half-center output at 10C and 20C

for i=1:numel(data.V)
    i
    for j=1:2 % temperature: 10,20C
        for ii=1:2 % neuron
            ERQ_{i}{j}(ii) = analysis.erq(data.V{i}{j}(ii,:),data.Vth(i));
            [hcostat{i}{j}{ii}] = analysis.hco_stat(data.V{i}{j}(ii,:), data.Fs(i));
            
            FR(i,j,ii) = 1./hcostat{i}{j}{ii}.T1_mean;
            dc(i,j,ii) = hcostat{i}{j}{ii}.dc_mean;
            A(i,j,ii) = hcostat{i}{j}{ii}.A;
            SpkFR(i,j,ii) = hcostat{i}{j}{ii}.SpkFreq_mean;
            Nspks(i,j,ii) = hcostat{i}{j}{ii}.nSpks_mean;
            ERQ(i,j,ii) = ERQ_{i}{j}(ii);
            
        end
        meanFR(i,j) = mean(FR(i,j,:)); % average cycle frequency based on calculations from both neurons
    end
    dFR(i) = diff(meanFR(i,:)); % change in cycle frequency with temperature
end

%% calculate mean change of each characteristic of half-center output 

ie=find(data.mode=="escape");
ir=find(data.mode=="release");
iq1=find(data.Q10==11);
iq2g=find(data.Q10==21);
iq2=find(data.Q10==22);

idxeq1=intersect(ie,iq1); idxeq2g=intersect(ie,iq2g); idxeq2=intersect(ie,iq2);  
idxrq1=intersect(ir,iq1); idxrq2g=intersect(ir,iq2g); idxrq2=intersect(ir,iq2);  

prop=[FR, SpkFR, Nspks, A, dc, ERQ]; % all the properties into a single array

propall = [prop(idxeq1,:,1);prop(idxeq1,:,2);prop(idxeq2g,:,1);prop(idxeq2g,:,2);...
    prop(idxeq2,:,1);prop(idxeq2,:,2); prop(idxrq1,:,1);prop(idxrq1,:,2);
    prop(idxrq2g,:,1);prop(idxrq2g,:,2); prop(idxrq2,:,1);prop(idxrq2,:,2)];

%dpropall=[];

for i=1:size(prop,2)/2
    dpropall(i,:) = diff(propall(:,(i-1)*2+1:2*i)');
end

propall =  propall([1:21,23:end],:); % exclude LG
dpropall =  dpropall(:,[1:21,23:end]);

% mean and standart deviation for each neuron

idxeq1all = 1:idxeq1(end)*2-1;
idxeq2gall = idxeq1all(end)+1:idxeq1all(end)+numel(idxeq2g)*2;
idxeq2all = idxeq2gall(end)+1:idxeq2gall(end)+numel(idxeq2)*2;
idxrq1all = idxeq2all(end)+1:idxeq2all(end)+numel(idxrq1)*2;
idxrq2gall = idxrq1all(end)+1:idxrq1all(end)+numel(idxrq2g)*2;
idxrq2all = idxrq2gall(end)+1:idxrq2gall(end)+numel(idxrq2)*2;

for i=1:size(dpropall,1)
    
    meandpropall(i,:)=[nanmean(dpropall(i,idxeq1all)), nanmean(dpropall(i,idxeq2gall)),...
        nanmean(dpropall(i,idxeq2all)), nanmean(dpropall(i,idxrq1all)),...
        nanmean(dpropall(i,idxrq2gall)),nanmean(dpropall(i,idxrq2all))];
    
    errdpropall(i,:)=[nanstd(dpropall(i,idxeq1all)), nanstd(dpropall(i,idxeq2gall)),...
        nanstd(dpropall(i,idxeq2all)), nanstd(dpropall(i,idxrq1all)),...
        nanstd(dpropall(i,idxrq2gall)),nanstd(dpropall(i,idxrq2all))];
end

groupIdx = [ones(size(idxeq1all))'; 2*ones(size(idxeq2gall))';...
    3*ones(size(idxeq2all))'; 4*ones(size(idxrq1all))';...
    5*ones(size(idxrq2gall))'; 6*ones(size(idxrq2all))'];

groupIdx1 = [ones(size(idxeq1))'; 2*ones(size(idxeq2g))';...
    3*ones(size(idxeq2))'; 4*ones(size(idxrq1))';...
    5*ones(size(idxrq2g))'; 6*ones(size(idxrq2))'];

%% plot absolute change in half-center characteristics with temperature
% (Figue 7 E-H + chnge in duty cycle and ERQ)

clf
color = {[0.36,0.42,0.6],[0.6,0.4,1],[0.36,0.23,0.6],[0.8, 0.6, 0],[0.8, 0.3, 0],[0.5, 0, 0]};
symbol = {'o','o','o','s','s','s'};

for i = 1:size(dpropall,1)
    if i>1
        display.bigsubplot(1,size(dpropall,1),1,i, 0.05, 0.05)
        display.plot_barerror(meandpropall(i,:),errdpropall(i,:),'color',color); hold on
        display.plot_scatter(dpropall(i,:), groupIdx,'color',{'k','k','k','k','k','k'},'symbol',symbol,'markersize',30);
        xticks([1,2,3,4,5,6]);
        xticklabels({'q_{10}=1','q_{10}=2 gs','q_{10}=2','q_{10}=1','q_{10}=2 gs','q_{10}=2'});
    end
    box off; axis square
    if i==1
        display.bigsubplot(1,size(dpropall,1),1,i, 0.05, 0.05)
        display.plot_barerror(meandpropall(i,:),errdpropall(i,:),'color',color); hold on
        display.plot_scatter(dFR,groupIdx1,'color',{'k','k','k','k','k','k'},'symbol',symbol,'markersize',30);
        xticks([1,2,3,4,5,6]);
        xticklabels({'q_{10}=1','q_{10}=2 gs','q_{10}=2','q_{10}=1','q_{10}=2 gs','q_{10}=2'});
        title('% change in cycle frequency')
        ylabel('% \Delta in cycle freq, Hz');
        ylim([-0.15 0.45])
        box off; axis square
 
    elseif i==2
        title('Change in spike frequency')
        ylabel('\Delta in spike freq, Hz');
        ylim([-1 20])
    elseif i==3
        title('Change in # spikes/burst')
        ylabel('\Delta in # of spikes/burst');
        ylim([-20 50])
    elseif i==4
        title('Change in amplitude')
        ylabel('\Delta in amplitude, mV');
        ylim([-15 20])
    elseif i==5
        title('Change in duty cycle')
        ylabel('\Delta in duty cycle, %');
        ylim([-20 20])
    elseif i==6
        title('Change in ERQ')
        ylabel('\Delta in ERQ');
        ylim([-0.08 0.08])
    end

    xlim([0.5 6.5])    
    set(gca,'Fontsize',14,'Fontname','Arial')
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
end

%% stats (compare each property of half-center output @ 10C and 20C)

% Wilcoxon signed rank test (non parametric for paired samples)

% also performed Wilcoxon signed rank test in SPSS 24, p-values are slightly different
% but significant pairs are the same in both SPSS and MATLAB

p_eq1(1) = signrank(meanFR(idxeq1,1),meanFR(idxeq1,2));
p_eq2g(1) = signrank(meanFR(idxeq2g,1),meanFR(idxeq2g,2));
p_eq2(1) = signrank(meanFR(idxeq2,1),meanFR(idxeq2,2));

p_rq1(1) = signrank(meanFR(idxrq1,1),meanFR(idxrq1,2));
p_rq2g(1) = signrank(meanFR(idxrq2g,1),meanFR(idxrq2g,2));
p_rq2(1) = signrank(meanFR(idxrq2,1),meanFR(idxrq2,2));

for i=2:size(propall,2)/2
    
    meanprop_eq1(i,:) = [mean(propall(idxeq1all,(i-1)*2+1)), mean(propall(idxeq1all,i*2))];
    stdprop_eq1(i,:) = [std(propall(idxeq1all,(i-1)*2+1)), std(propall(idxeq1all,i*2))];
    
    meanprop_eq2g(i,:) = [mean(propall(idxeq2gall,(i-1)*2+1)), mean(propall(idxeq2gall,i*2))];
    stdprop_eq2g(i,:) = [std(propall(idxeq2gall,(i-1)*2+1)), std(propall(idxeq2gall,i*2))];
    
    meanprop_eq2(i,:) = [mean(propall(idxeq2all,(i-1)*2+1)), mean(propall(idxeq2all,i*2))];
    stdprop_eq2(i,:) = [std(propall(idxeq2all,(i-1)*2+1)), std(propall(idxeq2all,i*2))];
    
    meanprop_rq1(i,:) = [mean(propall(idxrq1all,(i-1)*2+1)), mean(propall(idxrq1all,i*2))];
    stdprop_rq1(i,:) = [std(propall(idxrq1all,(i-1)*2+1)), std(propall(idxrq1all,i*2))];
    
    meanprop_rq2g(i,:) = [mean(propall(idxrq2gall,(i-1)*2+1)), mean(propall(idxrq2gall,i*2))];
    stdprop_rq2g(i,:) = [std(propall(idxrq2gall,(i-1)*2+1)), std(propall(idxrq2gall,i*2))];
    
    meanprop_rq2(i,:) = [mean(propall(idxrq2all,(i-1)*2+1)), mean(propall(idxrq2all,i*2))];
    stdprop_rq2(i,:) = [std(propall(idxrq2all,(i-1)*2+1)), std(propall(idxrq2all,i*2))];    
    
    p_eq1(i) = signrank(propall(idxeq1all,(i-1)*2+1),propall(idxeq1all,i*2));
    p_eq2g(i) = signrank(propall(idxeq2gall,(i-1)*2+1),propall(idxeq2gall,i*2));
    p_eq2(i) = signrank(propall(idxeq2all,(i-1)*2+1),propall(idxeq2all,i*2));
    
    p_rq1(i) = signrank(propall(idxrq1all,(i-1)*2+1),propall(idxrq1all,i*2));
    p_rq2g(i) = signrank(propall(idxrq2gall,(i-1)*2+1),propall(idxrq2gall,i*2));
    p_rq2(i) = signrank(propall(idxrq2all,(i-1)*2+1),propall(idxrq2all,i*2));
    
end
