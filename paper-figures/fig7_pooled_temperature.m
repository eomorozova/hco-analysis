%%
load('C:\Users\moroz\Documents\code\hco_analysis\data\data_temp_pooled.mat')

%%
Fs=10000; % Sampling frequency
ERQ=[]; hcostat=[];
FR=[]; dc=[]; A=[]; SpkFR=[]; Nspks=[]; ERQ1=[];
meanFR=[]; meandc=[]; meanA=[]; meanSpkFR=[]; meanNspks=[]; meanERQ=[];
dperFR=[]; dperdc=[]; dperA=[]; dperSpkFR=[]; dperNspks=[]; dperERQ1=[];

k=1;
for i=1:numel(data.V)
    i
    for j=1:2 % temperature: 10,20 C
        for ii=1:2 % neuron
            ERQ{i}{j}(ii) = analysis.erq(data.V{i}{j}(ii,:),data.Vth(i));
            [hcostat{i}{j}{ii}] = analysis.hco_stat(data.V{i}{j}(ii,:), Fs);
            
            FR(k,j,ii) = 1./hcostat{i}{j}{ii}.T1_mean;
            dc(k,j,ii) = hcostat{i}{j}{ii}.dc_mean;
            A(k,j,ii) = hcostat{i}{j}{ii}.A;
            SpkFR(k,j,ii) = hcostat{i}{j}{ii}.SpkFreq_mean;
            Nspks(k,j,ii) = hcostat{i}{j}{ii}.nSpks_mean;
            ERQ1(k,j,ii) = ERQ{i}{j}(ii);
        end
        
        % mean characteristics across two neurons
        if i~= 11
            meanFR(k,j) = mean(FR(k,j,:),3);
            meandc(k,j) = mean(dc(k,j,:),3);
            meanA(k,j) = mean(A(k,j,:),3);
            meanSpkFR(k,j) = mean(SpkFR(k,j,:),3);
            meanNspks(k,j) = mean(Nspks(k,j,:),3);
            meanERQ(k,j) = mean(ERQ1(k,j,:),3);
        else  % 2nd neuron in 11 is LG (exclude)
            meanFR(k,j) = FR(k,j,1);
            meandc(k,j) = dc(k,j,1);
            meanA(k,j) = A(k,j,1);
            meanSpkFR(k,j) = SpkFR(k,j,1);
            meanNspks(k,j) = Nspks(k,j,1);
            meanERQ(k,j) = ERQ1(k,j,1);
        end
        
    end
    
    % change in characteristics with temperature
    
    dFR(k) = diff(meanFR(k,:));
    ddc(k) = diff(meandc(k,:));
    dA(k) = diff(meanA(k,:));
    dSpkFR(k) = diff(meanSpkFR(k,:));
    dNspks(k) = diff(meanNspks(k,:));
    dERQ1(k) = diff(meanERQ(k,:));
    
    % percent change in characteristics with temperature
    
    dperFR(k) = dFR(k)/abs(meanFR(k,1))*100;
    dperdc(k) = ddc(k); % already in %
    dperA(k) = dA(k)/abs(meanA(k,1))*100;
    dperSpkFR(k) = dSpkFR(k)/abs(meanSpkFR(k,1))*100;
    dperNspks(k) = dNspks(k)/abs(meanNspks(k,1))*100;
    dperERQ1(k) = dERQ1(k)/abs(meanERQ(k,1))*100;
    
    k=k+1;
end


%% calculate mean and std for each HCO property

ie=find(data.mode=="escape");
ir=find(data.mode=="release");
iq1=find(data.Q10==11);
iq2g=find(data.Q10==21);
iq2=find(data.Q10==22);

dperprop=[dperFR; dperSpkFR; dperA; dperNspks; dperdc; dperERQ1];

idxeq1=intersect(ie,iq1); idxeq2g=intersect(ie,iq2g); idxeq2=intersect(ie,iq2);  
idxrq1=intersect(ir,iq1); idxrq2g=intersect(ir,iq2g); idxrq2=intersect(ir,iq2);  

meanprop=[]; errorprop=[];

for i=1:size(dperprop,1)
    meanprop(i,:)=[mean(dperprop(i,idxeq1)), mean(dperprop(i,idxeq2g)),...
        mean(dperprop(i,idxeq2)), mean(dperprop(i,idxrq1)),...
        mean(dperprop(i,idxrq2g)),mean(dperprop(i,idxrq2)),];
    
    errprop(i,:)=[std(dperprop(i,idxeq1)), std(dperprop(i,idxeq2g)),...
        std(dperprop(i,idxeq2)), std(dperprop(i,idxrq1)),...
        std(dperprop(i,idxrq2g)),std(dperprop(i,idxrq2)),];    
end

groupIdx = [ones(size(idxeq1))'; 2*ones(size(idxeq2g))';...
    3*ones(size(idxeq2))'; 4*ones(size(idxrq1))';...
    5*ones(size(idxrq2g))'; 6*ones(size(idxrq2))'];

%% anova to compare percent change in hco propperties across conditions

stats=[]; pvalsig=[]; sigpoints=[];

for i=1%:size(dperprop,1)
    
    [p(i),tbl{i},stats{i}]=anova1(dperprop(i,:)',groupIdx');
    comp=multcompare(stats{i})
    group=(comp(:,1:2));
    hcomp=comp(:,6)<0.05;
    pvalsig{i} = comp(hcomp,6);
    sigpoints{i}=([group(hcomp==1,1), group(hcomp==1,2)])
    groups{i}=(num2cell(sigpoints{i},2));
    
end

%% plot percent change in all the measures
clf

color = {[0.36,0.42,0.6],[0.6,0.4,1],[0.36,0.23,0.6],[0.8, 0.6, 0],[0.8, 0.3, 0],[0.5, 0, 0]};
symbol = {'o','o','o','s','s','s'};

for i = 1:3 %size(dperall,1)
    %bigsubplot(1,size(dperall,1),1,i, 0.03, 0.05)
    bigsubplot(2,3,1,i, 0.03, 0.06)
    display.plot_barerror(meanprop(i,:),errprop(i,:),'color',color); hold on
    display.plot_scatter(dperprop(i,[1:33]), groupIdx([1:33]),'color',{'k','k','k','k','k','k'},'symbol',symbol,'markersize',30);
    xticks([1,2,3,4,5,6]);
    xticklabels({'q_{10}=1','q_{10}=2 gs','q_{10}=2','q_{10}=1','q_{10}=2 gs','q_{10}=2'});
    
    box off; axis square
    if i==1
        title('% change in cycle frequency')
        ylabel('% \Delta in cycle freq');
        ylim([-60 160])
        % uncomment to plot asterisks
        %bar(meanprop(i,:),0.7,'FaceALpha',0); hold on
        %display.sigstar(groups{i}',pvalsig{i}')
    elseif i==2
        title('% change in spike frequency')
        ylabel('% \Delta in spike freq');
        ylim([-10 200])
        %bar(meanprop(i,:),0.7,'FaceALpha',0); hold on
        %display.sigstar(groups{i}',pvalsig{i}')
    elseif i==3
        title('% change in amplitude')
        ylabel('% \Delta in amplitude');
        ylim([-40 90])
    elseif i==4
        title('% change in # spikes/burst')
        ylabel('% \Delta in # of spikes/burst');
        ylim([-40 360])
    elseif i==5
        title('% change in duty cycle')
        ylabel('% \Delta in duty cycle');
        ylim([-20 20])
    elseif i==6
        title('% change in ERQ')
        ylabel('% \Delta in ERQ');
        ylim([-40 40])
    end
    
    set(gca,'Fontsize',14,'Fontname','Arial')
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
end

for i = 4:size(dperprop,1)
    %bigsubplot(1,size(dperall,1),1,i, 0.03, 0.05)
    bigsubplot(2,3,2,i-3, 0.03, 0.06)
    display.plot_barerror(meanprop(i,:),errprop(i,:),'color',color); hold on
    display.plot_scatter(dperprop(i,[1:33]), groupIdx([1:33]),'color',{'k','k','k','k','k','k'},'symbol',symbol,'markersize',30);
    xticks([1,2,3,4,5,6]);
    xticklabels({'q_{10}=1','q_{10}=2 gs','q_{10}=2','q_{10}=1','q_{10}=2 gs','q_{10}=2'});
    
    box off; axis square
    if i==1
        title('% change in cycle frequency')
        ylabel('% \Delta in cycle freq');
        ylim([-60 160])
    elseif i==2
        title('% change in spike frequency')
        ylabel('% \Delta in spike freq');
        ylim([-10 200])
    elseif i==3
        title('% change in amplitude')
        ylabel('% \Delta in amplitude');
        ylim([-40 90])
    elseif i==4
        title('% change in # spikes/burst')
        ylabel('% \Delta in # of spikes/burst');
        ylim([-55 360])
    elseif i==5
        title('% change in duty cycle')
        ylabel('% \Delta in duty cycle');
        ylim([-20 20])
    elseif i==6
        title('% change in ERQ')
        ylabel('% \Delta in ERQ');
        ylim([-40 40])
    end
    
    set(gca,'Fontsize',14,'Fontname','Arial')
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
end

%% stats (compare each hco property @ 10C and 20C)
% Wilcoxon signed rank test (non parametric for paired samples)

FRall = permute(FR,[1 3 2]); FRall = reshape(FRall,[],size(A,2),1);
SpkFRall = permute(SpkFR,[1 3 2]); SpkFRall = reshape(SpkFRall,[],size(A,2),1);
Nspksall = permute(Nspks,[1 3 2]); Nspksall = reshape(Nspksall,[],size(A,2),1);
Aall = permute(A,[1 3 2]); Aall = reshape(Aall,[],size(A,2),1);
dcall = permute(dc,[1 3 2]); dcall = reshape(dcall,[],size(A,2),1);
ERQall = permute(ERQ1,[1 3 2]); ERQall = reshape(ERQall,[],size(A,2),1);

prop=[FRall, SpkFRall, Aall, Nspksall, dcall, ERQall];

for i=1:size(prop,2)/2
    
    p_eq1(i) = signrank(prop(idxeq1(1):idxeq1(end)+numel(idxeq1),(i-1)*2+1),prop(idxeq1(1):idxeq1(end)+numel(idxeq1),i*2));
    p_eq2g(i) = signrank(prop(idxeq2g(1):idxeq2g(end)+numel(idxeq2g),(i-1)*2+1),prop(idxeq2g(1):idxeq2g(end)+numel(idxeq2g),i*2));
    p_eq2(i) = signrank(prop(idxeq2(1):idxeq2(end)+numel(idxeq2),(i-1)*2+1),prop(idxeq2(1):idxeq2(end)+numel(idxeq2),i*2));
    
    p_rq1(i) = signrank(prop(idxrq1(1):idxrq1(end)+numel(idxrq1),(i-1)*2+1),prop(idxrq1(1):idxrq1(end)+numel(idxrq1),i*2));
    p_rq2g(i) = signrank(prop(idxrq2g(1):idxrq2g(end)+numel(idxrq2g),(i-1)*2+1),prop(idxrq2g(1):idxrq2g(end)+numel(idxrq2g),i*2));
    p_rq2(i) = signrank(prop(idxrq2(1):idxrq2(end)+numel(idxrq2),(i-1)*2+1),prop(idxrq2(1):idxrq2(end)+numel(idxrq2),i*2));
    
end
