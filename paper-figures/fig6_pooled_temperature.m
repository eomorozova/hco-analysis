%%
load('C:\Users\moroz\Documents\code\hco_analysis\data\data_temp_pooled.mat')

%%
Fs=10000; % Sampling frequency
ERQ=[]; hcostat=[];
FR=[]; dc=[]; A=[]; SpkFR=[]; Nspks=[]; ERQ1=[];
k=1;
for i=[ie,ir]
    i
    for j=1:2 % temperature: 10,20 C
        for ii=1:2 % neuron
            ERQ{i}{j}(ii) = analysis.erq(data.V{i}{j}(ii,:),data.Vth(i));
            [hcostat{i}{j}{ii}] = analysis.hco_stat(data.V{i}{j}(ii,:), Fs);
            
        end
        
        % mean characteristics across two neurons
        
        FR(k,j) = mean([1./hcostat{i}{j}{1}.T1_mean,1./hcostat{i}{j}{2}.T1_mean]);
        dc(k,j) = mean([hcostat{i}{j}{1}.dc_mean,hcostat{i}{j}{2}.dc_mean]);
        A(k,j) = mean([hcostat{i}{j}{1}.A,hcostat{i}{j}{2}.A]);
        SpkFR(k,j) = mean([hcostat{i}{j}{1}.SpkFreq_mean,hcostat{i}{j}{2}.SpkFreq_mean]);
        Nspks(k,j) = mean([hcostat{i}{j}{1}.nSpks_mean,hcostat{i}{j}{2}.nSpks_mean]);
        ERQ1(k,j) = mean([ERQ{i}{j}(1),ERQ{i}{j}(2)]);
        
    end

    % change in characteristics with temperature

    dFR(k) = diff(FR(k,:)); 
    ddc(k) = diff(dc(k,:));
    dA(k) = diff(A(k,:));
    dSpkFR(k) = diff(SpkFR(k,:));
    dNspks(k) = diff(Nspks(k,:));
    dERQ1(k) = diff(ERQ1(k,:));

    % percent change in characteristics with temperature

    dperFR(k) = dFR(k)/FR(k,1)*100; 
    dperdc(k) = ddc(k)/dc(k,1)*100;
    dperA(k) = ddc(k)/dc(k,1)*100;
    dperSpkFR(k) = dSpkFR(k)/SpkFR(k,1)*100;
    dperNspks(k) = dNspks(k)/Nspks(k,1)*100;
    dperERQ1(k) = dERQ1(k)/ERQ1(k,1)*100;
    
    k=k+1;
end

%%
ie=find(data.mode=="escape");
ir=find(data.mode=="release");
iq1=find(data.Q10==11);
iq2g=find(data.Q10==21);
iq2=find(data.Q10==22);

%% plot percent change in frequency

color = {[0.36,0.42,0.6],[0.6,0.4,1],[0.36,0.23,0.6],[0.8, 0.6, 0],[0.8, 0.3, 0],[0.5, 0, 0]};
symbol = {'o','o','o','s','s','s'};
clf; subplot(2,4,1)

meanFR=[median(dperFR(intersect(ie,iq1))),median(dperFR(intersect(ie,iq2g))),median(dperFR(intersect(ie,iq2))),...
median(dperFR(intersect(ir,iq1))),median(dperFR(intersect(ir,iq2g))),median(dperFR(intersect(ir,iq2))),];

errFR=[std(dperFR(intersect(ie,iq1))),std(dperFR(intersect(ie,iq2g))),std(dperFR(intersect(ie,iq2))),...
std(dperFR(intersect(ir,iq1))),std(dperFR(intersect(ir,iq2g))),std(dperFR(intersect(ir,iq2))),];

display.plot_barerror(meanFR,errFR,'color',color); hold on

groupIdx = [ones(size(intersect(ie,iq1)))'; 2*ones(size(intersect(ie,iq2g)))';...
 3*ones(size(intersect(ie,iq2)))'; 4*ones(size(intersect(ir,iq1)))';...
 5*ones(size(intersect(ir,iq2g)))'; 6*ones(size(intersect(ir,iq2)))'];

display.plot_scatter(dperFR, groupIdx,'color',{'k','k','k','k','k','k'},'symbol',symbol,'markersize',50)

h=display.plothorzline(0); set(h,'color','k','linewidth',1)

xlim([0.5 6.5]); ylim([-60 180])
xticks([1,2,3,4,5,6]);
xticklabels({'q_{10}=1','q_{10}=2 gs','q_{10}=2','q_{10}=1','q_{10}=2 gs','q_{10}=2'});
ylabel('% \Delta in cycle frequency');
title('Percent change in cycle frequency')
set(gca,'Fontsize',14,'Fontname','Arial')
box off; axis square
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

