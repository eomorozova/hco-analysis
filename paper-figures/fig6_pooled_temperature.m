%%
load('C:\Users\moroz\Documents\code\hco_analysis\data\data_temp_pooled.mat')

%%
Fs=10000; % Sampling frequency
ERQ=[]; hcostat=[];
meanFR=[]; meandc=[]; meanA=[]; meanSpkFR=[]; meanNspks=[]; meanERQ1=[];
k=1;
for i=[ie,ir]
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
        
        meanFR(k,j) = mean(FR(k,j,:),3);
        meandc(k,j) = mean(dc(k,j,:),3);
        meanA(k,j) = mean(A(k,j,:),3);
        meanSpkFR(k,j) = mean(SpkFR(k,j,:),3);
        meanNspks(k,j) = mean(Nspks(k,j,:),3);
        meanERQ(k,j) = mean(ERQ1(k,j,:),3);
        
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


%%
ie=find(data.mode=="escape");
ir=find(data.mode=="release");
iq1=find(data.Q10==11);
iq2g=find(data.Q10==21);
iq2=find(data.Q10==22);
%% calculate mean and std for each HCO property

dperprop=[dperFR; dperSpkFR; dperA; dperNspks; dperdc; dperERQ1];

meanprop=[]; errorprop=[];

for i=1:size(dperall,1)
    meanprop(i,:)=[mean(dperprop(i,intersect(ie,iq1))), mean(dperprop(i,intersect(ie,iq2g))),...
        mean(dperprop(i,intersect(ie,iq2))), mean(dperprop(i,intersect(ir,iq1))),...
        mean(dperprop(i,intersect(ir,iq2g))),mean(dperprop(i,intersect(ir,iq2))),];
    
    errprop(i,:)=[std(dperprop(i,intersect(ie,iq1))), std(dperprop(i,intersect(ie,iq2g))),...
        std(dperprop(i,intersect(ie,iq2))), std(dperprop(i,intersect(ir,iq1))),...
        std(dperprop(i,intersect(ir,iq2g))),std(dperprop(i,intersect(ir,iq2))),];    
end

groupIdx = [ones(size(intersect(ie,iq1)))'; 2*ones(size(intersect(ie,iq2g)))';...
    3*ones(size(intersect(ie,iq2)))'; 4*ones(size(intersect(ir,iq1)))';...
    5*ones(size(intersect(ir,iq2g)))'; 6*ones(size(intersect(ir,iq2)))'];
%% plot perent change in all the measures
clf
color = {[0.36,0.42,0.6],[0.6,0.4,1],[0.36,0.23,0.6],[0.8, 0.6, 0],[0.8, 0.3, 0],[0.5, 0, 0]};

for i = 1:size(dperall,1)
    subplot(1,size(dperall,1),i)
    display.plot_barerror(meanprop(i,:),errprop(i,:),'color',color); hold on
    display.plot_scatter(dperprop(i,1:5), groupIdx(1:5),'color',{'k','k','k','k','k','k'},'symbol',symbol,'markersize',50)
    
    xticks([1,2,3,4,5,6]);
    xticklabels({'q_{10}=1','q_{10}=2 gs','q_{10}=2','q_{10}=1','q_{10}=2 gs','q_{10}=2'});
    
    box off; %axis square
    if i==1
        title('% change in cycle frequency')
        ylabel('% \Delta in cycle freq');
    elseif i==2
        title('% change in spike frequency')
        ylabel('% \Delta in spike freq');
    elseif i==3
        title('% change in amplitude')
        ylabel('% \Delta in amplitude');
    elseif i==4
        title('% change in # spikes/burst')
        ylabel('% \Delta in # of spikes/burst');
    elseif i==5
        title('% change in duty cycle')
        ylabel('% \Delta in duty cycle');
    elseif i==6
        title('% change in ERQ')
        ylabel('% \Delta in ERQ');
    end
    
    set(gca,'Fontsize',14,'Fontname','Arial')
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
end
