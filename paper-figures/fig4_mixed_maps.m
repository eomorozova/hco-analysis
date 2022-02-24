
%% load all the data with the mixed maps

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig4_fig4suppl.mat')

%% calculate half-center output characteristcs as a function of gH and gSyn

hcostat=cell(1,2); ERQ=cell(1,2);

for jj = 1:2 % experiment
    jj
    for j = 1:numel(data.V{jj}) % map number
        for i = 1:length(data.V{jj}{j}) % (gh,gsyn) combination
            for ii = 1:2 % neuron number
                if ~isempty(data.V{jj}{j}{i})==1 & ~isnan(data.V{jj}{j}{i})==1
                    
                    Fs(jj) = data.Fs(jj);
                    V = data.V{jj}{j}{i}(ii,5*Fs(jj):end);                  
                    Vth = data.Vth{jj}(j);
                    
                    ERQ{jj}{j}(ii,i) = analysis.erq(V,Vth);
                    [hcostat{jj}{j}{i}{ii}] = analysis.hco_stat(V, Fs(jj), Vth);
                    
                else
                    hcostat{jj}{j}{i}{ii} = [];  ERQ{jj}{j}(ii,i) = NaN;
                end
                
                FR{jj}{j}(ii,i) = 1./hcostat{jj}{j}{i}{ii}.T1_mean; % cycle frequency
                A{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.A; % slow-wave amplitude
                SpkFR{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.SpkFreq_mean; % spike frequency
                Nspk{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.nSpks_mean; % # spikes/burst
                DC{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.dc_mean; % #duty cycle (based on spikes)
                dc{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.DC_mean; % #duty cycle (based on slow-wave)
                
                stdFR{jj}{j}(ii,i) = nanstd(1./hcostat{jj}{j}{i}{ii}.T1); % standatd deviation of cycle frequency
                CVFR{jj}{j}(ii,i) = stdFR{jj}{j}(ii,i)/FR{jj}{j}(ii,i); % coefficient of variation of cycle frequency
                
                BurDur{jj}{j}{ii,i} = hcostat{jj}{j}{i}{ii}.BurDur; % burst durations based on slow-wave
                meanBurDur{jj}{j}(ii,i) =  mean(BurDur{jj}{j}{ii,i}); % mean burst duration
                
            end
        end
        
        % mean characteristics across two neurons
        FR_{jj}(j,:) = mean(FR{jj}{j});
        dc_{jj}(j,:) = mean(dc{jj}{j});
        DC_{jj}(j,:) = mean(DC{jj}{j});
        A_{jj}(j,:) = mean(A{jj}{j});
        SpkFR_{jj}(j,:) = mean(SpkFR{jj}{j});
        Nspks_{jj}(j,:) = mean(Nspk{jj}{j});
        ERQ_{jj}(j,:) = mean(ERQ{jj}{j});
        CVFR_{jj}(j,:) = mean(CVFR{jj}{j});
        
        %reshape
        ERQmean{jj,j} = reshape(ERQ_{jj}(j,:),[7,7]);
        ERQ1{jj,j}{1} = reshape(ERQ{jj}{j}(1,:),[7,7]);
        ERQ1{jj,j}{2} = reshape(ERQ{jj}{j}(2,:),[7,7]);
        
        FRmean{jj,j} = reshape(FR_{jj}(j,:),[7,7]);
        Amean{jj,j} = reshape(A_{jj}(j,:),[7,7]);
        SpkFRmean{jj,j} = reshape(SpkFR_{jj}(j,:),[7,7]);
        Nspkmean{jj,j} = reshape(Nspks_{jj}(j,:),[7,7]);
        DCmean{jj,j} = reshape(DC_{jj}(j,:),[7,7]);
        CVFRmean{jj,j} = reshape(CVFR_{jj}(j,:),[7,7]);
        
        BrstDur{jj,j}{1} = reshape(meanBurDur{jj}{j}(1,:),[7,7]);
        BrstDur{jj,j}{2} = reshape(meanBurDur{jj}{j}(2,:),[7,7]);
        asymm{jj,j} = abs(BrstDur{jj,j}{1}-BrstDur{jj,j}{2});
        
    end
end

%% classify the network state (silent/asymmetric/irregular/antiphase spiking/half-center)

state=cell(1,2);

for jj = 1:2 % experiment
    for j = 1:numel(data.V{jj}) % map
        for i = 1:numel(hcostat{jj}{j}) % (gh,gsyn) combination
            networkstate{jj}{j}(i)=analysis.classify_networkstate(hcostat{jj}{j}{i},Fs(jj));
        end
        state{jj}{j} = reshape(networkstate{jj}{j},[7,7]);
    end
end

%% if identified state is not a half-center, set output characteristic to NaN

for jj=1:2
    for j=1:5
        ERQmean{jj,j}(state{jj}{j}~=5)=NaN;
        ERQ1{jj,j}{1}(state{jj}{j}~=5)=NaN;
        ERQ1{jj,j}{2}(state{jj}{j}~=5)=NaN;
        
        FRmean{jj,j}(state{jj}{j}~=5)=NaN;
        Amean{jj,j}(state{jj}{j}~=5)=NaN;
        SpkFRmean{jj,j}(state{jj}{j}~=5)=NaN;
        Nspkmean{jj,j}(state{jj}{j}~=5)=NaN;
        DCmean{jj,j}(state{jj}{j}~=5)=NaN;
        
        CVFRmean{jj,j}(state{jj}{j}~=5)=NaN;
        asymm{jj,j}(state{jj}{j}~=5)=NaN;
    end
end

%% min and max limits for the maps

for jj=1:2
    for j=1:5
        minERQ_(jj,j) = min(ERQmean{jj,j}(:)); maxERQ_(jj,j) = max(ERQmean{jj,j}(:));
        minFR_(jj,j) = min(FRmean{jj,j}(:)); maxFR_(jj,j) = max(FRmean{jj,j}(:));
        minA_(jj,j) = min(Amean{jj,j}(:)); maxA_(jj,j) = max(Amean{jj,j}(:));
        minNspk_(jj,j) = min(Nspkmean{jj,j}(:)); maxNspk_(jj,j) = max(Nspkmean{jj,j}(:));
        minSpkFR_(jj,j) = min(SpkFRmean{jj,j}(:)); maxSpkFR_(jj,j) = max(SpkFRmean{jj,j}(:));
        minDC_(jj,j) = min(DCmean{jj,j}(:)); maxDC_(jj,j) = max(DCmean{jj,j}(:));
        %minCVFR_(jj,j) = min(CVFRmean{jj,j}(:)); maxCVFR_(jj,j) = max(CVFRmean{jj,j}(:));
        minCVFR_(jj,j) = min(CVFRmean{jj,j}(:)); maxCVFR_(jj,j) = max(CVFRmean{jj,j}(CVFRmean{jj,j}(:)<0.3));
        minasymm_(jj,j) = min(asymm{jj,j}(:)); maxasymm_(jj,j) = max(asymm{jj,j}(:));
    end
end

minERQ = min(minERQ_(:)); maxERQ = max(maxERQ_(:));
minFR = min(minFR_(:)); maxFR = max(maxFR_(:));
minA = min(minA_(:)); maxA = max(maxA_(:));
minNspk = min(minNspk_(:)); maxNspk = max(maxNspk_(:));
minSpkFR = min(minSpkFR_(:)); maxSpkFR = max(maxSpkFR_(:));
minDC = min(minDC_(:)); maxDC = max(maxDC_(:));
minCVFR = min(minCVFR_(:)); maxCVFR = max(maxCVFR_(:));
minasymm = min(minasymm_(:)); maxasymm = max(maxasymm_(:));

%% determine the mechanism of oscillation (escape, release or mixed)
% based on the ERQ threhold of -0.038 for escape and 0.105 for release

regime=[]; regime1=[]; regime2=[];
regime_=NaN(7,7); regime_1=NaN(7,7); regime_2=NaN(7,7);

for i = 1:2 % experiment
    for j=1:numel(data.V{i}) % map
        % for the mean across two neurons
        regime_(ERQmean{i,j}<-0.038) = 1; % escape
        regime_(ERQmean{i,j}>-0.038 & ERQmean{i,j}<0.105) = 2; % mixed
        regime_(ERQmean{i,j}>0.105) = 3; % release
        regime_(isnan(ERQmean{i,j})) = 0;
        regime{i}{j}=regime_;
        
        % for the first neuron
        regime_1(ERQ1{i,j}{1}<-0.038) = 1; % escape
        regime_1(ERQ1{i,j}{1}>-0.038 & ERQ1{i,j}{1}<0.105) = 2; % mixed
        regime_1(ERQ1{i,j}{1}>0.105) = 3; % release
        regime_1(isnan(ERQ1{i,j}{1})) = 0;
        regime1{i}{j}=regime_1;
        
        % for the second neuron
        regime_2(ERQ1{i,j}{2}<-0.038) = 1; % escape
        regime_2(ERQ1{i,j}{2}>-0.038 & ERQ1{i,j}{2}<0.105) = 2; % mixed
        regime_2(ERQ1{i,j}{2}>0.105) = 3; % release
        regime_2(isnan(ERQ1{jj,j}{2})) = 0;
        regime2{i}{j}=regime_2;

    end
end

%% reproduce Figure 4A (ERQ map)

clf
jj=2; j=5;

gH = data.gH{jj}{1}; gSyn = data.gSyn{jj}{1};

display.genimagesc(ERQmean{jj,j}',gH,gSyn);
myColorMap = display.linspecer(128); myColorMap(1,:) = 1; colormap(myColorMap);
if j ==5; h=colorbar; ylabel(h,'ERQ'); end
caxis([minERQ-0.02 maxERQ+0.02])
title(['Vth=',num2str(data.Vth{jj}(j))])

xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
xlabel('gSyn, nS'); ylabel('gH, nS'); axis square
set(gca,'Fontsize',14,'FontName','Arial');

for i=1:10 % create grid
    line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
    line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
end

%% Figure 4B (example traces)

jj=2; j=5;
ii=[43,30,17,12,7]; k=1;
for i=1:5 % 5 traces
    for n=1:2 % neuron
    x=1/Fs(jj):1/Fs(jj):length(data.V{jj}{j}{ii(k)}(n,Fs(jj)*5:end))/Fs(jj);
    subplot(4,5,1+5*(n-1)+(k-1))
    plot(x,data.V{jj}{j}{ii(k)}(n,Fs(jj)*5:end),'color','k','linewidth',1)
    ylim([-70 -15]); xlim([10 30]);
    xlabel('Time, sec'); ylabel('V_M, mV')
    set(gca, 'Fontsize',14,'FontName','Arial'); box off
    h=display.plothorzline(data.Vth{jj}(j)); hold on
    end
    k=k+1;
end

%% Figure 4 C-F
% ERQs, cycle frequency, spike frequency

jj=1;
gH = data.gH{jj}{1}; gSyn = data.gSyn{jj}{1};

g=0.04; l=0.07;
clf,

for j=1:5 % map 
    for i=1:4 % number of features to display
        
        h1=display.bigsubplot(4,5,i,j,g,l);
        
        switch i
            
            case 1
                display.genimagesc(ERQ1{jj,j}{2}',gH,gSyn);
                myColorMap = display.linspecer(128); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j == numel(data.V{jj}); h=colorbar; ylabel(h,'ERQ'); end               
                minERQ1 = min([ERQ1{jj,1}{1}(:);ERQ1{jj,1}{2}(:)]);
                maxERQ1 = max([ERQ1{jj,5}{1}(:);ERQ1{jj,5}{2}(:)]); 
                caxis([minERQ1-0.01 maxERQ1+0.01])
                
            case 2
                display.genimagesc(ERQ1{jj,j}{1}',gH,gSyn);
                myColorMap = display.linspecer(128); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j == numel(data.V{jj}); h=colorbar; ylabel(h,'ERQ'); end
                caxis([minERQ1-0.01 maxERQ1+0.01])
                
            case 3
                display.genimagesc(FRmean{jj,j}',gH,gSyn);
                myColorMap = colormaps.viridis; myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j==5; h=colorbar; ylabel(h,'Cycle frequency, Hz'); end
                caxis([minFR-0.02 maxFR+0.02])
                
            case 4
                display.genimagesc(SpkFRmean{jj,j}',gH,gSyn);
                myColorMap = spring(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j==5; h=colorbar; ylabel(h,'Spike Frequency, Hz'); end
                caxis([minSpkFR-0.2 maxSpkFR+0.2])                          
        end
        
        xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
        xlabel('gSyn, nS'); ylabel('gH, nS'); axis square
        set(gca,'Fontsize',12,'FontName','Arial');
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
        
        for i=1:10 % create grid
            line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
            line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
        end
    end
end

%% reproduce Figure 4 - figure supplement panel A
% plot asymmerty in burst durations and ireggularity in cycle frequency

clf
jj=1;

for j = 1:5 % map
    for i=1:2 % number of features to display
        
        h1=display.bigsubplot(4,5,i,j,g,l);
        
        switch i
                
            case 1 % CV of cycle frequency
                display.genimagesc(CVFRmean{jj,j}',gH,gSyn);
                myColorMap = display.linspecer(128); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j ==5; h=colorbar; ylabel(h,'CV of cycle frequency'); end
                caxis([minCVFR-0.01 maxCVFR+0.01])
                title(['Vth=',num2str(data.Vth{jj}(j))])
                
            case 2 % assymetry of burst durations
                display.genimagesc(squeeze(asymm{jj,j})',gH,gSyn);
                myColorMap = display.linspecer(128); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                colormap(h1,myColorMap);
                if j ==5; h=colorbar; ylabel(h,'\Delta in burst duration, s'); end             
        end
        
        xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
        xlabel('gSyn'); ylabel('gH'); axis square
        set(gca,'Fontsize',10,'FontName','Arial');
        
        for i=1:10 % create grid
            line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
            line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
        end
    end
end

%% Figure 4 - figure supplement (C-F)
% slow-wave amplitude, # of spikes/burst, duty cycle, state

clf
jj=1;

for j=1:5 % map 
    for i=1:4 % number of features to display
        
        h1=display.bigsubplot(4,5,i,j,g,l);
        
        switch i
                              
            case 1
                display.genimagesc(Amean{jj,j}',gH, gSyn);
                myColorMap = colormaps.magma; myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j==5; h=colorbar; ylabel(h,'Amplitude, mV'); end
                caxis([minA-0.2 maxA+0.2])
                
            case 2
                display.genimagesc(Nspkmean{jj,j}',gH,gSyn);
                myColorMap = summer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j==5; h=colorbar; ylabel(h,'# spikes/burst'); end
                %caxis([min(Nspksmean(:))-0.2 max(Nspksescmean(:))+0.2])
                caxis([minNspk-0.3 maxNspk+0.3])
                
            case 3
                display.genimagesc(DCmean{jj,j}',gH,gSyn);
                myColorMap = copper(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j==5; h=colorbar; ylabel(h,'Duty cycle, %'); end
                %caxis([minDC-0.2 maxDC+0.2])
                caxis([0 maxDC])
                
            case 4
                display.genimagesc(state{jj}{j}',gH,gSyn);
                myColorMap = [.65 .65 .65; .25 .6 .3; 0 .6 .7; .8 .5 .3; .8 .2 0;];
                colormap(h1,myColorMap);
                if j==5
                    h=colorbar; h.Ticks=[1:5];
                    h.TickLabels={'silent','asymmetric','irregular\newlinespiking','antiphase\newlinespiking','HCO'}; caxis([1 5]);
                end
                
        end
        
        xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
        xlabel('gSyn, nS'); ylabel('gH, nS'); axis square
        set(gca,'Fontsize',12,'FontName','Arial');
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
        
        for i=1:10 % create grid
            line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
            line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
        end
    end
end

%% mechanisms of oscillatinos defined based on ERQ of individual neurons
% for the boundaries in figure 4 C,D

clf
jj=1;

for j = 1:5 % map
    for i=1:2 % number of features to display
        
        h1=display.bigsubplot(4,5,i,j,g,l);
        
        switch i

            case 1 % classified mechanism for neuron #`
                display.genimagesc(regime1{jj}{j}',gH,gSyn);
                myColorMap = [1,1,1; 0.36,0.23,0.6; 0 0.54 0.54; 0.8,0.3,0;];
                colormap(h1,myColorMap);
                if j ==5
                    h=colorbar; h.Ticks=[0:3]; h.TickLabels={'none','escape','mixed','release'};
                end
                caxis([0 3]);
                ylabel('gH');
                
            case 2 % classified mechanism for neuron #2
                display.genimagesc(regime2{jj}{j}',gH,gSyn);
                myColorMap = [1,1,1; 0.36,0.23,0.6; 0 0.54 0.54; 0.8,0.3,0;];
                colormap(h1,myColorMap);
                if j ==5
                    h=colorbar; h.Ticks=[0:3]; h.TickLabels={'none','escape','mixed','release'};
                end
                caxis([0 3]);                            
        end
        
        xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
        xlabel('gH'); xlabel('gSyn'); axis square
        set(gca,'Fontsize',10,'FontName','Arial');
        
        for i=1:10 % create grid
            line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
            line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
        end
    end
end
