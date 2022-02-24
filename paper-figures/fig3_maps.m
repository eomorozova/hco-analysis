%% load the data for "escape" and "release" (gH,gSyn) maps

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig3_fig3suppl.mat')

%% calculate half-center output characteristcs as a function of gH and gSyn

for jj = 1:2 % escape/release
    
    if jj==1; Data = data.e; else; Data = data.r; end
    
    for j = 1:numel(Data.V) % map number
        j
        for i = 1:numel(Data.V{j}) % (gh,gsyn) combination
            for ii = 1:2 % neuron number
                if ~isempty(Data.V{j}{i})==1 & ~isnan(Data.V{j}{i})==1
                    
                    Fs(j) = data.e.Fs{j};
                    V = Data.V{j}{i}(ii,5*Fs(j):end);
                    Vth = Data.Vth{j};
                    
                    ERQ{jj}{j}(ii,i) = analysis.erq(V,Vth);
                    [hcostat{jj}{j}{i}{ii}] = analysis.hco_stat(V, Fs(j), Vth);
                    
                    FR{jj}{j}(ii,i) = 1./hcostat{jj}{j}{i}{ii}.T1_mean; % cycle frequency
                    A{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.A; % slow-wave amplitude
                    SpkFR{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.SpkFreq_mean; % spike frequency
                    Nspk{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.nSpks_mean; % # spikes/burst
                    DC{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.dc_mean; % #duty cycle (based on spikes)
                    dc{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.DC_mean; % #duty cycle (based on slow-wave)
                    DC1{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.DC1_mean; % #duty cycle (based on time above synaptic threshold)
                    
                else
                    hcostat{jj}{j}{i}{ii} = []; ERQ{jj}{j}(ii,i) = NaN;
                    FR{jj}{j}(ii,i) = NaN; A{jj}{j}(ii,i) = NaN; SpkFR{jj}{j}(ii,i) = NaN;
                    Nspk{jj}{j}(ii,i) = NaN; DC{jj}{j}(ii,i) = NaN; dc{jj}{j}(ii,i) = NaN;
                    DC1{jj}{j}(ii,i) = NaN;
                end
            end
        end
        
        % mean characteristics across two neurons
        FR_{jj}(j,:) = mean(FR{jj}{j});
        dc_{jj}(j,:) = mean(dc{jj}{j});
        DC_{jj}(j,:) = mean(DC{jj}{j});
        DC1_{jj}(j,:) = mean(DC1{jj}{j});        
        A_{jj}(j,:) = mean(A{jj}{j});
        SpkFR_{jj}(j,:) = mean(SpkFR{jj}{j});
        Nspks_{jj}(j,:) = mean(Nspk{jj}{j});
        ERQ_{jj}(j,:) = mean(ERQ{jj}{j});
        
        %reshape
        ERQmean{jj,j} = reshape(ERQ_{jj}(j,:),[7,7]);
        ERQ1{jj,j}{1} = reshape(ERQ{jj}{j}(1,:),[7,7]);
        ERQ1{jj,j}{2} = reshape(ERQ{jj}{j}(2,:),[7,7]);
        
        FRmean{jj,j} = reshape(FR_{jj}(j,:),[7,7]);
        Amean{jj,j} = reshape(A_{jj}(j,:),[7,7]);
        SpkFRmean{jj,j} = reshape(SpkFR_{jj}(j,:),[7,7]);
        Nspkmean{jj,j} = reshape(Nspks_{jj}(j,:),[7,7]);
        DCmean{jj,j} = reshape(DC_{jj}(j,:),[7,7]);
        DC1mean{jj,j} = reshape(DC1_{jj}(j,:),[7,7]);        
        
    end
end

%% classify the network state (silent/asymmetric/irregular/antiphase spiking/half-center)

state=cell(1,2);

for jj = 1:2 % escape/release
    
    if jj==1; Data = data.e; else; Data = data.r; end
    
    for j = 1:numel(Data.V) % map
        for i = 1:numel(hcostat{jj}{j}) % (gh,gsyn) combination
            networkstate{jj}{j}(i)=analysis.classify_networkstate(hcostat{jj}{j}{i},Fs(j));
        end
        state{jj}{j} = reshape(networkstate{jj}{j},[7,7]);
    end
end

%% if identified state is not a half-center, set output characteristic to NaN

for jj=1:2
    for j=1:numel(state{jj})
        ERQmean{jj,j}(state{jj}{j}~=5)=NaN;
        ERQ1{jj,j}{1}(state{jj}{j}~=5)=NaN;
        ERQ1{jj,j}{2}(state{jj}{j}~=5)=NaN;
        
        FRmean{jj,j}(state{jj}{j}~=5)=NaN;
        Amean{jj,j}(state{jj}{j}~=5)=NaN;
        SpkFRmean{jj,j}(state{jj}{j}~=5)=NaN;
        Nspkmean{jj,j}(state{jj}{j}~=5)=NaN;
        DCmean{jj,j}(state{jj}{j}~=5)=NaN;
        DC1mean{jj,j}(state{jj}{j}~=5)=NaN;        
    end
end

%% calculate mean maps across all the experiments

for jj=1:2
    ERQall{jj} = nanmean(cat(3,ERQmean{jj,:}),3);
    FRall{jj} = nanmean(cat(3,FRmean{jj,:}),3);
    Aall{jj} = nanmean(cat(3,Amean{jj,:}),3);
    SpkFRall{jj} = nanmean(cat(3,SpkFRmean{jj,:}),3);
    Nspkall{jj} = nanmean(cat(3,Nspkmean{jj,:}),3);
    DCall{jj} = nanmean(cat(3,DCmean{jj,:}),3);
    DC1all{jj} = nanmean(cat(3,DC1mean{jj,:}),3);    
end

%% min and max limits for the maps

for jj=1:2 % escape/release
    minERQ(jj) = min(ERQall{jj}(:)); maxERQ(jj) = max(ERQall{jj}(:));
    minFR(jj) = min(FRall{jj}(:)); maxFR(jj) = max(FRall{jj}(:));
    minA(jj) = min(Aall{jj}(:)); maxA(jj) = max(Aall{jj}(:));
    minNspk(jj) = min(Nspkall{jj}(:)); maxNspk(jj) = max(Nspkall{jj}(:));
    minSpkFR(jj) = min(SpkFRall{jj}(:)); maxSpkFR(jj) = max(SpkFRall{jj}(:));
    minDC(jj) = min(DCall{jj}(:)); maxDC(jj) = max(DCall{jj}(:));
    minDC1(jj) = min(DC1all{jj}(:)); maxDC1(jj) = max(DC1all{jj}(:));    
end

%% % of HCOs in escape maps and average maps

for jj = 1:2 % escape/release
    k=1;
    percent_state{jj} = zeros(7,7);
    
    for j = 1:numel(Data.file) % map number
        state1{jj}=state{jj}{j};
        state1{jj}(state1{jj}~=5)=0; % not a half-center, set to 0
        percent_state{jj} = percent_state{jj}+state1{jj};
        k=k+1;
    end
    percent_state{jj}=percent_state{jj}/(numel(Data.file{jj})*5)*100; % convert to percents
end

%% plot the maps of network output in escape and release

g=0.04; l=0.07;
clf,

gH = data.e.gH{1}; gSyn = data.e.gSyn{1};

for jj=1:2 % escape/release
    for i=1:6 % number of features to display
        
        h1=display.bigsubplot(2,6,jj,i,g,l);
        
        switch i
            case 1
                display.genimagesc(percent_state{jj},gH,gSyn);
                colormap(h1,flipud(gray(10))); h=colorbar('YTick',0:20:100); ylabel(h,'% HCOs')
                %h=colorbar('Position', [0.955  0.1  0.01  0.2]);
                title('% of HCOs for each (gH, gSyn)');  ylabel('gH');
                caxis([0 100])
                
            case 2
                display.genimagesc(FRall{jj},gH,gSyn);
                myColorMap = colormaps.viridis; myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Cycle frequency, Hz'),
                title('Cycle frequency');
                caxis([min(minFR)-0.02 max(maxFR)+0.02])
                
            case 3
                display.genimagesc(Aall{jj},gH, gSyn);
                myColorMap = colormaps.magma; myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Amplitude, mV'),
                title('Slow-wave amplitude');
                caxis([min(minA)-0.2 max(maxA)+0.2])
                
            case 4
                display.genimagesc(Nspkall{jj},gH,gSyn);
                myColorMap = summer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'# spikes/burst'),
                title('# spikes/burst');
                caxis([min(minNspk)-0.3 max(maxNspk)+0.3])
                
            case 5
                display.genimagesc(SpkFRall{jj},gH,gSyn);
                myColorMap = spring(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Spike Frequency, Hz'),
                title('Spike Frequency');
                caxis([min(minSpkFR)-0.2 max(maxSpkFR)+0.2])
                %caxis([min(minSpkFR)-0.2 25])
                
            case 6
                display.genimagesc(DCall{jj},gH,gSyn);
                myColorMap = copper(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Duty cycle, %'),
                title('Duty cycle');
                %caxis([minDC-0.2 maxDC+0.2])
                caxis([0 max(maxDC)])               
        end
        
        xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
        xlabel('gSyn'); axis square
        set(gca,'Fontsize',10,'FontName','Arial');
        
        for i=1:10 % create grid
            line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
            line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
        end
    end
end

%% plot maps for the duty cycle based on the slow wave (for the supplementary figure)

g=0.05; l=0.1;
clf,

for jj=1:2 % escape/release
    
    h1=display.bigsubplot(1,2,1,jj,g,l);
    display.genimagesc(DC1all{jj},gH,gSyn);
    myColorMap = copper(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
    h=colorbar; ylabel(h,'Duty cycle, %'),
    if jj==1; title('Duty cycle (escape)'); else; title('Duty cycle (release)'); end
    %caxis([minDC-0.2 maxDC+0.2])
    %caxis([0 max(dcmean{jj}(:))])
    %caxis([0 62])
    
    xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
    xlabel('gSyn'); axis square
    set(gca,'Fontsize',14,'FontName','Arial');
    
    for i=1:10 % create grid
        line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
        line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
    end
end

% calculate mean and std for duty cycles in release and escape
dc_mean_e = nanmean(dcmean{1}(:)*100);
dc_mean_r = nanmean(dcmean{2}(:)*100);

dc_std_e = nanstd(dcmean{1}(:)*100);
dc_std_r = nanstd(dcmean{2}(:)*100);

%% plot ERQ maps
g=0.07; l=0.07;

for jj=1:2 % escape/release
    
    h1=display.bigsubplot(1,2,1,jj,g,l);
    display.genimagesc(ERQall{jj},gH,gSyn);
    myColorMap = colormaps.linspecer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
    h=colorbar; ylabel(h,'ERQ'),
    if jj==1; title('Escape'); else; title('Release'); end
    caxis([min(minERQ)-0.02 max(maxERQ)+0.02])
    
    xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
    xlabel('gSyn'); ylabel('gH'); axis square
    set(gca,'Fontsize',14,'FontName','Arial');
    
    for i=1:10 % create grid
        line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
        line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
    end
end

%% figure 3 - figure supplement
%duty cycle based on the slow wave above the synaptic threhold for escape circuits

clf

gH = data.e.gH{1}; gSyn = data.e.gSyn{1};

jj=1 % escape

display.genimagesc(DC1all{jj},gH,gSyn);
myColorMap = copper(256); myColorMap(1,:) = 1; colormap(myColorMap);
h=colorbar; ylabel(h,'Duty cycle, %'),
title('Duty cycle');
caxis([minDC1(1)-0.2 maxDC1(1)+0.2])

xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
xlabel('gSyn, nS'); ylabel('gH, nS'); axis square
set(gca,'Fontsize',10,'FontName','Arial');

for i=1:10 % create grid
    line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
    line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
end
