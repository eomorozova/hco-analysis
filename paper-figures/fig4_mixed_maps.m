
%% load all the data with the mixed maps

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig4.mat')

%% calculate half-center output characteristcs as a function of gH and gSyn

hcostat=[]; ERQ=[];

for jj = 1 %:numel(data) % experiment
    jj
    for j = 1:numel(data{jj}.V) % map number
        for i = 1:length(data{jj}.V{j}) % (gh,gsyn) combination
            for ii = 1:2 % neuron number
                if ~isempty(data{jj}.V{j}{i})==1 & ~isnan(data{jj}.V{j}{i})==1
                    
                    Fs(j) = data{jj}.Fs;
                    V = data{jj}.V{j}{i}(ii,5*Fs(j):end);
                    Vth = data{jj}.Vth(j);
                    
                    ERQ{jj}{j}(ii,i) = analysis.erq(V,Vth);
                    [hcostat{jj}{j}{i}{ii}] = analysis.hco_stat(V, Fs(j));
                    
                else
                    hcostat{jj}{j}{i}{ii} = [];  ERQ{jj}{j}(ii,i) = NaN;
                end
            end
        end
    end
end

%% extract burst features

FR =[]; A=[]; SpkFR=[]; Nspk=[]; DC=[]; dc=[];

meanFR=[]; stfFR=[]; CVFR=[];

for jj = 1 %1:numel(data) % experiment
    for j = 1:numel(data{jj}.V) % map number
        for i = 1:length(data{jj}.V{j}) % (gh,gsyn) combination
            for ii = 1:2 % neuron number
                if ~isempty(hcostat{jj}{j}{i}{ii})
                    FR{jj}{j}(ii,i) = 1./hcostat{jj}{j}{i}{ii}.T1_mean; % cycle frequency
                    A{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.A; % slow-wave amplitude
                    SpkFR{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.SpkFreq_mean; % spike frequency
                    Nspk{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.nSpks_mean; % # spikes/burst
                    DC{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.dc_mean; % #duty cycle (based on spikes)
                    dc{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.DC_mean; % #duty cycle (based on slow-wave)

                    % evaluate asymmetry
                    BurDur{jj}{j}{ii,i} = hcostat{jj}{j}{i}{ii}.BurDur; % based on slow-wave
                    active{jj}{j}(ii,i) = sum(hcostat{jj}{j}{i}{ii}.BurDur(2:end-1));
                    
                    BrstDur{jj}{j}{ii,i} =  hcostat{jj}{j}{i}{ii}.burststat.BuDur;
                    meanBrstDur{jj}{j}(ii,i) =  mean(BrstDur{jj}{j}{ii,i});
                    
                    % coefficient of variation of burst characteristics
                    meanFR{jj}{j}(ii,i) = mean(hcostat{jj}{j}{i}{ii}.burststat.BuFreq(2:end-2));
                    stdFR{jj}{j}(ii,i) = std(hcostat{jj}{j}{i}{ii}.burststat.BuFreq(2:end-2));
                    CVFR{jj}{j}(ii,i) = stdFR{jj}{j}(ii,i)/meanFR{jj}{j}(ii,i);
                    
                    meanNspk{jj}{j}(ii,i) = mean(hcostat{jj}{j}{i}{ii}.burststat.nSp(2:end-2));
                    stdNspk{jj}{j}(ii,i) = std(hcostat{jj}{j}{i}{ii}.burststat.nSp(2:end-2));
                    CVNspk{jj}{j}(ii,i) = stdNspk{jj}{j}(ii,i)/meanNspk{jj}{j}(ii,i);
                    
                    meanSpkFR{jj}{j}(ii,i) = mean(hcostat{jj}{j}{i}{ii}.burststat.SpFreq(2:end-2));
                    stdSpkFR{jj}{j}(ii,i) = std(hcostat{jj}{j}{i}{ii}.burststat.SpFreq(2:end-2));
                    CVSpkFR{jj}{j}(ii,i) = stdSpkFR{jj}{j}(ii,i)/meanSpkFR{jj}{j}(ii,i);
                    
                else
                    FR{jj}{j}(ii,i)=NaN; A{jj}{j}(ii,i)=NaN;
                    SpkFR{jj}{j}(ii,i)=NaN; Nspk{jj}{j}(ii,i)=NaN;
                    DC{jj}{j}(ii,i)=NaN; dc{jj}{j}(ii,i)=NaN;
                                      
                end
            end
        end
    end
end

%% quantify irregularity in bursting
hcostat=[]; 
meanNspk=[];  stdNspk=[]; CVNspk=[];
meanFR=[];  stdFR=[]; CVFR=[];

i=1;

for j=1:2
    
    meanNspk(i,j) = mean(hcostat{j}.burststat.nSp(2:end-2));
    stdNspk(i,j) = std(hcostat{j}.burststat.nSp(2:end-2));
    CVNspk(i,j) = stdNspk(j)/meanNspk(j);

    meanFR(i,j) = mean(hcostat{j}.burststat.BuFreq(2:end-2));
    stdFR(i,j) = std(hcostat{j}.burststat.BuFreq(2:end-2));
    CVFR(i,j) = stdFR(j)/meanFR(j);


end

%% classify the network state (silent, spiking, half-center,..)

networkstate=[]; state=[]; percentSingleSpikeBursts=[];

for jj = 1 %1:numel(data) % experiment
    
    for j = 1:numel(data{jj}.V) % map
        
        for i = 1:length(data{jj}.V{j}) % (gh,gsyn) combination
            
            if ~isempty(hcostat{jj}{j}{i}{1})==1
                st1 = hcostat{jj}{j}{i}{1}.st/Fs(j); st2 = hcostat{jj}{j}{i}{2}.st/Fs(j); % spike times in s
                
                percentSingleSpikeBursts{jj}{j}(i) = calcPercentSingleSpikeBursts(st1,st2);
                
                % both neurons spike less than 5 times in a minute
                if length(st1(st1<60))<=5  && length(st2(st2<60))<=5
                    
                    networkstate{jj}{j}(i) = 1; % silent
                    
                    % one neuron spikes, the other one is silent (<5 spikes/minute)
                elseif (length(st1(st1<60))<5  && length(st2(st2<60))>5 || length(st1(st1<60))>5  && length(st2(st2<60))<5)
                    
                    networkstate{jj}{j}(i) = 2; % asymmetric
                    
                    % both neurons spike more than 5 spikes per minute
                elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 ...
                        && hcostat{jj}{j}{i}{1}.T1_mean==0 && hcostat{jj}{j}{i}{2}.T1_mean==0)
                    
                    % if more than 80% spikes are alternating
                    if percentSingleSpikeBursts{jj}{j}(i)>0.8
                        
                        networkstate{jj}{j}(i) = 4; % antiphase spiking
                        
                    else
                        
                        networkstate{jj}{j}(i) = 3; % irregular spiking
                        
                    end
                    
                    % alternating bursting pattern of activity
                elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 && ...
                        hcostat{jj}{j}{i}{1}.T1_mean>0 && hcostat{jj}{j}{i}{2}.T1_mean>0)
                    
                    networkstate{jj}{j}(i) = 5; % half-center oscillator
                    
                    % single spike half-center oscillator
                elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 && ...
                        hcostat{jj}{j}{i}{1}.T1_mean>0 && hcostat{jj}{j}{i}{2}.T1_mean>0 && ...
                        hcostatr{jj}{j}{i}{1}.nSpks_mean==0 && hcostat{jj}{j}{i}{2}.nSpks_mean==0)
                    
                    networkstate{jj}{j}(i) = 4; % antiphase spiking
                    
                else
                    
                    networkstate{jj}{j}(i) = 6; % misclassified
                    
                end
            end
        end
        
        state{jj}{j} = reshape(networkstate{jj}{j},[7,7]);
    end
end

 %% if identified state is not a half-center, make sure all the half-center characteristics are NaN

 for jj=1 % experiment
     for j=1:numel(state{jj}) % map number
         for i=1:length(networkstate{jj}{j}) % (gh,gsyn) combination
             if state{jj}{j}(i)~=5 % if not a half-center
                 ERQ{jj}{j}(:,i)=NaN;
                 FR{jj}{j}(:,i)=NaN;
                 A{jj}{j}(:,i)=NaN;
                 SpkFR{jj}{j}(:,i)=NaN;
                 Nspk{jj}{j}(:,i)=NaN;
                 DC{jj}{j}(:,i)=NaN;
                 dc{jj}{j}(:,i)=NaN;
                 
                 CVFR{jj}{j}(:,i)=NaN;
             end
         end
     end
 end

%% calculate mean characteristics across two neurons in a circuit

for jj = 1 %:numel(data)
    for j = 1:numel(data{jj}.V)
        ERQall{jj}(j,:) = mean(ERQ{jj}{j});
        FRall{jj}(j,:) = mean(FR{jj}{j});
        Aall{jj}(j,:) = mean(A{jj}{j});
        SpkFRall{jj}(j,:) = mean(SpkFR{jj}{j});
        Nspkall{jj}(j,:) = mean(Nspk{jj}{j});
        DCall{jj}(j,:) = mean(DC{jj}{j});
        dcall{jj}(j,:) = mean(dc{jj}{j});
        
        CVFRall{jj}(j,:) = mean(CVFR{jj}{j});
        
        %reshape
        ERQmean{jj}{j} = reshape(ERQall{jj}(j,:),[7,7])';
        for i=1:2
            ERQ1{jj}{j} = reshape(ERQ{jj}{j}(i,:),[7,7]);
        end
        FRmean{jj}{j} = reshape(FRall{jj}(j,:),[7,7])';
        Amean{jj}{j} = reshape(Aall{jj}(j,:),[7,7])';
        SpkFRmean{jj}{j} = reshape(SpkFRall{jj}(j,:),[7,7])';
        Nspkmean{jj}{j} = reshape(Nspkall{jj}(j,:),[7,7])';
        DCmean{jj}{j} = reshape(DCall{jj}(j,:),[7,7])';
        dcmean{jj}{j} = reshape(dcall{jj}(j,:),[7,7])';
        
        CVFRmean{jj}{j} = reshape(CVFRall{jj}(j,:),[7,7])';
    end
end

%% max and min limits for the maps

minERQ = min(ERQall{jj}(:)); maxERQ = max(ERQall{jj}(:)); 
minFR = min(FRall{jj}(:)); maxFR = max(FRall{jj}(:)); 
minA = min(Aall{jj}(:)); maxA = max(Aall{jj}(:)); 
minNspk = min(Nspkall{jj}(:)); maxNspk = max(Nspkall{jj}(:)); 
minSpkFR = min(SpkFRall{jj}(:)); maxSpkFR = max(SpkFRall{jj}(:)); 
minDC = min(DCall{jj}(:)); maxDC = max(DCall{jj}(:)); 
mindc = min(dcall{jj}(:)); maxdc = max(dcall{jj}(:)); 

minCVFR = min(CVFRall{jj}(:)); maxCVFR = max(CVFRall{jj}(:)); 

%% reproduce Figure 4 C-E

jj=1;
gH = data{jj}.gH{1}; gSyn = data{jj}.gSyn{1};

g=0.04; l=0.07;
clf,

for j=1:5 % map 
    for i=1:3 % number of features to display
        
        h1=display.bigsubplot(3,5,i,j,g,l);
        
        switch i
            
            case 1
                display.genimagesc(ERQmean{jj}{j},gH,gSyn);
                myColorMap = display.linspecer(128); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j ==5; h=colorbar; ylabel(h,'ERQ'); end
                caxis([minERQ-0.02 maxERQ+0.02])
                title(['Vth=',num2str(data{jj}.Vth(j))])
                
            case 2
                display.genimagesc(FRmean{jj}{j},gH,gSyn);
                myColorMap = colormaps.viridis; myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j==5; h=colorbar; ylabel(h,'Cycle frequency, Hz'); end
                caxis([minFR-0.02 maxFR+0.02])
                
            case 3
                display.genimagesc(SpkFRmean{jj}{j},gH,gSyn);
                myColorMap = spring(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j==5; h=colorbar; ylabel(h,'Spike Frequency, Hz'); end
                caxis([minSpkFR-0.2 maxSpkFR+0.2])
                %caxis([minSpkFR-0.2 25])
                
                
        end
        
        xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
        xlabel('gSyn, nS'); ylabel('gH, nS'); axis square
        set(gca,'Fontsize',14,'FontName','Arial');
        
        for i=1:10 % create grid
            line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
            line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
        end
    end
end

%% reproduce supplementaty figure for figure 4

jj=1;
gH = data{jj}.gH{1}; gSyn = data{jj}.gSyn{1};

g=0.04; l=0.07;
clf,

for j=1:5 % map 
    for i=1:4 % number of features to display
        
        h1=display.bigsubplot(4,5,i,j,g,l);
        
        switch i
            case 1
                display.genimagesc(state{jj}{j}',gH,gSyn);
                myColorMap = [0.3,0.33,0.35; 0.25,0.4,0.3; 0,0.4,0.5; 0.8,0.5,0.3; 0.6,0.2,0;];
                colormap(h1,myColorMap);
                h=colorbar;
                h.Ticks=[1:5]; h.TickLabels={'silent','asymmetric','irregular\newlinespiking','antiphase\newlinespiking','HCO'}; caxis([1 5]);
                ylabel('gH');
                axis square              
                title(['Vth=', num2str(data{jj}.Vth(j)),'mV'])
                              
            case 2
                display.genimagesc(Amean{jj}{j},gH, gSyn);
                myColorMap = colormaps.magma; myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Amplitude, mV'),
                title('Slow-wave amplitude');
                caxis([minA-0.2 maxA+0.2])
                
            case 3
                display.genimagesc(Nspkmean{jj}{j},gH,gSyn);
                myColorMap = summer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'# spikes/burst'),
                title('# spikes/burst');
                %caxis([min(Nspksmean(:))-0.2 max(Nspksescmean(:))+0.2])
                caxis([minNspk-0.3 maxNspk+0.3])
                
            case 4
                %display.genimagesc(DCmean{jj},gH,gSyn);
                display.genimagesc(DCmean{jj}{j},gH,gSyn);
                myColorMap = copper(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Duty cycle, %'),
                title('Duty cycle');
                %caxis([minDC-0.2 maxDC+0.2])
                caxis([0 50])
                
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

%% plot CV of the characteristics of the circuit output

jj=1;
gH = data{jj}.gH{1}; gSyn = data{jj}.gSyn{1};

g=0.04; l=0.07;
clf,

for j=1:5 % map
    for i=1%:3 % number of features to display
        
        h1=display.bigsubplot(3,5,i,j,g,l);
        
        switch i
            
            case 1
                display.genimagesc(CVFRmean{jj}{j},gH,gSyn);
                myColorMap = display.linspecer(128); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j ==5; h=colorbar; ylabel(h,'CV of cycle frequency'); end
                %caxis([minCVFR-0.01 maxCVFR+0.01])
                caxis([minCVFR-0.01 0.15])
                title(['Vth=',num2str(data{jj}.Vth(j))])
        end
    end
end

%% asymmerty

jj=1;
gH = data{jj}.gH{1}; gSyn = data{jj}.gSyn{1};

g=0.04; l=0.07;
clf,

for j=1:5 % map
    for i=1%:3 % number of features to display
        
        active1 = reshape(active{jj}{j}(1,:),[7,7]);
        active2 = reshape(active{jj}{j}(2,:),[7,7])
        h1=display.bigsubplot(3,5,i,j,g,l);
        
        switch i
            
            case 1
                display.genimagesc(abs(active1-active2)',gH,gSyn);
                myColorMap = display.linspecer(128); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                if j ==5; h=colorbar; ylabel(h,'CV of cycle frequency'); end
                %caxis([minCVFR-0.01 maxCVFR+0.01])
                %caxis([minCVFR-0.01 0.15])
                title(['Vth=',num2str(data{jj}.Vth(j))])
        end
    end
end

%%
i=8; %experiment

clf

for i=1:8
    subtightplot(2,4,i,g,l,l)
    
    j=1; %gmi
    active{i,j}(active{i,j}==0)=NaN;
    %plot(data.Vth{i}{j}, abs(active{i,j}(:,2) - active{i,j}(:,1)),'.-','linewidth',1)
    plot(data.Vth{i}{j}, abs(meanBrstDur{i,j}(:,2) - meanBrstDur{i,j}(:,1)),'.-','linewidth',1)
    
    j=2;
    active{i,j}(active{i,j}==0)=NaN;
    %hold on, plot(data.Vth{i}{j}, abs(active{i,j}(:,2) - active{i,j}(:,1)),'.-','linewidth',1)
    hold on, plot(data.Vth{i}{j}, abs(meanBrstDur{i,j}(:,2) - meanBrstDur{i,j}(:,1)),'.-','linewidth',1)
    %hold on, plot(data.Vth{i}{j},active{i,j}(:,2),'.-','linewidth',1)
   % ylim([0 35])
    xlim([-53 -33]); xticks([-54:2:-34])
    legend('ctrl','gmi=150nS')
end
