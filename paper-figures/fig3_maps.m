%% load the data for "escape" and "release" (gH,gSyn) maps

load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig3_maps.mat')

%% calculate half-center output characteristcs as a function of gH and gSyn

%hcostat=[];

data = [data_escape, data_release];

for jj = 1:numel(data) % escape/release    
    for j = 1:numel(data(jj).file) % map number
        j
        for i = 1:length(data(jj).V{j}) % (gh,gsyn) combination           
            for ii = 1:2 % neuron number              
                if ~isempty(data(jj).V{j}{i})==1 & ~isnan(data(jj).V{j}{i})==1
                    
                    Fs(j) = data(jj).Fs{j};
                    V = data(jj).V{j}{i}(ii,5*Fs(j):end);
                    Vth = data(jj).Vth{j};
                    
                    ERQ{jj}{j}{ii}(i) = analysis.erq(V,Vth);
                    [hcostat{jj}{j}{i}{ii}] = analysis.hco_stat(V, Fs(j));
                else
                    hcostat{jj}{j}{i}{ii} = [];  ERQ{jj}{j}{i}(ii) = NaN;
                end
            end
        end       
    end
end
%% extract burst features

FR =[]; A=[]; SpkFR=[]; Nspk=[]; DC=[]; dc=[];

for jj = 1:numel(data) % escape/release
    for j = 1:numel(data(jj).file) % map number
        for i = 1:length(data(jj).V{j}) % (gh,gsyn) combination
            for ii = 1:2 % neuron number
                if ~isempty(hcostat{jj}{j}{i}{ii})
                    FR{jj}{j}(ii,i) = 1./hcostat{jj}{j}{i}{ii}.T1_mean; % cycle frequency
                    A{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.A; % slow-wave amplitude
                    SpkFR{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.SpkFreq_mean; % spike frequency
                    Nspk{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.nSpks_mean; % # spikes/burst
                    DC{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.dc_mean; % #duty cycle (based on spikes)
                    dc{jj}{j}(ii,i) = hcostat{jj}{j}{i}{ii}.DC_mean; % #duty cycle (based on slow-wave)
                else
                    FR{jj}{j}(ii,i)=NaN; A{jj}{j}(ii,i)=NaN;
                    SpkFR{jj}{j}(ii,i)=NaN; Nspk{jj}{j}(ii,i)=NaN;
                    DC{jj}{j}(ii,i)=NaN; dc{jj}{j}(ii,i)=NaN;
                end
            end
        end
    end
end
%% calculate mean characteristics across two neurons in a circuit
FRall=[]; Aall=[]; SpkFRall=[]; Nspkall=[]; DCall=[]; dcall=[];

for jj = 1:numel(data)
    for j = 1:numel(data(1).file)
        FRall{jj}(j,:) = mean(FR{jj}{j}); 
        Aall{jj}(j,:) = mean(A{jj}{j}); 
        SpkFRall{jj}(j,:) = mean(SpkFR{jj}{j}); 
        Nspkall{jj}(j,:) = mean(Nspk{jj}{j}); 
        DCall{jj}(j,:) = mean(DC{jj}{j}); 
        dcall{jj}(j,:) = mean(dc{jj}{j});         
    end
end

%% mean characteristics across the experiments

FRmean=[]; Amean=[]; SpkFRmean=[]; Nspkmean=[]; DCmean=[]; dcmean=[];

for jj = 1:numel(data)
    FRmean{jj} = nanmean(FRall{jj}); FRstd(jj,:) = nanstd(FRall{jj});
    Amean{jj} = nanmean(Aall{jj}); Astd(jj,:) = nanstd(Aall{jj});
    SpkFRmean{jj} = mean(SpkFRall{jj}); SpkFRstd(jj,:) = nanstd(SpkFRall{jj});
    Nspkmean{jj} = nanmean(Nspkall{jj}); Nspksstd(jj,:) = nanstd(Nspkall{jj});
    DCmean{jj} = nanmean(DCall{jj}); DCstd(jj,:) = nanstd(DCall{jj});
    dcmean{jj} = nanmean(dcall{jj}); dcstd(jj,:) = nanstd(dcall{jj});  
    
    FRmean{jj} = reshape(FRmean{jj},[7,7]);
    Amean{jj} = reshape(Amean{jj},[7,7]);
    SpkFRmean{jj} = reshape(SpkFRmean{jj},[7,7]);
    Nspkmean{jj} = reshape(Nspkmean{jj},[7,7]);
    DCmean{jj} = reshape(DCmean{jj},[7,7]);
    dcmean{jj} = reshape(dcmean{jj},[7,7]);
end

%% max and min limits for the maps
minFR=min([FRmean{1}(:); FRmean{2}(:)]);
minA=min([Amean{1}(:); Amean{2}(:)]); 
minNspk=min([Nspkmean{1}(:); Nspkmean{2}(:)]); 
minSpkFR=min([SpkFRmean{1}(:); Nspkmean{2}(:)]); 
minDC=min([DCmean{1}(:); Nspkmean{2}(:)]); 
mindc=min([dcmean{1}(:); Nspkmean{2}(:)]); 

maxFR=max([FRmean{1}(~isinf(FRmean{1}(:)));FRmean{2}(~isinf(FRmean{2}(:)))]);
maxA=max([Amean{1}(~isinf(Amean{1}(:)));Amean{2}(~isinf(Amean{2}(:)))]);
maxNspk=max([Nspkmean{1}(~isinf(Nspkmean{1}(:)));Nspkmean{2}(~isinf(Nspkmean{2}(:)))]);
maxSpkFR=max([SpkFRmean{1}(~isinf(SpkFRmean{1}(:)));SpkFRmean{2}(~isinf(SpkFRmean{2}(:)))]);
maxDC=max([DCmean{1}(~isinf(DCmean{1}(:)));DCmean{2}(~isinf(DCmean{2}(:)))]);
maxdc=max([dcmean{1}(~isinf(dcmean{1}(:)));dcmean{2}(~isinf(dcmean{2}(:)))]);

%% plot an example map

jj=2; % escape
j = 10; % 1st map

%dc1 = reshape(dc{jj}{j}(1,:),[sqrt(numel(dc{jj}{j}(1,:))),sqrt(numel(dc{jj}{j}(1,:)))]);

FR1=[];
%FR1 = reshape(FR{jj}{j}(2,:),[7,7]);
FR1 = reshape(FRall{jj}(j,:),[7,7]);
    
gH = data(jj).gH{1}; gSyn = data(jj).gSyn{1};

clf
subplot(1,3,1)
%plot(data(jj).V{j}{2}(2,:))
title('gH=150, gsyn=150')
set(gca,'Fontsize',16,'FontName','Arial');

subplot(1,3,2)
display.genimagesc(FR1,gH, gSyn);
colormap([1 1 1; colormaps.viridis]);
caxis([min(FR1(:))-0.02 max(FR1(:))+0.02])
%caxis([min(dc1(:))-0.02 max(dc1(:)*100)+0.02])  
h=colorbar; ylabel(h,'Cycle frequency, Hz')
xlabel('gSyn, nS'); ylabel('gH, nS')
%title('Cycle frequency');
axis square
xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
set(gca,'Fontsize',16,'FontName','Arial');

% compare with previously generated map
subplot(1,3,3)
%display.genimagesc(maps_escape{j}.FR{1},gH, gSyn);
display.genimagesc(maps_release{j}.FR{1},gH, gSyn);
colormap([1 1 1; colormaps.viridis]);
caxis([min(FR1(:))-0.02 max(FR1(:))+0.02])
%caxis([min(dc1(:))-0.02 max(dc1(:)*100)+0.02])  
h=colorbar; ylabel(h,'Cycle frequency, Hz')
xlabel('gSyn, nS'); ylabel('gH, nS')
%title('Cycle frequency');
axis square
xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
set(gca,'Fontsize',16,'FontName','Arial');

 %% classify the state (silent, spiking, half-center,..)
 
 networkstate=[]; state=[]; percentSingleSpikeBursts=[];
  
 for jj=1:numel(data) % escape/release
     
     for j = 1:numel(data(jj).file) % map
         
         for i = 1:length(data(jj).V{j}) % (gh,gsyn) combination
             
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
         
         networkstate{1}{10}(43:49)=1;% fix manually (incomplete map)
         networkstate{2}{5}(49)=1;% fix manually (incomplete map)
         
         state{jj}{j} = reshape(networkstate{jj}{j},[7,7]);
     end
 end
 
 %% % of HCOs in escape maps and average maps
%gH_1=gH{1}{1}; gSyn_1=gSyn{1}{1};

for jj = 1:numel(data)
    k=1;
    percent_state{jj} = zeros(7,7);
    for j = 1:numel(data(jj).file) % map number
        state1{jj}=state{jj}{j};
        state1{jj}(state1{jj}~=5)=0; % not a half-center, set to 0
        percent_state{jj} = percent_state{jj}+state1{jj};
        k=k+1;
    end
    percent_state{jj}=percent_state{jj}/(numel(data(jj).file)*5)*100; % convert to percents
end

%%
jj =1; % escape

g=0.04; l=0.07;
clf,
for i=1:6 % number of features to display
    
    h1=display.bigsubplot(2,6,1,i,g,l);
    
    switch i
        case 1
            display.genimagesc(percent_state{jj},gH,gSyn);
            colormap(h1,flipud(gray(10))); h=colorbar('YTick',0:20:100); ylabel(h,'% HCOs')
            %h=colorbar('Position', [0.955  0.1  0.01  0.2]);
            title('% of HCOs for each (gH, gSyn)');  ylabel('gH');
            caxis([0 100])
            
        case 2
            display.genimagesc(FRmean{jj},gH,gSyn);
            myColorMap = colormaps.viridis; myColorMap(1,:) = 1; colormap(h1,myColorMap);
            h=colorbar; ylabel(h,'Cycle frequency, Hz'),
            title('Cycle frequency');
            %caxis([min(FRescmean(:))-0.02 max(FRescmean(:))+0.02])
            caxis([minFR-0.02 maxFR+0.02])
            
        case 3
            display.genimagesc(Amean{jj},gH, gSyn);
            myColorMap = colormaps.magma; myColorMap(1,:) = 1; colormap(h1,myColorMap);
            h=colorbar; ylabel(h,'Amplitude, mV'),
            title('Slow-wave amplitude');
            caxis([minA-0.2 maxA+0.2])
            
        case 4
            display.genimagesc(Nspkmean{jj},gH,gSyn);
            myColorMap = summer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
            h=colorbar; ylabel(h,'# spikes/burst'),
            title('# spikes/burst');
            caxis([min(Nspksescmean(:))-0.2 max(Nspksescmean(:))+0.2])
            caxis([minNspks-0.3 maxNspks+0.3])
            
        case 5
            display.genimagesc(SpkFRmean{jj},gH,gSyn);
            myColorMap = spring(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
            h=colorbar; ylabel(h,'Spike Frequency, Hz'),
            title('Spike Frequency');
            caxis([minSpkFR-0.2 maxSpkFR+0.2])
            
        case 6
            display.genimagesc(DCmean{jj},gH,gSyn);
            myColorMap = copper(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
            h=colorbar; ylabel(h,'Duty cycle, %'),
            title('Duty cycle');
            caxis([minDC-0.2 maxDC+0.2])
            
    end
    
    xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
    xlabel('gSyn'); axis square
    set(gca,'Fontsize',10,'FontName','Arial');
    
    for i=1:10 % create grid
        line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
        line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
    end
end

%% release
for i=1:6
    h1=display.bigsubplot(2,6,2,i,g,l);
    switch i
        case 1
   % display.genimagesc(state_release,gH,gSyn);
    colormap(h1,flipud(gray(10))); h=colorbar('YTick',0:20:100); ylabel(h,'% HCOs')
      %h=colorbar('Position', [0.955  0.1  0.01  0.2]);
    title('% of HCOs for each (gH, gSyn)');  ylabel('gH');
    
    case 2
        display.genimagesc(FRrelmean,gH,gSyn);
        myColorMap = colormaps.viridis; myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Cycle frequency, Hz'),
        title('Cycle frequency');
        %caxis([min(FRrelmean(:))-0.02 max(FRrelmean(:))+0.02])
        caxis([minFR-0.02 maxFR+0.02])
        
    case 3
        display.genimagesc(Arelmean,gH,gSyn);
        myColorMap = colormaps.magma; myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Amplitude, mV'),
        title('Slow-wave amplitude');
        %caxis([min(Arelmean(:))-0.2 max(Arelmean(:))+0.2])
        caxis([minA-0.2 maxA+0.2])          
        
    case 4
        display.genimagesc(Nspksrelmean,gH,gSyn);
        myColorMap = summer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'# spikes/burst'),
        title('# spikes/burst');
        %caxis([min(Nspksrelmean(:))-0.2 max(Nspksrelmean(:))+0.2])  
        caxis([minNspks-0.3 maxNspks+0.3])         
        
    case 5
        display.genimagesc(SpkFRrelmean,gH,gSyn);
        myColorMap = spring(256); myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Spike Frequency, Hz'),
        title('Spike Frequency');
        %caxis([min(SpkFRrelmean(:))-0.1 max(SpkFRrelmean(:))+0.1])
        caxis([minSpkFR-0.2 maxSpkFR+0.2])           
        
    case 6
        display.genimagesc(DCrelmean,gH,gSyn);
        myColorMap = copper(256); myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Duty cycle, %'),
        title('Duty cycle');
        %caxis([min(DCrelmean(:))-0.1 max(DCrelmean(:))+0.1]) 
        caxis([minDC-0.2 maxDC+0.2])           
        
    end
    
    xticks([150,300,450,600,750,900,1050]); yticks([150,300,450,600,750,900,1050]);
    xlabel('gSyn'); axis square
    set(gca,'Fontsize',10,'FontName','Myriard pro');
    
    for i=1:10 % create grid
        line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
        line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
    end
end
