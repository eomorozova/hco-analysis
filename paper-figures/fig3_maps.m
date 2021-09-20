%%
load('C:\Users\moroz\Documents\Katya_DynamicClamp\hco_related_codes\HCO_maps.mat')
load('C:\Users\moroz\Documents\Katya_DynamicClamp\hco_related_codes\gH_gSyn.mat')

%%
load('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig3_maps.mat')

%% calculate half-center output characteristcs as a function of the synaptic threhold
 
hcostatr = cell(1,numel(data_release.V)); 
hcostate = cell(1,numel(data_escape.V));

ERQr = cell(1,numel(data_release.V));
ERQe = cell(1,numel(data_escape.V));

for j = 1:numel(data_release.file) % experiment
    j
    for i = 1:length(data_release.V{j}) % (gh,gsyn) combination
        for ii = 1:2 % neuron number
            if ~isempty(data_release.V{j}{i})==1
                
                Vr = data_release.V{j}{i}(ii,5*Fs(j):end);
                Vth = data_release.Vth{j};
                
                ERQr{j}{ii}(i) = analysis.erq(Vr,Vth);
                
                [hcostatr{j}{i}{ii}] = hco_stat(Vr, Fs(j));
            else
                hcostatr{j}{i}{ii} = [];  ERQr{j}{i}(ii) = NaN;
            end
        end
    end
end

%%
 FR = cell(1,numel(data_release.V)); A = cell(1,numel(data_release.V));
 SpkFR = cell(1,numel(data_release.V)); Nspks = cell(1,numel(data_release.V));
 DC = cell(1,numel(data_release.V));
 
 for j = 1:numel(data_release.V) % experiment
     for  i = 1:length(data_release.V{j}) % (gh,gsyn) combination 
         for ii = 1:2 % neuron number
             if ~isempty(hcostatr{j}{i}{ii})
                 FR{j}{ii}(i) = 1./hcostatr{j}{i}{ii}.T1_mean; % cycle frequency
                 A{j}{ii}(i) = hcostatr{j}{i}{ii}.A; % slow-wave amplitude
                 SpkFR{j}{ii}(i) = hcostatr{j}{i}{ii}.SpkFreq_mean; % spike frequency
                 Nspks{j}{ii}(i) = hcostatr{j}{i}{ii}.nSpks_mean; % # spikes/burst
                 DC{j}{ii}(i) = hcostatr{j}{i}{ii}.dc_mean; % #duty cycle
             else
                 FR{j}{i}(ii)=NaN; A{j}{i}(ii)=NaN;
                 SpkFR{j}{i}(ii)=NaN; Nspks{j}{i}(ii)=NaN; DC{j}{i}(ii)=NaN;
             end
         end
     end
 end
 
%%
FR1=[];
j =1
FR{j}{1}(FR{j}{1}==Inf)=NaN;
FR1{j} = reshape(FR{j}{1},[sqrt(numel(FR{j}{1})),sqrt(numel(FR{j}{1}))]);

%% plot example map
gH = data_release.gH{1}; gSyn = data_release.gSyn{1};
clf,
h1=display.bigsubplot(2,6,2,1,g,l);
display.genimagesc(FR1{1},gH, gSyn);
myColorMap = [1 1 1; colormaps.viridis]; colormap(h1,myColorMap);
caxis([min(FR1{1}(:))-0.01 max(FR1{1}(:))+0.01])  
h=colorbar; ylabel(h,'Cycle frequency, Hz'),
title('Cycle frequency');

 %% classify the state (silent, spiking, half-center,..)
 %networkstate=[]; state=[]; percentSingleSpikeBursts=[];
 
 for j = 1:numel(data_release.file) % experiment
     
     for i=length(data_release.V{j}) % gh,gsyn combination
         
         if ~isempty(hcostat{j}{i}{1})==1
             st1 = hcostat{j}{i}{1}.st/Fs; st2 = hcostat{j}{i}{2}.st/Fs; % spike times in s
             
             percentSingleSpikeBursts{j}{k}(i) = calcPercentSingleSpikeBursts(st1,st2);
             
             % both neurons spike less than 5 times in a minute
             if length(st1(st1<60))<=5  && length(st2(st2<60))<=5 
                 
                 networkstate{j}(i) = 1; % silent
             
             % one neuron spikes, the other one is silent (<5 spikes/minute)     
             elseif (length(st1(st1<60))<5  && length(st2(st2<60))>5 || (st1(st1<60))>5  && length(st2(st2<60))<5)
                 
                 networkstate{j}{k}(i) = 2; % asymmetric
             
             % both neurons spike more than 5 spikes per minute    
             elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 ...
                     && hcostat{j}{k}{i}{1}.T1_mean==0 && hcostat{j}{k}{i}{2}.T1_mean==0)
                 
                 % if more than 80% spikes are alternating
                 if percentSingleSpikeBursts{j}{k}(i)>0.8 
                     
                     networkstate{j}{k}(i) = 4; % antiphase spiking
                    
                 else
                     
                     networkstate{j}{k}(i) = 3; % irregular spiking
                     
                 end
                 
             % alternating bursting pattern of activity    
             elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 && ...
                     hcostat{j}{k}{i}{1}.T1_mean>0 && hcostat{j}{k}{i}{2}.T1_mean>0)
                 
                 networkstate{j}{k}(i) = 5; % half-center oscillator
              
             % single spike half-center oscillator    
             elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 && ...
                     hcostat{j}{k}{i}{1}.T1_mean>0 && hcostat{j}{k}{i}{2}.T1_mean>0 && ...
                     hcostat{j}{k}{i}{1}.nSpks_mean==0 && hcostat{j}{k}{i}{2}.nSpks_mean==0)
                
                 networkstate{j}{k}(i) = 4; % antiphase spiking
                 
             else
                 
                 networkstate{j}{k}(i) = 6; % misclassified
                 
             end
         end
     end
     
     %networkstate{17}{1}(43:49)=1;% fix manually (incomplete map)
     state{j} = reshape(networkstate{j}{k},[7,7]);
 end
 
%% % of HCOs in escape maps and average maps
%gH_1=gH{1}{1}; gSyn_1=gSyn{1}{1};
state_escape=zeros(7,7); k=1;
FRescall=[]; Aescall=[]; Nspksescall=[]; SpkFRescall=[]; DCescall=[];
for j = [1:length(maps_escape)]

    state1=maps_escape{j}.state;
    state1(state1~=5)=0; % not an HCO st to 0
    state_escape=state_escape+state1;   
    if  j==7 % second neuron
        ii=2;
    FR_1=maps_escape{j}.FR{ii}; FR_1(FR_1==-Inf)=NaN; FRescall(:,:,k)=FR_1;
    A_1=maps_escape{j}.A{ii}; A_1(A_1==0)=NaN; Aescall(:,:,k)=A_1;
    Nspks_1=maps_escape{j}.Nspks{ii}; Nspks_1(Nspks_1==0)=NaN; Nspksescall(:,:,k)=Nspks_1;    
    SpkFR_1=maps_escape{j}.SpkFR{ii}; SpkFR_1(SpkFR_1==0)=NaN; SpkFRescall(:,:,k)=SpkFR_1;       
    DC_1=maps_escape{j}.DC{ii}; DC_1(DC_1==0)=NaN; DCescall(:,:,k)=DC_1;   
    else
    FR_1=maps_escape{j}.FR{1}; FR_1(FR_1==-Inf)=NaN; FRescall(:,:,k)=FR_1;
    A_1=maps_escape{j}.A{1}; A_1(A_1==0)=NaN; Aescall(:,:,k)=A_1;
    Nspks_1=maps_escape{j}.Nspks{1}; Nspks_1(Nspks_1==0)=NaN; Nspksescall(:,:,k)=Nspks_1;    
    SpkFR_1=maps_escape{j}.SpkFR{1}; SpkFR_1(SpkFR_1==0)=NaN; SpkFRescall(:,:,k)=SpkFR_1;       
    DC_1=maps_escape{j}.DC{1}; DC_1(DC_1==0)=NaN; DCescall(:,:,k)=DC_1;  
    end
    k=k+1;
end

 %SpkFRescall(:,:,7)=NaN;

FRescmean=nanmean(FRescall,3); FRescstd=nanstd(FRescall,[],3);
Aescmean=nanmean(Aescall,3); Aescstd=nanstd(Aescall,[],3);
Nspksescmean=nanmean(Nspksescall,3); Nspksescstd=nanstd(Nspksescall,[],3);
SpkFRescmean=nanmean(SpkFRescall,3); SpkFRescstd=nanstd(SpkFRescall,[],3);
DCescmean=nanmean(DCescall,3); DCescstd=nanstd(DCescall,[],3);
state_escape=state_escape/(length(maps_escape)*5)*100;

state_release=zeros(7,7); k=1;
FRrelall=[]; Arelall=[]; Nspksrelall=[]; SpkFRrelall=[]; DCrelall=[];
for j = 1:length(maps_release)

    state1=maps_release{j}.state;
    state1(state1~=5)=0; % not an HCO st to 0
    state_release=state_release+state1;   
    if j==2 || j==6 || j==8 % neuron #2
        ii=2;
        FR_1=maps_release{j}.FR{ii}; FR_1(FR_1==-Inf)=NaN; FRrelall(:,:,k)=FR_1;
        A_1=maps_release{j}.A{ii}; A_1(A_1==0)=NaN; Arelall(:,:,k)=A_1;
        Nspks_1=maps_release{j}.Nspks{ii}; Nspks_1(Nspks_1==0)=NaN; Nspksrelall(:,:,k)=Nspks_1;
        SpkFR_1=maps_release{j}.SpkFR{ii}; SpkFR_1(SpkFR_1==0)=NaN; SpkFRrelall(:,:,k)=SpkFR_1;
        DC_1=maps_release{j}.DC{ii}; DC_1(DC_1==0)=NaN; DCrelall(:,:,k)=DC_1;
    else
        
        FR_1=maps_release{j}.FR{1}; FR_1(FR_1==-Inf)=NaN; FRrelall(:,:,k)=FR_1;
        A_1=maps_release{j}.A{1}; A_1(A_1==0)=NaN; Arelall(:,:,k)=A_1;
        Nspks_1=maps_release{j}.Nspks{1}; Nspks_1(Nspks_1==0)=NaN; Nspksrelall(:,:,k)=Nspks_1;
        SpkFR_1=maps_release{j}.SpkFR{1}; SpkFR_1(SpkFR_1==0)=NaN; SpkFRrelall(:,:,k)=SpkFR_1;
        DC_1=maps_release{j}.DC{1}; DC_1(DC_1==0)=NaN; DCrelall(:,:,k)=DC_1;
    end
    k=k+1;
end

FRrelmean=nanmean(FRrelall,3); FRrelstd=nanstd(FRrelall,[],3);
Arelmean=nanmean(Arelall,3); Arelstd=nanstd(Arelall,[],3);
Nspksrelmean=nanmean(Nspksrelall,3); Nspksrelstd=nanstd(Nspksrelall,[],3);
SpkFRrelmean=nanmean(SpkFRrelall,3); SpkFRrelstd=nanstd(SpkFRrelall,[],3);
DCrelmean=nanmean(DCrelall,3); DCrelstd=nanstd(DCrelall,[],3);
state_release=state_release/(length(maps_release)*5)*100;
%
minFR=min([min(FRescmean(:)),min(FRrelmean(:))]); maxFR=max([max(FRescmean(:)),max(FRrelmean(:))]);
minA=min([min(Aescmean(:)),min(Arelmean(:))]); maxA=max([max(Aescmean(:)),max(Arelmean(:))]);
minNspks=min([min(Nspksescmean(:)),min(Nspksrelmean(:))]); maxNspks=max([max(Nspksescmean(:)),max(Nspksrelmean(:))]);
minSpkFR=min([min(SpkFRescmean(:)),min(SpkFRrelmean(:))]); maxSpkFR=max([max(SpkFRescmean(:)),max(SpkFRrelmean(:))]);
minDC=min([min(DCescmean(:)),min(DCrelmean(:))]); maxDC=max([max(DCescmean(:)),max(DCrelmean(:))]);
%

%%

g=0.04; l=0.07;
clf, 
for i=1:6
    h1=display.bigsubplot(2,6,1,i,g,l);
    switch i
        case 1
    display.genimagesc(state_escape,gH,gSyn);
    colormap(h1,flipud(gray(256))); h=colorbar('YTick',0:20:100); ylabel(h,'% HCOs')
      %h=colorbar('Position', [0.955  0.1  0.01  0.2]);
    title('% of HCOs for each (gH, gSyn)');  ylabel('gH');
    caxis([0 100])
    
    case 2
        display.genimagesc(FRescmean,gH,gSyn);
        myColorMap = colormaps.viridis; myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Cycle frequency, Hz'),
        title('Cycle frequency');
        %caxis([min(FRescmean(:))-0.02 max(FRescmean(:))+0.02])
        caxis([minFR-0.02 maxFR+0.02])
        
    case 3
        display.genimagesc(Aescmean,gH, gSyn);
        myColorMap = colormaps.magma; myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Amplitude, mV'),
        title('Slow-wave amplitude');
        %caxis([min(Aescmean(:))-0.2 max(Aescmean(:))+0.2])
        caxis([minA-0.2 maxA+0.2])        
        
    case 4
        display.genimagesc(Nspksescmean,gH,gSyn);
        myColorMap = summer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'# spikes/burst'),
        title('# spikes/burst');
        caxis([min(Nspksescmean(:))-0.2 max(Nspksescmean(:))+0.2])  
        caxis([minNspks-0.3 maxNspks+0.3])                
        
    case 5
        display.genimagesc(SpkFRescmean,gH,gSyn);
        myColorMap = spring(256); myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Spike Frequency, Hz'),
        title('Spike Frequency');
        %caxis([min(SpkFRescmean(:))-0.1 max(SpkFRescmean(:))+0.1])
        caxis([minSpkFR-0.2 maxSpkFR+0.2])           
        
    case 6
        display.genimagesc(DCescmean,gH,gSyn);
        myColorMap = copper(256); myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Duty cycle, %'),
        title('Duty cycle');
        %caxis([min(DCescmean(:))-0.1 max(DCescmean(:))+0.1]) 
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

% release
for i=1:6
    h1=display.bigsubplot(2,6,2,i,g,l);
    switch i
        case 1
    eogenimagesc(state_release,gH,gSyn);
    colormap(h1,flipud(gray(256))); h=colorbar('YTick',0:20:100); ylabel(h,'% HCOs')
      %h=colorbar('Position', [0.955  0.1  0.01  0.2]);
    title('% of HCOs for each (gH, gSyn)');  ylabel('gH');
    
    case 2
        eogenimagesc(FRrelmean,gH,gSyn);
        myColorMap = colormaps.viridis; myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Cycle frequency, Hz'),
        title('Cycle frequency');
        %caxis([min(FRrelmean(:))-0.02 max(FRrelmean(:))+0.02])
        caxis([minFR-0.02 maxFR+0.02])
        
    case 3
        eogenimagesc(Arelmean,gH,gSyn);
        myColorMap = colormaps.magma; myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Amplitude, mV'),
        title('Slow-wave amplitude');
        %caxis([min(Arelmean(:))-0.2 max(Arelmean(:))+0.2])
        caxis([minA-0.2 maxA+0.2])          
        
    case 4
        eogenimagesc(Nspksrelmean,gH,gSyn);
        myColorMap = summer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'# spikes/burst'),
        title('# spikes/burst');
        %caxis([min(Nspksrelmean(:))-0.2 max(Nspksrelmean(:))+0.2])  
        caxis([minNspks-0.3 maxNspks+0.3])         
        
    case 5
        eogenimagesc(SpkFRrelmean,gH,gSyn);
        myColorMap = spring(256); myColorMap(1,:) = 1; colormap(h1,myColorMap); 
        h=colorbar; ylabel(h,'Spike Frequency, Hz'),
        title('Spike Frequency');
        %caxis([min(SpkFRrelmean(:))-0.1 max(SpkFRrelmean(:))+0.1])
        caxis([minSpkFR-0.2 maxSpkFR+0.2])           
        
    case 6
        eogenimagesc(DCrelmean,gH,gSyn);
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
