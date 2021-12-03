%% all the map data
addpath('C:\Users\moroz\Desktop\HCO_git')
addpath('C:\Users\moroz\Documents\Katya_DynamicClamp\matlab-tools-master')
addpath(genpath('C:\Users\moroz\Documents\Katya_DynamicClamp\hco_related_codes'))

folder=[]; ifiles=[];
for i=1:7
    switch i
        case 1
            folder{i}= 'D:\Katya_DynamicClamp\HCO_data\951_060\';
            ifiles{i}=[14:18] % [-30, -35, -40, -45, -50]
        case 2
            folder{i}= 'D:\Katya_DynamicClamp\HCO_data\951_036\';
            ifiles{i}=[8:13] %  [-35, -40,-45,-50]
    end
end

%% load data (h5)
%data=[]; param=[]; key=[]; datalength=[];
for j=1:2
    for i = ifiles{j}
        file_search = strcat(folder{j}, '\', '*.h5');  files=dir(file_search);
        [param{j}{i},data{j}{i},key{j}{i}] = rtxih5load([folder{j},files(i).name]);
    end
    datalength{j}(i)=length(data{j}{i});
    Fs(j)=1./param{j}{ifiles{j}(1)}.dt;
end

%% plot the data
j=1; i=ifiles{j}(1); clf;
x=1/Fs(j)/60:1/Fs(j)/60:length(data{j}{i}(1,:))/Fs(j)/60;
for ii=1:size(data{j}{i},1)
subplot(size(data{j}{i},1),1,ii)
plot(x,data{j}{i}(ii,:),'k','linewidth',1)
end
%% assign variables
%V1=[]; V2=[]; I1=[]; I2=[]; gh=[]; gsyn=[]; time=[]; gmi=[];
for j= 1% experiment
    k=1;
    for l=ifiles{j}
        V1{j}{k}=data{j}{l}(3,1:end); V2{j}{k}=data{j}{l}(4,1:end);
        I1{j}{k}=data{j}{l}(1,1:end); I2{j}{k}=data{j}{l}(2,1:end);

        try; Vth{j}(k)=param{j}{l}.V_half_12_HCO_parameter_loop; end
        try; Vth{j}(k)=param{j}{l}.V_half_11_HCO_parameter_loop; end
        
        gh{j}{k}=data{j}{l}(6,1:end); gsyn{j}{k}=data{j}{l}(7,1:end);    
        time{j}{k}=data{j}{l}(9,1:end);
        k=k+1;
    end
end

%% find gH, gSyn values and break into segments
%gH=[]; gSyn=[]; gMI=[]; time_start=[]; Vhco=[]; Ihco=[];
for j=1
    for k=1:length(ifiles{j})
        [time1, time_start{j}{k}]=find(time{j}{k}==0.001);
        time_start{j}{k}=time_start{j}{k}(1:1:end);
        
        for ii=1:length(time_start{j}{k})    
            gH{j}{k}(ii)=mean(gh{j}{k}(time_start{j}{k}(ii):time_start{j}{k}(ii)+40*Fs(j))); gH{j}{k}(ii)=round(gH{j}{k}(ii)*10^4);
            gSyn{j}{k}(ii)=mean(gsyn{j}{k}(time_start{j}{k}(ii):time_start{j}{k}(ii)+40*Fs(j))); gSyn{j}{k}(ii)=round(gSyn{j}{k}(ii)*10^4);
            
            Vhco1{j}{k}{ii}=V1{j}{k}(time_start{j}{k}(ii):time_start{j}{k}(ii)+40*Fs(j))*10^2; % mV
            Vhco2{j}{k}{ii}=V2{j}{k}(time_start{j}{k}(ii):time_start{j}{k}(ii)+40*Fs(j))*10^2; % mV
            
            Ihco1{j}{k}{ii}=I1{j}{k}(time_start{j}{k}(ii):time_start{j}{k}(ii)+40*Fs(j));
            Ihco2{j}{k}{ii}=I2{j}{k}(time_start{j}{k}(ii):time_start{j}{k}(ii)+40*Fs(j));
            
            Vhco{j}{k}{ii}=[Vhco1{j}{k}{ii}; Vhco2{j}{k}{ii}];
            Ihco{j}{k}{ii}=[Ihco1{j}{k}{ii}; Ihco2{j}{k}{ii}];
        end
    end
end

%% here save data to structure

data=[];

j=1;

data{j}.file = '951_060';
data{j}.Vth = [-30, -35, -40, -45, -50];
data{j}.V = Vhco{j};
data{j}.I = Ihco{j};
data{j}.gH = gH{j};
data{j}.gSyn = gSyn{j};
data{j}.Fs = Fs(j);

%% save structure

save('C:\Users\moroz\Documents\code\hco_analysis\data\data_fig4.mat','data','-v7.3')

%% calculate half-center output characteristcs as a function of gH and gSyn

hcostat=[]; ERQ=[];

for jj = 1:numel(data) % experiment
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

for jj = 1:numel(data) % experiment
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
                else
                    FR{jj}{j}(ii,i)=NaN; A{jj}{j}(ii,i)=NaN;
                    SpkFR{jj}{j}(ii,i)=NaN; Nspk{jj}{j}(ii,i)=NaN;
                    DC{jj}{j}(ii,i)=NaN; dc{jj}{j}(ii,i)=NaN;
                end
            end
        end
    end
end

%% classify the network state (silent, spiking, half-center,..)

networkstate=[]; state=[]; percentSingleSpikeBursts=[];

for jj=1:numel(data) % experiment
    
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

%% if state is not HCO,  set all the measures to 0)

for jj=1 % experiment
    for j=1:numel(state{jj}) % map number
        for ii=1:2 % neurons number
            for i=1:length(networkstate{jj}{j})
                
                if networkstate{jj}{j}(i)~=5 % if not a half-center
                    ERQ{jj}{j}(ii,i)=NaN;
                    FR{jj}{j}(ii,i)=NaN;
                    A{jj}{j}(ii,i)=NaN;
                    SpkFR{jj}{j}(ii,i)=NaN;
                    Nspk{jj}{j}(ii,i)=NaN;
                    DC{jj}{j}(ii,i)=NaN;
                    dc{jj}{j}(ii,i)=NaN;
                end
            end
        end
    end
end

%% calculate mean characteristics across two neurons in a circuit

for jj = 1:numel(data)
    for j = 1:numel(data{jj}.V)
        ERQall{jj}(j,:) = mean(ERQ{jj}{j});
        FRall{jj}(j,:) = mean(FR{jj}{j}); 
        Aall{jj}(j,:) = mean(A{jj}{j}); 
        SpkFRall{jj}(j,:) = mean(SpkFR{jj}{j}); 
        Nspkall{jj}(j,:) = mean(Nspk{jj}{j}); 
        DCall{jj}(j,:) = mean(DC{jj}{j}); 
        dcall{jj}(j,:) = mean(dc{jj}{j});         
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

%% reshpe

FRmean=[]; Amean=[]; SpkFRmean=[]; Nspkmean=[]; DCmean=[]; dcmean=[]; ERQ1=[];

for jj = 1:numel(data)
    for j=1:numel(data{jj}.V)
        
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
    end
end

%% Figure 3 C-E

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
        xlabel('gSyn'); axis square
        set(gca,'Fontsize',14,'FontName','Arial');
        
        for i=1:10 % create grid
            line([0 1200], [150*i+75 150*i+75],'color',[0.5 0.5 0.5])
            line([150*i+75 150*i+75],[0 1200],'color',[0.5 0.5 0.5])
        end
    end
end



%% plot all the maps

jj=1;
gH = data{jj}.gH{1}; gSyn = data{jj}.gSyn{1};

g=0.04; l=0.07;
clf,

for j=1:5 % map 
    for i=1:6 % number of features to display
        
        h1=display.bigsubplot(6,5,i,j,g,l);
        
        switch i
            case 1
                display.genimagesc(state{jj}{j},gH,gSyn);
                myColorMap = [0.3,0.33,0.35; 0.25,0.4,0.3; 0,0.4,0.5; 0.8,0.5,0.3; 0.6,0.2,0;];
                colormap(h1,myColorMap);
                h=colorbar;
                h.Ticks=[1:5]; h.TickLabels={'silent','asymmetric','irregular\newlinespiking','antiphase\newlinespiking','HCO'}; caxis([1 5]);
                title('State');  ylabel('gH');
                axis square
                
                %title(['Vth=', num2str(Vth{j}(k)),'mV'])
                
            case 2
                display.genimagesc(FRmean{jj}{j},gH,gSyn);
                myColorMap = colormaps.viridis; myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Cycle frequency, Hz'),
                title('Cycle frequency');
                %caxis([min(FRescmean(:))-0.02 max(FRescmean(:))+0.02])
                caxis([minFR-0.02 maxFR+0.02])
                
            case 3
                display.genimagesc(Amean{jj}{j},gH, gSyn);
                myColorMap = colormaps.magma; myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Amplitude, mV'),
                title('Slow-wave amplitude');
                caxis([minA-0.2 maxA+0.2])
                
            case 4
                display.genimagesc(Nspkmean{jj}{j},gH,gSyn);
                myColorMap = summer(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'# spikes/burst'),
                title('# spikes/burst');
                %caxis([min(Nspksmean(:))-0.2 max(Nspksescmean(:))+0.2])
                caxis([minNspk-0.3 maxNspk+0.3])
                
            case 5
                display.genimagesc(SpkFRmean{jj}{j},gH,gSyn);
                myColorMap = spring(256); myColorMap(1,:) = 1; colormap(h1,myColorMap);
                h=colorbar; ylabel(h,'Spike Frequency, Hz'),
                title('Spike Frequency');
                %caxis([minSpkFR-0.2 maxSpkFR+0.2])
                caxis([minSpkFR-0.2 25])
                
            case 6
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


