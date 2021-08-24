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
    
    k=k+1;
end

%%
ie=find(data.mode=="escape");
ir=find(data.mode=="release");
iq1=find(data.Q10==11);
iq2g=find(data.Q10==21);
iq2=find(data.Q10==22);
