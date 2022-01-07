
function state=classify_networkstate(hcostat,Fs)

% classify the activity of the network into silent/asymmetric/irregular/half-center/antiphase

%inputs:
% V         membrane potential
% hcostat   characteristics of the circuit output calculated by analysis.hco_stat

% output:
%state      classified into silent/asymmetric/irregular/half-center/antiphase

if ~isempty(hcostat{1})
    
    st1 = hcostat{1}.st/Fs; st2 = hcostat{2}.st/Fs; % spike times in s
    
    percentSingleSpikeBursts = calcPercentSingleSpikeBursts(st1,st2);
    
    % both neurons spike less than 5 times in a minute
    if length(st1(st1<60))<=5  && length(st2(st2<60))<=5
        
        state = 1; % silent
        
    % one neuron spikes, the other one is silent (<5 spikes/minute)
    elseif (length(st1(st1<60))<5  && length(st2(st2<60))>5 || ...
            length(st1(st1<60))>5  && length(st2(st2<60))<5)
        
        state = 2; % asymmetric
        
    % both neurons spike more than 5 spikes per minute
    elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 ...
            && isnan(hcostat{1}.T1_mean) && isnan(hcostat{2}.T1_mean))
        
        % if more than 80% spikes are alternating
        if percentSingleSpikeBursts>0.8
            
            state = 4; % antiphase spiking
            
        else
            
            state = 3; % irregular spiking
            
        end
        
    % alternating bursting pattern of activity
    elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 && ...
            ~isnan(hcostat{1}.T1_mean) && ~isnan(hcostat{2}.T1_mean))
        
        state = 5; % half-center oscillator
        
    % single spike half-center oscillator
    elseif (length(st1(st1<60))>5  && length(st2(st2<60))>5 && ...
            ~isnan(hcostat{1}.T1_mean) && ~isnan(hcostat{2}.T1_mean) && ...
            isnan(hcostat{1}.nSpks_mean) && isnan(hcostat{2}.nSpks_mean))
        
        state = 4; % antiphase spiking
        
    else
        
        state = 6; % misclassified
        
    end
else
    
    state = 1; %silent
    
end

