function bursts=detect_bursts(st,start_ISI,stop_ISI)

% Inputs:
%   st              vector with spike times (sec)
%   start_ISI       inter-spike interval to start a burst (sec)
%   stop_ISI        inter-spike interval to terminate a burst (sec)
%   min_nspikes     minimum number of spikes in a burst
%
% Output:
%   burst        struct with the following fields:
%     .N         burst number
%     .nSp       number of spikes in each burst
%     .firstSp   burst onsets
%     .lastSp    burst offsets
%     .center    burst center
%     .BuDur     burst durations
%     .BuFreq    instanteneous inter-burst frequency
%     .SpFreq    mean frequency of spikes within bursts


% defaults
if nargin<2, start_ISI=mean(diff(st)); end
if nargin<3, stop_ISI=mean(diff(st))+0.5; end
%if nargin<4, min_nspks=2; end

if isrow(st); st=st'; end

in_burst=false;
k=0; % burst count

for i=1:numel(st)
    
    if ~in_burst && i~=nspikes && st(i+1)-st(i)<start_ISI
        
        in_burst=true;
        
        k=k+1; % count bursts
  
        first_spike=st(i); % first spike in a burst
        burst_onsets(k)=first_spike; % store the onset of this burst
        
        nspks=2;
        
    elseif in_burst && i~=nspikes && st(i+1)-st(i)<stop_ISI
        
        nspks=nspks+1; % count number of spikes in a burst
        
    else
        
        last_spike=st(i); % last spike in a burst
        
        burst_durs(k)=last_spike-first_spike; %burst duration
        Nspks(k)=nspks; % tota number of spikes per burst
        
        in_burst=false;
    end
end

% output
bursts=struct();
bursts.N=(1:k)';
bursts.nSp=Nspks;
bursts.firstSp=burst_onsets;
bursts.lastSp=burst_onsets+burst_durs;
bursts.center=(bursts.lastSp+bursts.firstSp)./2;
bursts.BuDur=burst_durs;
bursts.BuFreq=1./(diff(bursts.center));
bursts.SpFreq=bursts.nSp./ bursts.BuDur;