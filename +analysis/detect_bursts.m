function bursts=detect_bursts(st,start_ISI,stop_ISI, min_nspks)

% Inputs:
%   st              vector with spike times (sec)
%   start_ISI       inter-spike interval to start a burst (sec)
%   stop_ISI        inter-spike interval to terminate a burst (sec)
%   min_nspks       minimum number of spikes in a burst
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
if nargin<4, min_nspks=2; end

if isrow(st); st=st'; end

Nspks=[]; burst_on=[]; burst_off=[]; burst_durs=[]; 

in_burst=false;
k=0; % burst count

for i=1:numel(st)
    
    if ~in_burst && i~=numel(st) && st(i+1)-st(i)<start_ISI
        
        in_burst=true;
        
        k=k+1; % count bursts
        
        first_spike=st(i); % first spike in a burst
        burst_on(k)=first_spike; % store the onset of this burst

        nspks=2;
        
    elseif in_burst && i~=numel(st) && st(i+1)-st(i)<stop_ISI
        
        nspks=nspks+1; % count number of spikes in a burst

    else
        % check if # spikes/burst is more than minimum specified number
        if in_burst
            if nspks>=min_nspks
                last_spike=st(i); % last spike in a burst
                burst_off(k)=last_spike;
                
                burst_durs(k)=last_spike-first_spike; %burst duration
                Nspks(k)=nspks; % tota number of spikes per burst
            end
        end
            in_burst=false;
        end
    end

    if isempty(burst_off); burst_on=[]; k=0; end

% output
bursts=struct();
bursts.N=(1:k)';
bursts.nSp=Nspks;
bursts.firstSp=burst_on;
bursts.lastSp=burst_off;
bursts.center=(bursts.lastSp+bursts.firstSp)./2;
bursts.BuDur=burst_durs;
bursts.BuFreq=1./(diff(bursts.center));
bursts.SpFreq=bursts.nSp./ bursts.BuDur;