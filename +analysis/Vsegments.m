function [Vseg,Vth1]=Vsegments(V,Vth)

% Inputs:
% V     membrane potential recordings
% Vth   synaptic threhold
%
% Outputs:
% Vseg  segments of memrane potential with a constant synaptic threshold
% Vth1 median synaptic synaptic threshold 

Vth_change=find(abs(diff(Vth))>1);
n=numel(Vth_change)+1; % number of synaptic threholds

    Vth1(1)=median(Vth(1:Vth_change(1)));
    Vseg{1}=V(1:Vth_change(1));

for i=2:numel(Vth_change)
    Vth1(i)=median(Vth(Vth_change(i-1)+1:Vth_change(i)));
    Vseg{i}=V(Vth_change(i-1)+1:Vth_change(i));
end

    Vth1(n)=median(Vth(Vth_change(end):length(Vth)));
    Vseg{n}=V(Vth_change(end):length(Vth));
