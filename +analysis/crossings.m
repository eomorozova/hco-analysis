function [stup, stdown] = crossings(V,Vth)

if ~exist('Vth')
    Vth=-35;
end

for i=1:size(V,1)
    crossings(i,:)=(V(i,:)>Vth);
    di(i,:)=diff(crossings(i,:));
    stup{i}=(find(di(i,:)>0));
    stdown{i}=(find(di(i,:)<0));
end
