function ax = plot_scatter(data, groupIdx, varargin)
%% plot_scatter(data, groupIdx, pos, color, symbol, opt)

% INPUTS:
%       - data, <vector>, M * 1, M number of total points.
%       - groupIdx, <vector>, M * 1, different groups of data

% OPTIONAL INPUTS:
%       - pos, <vector>, position of each boxplot
%                         Default: 1 : N, N number of groups
%       - color, <cellarray>, cell, N * 1, scatter color of each group
%                         Default: random
%       - symbol, <cellarray>, cell, N * 1, scatter symbol of each group
%                         Default: random
%       - markersize, 1, Default: 50
%       - opt, <scalar>, scatter along center of box (0) or fullfill box (1)
%                         Default: 1


color = {'r','g','b','y','m','c'};
symbol = {'o','s','d','^','v','<','>','p','h'};
markersize = 50;
n = numel(unique(groupIdx));

p = inputParser;
% Required input
addRequired(p,'data');
addRequired(p,'groupIdx');
% Optional input
addOptional(p,'pos', 1:n);
addOptional(p,'color', color(randi(6,1, n)));
addOptional(p,'symbol',symbol(randi(7,1, n)));
addOptional(p,'markersize',markersize);
addOptional(p,'opt',1);
parse(p,data, groupIdx, varargin{:})


data = p.Results.data;
group = p.Results.groupIdx;
color = p.Results.color;
symbol = p.Results.symbol;
markersize = p.Results.markersize;

for i = 1 : n
    position = p.Results.pos(i);
    datai = data(group == i);
    if p.Results.opt
        s = scatter(ones(size(datai)).*(position +(rand(size(datai))-0.5)/4),datai,markersize);
    else
        s = scatter(ones(size(datai))* position ,datai,markersize);
    end
    s.MarkerFaceColor = color{i};
    s.MarkerEdgeColor = color{i};
    s.MarkerEdgeColor = 'w';    
    s.Marker = symbol{i};
    
    %line([i-0.3 i+0.3],[median(datai) median(datai)],'color','k','linewidth',2)
    hold on;
end