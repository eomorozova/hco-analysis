function ax = plot_barerror(data,err, varargin)
%% plot_scatter(data, groupIdx, pos, color, symbol, opt)
% Plot boxplot and scatter overlaid figure
% INPUT:
%       - data, <vector>, M * 1, M number of total points.
% OPTIONAL INPUTS:
%       - color, <cellarray>, cell, N * 1, scatter color of each group
%                         Default: random

%       - opt, <scalar>, scatter along center of box (0) or fullfill box (1)
%                         Default: 1

%% Input parsing
    color = {'r','g','b','y','m','c'};

    n = numel(data);
    
    p = inputParser;
    % Required input
    addRequired(p,'data');
    % Optional input
    addOptional(p,'color', color(randi(6,1, n)));
    addOptional(p,'opt',1);
    parse(p,data, varargin{:})
%%
data = p.Results.data;
color = p.Results.color;

for i = 1 : n
    datai = data(i);
    h = bar(i, datai,0.7);
    h.FaceColor = color{i};
    h.EdgeColor = 'k';
    h.FaceAlpha = 0.5
    h.LineWidth = 1.5
    hold on;
end

er = errorbar([1:n],data,err); 
er.Color ='k';
er.LineWidth = 1.5;
er.LineStyle = 'none';  
