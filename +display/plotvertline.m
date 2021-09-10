function h=plotvertline(x,dastring,linewidth)

dalim=get(gca,'YLim');
if nargin<2
    h=line([x x], dalim, 'Color','k','LineStyle','--','linewidth',1);
elseif nargin<3
    h=line([x x], dalim, 'Color',dastring);
else
    h=line([x x], dalim, 'Color',dastring,'linewidth',linewidth);
end
