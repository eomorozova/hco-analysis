function h=plothorzline(y,dastring,linewidth)

dalim=get(gca,'XLim');
if nargin<2
    h=line(dalim,[y y],'Color',[0 1 0],'LineStyle','--','linewidth',3);
elseif nargin<3
    h=line(dalim,[y y],'Color',dastring);
else
    h=line(dalim,[y y],'Color',dastring,'linewidth',linewidth);
end
