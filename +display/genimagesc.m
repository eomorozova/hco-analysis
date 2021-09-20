function h=genimagesc(x,xx,yy,dogray,docolormapNan,docolorbar)
%% imagesc with sane y direction
%genimagesc(x,dogray,docolormapNan,docolorbar)
if isempty(whos('dogray'))
    dogray=0;
end
if isempty(whos('docolormapNan'))
    docolormapNan=0;
end
if isempty(whos('docolorbar'))
    docolorbar=0;
end

if isempty(whos('dogray'))
    dogray=0;
end
if docolormapNan
    minv = min(min(x));
    maxv = max(max(x));
    x(isnan(x)) = minv-((maxv-minv)/5);
    ddd=[0 0 0;jet(10)];
    colormap(ddd);
% else
%     colormap(jet);
end

if (isempty(whos('xx')) || isempty(whos('yy')))
  h=imagesc(x);
else
    h=imagesc(xx,yy,x);
    set(gca,'YDir','normal')
end

    

set(gca,'YDir','normal')
if dogray
    colormap(gcf,gray)
end
if docolorbar
    colorbar
end