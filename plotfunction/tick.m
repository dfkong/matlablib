function tick(stat)

h_axes=findobj(gcf,'type','axes');

    for i=1:length(h_axes)
        set(h_axes(i),'XMinorTick',stat);
        set(h_axes(i),'YMinorTick',stat);
    end