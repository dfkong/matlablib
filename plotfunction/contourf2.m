function contourf2(varargin)
figure;

[temp,h]=contourf(varargin{:});
 h_axes=findobj(gcf,'type','axes');
    set(h_axes,'LineWidth',1.5);
    set(h_axes,'fontsize',14)
    set(h_axes,'fontweight','bold')
    set(h_axes,'XMinorTick','on')
    set(h_axes,'YMinorTick','on')
    set(h_axes,'layer','top')
%     xlabel2('f(kHz)')
%     ylabel2('PSD(a.u.)')

set(h,'linestyle','none');
colorbar4;