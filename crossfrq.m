function [cxy,coxy,coan,errorcoan,f] = crossfrq(V1,V2,nfft,Fs,overlap)
% CROSSPOWER used to calculated the spectrum of cross power, coherence, cross phase and
% errorbar of cross phase.
% SYNTAX
% [cxy,coxy,coan,errorcoan,f] = CROSSPOWER(V1,V2,nfft,Fs,overlaps)
%
% INPUT PARAMETERS
% A: mesh in patch format with the following additional subfields
% - A.FA:
% - A.FN:
%
% OUTPUT PARAMETERS
% output: output description
%
% DESCRIPTION
% Exhaustive and long description of functionality of crosspower.
%
% Examples:
% description of example for crosspower
% >> [output] = crosspower(A);
%
% See also:
% csd, cohere
%
% References:
% [1] First coded by Adi Liu in USTC

% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 17-Jul-2015 12:25:50


if nargin<3, nfft=1000;end
if nargin<4,Fs=1e6; end
if nargin<5,overlap=nfft/2;end
warning off
cxy=csd(V1,V2,nfft,Fs,hanning(nfft),overlap,'linear');
[coxy,f]=cohere(V1,V2,nfft,Fs,hanning(nfft),overlap,'linear');
coan=angle(cxy);

m=floor((length(V1)-overlap)/(nfft-overlap));
errorcoan=0.5*sqrt((1-coxy.^2)./(m*coxy.^2));
errorcoan=errorcoan/pi;

    figure('position',[ 117   0   560   642])
    subplot(3,1,1)
    %semilogx(f,abs(cxy),'linewidth',2.5)
    plot(f,abs(cxy),'linewidth',2.5)
   % turny;
   h_axes=findobj(gcf,'type','axes');
   set(h_axes,'LineWidth',1.5);
    set(h_axes,'fontsize',14)
    set(h_axes,'fontweight','bold')
    set(h_axes,'XMinorTick','on')
    set(h_axes,'YMinorTick','on')
    set(h_axes,'layer','top')
    ylabel2('Cross Power (a.u.)')

%     ylim([-40 0])
    myplot(3,1,1)
set(gca,'fontsize',13)
    
    subplot(3,1,2)
    %semilogx(f,sqrt(coxy),'linewidth',2.5);
    plot(f,sqrt(coxy),'linewidth',2.5);
    set(h_axes,'LineWidth',1.5);
    set(h_axes,'fontsize',14)
    set(h_axes,'fontweight','bold')
    set(h_axes,'XMinorTick','on')
    set(h_axes,'YMinorTick','on')
    set(h_axes,'layer','top')
    ylabel2('Coherence') ;
%     ylim([0 0.4])
    myplot(3,1,2)
set(gca,'fontsize',13)
    
     line([Fs/nfft Fs/2],[sqrt(1/m) sqrt(1/m)],'color',[0.502 0.502 0.502],'linestyle','--','linewidth',1.5)
     
    subplot(3,1,3)
    %semilogx(f,coan,'linestyle','none','marker','o','markersize',6,'linewidth',2.5,'markerfacecolor','w');
    %semilogx(f,coan,'linewidth',2.5);
    plot(f,coan,'*');
    set(h_axes,'LineWidth',1.5);
    set(h_axes,'fontsize',14)
    set(h_axes,'fontweight','bold')
    set(h_axes,'XMinorTick','on')
    set(h_axes,'YMinorTick','on')
    set(h_axes,'layer','top')
    ylabel2('Cross Phase (rad)');
    ylim([-pi,pi]);
    myplot(3,1,3)
    set(gca,'fontsize',13)
    xtick2
    xlabel2('f (kHz)')
    turnx
