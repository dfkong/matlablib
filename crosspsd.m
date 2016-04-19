function [pxy,poxy,coan,n] = crosspsd(x,y,nfft,Fs,overlap)
% CROSSPSD used to calculate the cross-power spectrum in f (x and y is 1D) or k (x and y is 2D) space.
% SYNTAX
% [output] = CROSSPSD(A)
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
% Exhaustive and long description of functionality of crosspsd.
%
% Examples:
% description of example for crosspsd
% >> [output] s crosspsd(A);
%
% See also:
% ALSO1, ALSO2
%
% References:
% [1] A. Einstein, Die Grundlage der allgemeinen Relativitaetstheorie,
%     Annalen der Physik 49, 769-822 (1916).

% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 19-Aug-2015 13:32:52

% write down your codes from here.
x=squeeze(x);
y=squeeze(y);
sizex=size(x);
sizey=size(y);
if sizex~=sizey
    error('x and y should be of the same length!');
end

if length(sizex)==2 && ((sizex(1)*sizex(2))<(sizex(1)+sizex(2))) %x is 1D array x(t)
    [pxy,n] = cpsd(x,y,hanning(nfft),overlap,nfft,Fs);
    poxy = mscohere(x,y,hanning(nfft),overlap,nfft,Fs);
    coan=angle(pxy); 
    m=floor((length(x)-overlap)/(nfft-overlap));
elseif length(sizex)==2 && ((sizex(1)*sizex(2)) > (sizex(1)+sizex(2))) %x is 2D array x(y,z)
    xtemp=reshape(x',sizex(1)*sizex(2),1);
    ytemp=reshape(y',sizex(1)*sizex(2),1);
    overlap=0;
    [pxy,n] = cpsd(xtemp,ytemp,hanning(nfft),overlap,nfft,Fs);
    poxy = mscohere(xtemp,ytemp,hanning(nfft),overlap,nfft,Fs);
    coan=angle(pxy);
    m=sizex(1);
elseif length(sizex)==3    %x is 3D array x(x,y,z)
    pxy=zeros(floor((nfft+1)/2),1);
    poxy=zeros(floor((nfft+1)/2),1);
    for ii=1:sizex(1)
        x1=squeeze(x(ii,:,:));
        y1=squeeze(y(ii,:,:));
        xtemp=reshape(x1',sizex(2)*sizex(3),1);
        ytemp=reshape(y1',sizex(2)*sizex(3),1);
        overlap=0;
        [pxyt,n] = cpsd(xtemp,ytemp,hanning(nfft),overlap,nfft,Fs);
        poxyt = mscohere(xtemp,ytemp,hanning(nfft),overlap,nfft,Fs);
        pxy=pxy+pxyt;
        poxy=poxy+poxyt;
    end
    m=sizex(1)*size(2);
    pxy=pxy/m;
    coan=angle(pxy);
end

% plot
if nargout==0
    figure('position',[ 117   0   560   642])
    subplot(3,1,1)
    plot(n,abs(pxy),'linewidth',2.0)
    set(gca,'LineWidth',1.5);
    set(gca,'fontsize',14)
    set(gca,'fontweight','bold')
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    set(gca,'layer','top')
    ylabel2('Cross Power (a.u.)')
    myplot(3,1,1)
    
    subplot(3,1,2)
    plot(n,sqrt(poxy),'linewidth',2.0);
    set(gca,'LineWidth',1.5);
    set(gca,'fontsize',14)
    set(gca,'fontweight','bold')
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    set(gca,'layer','top')
    ylabel2('Coherence') ;
    myplot(3,1,2)
    line([0 Fs/2],[sqrt(1/m) sqrt(1/m)],'color',[0.502 0.502 0.502],'linestyle','--','linewidth',2.0)
     
    subplot(3,1,3)
    plot(n,coan,'o-','linewidth',2.0);
    set(gca,'LineWidth',1.5);
    set(gca,'fontsize',14)
    set(gca,'fontweight','bold')
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    set(gca,'layer','top')
    ylabel2('Cross Phase (rad)');
    ylim([-pi,pi]);
    myplot(3,1,3)
    set(gca,'fontsize',13)
    if length(sizex)==2 && sizex(1)*sizex(2)<sizex(1)+sizex(2) 
        xlabel2('f (Hz)');
    elseif (length(sizex)==2 && sizex(1)*sizex(2) > sizex(1)+sizex(2))||length(sizex)==3  
        xlabel2('Mode Number');
    end
end
