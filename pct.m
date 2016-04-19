function [pct,pctall,pxy] = pct(pre,phi,nfft,Fs,n)
% PCT: Calculate the phase coherence time of a mode;
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
% Exhaustive and long description of functionality of pct.
%
% Examples:
% description of example for pct
% >> [output] = pct(A);
%
% See also:
% cpsd, crosspsd
%
% References:
% [1] P.W. Xi, X.Q. Xu and P. H. Diamond, Phase Dynamics Criterion for Fast 
%     Relaxation of High-Confinement-Mode Plasmas
%     PRL 112, 085001 (2014).

% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 27-Aug-2015 15:18:59

% write down your codes from here.
pre=squeeze(pre);
phi=squeeze(phi);
sizep=size(pre);
if length(sizep)>4&&length(sizep)==1
    error('pre and phi should be 2D, 3D or 4D array!')
end
if nargin<3
    nfft=sizep(end-1);
end
if nargin<4
    Fs=nfft;
end
if nargin<5
    n=20;
end

pxy=zeros(floor((nfft+1)/2),sizep(end));
poxy=zeros(floor((nfft+1)/2),sizep(end));
coan=zeros(floor((nfft+1)/2),sizep(end));
if length(sizep)==2
    for ii=1:sizep(end)
        pret=squeeze(pre(:,ii));
        phit=squeeze(phi(:,ii));
        [pxy(:,ii),poxy(:,ii),coan(:,ii),nall] = crosspsd(pret,phit,nfft,Fs,0);
    end
end
if length(sizep)==3
    for ii=1:sizep(end)
        pret=squeeze(pre(:,:,ii));
        phit=squeeze(phi(:,:,ii));
        [pxy(:,ii),poxy(:,ii),coan(:,ii),nall] = crosspsd(pret,phit,nfft,Fs,0);
    end
end
if length(sizep)==4
    for ii=1:sizep(end)
        pret=squeeze(pre(:,:,:,ii));
        phit=squeeze(phi(:,:,:,ii));
        [pxy(:,ii),poxy(:,ii),coan(:,ii),nall] = crosspsd(pret,phit,nfft,Fs,0);
    end
end
[~,nindex]=min(abs(nall-n));
pct=coan(nindex,:);
pctall=coan;

if nargout==0
    t=(0:length(pct)-1)*5;  %%%%%%dependent on the timestep in BOUT++
    figure;
    subplot(2,1,1)
    plot(t,abs(pxy(nindex,:)),'linewidth',2.0);
    set(gca,'yscale','log');
    set(gca,'LineWidth',1.5);
    set(gca,'fontsize',14)
    set(gca,'fontweight','bold')
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    set(gca,'layer','top')
    ylabel3('$Cross Power A_{\tilde{P}\tilde{\phi}} (a.u.)$')
    myplot(2,1,1)
    
    subplot(2,1,2)
    plot(t,pct,'linewidth',2.0);
    set(gca,'LineWidth',1.5);
    set(gca,'fontsize',14)
    set(gca,'fontweight','bold')
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    set(gca,'layer','top')
    xlabel2('t/\tau_A');
    ylabel3('$\delta \phi_{\tilde{P}\tilde{\phi}} (rad)$');
    ylim([-pi,pi]);
    myplot(2,1,2)
    set(gca,'fontsize',13)        
end
