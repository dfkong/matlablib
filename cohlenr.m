function [cohnx,psi,n,r,coxy] = cohlenr(p,uedge,position,nnum,tindex,nfft,maxgrid)
% COHLENR short description
% SYNTAX
% [output] = COHLENR(A)
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
% Exhaustive and long description of functionality of cohlenr.
%
% Examples:
% description of example for cohlenr
% >> [output] = cohlenr(A);
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
% $Revision: 1.0$ Created on: 16-Aug-2015 14:25:57

% write down your codes from here.
xindex=position(1); %
yindex=position(2);
%zindex=position(3);
g=uedge;


if xindex>=280&&xindex<=340 %xgrids=516
    xin=xindex-80;
    xout=xindex+80;
elseif xindex>=155&&xindex<=195 %xgrids=256
    xin=xindex-25;
    xout=xindex+25;
end


%sizep=size(p);
coxy=zeros(length(xin:xout),floor(nfft/2)+1);
for ii=xin:xout
    arr1=squeeze(p(ii,(yindex-10):(yindex+10),:,tindex));
    arr2=squeeze(p(xindex,(yindex-10):(yindex+10),:,tindex));
    sizearr1=size(arr1);
    sizearr2=size(arr2);
    arr1=reshape(arr1',sizearr1(1)*sizearr1(2),1);
    arr2=reshape(arr2',sizearr2(1)*sizearr2(2),1);
    arr1=arr1+max(abs(arr1)*5e-2)*rand(sizearr1(1)*sizearr1(2),1);
    arr2=arr2+max(abs(arr2)*5e-2)*rand(sizearr2(1)*sizearr2(2),1);
    %[coxy((ii-xin+1),:),n]=mscohere(arr1,arr2,hanning(nfft),0,nfft,maxgrid); %coherence
    [coxy((ii-xin+1),:),n]=cpsd(arr1,arr2,hanning(nfft),0,nfft,maxgrid); %cross power
end
coxy=abs(coxy);
[~,nindex]=min(abs(n-nnum));
cohnx=coxy(:,nindex);

r=g.Rxy(xin:xout,32);
psi=(g.psixy(xin:xout,32)-g.psi_axis)/(g.psi_bndry-g.psi_axis);

if nargout==0
    figure;plot(psi,cohnx);
    xlabel2('\Psi_N');
    ylabel2('Coherence: \gamma^2');
    set(gca, 'XMinorTick','on');set(gca, 'YMinorTick','on');
end

