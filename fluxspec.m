function [flx,n] = fluxspec(x,y,nfft,maxgrid,overlap)
% FLUXSPEC: This function is used to calculated the flux spectrum at one
% moment. If x is density fluctuation and y is potential fluctuation, then
% the output flx is particle flux. If x is tempareture fluctuation, then
% flx is heat flux. If x is pressure fluctuation, then flx is energy flux?
% SYNTAX
% [flx,n] = FLUXSPEC(x,y,nfft,maxgrid,overlap)
%
% INPUT PARAMETERS
% - x: 2D density fluctuation ne(y,z); 2D tempareture fluctuation te(y,z);
%       2D pressure fluctuation Pe(y,z);
% - y: 2D potential fluctuation phi(y,z), y is used to calculate the radial
%       velocity fluctuations as vr=(dphi/dy)/B. Here B is magnetic field.
% - nfft: corresponding MZ (number of points in z direction (2^n + 1)) in
%           BOUT.inp
% - maxgrid: The total maxgrid in z direction, maxgrid=MZ*ZPERIOD;
% - overlap: overlap=0;
%
% OUTPUT PARAMETERS
% - flx: flux spectrum
% - n: toroidal mode number index
%
% Examples:
% description of example for fluxspec
% >> [output] = fluxspec(A);
%
% See also:
% ALSO1, ALSO2
%
% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 13-Jul-2015 12:19:15

% write down your codes from here.

narginchk(4,6);
indsize=size(x);

if length(indsize)~=2,
error('x should be a 2D array with m=poloidal and n=torodial');
end

if nargin < 3, 
    nfft = indsize(2); 
end
if nargin < 4, 
    maxgrid = nfft; 
    warning('ZPERIOD=1, #Fraction of a torus to simulate in BOUT.inp');
end
if nargin < 5, 
    overlap = 0; 
    warning('In most of the time , overlap should be set 0');
end
Fs=maxgrid;
if mod(indsize(2),2)==1;
    x=x(:,1:end-1);
    y=y(:,1:end-1);
    Fs=maxgrid-maxgrid/nfft;
    nfft=nfft-1;  
end

xind=size(x);
x0=reshape(x',xind(1)*xind(2),1);
temp=fftmat(divmat(x0,nfft,hanning(nfft),overlap));
xfft=temp(1:nfft/2+1,:);

yind=size(y);
y0=reshape(y',yind(1)*yind(2),1);
temp=fftmat(divmat(y0,nfft,hanning(nfft),overlap));
yfft=temp(1:nfft/2+1,:);

% temp=conj(xfft).*yfft;
% temp=2*real(temp);
% temp=temp(1:nfft/2,:);
% flx=mean(temp,2);
n=linspace(0,Fs/2,nfft/2+1);
index=size(xfft);
Kmat=-n(ones(index(2),1),:);

temp=conj(xfft).*Kmat'.*yfft;
%temp=2*real(temp);
temp=2*imag(temp);
temp=temp(1:nfft/2+1,:);
flx=mean(temp,2);
% temp=sum(temp);
% tau_sum=mean(temp);


