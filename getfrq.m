function [fval] = getfrq(pzt,fs,passfrq)
% GETFRQ: Used to get a temporal frequncy (fval) of a coherent mode. The result is not
% that reliable, please double check with the fft result. If you want to
% make sure this result reliable, the input pzt should be a Instrinsic Mode
% Functions (IMF) in Hilbert-Huang Transform (HHT).
% SYNTAX
% [output] = GETFRQ(A)
%
% INPUT PARAMETERS
% pzt: pt=squeeze(p05(300,32,:,:));
% fs: sample frequency of simulations
% passfrq: if you want to filter the pzt 
%
% DESCRIPTION
% Exhaustive and long description of functionality of getfrq.
%
% Examples:
% description of example for getfrq
% >> pt=squeeze(p(300,32,:,:));
% >> [fval] = getfrq(pt);
%
% See also:
% flt_b, filter_b,cheb2ord,
% 

% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 28-Jul-2015 14:25:17

% write down your codes from here.
if nargin==1
    fs=1/(5*3.38*1e-7); 
    x=hilbert(pzt);
    theta=angle(x);
    fval=abs(fs*diff(theta)/(2*pi));
    fval=mean(fval);
elseif nargin==2; 
    passfrq=[5e3,200e3];
    temp=flt_b(pzt',passfrq,20,1,fs);
    x=hilbert(temp');
    theta=angle(x);
    fval=abs(fs*diff(theta)/(2*pi));
    fval=mean(fval);
elseif nargin==3    
    temp=flt_b(pzt',passfrq,20,1,fs);
    x=hilbert(temp');
    theta=angle(x);
    fval=abs(fs*diff(theta)/(2*pi));
    fval=mean(fval);
end


function [Y,a,b,W]=flt_b(X,passfre,Rs,Rp,Fs)
% cheb2ord      
% cheby2

if nargin<3,Rs=40;end
if nargin<4,Rp=1;end
if nargin<5,Fs=1e6;end

[b,a,W]=filter_b(passfre,Rs,Rp,Fs);

Y = filtfilt(b, a, X);
Y = filtfilt(b, a, Y);

function [b,a,W]=filter_b(frerange,Rs,Rp,Fs)
%bandpass

if nargin<2
    Rs=50;
end

if nargin<3
    Rp=1;
end

if nargin<4
    Fs=1e6;
end

Fs=Fs/2; %Nyquist frequency
Wp=[frerange(1)/Fs  frerange(2)/Fs];

step=0;
flag=0;
while flag==0&&(frerange(1)-step-0.5e3)>0
    step=step+0.5e3;
    Ws=[(frerange(1)-step)/Fs  (frerange(2)+step)/Fs];
    [n,ws]=cheb2ord(Wp,Ws,Rp,Rs,'s');
    [b,a]=cheby2(n,Rs,ws);   
    hd=dfilt.df2(b,a);
    flag=isstable(hd);
end

if flag==0
    step=0;
    while flag==0
        step=step+0.5e3;
        Wp=[(frerange(1)+step)/Fs  (frerange(2)+step)/Fs];
         [n,ws]=cheb2ord(Wp,Ws,Rp,Rs);
         [b,a]=cheby2(n,Rs,ws);
         hd=dfilt.df2(b,a);
         flag=isstable(hd);
    end
end

W=[Ws(1) Wp(1) Wp(2) Ws(2)]*Fs/1e3;       %kHz

