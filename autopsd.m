function [pxx,n] = autopsd(x,nfft,Fs,overlap)
% AUTOPSD: used to calculate the auto-power spectrum in f (x is 1D) or k (x is 2D) space.
% SYNTAX
% [output] = AUTOPSD(A)
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
% Exhaustive and long description of functionality of autopsd.
%
% Examples:
% description of example for autopsd
% >> [output] = autopsd(A);
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
% $Revision: 1.0$ Created on: 16-Aug-2015 21:58:01

% write down your codes from here.
sizex=size(x);
if length(sizex)>2
    error('x should be 1D or 2D array');
end
if sizex(1)*sizex(2) < size(1)+size(2)  %x is 1D
    [pxx,n] = pwelch(x,hanning(nfft),overlap,nfft,Fs);
elseif sizex(1)*sizex(2) > size(1)+size(2)
    xtemp=reshape(x',sizex(1)*sizex(2),1);
    overlap=0;
    [pxx,n] = pwelch(xtemp,hanning(nfft),overlap,nfft,Fs);
end

if nargout==0  
    figure;plot(n(2:end),pxx(2:end));
    if sizex(1)*sizex(2) < size(1)+size(2)
        xlabel2('Frequency (Hz)');
    elseif sizex(1)*sizex(2) > size(1)+size(2)
        xlabel2('Mode number');
    end
    ylabel2('Auto-power density (a.u.)');
    set(gca, 'XMinorTick','on');set(gca, 'YMinorTick','on');
    turny
end
