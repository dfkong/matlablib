function elmsize = elm_size(dcp,p0,uedge,xmin,xmax,yind,Bbar)
% ELM_SIZE used to calculate the loss energy during the elm crash. 
% The function is just translated from the "elm_size.pro" idl routine in
% bout/tools/idllib.
% SYNTAX
% [output] = ELM_SIZE(dcp,p0,uedge,xmin,xmax,yind,Bbar)
%
% INPUT PARAMETERS
% - dcp: 3D array from moment.m
% - p0: 2D pressure profile, can be collected from BOUT output
% - uedge: a struct data contents the profile and grid information in nc
% file
%
% OUTPUT PARAMETERS
% elmsize: a struct type
% elmsize.eloss: the energy loss calculated from the change of pressure
% profile
%
% DESCRIPTION
% Exhaustive and long description of functionality of elm_size.
%
% Examples:
% description of example for elm_size
% >> elmsize = elm_size(A);
%
% See also:
% moment, import_dmp, readgridnc
%
% References:
% [1] elm_size.pro in BOUT/tools/idllib/elm_size.pro

% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 01-Jul-2015 15:25:38

% write down your codes from here.
%narginchk(4,7);

if nargin < 5, 
    xmin = 1; 
end
if nargin < 6, 
    xmax = 328; 
end
if nargin < 7, 
    yind = 32;  %choose the poloidal location for 1D size
end
if nargin < 8, 
    Bbar = 1.992982; %the normalized magnetic field
end

mydcp=dcp;
myp0=p0;
g=uedge;
MU0=4.0e-7*pi;

s=size(mydcp);

if length(s)~=3
    error('dcp should be 3D (x,y,t)');
end
nx=s(1);
ny=s(2);
nt=s(3);

Dtheta=g.dy;     %using correct poloidal angle
psixy=g.psixy;
R=g.Rxy;
Bp=g.Bpxy;
hthe=g.hthe;

Dpsi=zeros(nx,ny);
Dpsi(1,:)=psixy(2,:)-psixy(1,:);
Dpsi(nx,:)=psixy(nx,:)-psixy(nx-1,:);
for i=2:nx-1
  Dpsi(i,:)=(psixy(i+1,:)-psixy(i-1,:))/2;
end

Ddcp1=zeros(1,nt);
Ddcp2=zeros(1,nt);
Ddcp3=zeros(1,nt);
% Tp01=0;
% Tp02=0;
% Tp03=0;

for t=1:nt
 Ddcp3(t)=2.0*pi*sum(sum(mydcp(xmin:xmax,:,t).*hthe(xmin:xmax,:).*Dtheta(xmin:xmax,:).*Dpsi(xmin:xmax,:)./Bp(xmin:xmax,:)));
 Ddcp2(t)=sum(sum(mydcp(xmin:xmax,:,t).*hthe(xmin:xmax,:).*Dtheta(xmin:xmax,:).*Dpsi(xmin:xmax,:)./(R(xmin:xmax,:).*Bp(xmin:xmax,:))));
 Ddcp1(t)=sum(sum(mydcp(xmin:xmax,yind,t).*Dpsi(xmin:xmax,yind)./(R(xmin:xmax,yind).*Bp(xmin:xmax,yind)))); 
end

Tp03=2.0*pi*sum(sum(myp0(xmin:xmax,:).*hthe(xmin:xmax,:).*Dtheta(xmin:xmax,:).*Dpsi(xmin:xmax,:)./Bp(xmin:xmax,:)));
Tp02=sum(sum(myp0(xmin:xmax,:).*hthe(xmin:xmax,:).*Dtheta(xmin:xmax,:).*Dpsi(xmin:xmax,:)./(R(xmin:xmax,:).*Bp(xmin:xmax,:))));
Tp01=sum(sum(myp0(xmin:xmax,yind).*Dpsi(xmin:xmax,yind)./(R(xmin:xmax,yind).*Bp(xmin:xmax,yind))));

% s1=zeros(1,nt);
% s2=zeros(1,nt);
% s3=zeros(1,nt);
% E_loss=zeros(1,nt);

s1=-Ddcp1/Tp01;   %1D elm size
s2=-Ddcp2/Tp02;   %2D elm size
s3=-Ddcp3/Tp03;   %3D elm size

E_loss=-Ddcp3*(0.5*Bbar*Bbar/MU0);    %energy loss, unit J
E_total=Tp03*(0.5*Bbar*Bbar/MU0);     %total energy, unit J

elmsize=struct('s1',s1,'s2',s2,'s3',s3,'eloss',E_loss,'etotal',E_total);