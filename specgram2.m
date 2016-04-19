function [s,T,F]=specgram2(sign,range,nfft,Fs,overlap)
% specgram2��Ҫ����matlab�Դ��?pecgram����Ļ��ϼ򵥵��޸���һ���������
% specgram���������������ʱ����Ҷ�任�ġ�?
if nargin<2, range=[0,1500];end
if nargin<3, nfft=1000;end
if nargin<4, Fs=1e6;end
if nargin<5, overlap=floor(1*nfft/2);end

[s,f,t1]=specgram(sign,nfft,Fs,hanning(nfft),overlap);
t=linspace(range(1),range(2),length(t1));
if sum(imag(sign))~=0
    f=[f(nfft/2+2:end)-Fs; f(1:nfft/2+1)];
    s=[s(nfft/2+2:end,:); s(1:nfft/2+1,:)];
end

[T,F]=meshgrid(t,f/1000);
s=abs(s);

if nargout==0
    if sum(imag(sign))~=0
        index=[2200:nfft/2-2,nfft/2+4:nfft-2200];
        %figure;pcolor(T(10:101,:),F(10:101,:),abs(s(10:101,:))); shading interp;
        figure;imagesc(t,f(index)/1e3,log(abs(s(index,:)))); 
        %figure;pcolor(T(3:end,:),F(3:end,:),abs(s(3:end,:)));shading interp;
        %contourf2(T(10:end,:),F(10:end,:),abs(s(10:end,:)),30);
        %contourf2(T(16:126,:),F(16:126,:),10*log(abs(s(16:126,:))),50);
        %contourf2(T(2:end,:),F(2:end,:),10*log(abs(s(2:end,:))),50);
    else
        index=3:nfft/2;
        %figure;pcolor(T(index,:),F(index,:),abs(s(index,:))); shading interp;
        %figure;pcolor(T(index,:),F(index,:),log(abs(s(index,:)))); shading interp;
        %figure;pcolor(T(3:end,:),F(3:end,:),abs(s(3:end,:)));shading interp;
        %contourf2(T(index,:),F(index,:),abs(s(index,:)),30);
        contourf2(T(index,:),F(index,:),10*log(abs(s(index,:))),30);
    end
    h_axes=findobj(gcf,'type','axes');
    set(h_axes,'LineWidth',1);
    set(h_axes,'fontsize',14)
    set(h_axes,'fontweight','bold')
    set(h_axes,'XMinorTick','on')
    set(h_axes,'YMinorTick','on')
    set(h_axes,'layer','top');
    xlabel2('t (s)');
    ylabel2('f (kHz)');
end
