
% path='qcflineargrd54_';
% p250=import_dmp(strcat(path,'n250/data'), 'P'); disp('p250');
% p200=import_dmp(strcat(path,'n200/data'), 'P'); disp('p200');
% p150=import_dmp(strcat(path,'n150/data'), 'P'); disp('p150');
% p100=import_dmp(strcat(path,'n100/data'), 'P'); disp('p100');
% p80=import_dmp(strcat(path,'n80/data'), 'P'); disp('p80');
% p75=import_dmp(strcat(path,'n75/data'), 'P'); disp('p75');
% p70=import_dmp(strcat(path,'n70/data'), 'P'); disp('p70');
% p65=import_dmp(strcat(path,'n65/data'), 'P'); disp('p65');
% p60=import_dmp(strcat(path,'n60/data'), 'P'); disp('p60');
% p55=import_dmp(strcat(path,'n55/data'), 'P'); disp('p55');
% p50=import_dmp(strcat(path,'n50/data'), 'P'); disp('p50');
% p45=import_dmp(strcat(path,'n45/data'), 'P'); disp('p45');
% p40=import_dmp(strcat(path,'n40/data'), 'P'); disp('p40');
% p35=import_dmp(strcat(path,'n35/data'), 'P'); disp('p35');
% p30=import_dmp(strcat(path,'n30/data'), 'P'); disp('p30');
% p25=import_dmp(strcat(path,'n25/data'), 'P'); disp('p25');
% p20=import_dmp(strcat(path,'n20/data'), 'P'); disp('p20');
% p15=import_dmp(strcat(path,'n15/data'), 'P'); disp('p15');
% p10=import_dmp(strcat(path,'n10/data'), 'P'); disp('p10');
% p05=import_dmp(strcat(path,'n5/data'), 'P');  disp('p05');


% pm05= moment(p05, 'rms');
% pm10= moment(p10, 'rms');
% pm15= moment(p15, 'rms');
% pm20= moment(p20, 'rms');
% pm25= moment(p25, 'rms');
% pm30= moment(p30, 'rms');
% pm35= moment(p35, 'rms');
% pm40= moment(p40, 'rms');
% pm45= moment(p45, 'rms');
% pm50= moment(p50, 'rms');
% pm55= moment(p55, 'rms');
% pm60= moment(p60, 'rms');
% pm65= moment(p65, 'rms');
% pm70= moment(p70, 'rms');
% pm75= moment(p75, 'rms');
% pm80= moment(p80, 'rms');

% [x,y,t]=size(pm);
% pmrt=zeros(x,t);
% for ii=1:t
%   
%        pmrt(:,ii)=pm(:,40,ii)/max(pm(:,40,ii));
%     
% end

% figure;
% for i=1:3:151
%     plot_pol_slice(p05(:,:,:,i),'cbm18_dens6_ne5_0.1_1.4.grid.nc',5);
%     text_time = ['t = ' num2str(i)];
%     title2(text_time);
%     M(i)=getframe(gcf);
% end
% movie2avi(M,'filename.avi');    

p0=p12phi0;
num=245;
Power=zeros(num,32);
% Power09=zeros(350,32);
% Power05=zeros(350,32);
for i=1:num
    Pxy0=zeros(1,32);
    for j=25:40
        x=squeeze(p0(300,j,1:64,i));
        Pxy = fft(x);
        Pxy = abs(Pxy(1:32));
        Pxy0=Pxy'+Pxy0;
    end
    Power(i,:)=Pxy0/16;
end
% for i=1:350
%     Pxy0=zeros(1,32);
%     for j=35:45
%         x=squeeze(p09(300,j,1:64,i));
%         Pxy = fft(x);
%         Pxy = abs(Pxy(1:32));
%         Pxy0=Pxy'+Pxy0;
%     end
%     Power09(i,:)=Pxy0/11;
% end
% for i=1:350
%     Pxy0=zeros(1,32);
%     for j=35:45
%         x=squeeze(p05(300,j,1:64,i));
%         Pxy = fft(x);
%         Pxy = abs(Pxy(1:32));
%         Pxy0=Pxy'+Pxy0;
%     end
%     Power05(i,:)=Pxy0/11;
% end

%y=squeeze(pm20(300,32,1:350));

bisum=zeros(32,num);
for i=1:num
    x=squeeze(p0(300,20:44,1:65,i));
    %maxp=max(max(x));
    maxp=0.01;
    xd=x+abs(maxp)*10e-2*rand([25,65]);
    %xd=x;
    [bisum(:,i),n0] = bispecmode(xd,xd,xd,65,65*5,[1,32],0,1);
    
end
