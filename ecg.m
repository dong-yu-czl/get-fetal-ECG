clc;clear;
[hdr,record]=edfRead('r01.edf'); %读取ECG信号
FECG0=record(1,1:5000); %原始的FECG信号
MECG0=record(2,1:5000);%原始的abdominal ECG信号

%FECG图像
figure(1)
plot(FECG0,'b')
detrend_FECG=detrend(FECG0);  %去趋势
d_record1=FECG0-detrend_FECG;
hold on
plot(detrend_FECG,'r')
fmaxd=3;
fs=1000;
fmaxn=fmaxd/(fs/2);
[b a]=butter(1,fmaxn,'low');
F_FECG=filtfilt(b,a,detrend_FECG);  %去漂移
FECG=detrend_FECG-F_FECG;
hold on
plot(FECG,'green')
l1=legend('original FECG','detrended FECG','filted FECG')
title('预处理的ECG信号','FontSize',20)
xlabel('length','FontSize',20);
ylabel('voltage(mV)','FontSize',20);
set(gcf,'unit','centimeters','position',[1 3 30 15]);
set(gca,'Position',[.1 .2 .8 .6]);
set(l1,'Fontname', 'Times New Roman','FontWeight','bold','FontSize',20)

%abdominal ECG 图像
figure(2)
plot(MECG0,'b')
detrend_MECG=detrend(MECG0);  %去趋势
hold on
plot(detrend_MECG,'r')
F_MECG=filtfilt(b,a,detrend_MECG);  %去漂移
MECG=detrend_MECG-F_MECG;
hold on
plot(MECG,'green')
l2=legend('original abdominal ECG','detrended abdominal ECG','filted abdominal ECG')
title('预处理的abdominal ECG信号','FontSize',20)
set(gcf,'unit','centimeters','position',[1 3 30 15]);
set(gca,'Position',[.1 .2 .8 .6]);
xlabel('length','FontSize',20);
ylabel('voltage(mV)','FontSize',20);
set(l2,'Fontname', 'Times New Roman','FontWeight','bold','FontSize',20)
axis([0 5000 -80 100])

ecgData = MECG;
ecgOneBeat = ecgData(2051:2200);  
ecgTemplate = ecgOneBeat;       
ecgOneBeat(length(ecgOneBeat)+1:length(ecgData)) = 0; 
figure(3)
[acor lag]=xcorr(MECG,ecgOneBeat,'coeff');
plot(acor)
set(gcf,'unit','centimeters','position',[1 3 30 15]);
set(gca,'Position',[.1 .2 .8 .6]);

lagROI = lag(-lag(1)+1:length(lag));
acorROI = acor(-lag(1)+1:length(lag));

figure(4)
plot(lagROI, acorROI, '-');
set(gcf,'unit','centimeters','position',[1 3 30 15]);
set(gca,'Position',[.1 .2 .8 .6]);
title('abdominal ECG信号与截取的部分abdominal ECG的互相关','FontSize',20)

meanROI = mean(acorROI);
stdROI = std(acorROI);
overTH = find(acorROI > 3*stdROI);
overTH = [overTH, 1000000];
diffTH = diff(overTH);
ind = find(diffTH > 2);

for iInd = 1:length(ind)
    rangeXCOR = 0;
    if iInd == 1
        rangeXCOR = overTH(1:ind(iInd));
    else
        rangeXCOR = overTH(ind(iInd-1) + 1:ind(iInd)); %
    end
    [v iTmp] = max(acorROI(rangeXCOR)); 
    indECGmatched(iInd) = iTmp + rangeXCOR(1) - 2;
    
end

% hold on ;plot(indECGmatched, 0.2* ones(1, length(indECGmatched)), '*');

ecgTemplateReplicate = zeros(1, length(ecgData));
eTlen = length(ecgTemplate);
for i = 1:length(indECGmatched)
    ecgTemplateReplicate(indECGmatched(i)+1:indECGmatched(i)+eTlen) = ecgTemplate;
end


diffECG = ecgData - ecgTemplateReplicate;
figure(5)
set(gcf,'unit','centimeters','position',[1 3 30 15]);
set(gca,'Position',[.1 .2 .8 .6]);
subplot(2, 1, 1); plot(ecgData); hold on; plot(ecgTemplateReplicate-50, 'r-');
axis([500 5000 -100 50]);title('abdominal ECG信号与MECG模板对比','FontSize',20)
subplot(2, 1, 2); plot(diffECG); axis([500 5000 -50 50]);title('去除母亲心电的胎儿心电','FontSize',20)




%  clc;clear;
% [hdr,record]=edfRead('r01.edf'); %读取ECG信号
% FECG0=record(1,1:5000);
% MECG0=record(2,1:5000);
% plot(FECG0,'b')
% detrend_FECG=detrend(FECG0);  %去趋势
% d_record1=FECG0-detrend_FECG;
% hold on
% plot(detrend_FECG,'r')
% [b a]=cheby1(4,0.5,0.3/180,'high');  %设计切比雪夫I型滤波器
% FECG=filter(b,a,detrend_FECG);  %去漂移
% hold on
% plot(FECG,'green')
% legend('原始的ECG信号','去趋势后的ECG','去漂移后的ECG')
% title('预处理的ECG信号')
% set(gcf,'unit','centimeters','position',[1 3 37 15]);
% set(gca,'Position',[.1 .2 .8 .6]);
% 
% figure (2)
% plot(MECG0,'b')
% detrend_MECG=detrend(MECG0);  %去趋势
% d_record2=MECG0-detrend_MECG;
% hold on
% plot(detrend_MECG,'r')
% MECG=filter(b,a,detrend_MECG);  %去漂移
% hold on
% plot(MECG,'green')
% legend('母亲的MECG信号','去趋势后的MECG','去漂移后的MECG')
% title('预处理的MECG信号')
% set(gcf,'unit','centimeters','position',[1 3 37 15]);
% set(gca,'Position',[.1 .2 .8 .6]);


% clc;
% clear all;
% close all;
% %loading input signal%
% load foetal_ecg.dat 
% t = foetal_ecg(:,1);%时间
% abdominal = foetal_ecg(:,2:6); %腹部信号
% thoraic = -foetal_ecg(:,7:9) ; %胸腔信号
% avg_abdominal= mean(abdominal,2); %腹部信号的平均值
% avg_thoraic= mean(thoraic,2); %胸腔信号的平均值
% reference=avg_thoraic;  %作为参考输入
% 
% % C=xcorr(avg_abdominal,reference);
% % plot(C)
% % axis([0 2500 -400000 100000])
% 
% 
% 
% figure
% subplot(211);
% plot(t,avg_abdominal);
% title('母亲腹部信号')
% xlabel('时间/s');
% ylabel('振幅/mV');
% subplot(212);%
% plot(t,reference,'r');
% title('母亲胸腔信号')
% xlabel('时间/s');
% ylabel('振幅/mV');
% 
% % 使用LLMS
% nord3 = 20; % 滤波器的阶数
% gamma = 0.0005; %LLMS中gamma的值
% ss3 = 0.000000075; % 步长
% [W3,Child_E3,Maternal_Y3] = llms(reference,avg_abdominal,ss3,gamma,nord3);
% figure
% plot(t,reference,'--k');
% hold on
% plot(t,Child_E3,'r')
% title('SISO System Filter Error & Thoraic L-LMS');
% xlabel('Time [Sec]');
% ylabel('Amplitude [mV]');
% legend('Mother HB','Child HB');
% figure(3)
% detrend_CB=detrend(Child_E3);
% [bbw abw]=cheby1(4,0.5,0.3/180,'high');  %设计切比雪夫I型滤波器
% detrend_CB1=filter(bbw,abw,detrend_CB);  %去漂移
% plot(t,detrend_CB1,'r')
% axis([1 2 -40 40])

