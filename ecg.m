clc;clear;
[hdr,record]=edfRead('r01.edf'); %��ȡECG�ź�
FECG0=record(1,1:5000); %ԭʼ��FECG�ź�
MECG0=record(2,1:5000);%ԭʼ��abdominal ECG�ź�

%FECGͼ��
figure(1)
plot(FECG0,'b')
detrend_FECG=detrend(FECG0);  %ȥ����
d_record1=FECG0-detrend_FECG;
hold on
plot(detrend_FECG,'r')
fmaxd=3;
fs=1000;
fmaxn=fmaxd/(fs/2);
[b a]=butter(1,fmaxn,'low');
F_FECG=filtfilt(b,a,detrend_FECG);  %ȥƯ��
FECG=detrend_FECG-F_FECG;
hold on
plot(FECG,'green')
l1=legend('original FECG','detrended FECG','filted FECG')
title('Ԥ�����ECG�ź�','FontSize',20)
xlabel('length','FontSize',20);
ylabel('voltage(mV)','FontSize',20);
set(gcf,'unit','centimeters','position',[1 3 30 15]);
set(gca,'Position',[.1 .2 .8 .6]);
set(l1,'Fontname', 'Times New Roman','FontWeight','bold','FontSize',20)

%abdominal ECG ͼ��
figure(2)
plot(MECG0,'b')
detrend_MECG=detrend(MECG0);  %ȥ����
hold on
plot(detrend_MECG,'r')
F_MECG=filtfilt(b,a,detrend_MECG);  %ȥƯ��
MECG=detrend_MECG-F_MECG;
hold on
plot(MECG,'green')
l2=legend('original abdominal ECG','detrended abdominal ECG','filted abdominal ECG')
title('Ԥ�����abdominal ECG�ź�','FontSize',20)
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
title('abdominal ECG�ź����ȡ�Ĳ���abdominal ECG�Ļ����','FontSize',20)

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
axis([500 5000 -100 50]);title('abdominal ECG�ź���MECGģ��Ա�','FontSize',20)
subplot(2, 1, 2); plot(diffECG); axis([500 5000 -50 50]);title('ȥ��ĸ���ĵ��̥���ĵ�','FontSize',20)




%  clc;clear;
% [hdr,record]=edfRead('r01.edf'); %��ȡECG�ź�
% FECG0=record(1,1:5000);
% MECG0=record(2,1:5000);
% plot(FECG0,'b')
% detrend_FECG=detrend(FECG0);  %ȥ����
% d_record1=FECG0-detrend_FECG;
% hold on
% plot(detrend_FECG,'r')
% [b a]=cheby1(4,0.5,0.3/180,'high');  %����б�ѩ��I���˲���
% FECG=filter(b,a,detrend_FECG);  %ȥƯ��
% hold on
% plot(FECG,'green')
% legend('ԭʼ��ECG�ź�','ȥ���ƺ��ECG','ȥƯ�ƺ��ECG')
% title('Ԥ�����ECG�ź�')
% set(gcf,'unit','centimeters','position',[1 3 37 15]);
% set(gca,'Position',[.1 .2 .8 .6]);
% 
% figure (2)
% plot(MECG0,'b')
% detrend_MECG=detrend(MECG0);  %ȥ����
% d_record2=MECG0-detrend_MECG;
% hold on
% plot(detrend_MECG,'r')
% MECG=filter(b,a,detrend_MECG);  %ȥƯ��
% hold on
% plot(MECG,'green')
% legend('ĸ�׵�MECG�ź�','ȥ���ƺ��MECG','ȥƯ�ƺ��MECG')
% title('Ԥ�����MECG�ź�')
% set(gcf,'unit','centimeters','position',[1 3 37 15]);
% set(gca,'Position',[.1 .2 .8 .6]);


% clc;
% clear all;
% close all;
% %loading input signal%
% load foetal_ecg.dat 
% t = foetal_ecg(:,1);%ʱ��
% abdominal = foetal_ecg(:,2:6); %�����ź�
% thoraic = -foetal_ecg(:,7:9) ; %��ǻ�ź�
% avg_abdominal= mean(abdominal,2); %�����źŵ�ƽ��ֵ
% avg_thoraic= mean(thoraic,2); %��ǻ�źŵ�ƽ��ֵ
% reference=avg_thoraic;  %��Ϊ�ο�����
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
% title('ĸ�׸����ź�')
% xlabel('ʱ��/s');
% ylabel('���/mV');
% subplot(212);%
% plot(t,reference,'r');
% title('ĸ����ǻ�ź�')
% xlabel('ʱ��/s');
% ylabel('���/mV');
% 
% % ʹ��LLMS
% nord3 = 20; % �˲����Ľ���
% gamma = 0.0005; %LLMS��gamma��ֵ
% ss3 = 0.000000075; % ����
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
% [bbw abw]=cheby1(4,0.5,0.3/180,'high');  %����б�ѩ��I���˲���
% detrend_CB1=filter(bbw,abw,detrend_CB);  %ȥƯ��
% plot(t,detrend_CB1,'r')
% axis([1 2 -40 40])

