clear
clc
%% mean ocean spectrogram distribution detector
load updata
f=256;
w=256;
h=128;
fs=2000;

Up=XC';
n=13201;
L=size(Up,1);
Up1=zeros(size(Up));
for ii=1:L
    Up1(ii,:)=Up(ii,:)/norm(Up(ii,:));
end

load hp_6th.mat
for ii=1:L
    Up_1(ii,:)=filter(hp_6th,1,Up1(ii,:));
end
[A1,F1,T1,~]=spectrogram(Up_1(ii,:),256,128,256,2000,'yaxis');
figure
imagesc(T1,F1, abs(A1).^2)

for ii=1:L/2
    [A1,F1,T1,~]=spectrogram(Up_1(ii,:),256,128,256,2000,'yaxis');
    S1=abs(A1).^2;
    S1k(ii,:)=S1(:);    
end
m_s=sum(S1k,1)/size(S1k,1); % mean of spectrogram elements "training"

Gaus1=zeros(size(Up));
for ii=1:L/2
    Gaus1(ii,:)=randn(1,13201);
    Gaus1(ii,:)=Gaus1(ii,:)/norm(Gaus1(ii,:));
    [C1,F3,T3,~]=spectrogram(Gaus1(ii,:),256,128,256,2000,'yaxis');
    [B1,F2,T2,~]=spectrogram(Up_1(ii+L/2,:),256,128,256,2000,'yaxis');
    S2=abs(B1).^2;
    S2k(ii,:)=S2(:);  
    S3=abs(C1).^2;
    S3k(ii,:)=S3(:);   

    zz1=2*sqrt(m_s.*S2k(ii,:))/f;
    lambda1=besseli(0,zz1);
    ln1=log(lambda1)-m_s/f;
    ln1_1(ii)=sum(ln1);
    
    zz0=2*sqrt(m_s.*S3k(ii,:))/f;
    lambda0=besseli(0,zz0);
    ln0=log(lambda0)-m_s/f;
    ln0_1(ii)=sum(ln0);
end    

inter=(max(max(ln0_1),max(ln1_1))-min(min(ln0_1),min(ln1_1)))/1000;
inter2=inter/10;
beta=min(min(ln0_1),min(ln1_1)):inter:max(max(ln0_1),max(ln1_1));
count=histc(ln1_1,beta);
n_p1=count/sum(count);
c_e1=cumsum(n_p1);

count2=histc(ln0_1,beta);
n_p2=count2/sum(count2);
c_e2=cumsum(n_p2);
 
vc=-abs(min(beta)):inter2:(max(beta)-inter);
t=1;
for j=1:length(vc)
    while(beta(t)<vc(j))
          t=t+1;
    end
    ww(j)=1-c_e1(t); %pd
    vv(j)=1-c_e2(t); %pf
    t=1;
end
KK1_1=vv; % spectrogram distribution
KK1_2=ww;
auc_n1=RocIntegral(KK1_2,KK1_1,length(KK1_1)-1); % 0.6039
save('roc_spec_dist_realdata.mat','KK1_1','KK1_2','auc_n1')

