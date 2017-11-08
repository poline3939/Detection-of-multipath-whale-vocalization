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
Up_11=Up_1(1:L/2,:);
s_mea=sum(Up_11,1)/(L/2);

%%% spectrogram of signal+noise and pure noise
ln1=zeros(1,L/2);
ln0=zeros(1,L/2);    
p_h0_f=ln1; p_h0_d=ln1; p_h1_f=ln1; p_h1_d=ln1;

% cat1=[]; cat0=[];
for ii=1:L/2
    X2=randn(1,size(Up_11,2));
    X2=X2/norm(X2);
    X1=X2+Up_1(ii+L/2,:);
    
    xx1=X1*s_mea';
    xx0=X2*s_mea';
%     p_h0_f(ii)=exp(-1/2*X2*X2');   
%     p_h0_d(ii)=exp(-1/2*X1*X1'); 
%     p_h1_f(ii)=exp(-1/2*(X2-s_mea)*(X2-s_mea)');
%     p_h1_d(ii)=exp(-1/2*(X1-s_mea)*(X1-s_mea)');
    ln1(ii)=xx1-1/2;
    ln0(ii)=xx0-1/2;    
end

inter=(max(max(ln1),max(ln0))-min(min(ln1),min(ln0)))/1000;
beta=min(min(ln1),min(ln0)):inter:max(max(ln1),max(ln0));
count=histc(ln1,beta);
n_p1=count/sum(count);
c_e1=cumsum(n_p1);
count2=histc(ln0,beta);
n_p2=count2/sum(count2);
c_e2=cumsum(n_p2);
inter2=inter/10;

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
KK2_1=vv; %time-STFT ROC
KK2_2=ww;
auc_n2=RocIntegral(KK2_1,KK2_2,length(KK2_1)-1); % 0.7677
save('roc_time_STFT_realdata.mat','KK2_1','KK2_2','auc_n2')

