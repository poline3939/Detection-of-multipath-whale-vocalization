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

s_mea=zeros(1,5500);
s_mea(1:length(s_m))=s_m;

ssu=su(1:50,:);
Su=zeros(size(ssu));
for ii=1:size(ssu,1)
    Su(ii,:)=ssu(ii,:)/sqrt(ssu(ii,:)*ssu(ii,:)');
end
S_mea=s_mea/sqrt(s_mea*s_mea');
dd=length(s_mea);

for kk=1:4
    es=Es(kk);
    %%% spectrogram of signal+noise and pure noise
    clear Ln1 Ln0
    Ln1=[];
    Ln0=[];
    
    Su1=sqrt(es)*Su;
    S_m=sqrt(es)*S_mea;
%     cat1=[]; cat0=[];
    for j=1:500
%         sk1=0; sk0=0;
        X2=randn(1,dd);
        Sk1=zeros(size(Su1)); Sk0=Sk1;
        for qq=1:size(Su1,1)
            X1=X2+Su1(qq,:);
            Sk1(qq,:)=X1*S_m'-es/2;
            Sk0(qq,:)=X2*S_m'-es/2;
        end
        t_sk1=Sk1';
        t_sk0=Sk0';
        Ln1=[Ln1;t_sk1(:)]; % signal+noise
        Ln0=[Ln0;t_sk0(:)]; % pure noise
    end
    inter=(max(max(Ln1),max(Ln0))-min(min(Ln1),min(Ln0)))/1000;
    beta=min(min(Ln1),min(Ln0)):inter:max(max(Ln1),max(Ln0));
    count=histc(Ln1,beta);
    n_p1=count/sum(count);
    c_e1=cumsum(n_p1);
    count2=histc(Ln0,beta);
    n_p2=count2/sum(count2);
    c_e2=cumsum(n_p2);
    switch kk
        case 1
            KK1_1=1-c_e2; %pf
            KK1_2=1-c_e1; %pd
        case 2
            KK2_1=1-c_e2;
            KK2_2=1-c_e1;
        case 3
            KK3_1=1-c_e2;
            KK3_2=1-c_e1;
        case 4
            KK4_1=1-c_e2;
            KK4_2=1-c_e1;
    end

end

N=100;
ii=round(linspace(1,numel(KK1_1),N));

figure
plot(KK1_1(ii),KK1_2(ii),'k^-')
hold on
plot(KK2_1(ii),KK2_2(ii),'k*-',KK3_1(ii),KK3_2(ii),'k+-',KK4_1(ii),KK4_2(ii),'ko-')
axis([0 1 0 1])
roc1_try(1);
xlabel('P_F')
ylabel('P_D')
legend('Es/\sigma^2=1','Es/\sigma^2=2','Es/\sigma^2=4','Es/\sigma^2=16','Theoretical SKE')
title('NARW propagated mismatched signal #1 amd signal #4')
roc1_try(2);roc1_try(4);roc1_try(16);
hold of