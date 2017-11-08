%% match of signal with spectrogram correaltion
clear
f=256;
w=256;
h=128;
fs=2000;

load updata
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
ln1=zeros(1,L/2);
ln0=zeros(1,L/2);   

[mu,T1,F1,~]=spectrogram(s_mea,256,128,256,2000,'yaxis');
mnu=mu/max(mu(:));
mu1=abs(mnu).^2;
[a1,b1]=size(mu1);
figure
imagesc(T1,F1,mu1)

%%% spectrogram of signal+noise and pure noise
for ii=1:L/2
    y=randn(1,n);
    zz=y+Up_1(L/2+ii,:);
    [sn2,Tn,Fn,~]=spectrogram(y,256,128,256,2000,'yaxis');
    [sz2,Tz,Fz,~]=spectrogram(zz,256,128,256,2000,'yaxis');
    ssn2=sn2/max(sn2(:));
    sn=abs(ssn2).^2;
    ssz2=sz2/max(sz2(:));
    sz=abs(ssz2).^2;

    con(ii)=sum(sum(sn.*mu1));
    co(ii)=sum(sum(sz.*mu1));    

    percent = ii/(L/2)*100;
    clc;display(sprintf('Completed: %.1f%%',percent));
end

inter=(max(max(con),max(co))-min(min(co),min(con)))/1000;
beta=min(min(co),min(con)):inter:max(max(co),max(con));
count=histc(co,beta);
n_p1=count/sum(count);
c_e1=cumsum(n_p1);
count2=histc(con,beta);
n_p2=count2/sum(count2);
c_e2=cumsum(n_p2);

vc=-abs(min(beta)):inter:(max(beta)-inter);
t=1;
for j=1:length(vc)
    while(beta(t)<vc(j))
          t=t+1;
    end
    ww(j)=1-c_e1(t); %pd
    vv(j)=1-c_e2(t); %pf
    t=1;
end
KK3_1=vv; % real data
KK3_2=ww; % real data
auc_n3=RocIntegral(KK3_1,KK3_2,length(KK3_1)-1); % 0.6734
save('roc_real_data.mat','KK3_1','KK3_2','auc_n3')
