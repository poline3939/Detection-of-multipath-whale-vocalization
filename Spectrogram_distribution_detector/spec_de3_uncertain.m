clear
% close all
% clc
load hudson_env_prop.mat;
f=256;
w=256;
h=128;
fs=2000;
dd=size(su,1);
Es=[1 2 4 16];

for ii=1:4
    es=Es(ii);
    ss=zeros(dd,size(su,2));
    for jj=1:dd
        ss(jj,:)=sqrt(es)*su(jj,:)/sqrt(su(jj,:)*su(jj,:)');
        [mu,t1,f1]=stft(ss(jj,:),f,w,h,fs);
%         mnu=mu/max(mu(:));
        mu1=abs(mu).^2;
        muk(jj,:)=mu1(:);
    end

    [a1,b1]=size(mu1);
    L2=length(muk);
    yy=zeros(length(muk),1);
    zz=yy; zz2=yy; yy2=yy;

for k=1:1000
    rn=randperm(dd);
    cr=rn(1);
    x0=randn(1,size(su,2));
    x1=ss(cr,:)+x0;
    [K0,T0,F0]=stft(x0,f,w,h,fs);
%     Kn0=K0/max(K0(:));
    S0=abs(K0).^2;
    S0k=S0(:);
    [K1,T1,F1]=stft(x1,f,w,h,fs);
%     Kn1=K1/max(K1(:));
    S1=abs(K1).^2;
    S1k=S1(:);
    sk1=0;
    sk0=0;
    
    for kk=1:dd   
        zz0=2*sqrt(muk(dd,:).*S0k')/f;
        zz1=2*sqrt(muk(dd,:).*S1k')/f;
        lambda0=besseli(0,zz0);
        lambda1=besseli(0,zz1);
        ln0=log(lambda0)-muk(dd,:)/f;
        ln1=log(lambda1)-muk(dd,:)/f;
        sk1=sk1+exp(sum(ln1))/dd;
        sk0=sk0+exp(sum(ln0))/dd;
    end
    
    Ln0(k)=sk0;
    Ln1(k)=sk1;   
    percent = k/1000*100;
    clc;display(sprintf('Completed: %.1f%%',percent));
end
    inter=(max(max(Ln0),max(Ln1))-min(min(Ln0),min(Ln1)))/1000;
    inter2=inter/10;
    beta=min(min(Ln0),min(Ln1)):inter:max(max(Ln0),max(Ln1));
    count=histc(Ln1,beta);
    n_p1=count/sum(count);
    c_e1=cumsum(n_p1);

    count2=histc(Ln0,beta);
    n_p2=count2/sum(count2);
    c_e2=cumsum(n_p2);

    switch ii
        case 1
            KK1_2=1-c_e2;
            KK1_1=1-c_e1;
        case 2
            KK2_2=1-c_e2;
            KK2_1=1-c_e1;
        case 3
            KK3_2=1-c_e2;
            KK3_1=1-c_e1;
        case 4
            KK4_2=1-c_e2;
            KK4_1=1-c_e1;
    end
end  
    % ES=1, vc_max=3; ES=sqrt(2), vc_max=4; Es=2, vc_max=9;  Es=4, vc_max=23;
%     vc=-abs(min(beta)):inter:(max(beta)-inter);
%     t=1;
%     for j=1:length(vc)
%         while(beta(t)<vc(j))
%             t=t+1;
%         end
%         ww(j)=1-c_e1(t); %pd
%         vv(j)=1-c_e2(t); %pf
%         t=1;
%     end
%     switch ii
%         case 1
%             KK1_1=vv;
%             KK1_2=ww;
%         case 2
%             KK2_1=vv;
%             KK2_2=ww;
%         case 3
%             KK3_1=vv;
%             KK3_2=ww;
%         case 4
%             KK4_1=vv;
%             KK4_2=ww;
%     end
% end

% figure
% plot(KK1_1,KK1_2,'b')
% hold on
% plot(KK2_1,KK2_2,'k:',KK3_1,KK3_2,'r-.',KK4_1,KK4_2,'m--')
% roc1_try(1);
% 
% xlabel('P_F')
% ylabel('P_D')
% axis([0 1 0 1])
% legend('Es/\sigma^2=1','Es/\sigma^2=2','Es/\sigma^2=4','Es/\sigma^2=16','Time domain SKE')
% title('Spectrogram detection of NARW propagated signal under approximation')
% roc1_try(2);roc1_try(4);roc1_try(16);
% hold off

N=100;
ii=round(linspace(1,numel(KK1_1),N));
figure;
plot(KK1_2(ii),KK1_1(ii),'m^-')
hold on
plot(KK2_2(ii),KK2_1(ii),'b*-',KK3_2(ii),KK3_1(ii),'g+-',KK4_2(ii),KK4_1(ii),'ro-')
roc1_try(1);
xlabel('P_F','FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
ylabel('P_D','FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
axis([0 1 0 1])
legend('Es/\sigma^2=1','Es/\sigma^2=2','Es/\sigma^2=4','Es/\sigma^2=16','Time domain SKE')
% title('Spectrogram detection of NARW propagated signal under approximation')
roc1_try(2);roc1_try(4);roc1_try(16);
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
hold off

spe_uncertain16_pf=KK4_2(ii);
spe_uncertain16_pd=KK4_1(ii);
spe_uncertain4_pf=KK3_2(ii);
spe_uncertain4_pd=KK3_1(ii);
spe_uncertain2_pf=KK2_2(ii);
spe_uncertain2_pd=KK2_1(ii);
spe_uncertain1_pf=KK1_2(ii);
spe_uncertain1_pd=KK1_1(ii);
save('roc_spe_uncertain.mat','spe_uncertain16_pf','spe_uncertain16_pd','spe_uncertain4_pf','spe_uncertain4_pd','spe_uncertain2_pf','spe_uncertain2_pd','spe_uncertain1_pf','spe_uncertain1_pd')
            
            
            
            
