clear
% close all
% clc
load NARWdata_prop.mat;
% load hudson_env_prop_mean_ocean.mat
s=su1;
f=256;
w=256;
h=128;
fs=2000;

Es=[1 2 4 16];

for ii=1:4
% ii=3;
s1=sqrt(Es(ii))*s/sqrt(s*s');
L=length(s1);

[mu,t1,f1]=stft(s1,f,w,h,fs);
mnu=mu/max(mu(:));
mu1=abs(mnu).^2;

[a1,b1]=size(mu1);
muk=mu1(:);
L2=length(muk);

yy=zeros(length(muk),1);
zz=yy; zz2=yy; yy2=yy;

for k=1:1000
    x0=randn(1,L);
    x1=s1+x0;
    
    [K0,T0,F0]=stft(x0,f,w,h,fs);
    Kn0=K0/max(K0(:));
    S0=abs(Kn0).^2;
    S0k=S0(:);
    [K1,T1,F1]=stft(x1,f,w,h,fs);
    Kn1=K1/max(K1(:));
    S1=abs(Kn1).^2;
    S1k=S1(:);
    
    zz0=2*sqrt(muk.*S0k)/f;
    zz1=2*sqrt(muk.*S1k)/f;
    lambda0=besseli(0,zz0);
    lambda1=besseli(0,zz1);
    ln0=log(lambda0)-muk/f;
    ln1=log(lambda1)-muk/f;    

   Ln0(k)=sum(ln0);
   Ln1(k)=sum(ln1);
   
   percent  = k/1000*100;
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

    % ES=1, vc_max=3; ES=sqrt(2), vc_max=4; Es=2, vc_max=9;  Es=4, vc_max=23;
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
    switch ii
        case 1
            KK1_1=vv;
            KK1_2=ww;
        case 2
            KK2_1=vv;
            KK2_2=ww;
        case 3
            KK3_1=vv;
            KK3_2=ww;
        case 4
            KK4_1=vv;
            KK4_2=ww;
    end
end

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
plot(KK1_1(ii),KK1_2(ii),'m^-')
hold on
plot(KK2_1(ii),KK2_2(ii),'b*-',KK3_1(ii),KK3_2(ii),'g+-',KK4_1(ii),KK4_2(ii),'ro-')
roc1_try(1);
xlabel('P_F','FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
ylabel('P_D','FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
axis([0 1 0 1])
legend('Es/\sigma^2=1','Es/\sigma^2=2','Es/\sigma^2=4','Es/\sigma^2=16','Time domain SKE')
% title('Spectrogram detection of NARW propagated signal under approximation')
roc1_try(2);roc1_try(4);roc1_try(16);
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
hold off

spe_ske16_pd=KK4_2(ii);
spe_ske16_pf=KK4_1(ii);
spe_ske4_pd=KK3_2(ii);
spe_ske4_pf=KK3_1(ii);
spe_ske2_pd=KK2_2(ii);
spe_ske2_pf=KK2_1(ii);
spe_ske1_pd=KK1_2(ii);
spe_ske1_pf=KK1_1(ii);
save('roc_spe_ske.mat','spe_ske16_pf','spe_ske16_pd','spe_ske4_pf','spe_ske4_pd','spe_ske2_pf','spe_ske2_pd','spe_ske1_pf','spe_ske1_pd')
             
%             
            
            
