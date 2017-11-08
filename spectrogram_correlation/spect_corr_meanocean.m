%% mismatch of signal with spectrogram correaltion
clear
% close all
% load NARWdata_prop.mat;
load hudson_env_prop.mat
load hudson_env_prop_mean_ocean.mat
s_mea=zeros(1,5500);
s_mea(1:length(s_m))=s_m;
S_mea=s_mea/sqrt(s_mea*s_mea');
dd=length(s_mea);
% s2=S_mea;

Su=zeros(size(su));
for ii=1:size(su,1)
    Su(ii,:)=su(ii,:)/sqrt(su(ii,:)*su(ii,:)');
end

% s=su1;
f=256;
w=256;
h=128;
fs=2000;

Es=[1 2 4 16];
for kk=1:4
    es=Es(kk);
    clear Ln1 Ln0 co con
    Ln1=[]; Ln0=[];
    
    s1=sqrt(es)*Su;
    sm2=sqrt(es)*S_mea;
    L=length(s1);

    [mu,t1,f1]=stft(sm2,f,w,h,fs);
    mnu=mu/max(mu(:));
    mu1=abs(mnu).^2;
    [a1,b1]=size(mu1);
%     figure
%     imagesc(t1,f1,mu1)
%     figure
%     imagesc(mu1)
%     axis([9, b1, 9,b1])

    f0=17;
    f1=21;
    d=27-22;
    sigma=2; % instant freq
    f2=1:a1;
    t=1:5;
    for i1=1:length(f2)
        for i2=1:length(t)
            x(i1,i2)=f2(i1)-(f0+(f1-f0)*t(i2)/d);
            ke(i1,i2)=(1-x(i1,i2)^2/sigma^2).*exp(-x(i1,i2)^2/2/sigma^2);
        end
    end
%     figure
%     imagesc(ke)
%     colormap(gray)
%     title('NARW call kernel')
%     xlabel('Time step')
%     ylabel('Frequency step')
%     set(gca, 'YDir', 'normal');
%     figure
%     imagesc(mu1)

    for i1=1:b1-length(t)+1
        alpha(i1)=sum(sum(mu1(f2,i1:(i1+length(t)-1)).*ke));
    end
% figure
% plot(alpha)
% title('Detection score for mismatch signal#1 and signal#4')
% xlabel('Time step')
% ylabel('Score')
alpha2=alpha;

for i1=1:length(alpha2)
    if alpha2(i1)<0
        alpha2(i1)=0;
    end
end
% figure
% plot(alpha2)
% title('Detection score for mismatch signal#1 and signal#4')
% xlabel('Time step')
% ylabel('Score')
    
    co=[]; con=co; co1=co; con1=co;
    clear Mu
    %%% spectrogram of signal+noise and pure noise
    for ii=1:1000
        y=randn(1,L);
        for kkk=1:size(Su,1)
            zz(kkk,:)=y+s1(kkk,:);
            [sz2,tz,fz]=stft(zz(kkk,:),f,w,h,fs);
            ssz2=sz2/max(sz2(:));
            sz=abs(ssz2).^2;
            Mu{kkk}=sz;
        end            
        [sn2,tn,fn]=stft(y,f,w,h,fs);
        ssn2=sn2/max(sn2(:));
        sn=abs(ssn2).^2;        
        
        for i1=1:b1-length(t)+1
            alpha_n(i1)=sum(sum(sn(f2,i1:(i1+length(t)-1)).*ke));
        end
        for i1=1:length(alpha_n)
            if alpha_n(i1)<0
                alpha_n(i1)=0;
            end
        end        
        con(ii)=max(alpha_n);
        clear ggg6 ggg5
        for kkk=1:size(Su,1)
            ggg5=Mu{kkk};
            for i1=1:b1-length(t)+1
                alpha_z(i1)=sum(sum(ggg5(f2,i1:(i1+length(t)-1)).*ke));
                if alpha_z(i1)<0
                    alpha_z(i1)=0;
                end
            end
            ggg6(kkk)=max(alpha_z);
        end           
        co=[co; ggg6(:)];        
        percent = ii/1000*100;
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
    switch kk
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

% figure
% plot(KK1_2,KK1_1,'r',KK2_2,KK2_1,'b',KK3_2,KK3_1,'m',KK4_2,KK4_1,'g')
% axis([0 1 0 1])
% hold on
% roc1_try(1);
% xlabel('P_F')
% ylabel('P_D')
% legend('Es/\sigma^2=1','Es/\sigma^2=2','Es/\sigma^2=4','Es/\sigma^2=16','SKE theoretical')
% title('Mismatched signal#1 with signal#4 using spectrogram correlation')
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
% title('Mismatched signal#1 with signal#4 using spectrogram correlation')
roc1_try(2);roc1_try(4);roc1_try(16);
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
hold off

%% su1 su4
speccor_mean16_pf=KK4_2(ii);
speccor_mean16_pd=KK4_1(ii);
speccor_mean4_pf=KK3_2(ii);
speccor_mean4_pd=KK3_1(ii);
speccor_mean2_pf=KK2_2(ii);
speccor_mean2_pd=KK2_1(ii);
speccor_mean1_pf=KK1_2(ii);
speccor_mean1_pd=KK1_1(ii);
save('roc_speccor_mean_ocean.mat','speccor_mean16_pf','speccor_mean16_pd','speccor_mean4_pf','speccor_mean4_pd','speccor_mean2_pf','speccor_mean2_pd','speccor_mean1_pf','speccor_mean1_pd')
            
            
            
            
