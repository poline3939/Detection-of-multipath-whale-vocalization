%% uncertain environmental parameters with spectrogram correaltion
% kernel_construct;
load NARWdata_prop.mat;
sss=[su1;su2;su3;su4;su5];
f=256;
w=256;
h=128;
fs=2000;

Es=[1 2 4 16];
t=1:5;
% load kernel.mat

for jj=1:4
%     kk=4 
s1=sqrt(Es(jj))*su1/sqrt(su1*su1');
L=length(s1);

[mu,t1,f1]=stft(s1,f,w,h,fs);
mnu=mu/max(mu(:));
mu1=abs(mnu).^2;
[a1,b1]=size(mu1);
% figure
% imagesc(t1,f1,mu1)
% figure
% imagesc(mu1)
% axis([9, b1, 9,b1])

f0=17;
f1=22;
d=27-22;
sigma=2; % instant freq
f2=1:a1;
tt1=1:5;
for i1=1:length(f2)
    for i2=1:length(tt1)
        x(i1,i2)=f2(i1)-(f0+(f1-f0)*tt1(i2)/d);
        ke1(i1,i2)=(1-x(i1,i2)^2/sigma^2).*exp(-x(i1,i2)^2/2/sigma^2);
    end
end
% figure
% imagesc(ke1)
% colormap(gray)
% title('NARW call kernel')
% xlabel('Time step')
% ylabel('Frequency step')
% set(gca, 'YDir', 'normal');

for i1=1:b1-length(tt1)+1
    alpha(i1)=sum(sum(mu1(f2,i1:(i1+length(tt1)-1)).*ke1));
end
% figure
% plot(alpha)
% title('Detection score for su1')
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
% title('Detection score for su1')
% xlabel('Time step')
% ylabel('Score')

%% su2
s2=sqrt(Es(jj))*su2/sqrt(su2*su2');
[mu2,t_2,f_2]=stft(s2,f,w,h,fs);
mnu2=mu2/max(mu2(:));
mu_2=abs(mnu2).^2;
% figure
% imagesc(t_2,f_2,mu_2)
% figure
% imagesc(mu_2)
% axis([9, b1, 9,b1])

f0=17;
f1=21;
d=27-22;
sigma=1.5; % instant freq
f2=1:a1;
tt2=1:5;
for i1=1:length(f2)
    for i2=1:length(tt2)
        x2(i1,i2)=f2(i1)-(f0+(f1-f0)*tt2(i2)/d);
        ke2(i1,i2)=(1-x2(i1,i2)^2/sigma^2).*exp(-x2(i1,i2)^2/2/sigma^2);
    end
end
% figure
% imagesc(ke2)
% colormap(gray)
% title('NARW call kernel for su2')
% xlabel('Time step')
% ylabel('Frequency step')
% set(gca, 'YDir', 'normal');

for i1=1:b1-length(tt2)+1
    alpha_s2(i1)=sum(sum(mu_2(f2,i1:(i1+length(tt2)-1)).*ke2));
end
% figure
% plot(alpha_s2)
% title('Detection score for s2')
% xlabel('Time step')
% ylabel('Score')
alpha2_s2=alpha_s2;

for i1=1:length(alpha2_s2)
    if alpha2_s2(i1)<0
        alpha2_s2(i1)=0;
    end
end
% figure
% plot(alpha2_s2)
% title('Detection score for s2')
% xlabel('Time step')
% ylabel('Score')

%% su3
s3=sqrt(Es(jj))*su3/sqrt(su3*su3');

[mu3,t_3,f_3]=stft(s3,f,w,h,fs);
mnu3=mu3/max(mu3(:));
mu_3=abs(mnu3).^2;
% figure
% imagesc(t_3,f_3,mu_3)
% figure
% imagesc(mu_3)
% axis([9, b1, 9,b1])

f0=17;
f1=21;
d=5;
sigma=1.5; % instant freq
f2=1:a1;
tt3=1:5;
for i1=1:length(f2)
    for i2=1:length(tt3)
        x3(i1,i2)=f2(i1)-(f0+(f1-f0)*tt3(i2)/d);
        ke3(i1,i2)=(1-x3(i1,i2)^2/sigma^2).*exp(-x3(i1,i2)^2/2/sigma^2);
    end
end
% figure
% imagesc(ke3)
% colormap(gray)
% title('NARW call kernel')
% xlabel('Time step')
% ylabel('Frequency step')
% set(gca, 'YDir', 'normal');

for i1=1:b1-length(tt3)+1
    alpha_s3(i1)=sum(sum(mu_3(f2,i1:(i1+length(tt3)-1)).*ke3));
end
% figure
% plot(alpha_s3)
% title('Detection score for s3')
% xlabel('Time step')
% ylabel('Score')
alpha2_s3=alpha_s3;

for i1=1:length(alpha2_s3)
    if alpha_s3(i1)<0
        alpha2_s3(i1)=0;
    end
end
% figure
% plot(alpha2_s3)
% title('Detection score for s3')
% xlabel('Time step')
% ylabel('Score')

%% su4 2nd segment
s4=sqrt(Es(jj))*su4/sqrt(su4*su4');

[mu4,t_4,f_4]=stft(s4,f,w,h,fs);
mnu4=mu4/max(mu4(:));
mu_4=abs(mnu4).^2;
% figure
% imagesc(t_4,f_4,mu_4)
% figure
% imagesc(mu_4)
% axis([9, b1, 9,b1])

f0=17;
f1=22;
d=27-22;
sigma=1.5; % instant freq
f2=1:a1;
tt4=1:5;
for i1=1:length(f2)
    for i2=1:length(tt4)
        x4(i1,i2)=f2(i1)-(f0+(f1-f0)*tt4(i2)/d);
        ke4_2(i1,i2)=(1-x4(i1,i2)^2/sigma^2).*exp(-x4(i1,i2)^2/2/sigma^2);
    end
end
% figure
% imagesc(ke4_2)
% colormap(gray)
% title('NARW call kernel for su4 2nd segment')
% xlabel('Time step')
% ylabel('Frequency step')
% set(gca, 'YDir', 'normal');

for i1=1:b1-length(tt4)+1
    alpha_s4(i1)=sum(sum(mu_4(f2,i1:(i1+length(tt4)-1)).*ke4_2));
end
% figure
% plot(alpha_s4)
% title('Detection score for su4')
% xlabel('Time step')
% ylabel('Score')
alpha2_s4=alpha_s4;

for i1=1:length(alpha2_s4)
    if alpha2_s4(i1)<0
        alpha2_s4(i1)=0;
    end
end
% figure
% plot(alpha2_s4)
% title('Detection score for su4 second segment')
% xlabel('Time step')

f0=14;
f1=14;
d=16-11;
sigma=0.5; % instant freq
f2=1:a1;
tt4_1=1:5;
for i1=1:length(f2)
    for i2=1:length(tt4_1)
        x4_1(i1,i2)=f2(i1)-(f0+(f1-f0)*tt4_1(i2)/d);
        ke4_1(i1,i2)=(1-x4_1(i1,i2)^2/sigma^2).*exp(-x4_1(i1,i2)^2/2/sigma^2);
    end
end
% figure
% imagesc(ke4_1)
% colormap(gray)
% title('NARW call kernel for 1st segment')
% xlabel('Time step')
% ylabel('Frequency step')
% set(gca, 'YDir', 'normal');

for i1=1:b1-length(tt4_1)+1
    alpha_s4_1(i1)=sum(sum(mu_4(f2,i1:(i1+length(tt4_1)-1)).*ke4_1));
end
% figure
% plot(alpha_s4_1)
% title('Detection score for su4 1st segment')
% xlabel('Time step')
% ylabel('Score')
alpha2_s4_1=alpha_s4_1;

for i1=1:length(alpha2_s4_1)
    if alpha2_s4_1(i1)<0
        alpha2_s4_1(i1)=0;
    end
end
% figure
% plot(alpha2_s4_1)
% title('Detection score for su4 1st segment')
% xlabel('Time step')
% ylabel('Score')

%% su5
s5=sqrt(Es(jj))*su5/sqrt(su5*su5');
[mu5,t_5,f_5]=stft(s5,f,w,h,fs);
mnu5=mu5/max(mu5(:));
mu_5=abs(mnu5).^2;
% figure
% imagesc(t_5,f_5,mu_5)
% figure
% imagesc(mu_5)
% axis([9, b1, 9,b1])

f0=17;
f1=21;
d=27-22;
sigma=1.5; % instant freq
f2=1:a1;
tt5=1:5;
for i1=1:length(f2)
    for i2=1:length(tt5)
        x5(i1,i2)=f2(i1)-(f0+(f1-f0)*tt5(i2)/d);
        ke5(i1,i2)=(1-x5(i1,i2)^2/sigma^2).*exp(-x5(i1,i2)^2/2/sigma^2);
    end
end
% figure
% imagesc(ke5)
% colormap(gray)
% title('NARW call kernel for su5')
% xlabel('Time step')
% ylabel('Frequency step')
% set(gca, 'YDir', 'normal');

for i1=1:b1-length(tt5)+1
    alpha_s5(i1)=sum(sum(mu_5(f2,i1:(i1+length(tt5)-1)).*ke5));
end
% figure
% plot(alpha_s5)
% title('Detection score for su5')
% xlabel('Time step')
% ylabel('Score')
alpha2_s5=alpha_s5;

for i1=1:length(alpha2_s5)
    if alpha2_s5(i1)<0
        alpha2_s5(i1)=0;
    end
end
% figure
% plot(alpha2_s5)
% title('Detection score for su5')
% xlabel('Time step')
% ylabel('Score')

ss=[s1;s2;s3;s4;s5];
t=1:5;

    co=zeros(1,10000); con=co; co1=co; con1=co;
    %%% spectrogram of signal+noise and pure noise
    for ii=1:10000
        y=randn(1,L);
        rn=randperm(5);
        cr=rn(1);        
        zz=y+ss(cr,:);
        
        [sn2,tn,fn]=stft(y,f,w,h,fs);
        [sz2,tz,fz]=stft(zz,f,w,h,fs);
        ssn2=sn2/max(sn2(:));
        sn=abs(ssn2).^2;
        ssz2=sz2/max(sz2(:));
        sz=abs(ssz2).^2;
        
        for i1=1:b1-length(t)+1
            alpha_n1(i1)=sum(sum(sn(f2,i1:(i1+length(t)-1)).*ke1));
            alpha_n2(i1)=sum(sum(sn(f2,i1:(i1+length(t)-1)).*ke2));
            alpha_n3(i1)=sum(sum(sn(f2,i1:(i1+length(t)-1)).*ke3));
            alpha_n41(i1)=sum(sum(sn(f2,i1:(i1+length(t)-1)).*ke4_1));
            alpha_n5(i1)=sum(sum(sn(f2,i1:(i1+length(t)-1)).*ke5));
            alpha_n42(i1)=sum(sum(sn(f2,i1:(i1+length(t)-1)).*ke4_2));
        end
        con(ii)=(max(alpha_n1)+max(alpha_n2)+max(alpha_n3)+(max(alpha_n41)+max(alpha_n42))/2+max(alpha_n5))/5;
        
        for i1=1:b1-length(t)+1
            alpha_z1(i1)=sum(sum(sz(f2,i1:(i1+length(t)-1)).*ke1));
            alpha_z2(i1)=sum(sum(sz(f2,i1:(i1+length(t)-1)).*ke2));
            alpha_z3(i1)=sum(sum(sz(f2,i1:(i1+length(t)-1)).*ke3));
            alpha_z41(i1)=sum(sum(sz(f2,i1:(i1+length(t)-1)).*ke4_1));
            alpha_z42(i1)=sum(sum(sz(f2,i1:(i1+length(t)-1)).*ke4_2));
            alpha_z5(i1)=sum(sum(sz(f2,i1:(i1+length(t)-1)).*ke5));
        end
        co(ii)=(max(alpha_z1)+max(alpha_z2)+max(alpha_z3)+(max(alpha_z41)+max(alpha_z42))/2+max(alpha_z5))/5;
        
        percent = ii/10000*100;
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
    switch jj
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
% title('LRT of NARW propagated signal using spectrogram correlation')
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
% title('LRT of NARW propagated signal using spectrogram correlation')
roc1_try(2);roc1_try(4);roc1_try(16);
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
hold off

speccor_uncer16_pf=KK4_2(ii);
speccor_uncer16_pd=KK4_1(ii);
speccor_uncer4_pf=KK3_2(ii);
speccor_uncer4_pd=KK3_1(ii);
speccor_uncer2_pf=KK2_2(ii);
speccor_uncer2_pd=KK2_1(ii);
speccor_uncer1_pf=KK1_2(ii);
speccor_uncer1_pd=KK1_1(ii);
save('speccor_uncer.mat','speccor_uncer16_pf','speccor_uncer16_pd','speccor_uncer4_pf','speccor_uncer4_pd','speccor_uncer2_pf','speccor_uncer2_pd','speccor_uncer1_pf','speccor_uncer1_pd')
            
            
            
            
