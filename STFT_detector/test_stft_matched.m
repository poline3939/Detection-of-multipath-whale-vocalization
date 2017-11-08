% test stft detector mismatch
clear
% close all
% clc
load hudson_env_prop.mat;

N=32; %N=length of test signal 
M=8; %M= size of the window
D=8; %D= move in each step 
B=1+fix(M/2); % upper bound of b
J=1+fix((N-M)/D); % upper bound of a
L_U=B*J;
L_X=2*B*J;

Es=[1 2 4 16];

HC_x=zeros(L_X,N+1);
tic;
for a=1:J
    for b=1:B
        for t=0:M-1
        HC_x((a-1)*B+b,(a-1)*D+t+1)=cos(-2*pi*(b-1)*t/M);
        HC_x(L_U+(a-1)*B+b,(a-1)*D+t+1)=sin(-2*pi*(b-1)*t/M);
        end
    end
end
C_x=HC_x*HC_x';

for i=1:L_X
    for j=1:L_X
        if abs(C_x(i,j))<10^(-5)
            C_x(i,j)=0;
        end
    end
end
toc;
[Q,Lambda]=eig(C_x); % Q\LambdaQ^T=C_x
[E_v,order]=sort(diag(Lambda),'descend');
Q=Q(:,order);
r=find(E_v>10^(-5));
Q_1=Q(:,1:r(end)); % Q=[Q_1^T,Q_2^T]^T Q_1 this the linearly idpendent part
Lambda_1=diag(E_v(1:r(end))); 
Pr=Q_1*Lambda_1^(-0.5);

% X=s1*HC_x';
% mu=X*100;
su1=su(37,:);
su4=su(37,:);
s1=su1(3001:3000+N);
s2=su4(3001:3000+N);

for jj=1:4
s1=sqrt(Es(jj))*s1/sqrt(s1*s1');
s2=sqrt(Es(jj))*s2/sqrt(s2*s2');

mu=stft3_new(s1, M, D, B, J);
mu=mu(:);
ne_mu=[real(mu);imag(mu)];

mu2=stft3_new(s2, M, D, B, J);
mu2=mu2(:);
ne_mu2=[real(mu2);imag(mu2)];

%---------------------------------------------------------------------
d_vecotr=Pr'*ne_mu;
d_square=d_vecotr'*d_vecotr;
dd_sq(jj)=d_square;

d_v2=Pr'*ne_mu2;
d_sq2=d_v2'*d_vecotr;
dd_sq2(jj)=d_sq2;

%---- Test the emperical and analytic pdf 
% F_vector=zeros(1,0);
T=2000;
L1=zeros(1,T);
L0=zeros(1,T);

for ii=1:T
    x_test=randn(1,N);
    x_test1=x_test+s2;
    X_0=stft3_new(x_test, M, D,B,J); %B*J matrix of FFT's
    X_1=stft3_new(x_test1, M, D,B,J); %B*J matrix of FFT's
    X_0=X_0(:);
    X_1=X_1(:);
    ne_X0=[real(X_0); imag(X_0)];
    ne_X1=[real(X_1); imag(X_1)];
    L0(ii)=ne_X0'*Pr*d_vecotr-d_square/2;
    L1(ii)=ne_X1'*Pr*d_vecotr-d_square/2;
    
    percent = ii/2000*100;
    clc;display(sprintf('Completed: %.1f%%',percent));   
end

inter=(max(max(L0),max(L1))-min(min(L0),min(L1)))/1000;
inter2=inter/10;
beta=min(min(L0),min(L1)):inter:max(max(L0),max(L1));
count=histc(L1,beta);
n_p1=count/sum(count);
c_e1=cumsum(n_p1); %pd

count2=histc(L0,beta);
n_p2=count2/sum(count2);
c_e2=cumsum(n_p2); %pf

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

N=100;
ii=round(linspace(1,numel(KK1_1),N));
figure;
plot(KK1_2(ii),KK1_1(ii),'m^-')
hold on
plot(KK2_2(ii),KK2_1(ii),'b*-',KK3_2(ii),KK3_1(ii),'g+-',KK4_2(ii),KK4_1(ii),'ro-')
roc1_try(dd_sq(1));
xlabel('P_F')
ylabel('P_D')
axis([0 1 0 1])
legend('Es/\sigma^2=1','Es/\sigma^2=2','Es/\sigma^2=4','Es/\sigma^2=16','Analytical STFT SKE')
% title('SKE based on STFT Distribution ')
roc1_try(dd_sq(2));roc1_try(dd_sq(3));roc1_try(dd_sq(4));
hold off

stft_snr16_pf=KK4_2(ii);
stft_snr16_pd=KK4_1(ii);
stft_snr4_pf=KK3_2(ii);
stft_snr4_pd=KK3_1(ii);
stft_snr2_pd=KK2_1(ii);
stft_snr2_pf=KK2_2(ii);
stft_snr1_pf=KK1_2(ii);
stft_snr1_pd=KK1_1(ii);
save('roc_stft_match.mat','stft_snr16_pd','stft_snr16_pf','stft_snr4_pd','stft_snr4_pf','stft_snr2_pd','stft_snr2_pf','stft_snr1_pd','stft_snr1_pf')
