% test stft
% clear
% close all
% clc
load NARWdata_prop.mat;

N=16; %N=length of test signal 
M=4; %M= size of the window
D=4; %D= move in each step 
B=1+fix(M/2); % upper bound of b
J=1+fix((N-M)/D); % upper bound of a
L_U=B*J;
L_X=2*B*J;

su1=Pp_35m;

test_signal=su1(3001:3001+N);

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

X=test_signal*HC_x';
mu=X*100;

%---------------------------------------------------------------------

[Q,Lambda]=eig(C_x); % Q\LambdaQ^T=C_x
[E_v,order]=sort(diag(Lambda),'descend');
Q=Q(:,order);
r=find(E_v>10^(-5));
Q_1=Q(:,1:r(end)); % Q=[Q_1^T,Q_2^T]^T Q_1 this the linearly idpendent part
Lambda_1=diag(E_v(1:r(end))); 
Pr=Q_1*Lambda_1^(-0.5);

d_vecotr=Pr'*mu';
d_square=d_vecotr'*d_vecotr;

%---- Test the emperical and analytic pdf 
F_vector=zeros(1,0);
T=200000;
for j=1:T
    x_test=randn(1,N);
    X_2=stft3(x_test,M, D,B,J); %B*J matrix of FFT's
    U_2=real(X_2);
    V_2=imag(X_2);
    U_vector2=U_2(:);
    V_vector2=V_2(:);
    x_test=[U_vector2' V_vector2'];
    temp_data=x_test*Pr*d_vecotr-d_square/2;
    F_vector=[F_vector temp_data];
    j
end
x=-3:0.005:3;
yy=hist(F_vector,x); 
yy=yy/length(F_vector)*200; %?????????
bar(x,yy) 
hold on
an_pdf=zeros(0,1);
for tt=1:length(x)
    an_pdf=[an_pdf normpdf(x(tt),-d_square/2,d_square^(0.5))];
end
plot(x,an_pdf,'r')


