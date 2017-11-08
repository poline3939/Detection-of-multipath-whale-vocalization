
%fft with w(x)\equiv 1
function d1= stft3_new(x, M, D, B, J)
%x= input signal 
%M= siaze of the window
%D= move in each step 

N = length(x);

d1 = zeros(B,J);

for a=1:J
    for b=1:B
        temp=0;
        for t=0:M-1
            temp=temp+x((a-1)*D+t+1)*exp(-i*2*pi*(b-1)*t/M);
        end
        d1(b,a)=temp;
    end
end
