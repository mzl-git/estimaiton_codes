function pa=proposed_complex(N,x1,L)%长度 点数，阶数 多频
n=0:N-1;%复信号的
%%
% clc;clear
% N=256;A1=1.1;A2=1.5;A3=0.3;A4=0.2;theta1=0.1;theta2=0.2;theta3=0.2;theta4=0.25;
% f1=0.07;f2=0.13;f3=0.26;f4=0.34;%w2=pi*0.223;
% w1=2*pi*f1;w2=2*pi*f2;w3=2*pi*f3;w4=2*pi*f4;%w2=pi*0.223;
% L=4;%谐波阶数
% n=0:N-1; %设置信号
% d1=0;d2=0;d3=0;d4=0;
% x1=A1*cos(w1*n+theta1)+A2*cos(w2*n+theta2)+A3*cos(w3*n+theta3)+A4*cos(w4*n+theta4);%;%;%
% Q=2;
n=0:N-1;
s=x1;
for k=1:L
    %%
    x=[s,zeros(1,N)];
    Y=fft(x);
    [~, m2]=max(abs(Y));
    delta=0;SC=x; 
    for Q1=1:2
        f0=(m2-1+delta)/2/N;     %粗估计w0
%         A=sum(SC(1:N).*exp(-1j*2*pi*f0*n))/N;figure;plot(abs(Y))
%         r=conj(A).*exp(-1j*f0*2*pi*n);
%         r=[r zeros(1,N)];    
%         SC=x-r;
        SC=x;
        x0=sum(SC(1:N).*exp(-1j*n*(m2-1+delta)*pi/N));
        xp=sum(SC(1:N).*exp(-1j*n*(m2-1+delta+0.5)*pi/N));
        xp2=sum(SC(1:N).*exp(-1j*n*(m2-1+delta-0.5)*pi/N));
        AAA=sqrt(2)*abs(x0)*(abs(xp)-abs(xp2));
        BBB=4*abs(xp)*abs(xp2)*cos(pi/4/N)-sqrt(2)*abs(x0)*(abs(xp)+abs(xp2));
        h=(2/pi)*atan(AAA/BBB);
        delta=delta+h;
    end
    fatf(1,k)=(m2-1+delta)/2/N;    %
    AA=sum(SC(1:N).*exp(-1j*fatf(1,k)*n*2*pi))/N;   %用f0效果会好
    %AA=2*A;
    fatf(2,k)=abs(AA);fatf(3,k)=angle(AA);
    fatf(4,k)=AA;
    %s=s-fatf(2,k)*cos(2*pi*fatf(1,k)*n+fatf(3,k));
    s=s-fatf(4,k)*exp(1j*2*pi*fatf(1,k)*n);
end
%精估计
for qq=1:2
    for k=1:L
        s=x1;
            for q=1:L
                if q==k
                else
                    %s=s-fatf(2,q)*cos(2*pi*fatf(1,q)*n+fatf(3,q));
                    s=s-fatf(4,q)*exp(1j*2*pi*fatf(1,q)*n);
                end
            end
            %r=conj(fatf(4,k)/2).*exp(-1j*fatf(1,k)*2*pi*n);
            %使用粗估计的复幅值fatf(4,k)，而不是A=sum(SC(1:N).*exp(-1j*2*pi*fatf(1,k)*n))/N;重新生成复幅值，
            %能提高精度，减少迭代次数
            %SC=s-r;
            SC=s;
            x0=sum(SC(1:N).*exp(-1j*n*2*pi*fatf(1,k)));
            xp=sum(SC(1:N).*exp(-1j*n*(fatf(1,k)+0.5/2/N)*2*pi));
            xp2=sum(SC(1:N).*exp(-1j*n*(fatf(1,k)-0.5/2/N)*2*pi));
            AAA=sqrt(2)*abs(x0)*(abs(xp)-abs(xp2));
            BBB=4*abs(xp)*abs(xp2)*cos(pi/4/N)-sqrt(2)*abs(x0)*(abs(xp)+abs(xp2));
            h=(2/pi)*atan(AAA/BBB);
            fatf(1,k)=fatf(1,k)+h/2/N;  %
        AA=sum(SC(1:N).*exp(-1j*fatf(1,k)*n*2*pi))/N;
        %AA=2*A;2*
        fatf(2,k)=abs(AA);
        fatf(3,k)=angle(AA);fatf(4,k)=AA;
    end
end

pa=sortrows(fatf',1)';