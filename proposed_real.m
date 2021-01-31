function pa=proposed_real(N,x1,L)%
%%
n=0:N-1;
s=x1;
for k=1:L% coarse estimation
    %%
    x=[s,zeros(1,N)];
    Y=fft(x);
    [~, m2]=max(abs(Y));
    delta=0;SC=x; 
    for Q1=1:2
        f0=(m2-1+delta)/2/N;     %
        A=sum(SC(1:N).*exp(-1j*2*pi*f0*n))/N;
        r=conj(A).*exp(-1j*f0*2*pi*n);
        r=[r zeros(1,N)];    
        SC=x-r;
        x0=sum(SC(1:N).*exp(-1j*n*(m2-1+delta)*pi/N));
        xp=sum(SC(1:N).*exp(-1j*n*(m2-1+delta+0.5)*pi/N));
        xp2=sum(SC(1:N).*exp(-1j*n*(m2-1+delta-0.5)*pi/N));
        AAA=sqrt(2)*abs(x0)*(abs(xp)-abs(xp2));
        BBB=4*abs(xp)*abs(xp2)*cos(pi/4/N)-sqrt(2)*abs(x0)*(abs(xp)+abs(xp2));
        h=(2/pi)*atan(AAA/BBB);
        delta=delta+h;
    end
    fatf(1,k)=(m2-1+delta)/2/N;    %
    AA=2*sum(SC(1:N).*exp(-1j*fatf(1,k)*n*2*pi))/N;   %用f0效果会好
    %AA=2*A;
    fatf(2,k)=abs(AA);fatf(3,k)=angle(AA);
    fatf(4,k)=AA;
    s=s-fatf(2,k)*cos(2*pi*fatf(1,k)*n+fatf(3,k));
end
%fine estimation
for qq=1:3
    for k=1:L
        s=x1;
            for q=1:L
                if q==k
                else
                    s=s-fatf(2,q)*cos(2*pi*fatf(1,q)*n+fatf(3,q));
                end
            end
            r=conj(fatf(4,k)/2).*exp(-1j*fatf(1,k)*2*pi*n);
            %使用粗估计的复幅值fatf(4,k)，而不是A=sum(SC(1:N).*exp(-1j*2*pi*fatf(1,k)*n))/N;重新生成复幅值，
            %能提高精度，减少迭代次数
            SC=s-r;
            x0=sum(SC(1:N).*exp(-1j*n*2*pi*fatf(1,k)));
            xp=sum(SC(1:N).*exp(-1j*n*(fatf(1,k)+0.5/2/N)*2*pi));
            xp2=sum(SC(1:N).*exp(-1j*n*(fatf(1,k)-0.5/2/N)*2*pi));
            AAA=sqrt(2)*abs(x0)*(abs(xp)-abs(xp2));
            BBB=4*abs(xp)*abs(xp2)*cos(pi/4/N)-sqrt(2)*abs(x0)*(abs(xp)+abs(xp2));
            h=(2/pi)*atan(AAA/BBB);
            fatf(1,k)=fatf(1,k)+h/2/N;  %
        AA=2*sum(SC(1:N).*exp(-1j*fatf(1,k)*n*2*pi))/N;
        %AA=2*A;
        fatf(2,k)=abs(AA);
        fatf(3,k)=angle(AA);fatf(4,k)=AA;
    end
end

pa=sortrows(fatf',1)';