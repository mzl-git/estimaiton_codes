function pa=proposed_real(N,x1,L,Qc,Qf)%length, signal, number_of_sinusoids, iteration_number_coarse_estimation,fine_estimation
%This is the code for multiple real sinusoids
%%
n=0:N-1;
s=x1;
for k=1:L
    %%
    x=[s,zeros(1,N)];
    Y=fft(x);
    [~, m2]=max(abs(Y));
    h=0;SC=x; 
    for Q1=1:Qc
        f0=(m2-1)/2/N;    
        A=sum(SC(1:N).*exp(-1j*2*pi*f0*n))/N;
        r=conj(A).*exp(-1j*f0*2*pi*n);
        r=[r zeros(1,N)];    
        SC=x-r;
        x0=sum(SC(1:N).*exp(-1j*n*(m2-1)*pi/N));
        xp=sum(SC(1:N).*exp(-1j*n*(m2-1+0.5)*pi/N));
        xp2=sum(SC(1:N).*exp(-1j*n*(m2-1-0.5)*pi/N));
        AAA=sqrt(2)*abs(x0)*(abs(xp)-abs(xp2));
        BBB=4*abs(xp)*abs(xp2)*cos(pi/4/N)-sqrt(2)*abs(x0)*(abs(xp)+abs(xp2));
        h=(2/pi)*atan(AAA/BBB);
    end
    fatf(1,k)=(m2-1+h)/2/N;    
    AA=2*sum(SC(1:N).*exp(-1j*fatf(1,k)*n*2*pi))/N;   
    fatf(2,k)=abs(AA);fatf(3,k)=angle(AA);
    fatf(4,k)=AA;
    s=s-fatf(2,k)*cos(2*pi*fatf(1,k)*n+fatf(3,k));
end
%fine_estimation
for qq=1:Qf
    for k=1:L
        s=x1;
            for q=1:L
                if q==k
                else
                    s=s-fatf(2,q)*cos(2*pi*fatf(1,q)*n+fatf(3,q));
                end
            end
            r=conj(fatf(4,k)/2).*exp(-1j*fatf(1,k)*2*pi*n);
            SC=s-r;
            x0=sum(SC(1:N).*exp(-1j*n*2*pi*fatf(1,k)));
            xp=sum(SC(1:N).*exp(-1j*n*(fatf(1,k)+0.5/2/N)*2*pi));
            xp2=sum(SC(1:N).*exp(-1j*n*(fatf(1,k)-0.5/2/N)*2*pi));
            AAA=sqrt(2)*abs(x0)*(abs(xp)-abs(xp2));
            BBB=4*abs(xp)*abs(xp2)*cos(pi/4/N)-sqrt(2)*abs(x0)*(abs(xp)+abs(xp2));
            h=(2/pi)*atan(AAA/BBB);
            fatf(1,k)=fatf(1,k)+h/2/N;  %
        AA=2*sum(SC(1:N).*exp(-1j*fatf(1,k)*n*2*pi))/N;
        fatf(2,k)=abs(AA);
        fatf(3,k)=angle(AA);fatf(4,k)=AA;
    end
end

pa=sortrows(fatf',1)';