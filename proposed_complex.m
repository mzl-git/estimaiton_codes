function pa=proposed_complex(N,x1,L,Qc,Qf)%length, signal, number_of_sinusoids, iteration_number_coarse_estimation,fine_estimation
%This is the code for multiple complex sinusoids
n=0:N-1;
SC=x1;fatf=zeros(2,L);
for k=1:L
    %%
    Y=fft(SC,2*N);
    [~, m2]=max(abs(Y));
    h=0; 
    for Q1=1:Qc
        x0=Y(m2);
        xp=sum(SC(1:N).*exp(-1j*n*(m2-1+h+0.5)*pi/N));
        xp2=sum(SC(1:N).*exp(-1j*n*(m2-1+h-0.5)*pi/N));
        AAA=sqrt(2)*abs(x0)*(abs(xp)-abs(xp2));
        BBB=4*abs(xp)*abs(xp2)*cos(pi/4/N)-sqrt(2)*abs(x0)*(abs(xp)+abs(xp2));
        h=(2/pi)*atan(AAA/BBB);
    end
    fatf(1,k)=(m2-1+h)/2/N;   
    AA=sum(SC(1:N).*exp(-1j*fatf(1,k)*n*2*pi))/N;
    fatf(2,k)=AA;
    SC=SC-fatf(2,k)*exp(1j*2*pi*fatf(1,k)*n);
end
% fine estimation
for qq=1:Qf
    for k=1:L
        s=x1;
        for q=1:L
            if q==k
            else
                s=s-fatf(2,q)*exp(1j*2*pi*fatf(1,q)*n);
            end
        end
        x0=sum(s(1:N).*exp(-1j*n*2*pi*fatf(1,k)));
        xp=sum(s(1:N).*exp(-1j*n*(fatf(1,k)+0.25/N)*2*pi));
        xp2=sum(s(1:N).*exp(-1j*n*(fatf(1,k)-0.25/N)*2*pi));
        AAA=sqrt(2)*abs(x0)*(abs(xp)-abs(xp2));
        BBB=4*abs(xp)*abs(xp2)*cos(pi/4/N)-sqrt(2)*abs(x0)*(abs(xp)+abs(xp2));
        h=(2/pi)*atan(AAA/BBB);
        fatf(1,k)=fatf(1,k)+h/2/N;  %
        AA=sum(s(1:N).*exp(-1j*fatf(1,k)*n*2*pi))/N;
        fatf(2,k)=AA;
    end
end

pa=sortrows(fatf',1)';