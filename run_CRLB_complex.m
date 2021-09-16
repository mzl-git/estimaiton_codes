clc
clear;
%close all
N=128;A1=1.5;A2=1.1;A3=0.3;A4=0.2;A5=0.5;theta1=0.2*pi;theta2=0.1*pi;theta3=0;theta4=0.25*pi;theta5=0.05*pi;
f1=0.05;f2=0.13;f3=0.26;f4=0.34;f5=0.44;
L=5;
n=0:N-1; 
idx=1;SNR1=-10;SNR2=70;step=1;

%% ·ÂÕæ
for snr=SNR1:step:SNR2
	SNR=10^(0.1*snr);% ¦Ñ
    CRB(idx)=6/(4*pi^2*(N*(N^2-1))*10^(0.1*snr));%CRLB for complex sinusoids
    theta1=(rand(1)-0.5)*2*pi;theta2=(rand(1)-0.5)*2*pi;theta3=(rand(1)-0.5)*2*pi;theta4=(rand(1)-0.5)*2*pi;theta5=(rand(1)-0.5)*2*pi;
    s=A1*exp(1i*((f1)*2*pi*n+theta1))+A2*exp(1i*((f2)*2*pi*n+theta2))+A3*exp(1i*((f3)*2*pi*n+theta3))...
        +A4*exp(1j*((f4)*2*pi*n+theta4))+A5*exp(1i*((f5)*2*pi*n+theta5));%;        
    for num_of_runs=1:100
        z=randn(1,N)+1j*randn(1,N);%noise for complex signal
        SIGMA=0.5*(A1^2)/(10^(snr/10));
        Z=sqrt(SIGMA)*z;
        x=s+Z;
        %%
        wfrealCTsubmf=proposed_complex(N,x,L,1,2);
        %%
        sum_mse_frealCTsubmf(num_of_runs)=(f1-wfrealCTsubmf(1,1)); 
    end
    %%
    mse_frealCTsubmf(idx)=(mean((sum_mse_frealCTsubmf).^2));
    idx=idx+1;snr
end
%% »æÍ¼
xs=SNR1:step:SNR2;
figure
set(gcf,'position',[200,300,800,400]);
plot(xs,10.*log10(CRB),'k-');
hold on
plot(xs,10.*log10(mse_frealCTsubmf),'^-','color',[235/256 1/256 14/256]);
legend('CRLB','Proposed_complex Qf=2','Location','southwest');%'2017SEA·¨','¸´ÐÅºÅCRLB','Nuttall windowing',
xlabel('SNR(dB)','FontName','Times New Roman','FontSize',10);
ylabel('MSE(dB)','FontName','Times New Roman','FontSize',10);
set(gca,'FontName','Times New Roman','FontSize',10);%,'LineWidth',1.5
