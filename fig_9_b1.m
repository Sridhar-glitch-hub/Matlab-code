clc;
clear;
gamma=0.05;
gamma_p=0.005;
phi=pi;
epsilon=2*gamma;
delta=0.5;
t=gamma;
w=linspace(-2,2,100000);
G_L=zeros(6); G_R=zeros(6); G_P=zeros(6);
G_L(1:3,1:3)=gamma; G_P(4,4)=gamma_p; G_R(4:6,4:6)=gamma;


for j=1:1:length(w)
    

% Retarded green function for 6QD(3,3)
                A = zeros(6);

                A(1,1)=w(j)-epsilon+delta+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(1,3)=1i*gamma/2;
                A(1,6)=-t*exp(-1i*phi/6);
                A(2,1)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(2,2)=w(j)-epsilon+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(3,1)=1i*gamma/2;
                A(3,2)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(3,3)=w(j)-epsilon-delta+1i*gamma/2;
                A(3,4)=-t*exp(1i*phi/6);
                A(4,3)=-t*exp(-1i*phi/6);
                A(4,4)=w(j)-epsilon+1i*gamma/2+1i*gamma_p/2;
                A(4,5)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(4,6)=1i*gamma/2;
                A(5,4)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(5,5)=w(j)-epsilon+1i*gamma/2;
                A(5,6)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(6,1)=-t*exp(1i*phi/6);
                A(6,4)=1i*gamma/2;
                A(6,5)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(6,6)=w(j)-epsilon+1i*gamma/2;
Gr=inv(A);
Ga=Gr';
% Transmission function for T_LR
T_LR(j)=trace(G_L*Gr*G_R*Ga);
T_LP(j)=trace(G_L*Gr*G_P*Ga);
T_PL(j)=trace(G_P*Gr*G_L*Ga);
T_PR(j)=trace(G_P*Gr*G_R*Ga);
end



x=w-epsilon;
%hold on
figure; 

plot(x,T_LR,'LineWidth',2)
xlabel('\omega-\epsilon')
ylabel('T_{LR}')
figure;
plot(x,T_LP,'LineWidth',2)
xlabel('\omega-\epsilon')
ylabel('T_{LP}')
figure;
plot(x,T_PL,'linewidth',2)
xlabel('\omega-\epsilon')
ylabel('T_{PL}')
figure; 
plot(x,T_LR+T_LP,'LineWidth',2)
xlabel('\omega-\epsilon')
ylabel('T_{LR}+T_{LP}')
figure; 
plot(x,T_PL+T_PR,'LineWidth',2)
xlabel('\omega-\epsilon')
ylabel('T_{PL}+T_{PR}')

