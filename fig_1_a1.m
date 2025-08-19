clc;
clear;
gamma=0.05;
E_0=0.3;
E=linspace(-5,5,100000);
T_s=8*gamma;
T_a=0.3 ; T_d=0.3;
mu=0.2; mu_d=0.2;

q=linspace(1,20,1000);

for i=1:1:length(q)
for j=1:1:length(E)
    epsilon(j)=(E(j)-E_0)/gamma;
    T(j,i)=(q(i)+epsilon(j)).^2/(1+epsilon(j).^2);
    T_n(j,i)=T(j,i)./(1+q(i).^2); % Normalized Fano transmission
    
   
    f(j)=1./(exp((E(j)-mu)/T_a)+1); % Fermi distribution
   
    n(j)=exp((E(j)-mu)/T_a);
    d(j)=T_a.*(exp((E(j)-mu)/T_a)+1).^2;
    Df(j)=-n(j)/d(j);   % derivative of Fermi distribution

    F(j,i)=-T_a.*T_n(j,i).*(E(j)-mu).*Df(j);
    F1(j,i)=-T_a.*T_n(j,i).*Df(j);
    F2(j,i)=-T_a.*T_n(j,i).*(E(j)-mu).^2.*Df(j);
end
L_11(i)=trapz(E,F1(:,i));
L_12(i)=trapz(E,F(:,i));
L_22(i)=trapz(E,F2(:,i));
L_21(i)=L_12(i);
ZT(i)=L_12(i).^2./(L_22(i).*L_11(i)-L_21(i).*L_12(i));


end

figure
plot(q,ZT,'LineWidth',2)
xlabel('q')
ylabel('ZT')






