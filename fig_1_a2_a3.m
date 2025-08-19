clc;
clear;
Gamma=0.05;
epsilon=0.4;
q=200;
%T=1;
%phi=pi/2;
T_s=0.6;
T_d=0.1;
E=linspace(-5,5,100000);
v=linspace(0,0.5532,1000);

mu_s=-v/2;
mu_d=v/2;
for i=1:1:length(v)
   for j=1:1:length(E)
       fs(j,i)=1./(exp((E(j)-mu_s(i))/T_s)+1);
       fd(j,i)=1./(exp((E(j)-mu_d(i))/T_d)+1);
      % T_BW(j)=Gamma.^2/4./((E(j)-epsilon).^2+Gamma.^2/4);
       %e(j)=(E(j)-epsilon)./Gamma;
       n(j)=(0.5*Gamma*q+(E(j)-epsilon)).^2;
       d(j)=0.25*Gamma^2+(E(j)-epsilon).^2;
       T_F(j)=n(j)./d(j);
       T_F1(j)=T_F(j)./(1+q^2);

       %T(j)=T_BW(j)+T_F1(j);
       
       f(j,i)=T_F1(j).*(fs(j,i)-fd(j,i));
       %g(j,i)=(mu_d(i)-mu_s(i)).*f(j,i);
       g1(j,i)=(E(j)-mu_s(i)).*f(j,i);
   end
   I(i)=trapz(E,f(:,i));
   Q(i)=trapz(E,g1(:,i));
   P(i)=(mu_d(i)-mu_s(i)).*I(i);
   eta(i)=P(i)./Q(i);
end
%hold on
%title('\Gamma=0.05 ,T_{S}=0.6,T_{D}=0.1,\epsilon=0.4')
%xlabel('P') 
%ylabel('\eta/\eta_{c}')
eta_c=1-(T_d/T_s);
y=eta/eta_c;
%x=w-epsilon;
%plot(x,T)
%plot(P,y,LineWidth=1.5);