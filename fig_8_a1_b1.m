clc;
clear;
gamma=0.5;
T_a=1;
mu=4;
epsilon=3;
t=2*gamma;
w=linspace(mu-50*gamma,mu+50*gamma,100001);
phi=linspace(0,8*pi,801);
for i=1:1:length(phi)
    for j=1:1:length(w)
        f(j)=1./(exp((w(j)-mu)/T_a)+1); % Fermi distribution
   
    n(j)=exp((w(j)-mu)/T_a);
    d(j)=T_a.*(exp((w(j)-mu)/T_a)+1).^2;
    Df(j)=-n(j)/d(j);  % derivative of Fermi distribution
    num(j)= 4.*gamma.^2.*t.^2.*((epsilon+(-1).*w(j)).*cos((1/4).*phi(i))+2.*t.* ...
  sin((1/4).*phi(i)).^2).^2;
    deno(j)=(gamma.^2+(epsilon+(-1).*w(j)).^2+2.*t.^2+4.*(epsilon+(-1).*w(j)) ...
  .*t.*cos((1/4).*phi(i))+2.*t.^2.*cos((1/2).*phi(i))).*(gamma.^2.*( ...
  epsilon+(-1).*w(j)).^2+(epsilon+(-1).*w(j)).^4+(-4).*(epsilon+(-1) ...
  .*w(j)).^2.*t.^2+6.*t.^4+4.*t.^2.*((epsilon+(-1).*w(j)).^2+(-2).* ...
  t.^2).*cos((1/2).*phi(i))+2.*t.^4.*cos(phi(i)));
    T(j)=num(j)./deno(j);

    F(j)=-T_a.*T(j).*(w(j)-mu).*Df(j);
    F1(j)=-T_a.*T(j).*Df(j);
    F2(j)=-T_a.*T(j).*(w(j)-mu).^2.*Df(j);
    end
     % Transport coefficients
    L_11 = trapz(w, F1);
    L_12 = trapz(w, F);
    L_22 = trapz(w, F2);
    L_21 = L_12;

    
    S(i)  = -(1/T_a)*(L_12 / L_11);
   
    G(i)  = L_11/T_a;
    PF(i)=S(i).^2.*G(i);
    S1(i)=86.17.*S(i);
    PF1(i)=287.*PF(i);
end
figure
plot(phi/pi,S1,'LineWidth',2.5)
xlabel('\phi/\pi')
ylabel('S')
figure
plot(phi/pi,PF1,'LineWidth',2.5)
xlabel('\phi/\pi')
ylabel('PF')
