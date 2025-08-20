clc;
clear;
gamma=0.05;
epsilon=8*gamma;
T_s=12*gamma;
T_d=2*gamma;
t=2*gamma;
phi=0;
w=linspace(-2,2,100000);
v=linspace(0,0.6667,1000);
mu_s=-v/2;
mu_d=v/2;
%for k=1:1:length(t)
    %for m=1:1:length(phi)
        for i=1:1:length(v)
           for j=1:1:length(w)
               
        
        fs(j,i)=1./(exp((w(j)-mu_s(i))/T_s)+1);
        fd(j,i)=1./(exp((w(j)-mu_d(i))/T_d)+1);
        num(j)= 4.*gamma.^2.*t.^2.*((epsilon+(-1).*w(j)).*cos((1/4).*phi)+2.*t.* ...
  sin((1/4).*phi).^2).^2;
    deno(j)=(gamma.^2+(epsilon+(-1).*w(j)).^2+2.*t.^2+4.*(epsilon+(-1).*w(j)) ...
  .*t.*cos((1/4).*phi)+2.*t.^2.*cos((1/2).*phi)).*(gamma.^2.*( ...
  epsilon+(-1).*w(j)).^2+(epsilon+(-1).*w(j)).^4+(-4).*(epsilon+(-1) ...
  .*w(j)).^2.*t.^2+6.*t.^4+4.*t.^2.*((epsilon+(-1).*w(j)).^2+(-2).* ...
  t.^2).*cos((1/2).*phi)+2.*t.^4.*cos(phi));
    T(j)=num(j)./deno(j);
    f(j,i)=T(j).*(fs(j,i)-fd(j,i));
    g(j,i)=(w(j)-mu_s(i)).*f(j,i);
                  end
               
    
    I(i)=trapz(w,f(:,i));
    Q(i)=trapz(w,g(:,i));
    P(i)=(mu_d(i)-mu_s(i)).*I(i);
    eta(i)=P(i)./Q(i);
        end
 eta_c=1-T_d/T_s;
 y=eta/eta_c;
figure       
plot(P,y,'LineWidth',2.5)
xlabel('P')
ylabel('\eta/\eta_c')
figure
plot(w-epsilon,T,'LineWidth',2.5)
xlabel('\omega-\epsilon')
ylabel('T_{4QD(2,2)}')