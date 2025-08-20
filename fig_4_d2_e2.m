clc;
clear;
T_s=0.6;
T_d=0.1;

gammad=0.05;
gammas=gammad;
t=0.05;
epsilon=0.4;
phi=5*pi/2;
w=linspace(-2,2,100000);
v=linspace(0,0.5637,1000);
mu_s=-v/2;
mu_d=v/2;
for i=1:1:length(v)
for j=1:1:length(w)
    fs(j,i)=1./(exp((w(j)-mu_s(i))/T_s)+1);
    fd(j,i)=1./(exp((w(j)-mu_d(i))/T_d)+1);
    num(j)=4.*gammad.*gammas.*t.^2.*(t.*((-1).*t.*cos((3/4).*phi)+cos((1/2).* ...
  phi).*(epsilon+(-1).*w(j)))+cos((1/4).*phi).*(epsilon+t+(-1).*w(j) ...
  ).*((-1).*epsilon+t+w(j))).^2;
    deno1(j)=(-2).*t.^4.*((-1)+cos(phi))+(1/4).*gammad.*gammas.*(2.*t.^2+(-2).* ...
  t.^2.*cos((1/2).*phi)+4.*t.*cos((1/4).*phi).*(epsilon+(-1).*w(j))+ ...
  (-3).*(epsilon+(-1).*w(j)).^2)+(-4).*t.^2.*(epsilon+(-1).*w(j)) ...
  .^2+(epsilon+(-1).*w(j)).^4;
    deno2(j)=2.*gammas.*t.^3.*cos((3/4).*phi)+(-1/2).*(gammad+3.*gammas).*((-2) ...
  .*t.^2+(epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w(j))+2.*gammas.* ...
  t.^2.*cos((1/2).*phi).*((-1).*epsilon+w(j))+(-2).*gammas.*t.*cos(( ...
  1/4).*phi).*(epsilon+t+(-1).*w(j)).*((-1).*epsilon+t+w(j));
    deno(j)=deno1(j).^2+deno2(j).^2;
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
ylabel('T_{4QD(3,1)}')