clc;
clear;
T_s=0.6;
T_d=0.1;

gammad=0.05;
%gammas=0.025:0.025:0.2;
%alpha=gammas./gammad;
gammas=gammad;
phi=5*pi/2;
epsilon=0.4;
t=0.01;
w=linspace(-2,2,100000);
v=linspace(0,0.5716,1000);
mu_s=-v/2;
mu_d=v/2;
%for k=1:1:length(gammas)
for i=1:1:length(v)
for j=1:1:length(w)
    fs(j,i)=1./(exp((w(j)-mu_s(i))/T_s)+1);
    fd(j,i)=1./(exp((w(j)-mu_d(i))/T_d)+1);
    num(j)=4.*gammad.*gammas.*t.^2.*(t.*(t.*cos((1/2).*phi)+2.*sin((1/6).* ...
  phi).^2.*(epsilon+(-1).*w(j)))+cos((1/6).*phi).*((-2).*t.^2+( ...
  epsilon+(-1).*w(j)).^2)).^2;
    deno(j)=2.*gammad.*gammas.*t.^2.*(t.*(t.*cos((1/2).*phi)+2.*sin((1/6).* ...
  phi).^2.*(epsilon+(-1).*w(j)))+cos((1/6).*phi).*((-2).*t.^2+( ...
  epsilon+(-1).*w(j)).^2)).^2+(1/4).*gammas.^2.*(4.*t.^3.*cos((1/2) ...
  .*phi)+2.*t.*cos((1/6).*phi).*(2.*t.^2+(-1).*(epsilon+(-1).*w(j)) ...
  .^2)+(7.*t.^2+2.*t.^2.*cos((1/3).*phi)+(-3).*(epsilon+(-1).*w(j)) ...
  .^2).*(epsilon+(-1).*w(j))).^2+((-2).*t.^3.*cos((1/2).*phi)+((-3) ...
  .*t.^2+(epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w(j))).^2.*( ...
  epsilon+2.*t.*cos((1/6).*phi)+(-1).*w(j)).^2+(1/16).*gammad.^2.*( ...
  4.*(4.*t.^3.*cos((1/2).*phi)+2.*t.*cos((1/6).*phi).*(2.*t.^2+(-1) ...
  .*(epsilon+(-1).*w(j)).^2)+(7.*t.^2+2.*t.^2.*cos((1/3).*phi)+(-3) ...
  .*(epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w(j))).^2+9.* ...
  gammas.^2.*(2.*t.*(t.*cos((1/3).*phi)+cos((1/6).*phi).*(epsilon+( ...
  -1).*w(j)))+3.*(epsilon+t+(-1).*w(j)).*((-1).*epsilon+t+w(j))).^2) ...
  ;
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