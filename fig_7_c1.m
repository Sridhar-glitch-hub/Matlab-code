clc;
clear;
gamma=0.05;
epsilon=0.4;
t=0.8*gamma;
phi=2.9*pi;
w=linspace(-2,2,100000);
for j=1:1:length(w)
    num(j)= 4.*gamma.^2.*t.^2.*((epsilon+(-1).*w(j)).*cos((1/4).*phi)+2.*t.* ...
  sin((1/4).*phi).^2).^2;
    deno(j)=(gamma.^2+(epsilon+(-1).*w(j)).^2+2.*t.^2+4.*(epsilon+(-1).*w(j)) ...
  .*t.*cos((1/4).*phi)+2.*t.^2.*cos((1/2).*phi)).*(gamma.^2.*( ...
  epsilon+(-1).*w(j)).^2+(epsilon+(-1).*w(j)).^4+(-4).*(epsilon+(-1) ...
  .*w(j)).^2.*t.^2+6.*t.^4+4.*t.^2.*((epsilon+(-1).*w(j)).^2+(-2).* ...
  t.^2).*cos((1/2).*phi)+2.*t.^4.*cos(phi));
    T(j)=num(j)./deno(j);
end
plot(w-epsilon,T,'LineWidth',2.5)