clc;
clear;
% y1 is maximum power and y2 is maximum eficiency for t less gamma
x=[2,3,4,5,6,8,10,12];
y1=[7.4e-3,1.2e-3,7.6247e-4,9.4e-4,6.7717e-4,3.96e-4,2.54e-4,1.769e-4];
y2=[0.5976,0.9626,0.87,0.9602,0.953,0.9712,0.9686,0.9669];
% y11 is maximum power and y22 is maximum eficiency for t eq gamma
y11=[7.4e-4,8.2e-3,10.5e-3,10.8e-3,10.2e-3,7.7e-3,5.2e-3,3.8e-3];
y22=[0.5976,0.8168,0.832,0.809,0.83,0.8581,0.8520,0.8472];
% y111 is maximum power and y222 is maximum eficiency for t gr gamma
y111=[7.4e-3,8.4e-3,16.8e-3,14e-3,17.6e-3,17.1e-3,13.5e-3,10.8e-3];
y222=[0.5976,0.6406,0.76,0.6306,0.685,0.72,0.7239,0.7283];
figure  % figure a1
plot(x,y1,'LineWidth',2.5)
xlabel('\Delta')
ylabel('P_{max}')
yyaxis right
plot(x,y2)
ylabel('\eta_{max}/\eta_c')
figure   % figure b1
plot(x,y11,'LineWidth',2.5)
xlabel('\Delta')
ylabel('P_{max}')
yyaxis right
plot(x,y22)
ylabel('\eta_{max}/\eta_c')
figure   % figure c1
plot(x,y111,'LineWidth',2.5)
xlabel('\Delta')
ylabel('P_{max}')
yyaxis right
plot(x,y222)
ylabel('\eta_{max}/\eta_c')