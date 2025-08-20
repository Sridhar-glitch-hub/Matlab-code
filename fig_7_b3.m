clc;
clear;


gamma = 0.05;

T_s=8*gamma;
T_d=6*gamma;
epsilon = 6 * gamma;

phi = linspace(0, 10 * pi, 50); % Grid for phi
w = linspace(-2, 2, 100000); % Grid for omega
t = linspace(0.04, 0.1, 50); % Grid for t
mu = 4*gamma;

T_a = 6 * gamma;

%eta_c=1-T_d/T_s;

% Preallocate ZT matrix

ZT = zeros(length(t), length(phi));
eta= zeros(length(t), length(phi));

for i = 1:length(t)
    for k = 1:length(phi)
        % Preallocate arrays for calculations
        f = zeros(1, length(w));
        Df = zeros(1, length(w));
        n1 = zeros(1, length(w));
        d1 = zeros(1, length(w));
        T = zeros(1, length(w));
        F = zeros(1, length(w));
        F1 = zeros(1, length(w));
        F2 = zeros(1, length(w));
        
        for j = 1:length(w)
            f(j) = 1 / (exp((w(j) - mu) / T_a) + 1);
            n = exp((w(j) - mu) / T_a);
            d = T_a * (exp((w(j) - mu) / T_a) + 1)^2;
            Df(j) = -n / d;
            
            % Use t(i) and phi(k) for current iteration
           num(j)=4.*gamma.^2.*t(i).^2.*(cos((1/5).*phi(k)).*((-2).*t(i).^2+(epsilon+(-1).*w( ...
  j)).^2).*(epsilon+(-1).*w(j))+t(i).*(t(i).*((-1).*t(i).*cos((4/5).*phi(k))+2.* ...
  cos((3/5).*phi(k)).*(epsilon+(-1).*w(j)))+2.*cos((2/5).*phi(k)).*( ...
  epsilon+t(i)+(-1).*w(j)).*((-1).*epsilon+t(i)+w(j)))).^2;
    deno(j)=(1/8).*(16.*t(i).^10+10.*t(i).^8.*(9.*gamma.^2+20.*(epsilon+(-1).*w(j)) ...
  .^2)+10.*t(i).^6.*(gamma.^4+(-8).*gamma.^2.*(epsilon+(-1).*w(j)).^2+( ...
  -40).*(epsilon+(-1).*w(j)).^4)+2.*gamma.^2.*t(i).*cos((1/5).*phi(k)).*( ...
  28.*t(i).^6+(-188).*t(i).^4.*(epsilon+(-1).*w(j)).^2+t(i).^2.*(25.* ...
  gamma.^2+136.*(epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w(j)).^2+( ...
  -2).*(21.*gamma.^2+16.*(epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w( ...
  j)).^4).*(epsilon+(-1).*w(j))+5.*t(i).^4.*((-5).*gamma.^4+42.* ...
  gamma.^2.*(epsilon+(-1).*w(j)).^2+56.*(epsilon+(-1).*w(j)).^4).*( ...
  epsilon+(-1).*w(j)).^2+t(i).^2.*(25.*gamma.^4+(-128).*gamma.^2.*( ...
  epsilon+(-1).*w(j)).^2+(-80).*(epsilon+(-1).*w(j)).^4).*(epsilon+( ...
  -1).*w(j)).^4+2.*(gamma.^2+(epsilon+(-1).*w(j)).^2).*(9.*gamma.^2+ ...
  4.*(epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w(j)).^6+t(i).^2.*(28.* ...
  gamma.^2.*t(i).^6.*cos((8/5).*phi(k))+16.*t(i).^8.*cos(2.*phi(k))+gamma.^2.* ...
  t(i).^4.*cos((6/5).*phi(k)).*(gamma.^2+(-48).*t(i).^2+84.*(epsilon+(-1).*w( ...
  j)).^2)+gamma.^2.*t(i).^2.*cos((4/5).*phi(k)).*(100.*t(i).^4+6.*t(i).^2.*( ...
  gamma.^2+(-60).*(epsilon+(-1).*w(j)).^2)+((-5).*gamma.^2+148.*( ...
  epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w(j)).^2)+gamma.^2.*cos(( ...
  2/5).*phi(k)).*((-72).*t(i).^6+5.*t(i).^4.*(3.*gamma.^2+32.*(epsilon+(-1).* ...
  w(j)).^2)+(-80).*t(i).^2.*(gamma.^2+2.*(epsilon+(-1).*w(j)).^2).*( ...
  epsilon+(-1).*w(j)).^2+(85.*gamma.^2+52.*(epsilon+(-1).*w(j)).^2) ...
  .*(epsilon+(-1).*w(j)).^4)+(-2).*gamma.^2.*t(i).*cos((3/5).*phi(k)).*( ...
  4.*t(i).^4+(-5).*t(i).^2.*(gamma.^2+16.*(epsilon+(-1).*w(j)).^2)+3.*(5.* ...
  gamma.^2+12.*(epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w(j)).^2).*( ...
  epsilon+(-1).*w(j))+2.*t(i).^3.*cos(phi(k)).*(3.*gamma.^4+80.*t(i).^4+4.* ...
  t(i).^2.*(23.*gamma.^2+(-20).*(epsilon+(-1).*w(j)).^2)+(-64).* ...
  gamma.^2.*(epsilon+(-1).*w(j)).^2+16.*(epsilon+(-1).*w(j)).^4).*( ...
  epsilon+(-1).*w(j))+48.*gamma.^2.*t(i).^5.*cos((7/5).*phi(k)).*((-1).* ...
  epsilon+w(j))));
    T(j)=num(j)./deno(j);
            F(j) = -T_a.*T(j) * (w(j) - mu) * Df(j);
            F1(j) = -T_a.*T(j) * Df(j);
            F2(j) = -T_a.*T(j) * (w(j) - mu)^2 * Df(j);
        end
        
        % Integrate over w
        L_11 = trapz(w, F1);
        L_12 = trapz(w, F);
        L_22 = trapz(w, F2);
        L_21 = L_12;

        % Compute ZT for this phi and t combination
        ZT(i,k) = L_12^2 / (L_22 * L_11 - L_21 * L_12);
    %    eta(i,k)=((ZT(i,k)+1)^0.5-1)/((ZT(i,k)+1)^0.5+1);
    end
end


%pcolor(phi/pi, t/gamma, eta);
pcolor(phi/pi, t/gamma, ZT);
shading interp
ylabel('t/\gamma');
xlabel('\phi/\pi');
%zlabel('ZT');
colormap jet
colorbar
% Find maximum ZT and corresponding indices
[maxZT, maxLinearIndex] = max(ZT(:)); % Maximum value and linear index
[maxRowIndex, maxColIndex] = ind2sub(size(ZT), maxLinearIndex); % Convert to row/column indices

% Find minimum ZT and corresponding indices
[minZT, minLinearIndex] = min(ZT(:)); % Minimum value and linear index
[minRowIndex, minColIndex] = ind2sub(size(ZT), minLinearIndex); % Convert to row/column indices

% Corresponding t and phi values for maximum and minimum
t_max = t(maxRowIndex);
phi_max = phi(maxColIndex);

t_min = t(minRowIndex);
phi_min = phi(minColIndex);

% Display the results
fprintf('Maximum ZT: %.4f\n', maxZT);
fprintf('Corresponding t (max): %.4f\n', t_max);
fprintf('Corresponding phi (max): %.4f\n', phi_max);

fprintf('Minimum ZT: %.4f\n', minZT);
fprintf('Corresponding t (min): %.4f\n', t_min);
fprintf('Corresponding phi (min): %.4f\n', phi_min);



