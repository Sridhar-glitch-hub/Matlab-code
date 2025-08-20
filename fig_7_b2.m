clc;
clear;

gamma = 0.05;
epsilon = 6 * gamma;

phi = linspace(0, 8 * pi, 50); % Grid for phi
w = linspace(-2, 2, 100000); % Grid for omega
%t = linspace(0.04, 0.125, 50); % Grid for t
t = linspace(0.01, 0.1, 50);
mu = 4 * gamma;

T_a = 6 * gamma;


% Preallocate ZT matrix
ZT = zeros(length(t), length(phi)); 

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
            n1(j)=4.*gamma.^2.*t(i).^2.*((epsilon+(-1).*w(j)+t(i)).*((-1).*epsilon+w(j)+t(i)) ...
         .*cos((1/4).*phi(k))+t(i).*((epsilon+(-1).*w(j)).*cos((1/2).*phi(k))+(-1).* ...
         t(i).*cos((3/4).*phi(k)))).^2;
          d1(j)=(1/16).*((gamma.^2+4.*(epsilon+(-1).*w(j)).^2).*(9.*gamma.^2+4.*( ...
  epsilon+(-1).*w(j)).^2).*(epsilon+(-1).*w(j)).^4+(-4).*(gamma.^4+ ...
  28.*gamma.^2.*(epsilon+(-1).*w(j)).^2+32.*(epsilon+(-1).*w(j)).^4) ...
  .*(epsilon+(-1).*w(j)).^2.*t(i).^2+2.*(3.*gamma.^4+56.*gamma.^2.*( ...
  epsilon+(-1).*w(j)).^2+160.*(epsilon+(-1).*w(j)).^4).*t(i).^4+32.*( ...
  3.*gamma.^2+(-8).*(epsilon+(-1).*w(j)).^2).*t(i).^6+96.*t(i).^8+(-8).* ...
  gamma.^2.*(epsilon+(-1).*w(j)).*t(i).*(3.*(gamma.^2+4.*(epsilon+(-1) ...
  .*w(j)).^2).*(epsilon+(-1).*w(j)).^2+(-1).*(gamma.^2+24.*(epsilon+ ...
  (-1).*w(j)).^2).*t(i).^2+24.*t(i).^4).*cos((1/4).*phi(k))+2.*t(i).^2.*((-2).* ...
  gamma.^2.*((-1).*(5.*gamma.^2+36.*(epsilon+(-1).*w(j)).^2).*( ...
  epsilon+(-1).*w(j)).^2+2.*(gamma.^2+24.*(epsilon+(-1).*w(j)).^2).* ...
  t(i).^2+12.*t(i).^4).*cos((1/2).*phi(k))+(-4).*gamma.^2.*(epsilon+(-1).*w( ...
  j)).*t(i).*(gamma.^2+24.*(epsilon+(-1).*w(j)).^2+(-36).*t(i).^2).*cos(( ...
  3/4).*phi(k))+t(i).^2.*(gamma.^4+72.*gamma.^2.*(epsilon+(-1).*w(j)).^2+( ...
  -32).*(epsilon+(-1).*w(j)).^4+(-16).*(3.*gamma.^2+(-8).*(epsilon+( ...
  -1).*w(j)).^2).*t(i).^2+(-64).*t(i).^4).*cos(phi(k))+48.*gamma.^2.*((-1).* ...
  epsilon+w(j)).*t(i).^3.*cos((5/4).*phi(k))+24.*gamma.^2.*t(i).^4.*cos((3/2) ...
  .*phi(k))+16.*t(i).^6.*cos(2.*phi(k))));
            
            T(j) = n1(j) / d1(j);
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
        ZT(i, k) = L_12^2 / (L_22 * L_11 - L_21 * L_12);
    end
end

pcolor(phi/pi, t/gamma, ZT);
shading interp
ylabel('t/\gamma');
xlabel('\phi/\pi');
%zlabel('ZT');
colormap jet
colorbar;
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


