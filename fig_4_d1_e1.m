clc;
clear;
T_s = 0.6;
T_d = 0.1;

gamma_d = 0.1;
gamma_s=3*gamma_d/4;
%gamma_s=0.0125:0.0125:0.2;
%alpha=gamma_s./gamma_d;
epsilon = 0.4;
phi = 2*pi;
t = 0.05;
w = linspace(-2, 2, 100000);
v = linspace(0, 0.5517, 1000);
mu_s = -v / 2;
mu_d = v / 2;

%I = zeros(1, length(v)); % Preallocate memory
%Q = zeros(1, length(v));
%P = zeros(1, length(v));
%eta = zeros(1, length(v));


  % for k=1:length(gamma_s)
   for i = 1:length(v)
   for j = 1:length(w)
    fs(j,i) = 1 ./ (exp((w(j) - mu_s(i)) / T_s) + 1);
    fd(j,i) = 1 ./ (exp((w(j) - mu_d(i)) / T_d) + 1);
    
   
    
        a = [(w(j) - epsilon + 1i/2 * gamma_s), -t * exp(1i/4 * phi) + 1i/2 * gamma_s, 0, -t * exp(-1i/4 * phi);
             -t * exp(-1i/4 * phi) + 1i/2 * gamma_s, (w(j) - epsilon + 1i/2 * gamma_s), -t * exp(1i/4 * phi), 0;
             0, -t * exp(-1i/4 * phi), (w(j) - epsilon + 1i/2 * gamma_d), -t * exp(1i/4 * phi) + 1i/2 * gamma_d;
             -t * exp(1i/4 * phi), 0, -t * exp(-1i * phi/4) + 1i/2 * gamma_d, (w(j) - epsilon + 1i/2 * gamma_d)];
        
        b = inv(a);
        c = b';
        x = [gamma_s, gamma_s, 0, 0;
             gamma_s, gamma_s, 0, 0;
             0, 0, 0, 0;
             0, 0, 0, 0];

        y = [0, 0, 0, 0;
             0, 0, 0, 0;
             0, 0, gamma_d, gamma_d;
             0, 0, gamma_d, gamma_d];

        n = x * b * y * c;
        T(j)=trace(n);
        %T(j,k) = trace(n);
       f(j,i)=T(j).* (fs(j,i) - fd(j,i));
        %f(j,k) = T(j,k) .* (fs(j,i) - fd(j,i));
        g(j,i) = (w(j) - mu_s(i)) .* f(j,i);
    end
    
    %I(i,k) = trapz(w, f(:,k));
    I(i)=trapz(w,f(:,i));
    Q(i)=trapz(w,g(:,i));
    P(i)=(mu_d(i) - mu_s(i)) .* I(i);
    eta(i)=P(i)./Q(i);

   % Q(i,k) = trapz(w, g);
    %P(i,k) = (mu_d(i) - mu_s(i)) .* I(i,k);

    %eta(i,k) = P(i,k) / Q(i,k);
 %  end
  % P_max(k)=max(P(:,k));
   end


eta_c = 1 - T_d / T_s;
y = eta / eta_c;
figure
plot(P,y,'LineWidth',2.5)
xlabel('P')
ylabel('\eta/\eta_c')
figure
plot(w-epsilon,T,'LineWidth',2.5)
xlabel('\omega-\epsilon')
ylabel('T_{4QD(2,2)}')
