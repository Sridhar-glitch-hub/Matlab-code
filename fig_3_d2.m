clc;
clear;

gamma   = 0.05;
t=linspace(0.01,0.1,100);

epsilon=2*gamma;
delta = 0.5;
phi     = pi;
w       = linspace(-2, 2, 1e5);  % 10^5 points for faster convergence
mu      = 0.01;

T_a=0.0334;

% Hybridzation matrices for 4QD(3,1)
     Gd=zeros(4); Gs=zeros(4);
     Gd(3,3)=gamma;
     Gs(1,1)=gamma;Gs(1,2)=gamma;Gs(1,4)=gamma;
     Gs(2,1)=gamma;Gs(2,2)=gamma;Gs(2,4)=gamma;
     Gs(4,1)=gamma;Gs(4,2)=gamma;Gs(4,4)=gamma;


nt = numel(t);




% Start parallel pool if needed
p = gcp('nocreate');
if isempty(p) || p.NumWorkers ~= 8
    delete(gcp('nocreate'));
    parpool('local', 8);
end



    ZT_k = zeros(1, nt);
    S_k  = zeros(1, nt);
    K_k  = zeros(1, nt);
    G_k  = zeros(1, nt);
    L_k  = zeros(1, nt);

    parfor k = 1:nt
       
        t_k = t(k);
        nw = length(w);

        F   = zeros(1, nw);
        F1  = zeros(1, nw);
        F2  = zeros(1, nw);
        Tvec = zeros(1, nw);

        for j = 1:nw
            E = w(j);
            f = 1 / (exp((E - mu) / T_a) + 1);
            n = exp((E - mu) / T_a);
            dfdT = T_a * (exp((E - mu) / T_a) + 1)^2;
            Df = -n / dfdT;

   
          % Transmission function for 4QD(3,1) 

     A=zeros(4);
           A(1,1)=w(j)-epsilon+1i*gamma/2;
           A(1,2)=-t_k*exp(1i*phi/4)+1i*gamma/2;
      %     A(1,3)=0;
           A(1,4)=-t_k*exp(-1i*phi/4)+1i*gamma/2;
           A(2,1)=-t_k*exp(-1i*phi/4)+1i*gamma/2;
           A(2,2)=w(j)-epsilon-delta+1i*gamma/2;
           A(2,3)=-t_k*exp(1i*phi/4);
           A(2,4)=1i*gamma/2;
      %     A(3,1)=0;
           A(3,2)=-t_k*exp(-1i*phi/4);
           A(3,3)=w(j)-epsilon+1i*gamma/2;
           A(3,4)=-t_k*exp(1i*phi/4);
           A(4,1)=-t_k*exp(1i*phi/4)+1i*gamma/2;
           A(4,2)=1i*gamma/2;
           A(4,3)=-t_k*exp(-1i*phi/4);
           A(4,4)=w(j)-epsilon+delta+1i*gamma/2;



     

            Gr = inv(A); Ga = Gr';
            Tj = real(trace(Gs * Gr * Gd * Ga));
            Tvec(j) = Tj;

            F(j)  = - Tj * (E - mu) * Df;
            F1(j) = - Tj * Df;
            F2(j) = - Tj * (E - mu)^2 * Df;
        end

        % Transport coefficients
        L_11 = trapz(w, F1);
        L_12 = trapz(w, F);
        L_22 = trapz(w, F2);
        L_21 = L_12;

        ZT_k(k) = (L_12^2) / (L_22 * L_11 - L_21 * L_12);
       % S_k(k)  = -(1/T_a)*(L_12 / L_11);
      %  K_k(k)  = (1/T_a)*((L_22 * L_11 - L_12^2) / L_11);
      %  G_k(k)  = L_11;
      %  L_k(k)  = (1/T_a)*(K_k(k) / G_k(k));
    end

 

%% ---- PLOT RESULTS ----
L0 = pi^2 / 3;

plot(t/gamma,ZT_k,'LineWidth',2)
xlabel('t/\gamma');
ylabel('ZT')