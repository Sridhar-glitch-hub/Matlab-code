%% -------------------------------------------------
% 1. Parameters (all broadcast to every worker)
% --------------------------------------------------
gamma   = 0.05;
gamma_p =0.025;
epsilon = linspace(0.05,0.4,50);
mu      = 0.01;
T_a     = linspace(0.004,0.2,50);
w       = linspace(-2,2,10000);
t       = gamma;
phi     = pi;
delta   = 0.9;
T_s     = 0.08;
T_p     = 0.06;

nw      = numel(w);
Neps    = numel(epsilon);
NTa     = numel(T_a);
G_L=zeros(6); G_R=zeros(6); G_P=zeros(6);
     G_L(1:3,1:3)=gamma; % source couples to dots 1,2,3
     G_R(4:6,4:6)=gamma; % drain couples to dots 4,5,6
     G_P(4,4)=gamma_p; % probe couples to dot 4


%% -------------------------------------------------
% 2. Spin‑up a pool with 8 workers (if one isn’t open)
% --------------------------------------------------
p = gcp('nocreate');  % Check for existing pool
if isempty(p)
    parpool('local', 8);  % Start new one
elseif p.NumWorkers ~= 8
    delete(p);            % Delete current pool if not using 8 workers
    parpool('local', 8);  % Restart with 8 workers
end


%% -------------------------------------------------
% 3. Pre‑allocate the output (shared)
% --------------------------------------------------
ZT  = zeros(Neps,NTa);



parfor jj = 1:Neps
    ZT_local = zeros(1,NTa);
    
    for k = 1:NTa
        % temporary vectors for integration
        F  = zeros(1,nw);
        F1 = zeros(1,nw);
        F2 = zeros(1,nw);
        F3 = zeros(1,nw);
        F4 = zeros(1,nw);
        F5 = zeros(1,nw);
        F6 = zeros(1,nw);
        F7 = zeros(1,nw);
        F8 = zeros(1,nw);
        F9 = zeros(1,nw);
        F10 = zeros(1,nw);
        F11 = zeros(1,nw);










        for j = 1:nw
            % --- Fermi function and its derivative
            n  = exp((w(j)-mu)/T_a(k));
            d  = T_a(k) * (n + 1)^2;
            Df = -n / d;
            
          

 % Retarded green function for 6QD(3,3) with voltage probe
                A = zeros(6);

                A(1,1)=w(j)-epsilon(jj)+delta+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(1,3)=1i*gamma/2;
                A(1,6)=-t*exp(-1i*phi/6);
                A(2,1)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(2,2)=w(j)-epsilon(jj)+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(3,1)=1i*gamma/2;
                A(3,2)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(3,3)=w(j)-epsilon(jj)-delta+1i*gamma/2;
                A(3,4)=-t*exp(1i*phi/6);
                A(4,3)=-t*exp(-1i*phi/6);
                A(4,4)=w(j)-epsilon(jj)+1i*gamma/2+1i*gamma_p/2;
                A(4,5)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(4,6)=1i*gamma/2;
                A(5,4)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(5,5)=w(j)-epsilon(jj)+1i*gamma/2;
                A(5,6)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(6,1)=-t*exp(1i*phi/6);
                A(6,4)=1i*gamma/2;
                A(6,5)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(6,6)=w(j)-epsilon(jj)+1i*gamma/2;
Gr=inv(A);
Ga=Gr';
% Transmission function for T_LR
T_LR=trace(G_L*Gr*G_R*Ga);
T_LP=trace(G_L*Gr*G_P*Ga);
T_PL=trace(G_P*Gr*G_L*Ga);
T_PR=trace(G_P*Gr*G_R*Ga);
     
                        
       
          
           
            
         

         
         F1(j) = -  (T_LR+T_LP) .* Df;                 % L11 integrand
    F(j)  = -(T_PL+T_PR) .*Df ;                    % L33 integrand
    F2(j) = -(-T_PL) .*Df;                 % L13 integrand
    F3(j)=-(-T_LP).*Df          % L31 integrand
    F4(j)=-(w(j)-mu).*(T_LR+T_LP).*Df;  % L12 or L21 integrand
    F5(j)=-(w(j)-mu).*(-T_LP).*Df;      % L32 or L41 integrand
    F6(j)=-(w(j)-mu).*(-T_PL).*Df;      % L14 or L23 integrand
    F7(j)=-(w(j)-mu).*(T_PL+T_PR).*Df;  % L43 or L34 integrand
    F8(j)=-(w(j)-mu).^2.*(T_PL+T_PR).*Df; % L44 integrand
    F9(j)=-(w(j)-mu).^2.*(-T_PL).*Df;     % L24 integrand
    F10(j)=-(w(j)-mu).^2.*(-T_LP).*Df;    % L42 integrand
    F11(j)=-(w(j)-mu).^2.*(T_LR+T_LP).*Df; % L22 integrand
        end

        % --- Landauer integrals
       
    L11 = trapz(w, F1);
    L33 = trapz(w, F);
    L13 = trapz(w, F2);
    L31 =trapz(w,F3);
    L12 =trapz(w,F4);
    L21=L12;
    L32=trapz(w,F5);
    L41=L32;
    L14=trapz(w,F6);
    L23=L14;
    L43=trapz(w,F7);
    L34=L43;
    L44=trapz(w,F8);
    L24=trapz(w,F9);
    L42=trapz(w,F10);
    L22=trapz(w,F11);
    X11=(L11.*L33-L13.*L31)./L33;
    X12=(L12.*L33-L13.*L32)./L33;
    X13=(L14.*L33-L13.*L34)./L33;
    X21=(L21.*L33-L23.*L31)./L33;
    X31=(L41.*L33-L43.*L31)./L33;
    X22=(L22.*L33-L23.*L32)./L33;
    X33=(L44.*L33-L43.*L34)./L33;
    X23=(L24.*L33-L23.*L34)./L33;
    X32=(L42.*L33-L43.*L32)./L33;
   
    ZT_local(k)=X12.*X21./(X11.*X22-X12.*X21);

   
    end
    ZT(jj,:) = ZT_local;
   
end

%% -------------------------------------------------
% 5. Plot
% --------------------------------------------------
figure
pcolor(T_a, epsilon, real(ZT));
shading interp
xlabel('T');
ylabel('\epsilon');
colormap jet
colorbar




