%% -------------------------------------------------
% 1. Parameters (all broadcast to every worker)
% --------------------------------------------------
gamma   = 0.05;
epsilon = linspace(0.05,0.4,50);
mu      = 0.01;
T_a     = linspace(0.004,0.2,50);
w       = linspace(-2,2,100000);
t       = gamma;
phi     = pi;
delta   = 0.5;

nw      = numel(w);
Neps    = numel(epsilon);
NTa     = numel(T_a);

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

%% -------------------------------------------------
% 4. Parallel region
%    Each worker handles one value of epsilon (jj),
%    loops serially over T_a (k) and w (j)
% --------------------------------------------------
parfor jj = 1:Neps
    ZT_local = zeros(1,NTa);         % private copy for this worker

    for k = 1:NTa
        % temporary vectors for integration
        F  = zeros(1,nw);
        F1 = zeros(1,nw);
        F2 = zeros(1,nw);

        for j = 1:nw
            % --- Fermi function and its derivative
            n  = exp((w(j)-mu)/T_a(k));
            d  = T_a(k) * (n + 1)^2;
            Df = -n / d;
            
           A =zeros(4)
           A(1,1)=w(j)-epsilon(jj)+1i*gamma/2;
           A(1,2)=-t*exp(1i*phi/4)+1i*gamma/2;
       %    A(1,3)=0;
           A(1,4)=-t*exp(-1i*phi/4);
           A(2,1)=-t*exp(-1i*phi/4)+1i*gamma/2;
           A(2,2)=w(j)-epsilon(jj)+delta+1i*gamma/2;
           A(2,3)=-t*exp(1i*phi/4);
       %    A(2,4)=0;
       %    A(3,1)=0;
           A(3,2)=-t*exp(-1i*phi/4);
           A(3,3)=w(j)-epsilon(jj)+1i*gamma/2;
           A(3,4)=-t*exp(1i*phi/4)+1i*gamma/2;
           A(4,1)=-t*exp(1i*phi/4);
       %    A(4,2)=0;
           A(4,3)=-t*exp(-1i*phi/4)+1i*gamma/2;
           A(4,4)=w(j)-epsilon(jj)-delta+1i*gamma/2;
           
       
          
           % Define hybridization matrices for 4QD(2,2)
            Gs = zeros(4); Gd = zeros(4);
            Gs(1:2,1:2) = gamma; Gd(3:4,3:4) = gamma;
            
         

            % --- Transmission
            B      = inv(A);        % Retarded GF
            C      = B';            % Advanced  GF
            T_val  = trace(Gs*B*Gd*C);

            % --- Integrands
            F(j)  = -T_a(k)*T_val * (w(j) - mu)    * Df;
            F1(j) = -T_a(k)*T_val *                Df;
            F2(j) = -T_a(k)*T_val * (w(j) - mu).^2 * Df;
        end

        % --- Landauer integrals
        L11 = trapz(w, F1);
        L12 = trapz(w, F );
        L22 = trapz(w, F2);
        L21 = L12;

        % --- Figure of merit ZT
        ZT_local(k) = L12^2 / (L22*L11 - L21*L12);
    end

    ZT(jj,:) = ZT_local;            % write results back for this ε
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



