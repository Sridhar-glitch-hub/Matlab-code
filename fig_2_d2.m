%% -------------  8‑core parallel ZT calculation  -----------------
clc;  clear;

% -----------------  parameters  ----------------------------------
gamma   = 0.05;
t       = gamma;

epsilon = 2*gamma;
mu      = 0.01;
phi     =   pi ;
delta_values   = 0.5:0.1:0.9;

w     = linspace(-2, 2, 10000);   % energy grid
T_a   = linspace(0, 0.2, 1000);    % temperature grid
Nw    = numel(w);
NT    = numel(T_a);
nd    = numel(delta_values);

%% -----------------  start / reuse 8‑worker pool  ----------------
p = gcp('nocreate');
if isempty(p)
    parpool('local', 8);
elseif p.NumWorkers ~= 8
    delete(p);
    parpool('local', 8);
end

%% -----------------  constant Γ‑matrices  ------------------------
% Hybridization matrices for 4QD(3,1)
Gd=zeros(4); Gs=zeros(4);
     Gd(3,3)=gamma;
     Gs(1,1)=gamma;Gs(1,2)=gamma;Gs(1,4)=gamma;
     Gs(2,1)=gamma;Gs(2,2)=gamma;Gs(2,4)=gamma;
     Gs(4,1)=gamma;Gs(4,2)=gamma;Gs(4,4)=gamma;


           

%% -----------------  preallocate result  -------------------------
ZT_all = zeros(nd, NT);

%% -----------------  main parallel loop over T_a -----------------
for jj =1:nd
    ZT=zeros(1,NT);
    delta=delta_values(jj);
parfor k = 1:NT
    Tk   = T_a(k);
    expo = exp((w - mu)./Tk);
    Df   = -expo ./ (Tk * (expo + 1).^2);   % row vector 1×Nw
    
    Tvec = zeros(1, Nw);                    % Transmission T(w)
    
    for j = 1:Nw
        wj = w(j);
        
         % Transmission function for 4QD(3,1) 
     
     A=zeros(4);
           A(1,1)=w(j)-epsilon+1i*gamma/2;
           A(1,2)=-t*exp(1i*phi/4)+1i*gamma/2;
      %     A(1,3)=0;
           A(1,4)=-t*exp(-1i*phi/4)+1i*gamma/2;
           A(2,1)=-t*exp(-1i*phi/4)+1i*gamma/2;
           A(2,2)=w(j)-epsilon-delta+1i*gamma/2;
           A(2,3)=-t*exp(1i*phi/4);
           A(2,4)=1i*gamma/2;
      %     A(3,1)=0;
           A(3,2)=-t*exp(-1i*phi/4);
           A(3,3)=w(j)-epsilon+1i*gamma/2;
           A(3,4)=-t*exp(1i*phi/4);
           A(4,1)=-t*exp(1i*phi/4)+1i*gamma/2;
           A(4,2)=1i*gamma/2;
           A(4,3)=-t*exp(-1i*phi/4);
           A(4,4)=w(j)-epsilon+delta+1i*gamma/2;
    
        Gr        = inv(A);
        Ga        = Gr';
        Tvec(j)   = real(trace(Gs * Gr * Gd * Ga));      % scalar
    end
    
    % -------- transport integrals (row arithmetic) -----------------
   % F1 = -Tk .* Tvec .* Df; % L11 integrand
    F1 = - Tvec .* Df;
    F  = F1 .* (w - mu);                    % L12 integrand
    F2 = F1 .* (w - mu).^2;                 % L22 integrand
    
    L11 = trapz(w, F1);
    L12 = trapz(w, F);
    L22 = trapz(w, F2);
    
    ZT(k) = L12^2 / (L22 * L11 - L12^2);
   % R(k)=(3/(pi*T_a(k)).^2).*(L22/L11-(L12/L11)^2)
end
ZT_all(jj,:)=ZT;
end


%% -----------------  plot result  --------------------------------

figure; hold on;
for jj = 1:nd
    plot(T_a, ZT_all(jj, :), 'LineWidth', 2, 'DisplayName', ['\Delta = ' num2str(delta_values(jj))])
end
xlabel('T'); ylabel('ZT'); legend show; title('ZT vs T');
%hold on
%plot(T_a, ZT, 'LineWidth', 2.5);
%xlabel('T');
%ylabel('ZT');


