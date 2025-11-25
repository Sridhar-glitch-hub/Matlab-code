%% -------------  8‑core parallel ZT calculation  -----------------
clc;  clear;

% -----------------  parameters  ----------------------------------
gamma   = 0.05;
gamma_p = 0.025;
t       = gamma;
T_s     = 0.08;
T_p     = 0.06;


epsilon = 2*gamma;
mu      = 0.01;
phi     =   pi ;
delta_values   = [0.5,0.6,0.7,0.8,0.9];

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
% Hybridization matrices for 5QD(3,2)
     G_L=zeros(6); G_R=zeros(6); G_P=zeros(6);
     G_L(1:3,1:3)=gamma; % source couples to dots 1,2,3
     G_R(4:6,4:6)=gamma; % drain couples to dots 4,5,6
     G_P(4,4)=gamma_p; % probe couples to dot 1,4

           

%% -----------------  preallocate result  -------------------------

ZT_all=zeros(nd,NT);

%% -----------------  main parallel loop over T_a -----------------
for jj=1:nd
    delta=delta_values(jj);
    ZT = zeros(1, NT);
parfor k = 1:NT
   
    Tk   = T_a(k);
    expo = exp((w - mu)./Tk);
    Df   = -expo ./ (Tk * (expo + 1).^2);   % row vector 1×Nw
    
                     % Transmission T(w)
    T_LR=zeros(1,Nw);
    T_LP=zeros(1,Nw);
    T_PL=zeros(1,Nw);
    T_PR=zeros(1,Nw);
    F   =zeros(1,Nw);
    F1   =zeros(1,Nw);
    F2   =zeros(1,Nw);
    F3   =zeros(1,Nw);
    F4   =zeros(1,Nw);
    F5   =zeros(1,Nw);
    F6   =zeros(1,Nw);
    F7   =zeros(1,Nw);
    F8   =zeros(1,Nw);
    F9   =zeros(1,Nw);
    F10   =zeros(1,Nw);
    F11   =zeros(1,Nw);
    for j = 1:Nw
        wj = w(j);
   
         
     
     % Retarded green function for 6QD(3,3)
                A = zeros(6);

                A(1,1)=w(j)-epsilon+delta+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(1,3)=1i*gamma/2;
                A(1,6)=-t*exp(-1i*phi/6);
                A(2,1)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(2,2)=w(j)-epsilon+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(3,1)=1i*gamma/2;
                A(3,2)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(3,3)=w(j)-epsilon-delta+1i*gamma/2;
                A(3,4)=-t*exp(1i*phi/6);
                A(4,3)=-t*exp(-1i*phi/6);
                A(4,4)=w(j)-epsilon+1i*gamma/2+1i*gamma_p/2;
                A(4,5)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(4,6)=1i*gamma/2;
                A(5,4)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(5,5)=w(j)-epsilon+1i*gamma/2;
                A(5,6)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(6,1)=-t*exp(1i*phi/6);
                A(6,4)=1i*gamma/2;
                A(6,5)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(6,6)=w(j)-epsilon+1i*gamma/2;
Gr=inv(A);
Ga=Gr';
% Transmission function for T_LR
T_LR(j)=trace(G_L*Gr*G_R*Ga);
T_LP(j)=trace(G_L*Gr*G_P*Ga);
T_PL(j)=trace(G_P*Gr*G_L*Ga);
T_PR(j)=trace(G_P*Gr*G_R*Ga);
    F1(j) = -  (T_LR(j)+T_LP(j)) .* Df(j);                 % L11 integrand
    F(j)  = -(T_PL(j)+T_PR(j)) .*Df(j) ;                    % L33 integrand
    F2(j) = -(-T_PL(j)) .*Df(j);                 % L13 integrand
    F3(j)=-(-T_LP(j)).*Df(j)  ;        % L31 integrand
    F4(j)=-(w(j)-mu).*(T_LR(j)+T_LP(j)).*Df(j);  % L12 or L21 integrand
    F5(j)=-(w(j)-mu).*(-T_LP(j)).*Df(j);      % L32 or L41 integrand
    F6(j)=-(w(j)-mu).*(-T_PL(j)).*Df(j);      % L14 or L23 integrand
    F7(j)=-(w(j)-mu).*(T_PL(j)+T_PR(j)).*Df(j);  % L43 or L34 integrand
    F8(j)=-(w(j)-mu).^2.*(T_PL(j)+T_PR(j)).*Df(j); % L44 integrand
    F9(j)=-(w(j)-mu).^2.*(-T_PL(j)).*Df(j);     % L24 integrand
    F10(j)=-(w(j)-mu).^2.*(-T_LP(j)).*Df(j);    % L42 integrand
    F11(j)=-(w(j)-mu).^2.*(T_LR(j)+T_LP(j)).*Df(j); % L22 integrand
    
      
    end
    
    % -------- transport integrals (row arithmetic) -----------------
   
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
   
    ZT(k)=X12.*X21./(X11.*X22-X12.*X21);
   
end
ZT_all(jj,:)=ZT;

end
figure; hold on;
for jj = 1:nd
    plot(T_a, ZT_all(jj, :), 'LineWidth', 2, 'DisplayName', ['\Delta = ' num2str(delta_values(jj))])
end
xlabel('T'); ylabel('ZT'); legend show; title('ZT vs T');



