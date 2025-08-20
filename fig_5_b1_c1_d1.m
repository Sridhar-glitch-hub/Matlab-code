%% -------------  8-core parallel I-Q-η calculation  --------------
clc;  clear;

%% ---------- parameters ----------
gamma   = 0.05;
t       = 2*gamma;
epsilon = 8 * gamma;
phi     = 0 ;
delta   = 0:0.01:0.1;                 % Nk = 11 values
T_s     = 12 * gamma;
T_d     =  2 * gamma;

v     = linspace(0, 0.6667, 1000);    % Nv = 1 000
mu_s  = -v/2;
mu_d  =  v/2;

b = -t * exp( 1i * phi / 4);
c = -t * exp(-1i * phi / 4);
d =  1i * gamma / 2;

%% ---------- Γ-matrices ----------
Gs          = zeros(4);      Gd = zeros(4);
Gs(1:2,1:2) = gamma;         % source couples to dots 1-3
Gd(3:4,3:4) = gamma;         % drain  couples to dots 4-6

%% ---------- energy grid ----------
w  = linspace(-2,  2, 100000);   % Nw = 100 000
Nw = numel(w);   Nv = numel(v);   Nk = numel(delta);

%% ---------- start / reuse 8-worker pool ----------
p = gcp('nocreate');
if isempty(p);               parpool('local', 8);
elseif p.NumWorkers ~= 8;     delete(p);  parpool('local', 8);
end

%% ---------- storage (sliceable) ----------
I   = zeros(Nv, Nk);
Q   = zeros(Nv, Nk);
P   = zeros(Nv, Nk);
eta = zeros(Nv, Nk);

P_max   = zeros(1, Nk);
eta_max = zeros(1, Nk);

%% ---------- SINGLE parfor over Δ -----------------
parfor k = 1:Nk
    Dk   = delta(k);
    Tvec = zeros(1, Nw);          % row for this Δ only

    % ----- loop over ω (cheap: 100k × 11 ≈ 1.1 M inverses) ----------
    for j = 1:Nw
        wj = w(j);

        a2 = wj - epsilon ;
        a3 = wj - epsilon+Dk;
        a4=  wj - epsilon ;
        a1 = wj - epsilon - Dk;
        A=[a1+d, b+d, 0, c;
           c+d, a2+d, b, 0;
           0, c, a3+d, b+d;
           b, 0, c+d, a4+d];

        %A = [ a1+d,  b+d,    d,      0, 0, c ;
        %      c+d,   a2+d,   b+d,    0, 0, 0 ;
         %     d,     c+d,    a3+d,   b, 0, 0 ;
         %     0,     0,      c,      a3+d, b+d, d ;
          %    0,     0,      0,      c+d, a2+d, b+d ;
          %    b,     0,      0,      d, c+d, a1+d ];

        Gr       = inv(A);
        Ga       = Gr';
        Tvec(j)  = real(trace(Gs * Gr * Gd * Ga));
    end

    % ----- loop over bias points v -------------------------------
    for ii = 1:Nv
        fs = 1 ./ (exp((w - mu_s(ii)) / T_s) + 1);
        fd = 1 ./ (exp((w - mu_d(ii)) / T_d) + 1);

        fvec      = Tvec .* (fs - fd);
        gvec      = (w - mu_s(ii)) .* fvec;

        Iii       = trapz(w, fvec);
        Qii       = trapz(w, gvec);
        Pii       = (mu_d(ii) - mu_s(ii)) * Iii;

        I(ii,k)   = Iii;
        Q(ii,k)   = Qii;
        P(ii,k)   = Pii;
        eta(ii,k) = Pii / Qii;
    end

   % P_max(k)   = max(P(:,k));
   % eta_max(k) = max(eta(:,k));
end
P_max   = max(P,   [], 1);   % 1-by-Nk    (do this on the client)
eta_max = max(eta, [], 1);

%% ---------- quick plot ------------------------------------------
eta_c = 1 - T_d / T_s;
%scatter(P_max, eta_max/eta_c, 60, delta, 'filled')
%xlabel('P_{max}');   ylabel('\eta_{max}/\eta_C')
%title('Peak efficiency vs. peak output power (colour = \delta)')
%colorbar
%grid on
yyaxis left
plot(delta,P_max)
ylabel('P_{max}')
xlabel('\Delta')
yyaxis right
plot(delta,eta_max/eta_c)
ylabel('\eta_{max}/\eta_c')
