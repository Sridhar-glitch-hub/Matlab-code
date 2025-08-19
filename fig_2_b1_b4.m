clc;
clear;

%% Parameters
gamma   = 0.05;
epsilon = 2*gamma;
phi     = pi;
t       = 0.05;
w       = linspace(-2,2,100000);

% Range of delta values
deltas = 0.5:0.1:0.9;

%% Configurations to compute
configs = {'3QD21','4QD22','4QD31','5QD22','5QD32','6QD33'};

for c = 1:length(configs)
    config = configs{c};

    figure; hold on;
    legends = cell(length(deltas),1);

    for d = 1:length(deltas)
        delta = deltas(d);

        % Compute transmission
        T = compute_transmission(w,config,gamma,epsilon,phi,t,delta);

        % Plot
        plot(w-epsilon,T,'LineWidth',1.5);
        legends{d} = ['\Delta = ',num2str(delta)];
    end

    xlabel('\omega - \epsilon');
    ylabel('Transmission T(\omega)');
    title(['Transmission for ', config],'FontWeight','bold');
    legend(legends,'Location','best');
   % grid on;
    hold off;
end

%% ================= Helper Function =================
function T = compute_transmission(w,config,gamma,epsilon,phi,t,delta)
    switch config
        case '3QD21', N=3;
        case {'4QD22','4QD31'}, N=4;
        case {'5QD22','5QD32'}, N=5;
        case '6QD33', N=6;
    end

    T = zeros(size(w));
    Gs = zeros(N); Gd = zeros(N);

    switch config
        case '3QD21'
            Gs([1 3],[1 3]) = gamma;
            Gd(2,2) = gamma;

        case '4QD22'
            Gs(1:2,1:2) = gamma;
            Gd(3:4,3:4) = gamma;

        case '4QD31'
            Gs([1 2 4],[1 2 4]) = gamma;
            Gd(3,3) = gamma;

        case '5QD22'
            Gs(1:2,1:2) = gamma;
            Gd(4:5,4:5) = gamma;

        case '5QD32'
            Gs(1:3,1:3) = gamma;
            Gd(4:5,4:5) = gamma;

        case '6QD33'
            Gs(1:3,1:3) = gamma;
            Gd(4:6,4:6) = gamma;
    end

    for j = 1:length(w)
        E = w(j);
        A = zeros(N);

        % Hamiltonians (your manual fill-in cases stay the same)
        switch config
            case '3QD21'
                A(1,1)=E-epsilon-delta+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/3);
                A(1,3)=-t*exp(-1i*phi/3)+1i*gamma/2;
                A(2,1)=-t*exp(-1i*phi/3);
                A(2,2)=E-epsilon+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/3);
                A(3,1)=-t*exp(1i*phi/3)+1i*gamma/2;
                A(3,2)=-t*exp(-1i*phi/3);
                A(3,3)=E-epsilon+delta+1i*gamma/2;

            case '4QD22'
                A(1,1)=E-epsilon-delta+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/4)+1i*gamma/2;
                A(1,4)=-t*exp(-1i*phi/4);
                A(2,1)=-t*exp(-1i*phi/4)+1i*gamma/2;
                A(2,2)=E-epsilon+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/4);
                A(3,2)=-t*exp(-1i*phi/4);
                A(3,3)=E-epsilon+delta+1i*gamma/2;
                A(3,4)=-t*exp(1i*phi/4)+1i*gamma/2;
                A(4,1)=-t*exp(1i*phi/4);
                A(4,3)=-t*exp(-1i*phi/4)+1i*gamma/2;
                A(4,4)=E-epsilon+1i*gamma/2;

            case '4QD31'
                A(1,1)=E-epsilon+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/4)+1i*gamma/2;
                A(1,4)=-t*exp(-1i*phi/4)+1i*gamma/2;
                A(2,1)=-t*exp(-1i*phi/4)+1i*gamma/2;
                A(2,2)=E-epsilon-delta+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/4);
                A(2,4)=1i*gamma/2;
                A(3,2)=-t*exp(-1i*phi/4);
                A(3,3)=E-epsilon+1i*gamma/2;
                A(3,4)=-t*exp(1i*phi/4);
                A(4,1)=-t*exp(1i*phi/4)+1i*gamma/2;
                A(4,2)=1i*gamma/2;
                A(4,3)=-t*exp(-1i*phi/4);
                A(4,4)=E-epsilon+delta+1i*gamma/2;

            case '5QD22'
                A(1,1)=E-epsilon+delta+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/5)+1i*gamma/2;
                A(1,5)=-t*exp(-1i*phi/5);
                A(2,1)=-t*exp(-1i*phi/5)+1i*gamma/2;
                A(2,2)=E-epsilon+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/5);
                A(3,2)=-t*exp(-1i*phi/5);
                A(3,3)=E-epsilon-delta;
                A(3,4)=-t*exp(1i*phi/5);
                A(4,3)=-t*exp(-1i*phi/5);
                A(4,4)=E-epsilon+1i*gamma/2;
                A(4,5)=-t*exp(1i*phi/5)+1i*gamma/2;
                A(5,1)=-t*exp(1i*phi/5);
                A(5,4)=-t*exp(-1i*phi/5)+1i*gamma/2;
                A(5,5)=E-epsilon+1i*gamma/2;

            case '5QD32'
                A(1,1)=E-epsilon+delta+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/5)+1i*gamma/2;
                A(1,3)=1i*gamma/2;
                A(1,5)=-t*exp(-1i*phi/5);
                A(2,1)=-t*exp(-1i*phi/5)+1i*gamma/2;
                A(2,2)=E-epsilon+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/5)+1i*gamma/2;
                A(3,1)=1i*gamma/2;
                A(3,2)=-t*exp(-1i*phi/5)+1i*gamma/2;
                A(3,3)=E-epsilon-delta+1i*gamma/2;
                A(3,4)=-t*exp(1i*phi/5);
                A(4,3)=-t*exp(-1i*phi/5);
                A(4,4)=E-epsilon+1i*gamma/2;
                A(4,5)=-t*exp(1i*phi/5)+1i*gamma/2;
                A(5,1)=-t*exp(1i*phi/5);
                A(5,4)=-t*exp(-1i*phi/5)+1i*gamma/2;
                A(5,5)=E-epsilon+1i*gamma/2;

            case '6QD33'
                A(1,1)=E-epsilon+delta+1i*gamma/2;
                A(1,2)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(1,3)=1i*gamma/2;
                A(1,6)=-t*exp(-1i*phi/6);
                A(2,1)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(2,2)=E-epsilon+1i*gamma/2;
                A(2,3)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(3,1)=1i*gamma/2;
                A(3,2)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(3,3)=E-epsilon-delta+1i*gamma/2;
                A(3,4)=-t*exp(1i*phi/6);
                A(4,3)=-t*exp(-1i*phi/6);
                A(4,4)=E-epsilon+1i*gamma/2;
                A(4,5)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(4,6)=1i*gamma/2;
                A(5,4)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(5,5)=E-epsilon+1i*gamma/2;
                A(5,6)=-t*exp(1i*phi/6)+1i*gamma/2;
                A(6,1)=-t*exp(1i*phi/6);
                A(6,4)=1i*gamma/2;
                A(6,5)=-t*exp(-1i*phi/6)+1i*gamma/2;
                A(6,6)=E-epsilon+1i*gamma/2;
        end

        % Transmission
        Gr = inv(A); Ga = Gr';
        T(j) = real(trace(Gs * Gr * Gd * Ga));
    end
end
