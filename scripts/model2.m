clear; close all; clc;

fs = 5e3;

%% --------------------------------------------------------
% 1) DATA ACQUISITION
%% --------------------------------------------------------
[uFull, yFull, ~, sig, ~, ~] = acquisition('fastMethod1/odd_5k');
%[uFull, yFull, ~, sig, ~, ~] = acquisition('fastMethodEven/even');

periodN = size(sig,1);

u = uFull(end-periodN+1:end,1);
y = yFull(end-periodN+1:end,1);

u = u(:); 
y = y(:);
N = length(u);

%% --------------------------------------------------------
% 2) MODEL ORDER
%% --------------------------------------------------------
Norder = 5;       % polynomial order
Nlag   = 10;       % memory depth

Nh = Nlag + 1;    
Nb = Nlag;       
Na = Norder - 1;   % nonlinear feedback powers y^2 ... y^P

%% --------------------------------------------------------
% 3) PARAMETER INITIALIZATION
%% --------------------------------------------------------
h0 = zeros(Nh,1);  h0(1) = 1;    % FIR feedforward
b0 = zeros(Nb,1);              % IIR poles
a0 = zeros(Na,1);              % nonlinear feedback

theta0 = [h0; b0; a0];

%% --------------------------------------------------------
% 4) PARAMETER ESTIMATION
%% --------------------------------------------------------
cost = @(theta) iir_nonlinear_cost(theta,u,y,Nh,Nb,Na);

opts = optimoptions('lsqnonlin',...
    'Algorithm','levenberg-marquardt',...
    'Display','iter');

theta = lsqnonlin(cost,theta0,[],[],opts);

%% --------------------------------------------------------
% 5) PARAMETER EXTRACTION
%% --------------------------------------------------------
h = theta(1:Nh);
b = theta(Nh+1:Nh+Nb);
a = [0; 0; theta(Nh+Nb+1:end)];     % map to y^0,y^1,y^2...

%% --------------------------------------------------------
% 6) MODEL SIMULATION
%% --------------------------------------------------------
y_est = simulate_iir_nonlinear(u,y,h,b,a);

%% --------------------------------------------------------
% 7) VALIDATION
%% --------------------------------------------------------
NL_model_comp(y,y_est,fs);






function y_est = simulate_iir_nonlinear(u,y,h,b,a)
    
    Nh = length(h);
    Nb = length(b);
    Na = length(a)-2;
    
    N = length(u);
    y_est = zeros(N,1);
    
    maxLag = max(Nh,Nb);
    y_est(1:maxLag) = y(1:maxLag);
    
    for t = maxLag+1:N
        
        % ---- nonlinear feedback f(y) ----
        f = 0;
        for p = 1:Na
            f = f + a(p+2)*y_est(t)^(p+1);
        end
        
        % ---- FIR path ----
        y_fir = 0;
        for l = 0:Nh-1
            f_past = 0;
            for p = 1:Na
                f_past = f_past + a(p+2)*y_est(t-l)^(p+1);
            end
            y_fir = y_fir + h(l+1)*(u(t-l) + f_past);
        end
        
        % ---- IIR path ----
        y_iir = 0;
        for l = 1:Nb
            y_iir = y_iir + b(l)*y_est(t-l);
        end
        
        y_est(t) = y_fir + y_iir;
    end
end

function err = iir_nonlinear_cost(theta,u,y,Nh,Nb,Na)

    h = theta(1:Nh);
    b = theta(Nh+1:Nh+Nb);
    a = theta(Nh+Nb+1:end);
    
    %% Nonlinear feedback
    f = zeros(size(y));
    for p = 1:Na
        f = f + a(p)*y.^(p+1);     % y^2, y^3, ...
    end
    
    %% LTI feedforward
    y_ff = filter(h,1,u + f);
    
    %% IIR feedback
    y_fb = filter([0; b],1,y);
    
    y_est = y_ff + y_fb;
    
    err = y - y_est;
end
