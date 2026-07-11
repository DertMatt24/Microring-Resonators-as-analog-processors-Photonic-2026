close all
clear all 
clc

k_eq = 62.5; % ns^-1
k = 62.5;  % ns^-1
A=1e9;  %time scaling parameter  [nano 10^9]
A_time = A * 1e3;

% ODE SOLVER: k = Bs -> A_time = A * 1e3;
% FWHM = 19.07 * 1e-3;
FWHM = 19.07 * 1e-3; 

% Input signal x(t)
%x = Model_utils.sin_function(A);
x = Model_utils.super_gaussian_function(FWHM, A);
%x = Model_utils.arbitrary_signal(A);

% Definition of the ODE
odefun = Model_utils.first_order_ode(A, k_eq, x);

% Initial condition
y0 = 0;

% Time span
N=1e5;
interval = 100/A_time;
tspan = linspace(-interval , interval, N);

% Solve using ode45
[t, y] = ode45(odefun, tspan, y0);


%% creatinng MRR
R=30e-6;  %radius of the MRR
neff=2.4;   %effective index of the MRR waveguide
MRR = mrr(R, neff, k, A, 0);

%% Computing input

time=linspace(min(t),max(t),N);
dt=time(2)-time(1);
in_ring = x(time);
IN_ring=fftshift(fft(in_ring));

Df=linspace(-1/(2*dt),1/(2*dt),N);

fsr = MRR.FSR(neff); % no chromatic dispersion assumption
b3db = MRR.B3dB(fsr);
N_det = 1e3;
delta_f_array = linspace(-2 * fsr, 2 * fsr, N_det);

% helpful vectors
phase_detuning = zeros(1, N_det);

% output vectors
P_out  = zeros(1, N_det);
P_in = zeros(1, N_det);
P_lost_db = zeros(1, N_det);

rmse_det = zeros(1, N_det);

for i = 1:N_det
    delta_f = delta_f_array(i);

    [H_drop, H_drop_norm] = MRR.h_drop_f(Df, delta_f);

    Out_ring = IN_ring .* H_drop;
    out_ring = real(ifft(ifftshift(Out_ring)));
    
    Out_ring_norm = IN_ring .* H_drop_norm;
    out_ring_norm = real(ifft(ifftshift(Out_ring_norm)));

    P_in(i) = MRR.power(in_ring, dt);
    P_out(i) = MRR.power(out_ring, dt);
    
    [~, P_lost_db(i)] = MRR.power_loss(in_ring, out_ring, dt);

    rmse_det(i) = Model_utils.computing_rmse(time, t, out_ring_norm, y, N);
    
end

delta_f_array = delta_f_array / fsr;


% Plot grph to see power loss depending from phase detuning
figure(401)
graph_drawer.power_loss(delta_f_array, P_lost_db);

figure(402)
graph_drawer.plot_rmse(delta_f_array, rmse_det);