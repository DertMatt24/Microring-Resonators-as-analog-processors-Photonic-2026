close all
clear all
clc

% =====================================================================
%  Cascaded N-th order ODE solver  (drop-port MRR cascade)
% =====================================================================

%% ---- Ring parameters (shared geometry) --------------------------------
A     = 1e9;      % time scaling parameter [ns -> s]
R     = 30e-6;    % ring radius [m]
neff  = 2.4;      % effective index
loss  = 0;        % loss factor
A_time = A * 1e3; % ODE-solver time regime (k ~ signal bandwidth)

max_order = 4;

k_base  = 62.5;
k_vec   = linspace(k_base*0.85, k_base*1.15, max_order);

fprintf('Ring parameters:\n');
fprintf('  R = %.0f um, neff = %.2f, loss = %g\n', R*1e6, neff, loss);
fprintf('  k_i [ns^-1] = %s\n', mat2str(round(k_vec*10)/10));
fprintf('  all within Yang tunable range [%.0f, %.0f]\n\n', ...
    min(mrr.k_Yang), max(mrr.k_Yang));

inputs = {'gaussian', 'super-gaussian'};

orders   = 1:max_order;
rmse_all = zeros(2, max_order);
loss_ins = zeros(2, max_order);
loss_shp = zeros(2, max_order);

N        = 1e5;
% tspan    = linspace(-interval, interval, N);
t_min = -50; t_max = 150;                 

for s = 1:2
    is_sg = strcmp(inputs{s}, 'super-gaussian');

    % Input signal x(t): Gaussian (m=1) or super-Gaussian (m=3)
    if is_sg
        FWHM = 41.54e-3;   % ps-scale FWHM for super-Gaussian
        x = Model_utils.super_gaussian_function(FWHM, A, 3);
    else
        FWHM = 19.07e-3;   % ps-scale FWHM for Gaussian
        x = Model_utils.super_gaussian_function(FWHM, A, 1);
    end

    ring0 = mrr(R, neff, k_base, A, loss);
    fsr   = ring0.FSR(neff);

    n_fsr_range = 1;
    dt   = 1/(2 * n_fsr_range * fsr);   
    T    = 800e-12;                     
    N    = round(T/dt);
    if mod(N,2) == 1, N = N + 1; end    

    time = (-(N/2):(N/2-1)) * dt;              
    Df   = (-(N/2):(N/2-1)) / (N*dt);          

    p = 1;
    for i = 1:max_order, p = conv(p, [1, k_vec(i)]); end
    a = p(2:end);
    if max_order == 1
        comp = -a;
    else
        comp = [zeros(max_order-1,1), eye(max_order-1); -fliplr(a)];
    end
    odefun_n = @(tt, z) A * (comp*z + [zeros(max_order-1,1); x(tt)]);
    [t_exact, Z] = ode45(odefun_n, time, zeros(max_order,1));
    y_exact = Z(:,1);

    in_ring = x(time);
    IN_ring = fftshift(fft(in_ring));


    % Collect transfer functions for the spectrum-by-order plot
    H_by_ord = cell(1, max_order);
    H_norm_by_ord = cell(1, max_order);

    for n = orders
        casc = mrr_cascade(R, neff, k_vec(1:n), A, loss, n);

        [H, H_norm] = casc.h_drop_cascade(Df, 0);   % real MRR cascade
        H_ideal     = casc.h_ode_cascade(Df, 0);    % ideal Lorentzian cascade
        H_by_ord{n} = H;
        H_norm_by_ord{n} = H_norm;


        out_real  = real(ifft(ifftshift(IN_ring .* H_norm)));
        out_raw   = real(ifft(ifftshift(IN_ring .* H)));
        out_ideal = real(ifft(ifftshift(IN_ring .* H_ideal)));

        % RMSE of shape (peak-normalized) vs ideal ODE output
        oN  = out_real  / max(abs(out_real));
        oiN = out_ideal / max(abs(out_ideal));
        rmse_all(s, n) = sqrt(mean(abs(oN - oiN).^2));

        % Raw insertion loss
        [~, loss_ins(s, n)] = casc.power_loss(in_ring, out_raw, dt);

        % Shape loss vs ideal ODE
        Preal  = casc.power(out_real,dt);
        Pideal = casc.power(out_ideal,dt);
        loss_shp(s, n) = 10*log10(Preal / Pideal);


    end

    out_ring = real(ifft(ifftshift(IN_ring .* H_norm_by_ord{max_order})));


    figure('Name', ['input vs output - ' inputs{s}]);
    utils = graph_drawer(t_min, t_max, time);
    utils.input_output_power(in_ring, out_ring, t_exact, y_exact, A_time);
    subplot(311); title(sprintf('%s | order %d - input x(t)', inputs{s}, max_order));

    % ---- Spectrum by order----
    figure('Name', ['Spectra by order - ' inputs{s}]);
    casc_top = mrr_cascade(R, neff, k_vec, A, loss, max_order);
    H_ideal_top = casc_top.h_ode_cascade(Df, 0);
    graph_drawer.spectrum_orders(Df, H_by_ord, H_ideal_top, IN_ring, fsr, 1);
    title(sprintf('Drop-port spectra vs order (%s)', inputs{s}));


    % ---------------- PHASE-DETUNING POWER ANALYSIS ----------------
    det_order = max_order;
    casc_d = mrr_cascade(R, neff, k_vec, A, loss, det_order);

    n_det   = 4001;
    delta_f_array = linspace(-2 * fsr, 2 * fsr, n_det);

    P_lost_db = zeros(1, n_det);
    rmse_det = zeros(1, n_det);
    [~, Hn0] = casc_d.h_drop_cascade(Df, zeros(1, det_order));
    out0 = ifft(ifftshift(IN_ring .* Hn0));
    casc_d.power(out0,dt);
    ref_peak = max(abs(out0));

    for i = 1:n_det
        delta_f = delta_f_array(i);
        [H_drop, H_drop_norm] = casc_d.h_drop_cascade(Df, delta_f * ones(1, det_order));

        Out_ring = IN_ring .* H_drop;
        out_ring = abs(ifft(ifftshift(Out_ring)));

        Out_ring_norm = IN_ring .* H_drop_norm;
        out_ring_norm = abs(ifft(ifftshift(Out_ring_norm)));

        [~, P_lost_db(i)] = casc_d.power_loss(in_ring, out_ring, dt);

        out_n = out_ring_norm(:) / ref_peak;
        y_n   = y_exact(:) / max(abs(y_exact));
        rmse_det(i) = Model_utils.computing_rmse(time, t_exact, out_n, y_n, N);

    end

    delta_f_array = delta_f_array / fsr;


    % Plot grph to see power loss depending from phase detuning
    figure("Name",['power loss vs phase detnuning' inputs{s}])
    graph_drawer.power_loss(delta_f_array, P_lost_db);

    figure("Name",['rmse vs phase detnuning' inputs{s}])
    graph_drawer.plot_rmse(delta_f_array, rmse_det);

    fprintf('\n');

end

%% ---- Power discussion over order --------------------------------------
graph_drawer.plot_power_vs_order(orders, loss_ins(1,:), loss_ins(2,:));
title('Insertion loss vs order');

graph_drawer.plot_power_vs_order(orders, loss_shp(1,:), loss_shp(2,:));
title('Shape loss vs order');

%% ---- RMSE vs order, both inputs ---------------------------------------
graph_drawer.plot_rmse_vs_order(orders, rmse_all(1,:), rmse_all(2,:));

%% ---- Summary table -----------------------------------------------------
fprintf('\n=== Summary ===\n');
fprintf('Order |  RMSE(G)  RMSE(SG) | InsLoss(G) ShapeLoss(G) | InsLoss(SG) ShapeLoss(SG)\n');
for n = orders
    fprintf('  %d   |  %.4f   %.4f  |  %6.1f dB    %5.2f dB|  %6.1f dB    %5.2f dB\n', ...
        n, rmse_all(1,n), rmse_all(2,n), loss_ins(1,n), loss_shp(1,n),loss_ins(2,n), loss_shp(2,n));
end
