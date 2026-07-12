classdef mrr_cascade < handle
    % MRR_CASCADE  N cascaded microring resonators solving an N-th order ODE.
    
    properties
        rings    % cell array of mrr objects, one per stage
        order    % number of cascaded rings N
        A        % time scale factor (shared by all rings)
    end

    methods
        % CONSTRUCTOR
        %   R, neff : ring geometry (shared by every stage)
        %   k_vec   : 1xN vector of coupling coefficients (one per ring).
        %             A scalar is broadcast to all N rings.
        %   A       : time scaling parameter
        %   loss    : loss factor (shared)
        %   order   : number of rings N
        function obj = mrr_cascade(R, neff, k_vec, A, loss, order)
            if nargin ~= 6
                error("mrr_cascade needs 6 args: R, neff, k_vec, A, loss, order");
            end
            if isscalar(k_vec)
                k_vec = k_vec * ones(1, order);
            end
            if numel(k_vec) ~= order
                error("k_vec must be scalar or have length = order (%d)", order);
            end

            obj.order = order;
            obj.A = A;
            obj.rings = cell(1, order);
            for i = 1:order
                obj.rings{i} = mrr(R, neff, k_vec(i), A, loss);
            end
        end

        % Retune ring i to a new k
        function tune_ring(obj, i, k_new)
            obj.rings{i}.tuning_k(k_new);
        end

        % Retune every ring to the same k
        function tune_all(obj, k_new)
            for i = 1:obj.order
                obj.rings{i}.tuning_k(k_new);
            end
        end

        %% Total drop-port transfer function
        %   delta_f : scalar (same detuning on all rings) OR 1xN vector
        %             (per-ring detuning, for the phase-detuning study).
        %   Returns H (raw) and H_norm (scaled by prod(1/k_i), DC gain ~ 1).
        function [H, H_norm] = h_drop_cascade(obj, Df, delta_f)
            if nargin < 3, delta_f = 0; end
            if isscalar(delta_f), delta_f = delta_f * ones(1, obj.order); end

            H = ones(size(Df));
            for i = 1:obj.order
                [h_i, ~] = obj.rings{i}.h_drop_f(Df, delta_f(i));
                H = H .* h_i;
            end
            H_norm = H;
            for i = 1:obj.order
                k_i = obj.rings{i}.k_ring / obj.A;
                H_norm = H_norm ./ k_i;   % remove per-stage DC gain 1/k_i
            end
        end

        %% Total ideal ODE transfer function (pure Lorentzian cascade)
        function H = h_ode_cascade(obj, Df, delta_f)
            if nargin < 3, delta_f = 0; end
            if isscalar(delta_f), delta_f = delta_f * ones(1, obj.order); end

            H = ones(size(Df));
            for i = 1:obj.order
                H = H .* obj.rings{i}.h_ode(Df, delta_f(i));
            end
        end

        %% Power helpers (delegate to a ring; identical geometry)
        function p = power(obj, sig, dt)
            p = obj.rings{1}.power(sig, dt);
        end

        function [P_loss, P_loss_dB] = power_loss(obj, x, y, dt)
            [P_loss, P_loss_dB] = obj.rings{1}.power_loss(x, y, dt);
        end
    end
end
