classdef SafetyFilterSystemDPC < matlab.System
   properties
        ts
        f_data
        G_data
        V_data
        B_data
        nV_data
        nB_data
        gamma_V_data
        gamma_B_data
        r
        l
        init_v_pcc_f = [1; 0];
    end

    % Pre-computed constants or internal states
    properties (Access = private)
        init_dm = [0; 0];
        % init_v_pcc_f = [1; 0];
        init_dm_corr = [0; 0];
        init_dv_corr = [0; 0];
        init_active = 0;
        init_lambda = [0; 0];
    end

    properties(DiscreteState)
        dm
        v_pcc_f
        dm_corr
        dv_corr
        active
        lambda
    end

    methods (Access = protected)

        function setupImpl(obj)
            obj.dm = obj.init_dm;
            obj.v_pcc_f = obj.init_v_pcc_f;
            obj.dm_corr = obj.init_dm_corr;
            obj.dv_corr = obj.init_dv_corr;
            obj.active = obj.init_active;
            obj.lambda = obj.init_lambda;
        end

        function resetImpl(obj)
            obj.dm = obj.init_dm;
            obj.v_pcc_f = obj.init_v_pcc_f;
            obj.dm_corr = obj.init_dm_corr;
            obj.dv_corr = obj.init_dv_corr;
            obj.active = obj.init_active;
            obj.lambda = obj.init_lambda;
        end

        function [dm, v_pcc_f, dm_corr, dv_corr, active, lambda] = outputImpl(obj, ~)
            dm = obj.dm;
            v_pcc_f = obj.v_pcc_f;
            dm_corr = obj.dm_corr;
            dv_corr = obj.dv_corr;
            active = obj.active;
            lambda = obj.lambda;
        end

        function updateImpl(obj, u)

            % read input arguments
            dm_ref = u(1:2);
            i_ref = u(3:4);
            i = u(5:6);
            v_pcc = u(7:8);
            enable = u(9);

            if ~enable

                obj.active = 0;

                obj.dm = dm_ref;
                obj.v_pcc_f = obj.init_v_pcc_f;

                obj.dm_corr = [0; 0];
                obj.dv_corr = [0; 0];
                obj.lambda = [0; 0];

                return

            end
                
            dv_pcc_f = obj.v_pcc_f - v_pcc;
            x = [i; dv_pcc_f; i_ref];

            xx = reshape(x * x', [], 1);

            x01 = [1;x];
            x02 = [1;xx];
            x012 = [1;x;xx];

            V = obj.V_data * x012;
            B = obj.B_data * x012;

            obj.active = 1;

            u_ref = [
                dm_ref + dv_pcc_f;
                -100 * dv_pcc_f
                % -10 * dv_pcc_f
            ];

            f_sys = obj.f_data * x;
            G_sys = reshape(obj.G_data, 6, 4);

            nV = obj.nV_data * x01;
            nB = obj.nB_data * x01;

            gamma_V = obj.gamma_V_data;
            gamma_B = obj.gamma_B_data;

            % cost = (1/2) * x' * H * x + f' * x
            H = eye(4);
            f = -u_ref;

            % A x <= b
            A = [ ...
                   nB' * G_sys;
                   nV' * G_sys ...
                ];       
            b = [ ...
                   -nB' * f_sys - gamma_B * B;
                   (-nV' * f_sys - gamma_V * V) ...
                 ];

            lb = [-2; -2; -500; -500];
            ub = [2; 2; 500; 500];

            x0 = u_ref;

            options = optimoptions('quadprog', 'Algorithm', 'active-set');
            [u, ~, ~, ~, lambda_val] = quadprog(H, f, A, b, [], [], lb, ub, x0, options);
    
            % next_dv_pcc_f = dv_pcc_f + obj.ts * u(3:4, 1);
            next_dv_pcc_f = dv_pcc_f + obj.ts * u_ref(3:4, 1);

            % update states
            obj.v_pcc_f = v_pcc + next_dv_pcc_f;
            % obj.dm = u(1:2, 1) - next_dv_pcc_f;
            obj.dm = u(1:2, 1) - dv_pcc_f;

            % debugging
            obj.dm_corr = u(1:2, 1) - u_ref(1:2, 1);
            obj.dv_corr = u(3:4, 1) - u_ref(3:4, 1);
            obj.lambda = lambda_val.ineqlin(1:2,1);

        end

        % function flag = isInputDirectFeedthroughImpl(~)
        %     flag = true;
        % end

    end
end
