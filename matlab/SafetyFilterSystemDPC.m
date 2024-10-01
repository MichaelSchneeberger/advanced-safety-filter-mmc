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
        init_vpcc_f
    end

    % Pre-computed constants or internal states
    properties (Access = private)
        init_m = [0; 0];
        init_dv_pcc_f = [0; 0];
        init_dm_corr = [0; 0];
        init_active = 0;
        init_lambda = [0; 0];
    end

    properties(DiscreteState)
        m
        dv_pcc_f
        dm_corr
        active
        lambda
    end

    methods (Access = protected)

        function setupImpl(obj)
            obj.m = obj.init_m;
            obj.dv_pcc_f = obj.init_dv_pcc_f;
            obj.dm_corr = obj.init_dm_corr;
            obj.active = obj.init_active;
            obj.lambda = obj.init_lambda;
        end

        function resetImpl(obj)
            obj.m = obj.init_m;
            obj.dv_pcc_f = obj.init_dv_pcc_f;
            obj.dm_corr = obj.init_dm_corr;
            obj.active = obj.init_active;
            obj.lambda = obj.init_lambda;
        end

        function [m, v_pcc_f, dm_corr, active, lambda] = outputImpl(obj, ~)
            m = obj.m;
            v_pcc_f = obj.dv_pcc_f + [1; 0];
            dm_corr = obj.dm_corr;
            active = obj.active;
            lambda = obj.lambda;
        end

        function updateImpl(obj, u)

            % read arguments
            ff = u(1:2);
            i_ref = u(3:4);
            i = u(5:6);
            v_pcc = u(7:8);
            w = u(9);
            enable = u(10);

            dv_pcc = v_pcc - [1; 0];

            % set state default values
            obj.dm_corr = [0; 0];
            obj.lambda = [0; 0];

            if ~enable

                obj.active = 0;

                obj.m = ff + obj.init_vpcc_f;
                obj.dv_pcc_f = obj.init_vpcc_f - [1; 0];

                return

            end
                
            x = [i; obj.dv_pcc_f; i_ref; dv_pcc];

            xx = reshape(x * x', [], 1);

            x01 = [1;x];
            x02 = [1;xx];
            x012 = [1;x;xx];

            V = obj.V_data * x012;
            B = obj.B_data * x012;

            V2 = 0.1 * V + 0.9 * B;

            J = [0 -1; 1 0];
            cross_term = (obj.r + J*w*obj.l) * i;

            obj.active = 1;

            u_ref = [
                % desired voltage applied to the converter (without
                % cross-terms) relative to v_c = [1; 0]
                ff - cross_term + obj.dv_pcc_f;

                100 * (dv_pcc - obj.dv_pcc_f)
            ];

            f_sys = obj.f_data * x;
            G_sys = reshape(obj.G_data, 8, 4);

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
                   -nV' * f_sys - gamma_V * V ...
                 ];

            lb = [-2.3; -1.3; -500; -500];
            ub = [0.3; 1.3; 500; 500];

            x0 = u_ref;

            options = optimoptions('quadprog', 'Algorithm', 'active-set');
            [u, ~, ~, ~, lambda_val] = quadprog(H, f, A, b, [], [], lb, ub, x0, options);
    
            obj.lambda = lambda_val.ineqlin(1:2,1);
            obj.dm_corr = u(1:2, 1) - u_ref(1:2, 1);

            obj.m = u(1:2, 1) + cross_term + [1; 0];

            obj.dv_pcc_f = obj.dv_pcc_f + obj.ts * u(3:4, 1);

        end

        function flag = isInputDirectFeedthroughImpl(~)
            flag = false;
        end

    end
end
