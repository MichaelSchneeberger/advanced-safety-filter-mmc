% 0 - MMC cell
% 1 - averaged model
ctrl.cell.sel = 0;
% ctrl.cell.sel = 1;

% 0: infinite bus
% 1: synchronous machine
% 2: GFL converter v_DC=const, 
% 3: passive load
% 4: GFL converter
for grid_sel = [1 4]
% for grid_sel = [4]
    machine.grid = grid_sel;
    
    if grid_sel == 1
        grid_name = 'sm';
        testScenario.p = 0.9;
    elseif grid_sel == 2
        grid_name = 'gfl';
        testScenario.p = -0.9;
    elseif grid_sel == 4
        grid_name = 'gfl';
        testScenario.p = -0.9;
    end

    for gfm_sel = [1 2]
    % for gfm_sel = [2]
        ctrl.machine.gfm_sel = gfm_sel;

        if gfm_sel == 1
            gfm_name = 'edpc';
        elseif gfm_sel == 2
            gfm_name = 'vsm';
        end

        for lim_sel = [1 3 4 7 9]
        % for lim_sel = [7 9]
            % lim_sel = 7;
            ctrl.machine.ictrl.sel = lim_sel;
    
            if lim_sel == 1
                lim_name = 'scc';
            elseif lim_sel == 3
                lim_name = 'ai';
            elseif lim_sel == 4
                lim_name = 'cl';
            elseif lim_sel == 7
                lim_name = 'sf';
            elseif lim_sel == 9
                lim_name = 'sfB';
            end

            testName = strcat(grid_name, '_', gfm_name, '_', lim_name);
            disp(testName)

            modelName = 'mmc_model';
    
            set_param(modelName,SimulationCommand='Update')
            simResult = sim(modelName);

            sim_time = simResult.simout.i_m.Time;
            sel = 1.2 <= sim_time;

            i = simResult.simout.i_m_abc.Data(sel,:);
            v = simResult.simout.v_m_abc.Data(sel,:);
            dm_corr = simResult.simout.dm_corr.Data(sel,:);

            dm_corr_avg = sum(sqrt(dm_corr(:,1).^2 + dm_corr(:,2).^2)) / length(sel);
            dm_corr_avg
            max(max(abs(i)))
            max([thd(i(:, 1)) thd(i(:, 2)) thd(i(:, 3))])
            max([thd(v(:, 1)) thd(v(:, 2)) thd(v(:, 3))])

            output(simResult, testName)

            % load(strcat(testName, '.mat'))
            % output(out, testName)
            % % return;

        end

    end

end