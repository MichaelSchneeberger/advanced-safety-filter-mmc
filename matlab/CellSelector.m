classdef CellSelector < matlab.System
    properties (Nontunable)
        ts
        f_carr
        n_cells
    end

    properties(DiscreteState)
        time
        pre_step_up
        pwm_counter
        selection_des
        selection
        delta_flux
    end

    methods (Access = protected)

        function resetImpl(obj)
            obj.time = 0;
            obj.pwm_counter = -1;
            obj.pre_step_up = 0;
            obj.selection_des = zeros(obj.n_cells, 1);
            obj.selection = zeros(obj.n_cells, 1);
            obj.delta_flux = 0;
        end

        function [g, delta_flux] = outputImpl(obj, ~, ~, ~, ~)
            g = kron(obj.selection==1, [1 0 0 1]') ...
                + kron(obj.selection==0, [1 0 1 0]') ...
                + kron(obj.selection==-1, [0 1 1 0]');
            delta_flux = obj.delta_flux;
        end

        function num = getNumInputsImpl(~)
            num = 4;
        end

        function updateImpl(obj, v_ref, v_dc, i, enable)

            if enable == 0
                obj.selection = zeros(obj.n_cells, 1);

                return
            end

            % Update flux
            % =====================

            delta_flux_der = sum(obj.selection .* v_dc) - v_ref;

            % update delta flux
            obj.delta_flux = obj.delta_flux + delta_flux_der * obj.ts;


            obj.time = obj.time + obj.ts;
            selector_des_sum = sum(obj.selection_des);
            step_up = 0;
            

            % Magnetic flux control
            % =====================

            hysteresis = 5e3*200e-6;

            if obj.delta_flux < -hysteresis % && obj.pre_step_up < 1
                step_up = 1;
                % obj.pre_step_up = 1;
                obj.pre_step_up = max(obj.pre_step_up + 1, 1);

            elseif hysteresis < obj.delta_flux % && -1 < obj.pre_step_up
                step_up = -1;
                % obj.pre_step_up = -1;
                obj.pre_step_up = min(obj.pre_step_up - 1, -1);

            else

                % PWM
                % ===
    
                uref = max(min(v_ref / sum(v_dc), 1), -1);
                curr_counter = fix(obj.time*obj.f_carr*2*obj.n_cells);
    
                % ensure that pwm steps up/down once per carrier cycle
                if obj.pwm_counter < curr_counter || 1.5 < abs(uref - selector_des_sum)
        
                    saw_tooth = 2 * rem(obj.time, 1/obj.f_carr) * obj.f_carr;
        
                    if saw_tooth > 1
                        carrier = 2 - saw_tooth;
                    else
                        carrier = saw_tooth;
                    end
        
                    uref_abc = abs(uref);
                    slot = floor(uref_abc * obj.n_cells);
                    uref_slot = uref_abc - slot / obj.n_cells;
        
                    if carrier < uref_slot * obj.n_cells
                        level = slot + 1;
                    else
                        level = slot;
                    end
        
                    level = sign(v_ref) * level;
        
                    step_up = level - selector_des_sum;
    
                end
    
                if step_up ~= 0
                    obj.pwm_counter = curr_counter;
                end

                if 0 < step_up
                    step_up = 1;
                elseif step_up < 0
                    step_up = -1;
                end

                % % Magnetic flux control
                % % =====================
                % 
                % hysteresis = 5e3*200e-6;
                % 
                % if step_up == 0
                %     if obj.delta_flux < -hysteresis
                %         step_up = 1;
                %         obj.pre_step_up = obj.pre_step_up + 1;
                % 
                %     elseif hysteresis < obj.delta_flux
                %         step_up = -1;
                %         obj.pre_step_up = obj.pre_step_up - 1;
                % 
                %     end
    
                if 0 < step_up && step_up <= obj.pre_step_up
                    obj.pre_step_up = obj.pre_step_up - step_up;
                    step_up = 0;

                elseif 0 < obj.pre_step_up && obj.pre_step_up < step_up
                    step_up = step_up - obj.pre_step_up;
                    obj.pre_step_up = 0;

                elseif obj.pre_step_up <= step_up && step_up < 0
                    obj.pre_step_up = obj.pre_step_up - step_up;
                    step_up = 0;

                elseif step_up < obj.pre_step_up && obj.pre_step_up < 0
                    step_up = step_up - obj.pre_step_up;
                    obj.pre_step_up = 0;

                end


                % if 0 < step_up && 0 < obj.pre_step_up
                %     % step_up = step_up - 1;
                %     % obj.pre_step_up = 0; %obj.pre_step_up - 1;
                % 
                % elseif step_up < 0 && obj.pre_step_up < 0
                %     % step_up = step_up + 1;
                %     % obj.pre_step_up = 0; %obj.pre_step_up + 1;
                % end
            end


            % Cell Selector
            % =============

            if 0 < step_up && selector_des_sum < obj.n_cells

                % `selection` is a vector of length 9, whose elements are
                % either taken from {-1, 0} or {0, 1}.

                % In case of {-1, 0}, step up -1 to 0.
                if selector_des_sum < 0
                    candidates = logical(obj.selection_des);

                % In case of {0, 1}, step up 0 to 1.
                else
                    candidates = logical(~obj.selection_des);
                end
                
                % while 0 < step_up
                if i >= 0
                    [~, sel_next] = min(v_dc(candidates));
                else
                    [~, sel_next] = max(v_dc(candidates));
                end
                
                true_indices = find(candidates);
                sel_next = true_indices(sel_next);

                obj.selection_des(sel_next) = obj.selection_des(sel_next) + 1;

                %     step_up = step_up - 1;
                % end
                
            elseif step_up < 0 && -obj.n_cells < selector_des_sum

                if 0 < selector_des_sum
                    candidates = logical(obj.selection_des);
                else
                    candidates = logical(~obj.selection_des);
                end
                
                % while step_up < 0
                if i >= 0
                    [~, sel_next] = max(v_dc(candidates));
                else
                    [~, sel_next] = min(v_dc(candidates));
                end

                true_indices = find(candidates);
                sel_next = true_indices(sel_next);

                obj.selection_des(sel_next) = obj.selection_des(sel_next) - 1;

                %     step_up = step_up + 1;
                % end

            end


            % Delay PWM
            % =========

            delta_flux_des_der = sum(obj.selection_des .* v_dc) - v_ref;

            % hysteresis = 0;
            hysteresis = 5e3*40e-6;

            if obj.delta_flux < -hysteresis && delta_flux_der < delta_flux_des_der
                obj.selection = obj.selection_des;

            elseif hysteresis < obj.delta_flux && delta_flux_des_der < delta_flux_der
                obj.selection = obj.selection_des;
            end
            % obj.selection = obj.selection_des;

        end

    end

end