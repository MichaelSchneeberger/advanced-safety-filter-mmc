function output(out, fileName)

sim_time = out.simout.i_m.Time;
v_m_abc = out.simout.v_m_abc.Data;
i_m = out.simout.i_m.Data;
i_m_abc = out.simout.i_m_abc.Data;
i_abs = sqrt(i_m(:,1).^2 + i_m(:,2).^2);
dm_corr = out.simout.dm_corr.Data;
w_m = out.simout.w_m.Data;

% tStart = 0.7;
tStart = 1.2;
tEnd = 3;

sel = tStart <= sim_time & sim_time <= tEnd;
% sel(~mod(1:length(sel), 10) == 0) = 0;

time = sim_time(sel);

fig = figure(2);
% set(fig, 'Position', [0, 0, 1600, 1000])
% set(fig, 'Position', [0, 0, 500, 800])
set(fig, 'Position', [0, 0, 300, 500])

subplot(5,1,1)
plot(time, v_m_abc(sel,:), 'LineWidth', 0.5)
% legend({'a', 'b', 'c'})
ylabel('v [p.u.]')
xlim([tStart tEnd])
ylim([-1.3, 1.3])

subplot(5,1,2)
plot(time, i_m_abc(sel,:), 'LineWidth', 0.5)
hold on;
plot(time, 1.3*ones(size(time)), '--r', 'LineWidth', 0.5)
plot(time, -1.3*ones(size(time)), '--r', 'LineWidth', 0.5)
hold off;
% legend({'a', 'b', 'c'})
ylabel('i_{abc} [p.u.]')
xlim([tStart tEnd])
ylim([-1.8, 1.8])

subplot(5,1,3)
plot(time, i_m(sel,1:2), 'LineWidth', 0.5)
% hold on;
% plot(time, i_abs(sel,:), 'LineWidth', 1.5)
% plot(time, 1.3*ones(size(time)), '--r', 'LineWidth', 1.5)
% hold off;
ylabel('i_d, i_q [p.u.]')
% legend({'d component', 'q component', '|| . ||'})
% ylabel('i, || i || [p.u.]')
xlim([tStart tEnd])
ylim([-1.5, 1.5])

subplot(5,1,4)
plot(time, w_m(sel,:), 'LineWidth', 0.5)
ylabel('\omega_{PLL} [p.u.]')
xlim([tStart tEnd])
% ylim([-0.1, 1.34])

subplot(5,1,5)
plot(time, dm_corr(sel,:), 'LineWidth', 0.5)
% plot(time, out_curr.simout.lambda.Data(sel,:))
% legend({'d', 'q'})
ylabel('\Deltau [p.u.]')
xlabel('time [s]')
xlim([tStart tEnd])
% ylim([-0.2, 0.15])

% saveas(fig, 'sim.jpg')
exportgraphics(fig, strcat(fileName, '.pdf'),'ContentType','vector')

save(strcat(fileName, '.mat'), 'out')
% matFileName = strcat(fileName, '.mat');
% save matFileName out

end