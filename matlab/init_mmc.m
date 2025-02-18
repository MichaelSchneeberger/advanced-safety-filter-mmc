clear

J = [0 -1; 1 0];


% General parameters
%%%%%%%%%%%%%%%%%%%%%


simulation.ts = 1e-6;
simulation.tsim = 2;
simulation.tsim = 3;

% cpu.ts = 100e-6;
cpu.ts = 200e-6;

% nominal converter power
converter.s = 40e6;


% Grid side parameters
%%%%%%%%%%%%%%%%%%%%%%

grid.fn = 50;
% grid.v = 220e3;
grid.v = 33e3;
grid.tr.l = 0.1;  % in p.u.
grid.tr.r = 0.01;  % in p.u.

grid.breaker_close = 0.04;


% Machine side parameters
%%%%%%%%%%%%%%%%%%%%%%%%%

machine.fn = 60;
machine.v = 33e3;
machine.tr.l = 0.16;  % in p.u.
machine.tr.r = 0.02;  % in p.u.
machine.filter.c = 0.006;
machine.filter.r = 10;
machine.filter.l = 0.2;

machine.breaker_close = 0.04;

% 0: infinite bus
% 1: SM
% 2: GFL converter (no DC-link), 
% 3: passive load
% 4: GFL converter
machine.grid = 0;
machine.grid = 1;
% machine.grid = 2;
% machine.grid = 3;
machine.grid = 4;


% MMC parameters
%%%%%%%%%%%%%%%%

mmc.l_leg_si = 5e-3;
mmc.r_leg_si = 10e-3;

mmc.cell.v_dc = 5e3;
mmc.cell.c_dc = 5e-3;
mmc.cell.n = 14;

% 0 - MMC cell
% 1 - averaged model
ctrl.cell.sel = 0;
ctrl.cell.sel = 1;


% MMC Control
%%%%%%%%%%%%%

ctrl.grid.pll.kp = 0.096;
ctrl.grid.pll.ti = 0.0849;
% ctrl.grid.pll.init_w = 1;
ctrl.grid.pll.w0 = 1;
ctrl.grid.pll.enable = 0.01;

ctrl.grid.ictrl.kp = 0.342;
ctrl.grid.ictrl.ti = 2e-3;
ctrl.grid.ictrl.enable = 0.02;
ctrl.machine.ictrl.ref_enable = 0.1;

% 0 - EDPC (modified)
% 1 - EDPC
% 2 - VSM
% 3 - dVOC
ctrl.machine.gfm_sel = 0;
ctrl.machine.gfm_sel = 1;
% ctrl.machine.gfm_sel = 2;
% ctrl.machine.gfm_sel = 3;
% ctrl.machine.gfm_sel = 4;
% ctrl.machine.gfm_sel = 5;

% 0 - open loop
% 1 - vector current control
% 2 - safety filter
% 3 - adaptive VI
% 4 - Current Limiting Control
% 5 - safety filter (CBF only)
% 6 - safety filter (no vPCC filter)
% 7 - safety filter (i0) (2 states)
% 8 - safety filter (i0)
% 9 - safety filter (i0) (2 states, B only)
ctrl.machine.ictrl.sel = 0;
% ctrl.machine.ictrl.sel = 1;
% ctrl.machine.ictrl.sel = 2;
ctrl.machine.ictrl.sel = 3;
% ctrl.machine.ictrl.sel = 4;
% ctrl.machine.ictrl.sel = 5;
% ctrl.machine.ictrl.sel = 6;
ctrl.machine.ictrl.sel = 7;
% ctrl.machine.ictrl.sel = 8;
% ctrl.machine.ictrl.sel = 9;

ctrl.vdcctrl.kp = 0.342*20;
ctrl.vdcctrl.ti = 0.1;

pq_ctrl.f_ref = 1.0;
pq_ctrl.vg_ref = 1.0;
pq_ctrl.fp_droop = 0.02;
pq_ctrl.uq_droop = 0.05;
pq_ctrl.kp = 0.45;
pq_ctrl.ti = 0.04*3;
pq_ctrl.enable = 0.4;

% soa.i_lim = 1.18;
soa.i_lim = 1;
soa.m_max = 1.23;

% Switched Current Vector Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% enable current control if i > i_lim
currentControl.i_lim = 1.24;
% currentControl.i_lim = 0;


% Current Limitation Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc.k_pcc = 0.1;
clc.k_pcc = 4;
clc.k_pcc = ctrl.grid.ictrl.kp;
clc.i_lim = 1.18;

% Synchronization voltage
%%%%%%%%%%%%%%%%%%%%%%%%%

% sync.breaker_open = 100;
sync.breaker_open = 0.2;

if ctrl.machine.gfm_sel == 2
    sync.breaker_open = 0.2;
end

% Adapted VI
%%%%%%%%%%%%

avi.i_lim = 1.18;
avi.k_x = 10;
% avi.k_x = 40;
% avi.k_x = 100;

% Test Scenario
%%%%%%%%%%%%%%%

% active power consumed by the grid
if machine.grid == 0
    testScenario.p = -1;

% synchronous machine
elseif machine.grid == 1
    % testScenario.p = -0.9;
    testScenario.p = 0.9;
    % testScenario.p = 0.5;
    % testScenario.p = 0;

elseif machine.grid == 4
    testScenario.p = -0.9;
    % testScenario.p = -0.5;
end


% Stiff Grid
%%%%%%%%%%%%

stiffGrid.l = 0.16;
stiffGrid.r = 0.02;
stiffGrid.v = 1;
% stiffGrid.v = 1.1;


% Grid Following Converter
%%%%%%%%%%%%%%%%%%%%%%%%%%

gfl.vdc.kp = 6;
gfl.vdc.ti = 0.1;
gfl.tr.l = 0.16;  % in p.u.
gfl.tr.r = 0.01;  % in p.u.
gfl.i_ref = [testScenario.p; 0];
gfl.pll_enable = 0.01;
gfl.breaker_close = 0.04;
gfl.ctrl_enable = gfl.breaker_close;
gfl.ref_step = 0.3;
gfl.ref_step = 0.3;
gfl.c_dc = 1.5;


% Load Step
%%%%%%%%%%%

% ratio = 0;
% ratio = 0.01;
ratio = 0.1;
% ratio = 1;
% ratio = 2;
% ratio = 5;
% ratio = 10;
% if machine.grid == 4
%     ratio = 5;
%     % ratio = 10;
% end
loadStep.r = machine.tr.r * ratio;
loadStep.l = machine.tr.l * ratio;
% loadStep.stepDown = 1 - 2.5*simulation.ts;
loadStep.stepDown = 1.4 - 2.5*simulation.ts;
% loadStep.stepDown = 4 - 2.5*simulation.ts;
loadStep.stepUp = loadStep.stepDown + 0.3;
% loadStep.stepUp = loadStep.stepDown + 0.1;


dVOC.kappa = atan(machine.tr.l/machine.tr.r);
dVOC.Rkappa = [cos(dVOC.kappa) -sin(dVOC.kappa); sin(dVOC.kappa) cos(dVOC.kappa)];


%%

grid.wn = 2*pi*grid.fn;
grid.vn = grid.v / sqrt(3);
grid.in = converter.s / 3 / grid.vn;
grid.zn = grid.vn / grid.in;
grid.ln = grid.zn / grid.wn;
grid.vpeak = grid.vn * sqrt(2);
grid.ipeak = grid.in * sqrt(2);
grid.tr.l_si = grid.tr.l * grid.ln;
grid.tr.r_si = grid.tr.r * grid.zn;

machine.f_ratio_50 = machine.fn / 50;
machine.wn = 2*pi*machine.fn;
% machine.w0 = machine.wn * machine.w0_pu;
machine.vn = machine.v / sqrt(3);
machine.sn = converter.s / 3;
machine.in = machine.sn / machine.vn;
machine.zn = machine.vn / machine.in;
machine.ln = machine.zn / machine.wn;
machine.cn = 1 / machine.zn / machine.wn;
machine.vpeak = machine.vn * sqrt(2);
machine.ipeak = machine.in * sqrt(2);
machine.tr.l_si = machine.tr.l * machine.ln;
machine.tr.r_si = machine.tr.r * machine.zn;
machine.filter.c_si = machine.filter.c * machine.cn;
machine.filter.r_si = machine.filter.r * machine.zn;
machine.filter.l_si = machine.filter.l * machine.ln;

% mmc.cell.c_n = converter.s / 3 / mmc.cell.v_dc^2;
mmc.cell.c_n = converter.s / 9 / mmc.cell.v_dc^2;
mmc.cell.c_dc_pu = mmc.cell.c_dc / mmc.cell.c_n;
mmc.cell.total_v_dc = mmc.cell.v_dc * mmc.cell.n;
mmc.l_leg_g = mmc.l_leg_si / grid.ln;
mmc.l_leg_m = mmc.l_leg_si / machine.ln;
mmc.r_leg_g = mmc.r_leg_si / grid.zn;
mmc.r_leg_m = mmc.r_leg_si / machine.zn;

ctrl.grid.pll.ki = ctrl.grid.pll.kp / ctrl.grid.pll.ti;
ctrl.grid.ictrl.ki = ctrl.grid.ictrl.kp / ctrl.grid.ictrl.ti;

ctrl.machine.pll.kp = ctrl.grid.pll.kp;
ctrl.machine.pll.ti = ctrl.grid.pll.ti / machine.f_ratio_50 * 10.8;
ctrl.machine.pll.ki = ctrl.machine.pll.kp / ctrl.machine.pll.ti;
ctrl.machine.pll.enable = ctrl.grid.pll.enable;
ctrl.machine.pll.w0 = 1;
ctrl.machine.ictrl.kp = ctrl.grid.ictrl.kp;
ctrl.machine.ictrl.ti = ctrl.grid.ictrl.ti / machine.f_ratio_50;
ctrl.machine.ictrl.enable = ctrl.grid.ictrl.enable;
ctrl.machine.ictrl.ki = ctrl.machine.ictrl.kp / ctrl.machine.ictrl.ti;

ctrl.vdcctrl.ki = ctrl.vdcctrl.kp / ctrl.vdcctrl.ti;

pq_ctrl.ki = pq_ctrl.kp / pq_ctrl.ti;


% Stiff Grid
%%%%%%%%%%%%

stiffGrid.l_si = stiffGrid.l * machine.ln;
stiffGrid.r_si = stiffGrid.r * machine.zn;

% stiff grid case
if machine.grid == 0
    stiffGrid.w_pu = 1 + testScenario.p * pq_ctrl.fp_droop;

else
    stiffGrid.w_pu = 1;
end

stiffGrid.w = stiffGrid.w_pu * machine.wn;
% stiffGrid.v = 1 - testScenario.q * machine.tr.l;
stiffGrid.vpeak = machine.vpeak * stiffGrid.v;


% Grid Following Converter
%%%%%%%%%%%%%%%%%%%%%%%%%%

gfl.tr.l_si = gfl.tr.l * machine.ln;
gfl.tr.r_si = gfl.tr.r * machine.zn;
gfl.v_dc = 2*machine.vpeak;
gfl.i_dc_n = converter.s / gfl.v_dc;
gfl.z_dc_n = gfl.v_dc^2 / converter.s;
gfl.c_dc_n = 1 / (machine.wn * gfl.z_dc_n);
gfl.c_dc_si = gfl.c_dc * gfl.c_dc_n;


% Load Step
%%%%%%%%%%%

loadStep.r_si = loadStep.r * machine.zn;
loadStep.l_si = loadStep.l * machine.ln;


%% MMC mappings

sqrt3 = sqrt(3);

abc_to_ab0 = [
    2 -1 -1;
    0 sqrt(3) -sqrt(3);
    1 1 1
]/3;

ab0_to_abc = inv(abc_to_ab0);

% two stage ab0 transformation
cell_to_2ab0 = kron(abc_to_ab0, abc_to_ab0);
ab0_to_cell = inv(cell_to_2ab0);

grid_to_cell = kron(eye(3), ones(3, 1));
machine_to_cell = kron(ones(3, 1), eye(3));

l1 = 1/3*mmc.l_leg_g/grid.tr.l;
l2 = 1/3*mmc.l_leg_m/machine.tr.l;

volt_to_cell = [ ...
    [ ...
        1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 0 -1 0 0; 
        0 0 0 0 -1 0;
        0 0 1 0 0 -1;
    ], [ ...
        l1+1 0 0 0;
        0 l1+1 0 0;
        0 0 -l2-1 0;
        0 0 0 -l2-1;
        0 0 0 0;
    ]
];
