% ** Minimum working exaple of RT-GSFit **
%
% As this is a Minimum working example we deliberately have all of the code
% in a single file
%
% Reminder: to run this type of file you need to click the big green "Run"
% button. This is so that the functions at the end of the file are found
%
% Muti-dimension array shape:
% * Matlabe/Fortran are column major and following MDSplus convention, time
% should be in the last index
% * C/Python/Rust are row major, with with time in the first index

clear;

%% Settings

% Choose the pulse and run where to read RT-GSFit setting from
rtgsfit_setup_pluse = 99000230;
rtgsfit_setup_run_name = 'RUN01';

% Choose the pulse to read the magnetic sensors from
% mag_data_pulse = 13343;
% mag_data_pulse = 13923;
mag_data_pulse = 13916;

% Time to start RT-GSFit
time_start = 15e-3;
time_end = 250e-3;
n_timestepping = 20;

% Path to RT-GSFit C library and header file
username = getenv('USER');
if strcmp(username, 'peter.buxton')
    lib_path = [getenv('HOME') '/github/rtgsfit/lib/librtgsfit.so'];
    header_path = [getenv('HOME') '/github/rtgsfit/src/rtgsfit.h'];
elseif strcmp(username, 'alex.prokopyszyn')
    lib_path = [getenv('HOME') '/GitHub/rtgsfit/lib/librtgsfit.so'];
    header_path = [getenv('HOME') '/GitHub/rtgsfit/src/rtgsfit.h'];
elseif strcmp(username, 'filip.janky')
    lib_path = [getenv('HOME') '/ops/rtgsfit/lib/librtgsfit.so'];
    header_path = [getenv('HOME') '/ops/rtgsfit/src/rtgsfit.h'];
end

%% Read RT-GSFit set-up data from MDSplus
mdsconnect('smaug');
mdsopen('RTGSFIT', rtgsfit_setup_pluse);

% sensor names PCS should read
sens_names_pcs = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:SENS_NAMES']);
n_sens_pcs = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:N_SENS_PCS']);

% coil signal names
coil_signal_names = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:COIL_SIGNALS']);
n_coil_signals = length(coil_signal_names);

% I_coil = coil_matrix * I_coil_pcs
% IMPORTANT: there is a matrix transpose because Matlab is column major,
% whereas Python/Rust/C are row major
coil_matrix = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:COIL_MATRIX'])';

initial_condition_flux_norm = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT.INITIAL_COND:FLUX_NORM']);
initial_condition_flux_mask = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT.INITIAL_COND:MASK']);
initial_condition_flux_psi_total = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT.INITIAL_COND:PSI_TOTAL']);

n_coef = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:N_COEF']);
n_lcfs_max = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:N_LCFS_MAX']);

% These variables are not needed to run RT-GSFit, but usefult for plotting 
% results
r_vec = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:R_VEC']);
z_vec = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:Z_VEC']);
n_r = length(r_vec);
n_z = length(z_vec);
coil_names = mdsvalue(['\RTGSFIT::TOP.' rtgsfit_setup_run_name '.PRESHOT:COIL_NAMES']);
n_pf_coils = length(coil_names);

% Check some array shapes
assert(n_sens_pcs == length(sens_names_pcs), 'Error n_sens_pcs should = length(.PRESHOT:SENS_NAMES)');
assert(length(initial_condition_flux_norm)==n_r*n_z, 'error in `initial_condition_flux_norm` shape');
assert(length(initial_condition_flux_mask)==n_r*n_z, 'error in `initial_condition_flux_mask` shape');
assert(length(initial_condition_flux_psi_total)==n_r*n_z, 'error in `initial_condition_flux_psi_total` shape');

mdsclose();

%% Read magnetic sensor data from MDSplus
mdsopen('ST40', mag_data_pulse);

% Get the time vector
time = mdsvalue('dim_of(\MAG::TOP.BEST.BPPROBE.P101:B)');
n_time = length(time);

% Loop over all sensors
sens_data = nan(n_sens_pcs, n_time);  % array ordered using MDSplus standard
for i_signal = 1: n_sens_pcs
    sens_name_pcs = sens_names_pcs(i_signal);
    mds_path = convert_pcs_signal_name_to_mdsplus_path(sens_name_pcs);
    sens_data(i_signal, :) = mdsvalue(mds_path);
end

%% Read PSU data from MDSplus
coil_curr_pcs = nan(n_coil_signals, n_time);  % array ordered using MDSplus standard
for i_signal = 1: n_coil_signals
    coil_signal_name = coil_signal_names(i_signal);
    mds_path = convert_pcs_signal_name_to_mdsplus_path(coil_signal_name);
    coil_curr_pcs(i_signal, :) = mdsvalue(mds_path);
end

mdsclose();

%% Read the plasma current from GSFit
mdsopen('GSFIT', mag_data_pulse);
gsfit_results.time = mdsvalue('\GSFIT::TOP.BEST:TIME');
gsfit_results.plasma_current = mdsvalue('\GSFIT::TOP.BEST.GLOBAL:IP');

mdsclose();
mdsdisconnect();

%% Load RT-GSFit C library
% Load the library
loadlibrary(lib_path, header_path);

% List all available functions in the library
libfunctions('librtgsfit');

%% Prepare input / outputs for RT-GSFit
flux_norm = libpointer('doublePtr', initial_condition_flux_norm);
mask = libpointer('int32Ptr', initial_condition_flux_mask);
flux_total = libpointer('doublePtr', initial_condition_flux_psi_total);
rtgsfit_error = libpointer('doublePtr', 0.0);
lcfs_r = libpointer('doublePtr', zeros(1, n_lcfs_max));
lcfs_z = libpointer('doublePtr', zeros(1, n_lcfs_max));
lcfs_n = libpointer('int32Ptr', zeros(1, n_lcfs_max));
coef = libpointer('doublePtr', zeros(1, n_coef));
flux_boundary = libpointer('doublePtr', 0.0);
plasma_current = libpointer('doublePtr', 0.0);

%% Loop over time
% Starting/ending time index
[~, i_time_start] = min(abs(time - time_start));
[~, i_time_end] = min(abs(time - time_end));

% Figure out the sie of the results
i_counter = 1;
for i_time = i_time_start : n_timestepping : i_time_end
    i_counter = i_counter + 1;
end
n_time_result = i_counter;

% Initialise empty results struct
rtgsfit_results.psi = nan(n_z, n_r, n_time_result);
rtgsfit_results.time = nan(1, n_time_result);
rtgsfit_results.plasma_current = nan(1, n_time_result);
rtgsfit_results.flux_boundary = nan(1, n_time_result);

i_counter = 1;
for i_time = i_time_start : n_timestepping : i_time_end
    % Magnetic sensor data
    meas_pcs = libpointer('doublePtr', sens_data(:, i_time)');

    % Coil currents
    coil_curr = coil_matrix * coil_curr_pcs(:, i_time);
    coil_curr = coil_curr(:)';
    % Manually set the PSH coil to zero
    % coil_curr(end) = 0.0;

    coil_curr = libpointer('doublePtr', coil_curr);

    % Call the function
    ret = calllib( ...
        'librtgsfit', ...
        'rtgsfit', ...
        meas_pcs, ...
        coil_curr, ...
        flux_norm, ...
        mask, ...
        flux_total, ...
        rtgsfit_error, ...
        lcfs_r, ...
        lcfs_z, ...
        lcfs_n, ...
        coef, ...
        flux_boundary, ...
        plasma_current ...
    );

    % Store results
    rtgsfit_results.time(i_counter) = time(i_time);
    rtgsfit_results.psi(:, :, i_counter) = reshape(flux_total.Value, n_r, n_z)';
    rtgsfit_results.plasma_current(i_counter) = plasma_current.Value;
    rtgsfit_results.flux_boundary(i_counter) = flux_boundary.Value;

    % Add to the results counter
    i_counter = i_counter + 1;
end

% Unload the library (this is useful, because if we make a change to the
% binary by recompiling the change only takes effect if we unload and then
% reload the library)
unloadlibrary('librtgsfit');

%% Plot the results

figure;
ax1 = subplot(2,1,1);
hold(ax1, 'on');
plot(rtgsfit_results.time*1e3, rtgsfit_results.plasma_current/1e3);
plot(gsfit_results.time*1e3, gsfit_results.plasma_current/1e3);
xlabel('Time [ms]');
ylabel('Plasma Current [kA]');
legend('RT-GSFit', 'GSFit');
title('Plasma Current Comparison');

ax2 = subplot(2,1,2);
hold(ax2, 'on');
axis(ax2, 'equal');
plot_results(ax1, ax2, r_vec, z_vec, rtgsfit_results);

%% Plot PF coil currents (for checking)
coil_curr_tmp = nan(n_pf_coils, n_time);
for i_time = 1: n_time
    coil_curr_tmp(:, i_time) = coil_matrix * coil_curr_pcs(:, i_time);
end
n_cols = 2;
n_rows = ceil(n_pf_coils / n_cols);
figure; hold on;
for i_pf_coil = 1:n_pf_coils
    subplot(n_rows, n_cols, i_pf_coil);
    plot(time*1e3, coil_curr_tmp(i_pf_coil, :)/1e3);
    title(coil_names(i_pf_coil));
    xlabel('Time [ms]');
    ylabel('Current [kA]');
end

%% Function to convert from PCS name format to MDSplus MAG path
function mds_path = convert_pcs_signal_name_to_mdsplus_path(pcs_signal_name)
    % CONVERT_PCS_SIGNAL_NAME_TO_MDSPLUS_PATH convert a signal using PCS formatting to an MDSplus path.
    %
    %   Example conversions:
    %   'I_BVL_PSU'      -> '\PSU::TOP.BEST.BVL:I'
    %   'I_ROG_INIVC000' -> '\MAG::TOP.BEST.INIVC000:I_ROG'
    %   'PSI_FLOOP_001'  -> '\MAG::TOP.BEST.FLOOP.L001:PSI'

    arguments (Input)
        pcs_signal_name char {mustBeNonempty}
    end

    quantity_system_channel = split(strip(pcs_signal_name, ' '), '_');

    if strcmp(quantity_system_channel(3), 'PSU')
        quantity = char(quantity_system_channel(1));
        system = '';
        channel = [char(quantity_system_channel(2)) ':'];
        tree = 'PSU';
        best = '';
    else
        quantity = char(quantity_system_channel(1));
        system = [char(quantity_system_channel(2)) '.'];
        channel = char(quantity_system_channel(3));
        tree = 'MAG';
        best = 'BEST.';
        if strcmp(system, 'BPPROBE.')
            channel = ['P' channel ':'];
        elseif strcmp(system, 'FLOOP.')
            channel = ['L' channel ':'];
        elseif strcmp(system, 'ROG.')
            channel = [channel ':'];
        end
    end

    mds_path = ['\' tree '::TOP.' best system channel quantity];
end

%% Interactive results plotting
function plot_results(ax1, ax2, r_vec, z_vec, results)
    n_time = length(results.time);
    [~, psi_contour_handle] = contour(ax2, r_vec, z_vec, results.psi(:, :, 1));
    [~, lcfs_contour_handle] = contour(ax2, r_vec, z_vec, results.psi(:, :, 1), [results.flux_boundary(1) results.flux_boundary(1)], Color='Black');
    plasma_current_handle = plot(ax1, [results.time(1)*1e3 results.time(1)*1e3], [0.0 400.0], LineStyle="--", Color='Black');
    uicontrol( ...
        Style='slider', ...
        Min=1, ...
        Max=n_time, ...
        Value=1, ...
        Position=[100 10 350 20], ...
        Callback=@react_to_slider ...
    );

    function react_to_slider(source, ~)
        i_time = round(get(source, 'Value'));
        psi_contour_handle.ZData = results.psi(:, :, i_time);
        lcfs_contour_handle.ZData = results.psi(:, :, i_time);
        lcfs_contour_handle.LevelList = results.flux_boundary(i_time);
        plasma_current_handle.XData = [results.time(i_time)*1e3 results.time(i_time)*1e3];
    end
end