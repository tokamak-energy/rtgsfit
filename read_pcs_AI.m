PCSIO_file = 'pss/PCS_IO.xlsx';
PCSIO_path = fullfile(pcs_home_path, PCSIO_file);
if not(isfile(PCSIO_path))
    error(['File not found: ', PCSIO_file]);
end

% Get the name of tabs in PCSIO_file
sheets = cellstr(sheetnames(PCSIO_path));

%% AI list
assert(numel(find(contains(sheets, 'PCS1 AI'))) == 1)
[AInum1, AItxt1] = xlsread(PCSIO_path, 'PCS1 AI', 'A1:N163');
List1 = AItxt1(4:end, 2);

assert(numel(find(contains(sheets, 'PCS2 AI'))) == 1)
[AInum2, AItxt2] = xlsread(PCSIO_path, 'PCS2 AI', 'A1:N131');
List2 = AItxt2(4:end, 2);

AI = [];
AI.List = [List1; List2];
AI.Nin1 = length(List1); % shorts
AI.Nin2 = length(List2); % shorts
AInum = [AInum1; AInum2];
% We are missing the Stop data
AI.Nin = AI.Nin1 + AI.Nin2; % shorts
% The following line caused an error, just check it
if (size(AInum, 1) < AI.Nin)
    % this can never happen, FAIL
    error([mfilename ' wrong sizes of PCS_IO data']);
end
rtgsfit_path = '/home/filip.janky/rt-gsfit/rtgsfit/data';
rtgsfit_sensor_list = 'included_sensors.txt';
rtgsfit_sens_file = fullfile(rtgsfit_path, rtgsfit_sensor_list);
fid = fopen(rtgsfit_sens_file,'r');
sens_list = textscan(fid, '%s');
sens_list = sens_list{:};
sens_list_idx = zeros(size(sens_list));
for rtgsfit_list = 1:length(sens_list)
    if ~contains(sens_list{rtgsfit_list},'REG_')
        sens_list_idx(rtgsfit_list) = find(strcmp(AI.List, sens_list{rtgsfit_list}));
    end
end
fclose(fid);
fid = fopen('~/rt-gsfit/rtgsfit/data/sensor_index.txt','w');
for ii=1:length(sens_list_idx)
    fprintf(fid, "%d\n", sens_list_idx(ii) - 1);
end
fclose(fid);