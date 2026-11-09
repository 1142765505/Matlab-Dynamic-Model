function cfg = ck8011d_case_config()
% 单工况入口配置。切换采集时一起修改工况标签、WAV 路径和切削参数。
repo = fileparts(mfilename('fullpath'));
data = fullfile(fileparts(repo),'实验三工况崩刀', ...
    '工况2：刀具11左右新刀-LNMX191940-RR94Tp0501','采集3轮1下刀2正常');
cfg.caseID = 'tool11_acq03';
cfg.outputDir = fullfile(repo,'matlab_result_tool11_acq03_calibrated');
cfg.toolID = '刀具11';
cfg.toolSide = '右';
cfg.acquisitionDate = '2025-11-05';
cfg.acquisition = 3;
cfg.recordedDiameterInfeed = 2; % 原记录“进刀量（轮径）”字段。
cfg.wavFiles = {fullfile(data,'Record1.wav'),fullfile(data,'Record_1_2.wav')};
cfg.recordingStartS = 60;       % 两台设备各自文件的起点偏移，单位 s。
cfg.envelopeChannel = 1;        % CH1 为轴向。原包络方案用它估计相对波动。
cfg.envelopeBlockS = 0.5;
cfg.accelerationPerFS = [10 10 10 NaN]; % 前三路为加速度，CH4 为声音，保留 FS。

% 这三个数值取本地车削计算器默认值，缺少采集3对应的机床日志。
cfg.model.rpm = 30;                   % 工件转速，r/min。
cfg.model.cutting.feedSpeedMmMin = 30; % 沿局部加工路径的进给速度，mm/min。
cfg.model.cutting.diameterMm = 800;   % 局部切削直径，mm。
% 工况表未注明单位。暂按直径减小量为 mm，径向切深为其一半。
cfg.model.cutting.depthOfCutMm = cfg.recordedDiameterInfeed/2;
cfg.model.cutting.enteringAngleDeg = 90; % 局部直线车削近似。
% 按径向、轴向、切向排列。系数为演示输入，待三向测力标定。
cfg.model.cutting.shearNPerMm2 = [600 350 1800];
cfg.model.cutting.edgeNPerMm = [0 0 0];
cfg.model.fs = 64000;
cfg.model.duration = 8;
cfg.model.discard = 2;
cfg.model.showFigures = false;
cfg.comparisonBandHz = [20 4000];
cfg.calibration.enabled = true;
cfg.calibration.channelByDirection = [2 1 3]; % 用户确认：径向、轴向、切向依次 CH2、CH1、CH3。
cfg.calibration.mappingConfirmed = true;
cfg.calibration.splitTimeS = 5;              % 2–5 s 拟合，5–8 s 留出比较。
cfg.calibration.gainFloorRatio = 3e-4;       % 按拟合段能量损失检查后选取的频响下限。
cfg.parameterBasis = '30 r/min、30 mm/min、800 mm 取车削计算器示例值。进刀量 2 暂按直径减小 2 mm 解释。';
end
