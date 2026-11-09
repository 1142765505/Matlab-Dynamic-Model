function result = run_ck8011d_case(cfg)
% 读取一个实际采集片段，用其相对波动包络调制切削力并仿真。
if nargin==0, cfg = ck8011d_case_config(); end
p = cfg.model;
assert(numel(cfg.wavFiles)==2,'请设置两个立体声 WAV 的路径。');
assert(all(cellfun(@isfile,cfg.wavFiles)),'找不到 WAV，请在 ck8011d_case_config.m 中修改路径。');
validateattributes(cfg.envelopeChannel,{'numeric'},{'scalar','integer','>=',1,'<=',4});
validateattributes(cfg.recordingStartS,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(cfg.envelopeBlockS,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.accelerationPerFS,{'numeric'},{'vector','real'});
scale = cfg.accelerationPerFS(:).';
assert(ismember(numel(scale),[1 4]),'请提供统一满量程系数，或按 CH1–CH4 提供四个系数。');
if isscalar(scale), scale = repmat(scale,1,4); end
assert(all(isnan(scale) | (isfinite(scale) & scale>0)) && any(isfinite(scale)), ...
    '加速度通道系数须为正数，非加速度通道用 NaN。');

% 两设备分别截取各自文件的时间段，仅使用各通道幅值和功率谱。
recordingDurationS = zeros(1,2);
for device = 1:2
    info = audioinfo(cfg.wavFiles{device});
    assert(info.NumChannels==2,'每个 WAV 必须包含两个通道。');
    recordingDurationS(device) = info.Duration;
    if device==1
        fs = info.SampleRate;
        count = round(p.duration*fs);
        measured = zeros(count,4);
    end
    assert(info.SampleRate==fs,'两台设备的采样率必须一致。');
    first = round(cfg.recordingStartS*fs)+1;
    last = first+count-1;
    assert(last<=info.TotalSamples,'所选时间段超出录音长度。');
    measured(:,2*device-1:2*device) = audioread(cfg.wavFiles{device},[first last]);
end
assert(all(isfinite(measured),'all'),'WAV 含非有限值。');

% 逐块去均值后算 RMS，再除以片段 RMS 中位数得到无量纲增益。
block = round(cfg.envelopeBlockS*fs);
assert(block>=2 && block<=count/2,'包络块长应至少两点，且片段至少包含两块。');
start = 1:block:count-1;
stop = min(start+block-1,count);
level = zeros(size(start));
for k=1:numel(start)
    x = measured(start(k):stop(k),cfg.envelopeChannel);
    level(k) = sqrt(mean((x-mean(x)).^2));
end
reference = median(level);
assert(reference>0,'所选包络通道的分块 RMS 中位数为零。');
timeS = (start+stop-2)/(2*fs);
gain = level/reference;
assert(cfg.comparisonBandHz(1)>0 && diff(cfg.comparisonBandHz)>0 && ...
    cfg.comparisonBandHz(2)<min(fs,p.fs)/2,'比较频带须位于两种采样率的奈奎斯特频率以内。');
p.loadEnvelope = struct('timeS',[0,timeS,p.duration], 'gain',[gain(1),gain,gain(end)]);
p.measurement = struct('caseID',cfg.caseID,'toolID',cfg.toolID,'side',cfg.toolSide, ...
    'date',cfg.acquisitionDate,'acquisition',cfg.acquisition,'sourceFiles',{cfg.wavFiles}, ...
    'recordingStartS',cfg.recordingStartS,'sampleRateHz',fs, ...
    'envelopeChannel',cfg.envelopeChannel,'parameterBasis',cfg.parameterBasis, ...
    'accelerationPerFS',scale,'calibratedUnit','m/s^2');
obsTime = (0:count-1)'/fs;
keep = obsTime>=p.discard;
obsTime = obsTime(keep);
obs = measured(keep,:);
obsAC = obs-mean(obs,1);
obsAcceleration = obsAC.*scale;
isCalibrated = isfield(cfg,'calibration') && cfg.calibration.enabled;
axesForPlot = 1:3;
channelsForPlot = find(isfinite(scale));
if isCalibrated
    c = cfg.calibration;
    validateattributes(c.channelByDirection,{'numeric'},{'size',[1 3],'integer','>=',0,'<=',4});
    validateattributes(c.gainFloorRatio,{'numeric'},{'scalar','real','finite','positive','<=',1});
    validateattributes(c.splitTimeS,{'numeric'},{'scalar','real','finite','>',p.discard,'<',p.duration});
    axesForPlot = find(c.channelByDirection);
    channelsForPlot = c.channelByDirection(axesForPlot);
    assert(~isempty(axesForPlot) && numel(unique(channelsForPlot))==numel(channelsForPlot), ...
        '至少指定一个通道。同一通道只能对应一个方向。');
    assert(all(isfinite(scale(channelsForPlot))),'标定映射只能引用有加速度满量程系数的通道。');
    train = obsTime<c.splitTimeS;
    % 逐窗去均值，预先去掉的全段均值不会改变训练 PSD。
    [trainF,trainP] = welch_psd(obs(train,channelsForPlot).*scale(channelsForPlot),fs,cfg.comparisonBandHz);
    c.frequencyHz = trainF;
    c.accelerationPsd = zeros(numel(trainF),3);
    c.accelerationPsd(:,axesForPlot) = trainP;
    p.forceCalibration = c;
end
result = wheel_lathe_lumped_demo(cfg.outputDir,p);

% 频谱按各通道的带内功率归一化，时域幅值按满量程换算。
simTime = result.time(result.time>=p.discard).';
sim = result.holderRearAcceleration(:,result.time>=p.discard).';
[f,obsP] = normalized_psd(obs,fs,cfg.comparisonBandHz);
[fSim,simP] = normalized_psd(sim,p.fs,cfg.comparisonBandHz);
assert(isequal(f,fSim),'实测与仿真 PSD 的频率网格不一致。');
[~,ix] = max(obsP,[],1);
names = {'CH1','CH2','CH3','CH4'};
summary = table(string(names(:)),sqrt(mean(obsAcceleration.^2,1)).', ...
    max(abs(obsAcceleration),[],1).',f(ix),sqrt(mean(obsAC.^2,1)).', ...
    max(abs(obsAC),[],1).', ...
    'VariableNames',{'channel','rms_m_s2','peak_m_s2','dominant_frequency_Hz','rms_FS','peak_FS'});
writetable(summary,fullfile(cfg.outputDir,'measured_channel_metrics.csv'));
writetable(array2table([f,obsP,simP], ...
    'VariableNames',{'frequency_Hz','CH1','CH2','CH3','CH4','m11_radial','m11_axial','m11_tangential'}), ...
    fullfile(cfg.outputDir,'normalized_psd_comparison.csv'));
writetable(table(timeS(:),level(:)*scale(cfg.envelopeChannel),level(:),gain(:), ...
    'VariableNames',{'case_time_s','block_rms_m_s2','block_rms_FS','gain'}), ...
    fullfile(cfg.outputDir,'measured_load_envelope.csv'));

if isCalibrated
    obsTest = obsTime>=c.splitTimeS;
    simTest = simTime>=c.splitTimeS & simTime<p.duration;
    [testF,testObs] = welch_psd(obsAcceleration(obsTest,channelsForPlot),fs,cfg.comparisonBandHz);
    [~,testSim] = welch_psd(sim(simTest,axesForPlot),p.fs,cfg.comparisonBandHz);
    obsRms = sqrt(trapz(testF,testObs,1));
    simRms = sqrt(trapz(testF,testSim,1));
    [~,measuredPeak] = max(testObs,[],1);
    [~,simulatedPeak] = max(testSim,[],1);
    dirs = ["radial";"axial";"tangential"];
    assessment = table(string(names(channelsForPlot)).',dirs(axesForPlot),obsRms.',simRms.', ...
        (simRms./obsRms-1).'*100,testF(measuredPeak),testF(simulatedPeak), ...
        'VariableNames',{'channel','direction','measured_band_rms_m_s2','simulated_band_rms_m_s2', ...
        'rms_relative_error_pct','measured_peak_Hz','simulated_peak_Hz'});
    writetable(assessment,fullfile(cfg.outputDir,'calibration_holdout_metrics.csv'));
    spectrum = array2table([testF,testObs,testSim], ...
        'VariableNames',[{'frequency_Hz'},strcat('measured_',names(channelsForPlot)), ...
        cellstr('simulated_'+dirs(axesForPlot)).']);
    writetable(spectrum,fullfile(cfg.outputDir,'calibration_holdout_psd.csv'));
    disp(assessment);
end

% 工况清单同时保留来源、计算量与待核实的输入。
info = struct('config',cfg,'wavSampleRateHz',fs,'recordingDurationS',recordingDurationS, ...
    'process',result.process,'simulatedDurationS',p.duration, ...
    'measurementUnit','m/s^2 for calibrated acceleration channels; FS for other channels', ...
    'rawMeasurementUnit','WAV full scale (FS)', ...
    'conditioning','median-normalized block RMS');
fid = fopen(fullfile(cfg.outputDir,'case_metadata.json'),'w','n','UTF-8');
assert(fid~=-1,'无法写入工况清单。');
closeFile = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(info,'PrettyPrint',true));
clear closeFile

for lang = 1:2
    if lang==1
        font='宋体'; suffix='';
        titles={'实测信号频谱（按通道归一化）','刀杆后端仿真频谱（按方向归一化）','实测相对波动包络'};
        freqLabel='频率 (Hz)'; psdLabel='归一化功率谱密度 (1/Hz)';
        timeLabel='工况时间 (s)'; gainLabel='随机力增益'; dirs={'径向','轴向','切向'};
        accelLabel='加速度 (m/s^2)'; timeTitles={'实测加速度（满量程换算）','刀杆后端仿真加速度'};
    else
        font='Times New Roman'; suffix='_en';
        titles={'Measured channel spectra','Simulated holder rear spectra','Measured load envelope'};
        freqLabel='Frequency (Hz)'; psdLabel='Normalized PSD (1/Hz)';
        timeLabel='Case time (s)'; gainLabel='Random force gain'; dirs={'Radial','Axial','Tangential'};
        accelLabel='Acceleration (m/s^2)'; timeTitles={'Measured acceleration after full-scale conversion','Simulated holder rear acceleration'};
    end
    fig = figure('Visible','off','Color','w','Position',[50 50 1100 780]);
    subplot(3,1,1); semilogy(f,max(obsP,eps)); grid on; xlim(cfg.comparisonBandHz);
    xlabel(freqLabel); ylabel(psdLabel); title(titles{1}); legend(names,'Location','northeast');
    subplot(3,1,2); semilogy(f,max(simP(:,axesForPlot),eps)); grid on; xlim(cfg.comparisonBandHz);
    xlabel(freqLabel); ylabel(psdLabel); title(titles{2}); legend(dirs(axesForPlot),'Location','northeast');
    subplot(3,1,3); plot(p.loadEnvelope.timeS,p.loadEnvelope.gain,'-o'); grid on;
    if isCalibrated
        if lang==1, gainLabel='分块相对 RMS'; else, gainLabel='Relative block RMS'; end
    end
    xlabel(timeLabel); ylabel(gainLabel); title(titles{3});
    set(findall(fig,'-property','FontName'),'FontName',font);
    print(fig,fullfile(cfg.outputDir,['measurement_comparison' suffix '.svg']),'-dsvg','-vector');
    close(fig);

    % 标定模式下只绘制映射方向，展示留出段且采用同一加速度刻度。
    firstTime = p.discard;
    if isCalibrated
        firstTime = c.splitTimeS;
        if lang==1
            timeTitles={'留出段实测加速度','留出段刀杆后端仿真加速度'};
        else
            timeTitles={'Measured acceleration, holdout','Simulated holder rear acceleration, holdout'};
        end
        if ~c.mappingConfirmed
            if lang==1
                timeTitles={'留出段实测加速度（暂定通道对应）','留出段刀杆后端仿真加速度'};
            else
                timeTitles={'Measured holdout, provisional channel mapping','Simulated holder rear, holdout'};
            end
        end
    end
    lastTime = min(firstTime+0.25,p.duration);
    obsShow = obsTime>=firstTime & obsTime<=lastTime;
    simShow = simTime>=firstTime & simTime<=lastTime;
    fig = figure('Visible','off','Color','w','Position',[50 50 1100 660]);
    axMeasured = subplot(2,1,1);
    plot(obsTime(obsShow),obsAcceleration(obsShow,channelsForPlot)); grid on;
    ylabel(accelLabel); title(timeTitles{1}); legend(names(channelsForPlot),'Location','northeast');
    axSim = subplot(2,1,2);
    plot(simTime(simShow),sim(simShow,axesForPlot)); grid on;
    xlabel(timeLabel); ylabel(accelLabel); title(timeTitles{2}); legend(dirs(axesForPlot),'Location','northeast');
    if isCalibrated
        limit = 1.05*max([max(abs(obsAcceleration(obsShow,channelsForPlot)),[],'all'), ...
            max(abs(sim(simShow,axesForPlot)),[],'all'),eps]);
        linkaxes([axMeasured axSim],'xy'); ylim([-limit limit]);
    else
        linkaxes([axMeasured axSim],'x');
    end
    xlim([firstTime lastTime]);
    set(findall(fig,'-property','FontName'),'FontName',font);
    print(fig,fullfile(cfg.outputDir,['measurement_time_comparison' suffix '.svg']),'-dsvg','-vector');
    close(fig);

    if isCalibrated
        fig = figure('Visible','off','Color','w','Position',[50 50 1100 270*numel(axesForPlot)]);
        if lang==1
            head='留出段绝对功率谱'; powerLabel='PSD ((m/s^2)^2/Hz)'; keys={'实测','仿真'};
            if ~c.mappingConfirmed, head=[head '（暂定通道对应）']; end
        else
            head='Holdout absolute power spectrum'; powerLabel='PSD ((m/s^2)^2/Hz)'; keys={'Measured','Simulated'};
            if ~c.mappingConfirmed, head=[head ' (provisional mapping)']; end
        end
        for k=1:numel(axesForPlot)
            subplot(numel(axesForPlot),1,k);
            semilogy(testF,max([testObs(:,k),testSim(:,k)],realmin)); grid on;
            xlim(cfg.comparisonBandHz); xlabel(freqLabel); ylabel(powerLabel);
            title([head '  ' names{channelsForPlot(k)} ' / ' dirs{axesForPlot(k)}]);
            legend(keys,'Location','northeast');
        end
        set(findall(fig,'-property','FontName'),'FontName',font);
        print(fig,fullfile(cfg.outputDir,['calibration_holdout' suffix '.svg']),'-dsvg','-vector');
        close(fig);
    end
end
result.measurement = info;
end

function [f,P] = normalized_psd(x,fs,band)
[f,P] = welch_psd(x,fs,band);
P = P./max(trapz(f,P,1),realmin);
end

function [f,P] = welch_psd(x,fs,band)
% Welch 分段估计：0.25 s Hann 窗、50% 重叠、逐窗去均值。
n = round(fs/4);
assert(size(x,1)>=n,'分析段不足 0.25 秒。');
w = 0.5-0.5*cos(2*pi*(0:n-1)'/n);
ids = 1:floor(n/2):size(x,1)-n+1;
P = zeros(floor(n/2)+1,size(x,2));
for first = ids
    segment = x(first:first+n-1,:);
    y = fft((segment-mean(segment,1)).*w,[],1);
    P = P+abs(y(1:size(P,1),:)).^2;
end
P = P/(numel(ids)*fs*sum(w.^2));
if mod(n,2)==0
    P(2:end-1,:) = 2*P(2:end-1,:);
else
    P(2:end,:) = 2*P(2:end,:);
end
nativeF = (0:size(P,1)-1)'*(fs/n);
f = (ceil(band(1)/4)*4:4:band(2)).';
assert(numel(f)>=2,'比较频带至少需要两个 4 Hz 频点。');
P = interp1(nativeF,P,f,'linear');
end
