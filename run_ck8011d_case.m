function result = run_ck8011d_case(cfg)
% 读取一个实际采集片段，用其相对波动包络调制切削力并仿真。
if nargin==0, cfg = ck8011d_case_config(); end
p = cfg.model;
assert(numel(cfg.wavFiles)==2,'请设置两个立体声 WAV 的路径。');
assert(all(cellfun(@isfile,cfg.wavFiles)),'找不到 WAV，请在 ck8011d_case_config.m 中修改路径。');
validateattributes(cfg.envelopeChannel,{'numeric'},{'scalar','integer','>=',1,'<=',4});
validateattributes(cfg.recordingStartS,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(cfg.envelopeBlockS,{'numeric'},{'scalar','finite','positive'});

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
    'envelopeChannel',cfg.envelopeChannel,'parameterBasis',cfg.parameterBasis);
result = wheel_lathe_lumped_demo(cfg.outputDir,p);

% 同工况比较频率分布。实测为 FS，仿真为 m/s²，分别归一化频带积分。
obs = measured((0:count-1)/fs>=p.discard,:);
sim = result.holderRearAcceleration(:,result.time>=p.discard).';
[f,obsP] = normalized_psd(obs,fs,cfg.comparisonBandHz);
[fSim,simP] = normalized_psd(sim,p.fs,cfg.comparisonBandHz);
assert(isequal(f,fSim),'实测与仿真 PSD 的频率网格不一致。');
[~,ix] = max(obsP,[],1);
obsAC = obs-mean(obs,1);
names = {'CH1','CH2','CH3','CH4'};
summary = table(string(names(:)),sqrt(mean(obsAC.^2,1)).',max(abs(obsAC),[],1).',f(ix), ...
    'VariableNames',{'channel','rms_FS','peak_FS','dominant_frequency_Hz'});
writetable(summary,fullfile(cfg.outputDir,'measured_channel_metrics.csv'));
writetable(array2table([f,obsP,simP], ...
    'VariableNames',{'frequency_Hz','CH1','CH2','CH3','CH4','m11_radial','m11_axial','m11_tangential'}), ...
    fullfile(cfg.outputDir,'normalized_psd_comparison.csv'));
writetable(table(timeS(:),level(:),gain(:),'VariableNames',{'case_time_s','block_rms_FS','gain'}), ...
    fullfile(cfg.outputDir,'measured_load_envelope.csv'));

% 工况清单同时保留来源、计算量与待核实的输入。
info = struct('config',cfg,'wavSampleRateHz',fs,'recordingDurationS',recordingDurationS, ...
    'process',result.process,'simulatedDurationS',p.duration, ...
    'measurementUnit','WAV full scale (FS)','conditioning','median-normalized block RMS');
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
    else
        font='Times New Roman'; suffix='_en';
        titles={'Measured channel spectra','Simulated holder rear spectra','Measured load envelope'};
        freqLabel='Frequency (Hz)'; psdLabel='Normalized PSD (1/Hz)';
        timeLabel='Case time (s)'; gainLabel='Random force gain'; dirs={'Radial','Axial','Tangential'};
    end
    fig = figure('Visible','off','Color','w','Position',[50 50 1100 780]);
    subplot(3,1,1); semilogy(f,max(obsP,eps)); grid on; xlim(cfg.comparisonBandHz);
    xlabel(freqLabel); ylabel(psdLabel); title(titles{1}); legend(names,'Location','northeast');
    subplot(3,1,2); semilogy(f,max(simP,eps)); grid on; xlim(cfg.comparisonBandHz);
    xlabel(freqLabel); ylabel(psdLabel); title(titles{2}); legend(dirs,'Location','northeast');
    subplot(3,1,3); plot(p.loadEnvelope.timeS,p.loadEnvelope.gain,'-o'); grid on;
    xlabel(timeLabel); ylabel(gainLabel); title(titles{3});
    set(findall(fig,'-property','FontName'),'FontName',font);
    print(fig,fullfile(cfg.outputDir,['measurement_comparison' suffix '.svg']),'-dsvg','-vector');
    close(fig);
end
result.measurement = info;
end

function [f,P] = normalized_psd(x,fs,band)
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
nativeF = (0:size(P,1)-1)'*(fs/n);
f = (ceil(band(1)/4)*4:4:band(2)).';
assert(numel(f)>=2,'比较频带至少需要两个 4 Hz 频点。');
P = interp1(nativeF,P,f,'linear');
% 单边因子、窗功率、窗数等常数在同通道带内归一化时相消。
power = trapz(f,P,1);
P = P./max(power,realmin);
end
