function result = wheel_lathe_lumped_demo(outdir, p)
% CK8011D 轮对车削 33 自由度线性受迫振动模型，仅需基础 MATLAB。
% 直接运行生成演示结果。实测包络辅助的单工况见 run_ck8011d_case。
% 设置 p.showFigures=false 仅导出 SVG。结构求解采用 kg、N、m、s。
% 默认结构和切削力参数用于演示，机床参数需实测标定。
% 每个质量体有径向、轴向、切向平动，顶紧预载用切线刚度等效。

if nargin < 1 || isempty(outdir)
    outdir = fullfile(fileparts(mfilename('fullpath')), 'matlab_result_ck8011d');
end
if nargin < 2, p = struct; end
p = defaults(p);
if ~exist(outdir,'dir'), mkdir(outdir); end
assert(p.duration > p.discard && p.discard >= 0);
assert(p.fs > 2*p.forceNoiseCornerHz && p.fs > 2*p.plotMaxHz);
assert(p.zetaLink > 0 && all(p.massKg > 0));
rng(p.seed,'twister');

% m6 为刀杆前端，m11 为远离刀片的刀杆后端测点。
names = {'left_spindle_plate','left_centre','right_centre','wheelset', ...
         'carriage','holder','tip','bed','right_spindle_plate','cassette','holder_rear'};
hardware = struct('machine','青重机床 CK8011D', ...
    'insert','Seco LNMX191940-RR94 TP0501','cassette','Seco TRWR175-CA-19', ...
    'holder','Seco R175.33-5055');
nbody = numel(names); nd = 3*nbody;
M = diag(kron(p.massKg(:),ones(3,1)));
K = zeros(nd); C = zeros(nd);
% 每行 [i,j,径向刚度,轴向刚度,切向刚度]，j=0 表示地基。
% 顶尖连接是预紧工作点附近的等效刚度，刀具链拆分保持原串联静刚度。
links = [1 8 2e8 3e8 2e8; ...
         2 8 1e8 1.2e8 1e8; ...
         3 8 1e8 1.2e8 1e8; ...
         4 2 1.4e8 1e8 1.4e8; ...
         4 3 1.4e8 1e8 1.4e8; ...
         4 1 8e6 3e6 3e7; ...
         5 8 1.5e8 1.2e8 1e8; ...
         7 10 2.4e7 1.8e7 3e7; ...
         10 6 2.4e7 1.8e7 3e7; ...
         6 11 1.2e8 1e8 8e7; ...
         11 5 4e7 (4/3)*2.5e7 (4/3)*2e7; ...
         8 0 8e8 8e8 6e8; ...
         9 8 2e8 3e8 2e8; ...
         4 9 8e6 3e6 3e7];
if isfield(p,'links'), links = p.links; end
if all(links(:,1:2) <= 9,'all')
    links = expand_legacy_links(links);
end
p.links = links;
for j = 1:size(links,1)
    i1 = links(j,1); i2 = links(j,2); k = links(j,3:5);
    assert(all(k > 0));
    meq = p.massKg(i1);
    if i2 ~= 0
        meq = p.massKg(i1)*p.massKg(i2)/(p.massKg(i1)+p.massKg(i2));
    end
    % 用两端约化质量估计连接阻尼，zetaLink 是连接参数。
    c = 2*p.zetaLink*sqrt(k*meq);
    [K,C] = add_link(K,C,i1,i2,k,c);
end
assert(norm(K-K','fro') < 1e-8*norm(K,'fro'));
chol(M); chol(K); chol(C);
fn = sort(sqrt(real(eig(K,M)))/(2*pi));
fprintf('MATLAB: %s\n',version);
fprintf('Linear lumped model: %d DOF, mode range %.2f--%.2f Hz\n', ...
        nd,fn(1),fn(end));
fprintf('Samples per highest undamped mode: %.1f\n',p.fs/fn(end));
if p.fs/fn(end) < 20
    warning('Increase fs and check time-step convergence for acceleration.');
end

tip = dofs(7); holder = dofs(6); wheel = dofs(4); holderRear = dofs(11);
B = zeros(nd,3); B(tip,:) = eye(3); B(wheel,:) = -eye(3);
assert(norm(reshape(sum(reshape(B,3,nbody,3),2),3,3),'fro') < 1e-12);
% 刀片与轮对承受同一切削力对，符号相反。
dt = 1/p.fs; t = 0:dt:p.duration; N = numel(t);
[meanForce,process] = cutting_load(p);
calibration = struct('enabled',false);
if isfield(p,'forceCalibration')
    [dynamicForce,calibration] = calibrated_force(p,M,C,K,B,holderRear,N);
    force = meanForce+dynamicForce;
    gain = ones(size(t)); % 此分支只用训练段频谱，留出段不参与载荷生成。
else
    fr = p.rpm/60;
    modulation = zeros(1,N);
    for j = 1:numel(p.orders)
        modulation = modulation + p.orderFraction(j)* ...
            sin(2*pi*p.orders(j)*fr*t + p.orderPhase(j));
    end
    rho = exp(-2*pi*p.forceNoiseCornerHz/p.fs);
    eta = filter(1-rho,[1 -rho],randn(3,N),[],2);
    eta = bsxfun(@minus,eta,mean(eta,2));
    eta = bsxfun(@rdivide,eta,std(eta,0,2));
    gain = ones(size(t));
    if isfield(p,'loadEnvelope')
        e = p.loadEnvelope;
        assert(numel(e.timeS)==numel(e.gain) && all(diff(e.timeS)>0) && ...
            e.timeS(1)<=t(1) && e.timeS(end)>=t(end) && ...
            all(isfinite(e.gain) & e.gain>=0),'载荷包络必须非负且覆盖仿真时间。');
        gain = interp1(e.timeS,e.gain,t,'linear');
    end
    % 无量纲包络只调制随机力分量，切削力仍为给定外载。
    force = meanForce.*(1+modulation+p.noiseFraction*eta.*gain);
end

% Newmark 平均加速度法。有效质量矩阵恒定，仅分解一次。
beta = 1/4; gamma = 1/2;
L = chol(M+gamma*dt*C+beta*dt^2*K,'lower');
q = K\(B*force(:,1)); v = zeros(nd,1); acc = zeros(nd,1);
tipA = zeros(3,N); holderA = zeros(3,N); holderRearA = zeros(3,N);
relativeU = zeros(3,N);
relativeU(:,1) = q(tip)-q(wheel);
maxResidual = 0; maxRhs = 0;
tic;
for j = 2:N
    qp = q+dt*v+dt^2*(0.5-beta)*acc;
    vp = v+dt*(1-gamma)*acc;
    f = B*force(:,j);
    an = L'\(L\(f-C*vp-K*qp));
    q = qp+beta*dt^2*an;
    v = vp+gamma*dt*an;
    acc = an;
    tipA(:,j) = acc(tip);
    holderA(:,j) = acc(holder);
    holderRearA(:,j) = acc(holderRear);
    relativeU(:,j) = q(tip)-q(wheel);
    if mod(j,1000)==0
        maxResidual = max(maxResidual,norm(M*acc+C*v+K*q-f,inf));
        maxRhs = max(maxRhs,norm(f,inf));
    end
end
elapsed = toc;
equilibriumResidual = maxResidual/max(maxRhs,1);
assert(equilibriumResidual < 1e-6,'Dynamic equilibrium residual too large.');

% 加速度频响的第 i 行为输出方向，第 j 列为切削力对的方向。
freq = linspace(0,p.plotMaxHz,p.frfPoints);
Htip = complex(zeros(3,3,numel(freq)));
Hrear = Htip;
for j = 1:numel(freq)
    w = 2*pi*freq(j);
    U = (K-w^2*M+1i*w*C)\B;
    Htip(:,:,j) = -w^2*U(tip,:);
    Hrear(:,:,j) = -w^2*U(holderRear,:);
end
keep = t >= p.discard;
[freqFFT,ampFFT] = amplitude_spectrum(tipA(:,keep),p.fs);
[~,rearAmpFFT] = amplitude_spectrum(holderRearA(:,keep),p.fs);
rmsTip = sqrt(mean(tipA(:,keep).^2,2));
peakTip = max(abs(tipA(:,keep)),[],2);
axesNames = {'radial';'axial';'tangential'};
metrics = table(axesNames,rmsTip,peakTip,rmsTip/9.80665, ...
                'VariableNames',{'direction','RMS_m_s2','Peak_m_s2','RMS_g'});
rearMetrics = metrics;
rearMetrics.RMS_m_s2 = sqrt(mean(holderRearA(:,keep).^2,2));
rearMetrics.Peak_m_s2 = max(abs(holderRearA(:,keep)),[],2);
rearMetrics.RMS_g = rearMetrics.RMS_m_s2/9.80665;
disp(metrics);
fprintf('Integration %.2f s; normalized equilibrium residual %.3g\n', ...
        elapsed,equilibriumResidual);

% 以独立算例检查时间积分和频谱幅值缩放。
validation = validate_integrator();
assert(validation.relativePhasorError < 0.01,'Harmonic validation failed.');
assert(validation.spectrumAmplitudeError < 1e-10,'FFT scaling check failed.');
fprintf('Analytic harmonic check error %.4g; FFT amplitude error %.3g\n', ...
        validation.relativePhasorError,validation.spectrumAmplitudeError);

result = struct('parameters',p,'hardware',hardware,'process',process,'calibration',calibration, ...
    'loadEnvelopeGain',gain,'bodyNames',{names},'links',links, ...
    'M',M,'C',C,'K',K,'Bcut',B,'naturalFrequencyHz',fn,'time',t, ...
    'cuttingForceN',force,'tipAcceleration',tipA, ...
    'holderAcceleration',holderA,'holderRearAcceleration',holderRearA, ...
    'relativeDisplacement',relativeU, ...
    'frequencyHz',freq,'tipAccelerance',Htip,'holderRearAccelerance',Hrear, ...
    'fftFrequencyHz',freqFFT,'tipAmplitudeSpectrum',ampFFT, ...
    'holderRearAmplitudeSpectrum',rearAmpFFT, ...
    'metrics',metrics,'holderRearMetrics',rearMetrics,'validation',validation, ...
    'equilibriumResidual',equilibriumResidual);
save(fullfile(outdir,'simulation.mat'),'result','-v7');
writetable(metrics,fullfile(outdir,'metrics.csv'));
writetable(rearMetrics,fullfile(outdir,'holder_rear_metrics.csv'));
writetable(array2table([t(:),holderRearA.'], ...
    'VariableNames',{'time_s','radial_m_s2','axial_m_s2','tangential_m_s2'}), ...
    fullfile(outdir,'holder_rear_acceleration.csv'));
writetable(array2table([freqFFT(:),rearAmpFFT.'], ...
    'VariableNames',{'frequency_Hz','radial_m_s2','axial_m_s2','tangential_m_s2'}), ...
    fullfile(outdir,'holder_rear_amplitude_spectrum.csv'));
writetable(table(t(:),tipA(1,:)',tipA(2,:)',tipA(3,:)', ...
    'VariableNames',{'time_s','radial_m_s2','axial_m_s2','tangential_m_s2'}), ...
    fullfile(outdir,'tip_acceleration.csv'));
writetable(table(freqFFT(:),ampFFT(1,:)',ampFFT(2,:)',ampFFT(3,:)', ...
    'VariableNames',{'frequency_Hz','radial_m_s2','axial_m_s2','tangential_m_s2'}), ...
    fullfile(outdir,'tip_amplitude_spectrum.csv'));
writetable(table((1:nd)',fn,'VariableNames',{'mode','frequency_Hz'}), ...
    fullfile(outdir,'natural_frequencies.csv'));

short = t >= p.discard & t <= p.discard+0.25;
directionTitles = {'径向','轴向','切向'};
visible = freqFFT <= p.plotMaxHz;
for point = 1:2
    if point == 1
        nameZh = '刀尖'; nameEn = 'tool-tip'; file = 'acceleration.svg';
        signal = tipA; spectrum = ampFFT;
    else
        nameZh = '刀杆后端'; nameEn = 'holder rear'; file = 'holder_rear_acceleration.svg';
        signal = holderRearA; spectrum = rearAmpFFT;
    end
    fig = figure('Name',[nameZh '加速度'],'NumberTitle','off', ...
        'Visible','off','Color','w','Position',[50 50 1100 850]);
    for j = 1:3
        subplot(3,2,2*j-1);
        plot(t(short),signal(j,short),'LineWidth',0.7); grid on;
        setappdata(gca,'PlotLabels', ...
            {[directionTitles{j} nameZh '加速度时程'], '时间 (s)', '加速度 (m/s^2)'; ...
             [axesNames{j} ' ' nameEn ' acceleration'], 'Time (s)', 'Acceleration (m/s^2)'});
        subplot(3,2,2*j);
        plot(freqFFT(visible),spectrum(j,visible),'LineWidth',0.7); grid on;
        xlim([0 p.plotMaxHz]);
        setappdata(gca,'PlotLabels', ...
            {[directionTitles{j} nameZh '加速度幅度谱'], '频率 (Hz)', '单边幅值 (m/s^2)'; ...
             [axesNames{j} ' ' nameEn ' acceleration spectrum'], 'Frequency (Hz)', 'Single-sided amplitude (m/s^2)'});
    end
    save_result_figure(fig,outdir,file,p.showFigures);
end

fig = figure('Name','单位力加速度频响','NumberTitle','off', ...
    'Visible','off','Color','w','Position',[50 50 950 650]);
for j = 1:3
    subplot(3,1,j);
    semilogy(freq,max(squeeze(abs(Htip(j,j,:))),1e-16)); grid on;
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '单位力加速度频响'], '频率 (Hz)', '加速度频响幅值 ((m/s^2)/N)'; ...
         [axesNames{j} ' unit-force accelerance'], 'Frequency (Hz)', 'Accelerance magnitude ((m/s^2)/N)'});
end
save_result_figure(fig,outdir,'unit_force_frf.svg',p.showFigures);

% 用较长时窗展示切削力的低频阶次调制。
forceWindow = t >= p.discard & t <= p.discard+2;
fig = figure('Name','三向切削力','NumberTitle','off', ...
    'Visible','off','Color','w','Position',[50 50 950 650]);
for j = 1:3
    subplot(3,1,j);
    plot(t(forceWindow),force(j,forceWindow),'LineWidth',0.7); grid on;
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '切削力时程'], '时间 (s)', '切削力 (N)'; ...
         [axesNames{j} ' cutting force'], 'Time (s)', 'Cutting force (N)'});
end
save_result_figure(fig,outdir,'cutting_force.svg',p.showFigures);

fig = figure('Name','刀具与轮对响应','NumberTitle','off', ...
    'Visible','off','Color','w','Position',[50 50 1100 850]);
for j = 1:3
    subplot(3,2,2*j-1);
    plot(t(short),tipA(j,short),t(short),holderA(j,short), ...
        t(short),holderRearA(j,short),'LineWidth',0.7); grid on;
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '刀尖与刀杆前后端加速度'], '时间 (s)', '加速度 (m/s^2)'; ...
         [axesNames{j} ' tool and holder acceleration'], 'Time (s)', 'Acceleration (m/s^2)'});
    if j == 1
        lg = legend('刀尖','刀杆前端','刀杆后端','Location','best');
        setappdata(lg,'PlotLabels', ...
            {'刀尖','刀杆前端','刀杆后端'; 'Tool tip','Holder front','Holder rear'});
    end
    subplot(3,2,2*j);
    plot(t(short),1e6*relativeU(j,short),'LineWidth',0.7); grid on;
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '刀尖相对轮对的位移'], '时间 (s)', '相对位移 ({\mu}m)'; ...
         [axesNames{j} ' displacement: tip minus wheelset'], 'Time (s)', 'Relative displacement ({\mu}m)'});
end
save_result_figure(fig,outdir,'tool_response.svg',p.showFigures);

fig = figure('Name','未阻尼固有频率','NumberTitle','off', ...
    'Visible','off','Color','w','Position',[50 50 950 500]);
stem(1:nd,fn,'filled'); grid on;
xlim([0.5 nd+0.5]); xticks(1:2:nd);
setappdata(gca,'PlotLabels', ...
    {'未阻尼固有频率', '模态阶次', '固有频率 (Hz)'; ...
     'Undamped natural frequencies', 'Mode index', 'Natural frequency (Hz)'});
save_result_figure(fig,outdir,'natural_frequencies.svg',p.showFigures);
fprintf('Saved illustrative results to: %s\n',outdir);
end

function [force,fit] = calibrated_force(p,M,C,K,B,holderRear,N)
% 在给定结构频响下，反推可复现实测谱量级的等效力谱。
c = p.forceCalibration;
f = c.frequencyHz(:);
forcePsd = zeros(numel(f),3);
force = zeros(3,N);
bins = 0:N-1;
fftFrequency = min(bins,N-bins)*p.fs/N;
% 平均加速度 Newmark 的离散频响采用双线性变换后的角频率。
omega = 2*p.fs*tan(pi*f/p.fs);
for d = find(c.channelByDirection(:).')
    H = zeros(size(f));
    for k = 1:numel(f)
        w = omega(k);
        q = (K-w^2*M+1i*w*C)\B(:,d);
        H(k) = abs(w^2*q(holderRear(d)));
    end
    lower = c.gainFloorRatio*max(H);
    assert(lower>0,'所选输入和输出之间没有可用于标定的频响。');
    forcePsd(:,d) = c.accelerationPsd(:,d)./max(H,lower).^2;
    density = interp1(f,forcePsd(:,d),fftFrequency,'linear',0);
    phase = fft(randn(1,N));
    phase = phase./max(abs(phase),eps);
    % 单边 PSD 对应共轭 FFT 两侧各分一半功率。随机相位由种子决定。
    force(d,:) = ifft(phase.*sqrt(density*p.fs*N/2),'symmetric');
end
fit = struct('enabled',true,'channelByDirection',c.channelByDirection, ...
    'mappingConfirmed',c.mappingConfirmed,'frequencyHz',f,'forcePsd_N2PerHz',forcePsd, ...
    'dynamicForceRmsN',sqrt(mean(force.^2,2)));
end

function [force,process] = cutting_load(p)
% 保留旧的给定力入口，提供工艺参数时改用剪切项与刃口项。
validateattributes(p.rpm,{'numeric'},{'scalar','real','finite','nonnegative'});
force = p.meanForceN(:);
process = struct('mode','prescribed_force','rpm',p.rpm,'meanForceN',force);
if ~isfield(p,'cutting'), return; end
assert(p.rpm>0,'工艺力模式下工件转速必须大于零。');
c = p.cutting;
validateattributes(c.depthOfCutMm,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(c.feedSpeedMmMin,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(c.diameterMm,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(c.enteringAngleDeg,{'numeric'},{'scalar','real','finite','>',0,'<=',90});
validateattributes(c.shearNPerMm2,{'numeric'},{'vector','numel',3,'real','finite','nonnegative'});
validateattributes(c.edgeNPerMm,{'numeric'},{'vector','numel',3,'real','finite','nonnegative'});
feed = c.feedSpeedMmMin/p.rpm;
area = c.depthOfCutMm*feed;
speed = pi*c.diameterMm*p.rpm/1000;
edge = c.depthOfCutMm/sind(c.enteringAngleDeg);
if area==0, edge = 0; end % 零切深或零进给用于本模型的空切近似。
force = c.shearNPerMm2(:)*area+c.edgeNPerMm(:)*edge;
process = struct('mode','mechanistic_force','rpm',p.rpm, ...
    'feedPerRevMm',feed,'chipAreaMm2',area,'cuttingSpeedMMin',speed, ...
    'engagedEdgeLengthMm',edge,'removalRateCm3Min',speed*area,'meanForceN',force);
end

function save_result_figure(fig,outdir,filename,showFigures)
% 同一组曲线导出中英文两版，桌面显示中文版。
[~,stemName,ext] = fileparts(filename);
apply_plot_language(fig,2,'Times New Roman');
print(fig,fullfile(outdir,[stemName '_en' ext]),'-dsvg','-vector');
apply_plot_language(fig,1,'宋体');
print(fig,fullfile(outdir,filename),'-dsvg','-vector');
if showFigures
    set(fig,'Visible','on');
else
    close(fig);
end
end

function apply_plot_language(fig,language,fontName)
for ax = reshape(findall(fig,'Type','axes'),1,[])
    labels = getappdata(ax,'PlotLabels');
    ax.FontName = fontName;
    set([ax.Title; ax.XLabel; ax.YLabel],{'String'},labels(language,:).', ...
        'FontName',fontName,'FontWeight','normal');
end
for lg = reshape(findall(fig,'Type','legend'),1,[])
    labels = getappdata(lg,'PlotLabels');
    set(lg,'String',labels(language,:),'FontName',fontName);
end
end

function p = defaults(p)
d.fs=64000; d.duration=8; d.discard=2; d.rpm=90; d.seed=42;
d.massKg=[250 35 35 1400 160 4.5 0.8 3500 250 0.162 4.5];
d.zetaLink=0.02;
d.meanForceN=[600 350 1800];
d.orders=[1 6 12]; d.orderFraction=[0.08 0.04 0.02];
d.orderPhase=[0 0.3 0.8]; d.noiseFraction=0.04;
d.forceNoiseCornerHz=600; d.plotMaxHz=4000; d.frfPoints=3001;
d.showFigures=usejava('desktop');
fields=fieldnames(d);
for j=1:numel(fields)
    if ~isfield(p,fields{j}), p.(fields{j})=d.(fields{j}); end
end
if numel(p.massKg)==9
    % 兼容旧版 9 项质量：均分刀杆质量，并追加刀盒。
    p.massKg = p.massKg(:).';
    p.massKg(6) = p.massKg(6)/2;
    p.massKg(10:11) = [d.massKg(10),p.massKg(6)];
end
assert(numel(p.massKg)==11 && numel(p.meanForceN)==3);
assert(numel(p.orders)==numel(p.orderFraction));
assert(numel(p.orders)==numel(p.orderPhase));
end

function links = expand_legacy_links(links)
% 展开旧连接，保持 7–6 和 6–5 的静态串联等效刚度。
ends = sort(links(:,1:2),2);
insertLink = all(ends==[6 7],2);
holderLink = all(ends==[5 6],2);
assert(any(insertLink) && any(holderLink), ...
    'Legacy links must include 7–6 and 6–5. Otherwise supply the 11-body links explicitly.');
ki = links(insertLink,3:5);
kh = links(holderLink,3:5);
ni = size(ki,1); nh = size(kh,1);
links = [links(~(insertLink | holderLink),:); ...
    repmat([7 10],ni,1),2*ki; repmat([10 6],ni,1),2*ki; ...
    repmat([6 11],nh,1),4*kh; repmat([11 5],nh,1),(4/3)*kh];
end

function ids = dofs(i)
ids=3*(i-1)+(1:3);
end

function [K,C] = add_link(K,C,i,j,k,c)
a=dofs(i); K(a,a)=K(a,a)+diag(k); C(a,a)=C(a,a)+diag(c);
if j~=0
    b=dofs(j);
    K(b,b)=K(b,b)+diag(k); C(b,b)=C(b,b)+diag(c);
    K(a,b)=K(a,b)-diag(k); K(b,a)=K(b,a)-diag(k);
    C(a,b)=C(a,b)-diag(c); C(b,a)=C(b,a)-diag(c);
end
end

function [f,A] = amplitude_spectrum(x,fs)
% 周期 Hann 窗并作相干增益修正，输出单边幅度谱。
n=size(x,2); x=bsxfun(@minus,x,mean(x,2));
w=0.5-0.5*cos(2*pi*(0:n-1)/n);
y=fft(bsxfun(@times,x,w),[],2)/sum(w);
last=floor(n/2)+1; A=abs(y(:,1:last));
if rem(n,2)==0
    A(:,2:end-1)=2*A(:,2:end-1);
else
    A(:,2:end)=2*A(:,2:end);
end
f=(0:last-1)*fs/n;
end

function val = validate_integrator()
% 独立单自由度算例，以解析正弦稳态结果为参考。
m=2; k=2e5; c=100; f=30; w=2*pi*f; F=10;
dt=1/16000; t=0:dt:1;
X=F/(k-w*w*m+1i*w*c);
q=real(X); v=real(1i*w*X); a=real(-w*w*X);
aa=zeros(size(t)); aa(1)=a; beta=0.25; gamma=0.5;
D=m+gamma*dt*c+beta*dt^2*k;
for j=2:numel(t)
    qp=q+dt*v+dt^2*(0.5-beta)*a; vp=v+dt*(1-gamma)*a;
    a=(F*cos(w*t(j))-c*vp-k*qp)/D;
    q=qp+beta*dt^2*a; v=vp+gamma*dt*a; aa(j)=a;
end
use=t>=0.2;
fit=[cos(w*t(use))' sin(w*t(use))' ones(sum(use),1)]\aa(use)';
Anumeric=fit(1)-1i*fit(2); Aexact=-w*w*X;
val.relativePhasorError=abs(Anumeric-Aexact)/abs(Aexact);
fs=1000; tt=(0:1999)/fs;
[ff,sp]=amplitude_spectrum(3*cos(2*pi*50*tt),fs);
[~,j]=min(abs(ff-50));
val.spectrumAmplitudeError=abs(sp(j)-3);
end
