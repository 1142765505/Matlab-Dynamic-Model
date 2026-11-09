function result = wheel_lathe_lumped_demo(outdir, p)
% WHEEL_LATHE_LUMPED_DEMO Pure MATLAB wheel-lathe vibration demonstration.
% Run: result = wheel_lathe_lumped_demo;
% Custom: p.fs=32000; p.duration=8; p.rpm=90;
%         result = wheel_lathe_lumped_demo('my_result',p);
% Set p.showFigures=false to export PNGs without opening figure windows.
% Base MATLAB only. SI units: kg, N, m, s.
% ALL DEFAULT PHYSICAL PARAMETERS AND FORCE FLUCTUATIONS ARE ILLUSTRATIVE.
% This is a linear 27-DOF lumped model, not a calibrated machine prediction.
% Each body has 3 translations: radial, axial, tangential. No rotations,
% contact opening, regenerative chatter or cutting-chip FE model is included.
% Static preload effects must be represented by identified tangent stiffnesses.

if nargin < 1 || isempty(outdir)
    outdir = fullfile(fileparts(mfilename('fullpath')), 'matlab_result');
end
if nargin < 2, p = struct; end
p = defaults(p);
if ~exist(outdir,'dir'), mkdir(outdir); end
assert(p.duration > p.discard && p.discard >= 0);
assert(p.fs > 2*p.forceNoiseCornerHz && p.fs > 2*p.plotMaxHz);
assert(p.zetaLink > 0 && all(p.massKg > 0));
rng(p.seed,'twister');

% Bodies: left spindle/faceplate, left centre, right centre, wheelset,
% carriage, holder, tip-equivalent mass, bed, right spindle/faceplate.
names = {'left_spindle_plate','left_centre','right_centre','wheelset', ...
         'carriage','holder','tip','bed','right_spindle_plate'};
nbody = numel(names); nd = 3*nbody;
M = diag(kron(p.massKg(:),ones(3,1)));
K = zeros(nd); C = zeros(nd);
% Row [i,j,k_radial,k_axial,k_tangential]; j=0 means ground.
% Connections to centres are tangent equivalents under the operating preload.
links = [1 8 2e8 3e8 2e8; ...
         2 8 1e8 1.2e8 1e8; ...
         3 8 1e8 1.2e8 1e8; ...
         4 2 1.4e8 1e8 1.4e8; ...
         4 3 1.4e8 1e8 1.4e8; ...
         4 1 8e6 3e6 3e7; ...
         5 8 1.5e8 1.2e8 1e8; ...
         6 5 3e7 2.5e7 2e7; ...
         7 6 1.2e7 9e6 1.5e7; ...
         8 0 8e8 8e8 6e8; ...
         9 8 2e8 3e8 2e8; ...
         4 9 8e6 3e6 3e7];
if isfield(p,'links'), links = p.links; end
for j = 1:size(links,1)
    i1 = links(j,1); i2 = links(j,2); k = links(j,3:5);
    assert(all(k > 0));
    meq = p.massKg(i1);
    if i2 ~= 0
        meq = p.massKg(i1)*p.massKg(i2)/(p.massKg(i1)+p.massKg(i2));
    end
    % Illustrative link damping based on equivalent two-body mass.
    % zetaLink is NOT the damping ratio of every assembled mode.
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

tip = dofs(7); holder = dofs(6); wheel = dofs(4);
B = zeros(nd,3); B(tip,:) = eye(3); B(wheel,:) = -eye(3);
assert(norm(reshape(sum(reshape(B,3,nbody,3),2),3,3),'fro') < 1e-12);
% Positive vector is the force on the tool; opposite force acts on wheelset.
dt = 1/p.fs; t = 0:dt:p.duration; N = numel(t);
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
force = bsxfun(@times,p.meanForceN(:), ...
        bsxfun(@plus,1+modulation,p.noiseFraction*eta));
% No intentional chatter: this force is prescribed, independent of q.

% Newmark average-acceleration integration (beta=1/4, gamma=1/2).
% Constant linear matrices allow one Cholesky factorization for all steps.
beta = 1/4; gamma = 1/2;
L = chol(M+gamma*dt*C+beta*dt^2*K,'lower');
q = K\(B*force(:,1)); v = zeros(nd,1); acc = zeros(nd,1);
tipA = zeros(3,N); holderA = zeros(3,N); relativeU = zeros(3,N);
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
    relativeU(:,j) = q(tip)-q(wheel);
    if mod(j,1000)==0
        maxResidual = max(maxResidual,norm(M*acc+C*v+K*q-f,inf));
        maxRhs = max(maxRhs,norm(f,inf));
    end
end
elapsed = toc;
equilibriumResidual = maxResidual/max(maxRhs,1);
assert(equilibriumResidual < 1e-6,'Dynamic equilibrium residual too large.');

% Unit-force accelerance: output direction i, input direction j.
freq = linspace(0,p.plotMaxHz,p.frfPoints);
Htip = complex(zeros(3,3,numel(freq)));
for j = 1:numel(freq)
    w = 2*pi*freq(j);
    U = (K-w^2*M+1i*w*C)\B;
    Htip(:,:,j) = -w^2*U(tip,:);
end
keep = t >= p.discard;
[freqFFT,ampFFT] = amplitude_spectrum(tipA(:,keep),p.fs);
rmsTip = sqrt(mean(tipA(:,keep).^2,2));
peakTip = max(abs(tipA(:,keep)),[],2);
axesNames = {'radial';'axial';'tangential'};
metrics = table(axesNames,rmsTip,peakTip,rmsTip/9.80665, ...
                'VariableNames',{'direction','RMS_m_s2','Peak_m_s2','RMS_g'});
disp(metrics);
fprintf('Integration %.2f s; normalized equilibrium residual %.3g\n', ...
        elapsed,equilibriumResidual);

% Numerical checks independent of machine calibration.
validation = validate_integrator();
assert(validation.relativePhasorError < 0.01,'Harmonic validation failed.');
assert(validation.spectrumAmplitudeError < 1e-10,'FFT scaling check failed.');
fprintf('Analytic harmonic check error %.4g; FFT amplitude error %.3g\n', ...
        validation.relativePhasorError,validation.spectrumAmplitudeError);

result = struct('parameters',p,'bodyNames',{names},'links',links, ...
    'M',M,'C',C,'K',K,'Bcut',B,'naturalFrequencyHz',fn,'time',t, ...
    'cuttingForceN',force,'tipAcceleration',tipA, ...
    'holderAcceleration',holderA,'relativeDisplacement',relativeU, ...
    'frequencyHz',freq,'tipAccelerance',Htip, ...
    'fftFrequencyHz',freqFFT,'tipAmplitudeSpectrum',ampFFT, ...
    'metrics',metrics,'validation',validation, ...
    'equilibriumResidual',equilibriumResidual);
save(fullfile(outdir,'simulation.mat'),'result','-v7');
writetable(metrics,fullfile(outdir,'metrics.csv'));
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
fig = figure('Name','刀尖加速度','NumberTitle','off', ...
    'Visible','off','Color','w','Position',[50 50 1100 850]);
for j = 1:3
    subplot(3,2,2*j-1);
    plot(t(short),tipA(j,short),'LineWidth',0.7); grid on;
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '刀尖加速度时程'], '时间 (s)', '加速度 (m/s^2)'; ...
         [axesNames{j} ' tool-tip acceleration'], 'Time (s)', 'Acceleration (m/s^2)'});
    subplot(3,2,2*j);
    plot(freqFFT,ampFFT(j,:),'LineWidth',0.7); grid on;
    xlim([0 p.plotMaxHz]);
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '刀尖加速度幅度谱'], '频率 (Hz)', '单边幅值 (m/s^2)'; ...
         [axesNames{j} ' acceleration spectrum'], 'Frequency (Hz)', 'Single-sided amplitude (m/s^2)'});
end
save_result_figure(fig,outdir,'acceleration.png',p.showFigures);

fig = figure('Name','单位力加速度频响','NumberTitle','off', ...
    'Visible','off','Color','w','Position',[50 50 950 650]);
for j = 1:3
    subplot(3,1,j);
    semilogy(freq,max(squeeze(abs(Htip(j,j,:))),1e-16)); grid on;
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '单位力加速度频响'], '频率 (Hz)', '加速度频响幅值 ((m/s^2)/N)'; ...
         [axesNames{j} ' unit-force accelerance'], 'Frequency (Hz)', 'Accelerance magnitude ((m/s^2)/N)'});
end
save_result_figure(fig,outdir,'unit_force_frf.png',p.showFigures);

% A longer time window exposes the low-frequency cutting-force modulation.
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
save_result_figure(fig,outdir,'cutting_force.png',p.showFigures);

fig = figure('Name','刀具与轮对响应','NumberTitle','off', ...
    'Visible','off','Color','w','Position',[50 50 1100 850]);
for j = 1:3
    subplot(3,2,2*j-1);
    plot(t(short),tipA(j,short),t(short),holderA(j,short),'LineWidth',0.7); grid on;
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '刀尖与刀杆加速度对比'], '时间 (s)', '加速度 (m/s^2)'; ...
         [axesNames{j} ' tool-tip and holder acceleration'], 'Time (s)', 'Acceleration (m/s^2)'});
    if j == 1
        lg = legend('刀尖','刀杆','Location','best');
        setappdata(lg,'PlotLabels',{'刀尖','刀杆'; 'Tool tip','Tool holder'});
    end
    subplot(3,2,2*j);
    plot(t(short),1e6*relativeU(j,short),'LineWidth',0.7); grid on;
    setappdata(gca,'PlotLabels', ...
        {[directionTitles{j} '刀尖相对轮对的位移'], '时间 (s)', '相对位移 ({\mu}m)'; ...
         [axesNames{j} ' displacement: tip minus wheelset'], 'Time (s)', 'Relative displacement ({\mu}m)'});
end
save_result_figure(fig,outdir,'tool_response.png',p.showFigures);

fig = figure('Name','未阻尼固有频率','NumberTitle','off', ...
    'Visible','off','Color','w','Position',[50 50 950 500]);
stem(1:nd,fn,'filled'); grid on;
xlim([0.5 nd+0.5]); xticks(1:nd);
setappdata(gca,'PlotLabels', ...
    {'未阻尼固有频率', '模态阶次', '固有频率 (Hz)'; ...
     'Undamped natural frequencies', 'Mode index', 'Natural frequency (Hz)'});
save_result_figure(fig,outdir,'natural_frequencies.png',p.showFigures);
fprintf('Saved illustrative results to: %s\n',outdir);
end

function save_result_figure(fig,outdir,filename,showFigures)
% Reuse the same plotted data. Leave desktop figures in Chinese.
[~,stemName,ext] = fileparts(filename);
apply_plot_language(fig,2,'Times New Roman');
print(fig,fullfile(outdir,[stemName '_en' ext]),'-dpng','-r140');
apply_plot_language(fig,1,'宋体');
print(fig,fullfile(outdir,filename),'-dpng','-r140');
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
d.fs=32000; d.duration=8; d.discard=2; d.rpm=90; d.seed=42;
d.massKg=[250 35 35 1400 160 9 0.8 3500 250];
d.zetaLink=0.02;
d.meanForceN=[600 350 1800];
d.orders=[1 6 12]; d.orderFraction=[0.08 0.04 0.02];
d.orderPhase=[0 0.3 0.8]; d.noiseFraction=0.04;
d.forceNoiseCornerHz=600; d.plotMaxHz=1500; d.frfPoints=3001;
d.showFigures=usejava('desktop');
fields=fieldnames(d);
for j=1:numel(fields)
    if ~isfield(p,fields{j}), p.(fields{j})=d.(fields{j}); end
end
assert(numel(p.massKg)==9 && numel(p.meanForceN)==3);
assert(numel(p.orders)==numel(p.orderFraction));
assert(numel(p.orders)==numel(p.orderPhase));
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
% Periodic Hann window, coherent-gain corrected; no Signal Toolbox.
% For broadband noise this is an amplitude spectrum, NOT a PSD estimate.
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
% Independent scalar oscillator, analytic sinusoidal steady-state initial data.
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
