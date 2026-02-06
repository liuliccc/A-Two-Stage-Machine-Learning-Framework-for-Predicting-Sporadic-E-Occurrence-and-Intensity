%% plot_11stations_obs_oneStage_twoStage_intensity_eventOnly_LOCALTIME.m
% Unified Font Style: Arial
% 3-column figure: Observed | One-stage | Two-stage
clear; clc;

%% ===== paths/settings =====
dataDir = pwd;
csvPath    = fullfile(dataDir,"kp_f10720_25.csv");
model1File = fullfile(dataDir,"models_oneStage_nextObs","Unified_oneStage_nextObs.mat");
model2File = fullfile(dataDir,"models_twoStage_nextObs","Unified_twoStage_nextObs.mat");
outPng = fullfile(dataDir,'Fig_11stations_DOYLT_LOCALTIME_obs_oneStage_twoStage_intensity_eventOnly.png');

assert(isfile(csvPath),    "CSV not found: " + csvPath);
assert(isfile(model1File), "One-stage model not found: " + model1File);
assert(isfile(model2File), "Two-stage model not found: " + model2File);

% ---- Graphic Settings (Unified) ----
fontName = 'Arial';
titleFS  = 14; 
labelFS  = 12; 
tickFS   = 10; 
textFS   = 10; 

% ---- stations ----
stationList = ["Canberra","Casey","Christchurch","Hobart","Norfolk Is", ...
               "Townsville","brisbane","darwin","ATHENS","GUAM","JULIUSRUH"];
aliasMap = containers.Map('KeyType','char','ValueType','char');
aliasMap('ATHENS') = 'AT138';

yearMap = containers.Map;
yearMap("Canberra")      = 2008;
yearMap("Casey")         = 2008;
yearMap("Christchurch")  = 2008;
yearMap("Hobart")        = 2008;
yearMap("Norfolk Is")    = 2010;
yearMap("Townsville")    = 2008;
yearMap("brisbane")      = 2008;
yearMap("darwin")        = 2008;
yearMap("ATHENS")        = 2009;
yearMap("GUAM")          = 2015;
yearMap("JULIUSRUH")     = 2010;
minSamplesPerStation = 50;

eventThresh = 4.0;   
vmin = 4; vmax = 10; 
doyEdges   = 0.5:1:365.5;     
ltEdges    = 0:1:24;          
doyCenters = 1:365;
ltCenters  = 0.5:1:23.5;

% ---- longitude for LOCAL TIME ----
stationLonDeg = containers.Map('KeyType','char','ValueType','double');
stationLonDeg('Canberra')      = 149.13;
stationLonDeg('Casey')         = 110.52;
stationLonDeg('Christchurch')  = 172.63;
stationLonDeg('Hobart')        = 147.33;
stationLonDeg('Norfolk Is')    = 167.95;
stationLonDeg('Townsville')    = 146.82;
stationLonDeg('brisbane')      = 153.03;
stationLonDeg('darwin')        = 130.84;
stationLonDeg('ATHENS')        = 23.72;
stationLonDeg('GUAM')          = 144.86;
stationLonDeg('JULIUSRUH')     = 13.40;

%% ===== load one-stage model =====
M1 = load(model1File, "net","muX","sdX","muY","sdY","stationNames","lookback","maxGapHoursForSample");
net1 = M1.net;
muX1 = M1.muX; sdX1 = M1.sdX; muY1 = M1.muY; sdY1 = M1.sdY;
stationNames_1 = strtrim(string(M1.stationNames));
lookback1 = M1.lookback;
maxGap1 = M1.maxGapHoursForSample;

%% ===== load two-stage model =====
M2 = load(model2File, ...
    "net_occ","net_reg_all","net_reg_event", ...
    "muX","sdX","muY","sdY", ...
    "stationNames","lookback","occThresh","pThresh","useSoftMix","maxGapHoursForSample");
net_occ  = M2.net_occ; net2_all = M2.net_reg_all; net2_evt = M2.net_reg_event;
muX2 = M2.muX; sdX2 = M2.sdX; muY2 = M2.muY; sdY2 = M2.sdY;
stationNames_2 = strtrim(string(M2.stationNames));
lookback2 = M2.lookback; occThresh = M2.occThresh; pThresh = M2.pThresh; useSoftMix= M2.useSoftMix;
maxGap2   = M2.maxGapHoursForSample;

%% ===== read Kp/F107 =====
Tcsv = readtable(csvPath);
dtCsv = datetime(Tcsv.Year,1,1) + days(Tcsv.DOY-1) + hours(Tcsv.LT); dtCsv.TimeZone = "";
TTcsv = timetable(dtCsv, Tcsv.Kp, Tcsv.F107, 'VariableNames', {'Kp','F107'});
TTcsv = sortrows(TTcsv); TTcsv.Properties.RowTimes.TimeZone = "";

%% ===== metrics accumulators =====
ae1_all=[]; se1_all=[]; re1_all=[];
ae2_all=[]; se2_all=[]; re2_all=[];
nEvent_all = 0;

%% ===== plot (3 columns) =====
nSt = numel(stationList);
figure('Color','w','Position',[30 30 2100 1100]);
tiledlayout(nSt,3,'TileSpacing','compact','Padding','compact');

for k = 1:nSt
    stFile = stationList(k);      
    stModel = stFile;             
    if isKey(aliasMap, char(stFile)), stModel = string(aliasMap(char(stFile))); end
    
    matFile = fullfile(dataDir, stFile + "_hourly.mat");
    if ~isfile(matFile), nexttile; axis off; nexttile; axis off; nexttile; axis off; continue; end
    
    S = load(matFile);
    if ~(isfield(S,'t_hourly') && isfield(S,'foes_hourly')), nexttile; axis off; nexttile; axis off; nexttile; axis off; continue; end
    
    t = S.t_hourly(:); y = double(S.foes_hourly(:)); y(y<0) = NaN;
    if ~isdatetime(t), t = datetime(string(t)); end
    t.TimeZone = ""; [t,idx] = sort(t); y = y(idx);
    good = isfinite(y) & ~isnat(t); t = t(good); y = y(good);
    
    if isempty(t), nexttile; axis off; nexttile; axis off; nexttile; axis off; continue; end
    
    if isKey(yearMap, stFile), targetYear = yearMap(stFile); else, targetYear = year(t(1)); end
    idxY = (year(t) == targetYear); t = t(idxY); y = y(idxY);
    
    if isempty(t) || numel(y) < (max(lookback1,lookback2)+2)
        nexttile; axis off; nexttile; axis off; nexttile; axis off; continue;
    end
    
    if ~isKey(stationLonDeg, char(stFile)), lonDeg = 0; else, lonDeg = stationLonDeg(char(stFile)); end
    
    s1 = find(stationNames_1 == strtrim(stModel), 1);
    s2 = find(stationNames_2 == strtrim(stModel), 1);
    if isempty(s1) || isempty(s2), nexttile; axis off; nexttile; axis off; nexttile; axis off; continue; end
    onehot1 = zeros(1,numel(stationNames_1)); onehot1(s1)=1;
    onehot2 = zeros(1,numel(stationNames_2)); onehot2(s2)=1;
    
    TTbase = timetable(t, zeros(size(t)), 'VariableNames', {'dummy'}); TTbase.Properties.RowTimes.TimeZone = "";
    TTall = synchronize(TTbase, TTcsv, "first", "previous");
    Kp = double(TTall.Kp); F107 = double(TTall.F107);
    
    [tT1, ypred1, ytrue1] = predict_nextObs_intensity_oneStage(t, y, Kp, F107, onehot1, net1, muX1, sdX1, muY1, sdY1, lookback1, maxGap1);
    [tT2, ypred2, ytrue2] = predict_nextObs_intensity_twoStage(t, y, Kp, F107, onehot2, net_occ, net2_all, net2_evt, muX2, sdX2, muY2, sdY2, lookback2, occThresh, maxGap2, pThresh, useSoftMix);
    
    [tCommon, ia, ib] = intersect(tT1, tT2);
    if isempty(tCommon), nexttile; axis off; nexttile; axis off; nexttile; axis off; continue; end
    
    ytrue = ytrue1(ia); y1p = ypred1(ia); y2p = ypred2(ib);
    
    ev = (ytrue >= eventThresh); nEv = sum(ev);
    if nEv > 0
        yT = ytrue(ev); y1 = y1p(ev); y2 = y2p(ev);
        ae1 = abs(y1 - yT); se1 = (y1 - yT).^2;
        ae2 = abs(y2 - yT); se2 = (y2 - yT).^2;
        okDen = isfinite(yT) & (yT > 0);
        re1 = abs(y1(okDen) - yT(okDen)) ./ yT(okDen);
        re2 = abs(y2(okDen) - yT(okDen)) ./ yT(okDen);
        ae1_all=[ae1_all; ae1(:)]; se1_all=[se1_all; se1(:)];
        ae2_all=[ae2_all; ae2(:)]; se2_all=[se2_all; se2(:)];
        re1_all=[re1_all; re1(:)]; re2_all=[re2_all; re2(:)];
        nEvent_all = nEvent_all + nEv;
    end
    
    [~, DOYloc, LTloc] = utc2lt(tCommon, lonDeg);
    obsGrid = accum2d_mean(DOYloc, LTloc, ytrue, doyEdges, ltEdges);
    oneGrid = accum2d_mean(DOYloc, LTloc, y1p,   doyEdges, ltEdges);
    twoGrid = accum2d_mean(DOYloc, LTloc, y2p,   doyEdges, ltEdges);
    
    % ---- Col 1: Observed ----
    ax1 = nexttile;
    imagesc(ax1, doyCenters, ltCenters, obsGrid'); axis(ax1,'xy');
    xlim(ax1,[1 365]); ylim(ax1,[0 24]);
    
    % Style ax1
    ax1.FontName = fontName;
    ax1.FontSize = tickFS;
    ax1.Box = 'on'; ax1.LineWidth=1; ax1.XColor='k'; ax1.YColor='k';
    
    ylabel(ax1, 'LT', 'FontSize', labelFS, 'FontName', fontName);
    if k==1
        title(ax1, 'Observation', 'FontSize', titleFS, 'FontWeight', 'bold', 'FontName', fontName);
    end
    text(ax1, 2, 22.5, sprintf('%s', stFile), ...
        'FontWeight','bold','BackgroundColor','w','Margin',2, ...
        'FontName', fontName, 'FontSize', textFS);
        
    if k ~= nSt
        ax1.XTick = []; ax1.XLabel.String = '';
    else
        xlabel(ax1, 'DOY', 'FontSize', labelFS, 'FontName', fontName);
    end
    colormap(ax1, parula); caxis(ax1,[vmin vmax]);
    
    % ---- Col 2: One-stage ----
    ax2 = nexttile;
    imagesc(ax2, doyCenters, ltCenters, oneGrid'); axis(ax2,'xy');
    xlim(ax2,[1 365]); ylim(ax2,[0 24]);
    
    % Style ax2
    ax2.FontName = fontName;
    ax2.FontSize = tickFS;
    ax2.Box = 'on'; ax2.LineWidth=1; ax2.XColor='k'; ax2.YColor='k';
    
    if k==1
        title(ax2, 'Single-baseline', 'FontSize', titleFS, 'FontWeight', 'bold', 'FontName', fontName);
    end
    ax2.YTick = []; 
    if k ~= nSt
        ax2.XTick = []; ax2.XLabel.String = '';
    else
        xlabel(ax2, 'DOY', 'FontSize', labelFS, 'FontName', fontName);
    end
    colormap(ax2, parula); caxis(ax2,[vmin vmax]);
    
    % ---- Col 3: Two-stage ----
    ax3 = nexttile;
    imagesc(ax3, doyCenters, ltCenters, twoGrid'); axis(ax3,'xy');
    xlim(ax3,[1 365]); ylim(ax3,[0 24]);
    
    % Style ax3
    ax3.FontName = fontName;
    ax3.FontSize = tickFS;
    ax3.Box = 'on'; ax3.LineWidth=1; ax3.XColor='k'; ax3.YColor='k';
    
    if k==1
        title(ax3, 'Two-stage', 'FontSize', titleFS, 'FontWeight', 'bold', 'FontName', fontName);
    end
    ax3.YTick = []; 
    if k ~= nSt
        ax3.XTick = []; ax3.XLabel.String = '';
    else
        xlabel(ax3, 'DOY', 'FontSize', labelFS, 'FontName', fontName);
    end
    colormap(ax3, parula); caxis(ax3,[vmin vmax]);
end

% Shared Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'foEs (MHz)';
cb.Limits = [vmin vmax];
cb.FontName = fontName;
cb.FontSize = tickFS;
cb.Label.FontSize = labelFS;

exportgraphics(gcf, outPng, 'Resolution', 600);
fprintf("Saved: %s\n", outPng);

%% ===== overall metrics (EVENT-ONLY) =====
MAE1  = mean(ae1_all,"omitnan"); RMSE1 = sqrt(mean(se1_all,"omitnan"));
MAE2  = mean(ae2_all,"omitnan"); RMSE2 = sqrt(mean(se2_all,"omitnan"));
MRE1  = mean(re1_all,"omitnan"); MRE2  = mean(re2_all,"omitnan");
fprintf("\nOVERALL EVENT-ONLY (Ytrue>=%.1f), Nevent=%d\n", eventThresh, nEvent_all);
fprintf("Baseline:  MAE=%.4f  MRE=%.4f  RMSE=%.4f\n", MAE1, MRE1, RMSE1);
fprintf("Two-stage: MAE=%.4f  MRE=%.4f  RMSE=%.4f\n", MAE2, MRE2, RMSE2);

%% ===== helper functions =====
function [tLocal, DOY, LT] = utc2lt(tUtc, lonDeg)
    tUtc = tUtc(:); offsetHours = lonDeg / 15.0; tLocal = tUtc + hours(offsetHours);
    DOY = day(tLocal,'dayofyear'); LT = mod(hour(tLocal) + minute(tLocal)/60 + second(tLocal)/3600, 24);
end
function G = accum2d_mean(DOY, LT, V, doyEdges, ltEdges)
    V = double(V); nbx = numel(doyEdges)-1; nby = numel(ltEdges)-1;
    ix = discretize(DOY, doyEdges); iy = discretize(LT, ltEdges);
    ok = isfinite(ix) & isfinite(iy) & isfinite(V); G = nan(nbx,nby);
    if ~any(ok), return; end
    lin  = sub2ind([nbx,nby], ix(ok), iy(ok));
    sumV = accumarray(lin, V(ok), [nbx*nby,1], @sum, 0.0);
    cnt  = accumarray(lin, 1,    [nbx*nby,1], @sum, 0);
    m = sumV ./ max(cnt,1); m(cnt==0) = NaN; G(:) = m;
end
function [tTarget, Ypred, Ytrue] = predict_nextObs_intensity_oneStage(t, y, Kp, F107, onehot, net, muX, sdX, muY, sdY, lookback, maxGapHours)
    t=t(:); y=double(y(:)); Kp=double(Kp(:)); F107=double(F107(:));
    DOY = day(t,'dayofyear'); LT = hour(t)+minute(t)/60+second(t)/3600;
    doySin = sin(2*pi*double(DOY)/365); doyCos = cos(2*pi*double(DOY)/365);
    ltSin = sin(2*pi*double(LT)/24); ltCos = cos(2*pi*double(LT)/24);
    n=numel(y); X=[]; tTarget=NaT(0,1); Ytrue=[];
    for i=lookback:n-1
        hist=y(i-lookback+1:i); y1=y(i+1);
        if any(~isfinite(hist)) || ~isfinite(y1) || any(~isfinite([Kp(i),F107(i)])), continue; end
        if hours(t(i+1)-t(i)) > maxGapHours, continue; end
        exo=[doySin(i),doyCos(i),ltSin(i),ltCos(i),Kp(i),F107(i)];
        X(end+1,:)=[hist.' hours(t(i)-t(i-lookback+1)) minutes(t(i)-t(i-1)) minutes(t(i+1)-t(i)) exo onehot]; %#ok<AGROW>
        tTarget(end+1,1)=t(i+1); Ytrue(end+1,1)=y1; %#ok<AGROW>
    end
    if isempty(X), Ypred=[]; return; end
    Xz=(X-muX)./sdX; Ypred = max(expm1(double(predict(net, double(Xz))).*sdY + muY), 0);
end
function [tTarget, Ypred, Ytrue] = predict_nextObs_intensity_twoStage(t, y, Kp, F107, onehot, net_occ, net_reg_all, net_reg_event, muX, sdX, muY, sdY, lookback, occThresh, maxGapHours, pThresh, useSoftMix)
    t=t(:); y=double(y(:)); Kp=double(Kp(:)); F107=double(F107(:));
    DOY = day(t,'dayofyear'); LT = hour(t)+minute(t)/60+second(t)/3600;
    doySin = sin(2*pi*double(DOY)/365); doyCos = cos(2*pi*double(DOY)/365);
    ltSin = sin(2*pi*double(LT)/24); ltCos = cos(2*pi*double(LT)/24);
    n=numel(y); X=[]; tTarget=NaT(0,1); Ytrue=[];
    for i=lookback:n-1
        hist=y(i-lookback+1:i); y1=y(i+1);
        if any(~isfinite(hist)) || ~isfinite(y1) || any(~isfinite([Kp(i),F107(i)])), continue; end
        if hours(t(i+1)-t(i)) > maxGapHours, continue; end
        exo=[doySin(i),doyCos(i),ltSin(i),ltCos(i),Kp(i),F107(i)];
        X(end+1,:)=[hist.' hours(t(i)-t(i-lookback+1)) minutes(t(i)-t(i-1)) minutes(t(i+1)-t(i)) exo onehot]; %#ok<AGROW>
        tTarget(end+1,1)=t(i+1); Ytrue(end+1,1)=y1; %#ok<AGROW>
    end
    if isempty(X), Ypred=[]; return; end
    Xz=(X-muX)./sdX;
    Pocc = max(min(double(predict(net_occ, double(Xz))),1),0);
    Yall = max(expm1(double(predict(net_reg_all, double(Xz))).*sdY + muY),0);
    Yevt = max(expm1(double(predict(net_reg_event, double(Xz))).*sdY + muY),0);
    if useSoftMix, Ypred = Pocc .* Yevt + (1-Pocc) .* Yall; else, Ypred = Yall; sel = (Pocc >= pThresh); Ypred(sel) = Yevt(sel); end
    Ypred = max(Ypred,0);
end