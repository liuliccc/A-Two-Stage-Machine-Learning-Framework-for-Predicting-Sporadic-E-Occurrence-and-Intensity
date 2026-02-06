%% plot_11stations_obs_vs_predProb_withPrecision_LOCALTIME.m
% Unified Font Style: Arial
clear; clc;
%% ===== paths/settings =====
dataDir = pwd;
csvPath   = fullfile(dataDir,"kp_f10720_25.csv");
modelFile = fullfile(dataDir,"models_twoStage_nextObs","Unified_twoStage_nextObs.mat");
outPng    = fullfile(dataDir,'Fig_11stations_DOYLT_LOCALTIME_obsOcc_left_predProb_right_withPrecision.png');
assert(isfile(csvPath),   "CSV not found: " + csvPath);
assert(isfile(modelFile), "Model not found: " + modelFile);

% ---- Graphic Settings (Unified) ----
fontName = 'Arial';
titleFS  = 14; % Title Font Size
labelFS  = 12; % Axis Label Font Size
tickFS   = 10; % Tick Font Size
textFS   = 10; % In-plot Text Font Size

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
yearMap("ATHENS")        = 2014;
yearMap("GUAM")          = 2014;
yearMap("JULIUSRUH")     = 2008;
minSamplesPerStation = 50;

pThresh = 0.95; 
doyEdges   = 0.5:1:365.5;    
ltEdges    = 0:1:24;         
doyCenters = 1:365;
ltCenters  = 0.5:1:23.5;

% ---- longitude for LOCAL TIME ----
stationLonDeg = containers.Map('KeyType','char','ValueType','double');
stationLonDeg('Canberra')     = 149.13;
stationLonDeg('Casey')        = 110.52;
stationLonDeg('Christchurch') = 172.63;
stationLonDeg('Hobart')       = 147.33;
stationLonDeg('Norfolk Is')   = 167.95;
stationLonDeg('Townsville')   = 146.82;
stationLonDeg('brisbane')     = 153.03;
stationLonDeg('darwin')       = 130.84;
stationLonDeg('ATHENS')       = 23.72;
stationLonDeg('GUAM')         = 144.86;
stationLonDeg('JULIUSRUH')    = 13.40;

%% ===== load model =====
M = load(modelFile, "net_occ","muX","sdX","stationNames","lookback","occThresh","maxGapHoursForSample");
net_occ = M.net_occ;
muX = M.muX; sdX = M.sdX;
stationNames_model_raw = string(M.stationNames);
lookback   = M.lookback;
occThresh  = M.occThresh;
maxGapHours = M.maxGapHoursForSample;
stationNames_model = strtrim(stationNames_model_raw);

%% ===== read Kp/F107 timetable =====
Tcsv = readtable(csvPath);
req = ["Year","DOY","LT","Kp","F107"];
assert(all(ismember(req, string(Tcsv.Properties.VariableNames))), "CSV missing cols");
dtCsv = datetime(Tcsv.Year,1,1) + days(Tcsv.DOY-1) + hours(Tcsv.LT);
dtCsv.TimeZone = "";
TTcsv = timetable(dtCsv, Tcsv.Kp, Tcsv.F107, 'VariableNames', {'Kp','F107'});
TTcsv = sortrows(TTcsv);
TTcsv.Properties.RowTimes.TimeZone = "";

%% ===== global precision counters =====
TP_all = 0; FP_all = 0; FN_all = 0; TN_all = 0;

%% ===== plot =====
nSt = numel(stationList);
figure('Color','w','Position',[50 50 1700 1100]);
t = tiledlayout(nSt,2,'TileSpacing','compact','Padding','loose');

for k = 1:nSt
    stFileName = stationList(k);     
    stModelName = stFileName;        
    if isKey(aliasMap, char(stFileName))
        stModelName = string(aliasMap(char(stFileName)));
    end
    matFile = fullfile(dataDir, stFileName + "_hourly.mat");
    if ~isfile(matFile)
        warning("Missing: %s", matFile); nexttile; axis off; nexttile; axis off; continue;
    end
    S = load(matFile);
    if ~(isfield(S,'t_hourly') && isfield(S,'foes_hourly'))
        nexttile; axis off; nexttile; axis off; continue;
    end
    
    t_time = S.t_hourly(:);
    y = double(S.foes_hourly(:));
    y(y<0) = NaN;
    if ~isdatetime(t_time), t_time = datetime(string(t_time)); end
    t_time.TimeZone = "";
    [t_time,idx] = sort(t_time); y = y(idx);
    good = isfinite(y) & ~isnat(t_time);
    t_time = t_time(good); y = y(good);
    
    if isKey(yearMap, stFileName), targetYear = yearMap(stFileName); else, targetYear = year(t_time(1)); end
    idxY = (year(t_time) == targetYear);
    t_time = t_time(idxY); y = y(idxY);
    
    if isempty(t_time) || numel(y) < (lookback+2) || numel(y) < minSamplesPerStation
        nexttile; axis off; nexttile; axis off; continue;
    end
    
    if ~isKey(stationLonDeg, char(stFileName)), lonDeg = 0; else, lonDeg = stationLonDeg(char(stFileName)); end
    
    %% ===== LEFT: Observed =====
    [~, DOYloc, LTloc] = utc2lt(t_time, lonDeg);
    isOcc = double(y >= occThresh);
    obsGrid = accum2d_max(DOYloc, LTloc, isOcc, doyEdges, ltEdges);
    
    ax1 = nexttile;
    imagesc(doyCenters, ltCenters, obsGrid'); axis(ax1,'xy');
    xlim(ax1,[1 365]); ylim(ax1,[0 24]);
    
    % --- Font Styling ax1 ---
    ax1.FontName = fontName;
    ax1.FontSize = tickFS;
    ax1.Box = 'on';
    ax1.LineWidth = 1;
    ax1.XColor = 'k'; ax1.YColor = 'k';
    
    if k==1
        title(ax1,'Observation','FontSize',titleFS,'FontWeight','bold','FontName',fontName);
    end
    text(ax1, 2, 22.5, sprintf('%s', stFileName), ...
        'FontWeight','bold','BackgroundColor','w','Margin',2, ...
        'FontName',fontName, 'FontSize', textFS);
        
    ylabel(ax1,'LT','FontSize',labelFS,'FontName',fontName);
    if k ~= nSt
        ax1.XTick = [];
        ax1.XLabel.String = '';
    else
        xlabel(ax1,'DOY','FontSize',labelFS,'FontName',fontName);
    end
    colormap(ax1, [0.25 0.25 0.75; 1.0 0.95 0.1]);
    caxis(ax1,[0 1]);
    
    %% ===== RIGHT: Predicted =====
    sIdx = find(stationNames_model == strtrim(stModelName), 1);
    if isempty(sIdx), nexttile; axis off; continue; end
    onehot = zeros(1, numel(stationNames_model)); onehot(sIdx)=1;
    
    TTbase = timetable(t_time, zeros(size(t_time)), 'VariableNames', {'dummy'});
    TTbase.Properties.RowTimes.TimeZone = "";
    TTall = synchronize(TTbase, TTcsv, "first", "previous");
    Kp_vals = double(TTall.Kp); F107_vals = double(TTall.F107);
    
    [tTarget, Pocc, YtrueOcc] = predict_nextObs_Pocc_andTruth( ...
        t_time, y, Kp_vals, F107_vals, onehot, net_occ, muX, sdX, lookback, occThresh, maxGapHours);
    
    if isempty(Pocc), nexttile; axis off; continue; end
    
    YpredOcc = (Pocc >= pThresh);
    TP = sum(YpredOcc==1 & YtrueOcc==1);
    FP = sum(YpredOcc==1 & YtrueOcc==0);
    FN = sum(YpredOcc==0 & YtrueOcc==1);
    TN = sum(YpredOcc==0 & YtrueOcc==0);
    TP_all = TP_all + TP; FP_all = FP_all + FP;
    FN_all = FN_all + FN; TN_all = TN_all + TN;
    
    [~, DOYploc, LTploc] = utc2lt(tTarget, lonDeg);
    probGrid = accum2d_mean(DOYploc, LTploc, Pocc, doyEdges, ltEdges);
    
    ax2 = nexttile;
    imagesc(doyCenters, ltCenters, probGrid'); axis(ax2,'xy');
    xlim(ax2,[1 365]); ylim(ax2,[0 24]);
    
    % --- Font Styling ax2 ---
    ax2.FontName = fontName;
    ax2.FontSize = tickFS;
    ax2.Box = 'on';
    ax2.LineWidth = 1;
    ax2.XColor = 'k'; ax2.YColor = 'k';
    
    if k==1
        title(ax2,'One-stage','FontSize',titleFS,'FontWeight','bold','FontName',fontName);
    end
    ax2.YTick = []; ax2.YLabel.String = '';
    if k ~= nSt
        ax2.XTick = [];
        ax2.XLabel.String = '';
    else
        xlabel(ax2,'DOY','FontSize',labelFS,'FontName',fontName);
    end
    colormap(ax2, parula);
    caxis(ax2,[0 1]);
end

% Colorbar Styling
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Occurrence probability';
cb.FontName = fontName;
cb.FontSize = tickFS;
cb.Label.FontSize = labelFS; 
cb.Limits = [0 1];

exportgraphics(gcf, outPng, 'Resolution', 600);
fprintf("Saved: %s\n", outPng);

%% ===== helper functions (Unchanged logic, just compacted) =====
function [tLocal, DOY, LT] = utc2lt(tUtc, lonDeg)
    tUtc = tUtc(:); offsetHours = lonDeg / 15.0; tLocal = tUtc + hours(offsetHours);
    DOY = day(tLocal,'dayofyear'); LT  = mod(hour(tLocal) + minute(tLocal)/60 + second(tLocal)/3600, 24);
end
function G = accum2d_max(DOY, LT, V, doyEdges, ltEdges)
    V = double(V); nbx = numel(doyEdges)-1; nby = numel(ltEdges)-1;
    ix = discretize(DOY, doyEdges); iy = discretize(LT, ltEdges);
    ok = isfinite(ix) & isfinite(iy) & isfinite(V); G = zeros(nbx,nby);
    if ~any(ok), return; end
    lin = sub2ind([nbx,nby], ix(ok), iy(ok));
    G(:) = accumarray(lin, V(ok), [nbx*nby,1], @max, 0.0);
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
function [tTarget, Pocc, YtrueOcc] = predict_nextObs_Pocc_andTruth(t, y, Kp, F107, onehot, net_occ, muX, sdX, lookback, occThresh, maxGapHours)
    t = t(:); y = double(y(:)); Kp = double(Kp(:)); F107 = double(F107(:));
    DOY = day(t,'dayofyear'); LT = hour(t)+minute(t)/60+second(t)/3600;
    doySin = sin(2*pi*double(DOY)/365); doyCos = cos(2*pi*double(DOY)/365);
    ltSin = sin(2*pi*double(LT)/24); ltCos = cos(2*pi*double(LT)/24);
    n = numel(y); X = []; tTarget = NaT(0,1); YtrueOcc = []; m=0;
    for i = lookback:n-1
        hist = y(i-lookback+1:i); y1 = y(i+1);
        if any(~isfinite(hist)) || ~isfinite(y1) || any(~isfinite([Kp(i),F107(i)])), continue; end
        if hours(t(i+1)-t(i)) > maxGapHours, continue; end
        exo = [doySin(i),doyCos(i),ltSin(i),ltCos(i), Kp(i), F107(i)];
        X(end+1,:) = [hist.' hours(t(i)-t(i-lookback+1)) minutes(t(i)-t(i-1)) minutes(t(i+1)-t(i)) exo onehot]; %#ok<AGROW>
        tTarget(end+1,1) = t(i+1); YtrueOcc(end+1,1) = double(y1 >= occThresh); %#ok<AGROW>
    end
    if isempty(X), Pocc=[]; return; end
    Xz = (X - muX) ./ sdX; Pocc = predict(net_occ, double(Xz)); Pocc = max(min(double(Pocc),1),0);
end