%% Two-stage model 
clear; clc;
try, gpuDevice([]); catch, end

%% ===== settings =====
dataDir = pwd;
csvPath = fullfile(dataDir,"kp_f10720_25.csv");
outDir  = fullfile(dataDir,"models_twoStage_nextObs");
if ~exist(outDir,"dir"), mkdir(outDir); end

lookback = 12;
trainRatio = 0.80;
occThresh  = 4.0;
pThresh    = 0.8;
useSoftMix = true;

% gap filtering (avoid crazy gaps)
maxGapHoursForSample = 12;   % if next gap too large, drop that sample (tune 6~24)
% network
hidden1=42; hidden2=42;hidden3=42; dropout=0.10;
hidden4=52; hidden5=52;hidden6=52; dropout=0.10;
maxEpochs_occ=40; maxEpochs_reg=60;
miniBatch=512; learnRate=1e-3;
showTrainingPlot=true;
posFactor=3;

%% ===== read Kp/F107 =====
assert(isfile(csvPath),"CSV not found: "+csvPath);
Tcsv = readtable(csvPath);
req = ["Year","DOY","LT","Kp","F107"];
assert(all(ismember(req, string(Tcsv.Properties.VariableNames))), ...
    "CSV must contain Year,DOY,LT,Kp,F107");

dtCsv = datetime(Tcsv.Year,1,1) + days(Tcsv.DOY-1) + hours(Tcsv.LT);
dtCsv.TimeZone = "";
TTcsv = timetable(dtCsv, Tcsv.Kp, Tcsv.F107, 'VariableNames', {'Kp','F107'});
TTcsv = sortrows(TTcsv);
TTcsv.Properties.RowTimes.TimeZone = "";

%% ===== station files =====
files = dir(fullfile(dataDir,"*_hourly.mat"));
assert(~isempty(files), "No *_hourly.mat found.");

stationNames = string(erase({files.name},"_hourly.mat"));
numStations = numel(stationNames);
fprintf("Found %d stations.\n", numStations);

%% ===== collect train/test =====
Xtr_all=[]; Yocc_tr_all=[]; Yreg_tr_all=[];
Xte_all=[]; Yocc_te_all=[]; Yreg_te_all=[];

for s=1:numStations
    stName = stationNames(s);
    S = load(fullfile(files(s).folder, files(s).name));

    if ~(isfield(S,"t_hourly") && isfield(S,"foes_hourly"))
        fprintf("Skip %s (missing t_hourly/foes_hourly)\n", stName); 
        continue;
    end

    t = S.t_hourly(:);
    y = double(S.foes_hourly(:));
    y(y<0)=NaN;

    % timezone unify
    if isdatetime(t), t.TimeZone=""; else, t=datetime(string(t)); t.TimeZone=""; end

    % sort + drop NaN
    [t, idx] = sort(t);
    y = y(idx);
    good = isfinite(y) & ~isnat(t);
    t = t(good); y = y(good);

    if numel(y) < (lookback + 200)
        fprintf("Skip %s (too short after NaN drop: %d)\n", stName, numel(y));
        continue;
    end

    % DOY/LT at each observation
    DOY = day(t,'dayofyear');
    LT  = hour(t);
    doySin = sin(2*pi*double(DOY)/365);
    doyCos = cos(2*pi*double(DOY)/365);
    ltSin  = sin(2*pi*double(LT)/24);
    ltCos  = cos(2*pi*double(LT)/24);

    % Kp/F107 align to observation time (previous, no leakage)
    TTbase = timetable(t, zeros(size(t)), 'VariableNames', {'dummy'});
    TTbase.Properties.RowTimes.TimeZone = "";
    TTall = synchronize(TTbase, TTcsv, "first", "previous");
    Kp = TTall.Kp;
    F107 = TTall.F107;

    % Build irregular samples: last 12 obs -> next obs
    [X, Yreg, Yocc] = buildXY_nextObs(y, t, doySin,doyCos,ltSin,ltCos, Kp,F107, lookback, occThresh, maxGapHoursForSample);
    if isempty(X)
        fprintf("Skip %s (no valid samples)\n", stName);
        continue;
    end

    % station one-hot
    onehot = zeros(1,numStations); onehot(s)=1;
    X = [X, repmat(onehot, size(X,1), 1)];

    n = size(X,1);
    nTrain = floor(trainRatio*n);

    Xtr_all=[Xtr_all; X(1:nTrain,:)]; %#ok<AGROW>
    Yreg_tr_all=[Yreg_tr_all; Yreg(1:nTrain,:)]; %#ok<AGROW>
    Yocc_tr_all=[Yocc_tr_all; Yocc(1:nTrain,:)]; %#ok<AGROW>

    Xte_all=[Xte_all; X(nTrain+1:end,:)]; %#ok<AGROW>
    Yreg_te_all=[Yreg_te_all; Yreg(nTrain+1:end,:)]; %#ok<AGROW>
    Yocc_te_all=[Yocc_te_all; Yocc(nTrain+1:end,:)]; %#ok<AGROW>

    fprintf("%s: train=%d test=%d (event train=%d)\n", ...
        stName, nTrain, n-nTrain, sum(Yocc(1:nTrain)==1));
end

assert(~isempty(Xtr_all), "No training data collected.");
fprintf("\nTOTAL train=%d test=%d inputDim=%d\n", size(Xtr_all,1), size(Xte_all,1), size(Xtr_all,2));

%% ===== standardize X =====
muX = mean(Xtr_all,1,"omitnan");
sdX = std(Xtr_all,0,1,"omitnan");
sdX(sdX==0 | isnan(sdX)) = 1;

XtrZ = (Xtr_all - muX) ./ sdX;
XteZ = (Xte_all - muX) ./ sdX;
inputDim = size(XtrZ,2);

%% ===== prepare regression Y: log1p + zscore =====
Yreg_tr_log = log1p(max(Yreg_tr_all,0));
muY = mean(Yreg_tr_log,1,"omitnan");
sdY = std(Yreg_tr_log,0,1,"omitnan");
sdY(sdY==0 | isnan(sdY)) = 1;
Yreg_trZ = (Yreg_tr_log - muY) ./ sdY;

Yreg_te_log = log1p(max(Yreg_te_all,0));
Yreg_teZ = (Yreg_te_log - muY) ./ sdY;

%% ===== options =====
plotOpt="none"; if showTrainingPlot, plotOpt="training-progress"; end
opt_occ = trainingOptions("adam","ExecutionEnvironment","cpu", ...
    "MaxEpochs",maxEpochs_occ,"MiniBatchSize",miniBatch, ...
    "InitialLearnRate",learnRate,"Shuffle","every-epoch", ...
    "GradientThreshold",1,"Verbose",true,"Plots",plotOpt);

opt_reg = trainingOptions("adam","ExecutionEnvironment","cpu", ...
    "MaxEpochs",maxEpochs_reg,"MiniBatchSize",miniBatch, ...
    "InitialLearnRate",learnRate,"Shuffle","every-epoch", ...
    "GradientThreshold",1,"Verbose",true,"Plots",plotOpt);

%% ===== stage-1 occ =====
layers_occ = [
    featureInputLayer(inputDim,"Normalization","none")
    fullyConnectedLayer(hidden1); reluLayer; dropoutLayer(dropout)
    fullyConnectedLayer(hidden2); reluLayer; dropoutLayer(dropout)
    fullyConnectedLayer(hidden3); reluLayer; dropoutLayer(dropout)
    fullyConnectedLayer(1)
    sigmoidLayer
    regressionLayer
];

bad1 = any(~isfinite(XtrZ),2) | ~isfinite(Yocc_tr_all);
X1 = double(XtrZ(~bad1,:));
Y1 = double(Yocc_tr_all(~bad1,:));

pos = find(Y1==1);
if ~isempty(pos) && posFactor>1
    repIdx = repmat(pos, posFactor-1, 1);
    X1=[X1; X1(repIdx,:)]; %#ok<AGROW>
    Y1=[Y1; Y1(repIdx,:)]; %#ok<AGROW>
end
net_occ = trainNetwork(X1,Y1,layers_occ,opt_occ);

%% ===== reg layers =====
layers_reg = [
    featureInputLayer(inputDim,"Normalization","none")
    fullyConnectedLayer(hidden4); reluLayer; dropoutLayer(dropout)
    fullyConnectedLayer(hidden5); reluLayer; dropoutLayer(dropout)
    fullyConnectedLayer(hidden6); reluLayer; dropoutLayer(dropout)
    fullyConnectedLayer(1)
    regressionLayer
];

% reg all
badR = any(~isfinite(XtrZ),2) | any(~isfinite(Yreg_trZ),2);
Xr = double(XtrZ(~badR,:));
Yr = double(Yreg_trZ(~badR,:));
net_reg_all = trainNetwork(Xr,Yr,layers_reg,opt_reg);

% reg event-only
eventMask = (Yocc_tr_all==1);
Xev = XtrZ(eventMask,:);
Yev = Yreg_trZ(eventMask,:);
badE = any(~isfinite(Xev),2) | any(~isfinite(Yev),2);
Xev = double(Xev(~badE,:));
Yev = double(Yev(~badE,:));
net_reg_event = trainNetwork(Xev,Yev,layers_reg,opt_reg);

%% ===== eval =====
Pocc = predict(net_occ, double(XteZ));
Pocc = max(min(Pocc,1),0);

YallZ = predict(net_reg_all, double(XteZ));
YevtZ = predict(net_reg_event, double(XteZ));

Yall = max(expm1(YallZ.*sdY + muY),0);
Yevt = max(expm1(YevtZ.*sdY + muY),0);

Ytrue = Yreg_te_all;

if useSoftMix
    Ypred = Pocc .* Yevt + (1-Pocc) .* Yall;
else
    Ypred = Yall;
    sel = (Pocc >= pThresh);
    Ypred(sel) = Yevt(sel);
end

mae  = mean(abs(Ypred-Ytrue),"omitnan");
rmse = sqrt(mean((Ypred-Ytrue).^2,"omitnan"));
fprintf("\nTwo-stage nextObs: MAE=%.4f  RMSE=%.4f\n", mae, rmse);

%% ===== save =====
modelFile = fullfile(outDir,"Unified_twoStage_nextObs.mat");
save(modelFile, ...
    "net_occ","net_reg_all","net_reg_event", ...
    "lookback","occThresh","pThresh","useSoftMix", ...
    "stationNames","muX","sdX","muY","sdY", ...
    "mae","rmse","maxGapHoursForSample", "-v7.3");

fprintf("Saved: %s\n", modelFile);

%% ===== helper =====
function [X, Yreg, Yocc] = buildXY_nextObs(y, t, doySin,doyCos,ltSin,ltCos, Kp,F107, lookback, occThresh, maxGapHours)
% X includes:
%   - last 12 foEs (12)
%   - spanHours: hours(t(i)-t(i-11))  (1)
%   - dtLastMin: minutes(t(i)-t(i-1)) (1)
%   - dtNextMin: minutes(t(i+1)-t(i)) (1)
%   - cyc DOY/LT (4)
%   - Kp F107 (2)
% Total dims before onehot: 12 + 3 + 6 = 21

    y=double(y(:));
    t=t(:);
    doySin=double(doySin(:)); doyCos=double(doyCos(:));
    ltSin=double(ltSin(:));   ltCos=double(ltCos(:));
    Kp=double(Kp(:)); F107=double(F107(:));

    n = numel(y);
    X=[]; Yreg=[]; Yocc=[];
    lastI = n-1;

    for i = lookback:lastI
        hist = y(i-lookback+1:i);
        y1 = y(i+1);

        if any(~isfinite(hist)) || ~isfinite(y1), continue; end
        if any(~isfinite([Kp(i),F107(i)])), continue; end

        spanHours = hours(t(i) - t(i-lookback+1));
        dtLastMin = minutes(t(i) - t(i-1));
        dtNextMin = minutes(t(i+1) - t(i));

        if ~isfinite(spanHours) || ~isfinite(dtLastMin) || ~isfinite(dtNextMin), continue; end
        if hours(t(i+1) - t(i)) > maxGapHours, continue; end  % drop huge gaps

        exo = [doySin(i),doyCos(i),ltSin(i),ltCos(i), Kp(i), F107(i)];
        X(end+1,:) = [hist.' spanHours dtLastMin dtNextMin exo]; %#ok<AGROW>
        Yreg(end+1,1) = y1; %#ok<AGROW>
        Yocc(end+1,1) = double(y1 >= occThresh); %#ok<AGROW>
    end
end
