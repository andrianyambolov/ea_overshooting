clear
close all
% clc

baseName = 'baseline';
load(['estimation/', baseName,'_mp_1999m01-2019m12_chol.mat']);
tic
selVars = (1:size(data.y,2))';
selPeriods = 1:48;
maxIdx = 12;
dumStrong = 0;
clear table
hypT = table;
hypT.h = selPeriods';
plotFig = 1;

whichModDum = contains(modname, 'baseline');
if whichModDum
    sPos = 6;
else
    sPos = 7;
end

%% Marginal
if plotFig
    figure;
    for vv = 2:size(data.names,2)
        aux = prctile(squeeze(irfs_draws(vv,1,selPeriods,:)), [qtoplot(2), qtoplot(1), qtoplot(3)]*100,2);
        subplot(ceil((size(data.y,2)-1)/2),2,vv-1)
        cbound(aux)
        addHorizontalLine;
        title(data.names{vv});
        xlim([selPeriods(1) selPeriods(end)]);
        xticks(selPeriods([1, 6:6:end]))
        xticklabels(arrayfun(@num2str, selPeriods([1, 6:6:end]), 'UniformOutput', 0))
    end
    
    figure;
    aux = prctile(squeeze(irfs_draws(sPos,1,selPeriods,:)), [qtoplot(2), qtoplot(1), qtoplot(3)]*100,2);
    cbound(aux)
    addHorizontalLine;
    xlim([selPeriods(1) selPeriods(end)]);
    xticks(selPeriods([1, 6:6:end]))
    xticklabels(arrayfun(@num2str, selPeriods([1, 6:6:end]), 'UniformOutput', 0))
    ylim([-0.6, 1]);
    
    figure;
    uip = (squeeze(irfs_draws(end,1,1:end-1,:))' - squeeze(irfs_draws(2,1,1:end-1,:))' - 12*(squeeze(irfs_draws(sPos,1,2:end,:))' - squeeze(irfs_draws(sPos,1,1:end-1,:))'))';
    uip = uip(selPeriods',:);
    aux = prctile(uip, [qtoplot(2), qtoplot(1), qtoplot(3)]*100,2);
    cbound(aux)
    addHorizontalLine;
    xlim([selPeriods(1) selPeriods(end)]);
    xticks(selPeriods([1, 6:6:end]))
    xticklabels(arrayfun(@num2str, selPeriods([1, 6:6:end]), 'UniformOutput', 0))
end

%% Joint complete

IRF = squeeze(irfs_draws(:,1,:,:));
IRF = IRF(selVars, selPeriods, :);
IRF = reshape(IRF, size(IRF,1) * size(IRF,2), size(IRF,3))';

[irfabs,IRFcred]=jointabsolute_c_par(IRF);

IRFcred = IRFcred';
IRFcred = reshape(IRFcred, length(selVars), length(selPeriods), []);

sIrf = squeeze(IRFcred(sPos,:,:))';

[share, selHump] = hypothesisHump(sIrf, maxIdx, dumStrong);

if plotFig
    figure;
    for vv = 2:length(selVars)
        aux = squeeze(IRFcred(vv,:,:));
        subplot(ceil((length(selVars)-1)/2),2,vv-1)
        hold on
        plot(aux(:,1:5:end), 'Color', [0.1412    0.5176    0.8902])
        plot(aux(:,1:5:end), 'Color', [0.1412    0.5176    0.8902])
        [~, idxMax] = max(aux,[],2);
        idxMax = unique(idxMax);
        [~, idxMin] = min(aux,[],2);
        idxMin = unique(idxMin);
        plot(aux(:,[idxMax, idxMax]), 'Color', [0.1412    0.5176    0.8902])
        plot(max(aux,[],2), 'k', 'LineWidth', 2.5, 'Color', [0    0.3020    0.6000]);
        plot(min(aux,[],2), 'k', 'LineWidth', 2.5, 'Color', [0    0.3020    0.6000]);
        plot(aux(:,selHump'), 'Color', [0.8706    0.0784    0.0784], 'LineWidth', 0.5)
        addHorizontalLine;
        title(data.names{selVars(vv)});
        xlim([selPeriods(1) selPeriods(end)]);
        xticks(selPeriods([1, 6:6:end]))
        xticklabels(arrayfun(@num2str, selPeriods([1, 6:6:end]), 'UniformOutput', 0))
    end
end

uip = (squeeze(IRFcred(end,1:end-1,:))' - squeeze(IRFcred(2,1:end-1,:))' - 12*(squeeze(IRFcred(sPos,2:end,:))' - squeeze(IRFcred(sPos,1:end-1,:))'))';
if plotFig
    figure;
    hold on
    plot(uip(:,1:5:end), 'Color', [0.1412    0.5176    0.8902])
    plot(max(uip,[],2), 'k', 'LineWidth', 2.5, 'Color', [0    0.3020    0.6000]);
    plot(min(uip,[],2), 'k', 'LineWidth', 2.5, 'Color', [0    0.3020    0.6000]);
    [~, idxMax] = max(uip,[],2);
    idxMax = unique(idxMax);
    [~, idxMin] = min(uip,[],2);
    idxMin = unique(idxMin);
    plot(uip(:,[idxMax, idxMax]), 'Color', [0.1412    0.5176    0.8902])
    plot(uip(:,selHump'), 'Color', [0.8706    0.0784    0.0784], 'LineWidth', 0.5)
    addHorizontalLine;
    xlim([selPeriods(1) selPeriods(end)]);
    xticks(selPeriods([1, 6:6:end]))
    xticklabels(arrayfun(@num2str, selPeriods([1, 6:6:end]), 'UniformOutput', 0))
end


hypT.complete_hump = NaN(size(hypT,1),1);
hypT.complete_max = NaN(size(hypT,1),1);
for tt = selPeriods
    hypT.complete_hump(tt) = hypothesisHump(sIrf, tt, dumStrong);
    hypT.complete_max(tt) = hypothesisMax(sIrf, tt);
end


%% Joint simple
selVars = (sPos)';
IRF = squeeze(irfs_draws(:,1,:,:));
IRF = IRF(selVars, selPeriods, :);
IRF = reshape(IRF, size(IRF,1) * size(IRF,2), size(IRF,3))';

[irfabs,IRFcred]=jointabsolute_c_par(IRF);
IRFcred = IRFcred';
IRFcred = reshape(IRFcred, length(selVars), length(selPeriods), []);

sIrf = squeeze(IRFcred(1,:,:))';

[share, selHump] = hypothesisHump(sIrf, maxIdx, dumStrong);
[share, selMax] = hypothesisMax(sIrf, 70);

if plotFig
    figure;
    aux = squeeze(IRFcred);
    hold on
    plot(aux(:,1:5:end), 'Color', [0.1412    0.5176    0.8902], 'LineWidth', 0.05)
    [~, idxMax] = max(aux,[],2);
    idxMax = unique(idxMax);
    [~, idxMin] = min(aux,[],2);
    idxMin = unique(idxMin);
    plot(aux(:,[idxMax, idxMax]), 'Color', [0.1412    0.5176    0.8902], 'LineWidth', 0.05)
    plot(max(aux,[],2), 'k', 'LineWidth', 2.5, 'Color', [0    0.3020    0.6000]);
    plot(min(aux,[],2), 'k', 'LineWidth', 2.5, 'Color', [0    0.3020    0.6000]);
    plot(aux(:,selHump'), 'Color', [0.8706    0.0784    0.0784], 'LineWidth', 0.05)
    addHorizontalLine;
    xlim([selPeriods(1) selPeriods(end)]); 
    xticks(selPeriods([1, 6:6:end]))
    xticklabels(arrayfun(@num2str, selPeriods([1, 6:6:end]), 'UniformOutput', 0))
    ylim([-0.6, 1]);
end


hypT.simple_hump = NaN(size(hypT,1),1);
hypT.simple_max = NaN(size(hypT,1),1);
for tt = selPeriods
    hypT.simple_hump(tt) = hypothesisHump(sIrf, tt, dumStrong);
    hypT.simple_max(tt) = hypothesisMax(sIrf, tt);
end
finalT = hypT([1, 6:6:end],:);
disp(finalT);
misc;
toc