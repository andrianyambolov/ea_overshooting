clear
close all
% clc

load('main_est.mat');
tic
selVars = (1:7)';
selPeriods = 1:36;
maxIdx = 24;
dumStrong = 0;
clear table
hypT = table;
hypT.h = selPeriods';
plotFig = 1;

%% Marginal
if plotFig
    figure;
    for vv = 1:size(data.names,2)
        aux = prctile(squeeze(irfs_draws(vv,1,selPeriods,:)), [qtoplot(2), qtoplot(1), qtoplot(3)]*100,2);
        subplot(ceil(size(data.y,2)/2),2,vv)
        cbound(aux)
        addHorizontalLine;
        title(data.names{vv});
    end
end

%% Joint complete

IRF = squeeze(irfs_draws(:,1,:,:));
IRF = IRF(selVars, selPeriods, :);
IRF = reshape(IRF, size(IRF,1) * size(IRF,2), size(IRF,3))';

[irfabs,IRFcred]=jointabsolute_c_par(IRF);

IRFcred = IRFcred';
IRFcred = reshape(IRFcred, length(selVars), length(selPeriods), []);

sIrf = squeeze(IRFcred(6,:,:))';

[share, selHump] = hypothesisHump(sIrf, maxIdx, dumStrong);

if plotFig
    figure;
    for vv = 1:length(selVars)
        aux = squeeze(IRFcred(vv,:,:));
        subplot(ceil(length(selVars)/2),2,vv)
        hold on
        plot(aux, 'r-')
        plot(aux(:,selHump'), 'g-')
        addHorizontalLine;
        title(data.names{selVars(vv)});
    end
end

hypT.complete_hump = NaN(size(hypT,1),1);
hypT.complete_max = NaN(size(hypT,1),1);
for tt = selPeriods
    hypT.complete_hump(tt) = hypothesisHump(sIrf, tt, dumStrong);
    hypT.complete_max(tt) = hypothesisMax(sIrf, tt);
end


%% Joint simple
selVars = (6)';
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
    for vv = 1:length(selVars)
        aux = squeeze(IRFcred(vv,:,:));
        subplot(ceil(length(selVars)/2),1,vv)
        hold on
        plot(aux, 'r-')
        plot(aux(:,selHump'), 'g-')
        %     plot(aux(:,logical((1-selMax)')), 'y-')
        addHorizontalLine;
        title(data.names{selVars(vv)});
    end
end

hypT.simple_hump = NaN(size(hypT,1),1);
hypT.simple_max = NaN(size(hypT,1),1);
for tt = selPeriods
    hypT.simple_hump(tt) = hypothesisHump(sIrf, tt, dumStrong);
    hypT.simple_max(tt) = hypothesisMax(sIrf, tt);
end
disp(hypT);
toc