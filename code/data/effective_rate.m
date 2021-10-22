clear
close all
clc

bloomberg = readtable('input/bloomberg.xlsx');
rates = table;
rates.time = bloomberg.time;
rates.r1y_au = bloomberg.GTAUD1YRCorp;
rates.r1y_ca = bloomberg.GTCAD1YRCorp;
rates.r1y_gb = NaN(size(bloomberg.GTAUD1YRCorp));
rates.r1y_jp = bloomberg.GTJPY1YRCorp;
rates.r1y_ko = bloomberg.GTKRW1YRCorp;
rates.r1y_no = NaN(size(bloomberg.GTAUD1YRCorp));
rates.r1y_se = bloomberg.GTSEK1YRCorp;

aux = readtable('input/gs1.xls');
rates.r1y_us = aux{:,2};
rates.r1y_eu = bloomberg.GTDEM1YRCorp;

sweden = readtable('input/central_banks/sweden.xlsx', 'Sheet', 'data');
rates.r1y_se(1:120) = sweden{109:228,2};

australia = readtable('input/central_banks/australia.xls', 'Sheet', 'data');
rates = mergeTimeTables(rates, australia);
rates.r1y_au(11) = NaN;
extr = extrapolateRate(rates.r3y_aud, rates.r1y_au, 1, 108);
rates.r1y_au(1:11)  = extr(1:11);
rates.r3y_aud = [];

gb = readtable('input/central_banks/gb.xlsx');
rates.r1y_gb = gb.r1y_gb;


w = readtable('output/weights.csv');
w.no = [];
rates.r1y_no = [];
rates.r1y_ko(1:2) = rates.r1y_ko(3);

rstar = rates{:,2:end-1}.*w{:,2:end};
rstar = sum(rstar,2);

rates.rstar = rstar;
% rates(:,2:7) = [];

writetable(rates, 'output/rates.xlsx', 'WriteMode', 'overwrite');



function Y = extrapolateRate(X, Y, smpl_start, smpl_end)

X = X(smpl_start:smpl_end);
Y = Y(smpl_start:smpl_end);

auxX = X(~isnan(Y));
auxV = Y(~isnan(Y));
% [~, idx] = unique(auxX, 'sorted');

p = polyfit(auxX, auxV, 2);

% aux = interp1(auxX(idx), auxV(idx),X(isnan(Y)), 'linear', 'extrap');

Y(isnan(Y)) = polyval(p, X(isnan(Y)));

end