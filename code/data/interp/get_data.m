clear
close all
clc

checkDum = 0;
t0 = 1999;
tend = 2019;

jk_m = readtable('input/jk.xls', 'Sheet', 'monthly');
jk_m = jk_m(1:264,:);
jk_m = [jk_m, table(generateTime(1995, 1, 2016, 12))];
jk_m(:,1) = [];
jk_m = [jk_m(:,end), jk_m(:,1:end-1)];

jk_q = readtable('input/jk.xls', 'Sheet', 'quarterly');
jk_q = jk_q(1:88,:);

quarterly = readtable('input/q_main.xlsx');
quarterly = [quarterly, table((1995:0.25:2019.75)')];
quarterly.time = quarterly.Var1;
quarterly.Var1 = [];

if checkDum == 1
    aux_q = quarterly(1:88,:);
    corr(aux_q{:,2:end}, jk_q{:,2:end})
    max(abs(aux_q{:,2:end} - jk_q{:,2:end})) ./ mean(aux_q{:,2:end},1) * 100
end

quarterly = quarterly(quarterly.time > t0 - 0.01 & quarterly.time < tend + 1.01,:);

data = table;
data.time = generateTime(t0, 1, tend, 12);

%% B
aux = getEcbData('input/B.xlsx');
aux.Properties.VariableNames{2} = 'ConsGoodsImports';
aux{:,2} = aux{:,2}/1000;
data = mergeTimeTables(data, aux);
if checkDum
    data = mergeTimeTables(data, jk_m(:,[1, 2]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% C
aux = getEcbData('input/C.xlsx');
aux.Properties.VariableNames{2} = 'RetailTrade';
data = mergeTimeTables(data, aux);
if checkDum
    data = mergeTimeTables(data, jk_m(:,[1, 3]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% E
aux = readtable('input/E.xlsx');
aux = array2table([generateTime(2007, 12, 2019, 12), aux.s1]);
aux.Properties.VariableNames = {'time', 'GovDebtSec'};
data = mergeTimeTables(data, aux);

if checkDum
    data = mergeTimeTables(data, jk_m(:,[1,5]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% H
aux = getEcbData('input/H.xlsx');
aux.Properties.VariableNames{2} = 'ConstructionOutput';
data = mergeTimeTables(data, aux);
if checkDum
    data = mergeTimeTables(data, jk_m(:,[1, 8]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% J
aux = getFredData(['input/J.xls']);
aux{:,2} = aux{:,2}/100;
data = mergeTimeTables(data, aux);
data.Properties.VariableNames{end} = 'StocksFinished';

if checkDum
    data = mergeTimeTables(data, jk_m(:,[1,10]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end


%% K
aux = getFredData(['input/K.xls']);
aux{:,2} = aux{:,2}/100;
data = mergeTimeTables(data, aux);
data.Properties.VariableNames{end} = 'StocksVolume';

if checkDum
    data = mergeTimeTables(data, jk_m(:,[1,11]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% L --- CHECK AGAIN
aux = readtable('input/L.xlsx');
aux{:,end+1} = generateTime(1995, 1, 2019, 12);
aux = aux(:,[3,2]);
aux{:,2} = aux{:,2}*1000;
data = mergeTimeTables(data, aux);
data.Properties.VariableNames{end} = 'NetTrade';

if checkDum
    data = mergeTimeTables(data, jk_m(:,[1,12]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% M 
aux = getFredData(['input/M.xls']);
aux{:,2} = aux{:,2}/100;
data = mergeTimeTables(data, aux);
data.Properties.VariableNames{end} = 'ExportOrders';

if checkDum
    data = mergeTimeTables(data, jk_m(:,[1,13]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% N
aux = getEcbData('input/N.xlsx');
aux.Properties.VariableNames{2} = 'ManNewOrders';
data = mergeTimeTables(data, aux);
if checkDum
    data = mergeTimeTables(data, jk_m(:,[1, 14]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% Q
aux = getEcbData('input/P.xlsx');
aux.Properties.VariableNames{2} = 'PPI';
data = mergeTimeTables(data, aux);
if checkDum
    data = mergeTimeTables(data, jk_m(:,[1, 16]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

%% Q
aux = getEcbData('input/Q.xlsx');
aux.Properties.VariableNames{2} = 'CPI';
data = mergeTimeTables(data, aux);
if checkDum
    data = mergeTimeTables(data, jk_m(:,[1, 17]));
    auxcorr = corr(data{:,end-1}, data{:,end}, 'rows', 'complete')
    auxd = max(abs(data{:,end-1}- data{:,end})) ./ mean(data{:,end-1}, 'omitnan') * 100
    data(:,end) = [];
end

monthly = table;
monthly.time = data.time;
monthly.ConsGoodsImports = data.ConsGoodsImports;
monthly.RetailTrade = data.RetailTrade;
monthly.ConstructionOutput = data.ConstructionOutput;
monthly.StocksFinished = data.StocksFinished;
monthly.StocksVolume = data.StocksVolume;
monthly.NetTrade = data.NetTrade;
monthly.ExportOrders = data.ExportOrders;
monthly.ManNewOrders = data.ManNewOrders;
monthly.PPI = data.PPI;
monthly.CPI = data.CPI;

writetable(monthly, 'output/interp_dataset.xlsx', 'Sheet', 'Monthly', 'WriteMode', 'overwrite');
writetable(quarterly, 'output/interp_dataset.xlsx', 'Sheet', 'Quarterly');

%% Load files

function [out] = getFredData(filename)
    raw = readtable(filename);
    out = table;
    out.time = year(raw{:,1}) + (month(raw{:,1})-1)/12;
    out.(raw.Properties.VariableNames{:,2}) = raw.(raw.Properties.VariableNames{:,2});
end

function [out] = getEcbData(filename)
    out = readtable(filename);
    out.date = year(out.date) + (month(out.date)-1)/12;
    colnames = out.Properties.VariableNames;
    for ii = 1:size(colnames,2)
        if ~(strcmp(string(colnames(ii)),'date') || strcmp(string(colnames(ii)),'s1'))
            out.(string(colnames(ii))) = [];
        end
    end
end



