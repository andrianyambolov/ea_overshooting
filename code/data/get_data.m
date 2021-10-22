clear
close all
clc

data = table;
data.time = generateTime(1999, 1, 2019, 12);

aux = getEcbData('iprod');
data = mergeTimeTables(data, aux);

aux = getEcbData('iprod_construction');
data = mergeTimeTables(data, aux);

aux = getEcbData('cpi');
data = mergeTimeTables(data, aux);

aux = getEcbData('cpi_energy');
data = mergeTimeTables(data, aux);

aux = getEcbData('neer');
data = mergeTimeTables(data, aux);

aux = getEcbData('eurostox');
data = mergeTimeTables(data, aux);

aux = readtable('input/gm_ea_credit_risk_indicators_0.xlsx', 'Sheet', 'Monthly');
aux.date = generateTime(1999, 1, 2021, 9);
aux = aux(:,[1, 3, 6]);
data = mergeTimeTables(data, aux);

aux = readtable('input/fred_spread.xls');
aux.date = generateTime(1999, 1, 2019, 12);
aux.ice = aux{:,2};
aux = aux(:,[end-1, end]);
data = mergeTimeTables(data, aux);

aux = readtable('input/unemp.xls');
aux.date = generateTime(1999, 1, 2019, 12);
aux.unemp = aux{:,2};
aux = aux(:,[end-1, end]);
data = mergeTimeTables(data, aux);

aux = readtable('output/rates.xlsx');
data = mergeTimeTables(data, aux);

aux = readtable('output/fx.csv');
data = mergeTimeTables(data, aux);

aux = readtable('output/neer_ex.xlsx');
data = mergeTimeTables(data, aux);

own_inerp = readtable('interp/output/interpolated.xlsx');
aux = table;
aux.time = own_inerp{:,1};
aux.gdp = own_inerp{:,3};
aux.def = own_inerp{:,4};
data = mergeTimeTables(data, aux);

% cs = readtable('input/cs.xlsx', 'Sheet', 'US');
% aux = table;
% aux.date = generateTime(1979, 7, 2019, 12);
% aux.gdp = cs.RGDP;
% data = mergeTimeTables(data, aux);

jk = readtable('input/jk.csv');
aux = table;
aux.date = generateTime(1979, 7, 2016, 12);
aux.jk_gdp = exp(jk.ea_rgdp/100);
aux.jk_def = exp(jk.ea_gdpdef/100);
data = mergeTimeTables(data, aux);

noLog = {'time', 'spr_nfc_bund_ea', 'spr_bk_bund_ea', 'ice', 'unemp','r1y_us', 'r1y_eu', 'rstar', 'r1y_au', 'r1y_ca', 'r1y_gb', 'r1y_jp', 'r1y_ko', 'r1y_se'};
for ii = 1:size(data,2)
    colnames = data.Properties.VariableNames;
    if ~any(strcmp(string(colnames(ii)), noLog))
        data{:,ii} = 100*log(data{:,ii});
    end
end

aux = readtable('output/shocks.xlsx');
data = mergeTimeTables(data, aux);

data.year = floor(data.time);
data.month = round((data.time - data.year)*12 + 1,0);
data = data(:,[end-1:end, 2:end-2]);

ydict = table;
% ydict.Properties.VariableNames = {'name', 'nice_name', 'nonstationary', 'lb', 'ub'};
for ii = 3:(size(data,2)-6)
    colnames = data.Properties.VariableNames;
    ydict.name(ii-2) = string(colnames(ii));
    ydict.nice_name(ii-2) = string(colnames(ii));
    
    if any(strcmp(string(colnames(ii)), noLog))
        ydict.nonstationary(ii-2) = 0;
    else
        ydict.nonstationary(ii-2) = 1;
    end
end
ydict.lb = cell(size(ydict,1),1);
ydict.ub = cell(size(ydict,1),1);
writetable(ydict, 'output/ydict.csv', 'WriteMode', 'overwrite');

writetable(data, 'output/dataset.csv', 'WriteMode', 'overwrite');

function [out] = getEcbData(filename)
    out = readtable(['input/ecb/', filename, '.xlsx']);
    out.date = year(out.date) + (month(out.date)-1)/12;
    colnames = out.Properties.VariableNames;
    for ii = 1:size(colnames,2)
        if ~(strcmp(string(colnames(ii)),'date') || strcmp(string(colnames(ii)),'s1'))
            out.(string(colnames(ii))) = [];
        end
    end
    out.Properties.VariableNames = {'time', filename};
end


