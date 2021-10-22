clear
close all
clc

weights = readtable('output/weights.csv');
fx = readtable('output/fx.csv');
fx.fx_us = fx.fx_us;
neer = getEcbData('neer');

out = table;
out.time = generateTime(1999, 1, 2019, 12);
out = mergeTimeTables(out, neer);

out = mergeTimeTables(out, fx);
out = mergeTimeTables(out, weights);

out.dI = diffX(log(out.neer),1);

for nn = 2:size(fx,2)
    country = erase(fx.Properties.VariableNames{nn}, 'fx_')
    current_w = out.(country);
    current_fx = out.(['fx_', country]);
    cn = current_w.*diffX(log(current_fx),1);
    neer_ex = (out.dI - cn)./(1-current_w);
    neer_ex = neer_ex * 100;
    neer_ex(1) = 100;
    neer_ex = cumsum(neer_ex);
    out.(['neer_ex_',country]) = neer_ex;
    corr(current_fx, neer_ex)
end

out = out(:,[1, end-size(fx,2)+2:end]);

writetable(out, 'output/neer_ex.xlsx', 'WriteMode', 'overwrite');


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