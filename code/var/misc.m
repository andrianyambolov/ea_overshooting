fileName = 'output/table.xlsx';
fileExists = isfile(fileName);

if fileExists
    aux = readtable(fileName);
    lastRow = size(aux,1);
    writetable(hypT(1:6:end,:), fileName, ...
        'Range', ['A',int2str(lastRow+2),':','E',int2str(lastRow+8)], ... 
        'WriteVariableNames', 0);
else
    writetable(hypT(1:6:end,:), fileName);
end