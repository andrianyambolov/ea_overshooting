fileName = 'output/table.mat';
fileExists = isfile(fileName);

if fileExists
    aux = load(fileName);
    finalT = [aux.finalT; finalT];
    save(fileName, 'finalT');
else
    save(fileName, 'finalT');
end