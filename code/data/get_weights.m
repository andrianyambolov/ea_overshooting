clear
close all
clc

fileNames = dir('input/weights/');
fileN = size(fileNames,1) - 2;
weights = table;
weights.year = (1999:1:2019)';
for ii = 3:size(fileNames,1)
    currentFile = fileNames(ii);
    aux = readtable(['input/weights/',currentFile.name]);
    aux.year = year(aux.date);
    aux = aux(aux.year > 1998 & aux.year < 2020,:);
    country = table;
    country.year = aux.year;
    country.(erase(currentFile.name, '.xlsx')) = aux.s1;
    if country.year(end) == 2018
       country.year(end+1) = 2019;
       country(end,2) = country(end-1,2);
    end
    weights = mergeTimeTables(weights, country);
end

weightsBig = table;
weightsBig.time = generateTime(1999, 1, 2019, 12);
weightsBig = [weightsBig, repelem(weights(:,2:end),12,1)];
weightsBig{:,2:end} = weightsBig{:,2:end}./sum(weightsBig{:,2:end},2);

writetable(weightsBig, 'output/weights.csv', 'WriteMode', 'overwrite');