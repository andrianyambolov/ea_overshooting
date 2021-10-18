clear
close all
clc

tend = datenum(2019, 12, 12); 
% tend = datenum(2019, 6, 6); 
tdel1 = datenum(2001, 9, 13);
tdel2 = datenum(2001, 9, 17);
tdel3 = datenum(2008, 10, 8);


raw = readtable('input/Dataset_EA-MPD.xlsx', 'Sheet', 'Monetary Event Window');
raw = raw(1:datefind(tend, raw.date),:);
raw(datefind(tdel1, raw.date),:) = [];
raw(datefind(tdel2, raw.date),:) = [];
raw(datefind(tdel3, raw.date),:) = [];

ois = raw(:,[1,3:6]);
[coef, score, ~, expl] = pca(ois{:,2:end});

i_tot = score(:,1)*(std(ois{:,end}) / std(score(:,1)));
s = raw.STOXX50;

i_mp_poor = (sign(i_tot.*s) == -1).*i_tot;
i_cbi_poor = (sign(i_tot.*s) == 1).*i_tot;

M = [i_tot, s];
[Q, R] = qr(M,0);
alfa = acos(sqrt(var(i_mp_poor)/var(i_tot)));
P = [cos(alfa), sin(alfa); 
    -sin(alfa), cos(alfa)];
D = [R(1,1)*cos(alfa), 0; 
    0, R(1,1)*sin(alfa)];
U = Q*P*D;

nonsum = table;
nonsum.time = (month(ois.date)-1)/12 + year(ois.date);
nonsum.mp_poor = i_mp_poor;
nonsum.cbi_poor = i_cbi_poor;
nonsum.mp = U(:,1);
nonsum.cbi = U(:,2);
nonsum.s = s;
nonsum.i_tot = ois.OIS_1Y;

timevec = (1999:1/12:2019 + 11/12)';
out = table;
out.time = timevec;
out.mp_poor = NaN(size(timevec));
out.cbi_poor = NaN(size(timevec));
out.mp = NaN(size(timevec));
out.cbi = NaN(size(timevec));
out.s = NaN(size(timevec));
out.i_tot = NaN(size(timevec));


for tt = 1:size(out.time,1)
    aux = nonsum{nonsum.time == out.time(tt),2:end};
    out{tt,2:end} = sum(aux,1);
end

writetable(out, 'output/shocks.xlsx', 'WriteMode', 'overwritesheet');