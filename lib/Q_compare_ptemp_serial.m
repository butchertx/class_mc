measures = [20, 50, 100, 250, 500, 1000, 10000];
Q_ptemp = [.6822, .6104, .5498, .5122, .4976, .4884, .4791];
Q_normal = [.5061, .5216, .4978, .471, .501, .4946, .4882];


%stats for 100 measurement runs
Q1000_ptemp = [.469859, .479328, .450215, .484375, .500902, .4877, .494562, .495433, .481303, .4884];
Q1000_normal = [.495681, .479774, .494117, .491413, .501017, .502694, .485387, .489786, .505303, .4946];
display('Q1000_ptemp average: ');
ptemp_avg = sum(Q1000_ptemp) / length(Q1000_ptemp)
display('plus/minus');
sqrt(sum((Q1000_ptemp - ptemp_avg).*(Q1000_ptemp - ptemp_avg)) / length(Q1000_ptemp))
display('Q1000_normal average: ');
normal_avg = sum(Q1000_normal) / length(Q1000_normal)
display('plus/minus');
sqrt(sum((Q1000_normal - normal_avg).*(Q1000_normal - normal_avg)) / length(Q1000_normal))