function [J, D, a, A0] = parse_sb_params(params1) 
occ = strfind(params1, '_');
for i = occ
	J = str2num(params1(2:occ(1) - 1));
	D = str2num(params1(occ(1) + 2:occ(2) - 1));
	a = str2num(params1(occ(4) + 2:occ(5) - 1));
	A0 = str2num(params1(occ(5) + 3:end));
end