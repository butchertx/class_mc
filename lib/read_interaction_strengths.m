function hi = read_interaction_strengths(filename)
data = csvread(filename);
Lx = data(2,2);
Ly = data(2,3);
a = data(3, 2);
tc = data(3, 3);
tcJ = data(4, 2);
gamma = data(4, 3);
beta = data(5, 2);
A0 = data(6, 2);
delta = data(7, 2);
v = data(8, 2);
hi = 1;
end