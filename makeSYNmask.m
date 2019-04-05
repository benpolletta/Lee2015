function mask = makeSYNmask(Npre, Npost, fanout)

UB = max(Npre,Npost);
Xpre = linspace(1,UB,Npre)'*ones(1,Npost);
Xpost = ones(Npre, 1)*linspace(1,UB,Npost);
mask = Xpre-Xpost >= 0 & Xpre-Xpost < fanout | Xpre-Xpost < 0 & UB-(Xpost-Xpre) < fanout;

leftperm = zeros(Npre);
rightperm = zeros(Npost);
left_indices = randperm(Npre);
right_indices = randperm(Npost);
for i = 1:Npre, leftperm(i, left_indices(i)) = 1; end
for i = 1:Npost, rightperm(i, right_indices(i)) = 1; end
mask = leftperm*mask*rightperm;