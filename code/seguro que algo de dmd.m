%%
samplingRate = 3.75*1000;
dt = 1/samplingRate;
[nn,mm,pp] = size(ux);
uxx = ux(:,:,pp-2999:pp);
Uxx = zeros(nn*mm,3000);
for kk = 1:3000
    tempvec = ux(:,:,kk);
    Uxx(:,kk) = tempvec(:);
end
t = (0:3000)*dt;
V = Uxx;
Time = t;
threshold = 0.7E-2;
[Vreconst,deltas,omegas,amplitude]=DMD1_SIADS(V,Time,threshold,threshold);
