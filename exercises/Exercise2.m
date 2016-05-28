mu=logspace(-8, -1);

nsteps = zeros(length(mu),1);

for i=1:length(mu);
    nsteps(i) = ex2(mu(i));
end

figure;
semilogx(mu, nsteps);