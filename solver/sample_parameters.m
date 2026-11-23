function z = sample_parameters(count, seed)
%SAMPLE_PARAMETERS Generate reproducible inputs without changing the caller RNG.
old=rng; cleanup=onCleanup(@() rng(old));
rng(seed,'twister');z=2*rand(5,count)-1;
end
