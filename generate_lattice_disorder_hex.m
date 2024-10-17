function onsite_energies = generate_lattice_disorder(L,epsilon,alat,seed,...
    regime,R0)
% Use a Gaussian distribution for random generator
% Use a fix seed
rng(seed,'philox')
ang2bohr = 1.889725989;
ainbohr = alat/sqrt(3)*ang2bohr;
max_onsite = exp(-R0*ainbohr)/(2*ainbohr*epsilon);
onsite_energies = zeros(1,2*L^2);
% Generate random onsite energies in the range [-1/(2eps . a), 1/(2eps . a)]
if(regime==1)
    fac = 1/5000000;
elseif(regime==2)
    fac = 500;
elseif(regime==3)
    fac = 2;
end
emax = fac * max_onsite;
onsite_energies = (randn(1,2*L^2) - 0.5)*emax;
save(join(['disorder_','hex_',num2str(2*L^2),'_',num2str(seed),'_',num2str(regime)],'onsite_energies'))
end
