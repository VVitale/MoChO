%% Input parameters
   % Moiré cell lattice constant
   alat = 50;% Ang
   tfac = 0.95; % factor used to decrease T (T^(i+1) = tfac*T^(i))
   max_iter = 200000; % Number of iterations in external loop
   L = 6; % Linear dimension of triangular lattice. Number of sites is L^2
   restart = false; % Whether it's a restart calculation
   nrep = 1; % Number of iterations in the Montecarlo step
   %rng(1);
   
   % Adjustable parameters
   epsilon = 4.5; % Value of dielectric constant in computation of V(R)
   % Electronic temperature
   Treal = 50000; % Kelvin
   seed = 0; % Seed for random generator
   add_disorder = false;
   filling = 2;
   regime = 1;
   R0 = 0.01;%5*1.889725989;
   lattice = 'hex';
   onsite_hex = -0.01; % Onsite potential in eV
   onsite_tri = 0.2; % Onsite potential in eV
   Ef = 0.0; % Electric field in V/Hartree
%% End Input parameters

nsites = L^2; % Number of moire cells, aka sites

if(strcmp(lattice,'tri'))
if(add_disorder)
    % Load disorder energies from file
    filename = join(['disorder_',lattice,'_',num2str(L^2),'_',num2str(seed),'_',num2str(regime)]);
    if(~exist(filename,'file'))
        %Generate disorder energies
        disorder = generate_lattice_disorder(L,epsilon,alat,seed,regime,R0);
    end
else
    disorder = zeros(1,nsites);
end
disorder

[charge_state_best,ics_best,E_min,potential,mcell] = simul_annealing(alat,tfac,max_iter,...
    L,restart,nrep,epsilon,Treal,nsites,filling,disorder,add_disorder,R0);
elseif(strcmp(lattice,'hex'))
    if(add_disorder)
    % Load disorder energies from file
    filename = join(['disorder_',lattice,'_',num2str(2*L^2),'_',num2str(seed),'_',num2str(regime)]);
    if(~exist(filename,'file'))
        %Generate disorder energies
        disorder = generate_lattice_disorder_hex(L,epsilon,alat,seed,regime,R0);
    end
else
    disorder = zeros(1,nsites);
end
disorder

[charge_state_best,ics_best,E_min,potential,mcell] = simul_annealing_hex(alat,tfac,max_iter,...
    L,restart,nrep,epsilon,Treal,nsites,filling,disorder,add_disorder,R0,onsite_hex,onsite_tri,Ef);
end
