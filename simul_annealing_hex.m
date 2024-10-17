%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulated Annealing Montecarlo code%%
%% to search for the ground state     %%
%% of a gas of classical charges      %%
%% interacting via a screened Coulomb %%
%% potential on a 2D triangular       %%
%% lattice at different fillings.     %%
%% We use a standard Metropolis       %%
%% algorithm for the acceptance rate. %%
%% The poisson equation is discretised%%
%% and solved in reciprocal space.    %%
%% There are two loops: an external   %%
%% loop that controls the temperature.%%
%% At each fix temperature we run a   %%
%% standard Montecarlo algorithm      %%
%% (internal loop)                    %%
%%                                    %%
%% Ref.                               %%
%% 1) Jong-Rim et al, PRB 46,6 (1992) %%
%%                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Written by Valerio Vitale          %%
%%                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [charge_state_best,ics_best,E_min,potential,mcell,...
    onsite_energies] = simul_annealing_hex(alat,tfac,max_iter,L,restart,...
    nrep,epsilon,Treal,ncells,filling,disorder,add_disorder,R0,onsite1,...
    onsite2,Ef)
   
   % Useful conversion factors
   ang2bohr = 1.889725989; 
   eV2Ha = 0.0367493;
   alat = alat*ang2bohr; % Bohr
   onsite1 = onsite1*eV2Ha
   onsite2 = onsite2*eV2Ha
   Ef = Ef*eV2Ha
   kb =  8.617333262145e-5 * eV2Ha ;% Hartree K^-1;
   unit2K = kb*(epsilon*alat); % Unit of energy Hartree
   Ti = Treal*unit2K % Temperature in Hartree
   
   % Initialise charges according to filling
   sites_per_cell = 3;
   nsites = sites_per_cell*ncells;
   ncharges = floor(ncells*filling);
   if ncharges > 3*ncells
       ncharges = 6*ncells - ncharges;
   end
   % Disp the number of charges
   ncharges
   
   % If this is not a restart initialise geometry
   if(~restart)
   % Set up unit cell geometry
   % Geometry 120° convention
   a1 = alat*[1,0,0];%[1/2,sqrt(3)/2,0];
   a2 = alat*[1/2,sqrt(3)/2,0];%[1/2,-sqrt(3)/2,0];
   %a3 = a2 - a1;
   % Reciprocal vectors
   %b1 = 2*pi/alat * [1, -1/sqrt(3)];
   %b2 = 2*pi/alat * [0, 2/sqrt(3)];
   
   % Supercell
   mcell = zeros(nsites,3);
   ma1 = L*a1;
   ma2 = L*a2;
   % Supercell in real space
   onsite_energies = zeros(nsites,1);
   il = 1;
   for in = 0 : L-1
      for jn = 0 : L-1
          mcell(il,:) = in*a1 + jn*a2 + [0,0,3.5]*ang2bohr;
          onsite_energies(il) = 0.0;
          mcell(il+1,:) = 1/3*(a1+a2) + in*a1+jn*a2 + [0,0,3.5]*ang2bohr;
          onsite_energies(il+1) = onsite1;
          mcell(il+2,:) = in*a1+jn*a2 ;
          onsite_energies(il+2) = onsite2 + Ef;
          il = il + 3;
%          il = il + 1;
      end
   end
   onsite_energies
   
   % Initialise arrays of nearest and next-nearest neighbors for each site
   % considering PBC.
   nn = zeros(nsites,6);
   nnn = zeros(nsites,6);
   intra_nn = zeros(nsites,10);
   for in = 1 : nsites
      iat = mcell(in,:);
      inn = 0;
      innn = 0;
      intra_inn  = 0;
      for jn = 1 : nsites
          for ia = -1 : 1
              for ib = -1 : 1
                  jat = mcell(jn,:) + ia*ma1 + ib*ma2;
                  dist_ij = norm(iat(1:2) - jat(1:2));
                  if(iat(3) == 0.0 && jat(3) == 0.0)
                      if(dist_ij~=0 && abs(dist_ij-alat/sqrt(3))<0.01)
                          inn = inn + 1;
                          if(inn > 3)
                              error('Too many nearest neighbours')
                          end
                          nn(in,inn) = jn;
                      elseif(dist_ij > alat/sqrt(3) && abs(dist_ij - alat)<0.1)
                          innn = innn+1;
                          if(innn > 6)
                              error('Too many next-nearest neighbors')
                          end
                          nnn(in,innn) = jn;
                      end
                  elseif(iat(3) > 0.0 && jat(3) > 0.0)
                      if(dist_ij~=0 && abs(dist_ij-alat)<0.01)
                          inn = inn + 1;
                          if(inn > 6)
                              error('Too many nearest neighbours')
                          end
                          nn(in,inn) = jn;
                      elseif(dist_ij > alat && abs(dist_ij - alat*sqrt(3))<0.1)
                          innn = innn+1;
                          if(innn > 6)
                              error('Too many next-nearest neighbors')
                          end
                          nnn(in,innn) = jn;
                      end
                  elseif((iat(3) > 0.0 && jat(3) == 0.0) || iat(3) == 0.0 && jat(3) > 0.0)
                      if(dist_ij==0.0 || abs(dist_ij-alat)<0.01)
                          intra_inn = intra_inn + 1;
                          if(intra_inn > 10)
                              error('Too many next-nearest neighbors')
                          end
                          intra_nn(in,intra_inn) = jn;
                      end
                  end         
              end
          end
      end
   end
   
   % Initialise charge state
   charge_state = zeros(nsites,1);
   charge_state_best = zeros(nsites,1);
   % Randomly distribute charges on the lattice
   ics = randsample(nsites,ncharges);
   charge_state(ics) = 1;
   
   % Compute on-site potential
   %potential = zeros(nsites);
   potential = real_pot2(nsites,epsilon,ma1,ma2,mcell,R0);
   potential = potential + potential';
   % Initial Energy
   E = 0.5*(charge_state'*potential*charge_state) + dot(charge_state,onsite_energies);
   ics_best = ics;
   if(add_disorder)
       % Add disorder energies
       onsite_disorder = dot(charge_state,disorder);
       E = E + onsite_disorder;
   end
   % Initialise variables
   %Ehist = zeros(max_iter*nrep,1);
   rejected = 0;
   itemp = 0;
   ineigh = [];
   E_min = 10^4;
   E_old = E;
   iter = 0;
   else
   % If this is a restart, read in data
       E = E_min;
       ics = ics_best;
       charge_state = charge_state_best;
       itemp = 0;
       ineigh = [];
       E_old = E;
       iter = 0;
       %Ti = Tf;
   end
   
   % Start simulated annealing
   % External loop
   while iter < max_iter
       iter = iter + 1;
       temp_iter = 0;
       if(mod(iter,1000)==0)
           Ti = Ti + 100;
        % Print the current iteration
        disp(join(['Iter: ',num2str(iter)]))
       end
       % Start Montecarlo loop at constant T
       while temp_iter < nrep
           temp_iter = temp_iter + 1;
           ineigh = [];
           % New state swap i-->j
           % Randomly choose an occupied site
           while isempty(ineigh)
               itemp = randsample(ics,1);
               crand = rand(1);
               if(crand <= 1/3)
                   % Randomly choose an empty nearest or next-nearest
                   % unoccupied site
                   ineigh = find(~charge_state(nonzeros(nn(itemp,:)))); 
               elseif(crand > 1/3 && crand <= 2/3)
                   ineigh = find(~charge_state(nonzeros(nnn(itemp,:))));
               else
                   ineigh = find(~charge_state(nonzeros(intra_nn(itemp,:))));
               end
           end
           if(crand <= 1/3)
               ineigh_temp = randsample(nonzeros(nn(itemp,ineigh)),1);
           elseif(crand > 1/3 && crand <= 2/3)
               ineigh_temp = randsample(nonzeros(nnn(itemp,ineigh)),1);
           else
               ineigh_temp = randsample(nonzeros(intra_nn(itemp,ineigh)),1);
           end
           % Create virtual state with swapped i,j elements
           virtual_state = charge_state;
           virtual_state([itemp,ineigh_temp]) = charge_state([ineigh_temp,itemp]);
           new_ics = find(virtual_state);
           % Check total charge is OK
           virtual_ncharges = size(new_ics);
           if(virtual_ncharges ~= ncharges)
               error('Number of charges is not preserved during swapping')
           end
           % Compute Energy of virtual state
           Ev = 0.5*(virtual_state'*potential*virtual_state) + dot(virtual_state,onsite_energies);
           if(add_disorder)
               % Add disorder energies
               onsite_disorder = dot(virtual_state,disorder);
               Ev = Ev + onsite_disorder;
           end
           % Acceptance rule
           if(Ev < E  || exp(-(Ev-E)/Ti)>rand(1))
               charge_state = virtual_state;
               ics = new_ics;
               E=Ev;
               % Save local minima
               if(E < E_min)
                   E_min = E;
                   ics_best = ics;
                   charge_state_best = charge_state;
               end
           else
               rejected = rejected + 1;
           end
           %Ehist(nrep*(iter-1) + temp_iter) = E;
       end
%        % If 0 accepted moves in the constant T loop break
%        if(E==E_old)
%            break
%        else
%            E_old = E;
%        end
       % Update temperature
       Ti = Ti*tfac;%*unit2K;
   end
   % Compute number of accepted Montecarlo swaps
   accpt = (iter+1)*(nrep) + (temp_iter+1) - rejected;
   total_ratio = accpt/((iter+1)*(nrep) + (temp_iter+1))
   Tf = Ti %*tfac^(iter)
   % Count average number of nearest-neighbor and second nearest neighbor
%    sum_nn = sum(charge_state_best(nn(ics_best(1),:)))
%    sum_nnn = sum(charge_state_best(nnn(ics_best(1),:)))
%    tot_avg_n = (sum_nn + sum_nnn)/12
   % Save state as a matlab file
   filename = join(['rec_filling_',num2str(filling),'.mat']);
   save(filename,'charge_state_best','ics_best','E_min')%,'Ehist')
   % Plot histogram
   %figure
   %histogram(Ehist);
   %set(gca,'FontSize',24)
   %saveas(gca,join([filename,'_hist.fig']))
   %close
   %figure
   %plot(Ehist);
   %set(gca,'FontSize',24)
   %saveas(gca,join([filename,'_fluctuations.fig']))
   %close
   
   % Plot supercell of minimum energy state
    plot_min_state_hex
    clear nn nnn onsite_pot accpt rejected virtual_state charge_state ics ...
        ri0 ri1 ineigh new_ics itemp
end
