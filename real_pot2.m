%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                    %%
%% Function to compute the potential  %%
%% energy in real space from the      %%
%% reciprocal space components V(g)   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                    %%
%% Written by Valerio Vitale          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pot = real_pot2(nsites,epsilon,ma1,ma2,mcell,alpha)
pot = zeros(nsites);
for ir = 1 : nsites
    for jr = ir + 1 : nsites
        for in = -3 : 3
            for jn = -3 : 3
                super_jr = mcell(jr,:) + in*ma1 + jn*ma2;
                norm_irjr = norm(mcell(ir,:)-super_jr);
                %if(norm_irjr < cutoff)
                    %pot(ir,jr) = 1/norm_irjr - 1/sqrt(norm_irjr^2 + alpha^2);
                    pot(ir,jr) = pot(ir,jr) + exp(-alpha*norm_irjr)/norm_irjr;
                %end
            end
        end
    end
end       
pot = pot./epsilon;
end
