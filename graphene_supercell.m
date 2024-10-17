a = 2.46 % Lattice constant

% Lattice vectors (60° convention)
a1 = a*[1,0];
a2 = a/2*[1,sqrt(3)];

% Rectangular supercell lattice vectors
ncell = 6;
ma1 = ncell*a1
ma2 = [0,(ncell-1)*sin(pi/3)*a]


% Populate rectangular supercell with atoms
is = 1;
tol = 10^(-5);
supercell = [];
for in = -ncell:ncell
    for jn = -ncell:ncell
        c1 = in*a1 + jn*a2;
        c2 = c1 + 1/3*(a1 + a2);
        fc1 = c1*inv([ma1;ma2]);
        fc2 = c2*inv([ma1;ma2]);
        if(fc1(1) >= -tol && fc1(1) < 1 - tol &&  fc1(2) >= -tol && fc1(2) < 1 - tol)
            supercell(is,:) = c1;
            is = is + 1;
        end
        if(fc2(1) >= -tol && fc2(1) < 1 - tol &&  fc2(2) >= -tol && fc2(2) < 1 - tol)
            supercell(is,:) = c2;
            is = is + 1;
        end
    end
end
supercell