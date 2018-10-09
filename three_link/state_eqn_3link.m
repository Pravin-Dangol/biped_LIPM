function dx = state_eqn_biped(t,x)

[D,C,G,B,fx,gx] = three_link_matrices(x(1:6));

u = 0; %no control input for now

dx(1:3) = x(4:6);
dx(4:6) = fx+gx*u;
dx = dx';

end