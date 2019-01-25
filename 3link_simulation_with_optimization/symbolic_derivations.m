%This .m files derives D C G matrix for the state equation, using a defined
%forward position kinematics, B matrix must be defined manually 

%De and E matric also derived for use in impact map, pe vector must be defined

%End of the file contains derivation of Lfh, L2fh, LgLfh for computing control action 
%and Lfb, dLfb for ZD

%% DCG matrices

syms q1 q2 q3 p_h p_v dp_h dp_v dq1 dq2 dq3 real
syms r m Mh Mt l g real

%Mh -  mass of hip, Mt - mass of torso, m - mass of legs
%l - length from hip to torso, r - length of legs

%[r,m,Mh,Mt,l,g] = model_params_3link;

q = [q1; q2; q3];
dq = [dq1; dq2; dq3];

qN = q1; dqN = dq1;
p_e = [p_h; p_v];

%       *               *torso
%       |
%       |<--.
%  'q1  *    ' q3
%  '-->/ \   '          *hip
%  '  /---\--
%  ' /---->\
%  ./  q2   \.
%
% CCW - positive
% q1 is cyclic, and negative pre-impact using this convention

%Position of masses, relative to q1
p_Mh = r*[sin(-q1); cos(-q1)] + p_e;
p_Mt = p_Mh + l*[sin(-q1 + pi - q3); cos(-q1 + pi - q3)];
p_m1 = r/2*[sin(-q1); cos(-q1)] + p_e;
p_m2 = p_Mh + r/2*[sin(q1 + q2); -cos(q1 + q2)];
P2 = p_Mh + r*[sin(q1 + q2); -cos(q1 + q2)];

dp_Mh = jacobian(p_Mh,q)*dq;
dp_Mt = jacobian(p_Mt,q)*dq;
dp_m1 = jacobian(p_m1,q)*dq;
dp_m2 = jacobian(p_m2,q)*dq;

K_Mh =  Mh/2*(dp_Mh)'*dp_Mh;
K_Mt = Mt/2*(dp_Mt)'*dp_Mt;
K_m1 = m/2*(dp_m1)'*dp_m1;
K_m2 = m/2*(dp_m2)'*dp_m2;

K = K_m1 + K_Mh + K_Mt + K_m2;

V_Mh = p_Mh(2)*Mh*g;
V_Mt = p_Mt(2)*Mt*g;
V_m1 = p_m1(2)*m*g;
V_m2 = p_m2(2)*m*g;

PE = V_m1 + V_Mh + V_Mt + V_m2;

k2 = jacobian(PE,qN);
D = jacobian(jacobian(K,dq).',dq); D = simplify(D);

D = jacobian(jacobian(K,dq).',dq); D = simplify(D);

N = max(size(q));
syms C
for k = 1:N
    for j = 1:N
        C(k,j) = 0*g;
        for i = 1:N
            C(k,j) = C(k,j) + 1/2*(diff(D(k,j),q(i)) + ...
                diff(D(k,i),q(j)) - ...
                diff(D(i,j),q(k)))*dq(i);
        end
    end
end
C = simplify(C);

G = jacobian(PE,q).'; G = simplify(G);

B = [0 0;
    1 0;
    0 1];


%% Impact map

%Using same psotion vectors as above, but taking partial with respect to qe
%instead

%Extended configuration variables
qe = [q; p_h; p_v];
dqe = [dq; dp_h; dp_v];

dp_Mh = jacobian(p_Mh,qe)*dqe;
dp_Mt = jacobian(p_Mt,qe)*dqe;
dp_m1 = jacobian(p_m1,qe)*dqe;
dp_m2 = jacobian(p_m2,qe)*dqe;

K_Mh =  Mh/2*(dp_Mh)'*dp_Mh;
K_Mt = Mt/2*(dp_Mt)'*dp_Mt;
K_m1 = m/2*(dp_m1)'*dp_m1;
K_m2 = m/2*(dp_m2)'*dp_m2;

K = K_m1 + K_Mh + K_Mt + K_m2;

De = jacobian(jacobian(K,dqe).',dqe); De = simplify(De);

E = jacobian(P2,qe);


%% Bezier poly
syms s q1_plus delq
%s = (q1 - q1_plus)/delq; 
%delq = q1_minus - q1_plus (always -ve if q1 goes from +ve to -ve)
%ds/dt = dq1/delq; ds/dq1 = 1/delq;

syms a21 a22 a23 a24 a25 
syms a31 a32 a33 a34 a35

a2 = [a21 a22 a23 a24 a25];
a3 = [a31 a32 a33 a34 a35];
M = 4;

b2 = 0; b3 = 0;
for k = 0:M
    b2 = b2 + a2(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

for k = 0:M
    b3 = b3 + a3(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

db_ds2 = 0;
for k = 0:M-1
    db_ds2 = db_ds2 + (a2(1,k+2)-a2(1,k+1))*(factorial(M)/(factorial(k)*factorial(M-k-1)))*s^k*(1-s)^(M-k-1);
end

db_ds3 = 0;
for k = 0:M-1
    db_ds3 = db_ds3 + (a3(1,k+2)-a3(1,k+1))*(factorial(M)/(factorial(k)*factorial(M-k-1)))*s^k*(1-s)^(M-k-1);
end

%
db_ds2 = jacobian(b2,q1);
db_ds3 = jacobian(b3,q1);

h = [q2 - b2; q3 - b3]; hd = [b2; b3];

x = [q1, q2, q3, dq1, dq2, dq3]';

dh_dx = jacobian(h,x);

fx = D\(-C*x(4:6)-G); Fx = [x(4:6); fx];
gx = D\B; Gx = [zeros(3,2); gx];

Lfh = dh_dx*Fx;

dLfh = jacobian(Lfh,x);
%}
%{
%% For Zero Dynamics

Lsb2 = jacobian(b2,s)*dq1/delq; %ds_dt = dq1/delta_theta

Lsb3 = jacobian(b3,s)*dq1/delq; %ds_dt = dq1/delta_theta

%jacbian is same as the sum, so using jacobian for 2nd partial

dLsb2 = jacobian(db_ds2*dq1/delq,s);

dLsb3 = jacobian(db_ds3*dq1/delq,s);

beta1 = [dLsb2; dLsb2]*dq1/delq;


%% For Controller

h = [q2 - b2; q3 - b3]; hd = [b2; b3];

x = [q1, q2, q3, dq1, dq2, dq3]';

ds_dt = dq1/delq; 

dh_dx = jacobian(h,[s;x(2:end)]);

fx = D\(-C*x(4:6)-G); Fx = [ds_dt; x(5:6); fx];
gx = D\B; Gx = [zeros(3,2); gx];

Lfh = dh_dx*[ds_dt;x(5:6); fx];

dLfh = jacobian(Lfh,[s;x(2:end)]);

L2fh = dLfh*Fx;

LgLfh = dLfh*Gx;

beta1 = jacobian(jacobian(hd,s)*(dq1/delq),s)*(dq1/delq);

beta2 = [db_ds2; db_ds3]/(-delq);

%u = LgLfh\L2fh;        %Takes a long time to compute, doesn't display well


ddq1 = (D(1,1)+D(1,2:3)*beta2)\(-D(1,2:3)*beta1-C(1,1)*dq1-G(1,1));

dd_q1 = -G(1,1);
%}
