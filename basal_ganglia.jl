using DifferentialEquations

delays = false;

# spikes/s
const S_D1 = 300.0; const S_D2 = 300.0; const S_GPe = 400.0; const S_GPi = 300.0; const S_STN = 500.0;

# mV^-1
const k_D1 = 0.3; const k_D2 = 0.3; const k_GPe = 0.2;const k_GPi = 0.2; const k_STN = 0.2;

# mV
const vh_D1 = 27; const vh_D2 = 27; const vh_GPe = 14; const vh_GPi = 12; const vh_STN = 18.5;

#A mV, D ms^-1,
const A_A = 2*e;   const D_A = 1/4.0*e;
const A_N = 0.1*e; const D_N = 1/100.0*e;
const A_G = 2*e;   const D_G = 1/6.0*e;

#              GPi    GPe  D1   D2  STN  Ctx
const C =  [[   0.0    5   22    0   3    0  ]
            [    0     1    0   33   3    0  ]
            [    0     0    0    0   0   50  ]
            [    0     0    0    0   0   50  ]
            [    0    10    0    0   2    2  ]
            [    0     0    0    0   0    0  ]]';

# delays ms
const τ_1 = 1.0; const τ_3 = 3.0; const τ_4 = 4.0;

# Cortical firing Hz
const μ = 3.0; const σ2 = 1.0;

function sigmoid(v, vh, k, S_max)
  S = S_max/(1+exp(k*(vh - v)));
end

function bg_model(t,u,du)
  du[1] = u[2];
  du[2] = (A_A*D_A + A_N*D_N)*C[5,1]*sigmoid(u[9], vh_STN, k_STN, S_STN) -
          A_G*D_G*C[2,1]*sigmoid(u[3], vh_GPe, k_GPe, S_GPe) -
          A_G*D_G*C[3,1]*sigmoid(u[5], vh_D1,  k_D1,  S_D1) -
          2*(D_A+D_N-2*D_G)*u[2] - (D_A^2+D_N^2-2*D_G^2)*u[1];
  du[3] = u[4];
  du[4] = (A_A*D_A + A_N*D_N)*C[5,2]*sigmoid(u[9], vh_STN, k_STN, S_STN) -
          A_G*D_G*C[4,2]*sigmoid(u[7], vh_D2,  k_D2,  S_D2) -
          A_G*D_G*C[2,2]*sigmoid(u[3], vh_GPe, k_GPe, S_GPe) -
          2*(D_A+D_N-2*D_G)*u[4] - (D_A^2+D_N^2-2*D_G^2)*u[3];
  du[5] = u[6];
  du[6] = (A_A*D_A + A_N*D_N)*C[6,3]*u[11] - 2*(D_A+D_N)*u[6] - (D_A^2+D_N^2)*u[5];
  du[7] = u[8];
  du[8] = (A_A*D_A + A_N*D_N)*C[6,4]*u[11] - 2*(D_A+D_N)*u[8] - (D_A^2+D_N^2)*u[7];
  du[9] = u[10];
  du[10] = (A_A*D_A + A_N*D_N)*C[6,5]*u[11] +
           (A_A*D_A + A_N*D_N)*C[5,5]*sigmoid(u[9], vh_STN, k_STN, S_STN) -
           A_G*D_G*C[2,5]*sigmoid(u[3], vh_GPe, k_GPe, S_GPe) -
           2*(2*D_A+2*D_N-D_G)*u[10] - (2*D_A^2+2*D_N^2-D_G^2)*u[9];
  du[11] = -u[11] + sqrt(σ2) * randn() + μ; #not sure if this is right for the cortical input
end
u0 = zeros(11);
tspan = (0.0,10.0);
prob = ODEProblem(bg_model,u0,tspan)

if delays
  function bg_model(t,u,h,du)
    du[1] = u[2];
    du[2] = (A_A*D_A + A_N*D_N)*C[5,1]*sigmoid(h(t-τ_1)[9], vh_STN, k_STN, S_STN) -
    A_G*D_G*C[2,1]*sigmoid(h(t-τ_1)[3], vh_GPe, k_GPe, S_GPe) -
    A_G*D_G*C[3,1]*sigmoid(h(t-τ_1)[5], vh_D1,  k_D1,  S_D1) -
    2*(D_A+D_N-2*D_G)*u[2] - (D_A^2+D_N^2-2*D_G^2)*u[1];
    du[3] = u[4];
    du[4] = (A_A*D_A + A_N*D_N)*C[5,2]*sigmoid(h(t-τ_1)[9], vh_STN, k_STN, S_STN) -
    A_G*D_G*C[4,2]*sigmoid(h(t-τ_3)[7], vh_D2,  k_D2,  S_D2) -
    A_G*D_G*C[2,2]*sigmoid(h(t-τ_1)[3], vh_GPe, k_GPe, S_GPe) -
    2*(D_A+D_N-2*D_G)*u[4] - (D_A^2+D_N^2-2*D_G^2)*u[3];
    du[5] = u[6];
    du[6] = (A_A*D_A + A_N*D_N)*C[6,3]*h(t-τ_4)[11] - 2*(D_A+D_N)*u[6] - (D_A^2+D_N^2)*u[5];
    du[7] = u[8];
    du[8] = (A_A*D_A + A_N*D_N)*C[6,4]*h(t-τ_4)[11] - 2*(D_A+D_N)*u[8] - (D_A^2+D_N^2)*u[7];
    du[9] = u[10];
    du[10] = (A_A*D_A + A_N*D_N)*C[6,5]*h(t-τ_1)[11] +
    (A_A*D_A + A_N*D_N)*C[5,5]*sigmoid(h(t-τ_1)[9], vh_STN, k_STN, S_STN) -
    A_G*D_G*C[2,5]*sigmoid(h(t-τ_1)[3], vh_GPe, k_GPe, S_GPe) -
    2*(2*D_A+2*D_N-D_G)*u[10] - (2*D_A^2+2*D_N^2-D_G^2)*u[9];
    du[11] = -u[11] + sqrt(σ2) * randn() + μ; #not sure if this is right for the cortical input
  end
  lags = [τ_1, τ_3, τ_4];
  h(t) = zeros(11);
  prob = ConstantLagDDEProblem(bg_model,h,u0,lags,tspan);
  alg = MethodOfSteps(Tsit5());
end

sol = solve(prob);
