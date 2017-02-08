using DifferentialEquations

delays = false

#cortical e, D1 d1, D2 d2, GPi p1, GPe p2, STN Ϛ
# spikes/s
const S_p2 = 400.0; const S_p1 = 300.0; const S_d1 = 300.0; const S_d2 = 300.0; const S_Ϛ = 500.0
# mV^-1
const k_p2 = 0.2;const k_p1 = 0.2; const k_d1 = 0.3; const k_d2 = 0.3; const k_Ϛ = 0.2
# mV
const θ_p2 = 14; const θ_p1 = 12; const θ_d1 = 27; const θ_d2 = 27; const θ_Ϛ = 18.5
# decay and rise time constants of membrane
const α = 160; const β = 640
#               p1     p2   d1   d2   Ϛ    c
const ν =  [[    0  -0.03  -0.1   0  0.3   0  ]
            [    0   -0.1    0  -0.3 0.3   0  ]
            [    0     0     0    0   0   1.0 ]
            [    0     0     0    0   0   0.7 ]
            [    0  -0.04    0    0   0   0.1 ]
            [    0     0     0    0   0    0  ]]'
# delays ms
const τ1 = 1.0; τ2 = 2.0; const τ3 = 3.0; const τ4 = 4.0
# Cortical firing Hz
const μ = 3.0; const σ2 = 1.0

function ζ(v, θ, k, S_max)
  S = S_max./(1+exp(k*(θ - v)))
end

function bg_model(t,u,du)
  du[1] = u[2]
  du[2] = ν[5,1]*ζ(u[9], θ_Ϛ, k_Ϛ, S_Ϛ) +
          ν[2,1]*ζ(u[3], θ_p2, k_p2, S_p2) +
          ν[3,1]*ζ(u[5], θ_d1,  k_d1,  S_d1) +
          (β + α)*u[2] + (α*β)*u[1]
  du[3] = u[4]
  du[4] = ν[5,2]*ζ(u[9], θ_Ϛ, k_Ϛ, S_Ϛ) +
          ν[4,2]*ζ(u[7], θ_d2,  k_d2,  S_d2) +
          ν[2,2]*ζ(u[3], θ_p2, k_p2, S_p2) +
          (β + α)*u[4] + (α*β)*u[3]
  du[5] = u[6]
  du[6] = ν[6,3]*u[11] + (β + α)*u[6] + (α*β)*u[5]
  du[7] = u[8]
  du[8] = ν[6,4]*u[11] + (β + α)*u[8] + (α*β)*u[7]
  du[9] = u[10]
  du[10] = ν[6,5]*u[11] +
           ν[2,5]*ζ(u[3], θ_p2, k_p2, S_p2) +
           (β + α)*u[10] + (α*β)*u[9]
  du[11] = -u[11] + sqrt(σ2) * randn() + μ; #not sure if this is right for the cortical input
end
u0 = zeros(11)
tspan = (0.0,10.0)
prob = ODEProblem(bg_model,u0,tspan)
alg = Tsit5()
sol = solve(prob,alg)
