using DifferentialEquations

function bg_sim(v_e, delays, T, ξ)
  #cortical - e, D1 - d1, D2 - d2, GPi - p1, GPe - p2, STN - Ϛ
  # spikes/s
  const S_p1 = 250.0; const S_p2 = 300.0; const S_d1 = 65.0; const S_d2 = 65.0; const S_Ϛ = 500.0
  # mV^-1
  const k_p1 = 0.2;const k_p2 = 0.2; const k_d1 = 0.3; const k_d2 = 0.3; const k_Ϛ = 0.2
  # mV
  const θ_p1 = 10.0; const θ_p2 = 9.0; const θ_d1 = 19.0; const θ_d2 = 19.0; const θ_Ϛ = 10.0
  # decay and rise time constants of membrane s^-1
  const α = 160.0; const β = 640.0
  # dopamine factor
  #const ξ = 5.0
  # mV s          p1     p2   d1    d2   Ϛ    e
  const ν =  [[    0  -0.03  -0.1   0   0.3   0  ]
              [    0   -0.1    0  -0.3  0.3   0  ]
              [    0     0     0    0    0   1.0/ξ ]
              [    0     0     0    0    0   0.7*ξ ]
              [    0  -0.04    0    0    0   0.1 ]
              [    0     0     0    0    0    0  ]]'
  # delays s
  if delays
    τ = zeros(6,6)
  else
    τ =  [[    0    1.0   1.0   0   1.0   0  ]
          [    0     0     0   1.0  1.0   0  ]
          [    0     0     0    0    0   2.0 ]
          [    0     0     0    0    0   2.0 ]
          [    0    1.0    0    0    0   1.0 ]
          [    0     0     0    0    0    0  ]]' * 0.001
  end
  # Cortical firing Hz
  #const μ = 14.0; const σ2 = 1.0

  function ζ(v, θ, k, S_max)
    S = S_max./(1+exp(k*(θ - v)))
  end

  function dζ(v, θ, k, S_max)
    S = S_max*k*exp(k*(θ - v))./(1+exp(k*(θ - v))).^2
  end

  function bg_model(t,u,du)
    du[1] = u[2]
    du[2] = (α*β)*ν[5,1]*(ζ(u[9], θ_Ϛ, k_Ϛ, S_Ϛ) - τ[5,1]*u[10]*dζ(u[9], θ_Ϛ, k_Ϛ, S_Ϛ)) +
            (α*β)*ν[2,1]*(ζ(u[3], θ_p2, k_p2, S_p2) - τ[2,1]*u[4]*dζ(u[3], θ_p2, k_p2, S_p2)) +
            (α*β)*ν[3,1]*(ζ(u[5], θ_d1,  k_d1,  S_d1) - τ[3,1]*u[6]*dζ(u[5], θ_d1,  k_d1,  S_d1)) -
            (β + α)*u[2] - (α*β)*u[1]
    du[3] = u[4]
    du[4] = (α*β)*ν[5,2]*(ζ(u[9], θ_Ϛ, k_Ϛ, S_Ϛ) - τ[5,2]*u[10]*dζ(u[9], θ_Ϛ, k_Ϛ, S_Ϛ)) +
            (α*β)*ν[4,2]*(ζ(u[7], θ_d2,  k_d2,  S_d2) - τ[4,2]*u[8]*dζ(u[7], θ_d2,  k_d2,  S_d2)) +
            (α*β)*ν[2,2]*(ζ(u[3], θ_p2, k_p2, S_p2) - τ[2,2]*u[4]*dζ(u[3], θ_p2, k_p2, S_p2)) -
            (β + α)*u[4] - (α*β)*u[3]
    du[5] = u[6]
    du[6] = (α*β)*ν[6,3]*v_e(t-τ[6,3]) - (β + α)*u[6] - (α*β)*u[5]
    du[7] = u[8]
    du[8] = (α*β)*ν[6,4]*v_e(t-τ[6,4]) - (β + α)*u[8] - (α*β)*u[7]
    du[9] = u[10]
    du[10] = (α*β)*ν[6,5]*v_e(t-τ[6,5]) +
             (α*β)*ν[2,5]*(ζ(u[3], θ_p2, k_p2, S_p2) - τ[2,5]*u[4]*dζ(u[3], θ_p2, k_p2, S_p2))-
             (β + α)*u[10] - (α*β)*u[9]
  end
  u0 = zeros(10)
  #u0[[1,3,5,7,9]]=[θ_p1, θ_p2, θ_d1, θ_d2, θ_Ϛ]
  tspan = (0.0,T)
  prob = ODEProblem(bg_model,u0,tspan)
  alg = Tsit5()
  sol = solve(prob, alg)
  rates = [ζ(sol[:,1], θ_p1, k_p1, S_p1) ζ(sol[:,3], θ_p2, k_p2, S_p2) ζ(sol[:,5], θ_d1,  k_d1,  S_d1) ζ(sol[:,7], θ_d2,  k_d2,  S_d2) ζ(sol[:,9], θ_Ϛ, k_Ϛ, S_Ϛ)]
end
