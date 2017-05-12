using DifferentialEquations

const p0 = 0.2; const q0 = 0.3; const v0 = 1; const d0 = 5
const p1 = 0.2; const q1 = 0.3; const v1 = 1; const d1 = 1
const d2 = 1; const beta0 = 1; const beta1 = 1; const tau = 1
function bc_model(t,u,h,du)
  du[1] = (v0/(1+beta0*(h(t-tau)[3]^2))) * (p0 - q0)*u[1] - d0*u[1]
  du[2] = (v0/(1+beta0*(h(t-tau)[3]^2))) * (1 - p0 + q0)*u[1] +
          (v1/(1+beta1*(h(t-tau)[3]^2))) * (p1 - q1)*u[2] - d1*u[2]
  du[3] = (v1/(1+beta1*(h(t-tau)[3]^2))) * (1 - p1 + q1)*u[2] - d2*u[3]
end

lags = [tau]
h(t) = ones(3)
u0 = ones(3)
tspan = (0.0,10.0)
prob = ConstantLagDDEProblem(bc_model,h,u0,lags,tspan)
alg = MethodOfSteps(BS3())
sol = solve(prob,alg)
using Plots; plot(sol)
