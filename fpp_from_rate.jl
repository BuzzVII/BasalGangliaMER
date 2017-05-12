#simulation paramters
ξ = 1.0
include("basal_ganglia.jl")
include("generate_timing.jl")

#MER paramters
N = 1
T = 10.0
fs = 24000.0

#BG parameters
global v_e = 14.0
global t_n = 0.0
global dt = 0.0
# function ve(t)
#   global t_n, v_e, dt
#   if t > t_n
#     dt = t - t_n
#   end
#   dv = 10*randn()*sqrt(dt)
#   t_n = t
#   v_e += dv
#   return v_e
# end
ve = (t) -> 14.0

neuron_superposition = [];
rates = [];
sol = [];

delays = true
U = [9.293; 0; 3.857; 0; 14.0; 0; 9.8; 0; -1.78; 0]
for ind = 1:1
  global v_e = 14.0
  global t_n = 0.0
  global dt = 0.0
  rates, sol = bg_sim(ve, delays, U, T)

  #Neuron parameters
  It = squeeze(readcsv("apcurrent24k.dat"),1)
  It = It / maximum(abs(It))

  #the simulation
  neuron_superposition = generate_timing(T, fs, sol, N)
  MER = conv(neuron_superposition, It)
  MER = MER + 0.01*randn(size(MER))
  MER = MER - mean(MER)

  #f = "C:\\Users\\Kristian\\Documents\\PhD\\BG Simulations\\"
  #writedlm(string(f,"times_wiener_",ind,".dat"), find(neuron_superposition))
  #writedlm(string(f,"MER_wiener_",ind,".dat"), MER)
end
