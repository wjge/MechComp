#/usr/bin/env python3

import cantera as ct
import numpy as np

# Simulation parameters
p = ct.one_atm # pressure [Pa]
#Tin = 830.0 #unburned gas temperature
Tin = 730.0
#reactants = 'C3H8:1, O2:1' # premixed gas composition
fuel = "C3H8"
#oxidizer = "O2:0.164205,N2:0.720716,CO2:0.045168,H2O:0.024653"
oxidizer = "O2:0.21894,N2:0.720716"

width = 0.03 # m
loglevel = 1 # amount of diagnostic output (0 to 8)

# Solution object used to compute mixture properties, set to the state of the 
# upstream fuel-air mixture
gas = ct.Solution("../mechanisms/pachler/pachler.yaml")
#gas.TPY = Tin, 8.51*ct.one_atm, 'C3H8:0.045258,O2:0.164205,N2:0.720716,CO2:0.045168,H2O:0.024653'
gas.TPY = Tin, 8.51*ct.one_atm, 'C3H8:0.060344,O2:0.21894,N2:0.720716'
#print (gas.T)
#print (gas.Y)
#print(gas.mix_diff_coeffs)
Z=gas.mixture_fraction(fuel, oxidizer)
print("Z={:1.3f}".format(Z))

# Set up flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
f.show_solution()

#print (f.T)
#print (f.Y)
#print(gas.mix_diff_coeffs)
#print(gas.thermal_diff_coeffs)

f.transport_model = 'Mix'
f.solve(loglevel=loglevel, auto=True)
f.show_solution()
print('mixture-average flamespeed = {0:7f} m/s'.format(f.velocity[0]))

#compute flame thickness
z=f.flame.grid
T=f.T
size=z.size-1
grad=np.zeros(size)
for i in range(size):
  grad[i] = (T[i+1]-T[i])/(z[i+1]-z[i])
thickness = (max(T)-min(T))/max(grad)
print('laminar flame thickness= {0:7f} m'.format(thickness))

#c vs z
c1 = np.zeros(z.size)
c2 = np.zeros(z.size)
index_CO = gas.species_index("CO")
index_CO2 = gas.species_index("CO2")
index_O2 = gas.species_index("O2")
index_OH = gas.species_index("OH")
index_H = gas.species_index("H")
index_HO2 = gas.species_index("HO2")
index_H2O2 = gas.species_index("H2O2")

Y_OH = f.Y[index_OH]
Y_H = f.Y[index_H]
Y_HO2 = f.Y[index_HO2]
Y_H2O2 = f.Y[index_H2O2]

omega_OH = f.net_production_rates[index_OH]*gas.molecular_weights[index_OH]
omega_H = f.net_production_rates[index_H]*gas.molecular_weights[index_H]
omega_HO2 = f.net_production_rates[index_HO2]*gas.molecular_weights[index_HO2]
omega_H2O2 = f.net_production_rates[index_H2O2]*gas.molecular_weights[index_H2O2]

#c1min = 4.33565e-10
c1min = f.Y[index_CO][0]+f.Y[index_CO2][0]
#c1max = 0.1707708
c1max = f.Y[index_CO][size]+f.Y[index_CO2][size]
#c2min = 0.009664
c2min = f.Y[index_O2][size]
#c2max = 0.2189
c2max = f.Y[index_O2][0]

for i in range(z.size):
  c1[i] = (f.Y[index_CO][i]+f.Y[index_CO2][i]-c1min)/(c1max-c1min)
  c2[i] = (c2max-f.Y[index_O2][i])/(c2max-c2min)
#with open('propane_egr_cvsY.npy','wb') as output:
with open('propane_cvsY.npy','wb') as output:
#  np.save(output, z)
  np.save(output, c1)
  np.save(output, Y_OH)
  np.save(output, Y_H)
  np.save(output, Y_HO2)
  np.save(output, Y_H2O2)
  np.save(output, omega_OH)
  np.save(output, omega_H)
  np.save(output, omega_HO2)
  np.save(output, omega_H2O2)
#  np.save(output, c2)
#  np.save(output, T)

#Calculate Lewis number
mix_lewis_avg=np.zeros(gas.n_total_species)

#average over Tin+10K to Tadibatic-10K
index_0 = 0
index_1 = f.grid.size-1

for j in range(f.grid.size):
  #print(f.T[j])
  if f.T[j]>Tin+10: 
    index_0=j
    #print(index_0)
    break

for j in reversed(range(f.grid.size)):
  #print(f.T[j])
  if f.T[j]<f.T[f.grid.size-1]-10:
    index_1=j
    #print(index_1)
    break

size_of_reaction_zone = index_1 - index_0 + 1

#the array size is n_total_species, f.grid.size
for i in range(gas.n_total_species):
  lewis = f.thermal_conductivity/f.density_mass/f.cp/f.mix_diff_coeffs_mass[i]
  sum_lewis=0.0
  for j in range(index_0,index_1):
    sum_lewis += lewis[j]
  mix_lewis_avg[i] = sum_lewis/size_of_reaction_zone


#print("===================Lewis number for full mechanism==================")
#print(gas.species_names)
#print(mix_lewis_avg)

print("===================Lewis number for the reduced mechanism==================")
#pass the Lewis numbers to the reduced mechanism
reduced_gas = ct.Solution("./propane/tianfengnonstiff/chem.red53.cti")
lewis_reduced = np.zeros(reduced_gas.n_species)

i=0
for S in reduced_gas.species_names:
  index_fullmech = gas.species_index(S)
  lewis_reduced[i]=mix_lewis_avg[index_fullmech]
  print(str(lewis_reduced[i]) + " " + S)
  i=i+1

print("===================generate clookup table==================")
print(f.T[0], end=' ')
i=0
for S in reduced_gas.species_names:
  index_fullmech = gas.species_index(S)
  if(f.X[index_fullmech][1] > 0.0):
    print(f.X[index_fullmech][1], end=' ')
  else:
    print(0.0, end=' ')
  i=i+1
print()
print(f.T[f.grid.size-2], end=' ')
i=0
for S in reduced_gas.species_names:
  index_fullmech = gas.species_index(S)
  if(f.X[index_fullmech][f.grid.size-2] > 0.0):
    print(f.X[index_fullmech][f.grid.size-2], end=' ')
  else:
    print(0.0, end=' ')
  i=i+1
print()
#print(reduced_gas.species_names)
#print(lewis_reduced)

#f.transport_model = 'Multi'
#f.solve(loglevel)
#f.show_solution()
#print('multicomponent flamespeed = {0:7f} m/s'.format(f.velocity[0]))

