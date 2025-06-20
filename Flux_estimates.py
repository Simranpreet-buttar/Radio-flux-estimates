#Contribution by Simranpreet Kaur and Rebeva Pirvu-Malanda.
#In order to have an estimate of the radio flux of the studied sub-stellar objects, three different flux estimates models are used here below: 
#1) FARELL AND DESCH MODEL 2) STEVENS MODEL and 3) GRIESSMEIER AND ZARKA MODEL. 
#These flux estimates can indicate us if these objects are observable or not and can be useful for future proposals. 

### References:
    ## 1.FARRELL DESCH MODEL: https://ui.adsabs.harvard.edu/abs/1999JGR...10414025F/abstract
    ## 2.STEVENS MODEL: Stevens, I. R. 2005, MNRAS, 356, 1053
    ## 3.GRIESSMEIER AND ZARKA: Grießmeier, J. M., Zarka, P., & Spreeuw, H. 2007, A&A, 475, 359

import math
from math import pi
import matplotlib.pyplot as plt
import sys
import numpy as np

# Loading the CSV file
#csv_filename = 'planetary_data.csv'
#df = pd.read_csv(csv_filename)

# Get the planetary system name and parameters
#system_name = config['Parameters']['PlanetarySystem']

# Filtering the DataFrame to get data for the specified system
#system_data = df[df['pl_name'] == system_name]

#if not system_data.empty: 
#    M_planet = system_data['pl_massj'].values[0]
 #   A = system_data['pl_orbsmax'].values[0]
  #  R_planet = system_data['pl_radj'].values[0]
  #  s = system_data['sy_dist'].values[0]
  #  M_star = system_data['st_mass'].values[0]
  #  R_star = system_data['st_rad'].values[0]
  #  t1 = system_data['stellar_age1'].values[0]    ### We should add these values by hand 
  #  t2 = system_data['stellar_age2'].values[0]
    
    
    # Access other columns as needed
#else:
 #   print(f"Data for '{system_name}' not found in the CSV file.")
 
#System = str(input('Planetary system: '))                                                                     #system name
#p_rot = float(input('planetary rotation period_in hours: '))                                                  #Planetary rotation period in hours
#M_planet = float(input('Mass of the planet in Jupiter masses, m_p: '))                                        #Mass of planet in terms of Jupiter masses
#A = float(input('Semi major axis of the planet in au, A: '))                                                  #Distance between the planet and the star in AU
#s = float(input('Distance of the system from Earth in pc, s: '))                                              #Distance between the system and the observer(Earth) in Parsecs
#R_planet = float(input('Radius of the planet in Jupiter radii, R_p: '))                                       #Radius of the planet in Jupiter radii
#M_star = float(input('Mass of the star in Solar masses: '))                                                   # Stellar mass in terms of solar mass
#t1 = float(input('Stellar age in Gyr, t1:'))                                                                  #lower limit on age
#t2 = float(input('Stellar age in Gyr, t2:'))                                                                  #upper limit on age

#WE ARE WORKING IN THE DIRECTION TO EXTRACT THE STELLAR AND PLANETARY PARAMETERS DIRECTLY FROM NASA EXOPLANETS ARCHIVE

System = "bd004475B "
p_rot = 10
M_planet = 50.9
A = 42.6
s = 1.48
R_star = 1
R_planet = 1
M_star = 0.81
t1 = 1
t2 = 5

d_s = A / (R_star*0.00465)     #Star-planet separation in terms of stellar radii, required while calcualting Alfven Mach number, R_star converted to AU

# Constants

p_rot_j = 10.0                   #Jovian period in hours  
m_j = 1.0                        #Jupiter mass, in Jupiter masses
A_j = 5.2                        #Separation Sun-Jupiter in AU
R_j = 1.0                        #Radius of jupiter, in Jupiter radii                          
U_j = 1.6*10**30                 #Magnetic moment of Jupiter in G·cm⁻³(Stevens 2005, after eq.14)
M_sun = 1.0                      #Solar mass in Solar masses                                    
MLR_sun = 2*10**(-14)*M_sun      #Solar Mass loss rate in kg/yr (Stevens 2005, after eq.11)
R_sun = 1.0                      #Solar Radii in Solar radius                                                                        
L_sun = 1.0                      #Solar Luminosity in Solar luminosity                                                                             
R_sj = 40*R_j                    #Magnetospheric radius of Jupiter in Jupiter radius, (Griessmeier-07, eq.31)                               
Tau = 0.0256                     #Time constant used to get the wind density and wind velocity in Gyr (Griessmeier-07, after eq. 16)       
v0 = 3971000                     #wind velocity assumed from the solar system in m/s (not exact values) (Griessmeier-07, see after eq. 16) 
n0 = 1.04*10**11                 #wind density assumed from the solar system in m³ (not exact values) (Griessmeier-07, see after eq. 16)
n_j = 2.0*10**5                  #Jupiter wind density in m³, (Griessmeier, after eq 31)   
V_j = 520000                     #Stellar wind velocity for Jupiter in m/s, (Griessmeier, after eq 31)                                
Flux_j_1AU = 10**11              #Radio flux of Jupiter at 1 AU in mJy (Stevens 2005, after eq. 14)
d0 = 1.0   	                 #1AU 						
G = 6.67*10**(-11)               #Gravitational constant in SI units, m^3kg^(-1)s^(-2)                                                                
B_r0 = 2.6*10**(-9)              #Radial Interplanetary magnetic field for the Solar System in Tesla (Griessmeier, after eq.18)  
B_phi0 = 2.4*10**(-9)            #Azimuthal Interplanetary magnetic field for the Solar System in Tesla  (Griessmeier, after eq.19) 
t_sun= 4.6                       #Age of the sun in Gyr
beam_rad = 1.6                   #Solid angle in steradians (Griessmeier, after eq.7)
B_j= 4
B_sun= 10                        #Solar magnetic field at the surface in gauss
  
#conversions

pc_to_m = 3.086E+16              #Conversion from parsec to meter
pc_to_AU = 206265                #Conversion from parsec to AU
AU_to_m = 1.49*10**11            #Conversion from AU to meter
Wm2Hz_to_mJy = 1.0e29            #Conversion from W/m2/Hz to mJy
u = 1.257*10**(-6)               #Magnetic permeability  of free space in N/A**2
P_ja= 4*(10**9)                  #Radio power of the Jovian hectometric component in Watts (Farrel and Desch eq 6)
P_jb= 4*(10**11)                 #Radio power of the Jovian decametric component in Watts (Farrel and Desch eq 7)
P_j= 2.1*(10**11)                #Radio power of the Jovian decametric component in Watts (Griessmeier 2007 after eq 6)
kHz_to_Hz= 10**3                 #Conversion from KiloHertz to Hertz
MHz_to_Hz= 10**6                 #Conversion from MegaHertz to Hertz
R_j_meters=71492000              #Radius of Jupiter in meters 
 
 
 
#Reshape

N_v=1000
N_B=500
N_nd=100
B_star = np.linspace(1,100,N_B)   # Assumed stellar magnetic field at the surface in gauss (vector)


# Magnetic Moment scaling laws
# Both definitions of the Magnetic moment are scaled with Uj, which is the Magnetic Moment of Jupiter.

 # a) Dipole moment in G cm^3 (Farrel, eq .4) 
Dipole_Magnetic_Moment = (p_rot_j/p_rot)*((M_planet/m_j)**(5/3))*U_j
   
 # b) Planetary Magnetic moment (Stevens, prior to Eq. 15)        
Planetary_Magnetic_Moment = (M_planet/m_j)*(R_planet/R_j)**2*U_j    

Magnetic_field_Equator =  B_j * (Planetary_Magnetic_Moment/U_j)*(R_planet/R_j)**(-3) #At the equator NOT USED 

# fcycl = 2.8 * B   # frequency in MHz, magnetic field in Gauss NOT USED

# b) This estimation of the frequency is used in the Stevens' Model taking into account the equation 8 from Griessmeier paper. 
#   Here we consider the Equatorial Magnetic field.

Freq_Equator = 12*(Planetary_Magnetic_Moment/U_j)*(R_planet/R_j)**(-3) #in MHz, assuming that for Jupiter we have 12 MHz corresponding to the equatorial emission, i.e. 4.3 G
                    

# c)  This estimation of the frequency is used in the Griessmeier and Zarka's Model and is also defined as the Maximum frequency needed to compute the Radio Flux. 
# Is considered the Magnetic field on the pole. That's why is double of the value of Freq_Equator.

Freq_Pole = 2*Freq_Equator        #in MHz. 
Magnetic_field_Pole=Freq_Pole/2.8      #in Gauss

print('Frequency in the equator in MHz: ',Freq_Equator)
print('Frequency in the pole in MHz: ',Freq_Pole)
print('Magnetic field in the equator in G: ',Magnetic_field_Equator)
print('Magnetic field in the pole in G: ',Magnetic_field_Pole)
# Here they are assuming a B=8.6 for Jupiter magnetic field at the pole

# Computing Stellar Wind Density and Stellar Wind Velocity that will be used in Stevens and Griessmeier models.
# The values are based on the Solar wind modeling
# We consider a range of wind and densities, relying on two ages and enlarging by one order of magnitude on each side

# Stellar wind density, Eq.16 Greissmeier 
# In case there is no value of the stellar age, we do a range of possible values. In the case we have a value we can also take ranges around that value in order to compute the stellar wind velocity and density.

n1 = n0*(1+(t1/Tau))**(-1.86)                                 #Wind density at 1AU (particle density) for Stellar Age 1 (Lower value of age)
n_d1 = n1*(A/d0)**(-2)                                        #Stellar wind density at the planetary position for Stellar Age 1 (Lower value of age)

n2 = n0*(1+(t2/Tau))**(-1.86)                                 #Wind density at 1AU (particle density) for Stellar Age 2 (Higher value of age)
n_d2 = n2*(A/d0)**(-2)                                        #Stellar wind density at the planetary position for Stellar Age 2 (Higher value of age)

n_d = np.linspace(n_d2,n_d1,N_nd)                             # range of wind density per cubic meter.

# Wind velocity, eq.15 Griessmeier

# It is supposed that in Stevens' model are using just the Stellar wind velocity and in Griessmeier it is used the Effective Stellar wind velocity (which takes into account the orbital velocity and the stellar wind velocity). We will approximate the Veff=V

v1 = v0*(1+(t1/Tau))**(-0.43)                                #Stellar wind velocity for Stellar Age 1 (Lower value of age)
v2 = v0*(1+(t2/Tau))**(-0.43)                                #Stellar wind velocity for Stellar Age 2 (Higher value of age)

V = np.linspace(v2,v1,N_v)                                   #Range of wind velocity in m/s. 


# 1) FARRELL DESCH MODEL

# Based on W.M. Farrell and M.D. Desch
# It systematically gives very high estimates
# To compute the flux is needed the del_f which is defined as the emission bandwidth. Here it is assumed the emission bandwidth to be 50% of the radiation frequency (Af =  0.5f), which is consistent with the solar system planets. In this model, the frequency used is the Dipole Frequency. 

del_f = 0.5*Freq_Pole             #after eq.8 from Farrel

# Model 1: Considers radio power from Earth, Saturn and Jovian hectometric component only - not visible from ground in w
Pa = ((p_rot_j/p_rot)**0.58)*((M_planet/m_j)**0.98)*((A_j/A)**(-1.17))*P_ja                               #Farrel, eq. 6

#Model 2: Includes radio power from Uranus, Neptune and the Jovian decametric component in W
Pb = ((p_rot_j/p_rot)**0.79)*((M_planet/m_j)**1.33)*((A_j/A)**(-1.60))*P_jb                               #Farrel, eq.7

#Converting power to flux in W m^-2 Hz⁻1, (Farrel eq. 8)
Fluxa = Pa/(4*pi*((s*pc_to_m)**2)*del_f*kHz_to_Hz)
Fluxb = Pb/(4*pi*((s*pc_to_m)**2)*del_f*kHz_to_Hz)

#Flux in mJy 
Fluxa_mJy = Fluxa*Wm2Hz_to_mJy 
Fluxb_mJy = Fluxb*Wm2Hz_to_mJy

# We could plot the flux as a function of the period of the planet (unknown)



# 2) STEVENS MODEL
                                                                                          
#Flux on the basis of Mass loss rate
##In this particular case the free parameters are the wind velocity and the Mass loss rate (MLR) of the star

#log_Ml_s_N = -7.93 + 1.64*math.log(L/L_sun,10) + 0.16*math.log(M_star/M_sun,10) - 1.61*math.log(T_eff,10)                   # H.Nieuwenhuijzen and C.de Jager, good for hotter star
#log_Ml_s_R = -4.74 + 1.50*math.log(L/L_sun,10) - math.log(M_star/M_sun,10) - 2*math.log(T_eff,10)                           # Reiners, good for cooler stars, can be useful to decide a range for SMLR

#Ml_s_N = 10**(log_Ml_s_N)
#Ml_s_R = 10**(log_Ml_s_R)

MLR = np.linspace(0.001,1,200)*MLR_sun                                       #typical mass loss for M dwarfs (ranges can be changed depending on the spectral type of the star)

#Flux in mJy for Stevens model (Stevens, eq.14)
Flux_R = Flux_j_1AU*((MLR.reshape(1,200)/MLR_sun)**(2/3))*((Planetary_Magnetic_Moment/U_j)**(2/3))*((A/5)**(-4/3))*((V.reshape(N_v,1)/V_j)**(5/3))*(((s*pc_to_AU)/1)**(-2))    

#Plot of the flux estimates for Stevens model 

fig = plt.figure(figsize=(11,9))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 


cp = plt.contourf(np.log10(MLR), np.log10(V/1000), np.log10(Flux_R.astype('float64')),30, cmap = 'RdGy_r')
cb = plt.colorbar(cp)

#levels=[-2, -1.3, -1]
ax.contour(np.log10(MLR), np.log10(V/1000), np.log10(Flux_R.astype('float64')), levels=[-1], colors='black', linewidths=2)

#ax.set_title(fontsize=28)
ax.set_xlabel('$Log_{10}(M_{\odot}/yr)$', fontsize=26)
ax.set_ylabel('$Log_{10}$($V_{eff}$[km/s])', fontsize=26)
cb.set_label('$Log_{10}$($\Phi_{mJy}$)', fontsize=26)
ax.tick_params(axis='both', which='major', labelsize=20)
cb.ax.tick_params(axis='both', which='major', labelsize=20)
plt.show()
fig.savefig("FluxEstimates_StevensModel.png", dpi=500)

# 3) GRIESSMEIER AND ZARKA
#In this model, the expected radio flux depends on the source of available energy. In the Kinetic Energy model the input power is proportional to the total kinetic energy of the solar wind impacting into the magnetosphere; and for the Magnetic Energy model is assumed to be proportional to the magnetic energy flux. 

#For both cases it will be needed the value of the Magnetospheric standoff distance. 

R_standoff = 40*R_j*((Planetary_Magnetic_Moment /U_j)**2/((n_d.reshape(1,N_nd)/n_j)*(V.reshape(N_v, 1)/V_j)**2))**(1/6)       #Griessmeier, Eq. 31

#Because the magnetospheric radius cannot be smaller than the planetary radius

R_standoff[R_standoff<R_planet]=R_planet	                                                         


#1. Kinetic Energy Case: 
#The estimation of the flux will be obtained using a range of the effective velocity and the effective density

#Input power for the kinetic energy case in Watts (Griessmeier, eq.1)
P_KE = (n_d.reshape(1,N_nd)/n_j)*(V.reshape(N_v,1)/V_j)**3*(R_standoff/R_sj)**2*P_j              

#Radio flux observed by an observer at a disrance s from the emitter for the KE case (Griessmeier, eq.7)

Flux_KE = P_KE/(beam_rad*(s*pc_to_m)**2*del_f*MHz_to_Hz)       
Flux_KE_mJy = Flux_KE*Wm2Hz_to_mJy			       #Conversion from SI to mJy 

#2. Magnetic Energy Case

# A range of both the effective velocity and the magnetic field of the solar wind is used to obtain an estimate of the flux.

#According to Zarka's SW-Magnetospheric interaction

#Calculation of interplanetary magnetic field for jupiter
# Orbital velocity is needed to compute the value of beta. 
v_orb_j = ((G*(M_sun))/(A_j*AU_to_m))**(0.5) 
v_j = v0*(1+(t_sun/Tau))**(-0.43)

Br_j  = B_r0*(A_j/d0)**(-2)                                 #Radial component of the Interplanetary magnetic field for jupiter in Tesla
B_phi_j = B_phi0*(A_j/d0)**(-1)                             #Azimuthal component of the Interplanetary magnetic field for jupiter in Tesla

alpha_j = math.atan(B_phi_j/Br_j)
beta_j = math.atan(v_orb_j/v_j)
B_j = (Br_j**2 + B_phi_j**2)**(0.5)*(math.sin(alpha_j-beta_j))    #Interplanetary magnetic field for jupiter, used below

print('Interplanetary magnetic field for Jupiter, B_j=',B_j)

                    
#Calculation of interplanetary magnetic field:    #Griessmeier, Eq. 18-22    
#To compute beta are needed the values of v_orb and v:

v_orb = ((G*M_star*(M_sun))/(A*AU_to_m))**(0.5)              #Orbital velocity by Kepler law: v_orb = (GM/A)^ 0.5 (m/s) 
# For the wind density we use the average of stellar ages.
v= v0*(1+(((t1+t2)/2)/Tau))**(-0.43)                         # average value for the velocity (m/s)
                           

# Parker solution of magnetic field (wind stretches and twists the field, which has a strong toroidal and radial component only)
   
B_r = B_r0*(B_star/B_sun)*(A/d0)**(-2)                                   #Radial component of the Interplanetary magnetic field in Tesla
B_phi = (B_star/B_sun) *B_phi0*(A/d0)**(-1)                              #Azimuthal component of the Interplanetary magnetic field in Tesla

alpha = math.atan(B_phi[0]/B_r[0])                        
beta = math.atan(v_orb/v)
Bint = (B_r**2 + B_phi**2)**(0.5)*(math.sin(alpha-beta))     # Interplanetary magnetic field in Tesla
theta= alpha-beta    


#Input power for the magnetic energy case in Watts (Griessmeier, eq.2)
P_ME = ((V.reshape(N_v,1,1)/V_j)*(Bint.reshape(1,1,N_B)/B_j)**2*(R_standoff.reshape(N_v, N_nd, 1)/R_sj)**2)*P_j

#Radio flux observed by an observer at a disrance s from the emitter for the ME case (Griessmeier, eq.7)
Flux_ME = P_ME/(beam_rad*(s*pc_to_m)**2*del_f*MHz_to_Hz)
Flux_ME_mJy = Flux_ME*Wm2Hz_to_mJy                             #Conversion from SI to mJy 
 
 
#Plots for the flux estimates: Griessmeier 
 
# Griessmeier, kinetic contribution            
                                                                      
fig = plt.figure(figsize=(11,9))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 

cp = plt.contourf(np.log10(n_d), np.log10(V/1000), np.log10(Flux_KE_mJy.astype('float64')),30, cmap = 'RdGy_r')
cb = plt.colorbar(cp)

contour_line = ax.contour(np.log10(n_d), np.log10(V/1000), np.log10(Flux_KE_mJy.astype('float64')),
                          levels=[-1], colors='black', linewidths=2)
#ax.clabel(contour_line, inline=True, fontsize=10)


#levels=[-2, -1.3, -1]
#ax.contour(np.log10(n_d), np.log10(V/1000), np.log10(Flux_KE_mJy.astype('float64')))

#ax.set_title('Flux vs $n_{sw}$ and $V_{sw}$')
ax.set_xlabel('$Log_{10}$($n_{sw}$ [$m^{-3}$])', fontsize=26)
ax.set_ylabel('$Log_{10}$($v_{sw}$ [$km/s$])', fontsize=26)
cb.set_label('$Log_{10}$($\Phi_{mJy}$)', fontsize=26)
ax.tick_params(axis='both', which='major', labelsize=20)
cb.ax.tick_params(axis='both', which='major', labelsize=20)

plt.show()
fig.savefig("FluxEstimates_KineticModel.png", dpi=500)      
                                              
# Griessmeier, magnetic contribution.

fig = plt.figure(figsize=(11,9))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 

cp = plt.contourf(B_star, np.log10(V/1000), np.log10(Flux_ME_mJy[:,math.floor(N_nd/2),:].astype('float64')),30,cmap = 'RdGy_r')
cb = plt.colorbar(cp)

#levels=[-2, -1.3, -1]
ax.contour(B_star, np.log10(V/1000), np.log10(Flux_ME_mJy[:,math.floor(N_nd/2),:].astype('float64')),levels=[-1], colors='black', linewidths=2)

#ax.set_title('Flux vs $B_{\perp}$ and $V_{sw}$')
ax.set_xlabel('$B_\star[G]$', fontsize=26) 
ax.set_ylabel('$Log_{10}$($v_{sw}$ [$km/s$])', fontsize=26)
cb.set_label('$Log_{10}$($\Phi_{mJy}$)', fontsize=26)
ax.tick_params(axis='both', which='major', labelsize=20)
cb.ax.tick_params(axis='both', which='major', labelsize=20)
plt.show()

fig.savefig("FluxEstimates_MagneticModel.png", dpi=500)

