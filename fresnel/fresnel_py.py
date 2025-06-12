import numpy as np

# Parameters
D = 40.0                 # Depth in km
velocity = 6           # km/s
frequency = .158         # .05 - .5 Hz..centre frequency: .158 hz (sqrt(fi*f2))
wavelength = velocity / frequency
n = 1

# 1. General Fresnel zone formula (source depth = scatterer depth = D)
d1 = D  # source to scatterer
d2 = D  # scatterer to receiver
r_general = np.sqrt(n * wavelength * (d1 * d2) / (d1 + d2))

# 2. from Gudmundsson: https://academic.oup.com/gji/article/124/1/304/568624
r_vertical = .5*np.sqrt(2 * wavelength * D + (wavelength**2) / 4)

# Print results
print(f"Depth: {D:.2f} km")
print(f"Wavelength (λ): {wavelength:.2f} km")
print(f"General Fresnel radius: {r_general:.2f} km")
print(f"Vertical-incidence Fresnel radius: {r_vertical:.2f} km")
