import argparse
import numpy as np
from astropy import constants as const


def rho(M,R):
	"""
	Function to calculate the density of a source assuming spherical symmetry

	Arguments:
	M (float) : mass of source, in Msun
	R (float) : radius of source, in parsec

	Returns:
	p (float) : density of source, in kg m^-3
	"""
	p = ((3./4.) * M*const.M_sun.value) / (np.pi*((R*const.pc.value)**3))
	return p


def volume_density(M,R):
	"""
	Function to calculate the volume density of a source assuming spherical symmetry


	"""
	mu = 2.8
	vd = (3.*M*const.M_sun.value)/((mu*const.m_p.value*4.*np.pi)*((R*const.pc.value)**3)*1e6)  #/ 1e6
	return vd


def sigma_therm(T):
	"""
	Function to calculate the thermal broadening in m/s.

	Arguments:
	T (float) : temperature of clump, in K

	Returns:
	sig (float) : thermal velocity dispersion, in m s^-1.
	"""
	mu = 2.8
	sig = np.sqrt( (const.k_B.value*T) / (mu*const.m_p.value) )
	return sig


def jeans_length(M,R,T):
	"""
	Function to calculate the Jeans length

	Arguments:
	M (float) : mass of source, in Msun
	R (float) : radius of source, in parsec
	T (float) : temperature of clump, in K

	Returns:
	lj_m (float) : Jeans length, in metres
	lj_pc (float) : Jeans length, in parsec
	"""
	lj_m = sigma_therm(T) * np.sqrt(np.pi / (const.G.value*rho(M,R)))
	lj_pc = lj_m/const.pc.value
	return lj_m, lj_pc


def jeans_mass(M,R,T):
	"""
	Function to calculate the Jeans mass

	Arguments:
	M (float) : mass of source, in Msun
	R (float) : radius of source, in parsec
	T (float) : temperature of clump, in K

	Returns:
	mj_kg (float) : Jeans mass, in kg
	mj_msun (float) : Jeans mass, in Msun
	"""
	mj_kg = (np.pi**(5./2.) * sigma_therm(T)**3) / (6.*np.sqrt(const.G.value**3 * rho(M,R)))
	mj_msun = mj_kg/const.M_sun.value
	return mj_kg, mj_msun


def ncores(M,MJ):
	"""
	Function to calculate the number of fragments predicted by the competitive accretion model.

	Arguments:
	M (float) : mass of source, in Msun
	MJ (float) : Jeans mass, in Msun

	Returns:
	n (float) : number of cores
	"""	
	CFE = 0.13 # from Palau et al. 2015
	n = (CFE*M)/MJ
	return n


def do_calcs(M,R,T):
	"""
	Main function to do all the calculations and return the results

	Arguments:
	M (float) : mass of source, in Msun
	R (float) : radius of source, in parsec
	T (float) : temperature of clump, in K

	Returns:
	density (float) : density of source, in kgm^-3
	s_th (float) : thermal width, in ms^-1
	L_jeans (float) : jeans length, in pc
	M_jeans (float) : jeans mass, in Msun
	number (float) : Number of predicted fragments
	"""

	###### Results ######
	density = rho(M,R)
	s_th = sigma_therm(T)
	L_jeans = jeans_length(M,R,T)
	M_jeans = jeans_mass(M,R,T)
	number = ncores(M,M_jeans[1])
	vol_dens = volume_density(M,R)

	print("-------  Density, in kgm^-3 ------- ")
	print(density)
	print("")


	print("------- Volume density, in cm^-3 ------- ")
	print(vol_dens)
	print("")


	print("-------  Thermal width, in ms^-1 ------- ")
	print(s_th)
	print("")


	print("-------  Jeans Length, in pc ------- ")
	print(L_jeans[1])
	print("")


	print("-------  Jeans Mass, in Msun ------- ")
	print(M_jeans[1])
	print("")


	print("-------  Predicted number of fragments (CFE=13%) ------- ")
	print(number)
	print("")

	print(" ---- (assumed mu=2.8 and CFE=13%) ----")




if __name__ == '__main__':

	parser=argparse.ArgumentParser(description='''Script to calculate Jeans parameters for any given Mass, Radius and Temperature of source, assuming spherical symmetry.''')
	parser.add_argument("-m","--mass",action='store',dest='mass',type=float,default=100.,help="Mass of source, in Msun")
	parser.add_argument("-r","--radius",action='store',dest='radius',type=float,default=0.1,help="Radius of source, in pc")
	parser.add_argument("-t","--temperature",action='store',dest='temperature',type=float,default=10.,help="Temperature of source, in K")
	args=parser.parse_args()

	print("")
	print("------- Parameters: M={0} Msun, R={1} pc, T={2} K  -------".format(args.mass,args.radius,args.temperature))
	print("")

	# Get results
	do_calcs(args.mass,args.radius,args.temperature)

