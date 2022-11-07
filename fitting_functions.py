from math import sqrt, pi, exp, log

def polinome(x: float, params: list) -> float:
	return params[0] + params[1]*x + params[2]*(x**2) + params[3]*(x**3)

def gaussian(x: float, params: list) -> float:
	return params[0]/(sqrt(2*pi)*params[2])*exp(-((x - params[1])**2)/(2*params[2]))

def gaussian_and_polinome(x: float, params: list) -> float:
	intergal, mu, sigma, *polinome_params = params
	return gaussian(x, [intergal, mu, sigma]) + polinome(x, polinome_params)

def get_rapidity(pt: float, pz: float, mass: float) -> float:
	energy = sqrt(mass**2 + pt**2 + pt**2)
	rapidity = 0.5*log((energy + pz)/(energy - pz))
	return rapidity