from fitting_functions import *
from ROOT import TFile

if __name__ == '__main__':
	real_data = TFile('dt.root', 'READ')
	mc_data = TFile('mc.root', 'READ')
	histos = mc_data.Get('MyTask/Histos')
	for histo in histos:
		locals()[histo.GetName()] = histo

	