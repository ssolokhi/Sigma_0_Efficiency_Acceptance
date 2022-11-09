from fitting_functions import *
from ROOT import TFile, TCanvas

if __name__ == '__main__':
	real_data = TFile('dt.root', 'READ')
	mc_data = TFile('mc.root', 'READ')
	histos = mc_data.Get('MyTask/Histos')
	for histo in histos:
		locals()[histo.GetName()] = histo
		print(histo)




	output_file = TFile('output.root', 'RECREATE')

	output_file.Close()

	print('Task Done!\n')