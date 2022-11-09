from fitting_functions import *
from visualization import *
import ROOT as rt


print('===== STAGE 1 =====\n', 'Reading data from files...')
real_data = rt.TFile('dt.root', 'READ')
mc_data = rt.TFile('mc.root', 'READ')
histos = mc_data.Get('MyTask/Histos')
for histo in histos:
	locals()[histo.GetName()] = histo
out = rt.TFile('out.root', 'RECREATE')
print(' done')

print('===== STAGE 2 =====\n', 'binning...')
pt_edges_sigma = [1.1, 1.6, 2.0, 2.4, 2.8, 3.4, 4.0, 5.0, 8.0, 24.0]
bin_index = [int(pt_edges_sigma[i] * 10) + 1 for i in range(len(pt_edges_sigma))]
bin_mass = []
nptbin = len(bin_index) - 1
print(*bin_index)
for i in range(len(bin_index) - 1):
	local_name = 'hSigmaMass' + '_' + str(pt_edges_sigma[i]) + '_' + str(pt_edges_sigma[i + 1])
	bin_mass.append(hLamGv0.ProjectionX(local_name, bin_index[i], bin_index[i + 1]))
visualize_n(bin_mass, 'binned___')
print(' done')

print('===== STAGE 3 =====\n', 'fisting...') # doesn't work, also we must rename logs, but not now
fitFcn = rt.TF1('fitFcn', gaussian_and_polynome, 1.35, 1.22, 7)
print(fitFcn)
fitFcn.SetParameters(1., 1., 1., 1., 1., 1., 1.)
for i in range(len(bin_index) - 1):
        bin_mass[i].Fit('fitFcn')
visualize_n(bin_mass, 'fitted___')
print(' done')



	
