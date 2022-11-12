from fitting_functions import *
from visualization import *
import ROOT as rt


print('========= STAGE 1 =========\n', 'Reading data from files...')
real_data = rt.TFile('dt.root', 'READ')
mc_data = rt.TFile('mc.root', 'READ')
histos = mc_data.Get('MyTask/Histos')
for histo in histos:
	locals()[histo.GetName()] = histo
out = rt.TFile('out.root', 'RECREATE')
print(' done')

hLamGv0.Write('hLamGv0')

print('========= STAGE 2 =========\n', 'binning...')
pt_edges_sigma = [1.1, 1.6, 2.0, 2.4, 2.8, 3.4, 4.0, 5.0, 8.0, 24.0]
for i in range(len(pt_edges_sigma)): 
	pt_edges_sigma[i] /= 2
bin_index = [int(pt_edges_sigma[i] * 10) + 1 for i in range(len(pt_edges_sigma))]
hMass = []
nptbin = len(bin_index) - 1
print(*bin_index)
local_names = []
for i in range(len(bin_index) - 1):
	local_name = 'hSigmaMass' + '_' + str(pt_edges_sigma[i] * 2) + '_' + str(pt_edges_sigma[i + 1] * 2)
	local_names.append(local_name)
	hMass.append(hLamGv0.ProjectionX(local_name, bin_index[i], bin_index[i + 1]))
	hMass[-1].Write(local_name)
c = rt.TCanvas()
in_raw = 3
in_column = len(hMass) // in_raw
if in_column * in_raw != len(hMass): in_column += 1
pad = rt.TPad.Divide(c, in_raw, in_column)
for i in range(len(hMass)): 
	c.cd(i + 1)
	hMass[i].Draw()
c.Draw()
c.Write('binned_pad')
print(' done')

print('========= STAGE 3 =========\n', 'bkg pre-fitting...') 
fit_range = []
for i in range(len(hMass)): fit_range.append([1.14, 1.22])
bkg_pars = []
c2 = rt.TCanvas()
pad = rt.TPad.Divide(c2, in_raw, in_column)

for i in range(len(hMass)):
	c2.cd(i + 1)
	bkg = rt.TF1( 'bkg', 'pol3', fit_range[i][0], fit_range[i][1], 4)
	bkg.SetParameters(6,0)
	bkg.SetLineColor(2)
	hMass[i].Draw()
	hMass[i].Fit(bkg, 'Q', '', fit_range[i][0], fit_range[i][1]) #set NQ in options to not draw and not show calculations, respectfully 
	hMass[i].Write(local_names[i] + '_bg_fit')
	bkg_pars.append(bkg.GetParameters())
print(*bkg_pars, sep = '\n')
c2.Draw()
c2.Write('bkg_fitted_pad')
print(' done')

print('========= STAGE 4 =========\n', 'fitting...') 
fit_pars = []
c3 = rt.TCanvas()
pad = rt.TPad.Divide(c3, in_raw, in_column)
fit = '[3] + [4] * x + [5] * x * x + [6] * x * x * x + [0] * exp(-0.5 * ((x - [1]) / [2]) * ((x - [1]) / [2]))'
for i in range(len(hMass)):
	c3.cd(i + 1)
	hMass[i].Draw()
	fullfit = rt.TF1( 'fullfit', fit, fit_range[i][0], fit_range[i][1])
	fullfit.Print()
	fullfit.SetLineColor(2)
	fullfit.SetParLimits(1,1.175,1.21)
	fullfit.SetParLimits(2,0.001, 0.03)
	for j in range(4):
		fullfit.SetParameter(j + 3, bkg_pars[i][j])
	hMass[i].Fit(fullfit, 'I', '', fit_range[i][0], fit_range[i][1]) #set NQ in options to not draw and not show calculations, respectfully 
	hMass[i].Write(local_names[i] + '_fit')
	fit_pars.append(fullfit.GetParameters())
c3.Draw()
c3.Write('fitted_pad')
print(' done')