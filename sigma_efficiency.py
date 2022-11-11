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
canvas = rt.TCanvas()
visualize_n(bin_mass, canvas)
save(canvas, 'binned___')
print(' done')

print('===== STAGE 3 =====\n', 'fitting...') # doesn't work, also we must rename logs, but not now
fit_range = []
for i in range(len(bin_mass)): fit_range.append([1.14, 1.22])

canvas2 = rt.TCanvas()
n = len(bin_mass)
in_raw = 3
in_column = n // in_raw
if in_column * in_raw != n: in_column += 1        
pad = rt.TPad.Divide(canvas2, in_raw, in_column)
for i in range(len(bin_index) - 1):
        my_bkg = rt.TF1('my_bkg', 'pol3', fit_range[i][0], fit_range[i][1], 4)
        my_bkg.SetParameters(6, 0)
        my_bkg.SetLineColor(2)
        canvas2.cd(i + 1)
        visualize(bin_mass[i])
        bin_mass[i].Fit(my_bkg, '', '', fit_range[i][0], fit_range[i][1])
#save(canvas2, 'fitted___1')
        gs = rt.TF1('gs', 'gaus')
        p3 = rt.TF1('p3', 'pol3')
        FullFit_t2 = rt.TF1("FullFit_t2", 'gs + p3', fit_range[i][0], fit_range[i][1], 7)
#        FullFit_t2.Print()
        FullFit_t2.SetParNames()
        FullFit_t2.Print()

       # FullFit_t2.SetParameter(0,1.0);// mean
       # FullFit_t2.SetParLimits(0,0.01,10.);// 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian      
        FullFit_t2.SetParameter(5, 1.19)
        FullFit_t2.SetParLimits(5,1.175,1.21)
        FullFit_t2.SetParameter(6, .003)                                                                                                
        FullFit_t2.SetParLimits(6, 0.0001, .015)
        #FullFit_t2.FixParameter(2, .0001)          
                                                                                                    
        #FullFit_t2.SetParName(0,"Mass")
        #FullFit_t2.SetParName(1,"#sigma")
        #FullFit_t2.SetParName(2,"#Gamma")
        FullFit_t2.SetParameter(0,0)
        FullFit_t2.SetParameter(1,0)
        FullFit_t2.SetParameter(2,0)
        FullFit_t2.SetParameter(3,0)
        #FullFit_t2.SetParameter(8,0)
        #FullFit_t2.SetParameter(9,0)
        #pars = []
        #pars_e = [];

        #for polbin in range(PARBINS-offset): bkg_params[polbin] = my_bkg.GetParameter(polbin)
        for polbin in range(4): 
                FullFit_t2.SetParameter(i, my_bkg.GetParameter(i))
                k = abs(my_bkg.GetParameter(i)) * 0.3
                FullFit_t2.SetParLimits(i, my_bkg.GetParameter(i) - k, my_bkg.GetParameter(i) + k)
        bin_mass[i].Fit(FullFit_t2, '+', '', fit_range[i][0],fit_range[i][1])
save(canvas2, 'fitted___2')

#for(int polbin=0; polbin<PARBINS; polbin++) pars[polbin] = FullFit_t2->GetParameter(polbin);
''' FullFit=new TF1("FullFit",VoigtplusPol,FitRange[0],FitRange[1],PARBINS);
      FullFit->SetParLimits(0,0.01,10.);// 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian                                                              
      FullFit->SetParLimits(1,1.175,1.21);
      if ( kk == 0) FullFit->SetParLimits(1,1.190,1.198);
      //      if ( kk == 6) { FullFit->SetParLimits(1,1.190,1.198);  printf( " kk %d \n",kk);      }                                                             

     // FullFit->SetParLimits(2,0.0001,.015);                                                                                                                    
      FullFit->SetParLimits(2,0.001,.005);

      //        FullFit->FixParameter(3,.0091);// BW width                                                                                                       
      FullFit->FixParameter(3,0.0);
      FullFit->SetParName(0,"Norm");
      FullFit->SetParName(1,"Mass");
      FullFit->SetParName(2,"#sigma");
      FullFit->SetParName(3,"#Gamma");
      FullFit->SetParameter(4,0);
      FullFit->SetParameter(5,0);
      FullFit->SetParameter(6,0);
      FullFit->SetParameter(7,0);
      FullFit->SetParameter(8,0);
      FullFit->SetParameter(9,0);

      for(int polbin=0; polbin<PARBINS; polbin++) FullFit->SetParameter(polbin, pars[polbin]);

      hMassPt[kk]->Fit(FullFit,"IMEQ","",FitRange[0],FitRange[1]);// 3rd fit //IMEQ+              
'''
print(' done')



	
