import ROOT as rt


def visualize(histo, option = ''):
	histo.SetMarkerStyle(21)
	histo.SetMarkerSize(.5)
	histo.SetLineColor(rt.kBlue)
	histo.Draw(option)

def visualize_n(histos, name, in_raw = 3, option = ''):
	c = rt.TCanvas()
	n = len(histos)
	in_column = n // in_raw
	if in_column * in_raw != n: in_column += 1
	pad = rt.TPad.Divide(c, in_raw, in_column)
	for i in range(n):
        	c.cd(i + 1)
        	visualize(histos[i])
	c.Write(name)

