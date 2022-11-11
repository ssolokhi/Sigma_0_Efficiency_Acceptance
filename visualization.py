import ROOT as rt


def visualize(histo, option = ''):
	histo.SetMarkerStyle(21)
	histo.SetMarkerSize(.5)
	histo.SetLineColor(rt.kBlue)
	histo.Draw(option)

def visualize_n(histos, canvas, in_raw = 3, option = ''):
	n = len(histos)
	in_column = n // in_raw
	if in_column * in_raw != n: in_column += 1
	pad = rt.TPad.Divide(canvas, in_raw, in_column)
	for i in range(n):
       		canvas.cd(i + 1)
       		visualize(histos[i])

def save(obj, name):
	obj.Draw()
	obj.Write(name)

