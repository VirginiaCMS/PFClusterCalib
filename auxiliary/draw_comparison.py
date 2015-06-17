#!/usr/bin/env python
"""Makes comparison of different trainings.
"""

# python-2 compatibility
from __future__ import division        # 1/2 = 0.5, not 0
from __future__ import print_function  # print() syntax from python-3

import os
import pickle
import fnmatch
import ROOT

# for keeping drawed ROOT objects in memory
saves = []

def main():
    """Steering function.
    """
    # configuration
    mva1 = 'prev_Feb15'
    mva2 = 'new_Apr15'
    blockSize = 3000
    eregions = [(0, 1), (1, 2), (2, 10), (10, 20), (20, 100), (100, 1000)]

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(False)

    # make output directory
    os.mkdir('plots_comparison')

    # ntuples
    ntuples = fnmatch.filter(os.listdir('/afs/cern.ch/work/k/konush/public/PFClusterCalib_Apr2015/input'), '*.root')
    ntuples = sorted('/afs/cern.ch/work/k/konush/public/PFClusterCalib_Apr2015/input/' + f for f in ntuples)

    # MVAs
    mvas = [f[f.rfind('/') + 1:].replace('.root', '') for f in ntuples]

    for det in ['EB', 'EE']:
        for mva in mvas:
            for (i, ntuple) in enumerate(ntuples):
                r  = [get_graphs(mva1, ntuple, det, '',                blockSize)]
                r += [get_graphs(mva1, ntuple, det, 'mva_mean_' + mva, blockSize)]
                r += [get_graphs(mva2, ntuple, det, 'mva_mean_' + mva, blockSize)]

                # repack graphs into per-parameter tuples
                r = list(zip(*r))

                ntuple = ntuple[ntuple.rfind('gun_') + 4:].replace('.root', '')

                # text in legends
                txts  = [ntuple + ', no correction', ntuple + ', ' + mva1, ntuple + ', ' + mva2]

                # text in captions
                cap = 'trained on {0}, {1}'.format(mva[mva.rfind('gun_') + 4:], det)

                # mean vs E
                title = 'mean_vs_e_{0}_{1}_{2}'.format(mva, ntuple, det)
                combine(r[0], txts, title,           cap, 'E^{gen}', 'Mean_{E^{rec}/E^{gen}}')
                combine(r[0], txts, title + '_zoom', cap, 'E^{gen}', 'Mean_{E^{rec}/E^{gen}}', False, 20)

                # sigma vs E
                title = 'sigma_vs_e_{0}_{1}_{2}'.format(mva, ntuple, det)
                combineR(r[1], txts, title,           cap, 'E^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean')
                combineR(r[1], txts, title + '_zoom', cap, 'E^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean', False, 20)

                # mean vs pT
                title = 'mean_vs_pt_{0}_{1}_{2}'.format(mva, ntuple, det)
                combine(r[2], txts, title,           cap, 'p_{T}^{gen}', 'Mean_{E^{rec}/E^{gen}}')
                combine(r[2], txts, title + '_zoom', cap, 'p_{T}^{gen}', 'Mean_{E^{rec}/E^{gen}}', False, 20)

                # sigma vs pT
                title = 'sigma_vs_pt_{0}_{1}_{2}'.format(mva, ntuple, det)
                combineR(r[3], txts, title,           cap, 'p_{T}^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean')
                combineR(r[3], txts, title + '_zoom', cap, 'p_{T}^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean', False, 20)

                for (i, (e1, e2)) in enumerate(eregions):
                    cap = 'trained on {0}, {1} < E^{{gen}} < {2} GeV/c'.format(mva[mva.rfind('gun_') + 4:], e1, e2)

                    # mean and sigma vs eta
                    if det == 'EB': # not necessary to repeat for EE
                        title = 'mean_vs_eta_{0}_{1}_E{2}-{3}'.format(mva, ntuple, e1, e2)
                        combine([x[i] for x in r[4]], txts, title, cap, '#eta^{gen}', 'Mean_{E^{rec}/E^{gen}}')

                        # sigma vs eta
                        title = 'sigma_vs_eta_{0}_{1}_E{2}-{3}'.format(mva, ntuple, e1, e2)
                        combineR([x[i] for x in r[5]], txts, title, cap, '#eta^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean')

    os.rename('plots_comparison', 'plots_comparison_{0}_{1}'.format(mva1, mva2))

def get_graphs(ttype, infile, det, mva_branch, blockSize):
    """Returns cached results.
    """
    fname = os.path.basename(infile).replace('.root', '')
    fmt = 'output_{0}/cache/draw_results_{1}_{2}_{3}_{4}.pkl'
    cachefile = fmt.format(ttype, fname, det, mva_branch, blockSize)

    with open(cachefile, 'rb') as f:
        return pickle.load(f)

def combine(grs, txts, cname, title, xtitle, ytitle, skipNoPU=False, xmax=-1, ymax=-1):
    """Visualization of several graphs on single canvas.
    """
    c = ROOT.TCanvas(cname, cname, 700, 700)
    saves.append(c)

    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.08)
    c.SetTopMargin(0.06)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()

    # draw dummy histogram
    xmin = min(gr.GetX()[i] for gr in grs for i in range(gr.GetN()))
    ymin = min(gr.GetY()[i] for gr in grs for i in range(gr.GetN()))
    if xmax < 0:
        xmax = max(gr.GetX()[i] for gr in grs for i in range(gr.GetN()))
    if ymax < 0:
        ymax = max(gr.GetY()[i] for gr in grs for i in range(gr.GetN()))

    frame = c.DrawFrame(0 if xmin >= 0 else xmin, ymin * 0.9, xmax, ymax * 1.1)
    frame.SetTitle(title)
    frame.SetXTitle(xtitle)
    frame.SetYTitle(ytitle)
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.95, 'Y')
    frame.Draw()

    # number of corrected and uncorrected graphs
    ncorr = len(grs)//3

    # legend
    leg = ROOT.TLegend(0.58, 0.12, 0.89, 0.12 + 0.033 * len(txts))

    clrs = [ROOT.kOrange] * ncorr + [ROOT.kBlack, ROOT.kRed, 8]

    for (i, (gr, clr, txt)) in enumerate(zip(grs, clrs, txts)):
        if skipNoPU and (i in [0, 1] or i in [ncorr, ncorr + 1]):
           continue

        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(0.4)
        gr.SetMarkerColor(clr)
        gr.SetLineColor(clr)

        gr.Draw('PZL')
        leg.AddEntry(gr, txt, 'p')

    leg.SetFillColor(0)
    leg.Draw('same')

    saves.append((frame, grs, leg))

    c.Update()
    c.SaveAs('plots_comparison/{0}.png'.format(c.GetTitle()))

def combineR(grs, txts, cname, title, xtitle, ytitle, skipNoPU=False, xmax=-1, ymax=-1):
    """Visualization of graphs along with corrected/uncorrected ratios.
    """
    c = ROOT.TCanvas(cname, cname, 700, 700)
    saves.append(c)

    pad1 = ROOT.TPad('top', 'top', 0, 0.3, 1, 1)
    pad2 = ROOT.TPad('bottom', 'bottom', 0, 0, 1, 0.3)
    pad1.Draw()
    pad2.Draw()
    saves.append(pad1)
    saves.append(pad2)

    # top pad
    pad1.cd()
    pad1.SetTopMargin(0.07)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.095)
    pad1.SetRightMargin(0.02)
    pad1.SetGridx()
    pad1.SetGridy()

    # draw dummy histogram
    xmin = min(gr.GetX()[i] for gr in grs for i in range(gr.GetN()))
    ymin = min(gr.GetY()[i] for gr in grs for i in range(gr.GetN()))
    if xmax < 0:
        xmax = max(gr.GetX()[i] for gr in grs for i in range(gr.GetN()))
    if ymax < 0:
        ymax = max(gr.GetY()[i] for gr in grs for i in range(gr.GetN()))

    frame = pad1.DrawFrame(0 if xmin >= 0 else xmin, ymin * 0.9, xmax, ymax * 1.1)
    frame.SetTitle(title)
    frame.SetXTitle(xtitle)
    frame.SetYTitle(ytitle)
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.35, 'Y')
    frame.Draw()

    # number of corrected and uncorrected graphs
    ncorr = len(grs)//3

    # legend
    leg = ROOT.TLegend(0.58, 0.58, 0.95, 0.58 + 0.04 * len(txts))

    clrs = [ROOT.kOrange] * ncorr + [ROOT.kBlack, ROOT.kRed, 8]

    for (i, (gr, clr, txt)) in enumerate(zip(grs, clrs, txts)):
        if skipNoPU and (i in [0, 1] or i in [ncorr, ncorr + 1]):
           continue

        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(0.4)
        gr.SetMarkerColor(clr)
        gr.SetLineColor(clr)

        gr.Draw('PZL')
        leg.AddEntry(gr, txt, 'p')

    leg.SetFillColor(0)
    leg.Draw('same')

    saves.append((frame, grs, leg))

    # bottom pad
    pad2.cd()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.25)
    pad2.SetLeftMargin(0.095)
    pad2.SetRightMargin(0.02)
    pad2.SetGridx()
    pad2.SetGridy()

    # make corrected/uncorrected ratios
    rats = [[ROOT.TGraphErrors() for _ in range(ncorr)] for _ in range(2)]

    for t in range(2):
        for n in range(ncorr):
            gru = grs[n]         # uncorrected
            grc = grs[ncorr * (t + 1) + n] # corrected

            for i in range(grc.GetN()):
                x  = grc.GetX()[i]
                y1 = grc.GetY()[i]
                ex  = grc.GetEX()[i]
                ey1 = grc.GetEY()[i]

                # NOTE: gru may have different x and ex, but we ignore this
                # difference. Instead, y for gru is estimated from linear
                # extrapolation between corresponding nearest points.

                # find nearest points
                for j in range(gru.GetN() - 1):
                    x1 = gru.GetX()[j]
                    x2 = gru.GetX()[j + 1]

                    if  x1 <= x < x2:
                        a = (gru.GetY()[j] - gru.GetY()[j + 1]) / (x1 - x2)
                        b = gru.GetY()[j] - a * x1
                        y2 = a * x + b

                        ea = (gru.GetEY()[j] - gru.GetEY()[j + 1]) / (x1 - x2)
                        eb = gru.GetEY()[j] - ea * x1
                        ey2 = ea * x + eb

                        break
                else:
                    if x < gru.GetX()[0]:
                        y2 = gru.GetY()[0]
                        ey2 = gru.GetEY()[0]
                    else:
                        y2 = gru.GetY()[gru.GetN() - 1]
                        ey2 = gru.GetEY()[gru.GetN() - 1]

                rats[t][n].SetPoint(i, x, y1/y2)
                rats[t][n].SetPointError(i, ex, y1/y2 * ((ey1/y1)**2 + (ey2/y2)**2)**0.5)

    frame = pad2.DrawFrame(xmin, 0.1, xmax, 1.25)
    frame.SetXTitle(xtitle)
    frame.SetYTitle('Corr/uncorr')
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(0.58, 'Y')
    frame.SetTitleSize(0.08, 'XY')
    frame.SetLabelSize(0.08, 'XY')
    frame.Draw()

    # draw corrected/uncorrected ratios
    for t in range(2):
        for (i, (gr, clr)) in enumerate(zip(rats[t], clrs[ncorr * (t + 1):])):
            if skipNoPU and (i in [0, 1] or i in [ncorr, ncorr + 1]):
                continue

            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.4)
            gr.SetMarkerColor(clr)
            gr.SetLineColor(clr)
            gr.Draw('PZL')

    saves.append((frame, rats))

    c.Update()
    c.SaveAs('plots_comparison/{0}.png'.format(c.GetTitle()))


if __name__ == '__main__':
    main()
