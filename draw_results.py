#!/usr/bin/env python
"""Visualizes achieved energy resolutions.
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
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(False)

    # keep current working directory clean from .d and .so files;
    # CPU-intensive part is written in C++ (python is too slow)
    ROOT.gSystem.SetBuildDir('output', True);
    ROOT.gROOT.LoadMacro('draw_results_helper.cc+')

    # ntuples to process
    infiles = fnmatch.filter(os.listdir('input'), '*.root')
    infiles = sorted('input/' + f for f in infiles)

    # names of MVAs
    mva_names = ['mva_mean_' + f[f.rfind('/') + 1:].replace('.root', '') for f in infiles]

    # text in legends
    txts = [f[f.rfind('gun_') + 4:].replace('.root', '') for f in infiles]
    txts = [t + ', no correction' for t in txts] + txts[:]

    # make output directories
    for d in ['output', 'output/cache', 'output/plots_results']:
        if not os.access(d, os.X_OK):
            os.mkdir(d)

    for det in ['EB', 'EE']:
        for mva_branch in mva_names:
            r  = [make_graphs(f, det, '',         blockSize=10000) for f in infiles]
            r += [make_graphs(f, det, mva_branch, blockSize=10000) for f in infiles]

            # repack graphs into per-parameter tuples
            r = list(zip(*r))

            # mean vs E
            title = 'mean_vs_e_{0}_{1}'.format(mva_branch, det)
            combine(r[0], txts, title,           det, 'E^{gen}', 'Mean_{E^{rec}/E^{gen}}')
            combine(r[0], txts, title + '_zoom', det, 'E^{gen}', 'Mean_{E^{rec}/E^{gen}}', False, 20)

            # sigma vs E
            title = 'sigma_vs_e_{0}_{1}'.format(mva_branch, det)
            combineR(r[1], txts, title,           det, 'E^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean')
            combineR(r[1], txts, title + '_zoom', det, 'E^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean', False, 20)

            # mean vs pT
            title = 'mean_vs_pt_{0}_{1}'.format(mva_branch, det)
            combine(r[2], txts, title,           det, 'p_{T}^{gen}', 'Mean_{E^{rec}/E^{gen}}')
            combine(r[2], txts, title + '_zoom', det, 'p_{T}^{gen}', 'Mean_{E^{rec}/E^{gen}}', False, 20)

            # sigma vs pT
            title = 'sigma_vs_pt_{0}_{1}'.format(mva_branch, det)
            combineR(r[3], txts, title,           det, 'p_{T}^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean')
            combineR(r[3], txts, title + '_zoom', det, 'p_{T}^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean', False, 20)

            # mean vs eta
            if det == 'EB': # not necessary to repeat for EE
                title = 'mean_vs_eta_{0}_pT0-1'.format(mva_branch)
                combine([x[0] for x in r[4]], txts, title, '0 < p_{T} < 1 GeV/c', '#eta^{gen}', 'Mean_{E^{rec}/E^{gen}}')

                title = 'mean_vs_eta_{0}_pT1-10'.format(mva_branch)
                combine([x[1] for x in r[4]], txts, title, '1 < p_{T} < 10 GeV/c', '#eta^{gen}', 'Mean_{E^{rec}/E^{gen}}')

                title = 'mean_vs_eta_{0}_pT10-100'.format(mva_branch)
                combine([x[2] for x in r[4]], txts, title, '10 < p_{T} < 100 GeV/c', '#eta^{gen}', 'Mean_{E^{rec}/E^{gen}}')

                # sigma vs eta
                title = 'sigma_vs_eta_{0}_pT0-1'.format(mva_branch)
                combineR([x[0] for x in r[5]], txts, title, '0 < p_{T} < 1 GeV/c', '#eta^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean')

                title = 'sigma_vs_eta_{0}_pT1-10'.format(mva_branch)
                combineR([x[1] for x in r[5]], txts, title, '1 < p_{T} < 10 GeV/c', '#eta^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean')

                title = 'sigma_vs_eta_{0}_pT10-100'.format(mva_branch)
                combineR([x[2] for x in r[5]], txts, title, '10 < p_{T} < 100 GeV/c', '#eta^{gen}', '#sigma_{E^{rec}/E^{gen}}/mean')

            # mean vs nVtx (nopu results excluded)
            title = 'mean_vs_nvtx_{0}_{1}_pT0-1'.format(mva_branch, det)
            combine([x[0] for x in r[6]], txts, title, det, 'nVtx', 'Mean_{E^{rec}/E^{gen}}', True)

            title = 'mean_vs_nvtx_{0}_{1}_pT1-10'.format(mva_branch, det)
            combine([x[1] for x in r[6]], txts, title, det, 'nVtx', 'Mean_{E^{rec}/E^{gen}}', True)

            title = 'mean_vs_nvtx_{0}_{1}_pT10-100'.format(mva_branch, det)
            combine([x[2] for x in r[6]], txts, title, det, 'nVtx', 'Mean_{E^{rec}/E^{gen}}', True)

            # sigma vs nVtx
            title = 'sigma_vs_nvtx_{0}_{1}_pT0-1'.format(mva_branch, det)
            combineR([x[0] for x in r[7]], txts, title, det, 'nVtx', '#sigma_{E^{rec}/E^{gen}}/mean', True)

            title = 'sigma_vs_nvtx_{0}_{1}_pT1-10'.format(mva_branch, det)
            combineR([x[1] for x in r[7]], txts, title, det, 'nVtx', '#sigma_{E^{rec}/E^{gen}}/mean', True)

            title = 'sigma_vs_nvtx_{0}_{1}_pT10-100'.format(mva_branch, det)
            combineR([x[2] for x in r[7]], txts, title, det, 'nVtx', '#sigma_{E^{rec}/E^{gen}}/mean', True)

def make_graphs(infile, det, mva_branch, blockSize):
    """Fills, fits and visualizes distributions of Etrue/Erec.

    Results are cached into file.
    """
    # return cached results, if any
    fname = os.path.basename(infile).replace('.root', '')
    fmt = 'output/cache/draw_results_{0}_{1}_{2}_{3}.pkl'
    cachefile = fmt.format(fname, det, mva_branch, blockSize)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            return pickle.load(f)

    # fill necessary arrays of points in C++
    friend = 'output/friend_{0}.root'.format(fname)
    ROOT.fill_arrays(infile, friend, mva_branch, True if det == 'EE' else False)

    # resolution vs mcE
    title = 'mcE_{0}_{1}_{2}'.format(fname, det, mva_branch)
    ROOT.fit_slices(0, blockSize, title, 'E^{gen}')
    grMeanE = ROOT.grMean.Clone()
    grSigmaE = ROOT.grSigma.Clone()

    # resolution vs mcPt
    title = 'mcPt_{0}_{1}_{2}'.format(fname, det, mva_branch)
    ROOT.fit_slices(1, blockSize, title, 'p_{T}^{gen}')
    grMeanPt = ROOT.grMean.Clone()
    grSigmaPt = ROOT.grSigma.Clone()

    # resolution vs mcEta in pT ranges [0, 1), [1, 10) and [10, 100) GeV
    grMeanEta = []
    grSigmaEta = []
    for (pt1, pt2) in [(0, 1), (1, 10), (10, 100)]:
        title = 'mcEta_{0}_{1}_{2}_pT{3}-{4}'.format(fname, det, mva_branch, pt1, pt2)
        ROOT.fit_slices(2, blockSize, title, '#eta^{gen}', pt1, pt2)
        grMeanEta.append(ROOT.grMean.Clone())
        grSigmaEta.append(ROOT.grSigma.Clone())

    # resolution vs nVtx in pT ranges [0, 1), [1, 10) and [10, 100) GeV
    grMeanVtx = []
    grSigmaVtx = []
    for (pt1, pt2) in [(0, 1), (1, 10), (10, 100)]:
        title = 'nVtx_{0}_{1}_{2}_pT{3}-{4}'.format(fname, det, mva_branch, pt1, pt2)
        ROOT.fit_slices(3, blockSize, title, 'nVtx', pt1, pt2)
        grMeanVtx.append(ROOT.grMean.Clone())
        grSigmaVtx.append(ROOT.grSigma.Clone())

    result = (grMeanE, grSigmaE, grMeanPt, grSigmaPt, grMeanEta, grSigmaEta,
              grMeanVtx, grSigmaVtx)

    # save cache
    with open(cachefile, 'wb') as f:
        pickle.dump(result, f)

    return result

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

    frame = c.DrawFrame(xmin, ymin * 0.9, xmax, ymax * 1.1)
    frame.SetTitle(title)
    frame.SetXTitle(xtitle)
    frame.SetYTitle(ytitle)
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.95, 'Y')
    frame.Draw()

    # number of corrected and uncorrected graphs
    ncorr = len(grs)//2

    # legend
    leg = ROOT.TLegend(0.58, 0.12, 0.89, 0.12 + 0.033 * len(txts))

    clrs = [ROOT.kOrange] * ncorr + [8, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]

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
    c.SaveAs('output/plots_results/{0}.png'.format(c.GetTitle()))

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

    frame = pad1.DrawFrame(xmin, ymin * 0.9, xmax, ymax * 1.1)
    frame.SetTitle(title)
    frame.SetXTitle(xtitle)
    frame.SetYTitle(ytitle)
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.35, 'Y')
    frame.Draw()

    # number of corrected and uncorrected graphs
    ncorr = len(grs)//2

    # legend
    leg = ROOT.TLegend(0.58, 0.58, 0.95, 0.58 + 0.04 * len(txts))

    clrs = [ROOT.kOrange] * ncorr + [8, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]

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
    rats = [ROOT.TGraphErrors() for _ in range(ncorr)]

    for n in range(ncorr):
        gru = grs[n]         # uncorrected
        grc = grs[ncorr + n] # corrected

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

            rats[n].SetPoint(i, x, y1/y2)
            rats[n].SetPointError(i, ex, y1/y2 * ((ey1/y1)**2 + (ey2/y2)**2)**0.5)

    frame = pad2.DrawFrame(xmin, 0.1, xmax, 1.25)
    frame.SetXTitle(xtitle)
    frame.SetYTitle('Corr/uncorr')
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(0.58, 'Y')
    frame.SetTitleSize(0.08, 'XY')
    frame.SetLabelSize(0.08, 'XY')
    frame.Draw()

    # draw corrected/uncorrected ratios
    for (i, (gr, clr)) in enumerate(zip(rats, clrs[ncorr:])):
        if skipNoPU and (i in [0, 1] or i in [ncorr, ncorr + 1]):
           continue

        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(0.4)
        gr.SetMarkerColor(clr)
        gr.SetLineColor(clr)
        gr.Draw('PZL')

    saves.append((frame, rats))

    c.Update()
    c.SaveAs('output/plots_results/{0}.png'.format(c.GetTitle()))


if __name__ == '__main__':
    main()
