#!/usr/bin/env python
"""Visualizes real vs estimated energy resolutions.
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
    ROOT.gROOT.LoadMacro('draw_mva_pars.cc+')

    # ntuples to process
    infiles = fnmatch.filter(os.listdir('input'), '*.root')
    infiles = sorted('input/' + f for f in infiles)

    # names of MVAs
    mva_names = [f[f.rfind('/') + 1:].replace('.root', '') for f in infiles]

    # text in legends
    txts = [f[f.rfind('gun_') + 4:].replace('.root', '') for f in infiles]

    # make output directories
    for d in ['output', 'output/cache', 'output/plots_mva_pars']:
        if not os.access(d, os.X_OK):
            os.mkdir(d)

    for mva_name in mva_names:
        r = [make_graphs(f, mva_name, blockSize=10000) for f in infiles]

        # repack graphs into per-type tuples
        r = list(zip(*r))

        # pT ranges [0, 1), [1, 10) and [10, 100) GeV
        for (i, (pt1, pt2)) in enumerate([(0, 1), (1, 10), (10, 100)]):
            mva = mva_name[mva_name.rfind('gun_') + 4:]
            title = 'Trained on {0}, sliced in {1} < p_{{T}} < {2} GeV/c'.format(mva, pt1, pt2)

            # position, EB
            cname = 'position_EB_{0}_pT{1}-{2}'.format(mva_name, pt1, pt2)
            combine([x[i] for x in r[0]], txts, cname, title, 'Width (expected)', 'Position (real)')

            # position, EE
            cname = 'position_EE_{0}_pT{1}-{2}'.format(mva, pt1, pt2)
            combine([x[i] for x in r[2]], txts, cname, title, 'Width (expected)', 'Position (real)')

            # width, EB
            cname = 'width_EB_{0}_pT{1}-{2}'.format(mva, pt1, pt2)
            combine([x[i] for x in r[1]], txts, cname, title, 'Width (expected)', 'Sigma/Mean (real)')

            # width, EE
            cname = 'width_EE_{0}_pT{1}-{2}'.format(mva, pt1, pt2)
            combine([x[i] for x in r[3]], txts, cname, title, 'Width (expected)', 'Sigma/Mean (real)')

def make_graphs(infile, mva_name, blockSize):
    """Fills, fits and visualizes distributions of Etrue/Erec.

    Results are cached into file.
    """
    # return cached results, if any
    fname = os.path.basename(infile).replace('.root', '')
    fmt = 'output/cache/draw_mva_pars_{0}_{1}_{2}.pkl'
    cachefile = fmt.format(fname, mva_name, blockSize)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            return pickle.load(f)

    # fill necessary arrays of points in C++
    friend = 'output/friend_{0}.root'.format(fname)
    ROOT.fill_arrays(infile, friend, mva_name)

    grMM_EB = []
    grSS_EB = []
    grMM_EE = []
    grSS_EE = []

    # pT ranges [0, 1), [1, 10) and [10, 100) GeV
    for (pt1, pt2) in [(0, 1), (1, 10), (10, 100)]:
        # EB
        title = '{0}_EB_{1}_pT{2}-{3}'.format(fname, mva_name, pt1, pt2)
        ROOT.fit_slices(0, blockSize, title, 'EB, expected width', pt1, pt2)
        grMM_EB.append(ROOT.grMeanVsMean.Clone())
        grSS_EB.append(ROOT.grSigmaVsSigma.Clone())

        # EE
        title = '{0}_EE_{1}_pT{2}-{3}'.format(fname, mva_name, pt1, pt2)
        ROOT.fit_slices(1, blockSize, title, 'EE, expected width', pt1, pt2)
        grMM_EE.append(ROOT.grMeanVsMean.Clone())
        grSS_EE.append(ROOT.grSigmaVsSigma.Clone())

    result = (grMM_EB, grSS_EB, grMM_EE, grSS_EE)

    # save cache
    with open(cachefile, 'wb') as f:
        pickle.dump(result, f)

    return result

def combine(grs, txts, cname, title, xtitle, ytitle, xmax=-1, ymax=-1):
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

    frame = c.DrawFrame(0, 0, xmax, ymax * 1.1)
    frame.SetTitle(title)
    frame.SetXTitle(xtitle)
    frame.SetYTitle(ytitle)
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.95, 'Y')
    frame.Draw()

    # legend
    leg = ROOT.TLegend(0.58, 0.12, 0.89, 0.12 + 0.033 * len(txts))

    clrs = [8, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]

    for (i, (gr, clr, txt)) in enumerate(zip(grs, clrs, txts)):
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
    c.SaveAs('output/plots_mva_pars/{0}.png'.format(c.GetTitle()))


if __name__ == '__main__':
    main()
