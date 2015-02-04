#!/usr/bin/env python
"""Overtraining tests for the semi-parametric MVAs.
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
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(False)

    # ntuples to process
    infiles = fnmatch.filter(os.listdir('input'), '*.root')
    infiles = sorted('input/' + f for f in infiles)

    # names of ntuples
    fnames = [f[f.rfind('/') + 1:].replace('.root', '') for f in infiles]

    # make output directories
    for d in ['output', 'output/cache', 'output/plots']:
        if not os.access(d, os.X_OK):
            os.mkdir(d)

    for det in ['EB', 'EE']:
        for (infile, fname) in zip(infiles, fnames):
            r = make_histos(infile, det)

            cname = 'overtraining_{0}_{1}'.format(fname, det)
            combine(r, cname, det)

    # save all open canvases as images
    canvases = ROOT.gROOT.GetListOfCanvases()
    for i in range(canvases.GetEntries()):
        c = canvases.At(i)
        c.SaveAs('output/plots/{0}.png'.format(c.GetTitle()))

def make_histos(infile, det):
    """Fills energy resolution histograms for train and test trees.

    Results are cached into file.
    """
    # return cached results, if any
    fname = os.path.basename(infile).replace('.root', '')
    cachefile = 'output/cache/draw_overtraining_{0}_{1}.pkl'.format(fname, det)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            return pickle.load(f)

    hTrain = ROOT.TH1D('h', '', 500, 0., 1.2)
    hTest  = ROOT.TH1D('h', '', 500, 0., 1.2)
    hOrig  = ROOT.TH1D('h', '', 500, 0., 1.2)

    # get TTree with PFClusters
    fi = ROOT.TFile(infile)
    tree = fi.Get('ntuplizer/PFClusterTree')
    if not tree:
        raise Exception('TTree not found')

    # add branches with outputs from MVAs
    tree.AddFriend('ntuplizer/PFClusterTree', 'output/friend_{0}.root'.format(fname))

    # fill histograms
    for ev in range(tree.GetEntriesFast()):
        if tree.GetEntry(ev) <= 0:
            raise Exception

        # barrel vs endcaps
        if det == 'EB':
            if abs(tree.pfEta) > 1.479:
                continue
        else:
            if abs(tree.pfEta) < 1.479:
                continue

        #if tree.pfSize5x5_ZS != 2:
            #continue

        orig = tree.pfE/tree.mcE
        corr = getattr(tree, 'mva_mean_' + fname)

        if ev % 2 == 0:
            hTrain.Fill(corr * orig)
        else:
            hTest.Fill(corr * orig)

        hOrig.Fill(orig)

    # normalization
    hOrig.Scale(0.5)

    result = (hTrain, hTest, hOrig)

    # save cache
    with open(cachefile, 'wb') as f:
        pickle.dump(result, f)

    return result

def combine(histos, cname, det):
    """Visualization of train/test/original distributions on single canvas.
    """
    c = ROOT.TCanvas(cname, cname, 700, 700)
    saves.append(c)

    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.08)
    c.SetTopMargin(0.06)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()

    # draw empty histogram
    xmin = histos[0].GetXaxis().GetXmin()
    xmax = histos[0].GetXaxis().GetXmax()
    ymax = max(h.GetMaximum() for h in histos)
    frame = c.DrawFrame(xmin, 0, xmax, ymax * 1.1)

    frame.SetTitle('Overtraining test, ' + det)
    frame.SetXTitle('correction * E^{PF}/E^{gen}')
    frame.SetYTitle('Entries')
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.95, 'Y')
    frame.Draw()

    # legend
    leg = ROOT.TLegend(0.6, 0.79, 0.91, 0.89)

    clrs = [ROOT.kBlack, ROOT.kBlue, ROOT.kOrange]
    txts = ['Corrections from train', 'Corrections from test',
            'No corrections, test+train']

    for (h, clr, txt) in zip(histos, clrs, txts):
        h.SetLineColor(clr)
        h.Draw('same')
        leg.AddEntry(h, txt, 'l')

    leg.SetFillColor(0)
    leg.Draw('same')

    saves.append((frame, histos, leg))

    c.Update()


if __name__ == '__main__':
    main()
