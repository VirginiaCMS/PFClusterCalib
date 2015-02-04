#!/usr/bin/env python
"""Visualization of some slices with achieved energy resolutions.
"""

# python-2 compatibility
from __future__ import division        # 1/2 = 0.5, not 0
from __future__ import print_function  # print() syntax from python-3

import os
import pickle
import fnmatch
import ROOT

# for keeping drawed objects in memory
saves = []

def main():
    """Steering function.
    """
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(False)

    # ntuples to process
    infiles = fnmatch.filter(os.listdir('input'), '*.root')
    infiles = sorted('input/' + f for f in infiles)

    # suffixes of ntuples
    sfxs = [f[f.rfind('gun_') + 4:].replace('.root', '') for f in infiles]

    # names of ntuples
    fnames = [f[f.rfind('/') + 1:].replace('.root', '') for f in infiles]

    # names of branches with MVA outputs to take
    mvas = [''] + ['mva_mean_' + f for f in fnames]  # '' = no correction

    # text in legends
    txts  = ['no correction'] + sfxs

    # custom [true pt]+[reco eta] slices (EB + EE)
    slicesEB = [(0.3, 0.6, 0, 1.479), (0.6, 1, 0, 1.479), (1, 10, 0, 1.479), (10, 100, 0, 1.479)]
    slicesEE = [(0.3, 0.6, 1.479, 5), (0.6, 1, 1.479, 5), (1, 10, 1.479, 5), (10, 100, 1.479, 5)]

    # (rebin, xmin, xmax) tuples for better drawing below
    customEB = [(10, 0, 2), (10, 0, 2), (4, 0.6, 1.5), (1, 0.85, 1.2)]
    customEE = [(25, 0, 2), (20, 0, 2), (4, 0.6, 1.4), (1, 0.75, 1.3)]

    # make output directories
    for d in ['output', 'output/cache', 'output/plots_slices']:
        if not os.access(d, os.X_OK):
            os.mkdir(d)

    # evaluate/draw shapes of pfE/mcE distributions per ntuple for EB/EE
    for (slices, det, custom) in zip([slicesEB, slicesEE], ['EB', 'EE'], [customEB, customEE]):
        for (infile, sfx) in zip(infiles, sfxs):
            r = [make_histos(infile, slices, mva) for mva in mvas]

            # repack histograms into per-slice tuples
            r = list(zip(*r))

            for (histos, (pt1, pt2, _, _), (rebin, xmin, xmax)) in zip(r, slices, custom):
                cname = '{0}_{1}_pT_{2:.1f}_{3:.1f}'.format(sfx, det, pt1, pt2)
                fmt = '{0}, {1}, {2:.1f} <= p_{{T}}^{{gen}} < {3:.1f} GeV/c^{{2}}'
                title = ' ' * 20 + 'test sample: ' + fmt.format(sfx, det, pt1, pt2)

                combine(histos, txts, cname, title, 'correction * E^{PF}/E^{true}', rebin, xmin, xmax)

    # save all open canvases as images
    canvases = ROOT.gROOT.GetListOfCanvases()
    for i in range(canvases.GetEntries()):
        c = canvases.At(i)
        c.SaveAs('output/plots_slices/slices_{0}.png'.format(c.GetTitle()))

def make_histos(infile, regions, mva_branch_name):
    """Fills requested energy resolution histograms.

    Results are cached into file.

    regions = list of (true pt1, true pt2, reco |eta1|, reco |eta2|) tuples.

    Empty mva_branch_name implies filling histograms without energy corrections.
    """
    cache = {}

    # load cache, if any
    fname = os.path.basename(infile).replace('.root', '')
    cachefile = 'output/cache/draw_slices_{0}.pkl'.format(fname)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            cache = pickle.load(f)

    if mva_branch_name not in cache:
        cache[mva_branch_name] = {}

    cache1 = cache[mva_branch_name]

    # return cached results if the cache exists for all pt+eta slices
    try:
        return [cache1[key] for key in regions]
    except KeyError:
        pass

    tofill = {}

    # prepare empty histograms
    for key in regions:
        if key not in cache1:
            tofill[key] = ROOT.TH1D('h', '', 1000, 0, 2)

    # get TTree with PFClusters
    fi = ROOT.TFile(infile)
    tree = fi.Get('ntuplizer/PFClusterTree')
    if not tree:
        raise Exception('TTree not found')

    # add branches with outputs from MVAs
    tree.AddFriend('ntuplizer/PFClusterTree', 'output/friend_{0}.root'.format(fname))

    # fill histograms
    for ev in range(1, tree.GetEntriesFast(), 2):  # test events only
        if tree.GetEntry(ev) <= 0:
            raise Exception

        # value of correction
        if mva_branch_name:
            corr = getattr(tree, mva_branch_name)
        else:
            corr = 1

        pt    = tree.mcPt  # NOTE: true generated pT
        aeta  = abs(tree.pfEta)
        ratio = corr * tree.pfE / tree.mcE

        for (pt1, pt2, aeta1, aeta2) in tofill:
            if pt1 <= pt < pt2 and aeta1 <= aeta < aeta2:
                tofill[(pt1, pt2, aeta1, aeta2)].Fill(ratio)

    # normalization
    for h in tofill.values():
        h.Scale(1/tree.GetEntriesFast())

    # update cache
    for (key, item) in tofill.items():
        cache1[key] = item

    # save updated cache
    with open(cachefile, 'wb') as f:
        pickle.dump(cache, f)

    return [cache1[key] for key in regions]

def combine(histos, txts, cname, title, xtitle, rebin=1, xmin=0, xmax=-1, logY=False):
    """Visualization of several histograms on single canvas.
    """
    # rebin histograms, if requested
    if rebin > 1:
        for h in histos:
            h.Rebin(rebin)

    c = ROOT.TCanvas(cname, cname, 700, 700)
    saves.append(c)

    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.06)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()

    if logY:
        c.SetLogy()

    if xmax < 0:
        xmax = histos[0].GetXaxis().GetXmax()

    # draw empty histogram
    ymax = max(h.GetMaximum() for h in histos)
    ymin = 0 if not logY else ymax * 10e-4
    frame = c.DrawFrame(xmin, ymin, xmax, ymax * 1.1)

    frame.SetTitle(title)
    frame.SetXTitle(xtitle)
    frame.SetYTitle('Entries/nPhotons')
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.95, 'Y')
    frame.Draw()

    # legend
    leg = ROOT.TLegend(0.63, 0.89 - 0.033 * len(txts), 0.94, 0.89)

    clrs = [ROOT.kOrange, ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen, 8, 42] * 10

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
