#!/usr/bin/env python
"""Visualization of distributions of:
    - input variables;
    - deltaR(MC photon, PFCluster) vs pT slices;
    - pfE/mcE vs deltaR slices;
    - MC truth variables.
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

    # ntuples to process
    infiles = fnmatch.filter(os.listdir('input'), '*.root')
    infiles = sorted('input/' + f for f in infiles)

    # text in legends
    txts = [f[f.rfind('gun_') + 4:].replace('.root', '') for f in infiles]

    # make output directories
    for d in ['output', 'output/cache', 'output/plots_inputs']:
        if not os.access(d, os.X_OK):
            os.mkdir(d)

    # draw distributions of inputs
    for det in ['EB', 'EE']:
        r = [make_histos(infile, det) for infile in infiles]

        # repack histograms into per-pileup tuples;
        # order: hE hEta hPhi hR13 hR22 hR25 hR33 hR55 hNVtx hPs1R hPs2R hPs1N hPs2N
        r = list(zip(*r))

        combine(r[0], txts, 'energy_' + det, det, 'Energy (GeV)', 700 if det == 'EE' else 300)
        combine(r[1], txts, 'eta_' + det, det, '#eta')
        combine(r[2], txts, 'phi_' + det, det, '#phi')
        combine(r[3], txts, 'r13_' + det, det, 'E_{1x3}/E')
        combine(r[4], txts, 'r22_' + det, det, 'E_{2x2}/E')
        combine(r[5], txts, 'r25_' + det, det, 'E_{2x5,max}/E')
        combine(r[6], txts, 'r33_' + det, det, 'E_{3x3}/E')
        combine(r[7], txts, 'r55_' + det, det, 'E_{5x5}/E')
        combine(r[8], txts, 'nVtx_' + det, det, 'nVtx')

        # preshower
        if det == 'EE':
            combine(r[9],  txts, 'ps1r_' + det, det, 'E_{ps1}/E', -1, topLegend=True)
            combine(r[10], txts, 'ps2r_' + det, det, 'E_{ps2}/E', -1, topLegend=True)
            combine(r[11], txts, 'ps1n_' + det, det, 'N_{ps1}',  -1, topLegend=True)
            combine(r[12], txts, 'ps2n_' + det, det, 'N_{ps2}',  -1, topLegend=True)

    # draw distributions of deltaR vs pT slices
    for det in ['EB', 'EE']:
        ptPairs = [(0, 1), (1, 10), (10, 20), (20, 100)]
        r = [make_histos_deltaR(infile, det, ptPairs) for infile in infiles]

        # repack histograms into per-pileup tuples
        r = list(zip(*r))

        for ((pt1, pt2), histos) in zip(ptPairs, r):
            cname = 'deltaR_{0}_pT_{1:.1f}_{2:.1f}'.format(det, pt1, pt2)
            fmt = '{0}, {1:.1f} <= p_{{T}}^{{gen}} < {2:.1f} GeV/c^{{2}}'
            title = fmt.format(det, pt1, pt2) + ' ' * 50
            combine(histos, txts, cname, title, '#Delta R(MC photon, PFCluster)', -1, topLegend=True, logY=True)

    # draw distributions of pfE/mcE vs deltaR slices
    for det in ['EB', 'EE']:
        dRPairs = [(0, 0.005), (0, 0.01), (0, 0.02), (0, 0.03), (0, 0.1)]
        r = [make_histos_pfEToMcE(infile, det, dRPairs) for infile in infiles]

        # repack histograms into per-pileup tuples
        r = list(zip(*r))

        for ((dR1, dR2), histos) in zip(dRPairs, r):
            cname = 'pfEToMcE_{0}_dR_{1:.3f}_{2:.3f}'.format(det, dR1, dR2)
            title = '{0}, {1:.3f} <= #Delta R < {2:.3f}'.format(det, dR1, dR2)
            combine(histos, txts, cname, title, 'E^{PF}/E^{gen}', -1, topLegend=True, logY=True)

    # MC truth
    r = [make_histos_mc(infile) for infile in infiles]

    # repack histograms into per-pileup tuples;
    # order: hPtEB hPtEE hEta hPhi hPtZoomEB hPtZoomEE
    r = list(zip(*r))

    # pT
    txtsMC = [t + ', EB' for t in txts] + [t + ', EE' for t in txts]
    combine(r[0] + r[1], txtsMC, 'mcPt',     'MC truth', 'p_{T}^{gen} (GeV/c^{2})')
    combine(r[4] + r[5], txtsMC, 'mcPtZoom', 'MC truth', 'p_{T}^{gen} (GeV/c^{2})', -1, topLegend=True)

    combine(r[2], txts, 'mcEta', 'MC truth', '#eta^{gen}')
    combine(r[3], txts, 'mcPhi', 'MC truth', '#phi^{gen}')

def make_histos(infile, det='EB'):
    """Fills histograms with distributions of PFCluster parameters.

    Results are cached into file.
    """
    # return cached results, if any
    fname = os.path.basename(infile).replace('.root', '')
    cachefile = 'output/cache/draw_inputs_{0}_{1}.pkl'.format(fname, det)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            return pickle.load(f)

    hE    = ROOT.TH1D('h', '', 250, 0, 1000)
    hEta  = ROOT.TH1D('h', '', 150, -3.2, 3.2)
    hPhi  = ROOT.TH1D('h', '', 150, -3.4, 3.4)
    hR13  = ROOT.TH1D('h', '', 200, 0.4, 1.05)
    hR22  = ROOT.TH1D('h', '', 250, 0.7, 1.05)
    hR25  = ROOT.TH1D('h', '', 250, 0.8, 1.05)
    hR33  = ROOT.TH1D('h', '', 250, 0.8, 1.05)
    hR55  = ROOT.TH1D('h', '', 250, 0.9, 1.05)
    hNVtx = ROOT.TH1D('h', '', 60, 0, 60)

    if det == 'EE':
        hPs1R = ROOT.TH1D('h', '', 250, 0, 0.3e-3)
        hPs2R = ROOT.TH1D('h', '', 250, 0, 0.5e-3)
        hPs1N = ROOT.TH1D('h', '', 50, 0, 50)
        hPs2N = ROOT.TH1D('h', '', 50, 0, 50)

    # get TTree with PFClusters
    fi = ROOT.TFile(infile)
    tree = fi.Get('ntuplizer/PFClusterTree')
    if not tree:
        raise Exception('TTree not found')

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

        # for brevity and better performance
        E = tree.pfE

        hE.Fill(E)
        hEta.Fill(tree.pfEta)
        hPhi.Fill(tree.pfPhi)
        hR13.Fill(tree.pfE1x3/E)
        hR22.Fill(tree.pfE2x2/E)
        hR25.Fill(tree.pfE2x5Max/E)
        hR33.Fill(tree.pfE3x3/E)
        hR55.Fill(tree.pfE5x5/E)
        hNVtx.Fill(tree.nVtx)

        if det == 'EE':
            hPs1R.Fill(tree.ps1E/E)
            hPs2R.Fill(tree.ps2E/E)
            hPs1N.Fill(tree.ps1N)
            hPs2N.Fill(tree.ps2N)


    if det == 'EB':
        result = (hE, hEta, hPhi, hR13, hR22, hR25, hR33, hR55, hNVtx)
    else:
        result = (hE, hEta, hPhi, hR13, hR22, hR25, hR33, hR55, hNVtx,
                  hPs1R, hPs2R, hPs1N, hPs2N)

    # normalization
    for h in result:
        h.Scale(1/tree.GetEntriesFast())

    # save cache
    with open(cachefile, 'wb') as f:
        pickle.dump(result, f)

    return result

def combine(histos, txts, cname, title, xtitle, xmax=-1, topLegend=False, logY=False):
    """Visualization of several histograms on single canvas.
    """
    c = ROOT.TCanvas(cname, cname, 700, 700)
    saves.append(c)

    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.08)
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
    frame = c.DrawFrame(histos[0].GetXaxis().GetXmin(), ymin, xmax, ymax * 1.1)

    frame.SetTitle(title)
    frame.SetXTitle(xtitle)
    frame.SetYTitle('Entries/nPhotons')
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.95, 'Y')
    frame.Draw()

    # legend
    if topLegend:
        leg = ROOT.TLegend(0.70, 0.89 - 0.033 * len(txts), 0.91, 0.89)
    else:
        leg = ROOT.TLegend(0.18, 0.18, 0.39, 0.18 + 0.033 * len(txts))

    clrs = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kOrange, ROOT.kGreen, 8] * 10

    for (h, clr, txt) in zip(histos, clrs, txts):
        h.SetLineColor(clr)
        h.Draw('same')
        leg.AddEntry(h, txt, 'l')

    leg.SetFillColor(0)
    leg.Draw('same')

    saves.append((frame, histos, leg))

    c.Update()
    c.SaveAs('output/plots_inputs/distrib_{0}.png'.format(c.GetTitle()))

def make_histos_deltaR(infile, det, ptPairs):
    """Fills histograms with distributions of deltaR vs pT slices.

    Results are cached into file.
    """
    cache = {}

    # load cache, if any
    fname = os.path.basename(infile).replace('.root', '')
    cachefile = 'output/cache/draw_inputs_{0}_deltaR_{1}.pkl'.format(fname, det)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            cache = pickle.load(f)

    # return cached results if the cache exists for all pT slices
    try:
        return [cache[key] for key in ptPairs]
    except KeyError:
        pass

    tofill = {}

    # prepare empty histograms
    for key in ptPairs:
        if key not in cache:
            tofill[key] = ROOT.TH1D('h', '', 200, 0, 0.1)

    # get TTree with PFClusters
    fi = ROOT.TFile(infile)
    tree = fi.Get('ntuplizer/PFClusterTree')
    if not tree:
        raise Exception('TTree not found')

    # fill histograms
    for ev in range(tree.GetEntriesFast()):
        if tree.GetEntry(ev) <= 0:
            raise Exception

        # barrel vs endcaps
        if det == 'EB':
            if abs(tree.pfEta) > 1.479:  # NOTE: pfEta, not mcEta
                continue
        else:
            if abs(tree.pfEta) < 1.479:  # NOTE: pfEta, not mcEta
                continue

        for (pt1, pt2) in tofill:
            if pt1 <= tree.mcPt < pt2:
                tofill[(pt1, pt2)].Fill(tree.pfPhoDeltaR)

    # normalization
    for h in tofill.values():
        h.Scale(1/tree.GetEntriesFast())

    # update cache
    for (key, item) in tofill.items():
        cache[key] = item

    # save updated cache
    with open(cachefile, 'wb') as f:
        pickle.dump(cache, f)

    return [cache[key] for key in ptPairs]

def make_histos_pfEToMcE(infile, det, dRPairs):
    """Fills histograms with distributions of pfE/mcE vs deltaR slices.

    Results are cached into file.
    """
    cache = {}

    # load cache, if any
    fname = os.path.basename(infile).replace('.root', '')
    cachefile = 'output/cache/draw_inputs_{0}_pfEToMcE_{1}.pkl'.format(fname, det)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            cache = pickle.load(f)

    # return cached results if the cache exists for all pT slices
    try:
        return [cache[key] for key in dRPairs]
    except KeyError:
        pass

    tofill = {}

    # prepare empty histograms
    for key in dRPairs:
        if key not in cache:
            tofill[key] = ROOT.TH1D('h', '', 200, 0, 1.5)

    # get TTree with PFClusters
    fi = ROOT.TFile(infile)
    tree = fi.Get('ntuplizer/PFClusterTree')
    if not tree:
        raise Exception('TTree not found')

    # fill histograms
    for ev in range(tree.GetEntriesFast()):
        if tree.GetEntry(ev) <= 0:
            raise Exception

        # barrel vs endcaps
        if det == 'EB':
            if abs(tree.pfEta) > 1.479:  # NOTE: pfEta, not mcEta
                continue
        else:
            if abs(tree.pfEta) < 1.479:  # NOTE: pfEta, not mcEta
                continue

        for (dR1, dR2) in tofill:
            if dR1 <= tree.pfPhoDeltaR < dR2:
                tofill[(dR1, dR2)].Fill(tree.pfE/tree.mcE)

    # normalization
    for h in tofill.values():
        h.Scale(1/tree.GetEntriesFast())

    # update cache
    for (key, item) in tofill.items():
        cache[key] = item

    # save updated cache
    with open(cachefile, 'wb') as f:
        pickle.dump(cache, f)

    return [cache[key] for key in dRPairs]

def make_histos_mc(infile):
    """Fills histograms with distributions of parameters of generated photons.

    Results are cached into file.
    """
    # return cached results, if any
    fname = os.path.basename(infile).replace('.root', '')
    cachefile = 'output/cache/draw_inputs_{0}_mc.pkl'.format(fname)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            return pickle.load(f)

    hPtEB = ROOT.TH1D('h', '', 110, 0, 110)
    hPtEE = ROOT.TH1D('h', '', 110, 0, 110)
    hEta  = ROOT.TH1D('h', '', 250, -3.2, -3.2)
    hPhi  = ROOT.TH1D('h', '', 250, -3.4, 3.4)
    hPtZoomEB = ROOT.TH1D('h', '', 60, 0, 3)
    hPtZoomEE = ROOT.TH1D('h', '', 60, 0, 3)

    # get TTree with PFClusters
    fi = ROOT.TFile(infile)
    tree = fi.Get('ntuplizer/PFClusterTree')
    if not tree:
        raise Exception('TTree not found')

    # fill histograms
    for ev in range(tree.GetEntriesFast()):
        if tree.GetEntry(ev) <= 0:
            raise Exception

        hEta.Fill(tree.mcEta)
        hPhi.Fill(tree.mcPhi)

        # barrel vs endcaps
        if abs(tree.mcEta) < 1.479:
            hPtEB.Fill(tree.mcPt)
            hPtZoomEB.Fill(tree.mcPt)
        else:
            hPtEE.Fill(tree.mcPt)
            hPtZoomEE.Fill(tree.mcPt)

    result = (hPtEB, hPtEE, hEta, hPhi, hPtZoomEB, hPtZoomEE)

    # normalization
    for h in result:
        h.Scale(1/tree.GetEntriesFast())

    # save cache
    with open(cachefile, 'wb') as f:
        pickle.dump(result, f)

    return result


if __name__ == '__main__':
    main()
