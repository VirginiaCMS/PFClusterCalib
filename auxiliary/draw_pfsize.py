#!/usr/bin/env python
"""Visualization of fractions of various PFCluster sizes vs pT.
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

    # make output directories
    for d in ['output', 'output/cache', 'output/plots']:
        if not os.access(d, os.X_OK):
            os.mkdir(d)

    for det in ['EB', 'EE']:
        for infile in infiles:
            r = make_histos(infile, det)

            txt = infile[infile.rfind('gun_') + 4:].replace('.root', '')
            cname = 'pfsize_{0}_{1}'.format(txt, det)
            combine(r, cname, '{0}, {1}'.format(txt, det))

def make_histos(infile, det):
    """Fills histograms with PFCluster sizes vs pT.

    Results are cached into file.
    """
    # return cached results, if any
    fname = os.path.basename(infile).replace('.root', '')
    cachefile = 'output/cache/draw_pfsize_{0}_{1}.pkl'.format(fname, det)
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            return pickle.load(f)

    # get TTree with PFClusters
    fi = ROOT.TFile(infile)
    tree = fi.Get('ntuplizer/PFClusterTree')
    if not tree:
        raise Exception('TTree not found')

    result = [ROOT.TH1D('h', '', 60, 0, 6) for _ in range(6)]

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

        # remove "fakes"; NOTE: cuts were evaluated with draw_inputs.py
        if tree.pfPhoDeltaR > 0.03 or tree.pfE/tree.mcE < 0.4:
            continue

        pfSize = min(len(result), tree.pfSize5x5_ZS)
        ptTrue = tree.mcPt

        result[pfSize - 1].Fill(ptTrue)

    # save cache
    with open(cachefile, 'wb') as f:
        pickle.dump(result, f)

    return result

def combine(histos, cname, title=''):
    """Visualization on single canvas.
    """
    c = ROOT.TCanvas(cname, cname, 700, 700)
    saves.append(c)

    c.SetLeftMargin(0.10)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.06)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()

    # draw dummy histogram
    xmax = max(h.GetXaxis().GetXmax() for h in histos)
    frame = c.DrawFrame(0, 0, xmax, 100)

    frame.SetTitle(title)
    frame.SetXTitle('p_{T}^{gen}')
    frame.SetYTitle('%')
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.4, 'Y')
    frame.Draw()

    # re-sum histograms
    hsums = [histos[0].Clone()]
    for h in histos[1:]:
        h2 = h.Clone()
        h2.Add(hsums[-1])
        hsums.append(h2)

    # scale histograms to 100% in Y axis
    last = hsums[-1]
    for h in hsums:
        for b in range(1, h.GetNbinsX() + 1):
            y1 = h.GetBinContent(b)
            y2 = last.GetBinContent(b)
            if y2 > 0:
                h.SetBinContent(b, y1/y2 * 100)

    hsums.reverse()

    clrs = [49, 40, ROOT.kBlue, 6, ROOT.kRed, 8, ROOT.kGreen, ROOT.kOrange]

    for (h, clr) in zip(hsums, clrs):
        h.SetLineColor(clr)
        h.SetFillColor(clr)
        h.SetFillStyle(3001)
        h.Draw('same')

    # legend
    leg = ROOT.TLegend(0.70, 0.89 - 0.033 * len(hsums), 0.91, 0.89)

    for (i, h) in enumerate(reversed(hsums)):
        sign = '#geq' if i == len(hsums) - 1 else '='
        leg.AddEntry(h, 'pfSize {0} {1}'.format(sign, i + 1), 'f')

    leg.SetFillColor(ROOT.kWhite)
    leg.Draw('same')

    saves.append((frame, hsums, leg))

    c.Update()
    c.SaveAs('output/plots/{0}.png'.format(c.GetTitle()))


if __name__ == '__main__':
    main()
