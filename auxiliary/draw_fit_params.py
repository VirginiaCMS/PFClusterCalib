#!/usr/bin/env python
"""Visualizes evolution of fitting function parameters with pT.
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
    for d in ['output', 'output/cache', 'output/plots_fit_params']:
        if not os.access(d, os.X_OK):
            os.mkdir(d)

    for det in ['EB', 'EE']:
        r = [make_graphs(f, det, blockSize=20000, pfSize=0) for f in infiles]

        # repack graphs into per-parameter tuples
        r = list(zip(*r))

        combine(r[0], txts, 'fit_amp_'     + det, det, 'amp')
        combine(r[1], txts, 'fit_mean_'    + det, det, 'mean')
        combine(r[2], txts, 'fit_sigma_'   + det, det, '#sigma')
        combine(r[3], txts, 'fit_alphaL_'  + det, det, '#alpha_{L}')
        combine(r[4], txts, 'fit_alphaR_'  + det, det, '#alpha_{R}')
        combine(r[5], txts, 'fit_powerR_'  + det, det, 'n_{R}')
        combine(r[6], txts, 'fit_chi2ndf_' + det, det, '#chi^{2}/ndf')

        # draw distributions of alphaR and powerR
        combine_distr(r[4], txts, 'fit_distr_alphaR_' + det, det, '#alpha_{R}', 11)
        combine_distr(r[5], txts, 'fit_distr_powerR_' + det, det, 'n_{R}', 110)

def make_graphs(infile, det, blockSize=10000, pfSize=0):
    """Fills, fits and visualizes distributions of Etrue/Erec.

    Results are cached into file.

    pfSize > 0: take PFClusters of only this size;
    pfSize = 0: take all PFClusters;
    pfSize < 0: take PFClusters of size |pfSize| or bigger.
    """
    # return cached results, if any
    fname = os.path.basename(infile).replace('.root', '')
    fmt = 'output/cache/draw_fit_params_{0}_{1}_{2}_{3}{4}.pkl'
    sign = 'p' if pfSize >= 0 else 'm'
    cachefile = fmt.format(fname, det, blockSize, sign, abs(pfSize))
    if os.access(cachefile, os.R_OK):
        with open(cachefile, 'rb') as f:
            return pickle.load(f)

    # get TTree with PFClusters
    fi = ROOT.TFile(infile)
    tree = fi.Get('ntuplizer/PFClusterTree')
    if not tree:
        raise Exception('TTree not found')

    data = []  # (pTtrue, Etrue/Erec) tuples

    # collect data
    for ev in range(tree.GetEntriesFast()):
        if tree.GetEntry(ev) <= 0:
            raise Exception

        # skip PFClusters of wrong size
        if pfSize > 0:
            if tree.pfSize5x5_ZS != pfSize:
                continue
        elif pfSize < 0:
            if tree.pfSize5x5_ZS < -pfSize:
                continue

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

        data.append((tree.mcPt, tree.mcE/tree.pfE))

    # sort by pT in increasing order
    data.sort(key=lambda x: x[0])

    grAmp    = ROOT.TGraphErrors()
    grMean   = ROOT.TGraphErrors()
    grSigma  = ROOT.TGraphErrors()
    grAlphaL = ROOT.TGraphErrors()
    grAlphaR = ROOT.TGraphErrors()
    grPowerR = ROOT.TGraphErrors()
    grChiNdf = ROOT.TGraphErrors()

    result = (grAmp, grMean, grSigma, grAlphaL, grAlphaR, grPowerR, grChiNdf)

    # loop over blocks of pT-ordered data
    # NOTE: last block is excluded if it has less than 0.5 * blockSize entries
    for b in range(int(round(len(data)/blockSize))):
        block = data[b * blockSize:(b + 1) * blockSize]

        # mean and sigma in the block of pT and Etrue/Erec
        (meanPt, sigmaPt)   = mean_sigma([x[0] for x in block])
        (meanTgt, sigmaTgt) = mean_sigma([x[1] for x in block])

        xmin = meanTgt - 5 * sigmaTgt
        xmax = meanTgt + 6 * sigmaTgt
        h = ROOT.TH1D('h', '', 100, xmin, xmax)

        for (pt, tgt) in block:
            h.Fill(tgt)

        # create new canvas, if necessary
        if b % 9 == 0:
            try:
                c.SaveAs('output/plots_fit_params/{0}.png'.format(c.GetTitle()))
            except UnboundLocalError:
                pass

            fmt = 'fits_{0}_{1}_blk{2:03d}to{3:03d}'
            cname = fmt.format(fname, det, b + 1, b + 9)
            c = ROOT.TCanvas(cname, cname, 1000, 700)
            saves.append(c)

            c.SetLeftMargin(0)
            c.SetRightMargin(0)
            c.SetTopMargin(0)
            c.SetBottomMargin(0)
            c.Divide(3, 3)

        c.cd(b % 9 + 1)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.02)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.08)

        h.SetTitle('p_{{T}}^{{gen}} = {0:.2f} #pm {1:.2f}'.format(meanPt, sigmaPt))
        h.SetXTitle('E^{gen}/E^{PF}')
        h.SetYTitle('Entries')
        h.SetTitleOffset(1.6, 'Y')

        h.Sumw2(True)
        h.SetLineColor(ROOT.kBlack)
        h.Draw()
        saves.append(h)

        # Gaussian + exponential left tail + power-law right tail
        gaus = 'exp(-(x-[1])^2/(2*[2]*[2]))'
        left = 'exp(0.5*[3]*[3] + [3]*(x-[1])/[2])'
        rght = '([5]/[4])^[5] * exp(-0.5*[4]^2) * ((x-[1])/[2]-[4]+[5]/[4])^(-[5])'

        expr = ('[0] * (' +
                    '(x-[1])/[2] > -[3] ? ' +
                        '( (x-[1])/[2] < [4] ? {0}:{1} ) :'.format(gaus, rght) +
                    '{0} )'.format(left))

        fit = ROOT.TF1('fit', expr, xmin, xmax)
        fit.SetLineWidth(1)
        fit.SetNpx(2000)

        fit.SetParameters(h.GetMaximum(), meanTgt, sigmaTgt, 1.5, 1.5, 20)
        fit.SetParLimits(0, 0.33 * h.GetMaximum(), 2 * h.GetMaximum())
        fit.SetParLimits(1, 0.9, 1.6)
        fit.SetParLimits(2, 0.33 * sigmaTgt, 3 * sigmaTgt)
        fit.SetParLimits(3, 0, 10)
        fit.SetParLimits(4, 0, 10)
        fit.SetParLimits(5, 1.01, 100)

        h.Fit(fit, 'QEML', 'same', xmin, xmax)

        ROOT.gPad.Update()

        for (i, gr) in enumerate(result[:-1]):
            gr.SetPoint(b, meanPt, fit.GetParameter(i))
            gr.SetPointError(b, sigmaPt, fit.GetParError(i))

        grChiNdf.SetPoint(b, meanPt, fit.GetChisquare()/fit.GetNDF())
        grChiNdf.SetPointError(b, sigmaPt, 0)

    # save the very last canvas
    c.SaveAs('output/plots_fit_params/{0}.png'.format(c.GetTitle()))

    # save cache
    with open(cachefile, 'wb') as f:
        pickle.dump(result, f)

    return result

def mean_sigma(numbers, nsigmas=3):
    """Returns average (mean) and dispersion (sigma) for an array of numbers.

    Mean and sigma are recalculated iteratively several times. During each
    calculation, a region [mean - nsigmas * sigma, mean + nsigmas * sigma] is
    used, where 'mean' and 'sigma' are taken from a previous iteration.
    """
    # zero-order iteration
    mean = sum(numbers)/len(numbers)
    sigma = (sum(x**2 for x in numbers)/len(numbers) - mean**2)**0.5

    # iterations
    for _ in range(1000):
        mean_prev = mean
        sigma_prev = sigma

        xmin = mean - nsigmas * sigma
        xmax = mean + nsigmas * sigma

        arr = [x for x in numbers if xmin <= x <= xmax]

        mean = sum(arr)/len(arr)
        sigma = (sum(x**2 for x in arr)/len(arr) - mean**2)**0.5

        # break when converged
        if (abs(mean - mean_prev) <= 1e-6 * abs(mean) and
            abs(sigma - sigma_prev) <= 1e-6 * abs(sigma)):
           break
    else:
        raise Exception('mean/sigma did not converged')

    return (mean, sigma)

def combine(grs, txts, cname, title, ytitle):
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
    xmax = max(gr.GetX()[i] for gr in grs for i in range(gr.GetN()))
    ymin = min(gr.GetY()[i] for gr in grs for i in range(gr.GetN()))
    ymax = max(gr.GetY()[i] for gr in grs for i in range(gr.GetN()))
    frame = c.DrawFrame(0, ymin * 0.9, xmax, ymax * 1.1)

    frame.SetTitle(title)
    frame.SetXTitle('p_{T}^{gen}')
    frame.SetYTitle(ytitle)
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.95, 'Y')
    frame.Draw()

    # legend
    leg = ROOT.TLegend(0.18, 0.18, 0.39, 0.18 + 0.033 * len(txts))

    clrs = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kOrange, ROOT.kGreen, 8] * 10

    for (gr, clr, txt) in zip(grs, clrs, txts):
        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(0.8)
        gr.SetMarkerColor(clr)
        gr.SetLineColor(clr)

        gr.Draw('PZL')
        leg.AddEntry(gr, txt, 'p')

    leg.SetFillColor(0)
    leg.Draw('same')

    saves.append((frame, grs, leg))

    c.Update()
    c.SaveAs('output/plots_fit_params/{0}.png'.format(c.GetTitle()))

def combine_distr(grs, txts, cname, title, xtitle, xmax):
    """Visualization of several distributions on single canvas.
    """
    histos = [ROOT.TH1D('h', '', 100, 0, xmax) for _ in range(len(grs))]

    # fill histograms
    for (h, gr) in zip(histos, grs):
        for i in range(gr.GetN()):
            h.Fill(gr.GetY()[i])

    c = ROOT.TCanvas(cname, cname, 700, 700)
    saves.append(c)

    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.08)
    c.SetTopMargin(0.06)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()

    # draw dummy histogram
    ymax = max(h.GetMaximum() for h in histos)
    frame = c.DrawFrame(0, 0, xmax, ymax * 1.1)

    frame.SetTitle(title)
    frame.SetXTitle(xtitle)
    frame.SetYTitle('Entries')
    frame.SetTitleOffset(1.2, 'X')
    frame.SetTitleOffset(1.95, 'Y')
    frame.Draw()

    # legend
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
    c.SaveAs('output/plots_fit_params/{0}.png'.format(c.GetTitle()))


if __name__ == '__main__':
    main()
