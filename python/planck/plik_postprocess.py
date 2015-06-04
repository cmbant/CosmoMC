#! /usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import numpy as nm
import clik
import clik.smicahlp as smh
import re
import sys


def parse_minimumfile(minfile):
    dc = {}
    f = open(minfile)
    for l in f:
        rr = re.search("\s+(\d+)\s+(.+?)\s+(.+?)\s+", l)
        if rr is None:
            continue
        vl = float(rr.group(2))
        pr = rr.group(3)
        dc[pr] = vl
    return dc


def parse_paraname(parfile):
    f = open(parfile)
    return [l.strip().split()[0] for l in f if l.strip()]


cache = dict()


def getLk(plikfile):
    lkl = cache.get(plikfile, None)
    if not lkl:
        lkl = clik.clik(plikfile)
        cache[plikfile] = lkl
    return lkl


def main_fgfile(argv):
    plikfile = argv[1]
    parfile = argv[2]
    minfile = argv[3]
    lmax = 2999
    where = argv[4]

    assert lmax < 3000, "lmax too high (2999 max)"
    # lkl = clik.clik(plikfile)
    lkl = getLk(plikfile)
    keys = lkl.get_extra_parameter_names()
    lmaxs = lkl.lmax
    # del(lkl)
    prms = smh.parametric_from_smica(plikfile, lmin=0, lmax=lmax)
    dc = parse_minimumfile(minfile)
    ppf = parse_paraname(parfile)
    rdc = dict([(keys[i], dc[ppf[i]]) for i in range(len(ppf))])

    res = nm.zeros((lmax + 1, 3, 3))
    for p in prms:
        ppars = [rdc[k] for k in p.varpar]
        pes = p(ppars)
        if pes.shape[1] > res.shape[1]:
            pes[:, :res.shape[1], :res.shape[1]] += res
            res = pes
        else:
            res = res + pes

    if res.shape[1] != 6:

        ges = nm.zeros((lmax + 1, 6, 6))

        if lmaxs[0] == -1:
            ges[:, 3:, 3:] = res
        else:
            ges[:, :3, :3] = res
        res = ges

    llp1 = (nm.arange(lmax + 1) * (nm.arange(lmax + 1) + 1.) / 2. / nm.pi)
    Dl = res * (llp1[:, nm.newaxis, nm.newaxis])
    print("save fg file to :", where)

    f = open(where, "w")
    print("#    l      TT100x100      TT143X143      TT143X217      TT217X217      EE100X100      EE100X143      EE100X217      EE143X143      EE143X217      EE217X217      TE100X100      TE100X143      TE100X217      TE143X143      TE143X217      TE217X217", file=f)
    print("#", file=f)
    for i in range(2, lmax + 1):
        print(("%6d" + " %14E" * 16) % (
            i, Dl[i, 0, 0], Dl[i, 1, 1], Dl[i, 1, 2], Dl[i, 2, 2], Dl[i, 3, 3], Dl[i, 3, 4], Dl[i, 3, 5], Dl[i, 4, 4],
            Dl[i, 4, 5], Dl[i, 5, 5], Dl[i, 0, 3], Dl[i, 0, 4], Dl[i, 0, 5], Dl[i, 1, 4], Dl[i, 1, 5], Dl[i, 2, 5]), file=f)

    f.close()


def main_fg2000(argv, tag=''):
    plikfile = argv[1]
    parfile = argv[2]
    minfile = argv[3]
    chainparfile = argv[4]
    chains = argv[5:]
    # lkl = clik.clik(plikfile)
    lkl = getLk(plikfile)

    keys = lkl.get_extra_parameter_names()
    lmaxs = lkl.lmax
    # del(lkl)
    assert lmaxs[0] != -1, "can only add T fg in T cases !"
    llp1 = 2000 * 2001 / 2. / nm.pi

    prms = smh.parametric_from_smica(plikfile, lmin=2000, lmax=2000)
    ppf = parse_paraname(parfile)
    cpf = parse_paraname(chainparfile)
    ldc = []
    for i in range(len(ppf)):
        if ppf[i] in cpf:
            ldc += [(keys[i], cpf.index(ppf[i]))]

    ldc = dict(ldc)

    nchainparfile = chainparfile[:-len(".paramnames")] + tag + ".paramnames"
    ls = open(chainparfile, "r").readlines()
    print("save paranames file to :", nchainparfile)
    f = open(nchainparfile, "w")
    for i, l in enumerate(ls):
        if re.search("^chi2_", l):
            break
        print(l.rstrip(), file=f)
    ib = i
    print("f2000_143*\tf_{2000}^{143}", file=f)
    print("f2000_x*\tf_{2000}^{143\\times217}", file=f)
    print("f2000_217*\tf_{2000}^{217}", file=f)
    for i, l in enumerate(ls[ib:]):
        print(l.rstrip(), file=f)

    f.close()

    if minfile.endswith("ranges"):
        dc = {}
        for l in open(minfile):
            a, b, c = l.split()
            if b == c:
                dc[a] = float(b)
        rdc = dict([(keys[i], dc.get(ppf[i], None)) for i in range(len(ppf))])
        print(rdc)
    else:
        dc = parse_minimumfile(minfile)
        rdc = dict([(keys[i], dc[ppf[i]]) for i in range(len(ppf))])
        nminimumfile = minfile[:-len(".minimum")] + tag + ".minimum"
        print("save minimum file to :", nminimumfile)
        fo = open(minfile, "r")
        ls = fo.readlines()
        fo.close()
        f = open(nminimumfile, "w")
        fnd = False
        for i, l in enumerate(ls):
            if re.search(r"\s+\d+\s+.+?\s+(chi2_.+?)\s+", l):
                fnd = True
                break
        assert fnd, "bad format for %s" % minfile
        ibreak = i
        for l in ls[:ibreak]:
            print(l.rstrip(), file=f)

        vl = int(ls[ibreak].strip().split()[0])

        res = 0
        for p in prms:
            ppars = [rdc[k] for k in p.varpar]
            res = res + p(ppars)[:, :3, :3]

        print("%5d %#14E   %-21s %s" % (vl, res[0, 1, 1] * llp1, "f2000_143", "f_{2000}^{143}"), file=f)
        print("%5d %#14E   %-21s %s" % (vl + 1, res[0, 1, 2] * llp1, "f2000_x", "f_{2000}^{143\\times217}"), file=f)
        print("%5d %#14E   %-21s %s" % (vl + 2, res[0, 2, 2] * llp1, "f2000_217", "f_{2000}^{217}"), file=f)

        cnt = vl + 3
        keep = True
        for l in ls[ibreak:]:
            rr = re.search("\s+(\d+)\s+(.+?)\s+(.+?)\s+", l)
            if rr and keep:
                print("%5d%s" % (cnt, l[5:].rstrip()), file=f)
                cnt += 1
            else:
                keep = False
                print(l.rstrip(), file=f)
        f.close()

    for chain in chains:
        ch = nm.loadtxt(chain)
        nmm = chain.split("_")
        if tag:
            nchain = "_".join(nmm[:-1] + [tag.replace('_', ''), nmm[-1]])
        else:
            nchain = chain
        f = open(nchain, "w")
        print(nchain, ":", end=' ')
        for i in range((ch.shape[0])):
            if i % 1000 == 0:
                print("%d/%d" % (i, ch.shape[0]), end=' ')
                sys.stdout.flush()
            res = 0
            for p in prms:
                ppars = []
                for k in p.varpar:
                    if k in ldc:
                        ppars += [ch[i, ldc[k] + 2]]
                    else:
                        ppars += [rdc[k]]
                res = res + p(ppars)[:, :3, :3]
            print("", end=' ', file=f)
            for v in ch[i][:ib + 2]:
                print(" %14E" % v, end=' ', file=f)
            print(" %14E  %14E  %14E" % (res[0, 1, 1] * llp1, res[0, 1, 2] * llp1, res[0, 2, 2] * llp1), end=' ', file=f)
            for v in ch[i][ib + 2:]:
                print(" %14E" % v, end=' ', file=f)
            print("", file=f)
        f.close()
        print("")


def main(argv):
    if sys.argv[1] == "fgfile":
        main_fgfile(argv[1:])
    if sys.argv[1] == "fg2000":
        main_fg2000(argv[1:])


if __name__ == "__main__":
    main(sys.argv)
