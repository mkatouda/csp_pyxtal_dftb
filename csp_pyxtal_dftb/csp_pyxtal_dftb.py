import os
import argparse
import yaml
from time import time
from .gen_random_molcrys import gen_random_molcrys
from .opt_molcrys import opt_molcrys


def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        formatter_class=customHelpFormatter,
        description='PyXtal DFTB crystal structure prediction code'
    )
    parser.add_argument(
        '-c', '--config', type=str, #required=True,
        help = "yaml style input file, overwriting argument values"
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help = "debug mode"
    )
    parser.add_argument(
        '-nogen', action='store_false'
    )
    parser.add_argument(
        '-noopt', action='store_false'
    )
    #parser.add_argument('-nstruc', type=int, default=100)
    #parser.add_argument('-factor', type=float, default=1.1)
    #parser.add_argument('-tfactor', type=float, default=1.0)
    #parser.add_argument('-diag', action='store_true')
    #parser.add_argument('-qc', default='DFTB')
    #parser.add_argument('-opt', default='LBFGS')
    #parser.add_argument('-fmax', type=float, default=0.5)
    #parser.add_argument('-optmaxcycle', type=int, default=50)
    #parser.add_argument('-optstepsize', type=float, default=0.01)
    #parser.add_argument('-symprec', type=float, default=0.1)
    parser.add_argument(
        '-istbgn', type=int, default=1
    )
    parser.add_argument(
        '-istend', type=int, default=100
    )
    parser.add_argument(
        '-nstsub', type=int, default=100
    )

    return parser.parse_args()


def set_default_config(conf):
    conf.setdefault('basename', 'benzene')
    conf.setdefault('mols', ['c1ccccc1.smi'])
    #conf.setdefault('mols', ['/home/usr4/q70204a/work/csp_kyoto/pyxtal_dftb_test/aspirin_laqaconf_test01/aspirin_laqa_gfn2-xtb_id0354.xyz'])
    conf.setdefault('nmols', [4])
    conf.setdefault('spg',  14) # 'P21/c'
    conf.setdefault('diag', False)

    conf.setdefault('nstruc', 100)
    conf.setdefault('factor', 1.0)
    conf.setdefault('tfactor', 1.0)
    #conf.setdefault('nstruc_try', 500)

    conf.setdefault('symprec', 0.1)

    conf.setdefault('qc_method', 'DFTB')
    conf.setdefault('opt_method', 'LBFGS')
    conf.setdefault('opt_maxcycle', 50)
    conf.setdefault('opt_fmax', 0.01)
    conf.setdefault('opt_maxstepsize', 0.01)

    #conf.setdefault('istbgn', 1)
    #conf.setdefault('istend', 100)


def main():
    args = get_parser()
    print(args)
    conf = {}
    with open(args.config, "r") as f:
        conf = yaml.load(f, Loader=yaml.SafeLoader)
    set_default_config(conf)
    print(conf)
    print(type(conf))

    basename = conf['basename']
    mols = conf['mols']
    nmols = conf['nmols']
    spg = conf['spg']
    diag = conf['diag']

    gen_crys = args.nogen
    nstruc = conf['nstruc']
    factor =  conf['factor']
    t_factor = conf['tfactor']
    #nstruc_try = conf['nstruc_try']

    symprec = conf['symprec']

    geom_opt = args.noopt
    qc_method = conf['qc_method']
    opt_method = conf['opt_method']
    opt_maxcycle = conf['opt_maxcycle']
    opt_fmax = conf['opt_fmax']
    opt_maxstepsize = conf['opt_maxstepsize']

    nstruc_try = nstruc * 5
    istruc_bgn = args.istbgn
    istruc_end = args.istend + 1
    nstruc_sub = args.nstsub
    nstbat = int(nstruc / nstruc_sub)

    if gen_crys:
        gen_random_molcrys(basename, mols, nmols, spg, diag, nstruc, 
                           factor, t_factor, nstruc_try, istruc_bgn, istruc_end)

    if geom_opt:
        opt_molcrys(basename, mols, nmols, spg, diag, nstruc, factor, t_factor, 
                    qc_method, opt_method, opt_fmax, opt_maxcycle, opt_maxstepsize,
                    symprec, istruc_bgn, istruc_end)


if __name__ == "__main__":
    main()
