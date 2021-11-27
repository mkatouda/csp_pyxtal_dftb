import os
from time import time
#import pandas as pd
from pyxtal import pyxtal
from ase.io import write


def gen_random_molcrys(basename, mols, nmols, spg, diag=False, nstruc=100, 
                       factor=1.0, t_factor=1.0, nstruc_try=500, istruc_bgn=1, istruc_end=100):

    if istruc_end > nstruc + 1:
        istruc_end = nstruc + 1        

    t1 = time()
    cwdir = os.getcwd()
    print('Genetate molecular crystal in spacegroup {}.\n'.format(spg))

    sdiag = ''
    if diag:
        sdiag = '-diag'
    spgdir = '{}_spg{}{}'.format(basename, spg, sdiag)
    if not os.path.exists(spgdir):
        os.makedirs(spgdir)
    os.chdir(spgdir)
    gendir = 'factor{:.2f}_tfactor{:.2f}'.format(factor, t_factor)
    if not os.path.exists(gendir):
        os.makedirs(gendir)
    os.chdir(gendir)
    initdir = 'init'
    if not os.path.exists(initdir):
        os.makedirs(initdir)
    os.chdir(initdir)

    basename1 = '{}_spg{}{}_factor{:.2f}_tfactor{:.2f}'.format(basename, spg, sdiag, factor, t_factor)
    istruc = 0
    for istruc_try in range(1, nstruc_try+1):

        print('Trial {}'.format(istruc_try))

        pyxtal_struc = pyxtal(molecular=True)
        pyxtal_struc.from_random(3, spg, mols, nmols, factor=factor, 
                                 t_factor=t_factor, diag=diag)
        if pyxtal_struc.valid:
            istruc += 1
            print('Trial {} is valid, and add structure {}\n'.format(istruc_try, istruc))

            energy = None; energy_per_mol = None
            density = pyxtal_struc.get_density()
            opt_converged = None; opt_steps = None; tg = None

            print('Initial lattice:')
            print(pyxtal_struc)
            print('Density:', density)

            outfile = '{}_{:06}_init.cif'.format(basename1, istruc)
            #pyxtal_struc.to_file(outfile)
            ase_struc = pyxtal_struc.to_ase(resort=False)
            write(outfile, ase_struc, format='cif')

            if istruc == nstruc:
                break

    os.chdir(cwdir)
    t2 = time()
    print('All of crystal generation took {:3f} seconds.\n'.format(t2 - t1))


def get_parser():
    parser = argparse.ArgumentParser(
        description='PyXtal DFTB crystal structure prediction code',
        usage=f"python {os.path.basename(__file__)} -c CONFIG_FILE"
    )
    parser.add_argument(
        "-c", "--config", type=str, #required=True,
        help="path to a config file"
    )
    parser.add_argument(
        "-d", "--debug", action='store_true',
        help="debug mode"
    )
    parser.add_argument('-nogen', action='store_false')
    parser.add_argument('-noopt', action='store_false')
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
    parser.add_argument('-istbgn', type=int, default=1)
    parser.add_argument('-istend', type=int, default=100)

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

    gen_random_molcrys(basename, mols, nmols, spg, diag, nstruc, 
                       factor, t_factor, nstruc_try, istruc_bgn, istruc_end)


if __name__ == "__main__":
    main()
