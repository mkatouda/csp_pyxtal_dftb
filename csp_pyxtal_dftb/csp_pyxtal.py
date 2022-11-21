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
        description='PyXtal GAFF & DFTB crystal structure prediction code'
    )
    parser = argparse.ArgumentParser(
        formatter_class=customHelpFormatter,
    )
    parser.add_argument(
        '-i', '--inp', type=str,
        help = 'yaml style input file, overwriting argument values',
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help = 'Verbose output.'
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

def set_config(args):
    # Read config yaml file
    if args.inp is not None and os.path.isfile(args.inp):
        with open(args.inp, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    return conf

def set_default_config(conf):
    conf.setdefault('basename', 'benzene')
    conf.setdefault('mols', ['c1ccccc1.smi'])
    conf.setdefault('nmols', [4])
    conf.setdefault('spg',  14) # 'P21/c'

    conf.setdefault('nstruc', 100)
    conf.setdefault('istbgn', 1)
    conf.setdefault('istend', 100)
    conf.setdefault('nst_sub', 100)
    conf.setdefault('factor', 1.3)
    conf.setdefault('tfactor', 1.5)
    conf.setdefault('use_hall', False)
    #conf.setdefault('nstruc_try', 500)

    conf.setdefault('nxyz', [1, 1, 1])
    conf.setdefault('ff', 'gaff2')
    conf.setdefault('charge_model', 'bcc')

    conf.setdefault('sim_method', 'xtb')
    conf.setdefault('kpts', [1, 1, 1])
    conf.setdefault('ps_path', None)
    conf.setdefault('qe_input', {'system': {'vdw_corr': 'DFT-D3', 'dftd3_version': 3}})
    conf.setdefault('cp2k_input', '')
    conf.setdefault('xtb_hamiltonian', 'GFN2-xTB')

    conf.setdefault('opt_method', 'LBFGS')
    conf.setdefault('opt_fmax', 0.01)
    conf.setdefault('opt_maxsteps', 50)
    conf.setdefault('opt_maxstepsize', 0.01)
    conf.setdefault('symprec', 0.1)

    conf.setdefault('gen_crys', True)
    conf.setdefault('geom_opt', True)
    conf.setdefault('free_energy', False)

    conf.setdefault('verbose', True)

    return conf

def csp_pyxtal_main(conf):

    basename = conf['basename']
    mols = conf['mols']
    nmols = conf['nmols']
    spg = conf['spg']

    gen_crys = conf['nogen']
    nstruc = conf['nstruc']
    factor =  conf['factor']
    t_factor = conf['tfactor']
    use_hall = conf['use_hall']
    #nstruc_try = conf['nstruc_try']

    nxyz = conf['nxyz']
    ff = conf['ff']
    charge_model = conf['charge_model']

    geom_opt = conf['geom_opt']
    sim_method = conf['sim_method']
    kpts = conf['kpts']
    ps_path = conf['ps_path']
    qe_input = conf['qe_input']
    cp2k_input = conf['cp2k_input']
    xtb_hamiltonian = conf['xtb_hamiltonian']

    opt_method = conf['opt_method']
    opt_fmax = conf['opt_fmax']
    opt_maxsteps = conf['opt_maxsteps']
    opt_maxstepsize = conf['opt_maxstepsize']
    symprec = conf['symprec']

    free_energy = conf['free_energy']

    nstruc_try = nstruc * 5
    istruc_bgn = conf['istbgn']
    istruc_end = conf['istend'] + 1
    nstruc_sub = conf['nstsub']
    nstbat = int(nstruc / nstruc_sub)

    verbose = conf['verbose']

    if gen_crys:
        gen_random_molcrys(basename, mols, nmols, spg,
                           nstruc, factor, t_factor, use_hall,
                           nstruc_try, istruc_bgn, istruc_end,
                           nxyz, ff, charge_model,
                           verbose)

    if geom_opt:
        opt_molcrys(basename, mols, nmols, spg, nstruc,
                    factor, t_factor, use_hall,
                    istruc_bgn, istruc_end,
                    sim_method, kpts, ps_path,
                    qe_input, cp2k_input, xtb_hamiltonian,
                    opt_method, opt_fmax, 
                    opt_maxsteps, opt_maxstepsize,
                    symprec, free_energy, verbose)


def main():
    args = get_parser()
    print(args)

    conf = set_config(args)
    conf = set_default_config(conf)

    print('======= Input configulations =======')
    for k, v in conf.items():
        print('{}: {}'.format(k, v))
    print('====================================')

    csp_pyxtal_main(conf)

if __name__ == '__main__':
    main()
