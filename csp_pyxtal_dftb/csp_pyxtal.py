import os
import argparse
import yaml
from time import time
from .gen_molcrys import gen_random_molcrys
from .opt_molcrys import opt_molcrys


def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        formatter_class=customHelpFormatter,
        description='PyXtal GAFF & DFTB crystal structure prediction code'
    )
    parser.add_argument(
        '-i', '--inp', type=str,
        help = 'yaml style input file, overwriting argument values',
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
    conf.setdefault('nstsub', 100)
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

    conf.setdefault('strucdiff_method', 'POWDER')
    conf.setdefault('refstrucfile', None)

    if os.getenv('AMBERHOME') is None:
        dftb_prefix = '~/slko/3ob-3-1/'
    else:
        dftb_prefix = os.getenv('AMBERHOME') + '/dat/slko/3ob-3-1/'
    conf.setdefault('dftb_command', 'dftb+ > PREFIX.out')
    conf.setdefault('dftb_prefix', dftb_prefix)

    conf.setdefault('cp2k_command', 'cp2k_shell.ssmp')
    conf.setdefault('cp2k_data_dir', 'BASIS_MOLOPT')

    if os.getenv('CONDA_PREFIX') is None:
        espresso_pseudo = '~/espresso/pseudo/'
    else:
        espresso_pseudo = os.getenv('CONDA_PREFIX') + '/share/sssp/efficiency/'
    conf.setdefault('espresso_command', 'pw.x -in PREFIX.pwi > PREFIX.pwo')
    conf.setdefault('espresso_pseudo', espresso_pseudo)

    conf.setdefault('critic2_command', 'critic2')

    conf.setdefault('verbose', True)

    return conf

def csp_pyxtal_main(conf):

    basename = conf['basename']
    mols = conf['mols']
    nmols = conf['nmols']
    spg = conf['spg']

    gen_crys = conf['gen_crys']
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

    strucdiff_method = conf['strucdiff_method']
    refstrucfile = conf['refstrucfile']

    verbose = conf['verbose']

    if os.getenv('ASE_DFTB_COMMAND') is None:
        os.environ['ASE_DFTB_COMMAND'] = conf['dftb_command']
    if os.getenv('DFTB_PREFIX') is None:
        os.environ['DFTB_PREFIX'] = conf['dftb_prefix']

    if os.getenv('ASE_ESPRESSO_COMMAND') is None:
        os.environ['ASE_ESPRESSO_COMMAND'] = conf['espresso_command']
    if os.getenv('ESPRESSO_PSEUDO') is None:
        os.environ['ESPRESSO_PSEUDO'] = conf['espresso_pseudo']

    if os.getenv('ASE_CP2K_COMMAND') is None:
        os.environ['ASE_CP2K_COMMAND'] = conf['cp2k_command']
    if os.getenv('CP2K_DATA_DIR') is None:
        os.environ['CP2K_DATA_DIR'] = conf['cp2k_data_dir']

    if os.getenv('CRITIC2_COMMAND') is None:
        os.environ['CRITIC2_COMMAND'] = conf['critic2_command']

    if gen_crys:
        gen_random_molcrys(basename, mols, nmols, spg,
                           nstruc, factor, t_factor, use_hall,
                           nstruc_try, istruc_bgn, istruc_end,
                           nxyz, ff, charge_model,
                           strucdiff_method, verbose)

    if geom_opt:
        opt_molcrys(basename, mols, nmols, spg, nstruc,
                    factor, t_factor, use_hall,
                    istruc_bgn, istruc_end,
                    sim_method, kpts, ps_path,
                    qe_input, cp2k_input, xtb_hamiltonian,
                    opt_method, opt_fmax, 
                    opt_maxsteps, opt_maxstepsize,
                    symprec, strucdiff_method, refstrucfile,
                    free_energy, verbose)


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
