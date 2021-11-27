import os
import shutil
import pandas as pd
import argparse
import yaml


def merge_csv(basename, spg, diag, nstruc, factor, t_factor, qc_method, 
              opt_method, opt_fmax, opt_maxcycle, opt_maxstepsize, nstruc_sub):

    nstbat = int(nstruc / nstruc_sub)

    column_titles_summary = ['ID', 'energy', 'energy_per_mol', 'density', 'volume', 
                             'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'opt_converged']

    lendfmax = 20
    pd.set_option('display.width', 150)
    pd.set_option('display.max_columns', None)
    pd.options.display.float_format = '{:.2f}'.format

    cwdir = os.getcwd()
    csvdir = '{}_csv'.format(basename)
    if not os.path.exists(csvdir):
        os.makedirs(csvdir)
    print(csvdir)

    sdiag = ''
    if diag:
        sdiag = '-diag'
    spgdir = '{}_spg{}{}'.format(basename, spg, sdiag)
    gendir = 'factor{:.2f}_tfactor{:.2f}'.format(factor, t_factor)
    basename1 = '{}_spg{}_factor{:.2f}_tfactor{:.2f}'.format(basename, spg, factor, t_factor)
    qcdir = '{}'.format(qc_method)
    optdir = '{}'.format(opt_method)
    optstepdir = 'fmax{:.4f}_maxsteps{:.5f}_optcycles{:04}'.format(opt_fmax, opt_maxstepsize, opt_maxcycle)
    wrkdir = '{}/{}/{}/{}/{}/{}'.format(cwdir, spgdir, gendir, qcdir, optdir, optstepdir)
    basename2 = '{}_{}_{}_fmax{:.4f}_maxsteps{:.5f}_optcycles{:04}'.format(basename1, qc_method, opt_method, opt_fmax, opt_maxstepsize, opt_maxcycle)
    print(basename2)
    #print(wrkdir)
    if os.path.exists(wrkdir):
        print(wrkdir)
        os.chdir(wrkdir)
        csvinfiles = []
        for istbat in range(nstbat):
            istruc_bgn = nstruc_sub * istbat + 1
            istruc_end = nstruc_sub * (istbat + 1)
            csvinfile = '{}_nstruc{:06}-{:06}.csv'.format(basename2, istruc_bgn, istruc_end)
            #print(csvinfile)

            if os.path.isfile(csvinfile):
                print(csvinfile)
                csvinfiles.append(pd.read_csv(csvinfile))

        df = pd.concat(csvinfiles, axis=0, ignore_index=True)
        df = df.drop('Unnamed: 0', axis=1)
        csvoutfile = '{}.csv'.format(basename2)
        print(csvoutfile)
        df.to_csv(csvoutfile)
                            
        shutil.copyfile(csvoutfile, cwdir + '/' + csvdir + '/' + csvoutfile)

        df_summary = df[column_titles_summary]
        lendf = min([len(df_summary), lendfmax])

        print('\n')
        print('Summary of generated structures: ', basename2)
        print(df_summary)
        print('\n')
                    
        print('Summary of generated structures (sorted in energies): ', basename2)
        print(df_summary.sort_values('energy', ascending=True).head(lendf))

        print('\n')

        print('Summary of generated structures (sorted in densities) ', basename2)
        print(df_summary.sort_values('density', ascending=False).head(lendf))
        print('\n')

        os.chdir(cwdir)


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
    #parser.add_argument('-nstruc', type=int, default=100)
    #parser.add_argument('-factor', type=float, default=1.1)
    #parser.add_argument('-tfactor', type=float, default=1.0)
    #parser.add_argument('-diag', action='store_true')
    #parser.add_argument('-qc', default='DFTB')
    #parser.add_argument('-opt', default='LBFGS')
    #parser.add_argument('-opt_fmax', type=float, default=0.5)
    #parser.add_argument('-optmaxcycle', type=int, default=50)
    #parser.add_argument('-optstepsize', type=float, default=0.01)
    #parser.add_argument('-symprec', type=float, default=0.1)
    parser.add_argument('-istbgn', type=int, default=1)
    parser.add_argument('-istend', type=int, default=100)
    parser.add_argument('-nstsub', type=int, default=100)

    return parser.parse_args()


def set_default_config(conf):

    conf.setdefault('basename', 'benzene')
    conf.setdefault('mols', ['c1ccccc1.smi'])
    conf.setdefault('nmols', [4])
    conf.setdefault('spg',  14) # 'P21/c'
    conf.setdefault('diag', False)

    conf.setdefault('nstruc', 100)
    conf.setdefault('factor', 1.0)
    conf.setdefault('tfactor', 1.0)
    conf.setdefault('nstruc_try', 500)

    conf.setdefault('symprec', 0.1)

    conf.setdefault('qc_method', 'DFTB')
    conf.setdefault('opt_method', 'LBFGS')
    conf.setdefault('opt_maxcycle', 50)
    conf.setdefault('opt_opt_fmax', 0.01)
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

    nstruc = conf['nstruc']
    factor =  conf['factor']
    t_factor = conf['tfactor']
    diag = conf['diag']

    nstruc_try = nstruc * 5

    symprec = conf['symprec']

    qc_method = conf['qc_method']
    opt_method = conf['opt_method']
    opt_maxcycle = conf['opt_maxcycle']
    opt_fmax = conf['opt_fmax']
    opt_maxstepsize = conf['opt_maxstepsize']

    istruc_bgn = 1
    istruc_end = nstruc

    nstruc_sub = args.nstsub

    merge_csv(basename, spg, diag, nstruc, factor, t_factor,
              qc_method, opt_method, opt_fmax, opt_maxcycle, opt_maxstepsize,
              nstruc_sub)


if __name__ == "__main__":
    main()
