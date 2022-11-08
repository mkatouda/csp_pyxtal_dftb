import os
import shutil
from time import time
import pandas as pd
from pyxtal import pyxtal
from .crys_geom_opt import get_ase_calculator, crys_geom_opt
from ase.io import read, write


def opt_molcrys(basename, mols, nmols, spg, diag=False, nstruc=100, 
                factor=1.0, t_factor=1.0, 
                qc_method='DFTB', opt_method='LBFGS', opt_fmax=5e-1, 
                opt_maxsteps=50, opt_maxstepsize=0.01, 
                symprec=0.1, istruc_bgn=1, istruc_end=101):

    if istruc_end > nstruc + 1:
        istruc_end = nstruc + 1        

    initdir = 'init'

    xtb_paramfile = 'param_gfn0-xtb.txt'
    calc = get_ase_calculator(qc_method)
    column_titles = ['ID', 'energy', 'energy_per_mol', 'density', 'volume', 
                     'a', 'b', 'c', 'alpha', 'beta', 'gamma', 
                     'qc_method', 'opt_method', 
                     'opt_fmax', 'opt_maxstepsize', 'opt_converged', 'opt_steps', 'time']
    column_titles_summary = ['ID', 'energy', 'energy_per_mol', 'density', 'volume', 
                             'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'opt_converged']

    lendfmax = 20
    pd.set_option('display.width', 150)
    pd.set_option('display.max_columns', None)
    pd.options.display.float_format = '{:.2f}'.format

    t1 = time()
    cwdir = os.getcwd()
    print('Optimize molecular crystal in spacegroup {}.\n'.format(spg))

    sdiag = ''
    if diag:
        sdiag = '-diag'
    spgdir = '{}_spg{}{}'.format(basename, spg, sdiag)
    #if not os.path.exists(spgdir):
    #    os.makedirs(spgdir)
    os.chdir(spgdir)
    gendir = 'factor{:.2f}_tfactor{:.2f}'.format(factor, t_factor)
    #if not os.path.exists(gendir):
    #    os.makedirs(gendir)
    os.chdir(gendir)

    sdiag = ''
    if diag:
        sdiag = '-diag'
    basename1 = '{}_spg{}{}_factor{:.2f}_tfactor{:.2f}'.format(basename, spg, sdiag, factor, t_factor)
    qcdir = '{}'.format(qc_method)
    if not os.path.exists(qcdir):
        os.makedirs(qcdir)
    os.chdir(qcdir)
    optdir = '{}'.format(opt_method)
    if not os.path.exists(optdir):
        os.makedirs(optdir)
    os.chdir(optdir)

    optstepdir = 'fmax{:.4f}_maxsteps{:.5f}_optcycles{:04}'.format(opt_fmax, opt_maxstepsize, opt_maxsteps)
    if not os.path.exists(optstepdir):
        os.makedirs(optstepdir)
    os.chdir(optstepdir)

    if qc_method == 'xTB':
        shutil.copyfile(cwdir + '/' + xtb_paramfile, xtb_paramfile)

    basename2 = '{}_{}_{}_fmax{:.4f}_maxsteps{:.5f}_optcycles{:04}'.format(basename1, qc_method, opt_method, 
                                                                               opt_fmax, opt_maxstepsize, opt_maxsteps)

    df = pd.DataFrame(index=[], columns=column_titles)
    for istruc in range(istruc_bgn, istruc_end):

        basename3 = '{}_{:06}'.format(basename2, istruc)
        print(basename3)
        logfile = '{}.log'.format(basename3)
        trajfile = '{}.traj'.format(basename3)
        infile = '{}_{:06}_init.cif'.format(basename1, istruc)
        wrkdir = basename3
        if not os.path.exists(wrkdir):
            os.makedirs(wrkdir)
        os.chdir(wrkdir)
        #pyxtal_struc = pyxtal()
        #pyxtal_struc.from_seed('{}/{}/{}/{}/{}'.format(cwdir, spgdir, gendir, initdir, infile))
        #ase_struc = pyxtal_struc.to_ase(resort=False)
        ase_struc = read('{}/{}/{}/{}/{}'.format(cwdir, spgdir, gendir, initdir, infile))
        ase_struc.calc = calc
        tg1 = time()
        ase_struc, energy, opt_steps, opt_converged = crys_geom_opt(ase_struc, basename, 
                                                                    opt_method, opt_fmax, opt_maxsteps,
                                                                    logfile, trajfile, opt_maxstepsize, symprec)

        tg = time() - tg1
        print('Relaxation took {:3f} seconds.\n'.format(tg))
        os.chdir('../')
        shutil.copyfile(wrkdir + '/' + logfile, logfile)
        shutil.copyfile(wrkdir + '/' + trajfile, trajfile)
        shutil.rmtree(wrkdir)

        energy_per_mol = energy / nmols[0]
        pyxtal_struc = pyxtal()
        pyxtal_struc.from_seed(ase_struc)
        density = pyxtal_struc.get_density()

        print('Final lattice:')
        print(pyxtal_struc)
        print('Energy:', energy)
        print('Density:', density)
              
        volume = pyxtal_struc.lattice.volume
        [a, b, c, alpha, beta, gamma] = pyxtal_struc.lattice.get_para(degree=True)
        record = pd.Series([istruc, energy, energy_per_mol, density, volume, a, b, c, 
                            alpha, beta, gamma, qc_method, opt_method,
                            opt_fmax, opt_maxstepsize, opt_converged, opt_steps, tg], 
                           index=df.columns)
                                   
        df = df.append(record, ignore_index=True)
        #print(df)
        outfile = '{}_opt.cif'.format(basename3)
        #pyxtal_struc.to_file(outfile)
        write(outfile, ase_struc, format='cif')

    outfile = '{}_nstruc{:06}-{:06}.csv'.format(basename2, istruc_bgn, istruc_end - 1)
    df['ID'] = df['ID'].astype('int')
    df.to_csv(outfile, float_format='%.3f')

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

    t2 = time()
    print('All of crystal generation took {:3f} seconds.\n'.format(t2 - t1))
