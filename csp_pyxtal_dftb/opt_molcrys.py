import os
import shutil
from time import time

import pandas as pd
from ase.io import read, write
from pyxtal import pyxtal

from .molcrystop import crys_geom_opt, write_proteindatabank
from .phonon import free_energy_calc
from .critic2 import compare_crys


def opt_molcrys(basename, mols, nmols, spg, nstruc=100, 
                factor=1.0, t_factor=1.0, use_hall=False,
                istruc_bgn=1, istruc_end=101,
                sim_method='LAMMPS', kpts=(1,1,1), ps_path=None,
                qe_input={'system': {'vdw_corr': 'DFT-D3', 'dftd3_version': 3}}, cp2k_input='',
                xtb_hamiltonian='GFN2-xTB',
                opt_method='LBFGS', opt_fmax=5e-2, opt_maxsteps=100, 
                opt_maxstepsize=0.01, symprec=0.1, 
                strucdiff_method='POWDER', refstrucfile=None, free_energy=False, verbose=False):

    if istruc_end > nstruc + 1:
        istruc_end = nstruc + 1        

    initdir = 'init'

    column_titles = ['ID', 'potential_energy', 'potential_energy_per_mol', 'Helmholtz_free_energy',
                     'density', 'volume', 
                     'a', 'b', 'c', 'alpha', 'beta', 'gamma', 
                     'sim_method', 'opt_method', 
                     'opt_fmax', 'opt_maxstepsize', 'opt_converged', 'opt_steps', 'time']
    column_titles_summary = ['ID', 'potential_energy', 'potential_energy_per_mol', 
                             'Helmholtz_free_energy', 'density', 'volume', 
                             'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'opt_converged']

    lendfmax = 20
    pd.set_option('display.width', 150)
    pd.set_option('display.max_columns', None)
    pd.options.display.float_format = '{:.2f}'.format

    t1 = time()
    cwdir = os.getcwd()
    print('Optimize molecular crystal in spacegroup {}.\n'.format(spg))
    print('Simulation method: {}.\n'.format(sim_method))

    spgdir = 'spg{}'.format(spg)
    #if not os.path.exists(spgdir):
    #    os.makedirs(spgdir)
    os.chdir(spgdir)
    gendir = 'factor{:.2f}_tfactor{:.2f}'.format(factor, t_factor)
    #if not os.path.exists(gendir):
    #    os.makedirs(gendir)
    os.chdir(gendir)

    basename1 = '{}_spg{}_factor{:.2f}_tfactor{:.2f}'.format(basename, spg, factor, t_factor)
    simdir = '{}'.format(sim_method.lower())
    if not os.path.exists(simdir):
        os.makedirs(simdir)
    os.chdir(simdir)
    optdir = '{}'.format(opt_method.lower())
    if not os.path.exists(optdir):
        os.makedirs(optdir)
    os.chdir(optdir)

    optstepdir = 'fmax{:.4f}_maxsteps{:.5f}_optcycles{:04}'.format(opt_fmax, opt_maxstepsize, opt_maxsteps)
    if not os.path.exists(optstepdir):
        os.makedirs(optstepdir)
    os.chdir(optstepdir)

    basename2 = '{}_{}_{}_fmax{:.4f}_maxsteps{:.5f}_optcycles{:04}'.format(basename1, sim_method, opt_method, 
                                                                           opt_fmax, opt_maxstepsize, opt_maxsteps)

    df = pd.DataFrame(index=[], columns=column_titles)
    for istruc in range(istruc_bgn, istruc_end):

        basename3 = '{}_{:06}'.format(basename2, istruc)
        print(basename3)
        logfile = '{}.log'.format(basename3)
        trajfile = '{}.traj'.format(basename3)
        infile = '{}_{:06}_init.pdb'.format(basename1, istruc)
        lmpinfile = '{}_{:06}_init.input'.format(basename1, istruc)
        lmpdatafile = '{}_{:06}_init.lmp'.format(basename1, istruc)
        pdbin_path = '{}/{}/{}/{}/{}'.format(cwdir, spgdir, gendir, initdir, infile)
        lmpin_path = '{}/{}/{}/{}/{}'.format(cwdir, spgdir, gendir, initdir, lmpinfile)
        lmpdata_path = '{}/{}/{}/{}/{}'.format(cwdir, spgdir, gendir, initdir, lmpdatafile)

        wrkdir = basename3
        if not os.path.exists(wrkdir):
            os.makedirs(wrkdir)
        os.chdir(wrkdir)

        ase_struc = read(pdbin_path, format='proteindatabank')
        if sim_method.upper() == 'LAMMPS':
            shutil.copyfile(lmpin_path, lmpinfile)
            shutil.copyfile(lmpdata_path, lmpdatafile)

        tg1 = time()
        energy, opt_steps, opt_converged = crys_geom_opt(
            ase_struc=ase_struc,
            sim_method=sim_method,
            top_path=lmpinfile,
            kpts=kpts,
            ps_path=ps_path,
            qe_input=qe_input,
            cp2k_input=cp2k_input,
            xtb_hamiltonian=xtb_hamiltonian,
            opt_method=opt_method,
            fmax=opt_fmax,
            opt_maxsteps=opt_maxsteps,
            maxstep=opt_maxstepsize,
            symprec=symprec,
            log_path=logfile,
            traj_path=trajfile,
        )
        tg = time() - tg1
        print('Geometry relaxation took {:3f} seconds.\n'.format(tg))

        #energy = ase_struc.get_potential_energy()
        energy_per_mol = energy / nmols[0]
        pyxtal_struc = pyxtal()
        pyxtal_struc.from_seed(ase_struc)
        density = pyxtal_struc.get_density()
        volume = pyxtal_struc.lattice.volume
        [a, b, c, alpha, beta, gamma] = pyxtal_struc.lattice.get_para(degree=True)

        print('Final lattice:')
        print(pyxtal_struc)
        print('Potential energy:', energy)
        print('Density:', density)

        fe = None
        if free_energy:
            tf1 = time()
            fe = free_energy_calc(ase_struc)
            tf = time() - tf1
            print('Free energy calculation took {:3f} seconds.\n'.format(tf))

        os.chdir('../')
        shutil.copyfile(wrkdir + '/' + logfile, logfile)
        shutil.copyfile(wrkdir + '/' + trajfile, trajfile)
        shutil.rmtree(wrkdir)
              
        record = pd.Series([istruc, energy, energy_per_mol, fe, density, volume, a, b, c, 
                            alpha, beta, gamma, sim_method, opt_method,
                            opt_fmax, opt_maxstepsize, opt_converged, opt_steps, tg],
                           index=df.columns)
                                   
        df = df.append(record, ignore_index=True)
        #print(df)
        outfile = '{}_opt.cif'.format(basename3)
        write(outfile, ase_struc, format='cif')
        outfile = '{}_opt.pdb'.format(basename3)
        with open(outfile, mode='w') as f:
            write_proteindatabank(f, ase_struc, write_arrays=True)
       
        del ase_struc.calc
        del ase_struc, pyxtal_struc

    diffmatfile = 'diffmat_{}_nstruc{:06}-{:06}.csv'.format(basename2, istruc_bgn, istruc_end - 1)
    struclistfile = 'strulist_{}_nstruc{:06}-{:06}.csv'.format(basename2, istruc_bgn, istruc_end - 1)
    compare_crys(infmt='cif', refstrucfile=refstrucfile, comparison=strucdiff_method.upper(),
                 diffmatfile=diffmatfile, struclistfile=struclistfile,
                 verbose=verbose)

    outfile = 'summary_{}_nstruc{:06}-{:06}.csv'.format(basename2, istruc_bgn, istruc_end - 1)
    df['ID'] = df['ID'].astype('int')
    df.to_csv(outfile, float_format='%.3f')

    df_summary = df[column_titles_summary]
    lendf = min([len(df_summary), lendfmax])
    
    print('\n')
    print('Summary of generated structures: ', basename2)
    print(df_summary)
    print('\n')
                    
    print('Summary of generated structures (sorted in energies): ', basename2)
    print(df_summary.sort_values('potential_energy', ascending=True).head(lendf))
    print('\n')

    if free_energy:
        print('Summary of generated structures (sorted in energies): ', basename2)
        print(df_summary.sort_values('Helmholtz_free_energy', ascending=True).head(lendf))
        print('\n')

    print('Summary of generated structures (sorted in densities) ', basename2)
    print(df_summary.sort_values('density', ascending=False).head(lendf))
    print('\n')

    os.chdir(cwdir)

    t2 = time()
    print('All of crystal generation took {:3f} seconds.\n'.format(t2 - t1))
