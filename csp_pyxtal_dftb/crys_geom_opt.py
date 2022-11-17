oimport sys
import os
from time import time
from ase.constraints import ExpCellFilter
from ase.spacegroup.symmetrize import FixSymmetry
from ase.io import read, write
from ase.io.trajectory import Trajectory

def get_ase_calculator(qc_method):
    if qc_method.upper() == 'DFTB':
        from ase.calculators.dftb import Dftb
        calc = Dftb(kpts=(1,1,1),
                Hamiltonian_='DFTB',
                Hamiltonian_SCC='Yes',
                Hamiltonian_SCCTolerance=1e-6,
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_P='d',
                Hamiltonian_MaxAngularMomentum_Br='d',
                Hamiltonian_MaxAngularMomentum_Cl='d',
                Hamiltonian_MaxAngularMomentum_F='p',
                Hamiltonian_MaxAngularMomentum_I='d',
                Hamiltonian_Dispersion_='DftD3',
                Hamiltonian_Dispersion_Damping_='BeckeJohnson',
                Hamiltonian_Dispersion_Damping_a1=0.746,
                Hamiltonian_Dispersion_Damping_a2=4.191,
                Hamiltonian_Dispersion_s6=1.0,
                Hamiltonian_Dispersion_s8=3.209
            )
    else:
        from xtb.ase.calculator import XTB
        calc = XTB(method="GFN0-xTB")

    return calc

def crys_geom_opt(ase_struc, basename, method='LBFGS', fmax=5e-2, opt_maxsteps=100, 
                  logfile='-', trajfile='crys_geom_opt.traj', maxstep=0.01, symprec=0.1):

    #t1 = time()
    ase_struc.set_constraint(FixSymmetry(ase_struc, symprec=symprec, verbose=True))
    if method == 'FIRE':
        from ase.optimize import FIRE
        opt = FIRE(ExpCellFilter(ase_struc), logfile=logfile, trajectory=trajfile, maxstep=maxstep)
    elif method == 'LBFGS':
        from ase.optimize import LBFGS
        opt = LBFGS(ExpCellFilter(ase_struc), logfile=logfile, trajectory=trajfile, maxstep=maxstep)
    elif method == 'LBFGSLineSearch':
        from ase.optimize import LBFGSLineSearch
        opt = LBFGSLineSearch(ExpCellFilter(ase_struc), logfile=logfile, trajectory=trajfile, maxstep=maxstep)
    elif method == 'BFGS':
        from ase.optimize import BFGS
        opt = BFGS(ExpCellFilter(ase_struc), logfile=logfile, trajectory=trajfile, maxstep=maxstep)
    elif method == 'GPMin':
        from ase.optimize import GPMin
        opt = GPMin(ExpCellFilter(ase_struc), logfile=logfile, trajectory=trajfile)
    elif method == 'MDMin':
        from ase.optimize import MDMin
        opt = MDMin(ase_struc, logfile=logfile, trajectory=trajfile)
        opt.run(fmax=fmax, steps=opt_steps)
        from ase.optimize import LBFGS
        opt = LBFGS(ExpCellFilter(ase_struc), logfile=logfile, trajectory=trajfile, maxstep=maxstep)
    elif method =='pyberny':
        from ase.optimize import Berny
        opt = Berny(ase_struc, logfile=logfile, trajectory=trajfile)
        opt.run(fmax=fmax, steps=opt_steps)
        from ase.optimize import LBFGS
        opt = LBFGS(ExpCellFilter(ase_struc), logfile=logfile, trajectory=trajfile, maxstep=maxstep)
        #from ase.optimize import FIRE
        #opt = FIRE(ExpCellFilter(ase_struc), logfile=logfile, trajectory=trajfile, maxstep=maxstep)
    opt.run(fmax=fmax, steps=opt_maxsteps)
    #t2 = time()
    #print('Relaxation took {:3f} seconds.\n'.format(t2 - t1))

    opt_converged = opt.converged()
    print('Opt converged: ', opt_converged)
    opt_steps = opt.nsteps #opt.get_number_of_steps()
    print('Number of optimized step', opt_steps)

    energy = ase_struc.get_potential_energy()
    ase_struc.set_constraint()

    return ase_struc, energy, opt_steps, opt_converged

def get_geompng(atoms, outfile):
    write(outfile, atoms, format='png')

def get_trajoutput(trajfile, outfile, format):
    traj = Trajectory(trajfile)
    write(outfile, traj, format=format)

def crys_geom_opt_test():

    infilepath = sys.argv[1]
    basename = os.path.splitext(os.path.basename(infilepath))[0]
    infile = basename + '.cif'

    qc_method = 'DFTB'
    #qc_method = 'xTB'
    #opt_methods = ['LBFGS', 'FIRE']
    opt_methods = ['LBFGS']
    #opt_methods = ['FIRE']
    #opt_methods = ['GPMin']
    #opt_methods = ['LBFGSLineSearch']
    fmax = 1.0e-2
    opt_maxsteps = 300
    #maxsteps = [0.2, 0.1, 0.05, 0.01, 0.005, 0.001]
    #maxsteps = [0.1, 0.01]
    maxsteps = [0.01]

    symprec = 0.1

    atoms = read(infile)
    outfile = '{}.png'.format(basename)
    write(outfile, atoms, format='png')

    calc = get_ase_calculator(qc_method)

    csvfile = '{}_{}.csv'.format(basename, qc_method)
    foutcsv = open(csvfile, 'w')
    #foutcsv.write('qc_method,opt_method,fmax,maxstep,time,opt_converged,opt_steps,energy\n')
    foutcsv.write('qc_method,opt_method,fmax,maxstep,time,opt_converged,opt_steps,energy,volume,a,b,c,alpha,beta,gamma\n')
    for opt_method in opt_methods:
        if opt_method == 'LBFGS' or opt_method == 'LBFGSLineSearch' or opt_method == 'FIRE':
            for maxstep in maxsteps:

                atoms = read(infile)
                print('Initial lattice:')
                bl = atoms.cell.get_bravais_lattice()
                print(bl)

                atoms.calc = calc

                logfile = '{}_{}_{}_fmax{:.4f}_maxstep{:.4f}.log'.format(basename, qc_method, opt_method, fmax, maxstep)
                trajfile = '{}_{}_{}_fmax{:.4f}_maxstep{:.4f}.traj'.format(basename, qc_method, opt_method, fmax, maxstep)
                t1 = time()
                atoms, energy, opt_steps, opt_converged = crys_geom_opt(atoms, basename, 
                                                                        opt_method, fmax, opt_maxsteps,
                                                                        logfile, trajfile, maxstep)
                t2 = time()
                t_opt = t2 -t1
                print('Relaxation took {:3f} seconds.\n'.format(t_opt))


                print('Traj file: ', trajfile)
                print('Final lattice:')
                bl = atoms.cell.get_bravais_lattice()
                print(bl)
                print('Total energy: ', atoms.get_potential_energy())

                [a, b, c] = atoms.cell.lengths()
                [alpha, beta, gamma] = atoms.cell.angles()
                volume = atoms.cell.volume
                #foutcsv.write('{},{},{},{},{},{},{},{}\n'
                #              .format(qc_method, opt_method, fmax, maxstep, opt_converged, t_opt, opt_steps, energy))
                foutcsv.write('{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'
                              .format(qc_method, opt_method, fmax, maxstep, opt_converged, t_opt, opt_steps, energy,
                                      volume, a, b, c, alpha, beta, gamma))

                outfile = '{}_{}_{}_fmax{:.4f}_maxstep{:.4f}_opt.cif'.format(basename, qc_method, opt_method, fmax, maxstep)
                write(outfile, atoms, format='cif')

                outfile = '{}_{}_{}_fmax{:.4f}_maxstep{:.4f}_opt.png'.format(basename, qc_method, opt_method, fmax, maxstep)
                write(outfile, atoms, format='png')

                #outfile = '{}_{}_{}_fmax{:.4f}_maxstep{:.4f}.pdb'.format(basename, qc_method, opt_method, fmax, maxstep)
                #get_trajoutput(trajfile, outfile, 'proteindatabank')
                #outfile = '{}_{}_{}_fmax{:.4f}_maxstep{:.4f}.mp4'.format(basename, qc_method, opt_method, fmax, maxstep)
                #get_trajoutput(trajfile, outfile, 'mp4')


        elif opt_method == 'GPMin':

            maxstep = None

            atoms = read(infile)
            print('Initial lattice:')
            bl = atoms.cell.get_bravais_lattice()
            print(bl)

            atoms.calc = calc

            logfile = '{}_{}_{}_fmax{:.4f}.log'.format(basename, qc_method, opt_method, fmax)
            trajfile = '{}_{}_{}_fmax{:.4f}.traj'.format(basename, qc_method, opt_method, fmax)
            t1 = time()
            atoms, energy, opt_steps, opt_converged = crys_geom_opt(atoms, basename, 
                                                                    opt_method, fmax, opt_maxsteps,
                                                                    logfile, trajfile, maxstep, symprec)
            t2 = time()
            t_opt = t2 -t1
            print('Relaxation took {:3f} seconds.\n'.format(t_opt))

            print('Traj file: ', trajfile)
            print('Final lattice:')
            bl = atoms.cell.get_bravais_lattice()
            print(bl)
            print('Total energy: ', atoms.get_potential_energy())
            
            [a, b, c] = atoms.cell.lengths()
            [alpha, beta, gamma] = atoms.cell.angles()
            volume = atoms.cell.volume
            #foutcsv.write('{},{},{},{},{},{},{},{}\n'
            #              .format(qc_method, opt_method, fmax, maxstep, opt_converged, t_opt, opt_steps, energy))
            foutcsv.write('{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'
                          .format(qc_method, opt_method, fmax, maxstep, opt_converged, t_opt, opt_steps, energy,
                                  volume, a, b, c, alpha, beta, gamma))

            outfile = '{}_{}_{}_fmax{:.4f}_opt.cif'.format(basename, qc_method, opt_method, fmax)
            write(outfile, atoms, format='cif')

            outfile = '{}_{}_{}_fmax{:.4f}_opt.png'.format(basename, qc_method, opt_method, fmax)
            write(outfile, atoms, format='png')

            #outfile = '{}_{}_{}_fmax{:.4f}.pdb'.format(basename, qc_method, opt_method, fmax)
            #get_trajoutput(trajfile, outfile, 'proteindatabank')
            #get_trajpdb(trajfile, outfile)
            #outfile = '{}_{}_{}_fmax{:.4f}.mp4'.format(basename, qc_method, opt_method, fmax)
            #get_trajoutput(trajfile, outfile, 'mp4')


    foutcsv.close()

if __name__ == "__main__":
    crys_geom_opt_test()
