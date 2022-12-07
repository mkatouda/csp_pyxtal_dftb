import sys
import os
import argparse
import subprocess
from time import time

import yaml
import numpy as np
from ase.io import read, write
from pyxtal import pyxtal
from pyxtal.molecule import generate_molecules

from .critic2 import compare_crys


def pyxtal_to_ase(cs, nxyz, mol_natoms, mol_atomtypes):
    #print('cs.numMols:', cs.numMols, len(cs.numMols), type(cs.numMols))
    #print('cs.numMols.shape:', cs.numMols.shape, len(cs.numMols.shape))

    atoms = cs.to_ase(resort=False)

    #ij = 0
    #ijk = 0
    #for i in range(cs.numMols.shape[0]):
    #    print(i, type(cs.numMols[i]), cs.numMols[i])
    #    print('mol_atomtypes:', mol_atomtypes[i], 'mol_natoms:', mol_natoms[i])
    #    for j in range(cs.numMols[i]):
    #        ij += 1
    #        for k in range(mol_natoms[i]):
    #            ijk += 1
    #            print(i, j, k, ijk, 'res:', ij, 'M{:02}'.format(i+1))

    at = np.concatenate([
        mol_atomtypes[i]
        for i in range(cs.numMols.shape[0])
        for j in range(cs.numMols[i])
    ])
    print(type(at), len(at), at)
    resname = np.array([
        'M{:02}'.format(i+1)
        for i in range(cs.numMols.shape[0])
        for j in range(cs.numMols[i])
        for k in range(mol_natoms[i])
    ])
    print(type(resname), len(resname), resname)

    atoms.set_array('atomtypes', at)
    atoms.set_array('residuenames', resname)
    atoms *= nxyz

    ii = 0
    tmp = []
    for l in range(nxyz[0] * nxyz[1] * nxyz[2]):
        for i in range(cs.numMols.shape[0]):
            tmp.extend([ii+j+1 for j in range(cs.numMols[i]) for k in range(mol_natoms[i])])
            ii += cs.numMols[i]
    resnum = np.array(tmp)
    print(type(resnum), len(resnum), resnum)

    atoms.set_array('residuenumbers', resnum)

    return atoms

def write_proteindatabank(fileobj, images, write_arrays=True):
    """Write images to PDB-file."""

    if hasattr(images, 'get_positions'):
        images = [images]

    rotation = None
    if images[0].get_pbc().any():
        from ase.geometry import cell_to_cellpar, cellpar_to_cell

        currentcell = images[0].get_cell()
        cellpar = cell_to_cellpar(currentcell)
        exportedcell = cellpar_to_cell(cellpar)
        rotation = np.linalg.solve(currentcell, exportedcell)
        # ignoring Z-value, using P1 since we have all atoms defined explicitly
        format = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n'
        fileobj.write(format % (cellpar[0], cellpar[1], cellpar[2],
                                cellpar[3], cellpar[4], cellpar[5]))

    #     1234567 123 6789012345678901   89   67   456789012345678901234567 890
    #format = ('ATOM  %5d %4s %3s     1    %8.3f%8.3f%8.3f%6.2f%6.2f'
    format = ('ATOM  %5d %4s %3s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f'
              '          %2s  \n')

    # RasMol complains if the atom index exceeds 100000. There might
    # be a limit of 5 digit numbers in this field.
    MAXNUM = 100000

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    atomtypes = symbols
    if 'atomtypes' in images[0].arrays:
        atomtypes = images[0].get_array('atomtypes')

    for n, atoms in enumerate(images):
        fileobj.write('MODEL     ' + str(n + 1) + '\n')
        p = atoms.get_positions()
        if 'residuenames' in atoms.arrays:
            resnames = atoms.get_array('residuenames')
        else:
            resnames = np.array(['MOL' for i in range(len(atoms))])
        if 'residuenumbers'  in atoms.arrays:
            resnum = atoms.get_array('residuenumbers')
        else:
            resnum = np.ones(len(atoms), dtype=np.int8)
        occupancy = np.ones(len(atoms))
        bfactor = np.zeros(len(atoms))
        if write_arrays:
            if 'occupancy' in atoms.arrays:
                occupancy = atoms.get_array('occupancy')
            if 'bfactor' in atoms.arrays:
                bfactor = atoms.get_array('bfactor')
        if rotation is not None:
            p = p.dot(rotation)
        for a in range(natoms):
            x, y, z = p[a]
            occ = occupancy[a]
            bf = bfactor[a]
            fileobj.write(format % ((a+1) % MAXNUM, atomtypes[a], resnames[a], resnum[a],
                                    x, y, z, occ, bf, symbols[a].upper()))
        fileobj.write('ENDMDL\n')

def run_cmd(cmd, verbose=False):
    print(' '.join(cmd))
    results = subprocess.run(cmd, capture_output=True, check=True, text=True)
    if verbose: print('stdout:\n', results.stdout, '\nstderr:\n', results.stderr)
    return results

def g16resp_run(mol_path, g16inp_path, charge=0, multi=1, mem=4, nprocshared=4):
    atoms = read(mol_path, format='mol')

    s = '%mem={}GB\n%nprocshared={}\n'.format(mem, nprocshared)
    s += '#p hf/6-31g(d) pop=mk iop(6/33=2,6/42=6) scf=tight\n\n'
    s += 'HF/6-31g(d) RESP Charge\n\n'
    s += '{:d} {:d}\n'.format(charge, multi)
    for cs, xyz in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        s += '{:<2s} {:>15.8f} {:>15.8f} {:>15.8f}\n'\
            .format(cs, xyz[0], xyz[1], xyz[2])
    s += '\n'    
    print(s)

    with open(g16inp_path, mode='w') as f:
        f.write(s)

    cmd = ['g16', g16inp_path]
    run_cmd(cmd)

def psi4resp_run(mol_path, respchg_path='resp_charge.out'):
    import psi4
    import resp

    atoms = read(mol_path, format='mol')
    s = ''
    for cs, xyz in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        s += '{:<2s} {:>15.8f} {:>15.8f} {:>15.8f}\n'\
            .format(cs, xyz[0], xyz[1], xyz[2])
    #print(s)

    mol = psi4.geometry(s)
    mol.update_geometry()

    options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
               'VDW_POINT_DENSITY'  : 1.0,
               'RESP_A'             : 0.0005,
               'RESP_B'             : 0.1,
           }

    # Call for RESP fit
    charges1 = resp.resp([mol], options)
    print('Electrostatic Potential Charges')
    print(charges1[0])
    print('Restrained Electrostatic Potential Charges')
    print(charges1[1])

    s = '\n'.join(['{:.6f}'.format(c) for c in charges1[1]])
    with open(respchg_path, mode='w') as f:
        f.write(s)

    return charges1[1]

def ac_run(pm, ff='gaff2', charge_model='bcc'):
    mol_atomtypes = []
    mol_natoms = []
    print('ac_run\n', pm)
    for i, p in enumerate(pm):
        resname = 'M{:02}'.format(i+1)
        molfile = resname + '.mol'
        mol2file = resname + '.mol2'
        pdbfile = resname + '.pdb'
        frcmodfile = resname + '.frcmod'
        g16infile = resname + '.com'
        g16outfile = resname + '.log'
        respfile = resname + '_resp_charge.out'
        with open(molfile, mode='w') as f:
            f.write(p.rdkit_mb)

        if charge_model.upper() == 'RESP':
            g16resp_run(molfile, g16infile)
            cmd = ['antechamber', '-fi', 'gout', '-fo', 'mol2',  '-at', ff, '-c', 'resp', '-i', g16outfile, '-o', mol2file]
            #psi4resp_run(molfile, respfile)
            #cmd = ['antechamber', '-fi', 'mdl', '-fo', 'mol2', '-at', ff, '-c', 'rc', '-cf', respfile, '-rn', resname,'-i', molfile, '-o', mol2file]
        else:
            cmd = ['antechamber', '-fi', 'mdl', '-fo', 'mol2', '-at', ff, '-c', charge_model, '-rn', resname,'-i', molfile, '-o', mol2file]
        run_cmd(cmd)

        cmd = ['antechamber', '-fi', 'mol2', '-fo', 'pdb', '-i', mol2file, '-o', pdbfile]
        run_cmd(cmd)

        cmd = ['parmchk2', '-i', mol2file, '-o', frcmodfile, '-f', 'mol2']
        run_cmd(cmd)

        atoms = read(pdbfile, format='proteindatabank')
        at = atoms.get_array('atomtypes')
        mol_atomtypes.append(at)
        mol_natoms.append(len(at))

    return mol_atomtypes, mol_natoms

def tleap_run(pm, cryspdb_path='model.pdb', leapin_path='leap.in'):
    root, ext = os.path.splitext(cryspdb_path)
    prmtop_path = root + '.prmtop'
    crd_path = root + '.crd'

    s = 'source leaprc.gaff\n\n'
    for i, p in enumerate(pm):
        resname = 'M{:02}'.format(i+1)
        frcmod_path = resname + '.frcmod'
        s += 'loadamberparams {}\n'.format(frcmod_path)
    s += '\n'
    for i, p in enumerate(pm):
        resname = 'M{:02}'.format(i+1)
        mol2_path = resname + '.mol2'
        s += '{} = loadmol2 {}\n'.format(resname, mol2_path)
        s += 'charge {}\n'.format(resname)
    s += '\n'
    s += 'crys = loadPDB {}\n'.format(cryspdb_path)
    s += 'charge crys\n\n'
    s += 'saveAmberParm crys {} {}\n'.format(prmtop_path, crd_path)
    s += 'quit\n'

    with open(leapin_path, mode='w') as f:
        f.write(s)

    cmd = ['tleap', '-f', leapin_path]
    print(' '.join(cmd))
    run_cmd(cmd)

    cell = []
    with open(cryspdb_path) as f:
        for l in f.readlines():
            if 'CRYST1' in l: 
                cell = l.split()[1:7]

    if len(cell) == 6:
        cmd = ['ChBox', '-c', crd_path, '-o', crd_path,
               '-X', cell[0], '-Y', cell[1], '-Z', cell[2],
               '-al', cell[3], '-bt', cell[4], '-gm', cell[5]]
        print(' '.join(cmd))
        run_cmd(cmd)

    return prmtop_path, crd_path

def parmed_run(amb_top_path, amb_crd_path):
    root, ext = os.path.splitext(amb_top_path)
    gmx_top_path = root + '.top'
    gmx_crd_path = root + '.gro'

    # convert topology from Amber to Gromacs
    import parmed as pmd
    amb_top = pmd.load_file(amb_top_path, xyz=amb_crd_path)
    amb_top.save(gmx_top_path, overwrite=True)
    amb_top.save(gmx_crd_path, overwrite=True)

    return gmx_top_path, gmx_crd_path

def intermol_run(amb_prmtop_path, amb_inpcrd_path):
    from intermol.convert import _load_amber, _save_lammps

    system, prefix, prmtop_in, crd_in, amb_structure \
        = _load_amber(amber_files=[amb_prmtop_path, amb_inpcrd_path])

    oname = '{0}'.format(prefix)
    output_status = dict()
    lmpin_path = '{0}.input'.format(oname)
    lmp_settings = {}
    lmp_settings['lmp_settings'] = 'pair_style lj/cut/coul/long 9.0 9.0\npair_modify tail yes\nkspace_style pppm 1e-8\n\n'
    _save_lammps(system, oname, output_status, lmp_settings)
    # Specify the output energies that we are interested in.
    energy_terms = " ".join(['ebond', 'eangle', 'edihed', 'eimp',
                             'epair', 'evdwl', 'ecoul', 'elong',
                             'etail', 'pe', 'pxx'])

    s = ''
    with open(lmpin_path) as f:
        for l in f.readlines():
            if 'thermo_style' in l:
                s += 'thermo_style custom {0}\n'.format(energy_terms)
            else:
                s += l

    with open(lmpin_path, mode='w') as f:
        f.write(s)

    return lmpin_path

def gen_random_molcrys(basename, mols, nmols, spg, nstruc=100, 
                       factor=1.0, t_factor=1.0, use_hall=False,
                       nstruc_try=500, istruc_bgn=1, istruc_end=100,
                       nxyz=[1, 1, 1], ff='gaff2', charge_model='bcc',
                       strucdiff_method='POWDER', verbose=False):

    cwdir = os.getcwd()
    print('Current working directory: {}'.format(cwdir))

    if istruc_end > nstruc + 1:
        istruc_end = nstruc + 1        

    t1 = time()
    print('Genetate molecular crystal in spacegroup {}.\n'.format(spg))

    spgdir = 'spg{}'.format(spg)
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

    pm = []
    for idx, mol in enumerate(mols):
        smile = mol.replace('.smi', '')
        print('mols', idx, 'smiles:', smile)
        pxmols = generate_molecules(smile, wps=None, N_iter=4, N_conf=10, tol=0.5)
        engs = [m.energy for m in pxmols]
        imin = engs.index(min(engs))
        print('minimun energy index:', imin, engs[imin])
        if verbose:
            for m in pxmols:
                print(m.energy, m.pga.sch_symbol)
        pm.append(pxmols[imin])

    if verbose:
        print('pyxtal_molecules:')
        print('lenpm', len(pm))
        for p in pm:
            print(p.name)
            print(vars(p))
            print(p.rdkit_mb)

    mol_atomtypes, mol_natoms = ac_run(pm, ff, charge_model)

    basename1 = '{}_spg{}_factor{:.2f}_tfactor{:.2f}'.format(basename, spg, factor, t_factor)
    istruc = 0
    for istruc_try in range(1, nstruc_try+1):

        print('Trial {}'.format(istruc_try))

        pyxtal_struc = pyxtal(molecular=True)
        pyxtal_struc.from_random(
            dim=3,
            group=spg,
            species=mols,
            numIons=nmols,
            factor=factor, 
            t_factor=t_factor,
            use_hall=use_hall
        )
        if pyxtal_struc.valid:
            istruc += 1
            print('Trial {} is valid, and add structure {}\n'.format(istruc_try, istruc))

            energy = None; energy_per_mol = None
            density = pyxtal_struc.get_density()
            opt_converged = None; opt_steps = None; tg = None

            print('Initial lattice:')
            print(pyxtal_struc)
            print('Density:', density)
            if verbose: 
                print(vars(pyxtal_struc))
                for pyxtal_mol in pyxtal_struc.molecules:
                    print(vars(pyxtal_mol))
                    print(pyxtal_mol.rdkit_mb)

            # Create ASE atoms object
            ase_struc = pyxtal_to_ase(pyxtal_struc, nxyz, mol_natoms, mol_atomtypes)
            if verbose: print('atoms:', ase_struc)
            outfile = '{}_{:06}_init.cif'.format(basename1, istruc)
            write(outfile, ase_struc, format='cif')
            outfile = '{}_{:06}_init.pdb'.format(basename1, istruc)
            with open(outfile, mode='w') as f:
                write_proteindatabank(f, ase_struc, write_arrays=True)

            # Create Amber topology
            amb_prmtop_path, amb_inpcrd_path = tleap_run(pm, outfile)

            # Create Gromacs topology
            gmx_top_path, gmx_crd_path = parmed_run(amb_prmtop_path, amb_inpcrd_path)

            # Create Lammps topology
            lmpin_path = intermol_run(amb_prmtop_path, amb_inpcrd_path)

            if istruc == nstruc:
                break

    # Compare crystal structure similarity
    if strucdiff_method.upper() in ['POWDER', 'RDF', 'AMD']:
        compare_crys(infmt='cif', refstrucfile=None, comparison=strucdiff_method.upper(),
                     diffmatfile='diffmat.csv', struclistfile='struclist.csv',
                     verbose=verbose)

    os.chdir(cwdir)
    t2 = time()
    print('All of crystal generation took {:3f} seconds.\n'.format(t2 - t1))

def get_parser():
    parser = argparse.ArgumentParser(
        description='PyXtal DFTB crystal structure prediction code',
        usage=f"python {os.path.basename(__file__)} -c CONFIG_FILE"
    )
    parser.add_argument(
        "-i", "--input", type=str, #required=True,
        help="path to a yaml config file"
    )

    return parser.parse_args()

def set_default_config(conf):

    conf.setdefault('basename', 'benzene')
    conf.setdefault('mols', ['c1ccccc1.smi'])
    conf.setdefault('nmols', [4])
    conf.setdefault('spg', 14) # 'P21/c'

    conf.setdefault('nstruc', 10)
    conf.setdefault('factor', 1.0)
    conf.setdefault('tfactor', 1.0)
    conf.setdefault('use_hall', False)

    #conf.setdefault('nstruc_try', 500)
    conf.setdefault('istbgn', 1)
    conf.setdefault('istend', 1)

    conf.setdefault('nxyz', [1, 1, 1])
    conf.setdefault('ff', 'gaff2')
    conf.setdefault('charge_model', 'bcc')

    strucdiff_method('strucdiff_method', 'POWDER')

    conf.setdefault('verbose', True)

    return conf

def main():
    args = get_parser()
    print(args)

    conf = {}
    if args.input is not None and os.path.isfile(args.input):
        with open(args.input, "r") as f:
            conf = yaml.load(f, Loader=yaml.SafeLoader)
    set_default_config(conf)
    print(conf)

    basename = conf['basename']
    mols = conf['mols']
    nmols = conf['nmols']
    spg = conf['spg']

    nstruc = conf['nstruc']
    factor =  conf['factor']
    t_factor = conf['tfactor']
    use_hall = conf['use_hall']

    #nstruc_try = conf['nstruc_try']
    nstruc_try = nstruc * 5
    istruc_bgn = conf['istbgn']
    istruc_end = conf['istend'] + 1

    nxyz = conf['nxyz']
    ff = conf['ff']
    charge_model = conf['charge_model']

    strucdiff_method = conf['strucdiff_method']

    verbose = conf['verbose']

    gen_random_molcrys(basename, mols, nmols, spg,
                       nstruc, factor, t_factor, use_hall,
                       nstruc_try, istruc_bgn, istruc_end,
                       nxyz, ff, charge_model,
                       strucdiff_method, verbose)
                       

if __name__ == '__main__':
    main()
