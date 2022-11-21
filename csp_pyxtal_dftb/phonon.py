from ase.phonons import Phonons
from ase.thermochemistry import CrystalThermo


def free_energy_calc(atoms, supercell=[1, 1, 1], ph_delta=0.05,
                     kpts=[1, 1, 1], npts=3000, dos_delta=5e-4, temp=298.15):

    pe = atoms.get_potential_energy()
    print('Potential energy:', pe)

    # Perform phonon analysis
    ph = Phonons(atoms, atoms.calc, supercell=supercell, delta=ph_delta)
    ph.run()
    ph.read(acoustic=True)
    phonon_energies, phonon_DOS = ph.dos(kpts=kpts, npts=npts, delta=dos_delta)
    #print('phonon_energies:', phonon_energies)
    #print('phonon_DOS:', phonon_DOS)

    # Calculate the Helmholtz free energy
    thermo = CrystalThermo(phonon_energies=phonon_energies,
                           phonon_DOS=phonon_DOS,
                           potentialenergy=pe,
                           formula_units=4)
    F = thermo.get_helmholtz_energy(temperature=temp)
    print('Helmholtz free energy({:.2f}): {:.3f}'.format(temp, F))

    return F
