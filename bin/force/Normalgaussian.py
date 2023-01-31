from rdkit import Chem
from rdkit.Chem import AllChem
import os
import json
import subprocess


def _solvemol(filename):

    mol = Chem.MolFromMolFile('./model/' + filename + '.mol')
    mol = Chem.MolToSmiles(mol)
    Mol = Chem.AddHs(Chem.MolFromSmiles(mol))

    return Mol


def _information(filename):

    Mol = _solvemol(filename)
    Numatom = Mol.GetNumAtoms()
    formcharge = 0
    for atom in Mol.GetAtoms():
        formcharge += atom.GetFormalCharge()

    return Numatom, formcharge


def gaussianfile(filename, mem=30):

    Mol = _solvemol(filename)
    _, formcharge = _information(filename)
    AllChem.EmbedMolecule(Mol, randomSeed=0xf00d, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(Mol)

    folder = os.path.exists('./topper')
    if not folder:
        os.makedirs('./topper')

    f_sys = open('./control/system.json', )  # Read nthread
    system = json.load(f_sys)
    nthread = 4
    if 'CPU' in system:
        nthread = system['CPU']

    f_out = open('./topper/' + filename + '.gjf', 'w')
    print('%mem=', mem, 'GB', sep='', end='', file=f_out)
    print(file=f_out)
    print('%nproc=', nthread, sep='', end='', file=f_out)
    print(file=f_out)
    print(r'%chk=', filename, 'em.chk', sep='', file=f_out)
    print('#p PM6D3  opt ', file=f_out)
    print(file=f_out)
    print(filename, file=f_out)
    print(file=f_out)
    print(formcharge, '1', file=f_out)
    for atom in Mol.GetAtoms():
        atom_idx = atom.GetIdx()
        atom_sym = atom.GetSymbol()
        x, y, z = Mol.GetConformer().GetAtomPosition(atom_idx)
        print(atom_sym.rjust(2),
              str('{:.8f}'.format(x)).rjust(27),
              str('{:.8f}'.format(y)).rjust(13),
              str('{:.8f}'.format(z)).rjust(13),
              file=f_out)
    print(file=f_out)
    print('--Link1--', file=f_out)
    print('%mem=', mem, 'GB', sep='', end='', file=f_out)
    print(file=f_out)
    print('%nproc=', nthread, sep='', end='', file=f_out)
    print(file=f_out)
    print(r'%oldchk=', filename, 'em.chk', sep='', file=f_out)
    print(r'%chk=', filename, 'step1.chk', sep='', file=f_out)
    print(
        '#p B3LYP/6-31G* em=GD3BJ opt(maxstep=20,notrust,GDIIS,loose,maxcyc=120) nosymm guess=read geom=allcheck int=fine',
        file=f_out)
    print(file=f_out)
    print('--Link1--', file=f_out)
    print('%mem=', mem, 'GB', sep='', end='', file=f_out)
    print(file=f_out)
    print('%nproc=', nthread, sep='', end='', file=f_out)
    print(file=f_out)
    print(r'%oldchk=', filename, 'step1.chk', sep='', file=f_out)
    print(r'%chk=', filename, 'step2.chk', sep='', file=f_out)
    print(
        '#p B3LYP/genecp em=GD3BJ opt(maxstep=5,readfc,notrust,GDIIS,maxcyc=120) nosymm guess=read geom=allcheck freq int=fine',
        file=f_out)
    print(file=f_out)
    print(
        '-H -Li -Be -B -C -N -O -F -Na -Mg -Al -Si -P -S -Cl',
        '6-311G**',
        '****',
        '-K -Ca -Sr -Ti -V -Cr -Mn -Fe -Co -Ni -Cu -Zn -Ga -Ge -As -Se -Br -Mo -Ru -Rh -Pd -Ag -Cd -I -Pt -Au -Hg',
        'SDD',
        '****',
        sep='\n',
        end='\n',
        file=f_out)
    print(file=f_out)
    print(
        '-K -Ca -Sr -Ti -V -Cr -Mn -Fe -Co -Ni -Cu -Zn -Ga -Ge -As -Se -Br -Mo -Ru -Rh -Pd -Ag -Cd -I -Pt -Au -Hg',
        'SDD',
        sep='\n',
        end='\n',
        file=f_out)
    print(file=f_out)
    print(file=f_out)


def gaussiancalculate(filename):
    step1 = subprocess.Popen('g16 ' + filename + '.gjf',
                             shell=True,
                             cwd='./topper')
    step1.wait()
    if step1.returncode != 0:
        print('Gaussian caculating is wrong (' + filename + ')')
        return False
    else:
        return True


def itpcalculate(filename):
    step2 = subprocess.Popen('formchk ' + filename + 'step2.chk',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True,
                             cwd='./topper')
    step2.wait()

    step3 = subprocess.Popen('obabel -ifchk ' + filename + 'step2.fchk' +
                             ' -omol2 -O ' + filename + '.mol2',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True,
                             cwd='./topper')
    step3.wait()

    sobtop = open('./sobtop/' + filename + '.sh', 'w')
    print('chmod 777 sobtop', 'chmod 777 atomtype', sep='\n', file=sobtop)
    print('./sobtop  <<EOF', file=sobtop)
    print('../topper/', filename + '.mol2', sep='', file=sobtop)
    print(1, 2, 7, sep='\n', file=sobtop)
    print('../topper/', filename + 'step2.fchk', sep='', file=sobtop)
    print('../topper/', filename + '.top', sep='', file=sobtop)
    print('../topper/', filename + '.itp', sep='', file=sobtop)
    print(0, file=sobtop)
    print('EOF', end='', file=sobtop)
    sobtop.close()

    step4 = subprocess.Popen('sh ' + filename + '.sh',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True,
                             cwd='./sobtop')
    step4.wait()
    if not os.access('./topper/' + filename + '.itp', os.F_OK):
        print('Interaction caculating is wrong (' + filename + ')')
        return False
    else:
        return True


def chargecalculate(filename):

    if os.access('./topper/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./topper/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_out.close()

        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./topper')
        step4.wait()
        if not os.access('./topper/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            return False
        else:
            return True
