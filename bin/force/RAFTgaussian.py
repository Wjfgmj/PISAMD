import os
from rdkit import Chem
import subprocess


def _solvemol():
    Mol = []
    model = ['I.mol', 'R.mol', 'MA.mol', 'MB.mol', 'C.mol']
    for i in range(len(model)):
        mol = Chem.MolFromMolFile('./model/' + model[i])
        mol = Chem.MolToSmiles(mol)
        mol = Chem.AddHs(Chem.MolFromSmiles(mol))
        Mol.append(mol)
    return Mol


def _Numatom():
    Numatom = []
    Mols = _solvemol()
    for mol in Mols:
        temp = mol.GetNumAtoms()
        Numatom.append(temp)
    return Numatom


def _Formcharge():
    Mols = _solvemol()
    formcharge = []
    for Mol in Mols:
        temp = 0
        for atom in Mol.GetAtoms():
            temp += atom.GetFormalCharge()
        formcharge.append(temp)
    return formcharge


def gaussianfile(Mol, filename, mem, nthread):

    folder = os.path.exists('./gaussian')
    if not folder:
        os.makedirs('./gaussian')
    f_out = open('./gaussian/' + filename + '.gjf', 'w')
    print('%mem=', mem, 'GB', sep='', end='', file=f_out)
    print(file=f_out)
    print('%nproc=', nthread, sep='', end='', file=f_out)
    print(file=f_out)
    print(r'%chk=', filename, 'em.chk', sep='', file=f_out)
    print('#p PM6D3  opt ', file=f_out)
    print(file=f_out)
    print(filename, file=f_out)
    print(file=f_out)

    if filename == 'MA':
        charge = _Formcharge()
        print(charge[2], '1', file=f_out)

    if filename == 'MB':
        charge = _Formcharge()
        print(charge[3], '1', file=f_out)

    if filename == 'IC':
        charge = _Formcharge()
        print(charge[0] + charge[4] + 2, '1', file=f_out)

    if filename == 'RC':
        charge = _Formcharge()
        print(charge[1] + charge[4] + 2, '1', file=f_out)

    if filename == 'IMAC':
        charge = _Formcharge()
        print(charge[0] + charge[2] + charge[4] + 2, '1', file=f_out)

    if filename == 'RMAC':
        charge = _Formcharge()
        print(charge[1] + charge[2] + charge[4] + 2, '1', file=f_out)

    if filename == 'IMBC':
        charge = _Formcharge()
        print(charge[0] + charge[3] + charge[4] + 2, '1', file=f_out)

    if filename == 'RMBC':
        charge = _Formcharge()
        print(charge[1] + charge[3] + charge[4] + 2, '1', file=f_out)

    if filename == 'IMA2C':
        charge = _Formcharge()
        print(charge[0] + charge[2] * 2 + charge[4] + 2, '1', file=f_out)

    if filename == 'IMB2C':
        charge = _Formcharge()
        print(charge[0] + charge[3] * 2 + charge[4] + 2, '1', file=f_out)

    if filename == 'RMA3C':
        charge = _Formcharge()
        print(charge[1] + charge[2] * 3 + charge[4] + 2, '1', file=f_out)

    if filename == 'RMB3C':
        charge = _Formcharge()
        print(charge[1] + charge[3] * 3 + charge[4] + 2, '1', file=f_out)

    if filename == 'RMAMB2C':
        charge = _Formcharge()
        print(charge[1] + charge[2] + charge[3] * 2 + charge[4] + 2,
              '1',
              file=f_out)

    if filename == 'RMA2MBC':
        charge = _Formcharge()
        print(charge[1] + charge[2] * 2 + charge[3] + charge[4] + 2,
              '1',
              file=f_out)

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
        '#p B3LYP/6-31G* em=GD3BJ opt(maxstep=20,notrust,GDIIS,loose) nosymm guess=read geom=allcheck',
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
        '#p B3LYP/genecp em=GD3BJ opt(maxstep=5,readfc,notrust,GDIIS) nosymm guess=read geom=allcheck freq',
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
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             cwd='./gaussian')
    step1.wait()
    if step1.returncode != 0:
        return False
    else:
        return True


def itpcalculate(filename):
    step2 = subprocess.Popen('formchk ' + filename + 'step2.chk',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True,
                             cwd='./gaussian')
    step2.wait()

    step3 = subprocess.Popen('obabel -ifchk ' + filename + 'step2.fchk' +
                             ' -omol2 -O ' + filename + '.mol2',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True,
                             cwd='./gaussian')
    step3.wait()

    sobtop = open('./sobtop/' + filename + '.sh', 'w')
    print('chmod 777 sobtop', 'chmod 777 atomtype', sep='\n', file=sobtop)
    print('./sobtop  <<EOF', file=sobtop)
    print('../gaussian/', filename + '.mol2', sep='', file=sobtop)
    print(1, 2, 3, sep='\n', file=sobtop)
    print('../gaussian/', filename + 'step2.fchk', sep='', file=sobtop)
    print('../gaussian/', filename + '.top', sep='', file=sobtop)
    print('../gaussian/', filename + '.itp', sep='', file=sobtop)
    print(0, file=sobtop)
    print('EOF', end='', file=sobtop)
    sobtop.close()

    step4 = subprocess.Popen('sh ' + filename + '.sh',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True,
                             cwd='./sobtop')
    step4.wait()
    if not os.access('./gaussian/' + filename + '.itp', os.F_OK):
        print('Interaction caculating is wrong (' + filename + ')')
        return False
    else:
        return True


def chargecalculate():
    Numatom = _Numatom()
    charge = _Formcharge()
    cout = 0
    total = 0

    # I* state_1
    filename = 'IMAC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print('1', '-', Numatom[0], sep='', end=' ', file=f_text)
        print(charge[0] + 1, file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/I_state1.chg')
            total += 1

    # I* state_2
    filename = 'IC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print('1', '-', Numatom[0], sep='', end=' ', file=f_text)
        print(charge[0] + 1, file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/I_state2.chg')
            total += 1

    # R* state_1
    filename = 'RMAC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print('1', '-', Numatom[1], sep='', end=' ', file=f_text)
        print(charge[1] + 1, file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/R_state1.chg')
            total += 1

    # R* state_2
    filename = 'RC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print('1', '-', Numatom[1], sep='', end=' ', file=f_text)
        print(charge[1] + 1, file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/R_state2.chg')
            total += 1

    # C* state_1
    filename = 'RC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + 1,
              '-',
              Numatom[1] + Numatom[4],
              sep='',
              end=' ',
              file=f_text)
        print(charge[4] + 1, end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/C_state1.chg')
            total += 1

    # C* state_2
    filename = 'RMAC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + Numatom[2] + 1,
              '-',
              Numatom[1] + Numatom[2] + Numatom[4],
              sep='',
              end=' ',
              file=f_text)
        print(charge[4] + 1, end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/C_state2.chg')
            total += 1

    # C* state_3
    filename = 'RMBC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + Numatom[3] + 1,
              '-',
              Numatom[1] + Numatom[3] + Numatom[4],
              sep='',
              end=' ',
              file=f_text)
        print(charge[4] + 1, end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/C_state3.chg')
            total += 1

    # MA state_1
    filename = 'MA'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
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
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MA_state1.chg')
            total += 1

    # MA state_2
    filename = 'RMA3C'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + 1,
              '-',
              Numatom[1] + Numatom[2],
              sep='',
              end=' ',
              file=f_text)
        print(charge[2], end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MA_state2.chg')
            total += 1

    # MA state_3
    filename = 'RMA3C'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + Numatom[2] + 1,
              '-',
              Numatom[1] + Numatom[2] + Numatom[2],
              sep='',
              end=' ',
              file=f_text)
        print(charge[2], end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MA_state3.chg')
            total += 1

    # MA state_4
    filename = 'RMA3C'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + Numatom[2] + Numatom[2] + 1,
              '-',
              Numatom[1] + Numatom[2] + Numatom[2] + Numatom[2],
              sep='',
              end=' ',
              file=f_text)
        print(charge[2], end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MA_state4.chg')
            total += 1

    # MA state_5
    filename = 'RMA2MBC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + Numatom[2] + 1,
              '-',
              Numatom[1] + Numatom[2] + Numatom[2],
              sep='',
              end=' ',
              file=f_text)
        print(charge[2], end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MA_state5.chg')
            total += 1

    # MB state_1
    filename = 'MB'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
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
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MB_state1.chg')
            total += 1

    # MB state_2
    filename = 'RMAMB2C'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + Numatom[2] + 1,
              '-',
              Numatom[1] + Numatom[2] + Numatom[3],
              sep='',
              end=' ',
              file=f_text)
        print(charge[3], end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MB_state2.chg')
            total += 1

    # MB state_3
    filename = 'RMB3C'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + Numatom[3] + 1,
              '-',
              Numatom[1] + Numatom[3] + Numatom[3],
              sep='',
              end=' ',
              file=f_text)
        print(charge[3], end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MB_state3.chg')
            total += 1

    # MB state_4
    filename = 'RMA2MBC'
    if os.access('./gaussian/' + filename + 'step2.fchk', os.F_OK):
        f_out = open('./gaussian/' + filename + '.sh', 'w')
        print('Multiwfn  <<EOF', file=f_out)
        print('./', filename, 'step2.fchk', sep='', file=f_out)
        print('7', file=f_out)
        print('18', file=f_out)
        print('6', file=f_out)
        print('1', file=f_out)
        print(filename, '.txt', sep='', file=f_out)
        print('2', file=f_out)
        print('y', file=f_out)
        print('q', file=f_out)
        print('EOF', end='', file=f_out)

        f_text = open('./gaussian/' + filename + '.txt', 'w')
        print(Numatom[1] + Numatom[2] + Numatom[2] + 1,
              '-',
              Numatom[1] + Numatom[2] + Numatom[2] + Numatom[3],
              sep='',
              end=' ',
              file=f_text)
        print(charge[3], end='', file=f_text)

        f_out.close()
        f_text.close()
        step4 = subprocess.Popen('sh ' + filename + '.sh',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True,
                                 cwd='./gaussian')
        step4.wait()
        if not os.access('./gaussian/' + filename + 'step2.chg', os.F_OK):
            print('Charge caculating is wrong (' + filename + ')')
            cout += 1
        else:
            os.rename('./gaussian/' + filename + 'step2.chg',
                      './gaussian/MB_state4.chg')
            total += 1

    if cout != 0:
        print('Charge caculating is wrong!')
        return False
    if total != 16:
        print('Molecular states are missing!')
        return None
    return True
