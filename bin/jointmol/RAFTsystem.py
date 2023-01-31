import json
import os
import glob
import subprocess
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def _solvemol():
    Mol = []
    model = ['I.mol', 'R.mol', 'MA.mol', 'MB.mol', 'C.mol']
    for i in range(len(model)):
        mol = Chem.MolFromMolFile('./model/' + model[i])
        mol = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        Mol.append(mol)
    return Mol


def _rectpoint():
    reactlist = []

    Mol = _solvemol()
    # I.mol
    list1 = [
        x.GetIdx() for x in Mol[0].GetAtoms()
        if x.GetSymbol() == 'C' and x.GetFormalCharge() == -1
    ]
    # R.mol
    list2 = [
        x.GetIdx() for x in Mol[1].GetAtoms()
        if x.GetSymbol() == 'C' and x.GetFormalCharge() == -1
    ]
    # MA.mol
    list3 = [0] * 2
    patt = Chem.MolFromSmiles('C=C')
    patt_lists = Mol[2].GetSubstructMatches(patt)
    AllChem.ComputeGasteigerCharges(Mol[2])
    ath = int(patt_lists[0][0])
    bth = int(patt_lists[0][1])
    atom1 = Mol[2].GetAtomWithIdx(ath)
    atom2 = Mol[2].GetAtomWithIdx(bth)
    charge1 = float(atom1.GetProp('_GasteigerCharge'))
    charge2 = float(atom2.GetProp('_GasteigerCharge'))
    if charge1 <= charge2:
        list3[0] = ath
        list3[1] = bth
    else:
        list3[0] = bth
        list3[1] = ath
    # MB.mol
    list4 = [0] * 2
    patt = Chem.MolFromSmiles('C=C')
    patt_lists = Mol[3].GetSubstructMatches(patt)
    AllChem.ComputeGasteigerCharges(Mol[3])
    ath = int(patt_lists[0][0])
    bth = int(patt_lists[0][1])
    atom1 = Mol[3].GetAtomWithIdx(ath)
    atom2 = Mol[3].GetAtomWithIdx(bth)
    charge1 = float(atom1.GetProp('_GasteigerCharge'))
    charge2 = float(atom2.GetProp('_GasteigerCharge'))
    if charge1 <= charge2:
        list4[0] = ath
        list4[1] = bth
    else:
        list4[0] = bth
        list4[1] = ath
    # C.mol
    list5 = [0] * 2
    temp = [
        x.GetIdx() for x in Mol[4].GetAtoms()
        if x.GetSymbol() == 'S' and x.GetFormalCharge() == -1
    ]
    list5[0] = temp[0]
    temp = [
        x.GetIdx() for x in Mol[4].GetAtoms()
        if x.GetSymbol() == 'S' and x.GetFormalCharge() == 0
    ]
    list5[1] = temp[0]
    reactlist.append(list1)
    reactlist.append(list2)
    reactlist.append(list3)
    reactlist.append(list4)
    reactlist.append(list5)
    return reactlist


def _Numatom():
    Numatom = []
    Mols = _solvemol()
    for mol in Mols:
        temp = mol.GetNumAtoms()
        Numatom.append(temp)
    return Numatom


def _Npolymer():
    f_sys = open('./control/system.json', )  # boxsize & PISA
    system = json.load(f_sys)
    boxsize = []
    polymer = []
    nCTA = 0
    nINI = 0
    nMA = 0
    nMB = 0

    if 'Boxsize' in system:
        boxsize = system['Boxsize']
    else:
        print('Boxsize(nm) is missing')
        print('Please add "Boxsize": [X(nm), Y(nm), Z(nm)]')
        return nCTA, nINI, nMA, nMB

    if 'PISA' in system:
        polymer = system['PISA']
    else:
        print('PISA system(mol/L) is missing')
        print(
            'Please add: "PISA": [{"CTA": C_CTA(Mol/L), "INI": C_INI(Mol/L), "MA": C_MA(Mol/L), "MB": C_MB(Mol/L)}]'
        )
        return nCTA, nINI, nMA, nMB

    if boxsize != [] and polymer != []:
        x, y, z = boxsize[0], boxsize[1], boxsize[2]
        mCTA, mINI, mMA, mMB = polymer[0]['CTA'], polymer[0]['INI'], polymer[
            0]['MA'], polymer[0]['MB']

        nINI = int(round(0.6023 * x * y * z * mINI + 0.5))
        if nINI == 0:
            nINI = 1
        nCTA = int(round(nINI * (mCTA / mINI)))
        if nCTA == 0:
            nCTA = 1
        nMA = int(round(nCTA * (mMA / mCTA)))
        if nMA == 0:
            nMA = 1
        nMB = int(round(nCTA * (mMB / mCTA)))
        if nMB == 0:
            nMB = 1

    f_sys.close()
    return nCTA, nINI, nMA, nMB


def _getcharge(ele, state):

    f_in = open('./topper/' + ele + '.prm', 'r')
    pos = ele + state
    while True:
        line = f_in.readline()
        if line == ';' + pos + '\n':
            break
    while True:
        line = f_in.readline()
        if line == '[ charges ]\n':
            line = f_in.readline()
            break
    charge = line.split()

    f_in.close()
    return charge


def _getstr(filename):

    Name = filename
    if filename == 'R':
        Name = 'RC'
    if filename == 'I':
        Name = 'IC'
    if filename == 'C':
        Name = 'RC'
    f_in = open('./gaussian/' + Name + '.itp', 'r')

    atom = []
    bond = []
    angle = []
    pair = []
    proper = []
    improper = []
    # atom
    while True:
        line = f_in.readline()
        if line == '[ atoms ]\n':
            line = f_in.readline()
            break
        if line == '':
            break
    while True:
        if line == '':
            break
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            lines = line.split()
            atom.append(lines)
    # bond
    while True:
        line = f_in.readline()
        if line == '[ bonds ]\n':
            line = f_in.readline()
            break
        if line == '':
            break
    while True:
        if line == '':
            break
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            lines = line.split()
            bond.append(lines)
    # angle
    while True:
        line = f_in.readline()
        if line == '[ angles ]\n':
            line = f_in.readline()
            break
        if line == '':
            break
    while True:
        if line == '':
            break
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            lines = line.split()
            angle.append(lines)
    # proper
    while True:
        line = f_in.readline()
        if line == '[ dihedrals ] ; propers\n':
            line = f_in.readline()
            line = f_in.readline()
            break
        if line == '':
            break
    while True:
        if line == '':
            break
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            lines = line.split()
            proper.append(lines)
    # pair
    while True:
        line = f_in.readline()
        if line == '[ pairs ] ; Yielded based on rotatable dihedrals\n':
            line = f_in.readline()
            break
        if line == '':
            break
    while True:
        if line == '':
            break
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            lines = line.split()
            pair.append(lines)
    # improper
    while True:
        line = f_in.readline()
        if line == '[ dihedrals ] ; impropers\n':
            line = f_in.readline()
            break
        if line == '':
            break
    while True:
        if line == '':
            break
        line = f_in.readline()
        if line == '':
            break
        else:
            lines = line.split()
            improper.append(lines)

    if filename == 'I':
        Numatom = _Numatom()
        sort_atom = []
        sort_bond = []
        sort_angle = []
        sort_pair = []
        sort_proper = []
        sort_improper = []
        if atom != []:
            for i in range(0, Numatom[0]):
                sort_atom.append(atom[i])
        if bond != []:
            for i in range(len(bond)):
                line = bond[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[0]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[0]):
                    sort_bond.append(bond[i])
        if angle != []:
            for i in range(len(angle)):
                line = angle[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[0]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[0]) and (
                            int(line[2]) >= 1 and int(line[2]) <= Numatom[0]):
                    sort_angle.append(angle[i])
        if pair != []:
            for i in range(len(pair)):
                line = pair[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[0]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[0]):
                    sort_pair.append(pair[i])
        if proper != []:
            for i in range(len(proper)):
                line = proper[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[0]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[0]
                ) and (int(line[2]) >= 1 and int(line[2]) <= Numatom[0]) and (
                        int(line[3]) >= 1 and int(line[3]) <= Numatom[0]):
                    sort_proper.append(proper[i])
        if improper != []:
            for i in range(len(improper)):
                line = improper[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[0]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[0]
                ) and (int(line[2]) >= 1 and int(line[2]) <= Numatom[0]) and (
                        int(line[3]) >= 1 and int(line[3]) <= Numatom[0]):
                    sort_improper.append(improper[i])

        f_in.close()
        return sort_atom, sort_bond, sort_angle, sort_proper, sort_pair, sort_improper

    if filename == 'R':
        Numatom = _Numatom()
        sort_atom = []
        sort_bond = []
        sort_angle = []
        sort_pair = []
        sort_proper = []
        sort_improper = []
        if atom != []:
            for i in range(0, Numatom[1]):
                sort_atom.append(atom[i])
        if bond != []:
            for i in range(len(bond)):
                line = bond[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[1]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[1]):
                    sort_bond.append(bond[i])
        if angle != []:
            for i in range(len(angle)):
                line = angle[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[1]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[1]) and (
                            int(line[2]) >= 1 and int(line[2]) <= Numatom[1]):
                    sort_angle.append(angle[i])
        if pair != []:
            for i in range(len(pair)):
                line = pair[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[1]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[1]):
                    sort_pair.append(pair[i])
        if proper != []:
            for i in range(len(proper)):
                line = proper[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[1]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[1]
                ) and (int(line[2]) >= 1 and int(line[2]) <= Numatom[1]) and (
                        int(line[3]) >= 1 and int(line[3]) <= Numatom[1]):
                    sort_proper.append(proper[i])
        if improper != []:
            for i in range(len(improper)):
                line = improper[i]
                if (int(line[0]) >= 1 and int(line[0]) <= Numatom[1]) and (
                        int(line[1]) >= 1 and int(line[1]) <= Numatom[1]
                ) and (int(line[2]) >= 1 and int(line[2]) <= Numatom[1]) and (
                        int(line[3]) >= 1 and int(line[3]) <= Numatom[1]):
                    sort_improper.append(improper[i])

        f_in.close()
        return sort_atom, sort_bond, sort_angle, sort_proper, sort_pair, sort_improper

    if filename == 'C':
        Numatom = _Numatom()
        sort_atom = []
        sort_bond = []
        sort_angle = []
        sort_pair = []
        sort_proper = []
        sort_improper = []
        if atom != []:
            for i in range(Numatom[1], Numatom[1] + Numatom[4]):
                line = atom[i]
                line[0] = str(int(line[0]) - Numatom[1])
                line[5] = str(int(line[5]) - Numatom[1])
                sort_atom.append(line)
        if bond != []:
            for i in range(len(bond)):
                line = bond[i]
                if (int(line[0]) >= 1 + Numatom[1] and int(line[0]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[1]) >= 1 + Numatom[1] and int(line[1]) <= Numatom[1] + Numatom[4]):
                    line[0] = str(int(line[0]) - Numatom[1])
                    line[1] = str(int(line[1]) - Numatom[1])
                    sort_bond.append(line)
        if angle != []:
            for i in range(len(angle)):
                line = angle[i]
                if (int(line[0]) >= 1 + Numatom[1] and int(line[0]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[1]) >= 1 + Numatom[1] and int(line[1]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[2]) >= 1 + Numatom[1] and int(line[2]) <= Numatom[1] + Numatom[4]):
                    line[0] = str(int(line[0]) - Numatom[1])
                    line[1] = str(int(line[1]) - Numatom[1])
                    line[2] = str(int(line[2]) - Numatom[1])
                    sort_angle.append(line)
        if pair != []:
            for i in range(len(pair)):
                line = pair[i]
                if (int(line[0]) >= 1 + Numatom[1] and int(line[0]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[1]) >= 1 + Numatom[1] and int(line[1]) <= Numatom[1] + Numatom[4]):
                    line[0] = str(int(line[0]) - Numatom[1])
                    line[1] = str(int(line[1]) - Numatom[1])
                    sort_pair.append(line)
        if proper != []:
            for i in range(len(proper)):
                line = proper[i]
                if (int(line[0]) >= 1 + Numatom[1] and int(line[0]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[1]) >= 1 + Numatom[1] and int(line[1]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[2]) >= 1 + Numatom[1] and int(line[2]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[3]) >= 1 + Numatom[1] and int(line[3]) <= Numatom[1] + Numatom[4]):
                    line[0] = str(int(line[0]) - Numatom[1])
                    line[1] = str(int(line[1]) - Numatom[1])
                    line[2] = str(int(line[2]) - Numatom[1])
                    line[3] = str(int(line[3]) - Numatom[1])
                    sort_proper.append(line)
        if improper != []:
            for i in range(len(improper)):
                line = improper[i]
                if (int(line[0]) >= 1 + Numatom[1] and int(line[0]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[1]) >= 1 + Numatom[1] and int(line[1]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[2]) >= 1 + Numatom[1] and int(line[2]) <= Numatom[1] + Numatom[4]) and \
                   (int(line[3]) >= 1 + Numatom[1] and int(line[3]) <= Numatom[1] + Numatom[4]):
                    line[0] = str(int(line[0]) - Numatom[1])
                    line[1] = str(int(line[1]) - Numatom[1])
                    line[2] = str(int(line[2]) - Numatom[1])
                    line[3] = str(int(line[3]) - Numatom[1])
                    sort_improper.append(line)

        f_in.close()
        return sort_atom, sort_bond, sort_angle, sort_proper, sort_pair, sort_improper

    f_in.close()
    return atom, bond, angle, proper, pair, improper


def _getpos():

    f1 = open('./gaussian/I_state2.chg', 'r')
    f2 = open('./gaussian/R_state2.chg', 'r')
    f3 = open('./gaussian/C_state1.chg', 'r')
    f4 = open('./gaussian/MA_state1.chg', 'r')
    f5 = open('./gaussian/MB_state1.chg', 'r')
    Numatom = _Numatom()
    Ipos = []
    Rpos = []
    Cpos = []
    MApos = []
    MBpos = []
    cout = 0
    while True:
        line = f1.readline()
        if cout == Numatom[0]:
            break
        Ipos.append(line)
        cout += 1
    cout = 0
    while True:
        line = f2.readline()
        if cout == Numatom[1]:
            break
        Rpos.append(line)
        cout += 1
    cout = 0
    while True:
        line = f3.readline()
        if line == '':
            break
        if cout >= Numatom[1]:
            Cpos.append(line)
        if cout == Numatom[1] + Numatom[4]:
            break
        cout += 1
    while True:
        line = f4.readline()
        if line == '':
            break
        MApos.append(line)
    while True:
        line = f5.readline()
        if line == '':
            break
        MBpos.append(line)

    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    return Ipos, Rpos, Cpos, MApos, MBpos


def _getconnectmatrix(Natom, bondlist):

    matrix = np.zeros((Natom, Natom), dtype=np.uint8)
    for i in range(len(bondlist)):
        lines = bondlist[i]
        matrix[int(lines[0]) - 1][int(lines[1]) - 1] = 1
        matrix[int(lines[1]) - 1][int(lines[0]) - 1] = 1

    return matrix


def _CTAspecialffbond(ffbondtype):

    ffbond = []
    for modelname in ['I', 'R', 'C']:
        f_prm = open('./topper/' + modelname + '.prm', 'r')
        for i in range(1, 3):
            while True:
                line = f_prm.readline()
                if line == ';' + modelname + str(i) + '\n':
                    break
            while True:
                line = f_prm.readline()
                if line == '[ ' + ffbondtype + ' ]\n':
                    break
            while True:
                line = f_prm.readline()
                if line == '\n' or line == ' \n' or line == '':
                    break
                else:
                    lines = line.split()
                    ffbond.append(lines)
        f_prm.close()

    f_prm = open('./topper/MA.prm', 'r')
    for i in range(1, 5):
        while True:
            line = f_prm.readline()
            if line == ';MA' + str(i) + '\n':
                break
        while True:
            line = f_prm.readline()
            if line == '[ ' + ffbondtype + ' ]\n':
                break
        while True:
            line = f_prm.readline()
            if line == '\n' or line == ' \n' or line == '':
                break
            else:
                lines = line.split()
                ffbond.append(lines)
    f_prm.close()

    return ffbond


def _CTAchangeMol(atom, Molrange, modelname, state):

    if modelname != 'C':
        f_prm = open('./topper/' + modelname + '.prm', 'r')
        type = []
        charge = []

        while True:
            line = f_prm.readline()
            if line == ';' + modelname + str(state) + '\n':
                break
        while True:
            line = f_prm.readline()
            if line == '[ atoms ]\n':
                line = f_prm.readline()
                type = line.split()
                break
        while True:
            line = f_prm.readline()
            if line == '[ charges ]\n':
                line = f_prm.readline()
                charge = line.split()
                break
        if type != [] and charge != []:
            for i in range(Molrange[0] - 1, Molrange[1]):
                atom[i][1] = type[i - Molrange[0] + 1]
                atom[i][6] = charge[i - Molrange[0] + 1]
    else:
        f_prm = open('./topper/C.prm', 'r')
        type = []
        charge = []
        s_charge = None
        ss_charge = None
        s_id = None
        ss_id = None

        while True:
            line = f_prm.readline()
            if line == ';' + 'C' + str(state) + '\n':
                break
        while True:
            line = f_prm.readline()
            if line == '[ atoms ]\n':
                line = f_prm.readline()
                type = line.split()
                break
        while True:
            line = f_prm.readline()
            if line == '[ charges ]\n':
                line = f_prm.readline()
                charge = line.split()
                break
        if type != [] and charge != []:
            for i in range(len(type)):
                if type[i] == 's':
                    s_charge = charge[i]
                if type[i] == 'ss':
                    ss_charge = charge[i]
                if s_charge is not None and ss_charge is not None:
                    break

            for i in range(Molrange[0] - 1, Molrange[1]):
                if atom[i][1] != 's' and atom[i][1] != 'ss':
                    atom[i][1] = type[i - Molrange[0] + 1]
                    atom[i][6] = charge[i - Molrange[0] + 1]
                if atom[i][1] == 's':
                    atom[i][6] = s_charge
                    s_id = i
                if atom[i][1] == 'ss':
                    atom[i][6] = ss_charge
                    atom[i][4] = 'c_b'
                    ss_id = i
            for i in range(len(atom[0])):
                if i != 0 and i != 5:
                    a = atom[s_id][i]
                    atom[s_id][i] = atom[ss_id][i]
                    atom[ss_id][i] = a
    f_prm.close()


def _CTAchangeffbond(Natom, adjacent, reac, atomtype):
    bond = []
    angle = []
    dihedral = []
    improper = []
    fin_itp = open('./gromacs/step3/PISA.itp', 'r')

    while True:
        line = fin_itp.readline()
        if line == '[ bonds ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            bond.append(lines)

    while True:
        line = fin_itp.readline()
        if line == '[ angles ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            angle.append(lines)

    while True:
        line = fin_itp.readline()
        if line == '[ dihedrals ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            dihedral.append(lines)

    while True:
        line = fin_itp.readline()
        if line == '[ dihedrals ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            improper.append(lines)
    fin_itp.close()

    ffbond = _CTAspecialffbond('bonds')
    ffangle = _CTAspecialffbond('angles')
    ffdihedral = _CTAspecialffbond('dihedrals')
    # bond
    for i in range(len(reac)):
        ffbond_M = ffbond
        special = False
        if ffbond != []:
            if int(reac[i][0]) < int(reac[i][1]):
                a = atomtype[reac[i][0]]
                b = atomtype[reac[i][1]]
                for j in range(len(ffbond)):
                    if ffbond_M[j][0] == a and ffbond_M[j][1] == b:
                        temp = ffbond_M[j]
                        temp[j][0] = reac[i][0]
                        temp[j][1] = reac[i][1]
                        bond.append(temp)
                        if (j + 1) < len(ffbond):
                            if ffbond_M[j][0] == ffbond_M[
                                    j +
                                    1][0] and ffbond_M[j][1] == ffbond_M[j +
                                                                         1][1]:
                                ffbond_M.pop(j)
                        special = True
                        break
                if not special:
                    temp = [''] * 3
                    temp[0] = reac[i][0]
                    temp[1] = reac[i][1]
                    temp[2] = '1'
                    bond.append(temp)
            if int(reac[i][0]) > int(reac[i][1]):
                a = atomtype[reac[i][1]]
                b = atomtype[reac[i][0]]
                for j in range(len(ffbond)):
                    if ffbond_M[j][0] == a and ffbond_M[j][1] == b:
                        temp = ffbond_M[j]
                        temp[j][0] = reac[i][1]
                        temp[j][1] = reac[i][0]
                        bond.append(temp)
                        if (j + 1) < len(ffbond):
                            if ffbond_M[j][0] == ffbond_M[
                                    j +
                                    1][0] and ffbond_M[j][1] == ffbond_M[j +
                                                                         1][1]:
                                ffbond_M.pop(j)
                        special = True
                        break
                    if not special:
                        temp = [''] * 3
                        temp[0] = reac[i][1]
                        temp[1] = reac[i][0]
                        temp[2] = '1'
                        bond.append(temp)
            ffbond_M = ffbond
        else:
            if int(reac[i][0]) < int(reac[i][1]):
                temp = [''] * 3
                temp[0] = reac[i][0]
                temp[1] = reac[i][1]
                temp[2] = '1'
                bond.append(temp)
            if int(reac[i][0]) > int(reac[i][1]):
                temp = [''] * 3
                temp[0] = reac[i][1]
                temp[1] = reac[i][0]
                temp[2] = '1'
                bond.append(temp)
    if ffbond != []:
        for i in range(len(bond)):
            special = True
            for j in range(len(ffbond)):
                if atomtype[bond[i][0]] == ffbond[j][0] and atomtype[
                        bond[i][1]] == ffbond[j][1]:
                    special = False
                    break
            if special and len(bond[i]) > 3:
                temp = [''] * 3
                temp[0] = bond[i][0]
                temp[1] = bond[i][1]
                temp[2] = '1'
                bond[i] = temp
    bond.sort(key=lambda x: int(x[0]))

    # angle
    for i in range(len(reac)):
        special = False
        ffangle_L = ffangle
        ffangle_R = ffangle
        if ffangle != []:
            for j in range(Natom):
                if adjacent.has_edge(int(reac[i][0]) - 1, j):
                    if (j + 1) < int(reac[i][1]):
                        a = atomtype[str(j + 1)]
                        b = atomtype[reac[i][0]]
                        c = atomtype[reac[i][1]]
                        for k in range(len(ffangle)):
                            if ffangle_L[k][0] == a and ffangle_L[k][
                                    1] == b and ffangle_L[k][2] == c:
                                temp = ffangle_L[k]
                                temp[0] = str(j + 1)
                                temp[1] = reac[i][0]
                                temp[2] = reac[i][1]
                                angle.append(temp)
                                if (k + 1) < len(ffangle):
                                    if ffangle_L[k][0] == ffangle_L[
                                            k + 1][0] and ffangle_L[k][
                                                1] == ffangle_L[
                                                    k + 1][1] and ffangle_L[k][
                                                        2] == ffangle_L[k +
                                                                        1][2]:
                                        ffangle_L.pop(k)
                                special = True
                                break
                        if not special:
                            temp = [''] * 4
                            temp[0] = str(j + 1)
                            temp[1] = reac[i][0]
                            temp[2] = reac[i][1]
                            temp[3] = '1'
                            angle.append(temp)
                    if (j + 1) > int(reac[i][1]):
                        a = atomtype[reac[i][1]]
                        b = atomtype[reac[i][0]]
                        c = atomtype[str(j)]
                        for k in range(len(ffangle)):
                            if ffangle_L[k][0] == a and ffangle_L[k][
                                    1] == b and ffangle_L[k][2] == c:
                                temp = ffangle_L[k]
                                temp[0] = reac[i][1]
                                temp[1] = reac[i][0]
                                temp[2] = str(j + 1)
                                angle.append(temp)
                                if (k + 1) < len(ffangle):
                                    if ffangle_L[k][0] == ffangle_L[
                                            k + 1][0] and ffangle_L[k][
                                                1] == ffangle_L[
                                                    k + 1][1] and ffangle_L[k][
                                                        2] == ffangle_L[k +
                                                                        1][2]:
                                        ffangle_L.pop(k)
                                special = True
                                break
                        if not special:
                            temp = [''] * 4
                            temp[0] = reac[i][1]
                            temp[1] = reac[i][0]
                            temp[2] = str(j + 1)
                            temp[3] = '1'
                            angle.append(temp)
                if adjacent.has_edge(int(reac[i][1]) - 1, j):
                    if (j + 1) < int(reac[i][0]):
                        a = atomtype[str(j + 1)]
                        b = atomtype[reac[i][1]]
                        c = atomtype[reac[i][0]]
                        for k in range(len(ffangle)):
                            if ffangle_R[k][0] == a and ffangle_R[k][
                                    1] == b and ffangle_R[k][2] == c:
                                temp = ffangle_R[k]
                                temp[0] = str(j + 1)
                                temp[1] = reac[i][1]
                                temp[2] = reac[i][0]
                                angle.append(temp)
                                if (k + 1) < len(ffangle):
                                    if ffangle_R[k][0] == ffangle_R[
                                            k + 1][0] and ffangle_R[k][
                                                1] == ffangle_R[
                                                    k + 1][1] and ffangle_R[k][
                                                        2] == ffangle_R[k +
                                                                        1][2]:
                                        ffangle_R.pop(k)
                                special = True
                                break
                        if not special:
                            temp = [''] * 4
                            temp[0] = str(j + 1)
                            temp[1] = reac[i][1]
                            temp[2] = reac[i][0]
                            temp[3] = '1'
                            angle.append(temp)
                    if (j + 1) > int(reac[i][0]):
                        a = atomtype[reac[i][1]]
                        b = atomtype[reac[i][0]]
                        c = atomtype[str(j)]
                        for k in range(len(ffangle)):
                            if ffangle_R[k][0] == a and ffangle_R[k][
                                    1] == b and ffangle_R[k][2] == c:
                                temp = ffangle_R[k]
                                temp[0] = reac[i][0]
                                temp[1] = reac[i][1]
                                temp[2] = str(j + 1)
                                angle.append(temp)
                                if (k + 1) < len(ffangle):
                                    if ffangle_R[k][0] == ffangle_R[
                                            k + 1][0] and ffangle_R[k][
                                                1] == ffangle_R[
                                                    k + 1][1] and ffangle_R[k][
                                                        2] == ffangle_R[k +
                                                                        1][2]:
                                        ffangle_R.pop(k)
                                special = True
                                break
                        if not special:
                            temp = [''] * 4
                            temp[0] = reac[i][0]
                            temp[1] = reac[i][1]
                            temp[2] = str(j + 1)
                            temp[3] = '1'
                            angle.append(temp)
            ffangle_R = ffangle
            ffangle_L = ffangle
        else:
            if adjacent.has_edge(int(reac[i][0]) - 1, j):
                if (j + 1) < int(reac[i][1]):
                    temp = [''] * 4
                    temp[0] = str(j + 1)
                    temp[1] = reac[i][0]
                    temp[2] = reac[i][1]
                    temp[3] = '1'
                    angle.append(temp)
                if (j + 1) > int(reac[i][1]):
                    temp = [''] * 4
                    temp[0] = reac[i][1]
                    temp[1] = reac[i][0]
                    temp[2] = str(j + 1)
                    temp[3] = '1'
                    angle.append(temp)
            if adjacent.has_edge(int(reac[i][1]) - 1, j):
                if (j + 1) < int(reac[i][0]):
                    temp = [''] * 4
                    temp[0] = str(j + 1)
                    temp[1] = reac[i][1]
                    temp[2] = reac[i][0]
                    temp[3] = '1'
                    angle.append(temp)
                if (j + 1) > int(reac[i][0]):
                    temp = [''] * 4
                    temp[0] = reac[i][0]
                    temp[1] = reac[i][1]
                    temp[2] = str(j + 1)
                    temp[3] = '1'
                    angle.append(temp)
    if ffangle != []:
        for i in range(len(angle)):
            special = True
            for j in range(len(ffangle)):
                if atomtype[angle[i][0]] == ffangle[j][0] and atomtype[
                        angle[i][1]] == ffangle[j][1] and atomtype[
                            angle[i][2]] == ffangle[j][2]:
                    special = False
                    break
            if special and len(angle[i]) > 4:
                temp = [''] * 4
                temp[0] = angle[i][0]
                temp[1] = angle[i][1]
                temp[2] = angle[i][2]
                temp[3] = '1'
                angle[i] = temp
    angle.sort(key=lambda x: int(x[0]))

    # dihedral
    for i in range(len(reac)):
        special = False
        ffdihedral_M = ffdihedral
        ffdihedral_L = ffdihedral
        ffdihedral_R = ffdihedral
        if ffdihedral != []:
            for j in range(Natom):
                if adjacent.has_edge(int(reac[i][0]) - 1,
                                     j) and (j + 1) != int(reac[i][1]):
                    for k in range(Natom):
                        if adjacent.has_edge(int(reac[i][1]) - 1,
                                             k) and (k + 1) != int(reac[i][0]):
                            if j < k:
                                a = atomtype[str(j + 1)]
                                b = atomtype[reac[i][0]]
                                c = atomtype[reac[i][1]]
                                d = atomtype[str(k + 1)]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_M[ff][0] == a and \
                                       ffdihedral_M[ff][1] == b and \
                                       ffdihedral_M[ff][2] == c and \
                                       ffdihedral_M[ff][3] == d:
                                        temp = ffdihedral_M[ff]
                                        temp[0] = str(j + 1)
                                        temp[1] = reac[i][0]
                                        temp[2] = reac[i][1]
                                        temp[3] = str(k + 1)
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_M[ff][0] == ffdihedral_M[ff+1][0] and \
                                               ffdihedral_M[ff][1] == ffdihedral_M[ff+1][1] and \
                                               ffdihedral_M[ff][2] == ffdihedral_M[ff+1][2] and \
                                               ffdihedral_M[ff][3] == ffdihedral_M[ff+1][3]:
                                                ffdihedral_M.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = str(j + 1)
                                    temp[1] = reac[i][0]
                                    temp[2] = reac[i][1]
                                    temp[3] = str(k + 1)
                                    temp[4] = '9'
                                    dihedral.append(temp)
                            if j > k:
                                a = atomtype[str(k + 1)]
                                b = atomtype[reac[i][1]]
                                c = atomtype[reac[i][0]]
                                d = atomtype[str(j + 1)]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_M[ff][0] == a and \
                                       ffdihedral_M[ff][1] == b and \
                                       ffdihedral_M[ff][2] == c and \
                                       ffdihedral_M[ff][3] == d:
                                        temp = ffdihedral_M[ff]
                                        temp[0] = str(k + 1)
                                        temp[1] = reac[i][1]
                                        temp[2] = reac[i][0]
                                        temp[3] = str(j + 1)
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_M[ff][0] == ffdihedral_M[ff+1][0] and \
                                               ffdihedral_M[ff][1] == ffdihedral_M[ff+1][1] and \
                                               ffdihedral_M[ff][2] == ffdihedral_M[ff+1][2] and \
                                               ffdihedral_M[ff][3] == ffdihedral_M[ff+1][3]:
                                                ffdihedral_M.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = str(k + 1)
                                    temp[1] = reac[i][1]
                                    temp[2] = reac[i][0]
                                    temp[3] = str(j + 1)
                                    temp[4] = '9'
                                    dihedral.append(temp)
                    for k in range(Natom):
                        if adjacent.has_edge(j,
                                             k) and (k + 1) != int(reac[i][0]):
                            if (k + 1) < int(reac[i][1]):
                                a = atomtype[str(k + 1)]
                                b = atomtype[str(j + 1)]
                                c = atomtype[reac[i][0]]
                                d = atomtype[reac[i][1]]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_L[ff][0] == a and \
                                       ffdihedral_L[ff][1] == b and \
                                       ffdihedral_L[ff][2] == c and \
                                       ffdihedral_L[ff][3] == d:
                                        temp = ffdihedral_L[ff]
                                        temp[0] = str(k + 1)
                                        temp[1] = str(j + 1)
                                        temp[2] = reac[i][0]
                                        temp[3] = reac[i][1]
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_L[ff][0] == ffdihedral_L[ff+1][0] and \
                                               ffdihedral_L[ff][1] == ffdihedral_L[ff+1][1] and \
                                               ffdihedral_L[ff][2] == ffdihedral_L[ff+1][2] and \
                                               ffdihedral_L[ff][3] == ffdihedral_L[ff+1][3]:
                                                ffdihedral_L.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = str(k + 1)
                                    temp[1] = str(j + 1)
                                    temp[2] = reac[i][0]
                                    temp[3] = reac[i][1]
                                    temp[4] = '9'
                                    dihedral.append(temp)
                            if (k + 1) > int(reac[i][1]):
                                a = atomtype[reac[i][1]]
                                b = atomtype[reac[i][0]]
                                c = atomtype[str(j + 1)]
                                d = atomtype[str(k + 1)]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_L[ff][0] == a and \
                                       ffdihedral_L[ff][1] == b and \
                                       ffdihedral_L[ff][2] == c and \
                                       ffdihedral_L[ff][3] == d:
                                        temp = ffdihedral_L[ff]
                                        temp[0] = reac[i][1]
                                        temp[1] = reac[i][0]
                                        temp[2] = str(j + 1)
                                        temp[3] = str(k + 1)
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_L[ff][0] == ffdihedral_L[ff+1][0] and \
                                               ffdihedral_L[ff][1] == ffdihedral_L[ff+1][1] and \
                                               ffdihedral_L[ff][2] == ffdihedral_L[ff+1][2] and \
                                               ffdihedral_L[ff][3] == ffdihedral_L[ff+1][3]:
                                                ffdihedral_L.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = reac[i][1]
                                    temp[1] = reac[i][0]
                                    temp[2] = str(j + 1)
                                    temp[3] = str(k + 1)
                                    temp[4] = '9'
                                    dihedral.append(temp)
                if adjacent.has_edge(int(reac[i][1]) - 1,
                                     j) and (j + 1) != int(reac[i][0]):
                    for k in range(Natom):
                        if adjacent.has_edge(j,
                                             k) and (k + 1) != int(reac[i][1]):
                            if (k + 1) < int(reac[i][0]):
                                a = atomtype[str(k + 1)]
                                b = atomtype[str(j + 1)]
                                c = atomtype[reac[i][1]]
                                d = atomtype[reac[i][0]]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_R[ff][0] == a and \
                                       ffdihedral_R[ff][1] == b and \
                                       ffdihedral_R[ff][2] == c and \
                                       ffdihedral_R[ff][3] == d:
                                        temp = ffdihedral_R[ff]
                                        temp[0] = str(k + 1)
                                        temp[1] = str(j + 1)
                                        temp[2] = reac[i][1]
                                        temp[3] = reac[i][0]
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_R[ff][0] == ffdihedral_R[ff+1][0] and \
                                               ffdihedral_R[ff][1] == ffdihedral_R[ff+1][1] and \
                                               ffdihedral_R[ff][2] == ffdihedral_R[ff+1][2] and \
                                               ffdihedral_R[ff][3] == ffdihedral_R[ff+1][3]:
                                                ffdihedral_R.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = str(k + 1)
                                    temp[1] = str(j + 1)
                                    temp[2] = reac[i][1]
                                    temp[3] = reac[i][0]
                                    temp[4] = '9'
                                    dihedral.append(temp)
                            if (k + 1) > int(reac[i][0]):
                                a = atomtype[reac[i][0]]
                                b = atomtype[reac[i][1]]
                                c = atomtype[str(j + 1)]
                                d = atomtype[str(k + 1)]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_R[ff][0] == a and \
                                       ffdihedral_R[ff][1] == b and \
                                       ffdihedral_R[ff][2] == c and \
                                       ffdihedral_R[ff][3] == d:
                                        temp = ffdihedral_R[ff]
                                        temp[0] = reac[i][0]
                                        temp[1] = reac[i][1]
                                        temp[2] = str(j + 1)
                                        temp[3] = str(k + 1)
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_R[ff][0] == ffdihedral_R[ff+1][0] and \
                                               ffdihedral_R[ff][1] == ffdihedral_R[ff+1][1] and \
                                               ffdihedral_R[ff][2] == ffdihedral_R[ff+1][2] and \
                                               ffdihedral_R[ff][3] == ffdihedral_R[ff+1][3]:
                                                ffdihedral_R.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = reac[i][0]
                                    temp[1] = reac[i][1]
                                    temp[2] = str(j + 1)
                                    temp[3] = str(k + 1)
                                    temp[4] = '9'
                                    dihedral.append(temp)
            ffdihedral_R = ffdihedral
            ffdihedral_L = ffdihedral
        else:
            if adjacent.has_edge(int(reac[i][0]) - 1,
                                 j) and (j + 1) != int(reac[i][1]):
                for k in range(Natom):
                    if adjacent.has_edge(j, k) and (k + 1) != int(reac[i][0]):
                        if (k + 1) < int(reac[i][1]):
                            temp = [''] * 5
                            temp[0] = str(k + 1)
                            temp[1] = str(j + 1)
                            temp[2] = reac[i][0]
                            temp[3] = reac[i][1]
                            temp[4] = '9'
                            dihedral.append(temp)
                        if (k + 1) > int(reac[i][1]):
                            temp = [''] * 5
                            temp[0] = reac[i][1]
                            temp[1] = reac[i][0]
                            temp[2] = str(j + 1)
                            temp[3] = str(k + 1)
                            temp[4] = '9'
                            dihedral.append(temp)
            if adjacent.has_edge(int(reac[i][1]) - 1,
                                 j) and (j + 1) != int(reac[i][0]):
                for k in range(Natom):
                    if adjacent.has_edge(j, k) and (k + 1) != int(reac[i][1]):
                        if (k + 1) < int(reac[i][0]):
                            temp = [''] * 5
                            temp[0] = str(k + 1)
                            temp[1] = str(j + 1)
                            temp[2] = reac[i][1]
                            temp[3] = reac[i][0]
                            temp[4] = '9'
                            dihedral.append(temp)
                        if (k + 1) > int(reac[i][0]):
                            temp = [''] * 5
                            temp[0] = reac[i][0]
                            temp[1] = reac[i][1]
                            temp[2] = str(j + 1)
                            temp[3] = str(k + 1)
                            temp[4] = '9'
                            dihedral.append(temp)
    if ffdihedral != []:
        for i in range(len(dihedral)):
            special = True
            for j in range(len(ffdihedral)):
                if atomtype[dihedral[i][0]] == ffdihedral[j][0] and atomtype[
                        dihedral[i][1]] == ffdihedral[j][1] and atomtype[
                            dihedral[i][2]] == ffdihedral[j][2] and atomtype[
                                dihedral[i][3]] == ffdihedral[j][3]:
                    special = False
                    break
            if special and len(dihedral[i]) > 5:
                temp = [''] * 5
                temp[0] = dihedral[i][0]
                temp[1] = dihedral[i][1]
                temp[2] = dihedral[i][2]
                temp[3] = dihedral[i][3]
                temp[4] = '9'
                dihedral[i] = temp
    dihedral.sort(key=lambda x: int(x[0]))

    # improper
    ffimproper = []
    fin_itp = open('./topper/ffbonded.itp', 'r')
    while True:
        line = fin_itp.readline()
        if line == '[ dihedraltypes ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '[ dihedraltypes ]\n':
            line = fin_itp.readline()
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            ffimproper.append(lines[:4])
    fin_itp.close()

    for i in range(len(reac)):
        bondatom_L = []
        bondatom_R = []
        for j in range(Natom):
            if adjacent.has_edge(int(reac[i][0]) - 1, j):
                bondatom_L.append(str(j + 1))
            if adjacent.has_edge(int(reac[i][1]) - 1, j):
                bondatom_R.append(str(j + 1))
        if len(bondatom_L) >= 2:
            for j in range(len(bondatom_L) - 1):
                for k in range(j + 1, len(bondatom_L)):
                    unsort = [0] * 3
                    unsort[0] = int(bondatom_L[j])
                    unsort[1] = int(bondatom_L[k])
                    unsort[2] = int(reac[i][1])
                    unsort.sort()
                    temp = [''] * 5
                    temp[0] = str(unsort[0])
                    temp[1] = str(unsort[1])
                    temp[2] = reac[i][0]
                    temp[3] = str(unsort[2])
                    temp[4] = '4'
                    improper.append(temp)
        if len(bondatom_R) >= 2:
            for j in range(len(bondatom_R) - 1):
                for k in range(j + 1, len(bondatom_R)):
                    unsort = [0] * 3
                    unsort[0] = int(bondatom_R[j])
                    unsort[1] = int(bondatom_R[k])
                    unsort[2] = int(reac[i][0])
                    unsort.sort()
                    temp = [''] * 5
                    temp[0] = str(unsort[0])
                    temp[1] = str(unsort[1])
                    temp[2] = reac[i][1]
                    temp[3] = str(unsort[2])
                    temp[4] = '4'
                    improper.append(temp)
    for imp in improper[:]:
        special = False
        for i in range(len(ffimproper)):
            if atomtype[imp[0]] == ffimproper[i][0] and atomtype[
                    imp[1]] == ffimproper[i][1] and atomtype[
                        imp[2]] == ffimproper[i][2] and atomtype[
                            imp[3]] == ffimproper[i][3]:
                special = True
                break
        if not special:
            index = improper.index(imp)
            del improper[index]
    improper.sort(key=lambda x: int(x[2]))
    return bond, angle, dihedral, improper


def _PISAspecialffbond(ffbondtype):

    ffbond = []
    for modelname in ['I', 'R']:
        f_prm = open('./topper/' + modelname + '.prm', 'r')
        for i in range(1, 3):
            while True:
                line = f_prm.readline()
                if line == ';' + modelname + str(i) + '\n':
                    break
            while True:
                line = f_prm.readline()
                if line == '[ ' + ffbondtype + ' ]\n':
                    break
            while True:
                line = f_prm.readline()
                if line == '\n' or line == ' \n' or line == '':
                    break
                else:
                    lines = line.split()
                    ffbond.append(lines)
        f_prm.close()

    f_prm = open('./topper/C.prm', 'r')
    for i in range(1, 4):
        while True:
            line = f_prm.readline()
            if line == ';C' + str(i) + '\n':
                break
        while True:
            line = f_prm.readline()
            if line == '[ ' + ffbondtype + ' ]\n':
                break
        while True:
            line = f_prm.readline()
            if line == '\n' or line == ' \n' or line == '':
                break
            else:
                lines = line.split()
                ffbond.append(lines)
    f_prm.close()

    f_prm = open('./topper/MA.prm', 'r')
    for i in range(1, 6):
        while True:
            line = f_prm.readline()
            if line == ';MA' + str(i) + '\n':
                break
        while True:
            line = f_prm.readline()
            if line == '[ ' + ffbondtype + ' ]\n':
                break
        while True:
            line = f_prm.readline()
            if line == '\n' or line == ' \n' or line == '':
                break
            else:
                lines = line.split()
                ffbond.append(lines)
    f_prm.close()

    f_prm = open('./topper/MB.prm', 'r')
    for i in range(1, 5):
        while True:
            line = f_prm.readline()
            if line == ';MB' + str(i) + '\n':
                break
        while True:
            line = f_prm.readline()
            if line == '[ ' + ffbondtype + ' ]\n':
                break
        while True:
            line = f_prm.readline()
            if line == '\n' or line == ' \n' or line == '':
                break
            else:
                lines = line.split()
                ffbond.append(lines)
    f_prm.close()

    return ffbond


def _PISAchangeMol(atom, Molrange, modelname, state):

    if modelname != 'C':
        f_prm = open('./topper/' + modelname + '.prm', 'r')
        type = []
        charge = []

        while True:
            line = f_prm.readline()
            if line == ';' + modelname + str(state) + '\n':
                break
        while True:
            line = f_prm.readline()
            if line == '[ atoms ]\n':
                line = f_prm.readline()
                type = line.split()
                break
        while True:
            line = f_prm.readline()
            if line == '[ charges ]\n':
                line = f_prm.readline()
                charge = line.split()
                break
        if type != [] and charge != []:
            for i in range(Molrange[0] - 1, Molrange[1]):
                atom[i][1] = type[i - Molrange[0] + 1]
                atom[i][6] = charge[i - Molrange[0] + 1]
    else:
        f_prm = open('./topper/C.prm', 'r')
        type = []
        charge = []
        s_charge = None
        ss_charge = None
        s_id = None
        ss_id = None

        while True:
            line = f_prm.readline()
            if line == ';' + 'C' + str(state) + '\n':
                break
        while True:
            line = f_prm.readline()
            if line == '[ atoms ]\n':
                line = f_prm.readline()
                type = line.split()
                break
        while True:
            line = f_prm.readline()
            if line == '[ charges ]\n':
                line = f_prm.readline()
                charge = line.split()
                break
        if type != [] and charge != []:
            for i in range(len(type)):
                if type[i] == 's':
                    s_charge = charge[i]
                if type[i] == 'ss':
                    ss_charge = charge[i]
                if s_charge is not None and ss_charge is not None:
                    break

            for i in range(Molrange[0] - 1, Molrange[1]):
                if atom[i][1] != 's' and atom[i][1] != 'ss':
                    atom[i][1] = type[i - Molrange[0] + 1]
                    atom[i][6] = charge[i - Molrange[0] + 1]
                if atom[i][1] == 's':
                    atom[i][6] = s_charge
                    s_id = i
                if atom[i][1] == 'ss':
                    atom[i][6] = ss_charge
                    ss_id = i
            for i in range(len(atom[0])):
                if i != 0 and i != 5:
                    a = atom[s_id][i]
                    atom[s_id][i] = atom[ss_id][i]
                    atom[ss_id][i] = a
    f_prm.close()


def _PISAchangeffbond(Natom, adjacent, reac, atomtype, breakbond):
    bond = []
    angle = []
    dihedral = []
    improper = []
    fin_itp = open('./gromacs/step6/PISA.itp', 'r')

    while True:
        line = fin_itp.readline()
        if line == '[ bonds ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            bond.append(lines)

    if breakbond != []:
        for bre in breakbond[:]:
            for bon in bond[:]:
                unbond = [False for b in bre if b not in bon]
                if not unbond:
                    index = bond.index(bon)
                    del bond[index]

    while True:
        line = fin_itp.readline()
        if line == '[ angles ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            angle.append(lines)

    if breakbond != []:
        for bre in breakbond[:]:
            for ang in angle[:]:
                unbond = [False for b in bre if b not in ang]
                if not unbond:
                    index = angle.index(ang)
                    del angle[index]

    while True:
        line = fin_itp.readline()
        if line == '[ dihedrals ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            dihedral.append(lines)

    if breakbond != []:
        for bre in breakbond[:]:
            for dih in dihedral[:]:
                unbond = [False for b in bre if b not in dih]
                if not unbond:
                    index = dihedral.index(dih)
                    del dihedral[index]

    while True:
        line = fin_itp.readline()
        if line == '[ dihedrals ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            improper.append(lines)
    fin_itp.close()

    ffbond = _PISAspecialffbond('bonds')
    ffangle = _PISAspecialffbond('angles')
    ffdihedral = _PISAspecialffbond('dihedrals')
    # bond
    for i in range(len(reac)):
        ffbond_M = ffbond
        special = False
        if ffbond != []:
            if int(reac[i][0]) < int(reac[i][1]):
                a = atomtype[reac[i][0]]
                b = atomtype[reac[i][1]]
                for j in range(len(ffbond)):
                    if ffbond_M[j][0] == a and ffbond_M[j][1] == b:
                        temp = ffbond_M[j]
                        temp[j][0] = reac[i][0]
                        temp[j][1] = reac[i][1]
                        bond.append(temp)
                        if (j + 1) < len(ffbond):
                            if ffbond_M[j][0] == ffbond_M[
                                    j +
                                    1][0] and ffbond_M[j][1] == ffbond_M[j +
                                                                         1][1]:
                                ffbond_M.pop(j)
                        special = True
                        break
                if not special:
                    temp = [''] * 3
                    temp[0] = reac[i][0]
                    temp[1] = reac[i][1]
                    temp[2] = '1'
                    bond.append(temp)
            if int(reac[i][0]) > int(reac[i][1]):
                a = atomtype[reac[i][1]]
                b = atomtype[reac[i][0]]
                for j in range(len(ffbond)):
                    if ffbond_M[j][0] == a and ffbond_M[j][1] == b:
                        temp = ffbond_M[j]
                        temp[j][0] = reac[i][1]
                        temp[j][1] = reac[i][0]
                        bond.append(temp)
                        if (j + 1) < len(ffbond):
                            if ffbond_M[j][0] == ffbond_M[
                                    j +
                                    1][0] and ffbond_M[j][1] == ffbond_M[j +
                                                                         1][1]:
                                ffbond_M.pop(j)
                        special = True
                        break
                    if not special:
                        temp = [''] * 3
                        temp[0] = reac[i][1]
                        temp[1] = reac[i][0]
                        temp[2] = '1'
                        bond.append(temp)
            ffbond_M = ffbond
        else:
            if int(reac[i][0]) < int(reac[i][1]):
                temp = [''] * 3
                temp[0] = reac[i][0]
                temp[1] = reac[i][1]
                temp[2] = '1'
                bond.append(temp)
            if int(reac[i][0]) > int(reac[i][1]):
                temp = [''] * 3
                temp[0] = reac[i][1]
                temp[1] = reac[i][0]
                temp[2] = '1'
                bond.append(temp)
    if ffbond != []:
        for i in range(len(bond)):
            special = True
            for j in range(len(ffbond)):
                if atomtype[bond[i][0]] == ffbond[j][0] and atomtype[
                        bond[i][1]] == ffbond[j][1]:
                    special = False
                    break
            if special and len(bond[i]) > 3:
                temp = [''] * 3
                temp[0] = bond[i][0]
                temp[1] = bond[i][1]
                temp[2] = '1'
                bond[i] = temp
    bond.sort(key=lambda x: int(x[0]))

    # angle
    for i in range(len(reac)):
        special = False
        ffangle_L = ffangle
        ffangle_R = ffangle
        if ffangle != []:
            for j in range(Natom):
                if adjacent.has_edge(int(reac[i][0]) - 1, j):
                    if (j + 1) < int(reac[i][1]):
                        a = atomtype[str(j + 1)]
                        b = atomtype[reac[i][0]]
                        c = atomtype[reac[i][1]]
                        for k in range(len(ffangle)):
                            if ffangle_L[k][0] == a and ffangle_L[k][
                                    1] == b and ffangle_L[k][2] == c:
                                temp = ffangle_L[k]
                                temp[0] = str(j + 1)
                                temp[1] = reac[i][0]
                                temp[2] = reac[i][1]
                                angle.append(temp)
                                if (k + 1) < len(ffangle):
                                    if ffangle_L[k][0] == ffangle_L[
                                            k + 1][0] and ffangle_L[k][
                                                1] == ffangle_L[
                                                    k + 1][1] and ffangle_L[k][
                                                        2] == ffangle_L[k +
                                                                        1][2]:
                                        ffangle_L.pop(k)
                                special = True
                                break
                        if not special:
                            temp = [''] * 4
                            temp[0] = str(j + 1)
                            temp[1] = reac[i][0]
                            temp[2] = reac[i][1]
                            temp[3] = '1'
                            angle.append(temp)
                    if (j + 1) > int(reac[i][1]):
                        a = atomtype[reac[i][1]]
                        b = atomtype[reac[i][0]]
                        c = atomtype[str(j)]
                        for k in range(len(ffangle)):
                            if ffangle_L[k][0] == a and ffangle_L[k][
                                    1] == b and ffangle_L[k][2] == c:
                                temp = ffangle_L[k]
                                temp[0] = reac[i][1]
                                temp[1] = reac[i][0]
                                temp[2] = str(j + 1)
                                angle.append(temp)
                                if (k + 1) < len(ffangle):
                                    if ffangle_L[k][0] == ffangle_L[
                                            k + 1][0] and ffangle_L[k][
                                                1] == ffangle_L[
                                                    k + 1][1] and ffangle_L[k][
                                                        2] == ffangle_L[k +
                                                                        1][2]:
                                        ffangle_L.pop(k)
                                special = True
                                break
                        if not special:
                            temp = [''] * 4
                            temp[0] = reac[i][1]
                            temp[1] = reac[i][0]
                            temp[2] = str(j + 1)
                            temp[3] = '1'
                            angle.append(temp)
                if adjacent.has_edge(int(reac[i][1]) - 1, j):
                    if (j + 1) < int(reac[i][0]):
                        a = atomtype[str(j + 1)]
                        b = atomtype[reac[i][1]]
                        c = atomtype[reac[i][0]]
                        for k in range(len(ffangle)):
                            if ffangle_R[k][0] == a and ffangle_R[k][
                                    1] == b and ffangle_R[k][2] == c:
                                temp = ffangle_R[k]
                                temp[0] = str(j + 1)
                                temp[1] = reac[i][1]
                                temp[2] = reac[i][0]
                                angle.append(temp)
                                if (k + 1) < len(ffangle):
                                    if ffangle_R[k][0] == ffangle_R[
                                            k + 1][0] and ffangle_R[k][
                                                1] == ffangle_R[
                                                    k + 1][1] and ffangle_R[k][
                                                        2] == ffangle_R[k +
                                                                        1][2]:
                                        ffangle_R.pop(k)
                                special = True
                                break
                        if not special:
                            temp = [''] * 4
                            temp[0] = str(j + 1)
                            temp[1] = reac[i][1]
                            temp[2] = reac[i][0]
                            temp[3] = '1'
                            angle.append(temp)
                    if (j + 1) > int(reac[i][0]):
                        a = atomtype[reac[i][1]]
                        b = atomtype[reac[i][0]]
                        c = atomtype[str(j)]
                        for k in range(len(ffangle)):
                            if ffangle_R[k][0] == a and ffangle_R[k][
                                    1] == b and ffangle_R[k][2] == c:
                                temp = ffangle_R[k]
                                temp[0] = reac[i][0]
                                temp[1] = reac[i][1]
                                temp[2] = str(j + 1)
                                angle.append(temp)
                                if (k + 1) < len(ffangle):
                                    if ffangle_R[k][0] == ffangle_R[
                                            k + 1][0] and ffangle_R[k][
                                                1] == ffangle_R[
                                                    k + 1][1] and ffangle_R[k][
                                                        2] == ffangle_R[k +
                                                                        1][2]:
                                        ffangle_R.pop(k)
                                special = True
                                break
                        if not special:
                            temp = [''] * 4
                            temp[0] = reac[i][0]
                            temp[1] = reac[i][1]
                            temp[2] = str(j + 1)
                            temp[3] = '1'
                            angle.append(temp)
            ffangle_R = ffangle
            ffangle_L = ffangle
        else:
            if adjacent.has_edge(int(reac[i][0]) - 1, j):
                if (j + 1) < int(reac[i][1]):
                    temp = [''] * 4
                    temp[0] = str(j + 1)
                    temp[1] = reac[i][0]
                    temp[2] = reac[i][1]
                    temp[3] = '1'
                    angle.append(temp)
                if (j + 1) > int(reac[i][1]):
                    temp = [''] * 4
                    temp[0] = reac[i][1]
                    temp[1] = reac[i][0]
                    temp[2] = str(j + 1)
                    temp[3] = '1'
                    angle.append(temp)
            if adjacent.has_edge(int(reac[i][1]) - 1, j):
                if (j + 1) < int(reac[i][0]):
                    temp = [''] * 4
                    temp[0] = str(j + 1)
                    temp[1] = reac[i][1]
                    temp[2] = reac[i][0]
                    temp[3] = '1'
                    angle.append(temp)
                if (j + 1) > int(reac[i][0]):
                    temp = [''] * 4
                    temp[0] = reac[i][0]
                    temp[1] = reac[i][1]
                    temp[2] = str(j + 1)
                    temp[3] = '1'
                    angle.append(temp)
    if ffangle != []:
        for i in range(len(angle)):
            special = True
            for j in range(len(ffangle)):
                if atomtype[angle[i][0]] == ffangle[j][0] and atomtype[
                        angle[i][1]] == ffangle[j][1] and atomtype[
                            angle[i][2]] == ffangle[j][2]:
                    special = False
                    break
            if special and len(angle[i]) > 4:
                temp = [''] * 4
                temp[0] = angle[i][0]
                temp[1] = angle[i][1]
                temp[2] = angle[i][2]
                temp[3] = '1'
                angle[i] = temp
    angle.sort(key=lambda x: int(x[0]))

    # dihedral
    for i in range(len(reac)):
        special = False
        ffdihedral_M = ffdihedral
        ffdihedral_L = ffdihedral
        ffdihedral_R = ffdihedral
        if ffdihedral != []:
            for j in range(Natom):
                if adjacent.has_edge(int(reac[i][0]) - 1,
                                     j) and (j + 1) != int(reac[i][1]):
                    for k in range(Natom):
                        if adjacent.has_edge(int(reac[i][1]) - 1,
                                             k) and (k + 1) != int(reac[i][0]):
                            if j < k:
                                a = atomtype[str(j + 1)]
                                b = atomtype[reac[i][0]]
                                c = atomtype[reac[i][1]]
                                d = atomtype[str(k + 1)]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_M[ff][
                                            0] == a and ffdihedral_M[ff][
                                                1] == b and ffdihedral_M[ff][
                                                    2] == c and ffdihedral_M[
                                                        ff][3] == d:
                                        temp = ffdihedral_M[ff]
                                        temp[0] = str(j + 1)
                                        temp[1] = reac[i][0]
                                        temp[2] = reac[i][1]
                                        temp[3] = str(k + 1)
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_M[ff][0] == ffdihedral_M[
                                                    ff +
                                                    1][0] and ffdihedral_M[ff][
                                                        1] == ffdihedral_M[
                                                            ff +
                                                            1][1] and ffdihedral_M[
                                                                ff][2] == ffdihedral_M[
                                                                    ff +
                                                                    1][2] and ffdihedral_M[
                                                                        ff][3] == ffdihedral_M[
                                                                            ff
                                                                            +
                                                                            1][3]:
                                                ffdihedral_M.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = str(j + 1)
                                    temp[1] = reac[i][0]
                                    temp[2] = reac[i][1]
                                    temp[3] = str(k + 1)
                                    temp[4] = '9'
                                    dihedral.append(temp)
                            if j > k:
                                a = atomtype[str(k + 1)]
                                b = atomtype[reac[i][1]]
                                c = atomtype[reac[i][0]]
                                d = atomtype[str(j + 1)]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_M[ff][
                                            0] == a and ffdihedral_M[ff][
                                                1] == b and ffdihedral_M[ff][
                                                    2] == c and ffdihedral_M[
                                                        ff][3] == d:
                                        temp = ffdihedral_M[ff]
                                        temp[0] = str(k + 1)
                                        temp[1] = reac[i][1]
                                        temp[2] = reac[i][0]
                                        temp[3] = str(j + 1)
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_M[ff][0] == ffdihedral_M[
                                                    ff +
                                                    1][0] and ffdihedral_M[ff][
                                                        1] == ffdihedral_M[
                                                            ff +
                                                            1][1] and ffdihedral_M[
                                                                ff][2] == ffdihedral_M[
                                                                    ff +
                                                                    1][2] and ffdihedral_M[
                                                                        ff][3] == ffdihedral_M[
                                                                            ff
                                                                            +
                                                                            1][3]:
                                                ffdihedral_M.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = str(k + 1)
                                    temp[1] = reac[i][1]
                                    temp[2] = reac[i][0]
                                    temp[3] = str(j + 1)
                                    temp[4] = '9'
                                    dihedral.append(temp)
                    for k in range(Natom):
                        if adjacent.has_edge(j,
                                             k) and (k + 1) != int(reac[i][0]):
                            if (k + 1) < int(reac[i][1]):
                                a = atomtype[str(k + 1)]
                                b = atomtype[str(j + 1)]
                                c = atomtype[reac[i][0]]
                                d = atomtype[reac[i][1]]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_L[ff][
                                            0] == a and ffdihedral_L[ff][
                                                1] == b and ffdihedral_L[ff][
                                                    2] == c and ffdihedral_L[
                                                        ff][3] == d:
                                        temp = ffdihedral_L[ff]
                                        temp[0] = str(k + 1)
                                        temp[1] = str(j + 1)
                                        temp[2] = reac[i][0]
                                        temp[3] = reac[i][1]
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_L[ff][0] == ffdihedral_L[
                                                    ff +
                                                    1][0] and ffdihedral_L[ff][
                                                        1] == ffdihedral_L[
                                                            ff +
                                                            1][1] and ffdihedral_L[
                                                                ff][2] == ffdihedral_L[
                                                                    ff +
                                                                    1][2] and ffdihedral_L[
                                                                        ff][3] == ffdihedral_L[
                                                                            ff
                                                                            +
                                                                            1][3]:
                                                ffdihedral_L.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = str(k + 1)
                                    temp[1] = str(j + 1)
                                    temp[2] = reac[i][0]
                                    temp[3] = reac[i][1]
                                    temp[4] = '9'
                                    dihedral.append(temp)
                            if (k + 1) > int(reac[i][1]):
                                a = atomtype[reac[i][1]]
                                b = atomtype[reac[i][0]]
                                c = atomtype[str(j + 1)]
                                d = atomtype[str(k + 1)]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_L[ff][
                                            0] == a and ffdihedral_L[ff][
                                                1] == b and ffdihedral_L[ff][
                                                    2] == c and ffdihedral_L[
                                                        ff][3] == d:
                                        temp = ffdihedral_L[ff]
                                        temp[0] = reac[i][1]
                                        temp[1] = reac[i][0]
                                        temp[2] = str(j + 1)
                                        temp[3] = str(k + 1)
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_L[ff][0] == ffdihedral_L[
                                                    ff +
                                                    1][0] and ffdihedral_L[ff][
                                                        1] == ffdihedral_L[
                                                            ff +
                                                            1][1] and ffdihedral_L[
                                                                ff][2] == ffdihedral_L[
                                                                    ff +
                                                                    1][2] and ffdihedral_L[
                                                                        ff][3] == ffdihedral_L[
                                                                            ff
                                                                            +
                                                                            1][3]:
                                                ffdihedral_L.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = reac[i][1]
                                    temp[1] = reac[i][0]
                                    temp[2] = str(j + 1)
                                    temp[3] = str(k + 1)
                                    temp[4] = '9'
                                    dihedral.append(temp)
                if adjacent.has_edge(int(reac[i][1]) - 1,
                                     j) and (j + 1) != int(reac[i][0]):
                    for k in range(Natom):
                        if adjacent.has_edge(j,
                                             k) and (k + 1) != int(reac[i][1]):
                            if (k + 1) < int(reac[i][0]):
                                a = atomtype[str(k + 1)]
                                b = atomtype[str(j + 1)]
                                c = atomtype[reac[i][1]]
                                d = atomtype[reac[i][0]]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_R[ff][
                                            0] == a and ffdihedral_R[ff][
                                                1] == b and ffdihedral_R[ff][
                                                    2] == c and ffdihedral_R[
                                                        ff][3] == d:
                                        temp = ffdihedral_R[ff]
                                        temp[0] = str(k + 1)
                                        temp[1] = str(j + 1)
                                        temp[2] = reac[i][1]
                                        temp[3] = reac[i][0]
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_R[ff][0] == ffdihedral_R[
                                                    ff +
                                                    1][0] and ffdihedral_R[ff][
                                                        1] == ffdihedral_R[
                                                            ff +
                                                            1][1] and ffdihedral_R[
                                                                ff][2] == ffdihedral_R[
                                                                    ff +
                                                                    1][2] and ffdihedral_R[
                                                                        ff][3] == ffdihedral_R[
                                                                            ff
                                                                            +
                                                                            1][3]:
                                                ffdihedral_R.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = str(k + 1)
                                    temp[1] = str(j + 1)
                                    temp[2] = reac[i][1]
                                    temp[3] = reac[i][0]
                                    temp[4] = '9'
                                    dihedral.append(temp)
                            if (k + 1) > int(reac[i][0]):
                                a = atomtype[reac[i][0]]
                                b = atomtype[reac[i][1]]
                                c = atomtype[str(j + 1)]
                                d = atomtype[str(k + 1)]
                                for ff in range(len(ffdihedral)):
                                    if ffdihedral_R[ff][
                                            0] == a and ffdihedral_R[ff][
                                                1] == b and ffdihedral_R[ff][
                                                    2] == c and ffdihedral_R[
                                                        ff][3] == d:
                                        temp = ffdihedral_R[ff]
                                        temp[0] = reac[i][0]
                                        temp[1] = reac[i][1]
                                        temp[2] = str(j + 1)
                                        temp[3] = str(k + 1)
                                        dihedral.append(temp)
                                        if (ff + 1) < len(ffdihedral):
                                            if ffdihedral_R[ff][0] == ffdihedral_R[
                                                    ff +
                                                    1][0] and ffdihedral_R[ff][
                                                        1] == ffdihedral_R[
                                                            ff +
                                                            1][1] and ffdihedral_R[
                                                                ff][2] == ffdihedral_R[
                                                                    ff +
                                                                    1][2] and ffdihedral_R[
                                                                        ff][3] == ffdihedral_R[
                                                                            ff
                                                                            +
                                                                            1][3]:
                                                ffdihedral_R.pop(ff)
                                        special = True
                                        break
                                if not special:
                                    temp = [''] * 5
                                    temp[0] = reac[i][0]
                                    temp[1] = reac[i][1]
                                    temp[2] = str(j + 1)
                                    temp[3] = str(k + 1)
                                    temp[4] = '9'
                                    dihedral.append(temp)
            ffdihedral_R = ffdihedral
            ffdihedral_L = ffdihedral
        else:
            if adjacent.has_edge(int(reac[i][0]) - 1,
                                 j) and (j + 1) != int(reac[i][1]):
                for k in range(Natom):
                    if adjacent.has_edge(j, k) and (k + 1) != int(reac[i][0]):
                        if (k + 1) < int(reac[i][1]):
                            temp = [''] * 5
                            temp[0] = str(k + 1)
                            temp[1] = str(j + 1)
                            temp[2] = reac[i][0]
                            temp[3] = reac[i][1]
                            temp[4] = '9'
                            dihedral.append(temp)
                        if (k + 1) > int(reac[i][1]):
                            temp = [''] * 5
                            temp[0] = reac[i][1]
                            temp[1] = reac[i][0]
                            temp[2] = str(j + 1)
                            temp[3] = str(k + 1)
                            temp[4] = '9'
                            dihedral.append(temp)
            if adjacent.has_edge(int(reac[i][1]) - 1,
                                 j) and (j + 1) != int(reac[i][0]):
                for k in range(Natom):
                    if adjacent.has_edge(j, k) and (k + 1) != int(reac[i][1]):
                        if (k + 1) < int(reac[i][0]):
                            temp = [''] * 5
                            temp[0] = str(k + 1)
                            temp[1] = str(j + 1)
                            temp[2] = reac[i][1]
                            temp[3] = reac[i][0]
                            temp[4] = '9'
                            dihedral.append(temp)
                        if (k + 1) > int(reac[i][0]):
                            temp = [''] * 5
                            temp[0] = reac[i][0]
                            temp[1] = reac[i][1]
                            temp[2] = str(j + 1)
                            temp[3] = str(k + 1)
                            temp[4] = '9'
                            dihedral.append(temp)
    if ffdihedral != []:
        for i in range(len(dihedral)):
            special = True
            for j in range(len(ffdihedral)):
                if atomtype[dihedral[i][0]] == ffdihedral[j][0] and atomtype[
                        dihedral[i][1]] == ffdihedral[j][1] and atomtype[
                            dihedral[i][2]] == ffdihedral[j][2] and atomtype[
                                dihedral[i][3]] == ffdihedral[j][3]:
                    special = False
                    break
            if special and len(dihedral[i]) > 5:
                temp = [''] * 5
                temp[0] = dihedral[i][0]
                temp[1] = dihedral[i][1]
                temp[2] = dihedral[i][2]
                temp[3] = dihedral[i][3]
                temp[4] = '9'
                dihedral[i] = temp
    dihedral.sort(key=lambda x: int(x[0]))

    # improper
    ffimproper = []
    fin_itp = open('./topper/ffbonded.itp', 'r')
    while True:
        line = fin_itp.readline()
        if line == '[ dihedraltypes ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '[ dihedraltypes ]\n':
            line = fin_itp.readline()
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            ffimproper.append(lines[:4])
    fin_itp.close()

    for i in range(len(reac)):
        bondatom_L = []
        bondatom_R = []
        for j in range(Natom):
            if adjacent.has_edge(int(reac[i][0]) - 1, j):
                bondatom_L.append(str(j + 1))
            if adjacent.has_edge(int(reac[i][1]) - 1, j):
                bondatom_R.append(str(j + 1))
        if len(bondatom_L) >= 2:
            for j in range(len(bondatom_L) - 1):
                for k in range(j + 1, len(bondatom_L)):
                    unsort = [0] * 3
                    unsort[0] = int(bondatom_L[j])
                    unsort[1] = int(bondatom_L[k])
                    unsort[2] = int(reac[i][1])
                    unsort.sort()
                    temp = [''] * 5
                    temp[0] = str(unsort[0])
                    temp[1] = str(unsort[1])
                    temp[2] = reac[i][0]
                    temp[3] = str(unsort[2])
                    temp[4] = '4'
                    improper.append(temp)
        if len(bondatom_R) >= 2:
            for j in range(len(bondatom_R) - 1):
                for k in range(j + 1, len(bondatom_R)):
                    unsort = [0] * 3
                    unsort[0] = int(bondatom_R[j])
                    unsort[1] = int(bondatom_R[k])
                    unsort[2] = int(reac[i][0])
                    unsort.sort()
                    temp = [''] * 5
                    temp[0] = str(unsort[0])
                    temp[1] = str(unsort[1])
                    temp[2] = reac[i][1]
                    temp[3] = str(unsort[2])
                    temp[4] = '4'
                    improper.append(temp)
    for imp in improper[:]:
        special = False
        for i in range(len(ffimproper)):
            if atomtype[imp[0]] == ffimproper[i][0] and atomtype[
                    imp[1]] == ffimproper[i][1] and atomtype[
                        imp[2]] == ffimproper[i][2] and atomtype[
                            imp[3]] == ffimproper[i][3]:
                special = True
                break
        if not special:
            index = improper.index(imp)
            del improper[index]
    improper.sort(key=lambda x: int(x[2]))
    return bond, angle, dihedral, improper


def RAFTsteponefile():
    # create step1 polymer system(simulation macro-CTA)
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    Pini = 0
    if 'Pini_in_macroCTA':
        Pini = system['Pini_in_macroCTA']
    else:
        print(
            'The proportion of INI in mactoCTA uses the default value! (Pini=0)'
        )
    if 'Version' in system:
        version = system['Version']

    nCTA, _, nMA, _ = _Npolymer()
    if nCTA == 0 and nMA == 0:
        return False
    nI = int(nCTA * Pini)
    nR = nCTA - nI
    nC = nCTA * 100
    if nR == 0:
        nR = 1
    if nI == 0:
        nI = 1
    Numatom = _Numatom()
    Rcharge = _getcharge('R', '2')
    Ccharge = _getcharge('C', '1')
    Icharge = _getcharge('I', '2')
    MAcharge = _getcharge('MA', '1')

    Ratom, Rbond, Rangle, Rproper, Rpair, Rimproper = _getstr('R')
    Iatom, Ibond, Iangle, Iproper, Ipair, Iimproper = _getstr('I')
    Catom, Cbond, Cangle, Cproper, Cpair, Cimproper = _getstr('C')
    MAatom, MAbond, MAangle, MAproper, MApair, MAimproper = _getstr('MA')
    reactpoint = _rectpoint()

    f_out = open('./topper/PISA.itp', 'w')
    print(';created with PISAMD version ' + version, file=f_out)
    print('\n[ moleculetype ]', file=f_out)
    print('PISA     3', file=f_out)

    print('\n[ atoms ]', file=f_out)

    for i in range(nI):
        for j in range(len(Iatom)):
            line = Iatom[j]
            if j == reactpoint[0][0]:
                print(str(int(line[0]) + i * Numatom[0]).rjust(6),
                      line[1].rjust(6),
                      str(i + 1).rjust(6),
                      'INI'.rjust(6),
                      'i_p'.rjust(6),
                      str(int(line[5]) + i * Numatom[0]).rjust(6),
                      Icharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)
            else:
                print(str(int(line[0]) + i * Numatom[0]).rjust(6),
                      line[1].rjust(6),
                      str(i + 1).rjust(6),
                      'INI'.rjust(6), (line[4] + line[0]).rjust(6),
                      str(int(line[5]) + i * Numatom[0]).rjust(6),
                      Icharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)

    for i in range(nR):
        for j in range(len(Ratom)):
            line = Ratom[j]
            if j == reactpoint[1][0]:
                print(str(int(line[0]) + nI * Numatom[0] +
                          i * Numatom[1]).rjust(6),
                      line[1].rjust(6),
                      str(nI + i + 1).rjust(6),
                      'CTA'.rjust(6),
                      'r_p'.rjust(6),
                      str(int(line[5]) + nI * Numatom[0] +
                          i * Numatom[1]).rjust(6),
                      Rcharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)
            else:
                print(str(int(line[0]) + nI * Numatom[0] +
                          i * Numatom[1]).rjust(6),
                      line[1].rjust(6),
                      str(nI + i + 1).rjust(6),
                      'CTA'.rjust(6), (line[4] + line[0]).rjust(6),
                      str(int(line[5]) + nI * Numatom[0] +
                          i * Numatom[1]).rjust(6),
                      Rcharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)

    for i in range(nC):
        for j in range(len(Catom)):
            line = Catom[j]
            if j == reactpoint[4][0]:
                print(str(
                    int(line[0]) + nI * Numatom[0] + nR * Numatom[1] +
                    i * Numatom[4]).rjust(6),
                      line[1].rjust(6),
                      str(nI + nR + i + 1).rjust(6),
                      'CTA'.rjust(6),
                      'c_x'.rjust(6),
                      str(
                          int(line[5]) + nI * Numatom[0] + nR * Numatom[1] +
                          i * Numatom[4]).rjust(6),
                      Ccharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)
            if j == reactpoint[4][1]:
                print(str(
                    int(line[0]) + nI * Numatom[0] + nR * Numatom[1] +
                    i * Numatom[4]).rjust(6),
                      line[1].rjust(6),
                      str(nI + nR + i + 1).rjust(6),
                      'CTA'.rjust(6),
                      'c_n'.rjust(6),
                      str(
                          int(line[5]) + nI * Numatom[0] + nR * Numatom[1] +
                          i * Numatom[4]).rjust(6),
                      Ccharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)
            if j != reactpoint[4][0] and j != reactpoint[4][1]:
                print(str(
                    int(line[0]) + nI * Numatom[0] + nR * Numatom[1] +
                    i * Numatom[4]).rjust(6),
                      line[1].rjust(6),
                      str(nI + nR + i + 1).rjust(6),
                      'CTA'.rjust(6), (line[4] + line[0]).rjust(6),
                      str(
                          int(line[5]) + nI * Numatom[0] + nR * Numatom[1] +
                          i * Numatom[4]).rjust(6),
                      Ccharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)

    for i in range(nMA):
        for j in range(len(MAatom)):
            line = MAatom[j]
            if j != reactpoint[2][0] and j != reactpoint[2][1]:
                print(str(
                    int(line[0]) + nI * Numatom[0] + nR * Numatom[1] +
                    nC * Numatom[4] + i * Numatom[2]).rjust(6),
                      line[1].rjust(6),
                      str(nI + nR + nC + i + 1).rjust(6),
                      'MOA'.rjust(6), (line[4] + line[0]).rjust(6),
                      str(
                          int(line[5]) + nI * Numatom[0] + nR * Numatom[1] +
                          nC * Numatom[4] + i * Numatom[2]).rjust(6),
                      MAcharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)
            if j == reactpoint[2][0]:
                print(str(
                    int(line[0]) + nI * Numatom[0] + nR * Numatom[1] +
                    nC * Numatom[4] + i * Numatom[2]).rjust(6),
                      line[1].rjust(6),
                      str(nI + nR + nC + i + 1).rjust(6),
                      'MOA'.rjust(6),
                      'm_n'.rjust(6),
                      str(
                          int(line[5]) + nI * Numatom[0] + nR * Numatom[1] +
                          nC * Numatom[4] + i * Numatom[2]).rjust(6),
                      MAcharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)
            if j == reactpoint[2][1]:
                print(str(
                    int(line[0]) + nI * Numatom[0] + nR * Numatom[1] +
                    nC * Numatom[4] + i * Numatom[2]).rjust(6),
                      line[1].rjust(6),
                      str(nI + nR + nC + i + 1).rjust(6),
                      'MOA'.rjust(6),
                      'm_p'.rjust(6),
                      str(
                          int(line[5]) + nI * Numatom[0] + nR * Numatom[1] +
                          nC * Numatom[4] + i * Numatom[2]).rjust(6),
                      MAcharge[j].rjust(14),
                      line[7].rjust(12),
                      file=f_out)

    if Ibond != [] or Rbond != [] or Cbond != [] or MAbond != []:
        print('\n[ bonds ]', file=f_out)
        if Ibond != []:
            for i in range(nI):
                for j in range(len(Ibond)):
                    if 'mSeminario' in Ibond[j]:
                        lines = Ibond[j]
                        print(str(int(lines[0]) + i * Numatom[0]).rjust(5),
                              str(int(lines[1]) + i * Numatom[0]).rjust(6),
                              lines[2].rjust(6),
                              lines[3].rjust(12),
                              str('{:.2f}'.format(float(lines[4]))).rjust(14),
                              file=f_out)
                    else:
                        lines = Ibond[j]
                        print(str(int(lines[0]) + i * Numatom[0]).rjust(5),
                              str(int(lines[1]) + i * Numatom[0]).rjust(6),
                              lines[2].rjust(6),
                              file=f_out)

        if Rbond != []:
            for i in range(nR):
                for j in range(len(Rbond)):
                    if 'mSeminario' in Rbond[j]:
                        lines = Rbond[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              lines[2].rjust(6),
                              lines[3].rjust(12),
                              str('{:.2f}'.format(float(lines[4]))).rjust(14),
                              file=f_out)
                    else:
                        lines = Rbond[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              lines[2].rjust(6),
                              file=f_out)
        if Cbond != []:
            for i in range(nC):
                for j in range(len(Cbond)):
                    if 'mSeminario' in Cbond[j]:
                        lines = Cbond[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            i * Numatom[4]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              lines[2].rjust(6),
                              lines[3].rjust(12),
                              str('{:.2f}'.format(float(lines[4]))).rjust(14),
                              file=f_out)
                    else:
                        lines = Cbond[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            i * Numatom[4]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              lines[2].rjust(6),
                              file=f_out)

        if MAbond != []:
            for i in range(nMA):
                for j in range(len(MAbond)):
                    if 'mSeminario' in MAbond[j]:
                        lines = MAbond[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              lines[2].rjust(6),
                              lines[3].rjust(12),
                              str('{:.2f}'.format(float(lines[4]))).rjust(14),
                              file=f_out)
                    else:
                        lines = MAbond[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              lines[2].rjust(6),
                              file=f_out)

    if Iangle != [] or Rangle != [] or Cangle != [] or MAangle != []:
        print('\n[ angles ]', file=f_out)
        if Iangle != []:
            for i in range(nI):
                for j in range(len(Iangle)):
                    if 'mSeminario' in Iangle[j]:
                        lines = Iangle[j]
                        print(str(int(lines[0]) + i * Numatom[0]).rjust(5),
                              str(int(lines[1]) + i * Numatom[0]).rjust(6),
                              str(int(lines[2]) + i * Numatom[0]).rjust(6),
                              lines[3].rjust(6),
                              lines[4].rjust(12),
                              str('{:.5f}'.format(float(lines[5]))).rjust(14),
                              file=f_out)
                    else:
                        lines = Iangle[j]
                        print(str(int(lines[0]) + i * Numatom[0]).rjust(5),
                              str(int(lines[1]) + i * Numatom[0]).rjust(6),
                              str(int(lines[2]) + i * Numatom[0]).rjust(6),
                              lines[3].rjust(6),
                              file=f_out)

        if Rangle != []:
            for i in range(nR):
                for j in range(len(Rangle)):
                    if 'mSeminario' in Rangle[j]:
                        lines = Rangle[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              lines[3].rjust(6),
                              lines[4].rjust(12),
                              str('{:.5f}'.format(float(lines[5]))).rjust(14),
                              file=f_out)
                    else:
                        lines = Rangle[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              lines[3].rjust(6),
                              file=f_out)

        if Cangle != []:
            for i in range(nC):
                for j in range(len(Cangle)):
                    if 'mSeminario' in Cangle[j]:
                        lines = Cangle[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            i * Numatom[4]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              lines[3].rjust(6),
                              lines[4].rjust(12),
                              str('{:.5f}'.format(float(lines[5]))).rjust(14),
                              file=f_out)
                    else:
                        lines = Cangle[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            i * Numatom[4]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              lines[3].rjust(6),
                              file=f_out)

        if MAangle != []:
            for i in range(nMA):
                for j in range(len(MAangle)):
                    if 'mSeminario' in MAangle[j]:
                        lines = MAangle[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              lines[3].rjust(6),
                              lines[4].rjust(12),
                              str('{:.5f}'.format(float(lines[5]))).rjust(14),
                              file=f_out)
                    else:
                        lines = MAangle[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              lines[3].rjust(6),
                              file=f_out)

    if Iproper != [] or Rproper != [] or Cproper != [] or MAproper != []:
        print('\n[ dihedrals ]', file=f_out)
        if Iproper != []:
            for i in range(nI):
                for j in range(len(Iproper)):
                    if 'mSeminario' in Iproper[j]:
                        lines = Iproper[j]
                        print(str(int(lines[0]) + i * Numatom[0]).rjust(5),
                              str(int(lines[1]) + i * Numatom[0]).rjust(6),
                              str(int(lines[2]) + i * Numatom[0]).rjust(6),
                              str(int(lines[3]) + i * Numatom[0]).rjust(6),
                              lines[4].rjust(6),
                              lines[5].rjust(12),
                              str('{:.5f}'.format(float(lines[6]))).rjust(14),
                              file=f_out)
                    else:
                        if j != 0:
                            line1 = Iproper[j - 1]
                            line2 = Iproper[j]
                            lines = line2
                            if line1[0] != line2[0] or \
                               line1[1] != line2[1] or \
                               line1[2] != line2[2] or \
                               line1[3] != line2[3] or \
                               line1[4] != line2[4]:
                                print(str(int(lines[0]) +
                                          i * Numatom[0]).rjust(5),
                                      str(int(lines[1]) +
                                          i * Numatom[0]).rjust(6),
                                      str(int(lines[2]) +
                                          i * Numatom[0]).rjust(6),
                                      str(int(lines[3]) +
                                          i * Numatom[0]).rjust(6),
                                      lines[4].rjust(6),
                                      file=f_out)
                        else:
                            lines = Iproper[j]
                            print(str(int(lines[0]) + i * Numatom[0]).rjust(5),
                                  str(int(lines[1]) + i * Numatom[0]).rjust(6),
                                  str(int(lines[2]) + i * Numatom[0]).rjust(6),
                                  str(int(lines[3]) + i * Numatom[0]).rjust(6),
                                  lines[4].rjust(6),
                                  file=f_out)

        if Rproper != []:
            for i in range(nR):
                for j in range(len(Rproper)):
                    if 'mSeminario' in Rproper[j]:
                        lines = Rproper[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              str(
                                  int(lines[3]) + nI * Numatom[0] +
                                  i * Numatom[1]).rjust(6),
                              lines[4].rjust(6),
                              lines[5].rjust(12),
                              str('{:.5f}'.format(float(lines[6]))).rjust(14),
                              file=f_out)
                    else:
                        if j != 0:
                            line1 = Rproper[j - 1]
                            line2 = Rproper[j]
                            lines = line2
                            if line1[0] != line2[0] or \
                               line1[1] != line2[1] or \
                               line1[2] != line2[2] or \
                               line1[3] != line2[3] or \
                               line1[4] != line2[4]:
                                print(str(
                                    int(lines[0]) + nI * Numatom[0] +
                                    i * Numatom[1]).rjust(5),
                                      str(
                                          int(lines[1]) + nI * Numatom[0] +
                                          i * Numatom[1]).rjust(6),
                                      str(
                                          int(lines[2]) + nI * Numatom[0] +
                                          i * Numatom[1]).rjust(6),
                                      str(
                                          int(lines[3]) + nI * Numatom[0] +
                                          i * Numatom[1]).rjust(6),
                                      lines[4].rjust(6),
                                      file=f_out)
                        else:
                            lines = Rproper[j]
                            print(str(
                                int(lines[0]) + nI * Numatom[0] +
                                i * Numatom[1]).rjust(5),
                                  str(
                                      int(lines[1]) + nI * Numatom[0] +
                                      i * Numatom[1]).rjust(6),
                                  str(
                                      int(lines[2]) + nI * Numatom[0] +
                                      i * Numatom[1]).rjust(6),
                                  str(
                                      int(lines[3]) + nI * Numatom[0] +
                                      i * Numatom[1]).rjust(6),
                                  lines[4].rjust(6),
                                  file=f_out)

        if Cproper != []:
            for i in range(nC):
                for j in range(len(Cproper)):
                    if 'mSeminario' in Cproper[j]:
                        lines = Cproper[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            i * Numatom[4]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              str(
                                  int(lines[3]) + nI * Numatom[0] +
                                  nR * Numatom[1] + i * Numatom[4]).rjust(6),
                              lines[4].rjust(6),
                              lines[5].rjust(12),
                              str('{:.5f}'.format(float(lines[6]))).rjust(14),
                              file=f_out)
                    else:
                        if j != 0:
                            line1 = Cproper[j - 1]
                            line2 = Cproper[j]
                            lines = line2
                            if line1[0] != line2[0] or \
                               line1[1] != line2[1] or \
                               line1[2] != line2[2] or \
                               line1[3] != line2[3] or \
                               line1[4] != line2[4]:
                                print(str(
                                    int(lines[0]) + nI * Numatom[0] +
                                    nR * Numatom[1] + i * Numatom[4]).rjust(5),
                                      str(
                                          int(lines[1]) + nI * Numatom[0] +
                                          nR * Numatom[1] +
                                          i * Numatom[4]).rjust(6),
                                      str(
                                          int(lines[2]) + nI * Numatom[0] +
                                          nR * Numatom[1] +
                                          i * Numatom[4]).rjust(6),
                                      str(
                                          int(lines[3]) + nI * Numatom[0] +
                                          nR * Numatom[1] +
                                          i * Numatom[4]).rjust(6),
                                      lines[4].rjust(6),
                                      file=f_out)
                        else:
                            lines = Cproper[j]
                            print(
                                str(
                                    int(lines[0]) + nI * Numatom[0] +
                                    nR * Numatom[1] + i * Numatom[4]).rjust(5),
                                str(
                                    int(lines[1]) + nI * Numatom[0] +
                                    nR * Numatom[1] + i * Numatom[4]).rjust(6),
                                str(
                                    int(lines[2]) + nI * Numatom[0] +
                                    nR * Numatom[1] + i * Numatom[4]).rjust(6),
                                str(
                                    int(lines[3]) + nI * Numatom[0] +
                                    nR * Numatom[1] + i * Numatom[4]).rjust(6),
                                lines[4].rjust(6),
                                file=f_out)

        if MAproper != []:
            for i in range(nMA):
                for j in range(len(MAproper)):
                    if 'mSeminario' in MAproper[j]:
                        lines = MAproper[j]
                        print(str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(5),
                              str(
                                  int(lines[1]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              str(
                                  int(lines[2]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              str(
                                  int(lines[3]) + nI * Numatom[0] +
                                  nR * Numatom[1] + nC * Numatom[4] +
                                  i * Numatom[2]).rjust(6),
                              lines[4].rjust(6),
                              lines[5].rjust(12),
                              str('{:.5f}'.format(float(lines[6]))).rjust(14),
                              file=f_out)
                    else:
                        if j != 0:
                            line1 = MAproper[j - 1]
                            line2 = MAproper[j]
                            lines = line2
                            if line1[0] != line2[0] or \
                               line1[1] != line2[1] or \
                               line1[2] != line2[2] or \
                               line1[3] != line2[3] or \
                               line1[4] != line2[4]:
                                print(str(
                                    int(lines[0]) + nI * Numatom[0] +
                                    nR * Numatom[1] + nC * Numatom[4] +
                                    i * Numatom[2]).rjust(5),
                                      str(
                                          int(lines[1]) + nI * Numatom[0] +
                                          nR * Numatom[1] + nC * Numatom[4] +
                                          i * Numatom[2]).rjust(6),
                                      str(
                                          int(lines[2]) + nI * Numatom[0] +
                                          nR * Numatom[1] + nC * Numatom[4] +
                                          i * Numatom[2]).rjust(6),
                                      str(
                                          int(lines[3]) + nI * Numatom[0] +
                                          nR * Numatom[1] + nC * Numatom[4] +
                                          i * Numatom[2]).rjust(6),
                                      lines[4].rjust(6),
                                      file=f_out)
                        else:
                            lines = MAproper[j]
                            print(str(
                                int(lines[0]) + nI * Numatom[0] +
                                nR * Numatom[1] + nC * Numatom[4] +
                                i * Numatom[2]).rjust(5),
                                  str(
                                      int(lines[1]) + nI * Numatom[0] +
                                      nR * Numatom[1] + nC * Numatom[4] +
                                      i * Numatom[2]).rjust(6),
                                  str(
                                      int(lines[2]) + nI * Numatom[0] +
                                      nR * Numatom[1] + nC * Numatom[4] +
                                      i * Numatom[2]).rjust(6),
                                  str(
                                      int(lines[3]) + nI * Numatom[0] +
                                      nR * Numatom[1] + nC * Numatom[4] +
                                      i * Numatom[2]).rjust(6),
                                  lines[4].rjust(6),
                                  file=f_out)

    if Ipair != [] or Rpair != [] or Cpair != [] or MApair != []:
        print('\n[ pairs ]', file=f_out)
        if Ipair != []:
            for i in range(nI):
                for j in range(len(Ipair)):
                    lines = Ipair[j]
                    print(str(int(lines[0]) + i * Numatom[0]).rjust(5),
                          str(int(lines[1]) + i * Numatom[0]).rjust(6),
                          lines[2].rjust(6),
                          file=f_out)

        if Rpair != []:
            for i in range(nR):
                for j in range(len(Rpair)):
                    lines = Rpair[j]
                    print(
                        str(int(lines[0]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(5),
                        str(int(lines[1]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(6),
                        lines[2].rjust(6),
                        file=f_out)

        if Cpair != []:
            for i in range(nC):
                for j in range(len(Cpair)):
                    lines = Cpair[j]
                    print(str(
                        int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                        i * Numatom[4]).rjust(5),
                          str(
                              int(lines[1]) + nI * Numatom[0] +
                              nR * Numatom[1] + i * Numatom[4]).rjust(6),
                          lines[2].rjust(6),
                          file=f_out)

        if MApair != []:
            for i in range(nMA):
                for j in range(len(MApair)):
                    lines = MApair[j]
                    print(
                        str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(5),
                        str(
                            int(lines[1]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(6),
                        lines[2].rjust(6),
                        file=f_out)

    if Iimproper != [] or Rimproper != [] or Cimproper != [] or MAimproper != []:
        print('\n[ dihedrals ]', file=f_out)
        if Iimproper != []:
            for i in range(nI):
                for j in range(len(Iimproper)):
                    lines = Iimproper[j]
                    print(str(int(lines[0]) + i * Numatom[0]).rjust(5),
                          str(int(lines[1]) + i * Numatom[0]).rjust(6),
                          str(int(lines[2]) + i * Numatom[0]).rjust(6),
                          str(int(lines[3]) + i * Numatom[0]).rjust(6),
                          lines[4].rjust(6),
                          file=f_out)

        if Rimproper != []:
            for i in range(nR):
                for j in range(len(Rimproper)):
                    lines = Rimproper[j]
                    print(
                        str(int(lines[0]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(5),
                        str(int(lines[1]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(6),
                        str(int(lines[2]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(6),
                        str(int(lines[3]) + nI * Numatom[0] +
                            i * Numatom[1]).rjust(6),
                        lines[4].rjust(6),
                        file=f_out)

        if Cimproper != []:
            for i in range(nC):
                for j in range(len(Cimproper)):
                    lines = Cimproper[j]
                    print(str(
                        int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                        i * Numatom[4]).rjust(5),
                          str(
                              int(lines[1]) + nI * Numatom[0] +
                              nR * Numatom[1] + i * Numatom[4]).rjust(6),
                          str(
                              int(lines[2]) + nI * Numatom[0] +
                              nR * Numatom[1] + i * Numatom[4]).rjust(6),
                          str(
                              int(lines[3]) + nI * Numatom[0] +
                              nR * Numatom[1] + i * Numatom[4]).rjust(6),
                          lines[4].rjust(6),
                          file=f_out)

        if MAimproper != []:
            for i in range(nMA):
                for j in range(len(MAimproper)):
                    lines = MAimproper[j]
                    print(
                        str(
                            int(lines[0]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(5),
                        str(
                            int(lines[1]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(6),
                        str(
                            int(lines[2]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(6),
                        str(
                            int(lines[3]) + nI * Numatom[0] + nR * Numatom[1] +
                            nC * Numatom[4] + i * Numatom[2]).rjust(6),
                        lines[4].rjust(6),
                        file=f_out)

    f_Ipdb = open('./topper/I.pdb', 'w')
    f_Rpdb = open('./topper/R.pdb', 'w')
    f_Cpdb = open('./topper/C.pdb', 'w')
    f_MApdb = open('./topper/MA.pdb', 'w')
    I_matrix = _getconnectmatrix(Numatom[0], Ibond)
    R_matrix = _getconnectmatrix(Numatom[1], Rbond)
    C_matrix = _getconnectmatrix(Numatom[4], Cbond)
    MA_matrix = _getconnectmatrix(Numatom[2], MAbond)
    Ipos, Rpos, Cpos, MApos, _ = _getpos()

    print('REMARK  PISAMD PDB file', file=f_Ipdb)
    for i in range(len(Ipos)):
        lines = Ipos[i].split()
        if i == reactpoint[0][0]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ',
                  'i_p'.ljust(4, ' '),
                  ' ',
                  'INI'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_Ipdb)
        else:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ', (lines[0] + str(i + 1)).ljust(4, ' '),
                  ' ',
                  'INI'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_Ipdb)
    print('TER'.ljust(4, ' '),
          '  ',
          str(len(Ipos) + 1).rjust(5, ' '),
          sep='',
          file=f_Ipdb)
    if Ibond != []:
        for i in range(len(I_matrix)):
            print('CONECT',
                  str(i + 1).rjust(5, ' '),
                  sep='',
                  end='',
                  file=f_Ipdb)
            for j in range(len(I_matrix[i])):
                if I_matrix[i][j] == 1:
                    print(str(j + 1).rjust(5, ' '),
                          sep='',
                          end='',
                          file=f_Ipdb)
            print(file=f_Ipdb)
    print('END', file=f_Ipdb)

    print('REMARK  PISAMD PDB file', file=f_Rpdb)
    for i in range(len(Rpos)):
        lines = Rpos[i].split()
        if i == reactpoint[1][0]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ',
                  'r_p'.ljust(4, ' '),
                  ' ',
                  'CTA'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_Rpdb)
        else:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ', (lines[0] + str(i + 1)).ljust(4, ' '),
                  ' ',
                  'CTA'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_Rpdb)
    print('TER'.ljust(4, ' '),
          '  ',
          str(len(Rpos) + 1).rjust(5, ' '),
          sep='',
          file=f_Rpdb)
    if Rbond != []:
        for i in range(len(R_matrix)):
            print('CONECT',
                  str(i + 1).rjust(5, ' '),
                  sep='',
                  end='',
                  file=f_Rpdb)
            for j in range(len(R_matrix[i])):
                if R_matrix[i][j] == 1:
                    print(str(j + 1).rjust(5, ' '),
                          sep='',
                          end='',
                          file=f_Rpdb)
            print(file=f_Rpdb)
    print('END', file=f_Rpdb)

    print('REMARK  PISAMD PDB file', file=f_Cpdb)
    for i in range(len(Cpos)):
        lines = Cpos[i].split()
        if i == reactpoint[4][0]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ',
                  'c_x'.ljust(4, ' '),
                  ' ',
                  'CTA'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_Cpdb)
        if i == reactpoint[4][1]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ',
                  'c_n'.ljust(4, ' '),
                  ' ',
                  'CTA'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_Cpdb)
        if i != reactpoint[4][0] and i != reactpoint[4][1]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ', (lines[0] + str(i + 1)).ljust(4, ' '),
                  ' ',
                  'CTA'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_Cpdb)
    print('TER'.ljust(4, ' '),
          '  ',
          str(len(Cpos) + 1).rjust(5, ' '),
          sep='',
          file=f_Cpdb)
    if Cbond != []:
        for i in range(len(C_matrix)):
            print('CONECT',
                  str(i + 1).rjust(5, ' '),
                  sep='',
                  end='',
                  file=f_Cpdb)
            for j in range(len(C_matrix[i])):
                if C_matrix[i][j] == 1:
                    print(str(j + 1).rjust(5, ' '),
                          sep='',
                          end='',
                          file=f_Cpdb)
            print(file=f_Cpdb)
    print('END', file=f_Cpdb)

    print('REMARK  PISAMD PDB file', file=f_MApdb)
    for i in range(len(MApos)):
        lines = MApos[i].split()
        if i != reactpoint[2][0] and i != reactpoint[2][1]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ', (lines[0] + str(i + 1)).ljust(4, ' '),
                  ' ',
                  'MOA'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_MApdb)
        if i == reactpoint[2][0]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ',
                  'm_n'.ljust(4, ' '),
                  ' ',
                  'MOA'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_MApdb)
        if i == reactpoint[2][1]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ',
                  'm_p'.ljust(4, ' '),
                  ' ',
                  'MOA'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_MApdb)

    print('TER'.ljust(4, ' '),
          '  ',
          str(len(MApos) + 1).rjust(5, ' '),
          sep='',
          file=f_MApdb)
    if MAbond != []:
        for i in range(len(MA_matrix)):
            print('CONECT',
                  str(i + 1).rjust(5, ' '),
                  sep='',
                  end='',
                  file=f_MApdb)
            for j in range(len(MA_matrix[i])):
                if MA_matrix[i][j] == 1:
                    print(str(j + 1).rjust(5, ' '),
                          sep='',
                          end='',
                          file=f_MApdb)
            print(file=f_MApdb)
    print('END', file=f_MApdb)

    f_Ipdb.close()
    f_Rpdb.close()
    f_Cpdb.close()
    f_MApdb.close()
    f_sys.close()
    f_out.close()
    return True


def CTA_changeitp(Ilist, Rlist, Clist, MAlist, reac):

    fin_itp = open('./gromacs/step3/PISA.itp', 'r')
    atom = []
    atomtype = {}  # aid to atype
    reac_molname = [['MOA', 'MOA'], ['CTA', 'MOA'], ['INI', 'MOA']]
    reac_list = [['r_p', 'm_n'], ['r_p', 'c_n'], ['i_p', 'm_n'],
                 ['i_p', 'c_n']]

    # atom: atom_id atom_type mol_id mol_name atom_name index charge mass (str)
    while True:
        line = fin_itp.readline()
        if line == '[ atoms ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            lines = line.split()
            atom.append(lines)
            atomtype[lines[0]] = lines[1]

    # adjacent
    adjacent = nx.Graph()
    while True:
        line = fin_itp.readline()
        if line == '[ bonds ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            adjacent.add_edge(int(lines[0]) - 1, int(lines[1]) - 1)
    fin_itp.close()

    # change atom & adjacent
    for i in range(len(reac)):
        reac_group = []
        molname_group = []
        reacindex = None
        molnameindex = None
        reac_group.append(reac[i][4])
        reac_group.append(reac[i][5])
        reac_group.sort()
        molname_group.append(atom[int(reac[i][0]) - 1][3])
        molname_group.append(atom[int(reac[i][1]) - 1][3])
        molname_group.sort()
        for j in range(4):
            temp = reac_list[j]
            temp.sort()
            if reac_group == temp:
                reacindex = j
                break
        for j in range(3):
            temp = reac_molname[j]
            temp.sort()
            if molname_group == temp:
                molnameindex = j
                break

        mola_id = None
        molb_id = None
        mola_range = [0, 0]
        molb_range = [0, 0]

        # MA: state1-->state3
        if molnameindex == 0 and (reacindex == 0 or reacindex == 2):
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'm_n':
                mola_id = int(reac[i][3])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            if reac[i][4] == 'm_n':
                mola_id = int(reac[i][2])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            for k in range(len(MAlist)):
                if MAlist[k][0] == mola_id:
                    mola_range[0] = MAlist[k][1]
                    mola_range[1] = MAlist[k][2]
                    break
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'm_p':
                    if reacindex == 0:
                        atom[k][4] = 'r_p'
                    if reacindex == 2:
                        atom[k][4] = 'i_p'
                    break
            _CTAchangeMol(atom, mola_range, 'MA', 3)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]

    # MA: state1-->state2
        if molnameindex == 1 and reacindex == 0:
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'm_n':
                mola_id = int(reac[i][3])
                molb_id = int(reac[i][2])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            if reac[i][4] == 'm_n':
                mola_id = int(reac[i][2])
                molb_id = int(reac[i][3])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            for k in range(len(MAlist)):
                if MAlist[k][0] == mola_id:
                    mola_range[0] = MAlist[k][1]
                    mola_range[1] = MAlist[k][2]
                    break
            for k in range(len(Rlist)):
                if Rlist[k][0] == molb_id:
                    molb_range[0] = Rlist[k][1]
                    molb_range[1] = Rlist[k][2]
                    break
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'm_p':
                    atom[k][4] = 'r_p'
                    break
            _CTAchangeMol(atom, mola_range, 'MA', 2)
            _CTAchangeMol(atom, molb_range, 'R', 1)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molb_range[0] - 1, molb_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
        if molnameindex == 2 and reacindex == 2:
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'm_n':
                mola_id = int(reac[i][3])
                molb_id = int(reac[i][2])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            if reac[i][4] == 'm_n':
                mola_id = int(reac[i][2])
                molb_id = int(reac[i][3])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            for k in range(len(MAlist)):
                if MAlist[k][0] == mola_id:
                    mola_range[0] = MAlist[k][1]
                    mola_range[1] = MAlist[k][2]
                    break
            for k in range(len(Ilist)):
                if Ilist[k][0] == molb_id:
                    molb_range[0] = Ilist[k][1]
                    molb_range[1] = Ilist[k][2]
                    break
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'm_p':
                    atom[k][4] = 'i_p'
                    break
            _CTAchangeMol(atom, mola_range, 'MA', 2)
            _CTAchangeMol(atom, molb_range, 'I', 1)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molb_range[0] - 1, molb_range[1]):
                atomtype[atom[k][0]] = atom[k][1]

    # MA: state3<-->state4
        if (molnameindex == 1 and reacindex == 1) or (molnameindex == 1
                                                      and reacindex == 3):
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'c_n':
                mola_id = int(reac[i][3])
                molb_id = int(reac[i][2])
                if reacindex == 1:
                    atom[int(reac[i][0]) - 1][4] = 'r_b'
                if reacindex == 3:
                    atom[int(reac[i][0]) - 1][4] = 'i_b'
            if reac[i][4] == 'c_n':
                mola_id = int(reac[i][2])
                molb_id = int(reac[i][3])
                if reacindex == 1:
                    atom[int(reac[i][1]) - 1][4] = 'r_b'
                if reacindex == 3:
                    atom[int(reac[i][1]) - 1][4] = 'i_b'
            for k in range(len(Clist)):
                if Clist[k][0] == mola_id:
                    mola_range[0] = Clist[k][1]
                    mola_range[1] = Clist[k][2]
                    break
            for k in range(len(MAlist)):
                if MAlist[k][0] == molb_id:
                    molb_range[0] = MAlist[k][1]
                    molb_range[1] = MAlist[k][2]
                    break
            _CTAchangeMol(atom, mola_range, 'C', 2)
            _CTAchangeMol(atom, molb_range, 'MA', 4)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molb_range[0] - 1, molb_range[1]):
                atomtype[atom[k][0]] = atom[k][1]

    # change ffbond
    bond, angle, dihedral, improper = _CTAchangeffbond(len(atom), adjacent,
                                                       reac, atomtype)

    # write PISA.itp
    fout_itp = open('./gromacs/step3/PISA.itp', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=fout_itp)
    print('\n[ moleculetype ]', file=fout_itp)
    print('PISA     3', file=fout_itp)

    print('\n[ atoms ]', file=fout_itp)
    for i in range(len(atom)):
        print(atom[i][0].rjust(6),
              atom[i][1].rjust(6),
              atom[i][2].rjust(6),
              atom[i][3].rjust(6),
              atom[i][4].rjust(6),
              atom[i][5].rjust(6),
              atom[i][6].rjust(14),
              atom[i][7].rjust(12),
              file=fout_itp)
    print('\n[ bonds ]', file=fout_itp)
    for i in range(len(bond)):
        if len(bond[i]) > 3:
            print(bond[i][0].rjust(5),
                  bond[i][1].rjust(6),
                  bond[i][2].rjust(6),
                  bond[i][3].rjust(12),
                  bond[i][4].rjust(14),
                  file=fout_itp)
        else:
            print(bond[i][0].rjust(5),
                  bond[i][1].rjust(6),
                  bond[i][2].rjust(6),
                  file=fout_itp)
    print('\n[ angles ]', file=fout_itp)
    for i in range(len(angle)):
        if len(angle[i]) > 4:
            print(angle[i][0].rjust(5),
                  angle[i][1].rjust(6),
                  angle[i][2].rjust(6),
                  angle[i][3].rjust(6),
                  angle[i][4].rjust(12),
                  angle[i][5].rjust(14),
                  file=fout_itp)
        else:
            print(angle[i][0].rjust(5),
                  angle[i][1].rjust(6),
                  angle[i][2].rjust(6),
                  angle[i][3].rjust(6),
                  file=fout_itp)
    print('\n[ dihedrals ]', file=fout_itp)
    for i in range(len(dihedral)):
        if len(dihedral[i]) > 5:
            print(dihedral[i][0].rjust(5),
                  dihedral[i][1].rjust(6),
                  dihedral[i][2].rjust(6),
                  dihedral[i][3].rjust(6),
                  dihedral[i][4].rjust(6),
                  dihedral[i][5].rjust(12),
                  dihedral[i][6].rjust(14),
                  file=fout_itp)
        else:
            print(dihedral[i][0].rjust(5),
                  dihedral[i][1].rjust(6),
                  dihedral[i][2].rjust(6),
                  dihedral[i][3].rjust(6),
                  dihedral[i][4].rjust(6),
                  file=fout_itp)
    print('\n[ pairs ]', file=fout_itp)
    for i in range(len(dihedral)):
        print(dihedral[i][0].rjust(5),
              dihedral[i][3].rjust(6),
              '1'.rjust(6),
              file=fout_itp)
    print('\n[ dihedrals ]', file=fout_itp)
    for i in range(len(improper)):
        print(improper[i][0].rjust(5),
              improper[i][1].rjust(6),
              improper[i][2].rjust(6),
              improper[i][3].rjust(6),
              improper[i][4].rjust(6),
              file=fout_itp)
    fout_itp.close()
    f_sys.close()
    return atom


def CTA_changegro(frame, atom):

    system_name = ''
    Natom = ''
    boxsize = ''
    old_atom = []
    fin_gro = open('./gromacs/step3/' + str(frame) + '.gro', 'r')
    system_name = fin_gro.readline()
    Natom = fin_gro.readline()
    for i in range(int(Natom)):
        line = fin_gro.readline()
        temp = [''] * 10
        temp[0] = line[0:5].strip()
        temp[1] = line[5:10].strip()
        temp[2] = line[10:15].strip()
        temp[3] = line[15:20].strip()
        temp[4] = line[20:28].strip()
        temp[5] = line[28:36].strip()
        temp[6] = line[36:44].strip()
        temp[7] = line[44:52].strip()
        temp[8] = line[52:60].strip()
        temp[9] = line[60:68].strip()
        old_atom.append(temp)
    boxsize = fin_gro.readline()
    fin_gro.close()

    for i in range(len(atom)):
        old_atom[i][2] = atom[i][4]

    fout_gro = open('./gromacs/step3/' + str(frame) + '.gro', 'w')
    print(system_name, end='', file=fout_gro)
    print(Natom, end='', file=fout_gro)
    for i in range(len(old_atom)):
        print(old_atom[i][0].rjust(5),
              old_atom[i][1].ljust(5),
              old_atom[i][2].rjust(5),
              old_atom[i][3].rjust(5),
              old_atom[i][4].rjust(8),
              old_atom[i][5].rjust(8),
              old_atom[i][6].rjust(8),
              old_atom[i][7].rjust(8),
              old_atom[i][8].rjust(8),
              old_atom[i][9].rjust(8),
              sep='',
              file=fout_gro)
    print(boxsize, file=fout_gro)
    fout_gro.close()


def CTA_grotodata(frame):

    PISA_atom = []
    PISA_bond = []
    typelist = []
    boxsize = [0] * 3
    Natom = 0
    Nmol = 0
    Nbondtype = 0

    # atomtype to [index, mass]
    atomtype = {}
    fin_itp = open('./topper/ffnonbonded.itp', 'r')
    while True:
        line = fin_itp.readline()
        if line == '[ atomtypes ]\n':
            line = fin_itp.readline()
            break
    index = 1
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            temp = [''] * 2
            lines = line.split()
            temp[0] = str(index)
            temp[1] = lines[2]
            atomtype[lines[0]] = temp
            index += 1
            typelist.append(lines[0])
    fin_itp.close()

    # polymer atom and bond
    fin_itp = open('./gromacs/step3/PISA.itp', 'r')
    while True:
        line = fin_itp.readline()
        if line == '[ atoms ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            temp = [''] * 7  # atomid molid atomtype charge X Y Z
            lines = line.split()
            temp[0] = lines[0]
            temp[1] = lines[2]
            temp[2] = atomtype[lines[1]][0]
            temp[3] = lines[6]
            PISA_atom.append(temp)
    Natom = Natom + len(PISA_atom)
    Nmol = int(PISA_atom[-1][1])
    while True:
        line = fin_itp.readline()
        if line == '[ bonds ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            temp = [''] * 3  # atomaid atombid bondtype
            lines = line.split()
            temp[0] = lines[0]
            temp[1] = lines[1]
            temp[2] = '1'
            PISA_bond.append(temp)
    fin_itp.close()
    Nbondtype += 1

    # Ion/solvent atom and bond
    mol_list = []
    fin_top = open('./gromacs/step3/system.top', 'r')
    while True:
        line = fin_top.readline()
        if line == '[ molecules ]\n':
            break
    while True:
        line = fin_top.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            mol_list.append(lines)
    fin_top.close()
    for i in range(1, len(mol_list)):
        atom_temp = []
        bond_temp = []
        fin_itp = open('./topper/' + mol_list[i][0] + '.itp', 'r')
        while True:
            line = fin_itp.readline()
            if line == '[ atoms ]\n':
                break
        while True:
            line = fin_itp.readline()
            if line == '\n' or line == ' \n' or line == '':
                break
            else:
                lines = line.split()
                atom_temp.append(lines)
        while True:
            if line == '':
                break
            line = fin_itp.readline()
            if line == '[ bonds ]\n':
                break
        while True:
            line = fin_itp.readline()
            if line == '\n' or line == ' \n' or line == '':
                break
            else:
                lines = line.split()
                bond_temp.append(lines)
        for j in range(int(mol_list[i][1])):
            for k in range(len(atom_temp)):
                temp = [''] * 7  # atomid molid atomtype charge X Y Z
                temp[0] = int(atom_temp[k][0]) + Natom + j * len(atom_temp)
                temp[1] = int(atom_temp[k][2]) + Nmol + j
                temp[2] = atomtype[atom_temp[k][1]][0]
                temp[3] = atom_temp[k][6]
                PISA_atom.append(temp)
        if bond_temp != []:
            Nbondtype += 1
            for j in range(int(mol_list[i][1])):
                for k in range(len(bond_temp)):
                    temp = [''] * 3  # atomaid atombid bondtype
                    temp[0] = int(bond_temp[k][0]) + Natom + j * len(atom_temp)
                    temp[1] = int(bond_temp[k][1]) + Natom + j * len(atom_temp)
                    temp[2] = Nbondtype
                    PISA_bond.append(temp)
        Natom = Natom + int(mol_list[i][1]) * len(atom_temp)
        Nmol = Nmol + int(mol_list[i][1])
        fin_top.close()
    # atom position
    fin_gro = open('./gromacs/step3/' + str(frame) + '.gro')
    line = fin_gro.readline()
    line = fin_gro.readline()
    Total_atom = int(line)
    for i in range(Total_atom):
        line = fin_gro.readline()
        lines = line.split()
        PISA_atom[i][4] = lines[-4]
        PISA_atom[i][5] = lines[-5]
        PISA_atom[i][6] = lines[-6]
    line = fin_gro.readline()
    lines = line.split()
    boxsize[0] = lines[0]
    boxsize[1] = lines[1]
    boxsize[2] = lines[2]
    fin_gro.close()

    folder = os.path.exists('./gromacs/step3/traj')
    if not folder:
        os.makedirs('./gromacs/step3/traj')

    fout_data = open('./gromacs/step3/traj/' + str(frame) + '.data', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print('# LAMMPS data file written by PISAMD version ' + version,
          file=fout_data)
    print(len(PISA_atom), 'atoms', sep=' ', file=fout_data)
    print(len(PISA_bond), 'bonds', sep=' ', file=fout_data)
    print(len(typelist), 'atom types', sep=' ', file=fout_data)
    print(Nbondtype, 'bond types', sep=' ', file=fout_data)
    print('0.0', boxsize[0], 'xlo xhi', sep=' ', file=fout_data)
    print('0.0', boxsize[1], 'ylo yhi', sep=' ', file=fout_data)
    print('0.0', boxsize[2], 'zlo zhi', sep=' ', file=fout_data)
    print('\nMasses\n', file=fout_data)
    for i in range(len(typelist)):
        print(atomtype[typelist[i]][0],
              atomtype[typelist[i]][1],
              ' # ',
              typelist[i],
              sep=' ',
              file=fout_data)
    print('\nAtoms  # full\n', file=fout_data)
    for i in range(len(PISA_atom)):
        print(PISA_atom[i][0],
              PISA_atom[i][1],
              PISA_atom[i][2],
              PISA_atom[i][3],
              PISA_atom[i][4],
              PISA_atom[i][5],
              PISA_atom[i][6],
              sep=' ',
              file=fout_data)
    print('\nBonds\n', file=fout_data)
    for i in range(len(PISA_bond)):
        print(str(i + 1),
              PISA_bond[i][2],
              PISA_bond[i][0],
              PISA_bond[i][1],
              sep=' ',
              file=fout_data)
    fout_data.close()


def RAFTsteptwofile():

    # PISA system .itp file
    fin_itp = open('./gromacs/step3/PISA.itp', 'r')
    # Atom information
    Atom = []
    CTAlist = []
    molid_to_new = {}
    atomid_to_new = {}
    temps = []
    while True:
        line = fin_itp.readline()
        if line == '[ atoms ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            temps.append(lines)
            if lines[4] == 'm_b' or lines[4] == 'c_b':
                if lines[2] not in CTAlist:
                    CTAlist.append(lines[2])
    for temp in temps:
        if temp[2] in CTAlist:
            Atom.append(temp)
    for i in range(len(CTAlist)):
        molid_to_new[CTAlist[i]] = str(i + 1)
    for i in range(len(Atom)):
        atomid_to_new[Atom[i][0]] = str(i + 1)
        Atom[i][0], Atom[i][2], Atom[i][5] = str(i + 1), molid_to_new[
            Atom[i][2]], str(i + 1)

    # bond information
    bond = []
    while True:
        line = fin_itp.readline()
        if line == '[ bonds ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            if lines[0] in atomid_to_new.keys() and \
               lines[1] in atomid_to_new.keys():
                lines[0] = atomid_to_new[lines[0]]
                lines[1] = atomid_to_new[lines[1]]
                bond.append(lines)

    # angle information
    angle = []
    while True:
        line = fin_itp.readline()
        if line == '[ angles ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            if lines[0] in atomid_to_new.keys() and \
               lines[1] in atomid_to_new.keys() and \
               lines[2] in atomid_to_new.keys():
                lines[0] = atomid_to_new[lines[0]]
                lines[1] = atomid_to_new[lines[1]]
                lines[2] = atomid_to_new[lines[2]]
                angle.append(lines)

    # dihedral information
    dihedral = []
    while True:
        line = fin_itp.readline()
        if line == '[ dihedrals ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            if lines[0] in atomid_to_new.keys() and \
               lines[1] in atomid_to_new.keys() and \
               lines[2] in atomid_to_new.keys() and \
               lines[3] in atomid_to_new.keys():
                lines[0] = atomid_to_new[lines[0]]
                lines[1] = atomid_to_new[lines[1]]
                lines[2] = atomid_to_new[lines[2]]
                lines[3] = atomid_to_new[lines[3]]
                dihedral.append(lines)

    # improper information
    improper = []
    while True:
        line = fin_itp.readline()
        if line == '[ dihedrals ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            if lines[0] in atomid_to_new.keys() and \
               lines[1] in atomid_to_new.keys() and \
               lines[2] in atomid_to_new.keys() and \
               lines[3] in atomid_to_new.keys():
                lines[0] = atomid_to_new[lines[0]]
                lines[1] = atomid_to_new[lines[1]]
                lines[2] = atomid_to_new[lines[2]]
                lines[3] = atomid_to_new[lines[3]]
                improper.append(lines)

    Numatom = _Numatom()
    reactpoint = _rectpoint()
    MBcharge = _getcharge('MB', '1')
    Icharge = _getcharge('I', '2')
    MBatom, MBbond, MBangle, MBproper, MBpair, MBimproper = _getstr('MB')
    Iatom, Ibond, Iangle, Iproper, Ipair, Iimproper = _getstr('I')
    MB_matrix = _getconnectmatrix(Numatom[3], MBbond)
    _, nINI, _, nMB = _Npolymer()
    Ipos, _, _, _, MBpos = _getpos()

    # PISA itp file
    fout_itp = open('./topper/PISA.itp', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=fout_itp)
    print('\n[ moleculetype ]', file=fout_itp)
    print('PISA     3', file=fout_itp)

    print('\n[ atoms ]', file=fout_itp)
    for i in range(len(Atom)):
        print(Atom[i][0].rjust(6),
              Atom[i][1].rjust(6),
              Atom[i][2].rjust(6),
              Atom[i][3].rjust(6),
              Atom[i][4].rjust(6),
              Atom[i][5].rjust(6),
              Atom[i][6].rjust(14),
              Atom[i][7].rjust(12),
              file=fout_itp)

    for i in range(nINI):
        for j in range(len(Iatom)):
            line = Iatom[j]
            if j != reactpoint[0][0]:
                print(str(int(line[0]) + len(Atom) + i * Numatom[0]).rjust(6),
                      line[1].rjust(6),
                      str(int(Atom[-1][2]) + i + 1).rjust(6),
                      'INI'.rjust(6), (line[4] + line[0]).rjust(6),
                      str(int(line[5]) + len(Atom) + i * Numatom[0]).rjust(6),
                      Icharge[j].rjust(14),
                      line[7].rjust(12),
                      file=fout_itp)
            if j == reactpoint[0][0]:
                print(str(int(line[0]) + len(Atom) + i * Numatom[0]).rjust(6),
                      line[1].rjust(6),
                      str(int(Atom[-1][2]) + i + 1).rjust(6),
                      'INI'.rjust(6),
                      'i_p'.rjust(6),
                      str(int(line[5]) + len(Atom) + i * Numatom[0]).rjust(6),
                      Icharge[j].rjust(14),
                      line[7].rjust(12),
                      file=fout_itp)

    for i in range(nMB):
        for j in range(len(MBatom)):
            line = MBatom[j]
            if j != reactpoint[3][0] and j != reactpoint[3][1]:
                print(str(
                    int(line[0]) + len(Atom) + nINI * len(Iatom) +
                    i * Numatom[3]).rjust(6),
                      line[1].rjust(6),
                      str(int(Atom[-1][2]) + nINI + i + 1).rjust(6),
                      'MOB'.rjust(6), (line[4] + line[0]).rjust(6),
                      str(
                          int(line[5]) + len(Atom) + nINI * len(Iatom) +
                          i * Numatom[3]).rjust(6),
                      MBcharge[j].rjust(14),
                      line[7].rjust(12),
                      file=fout_itp)
            if j == reactpoint[3][0]:
                print(str(
                    int(line[0]) + len(Atom) + nINI * len(Iatom) +
                    i * Numatom[3]).rjust(6),
                      line[1].rjust(6),
                      str(int(Atom[-1][2]) + nINI + i + 1).rjust(6),
                      'MOB'.rjust(6),
                      'm_n'.rjust(6),
                      str(
                          int(line[5]) + len(Atom) + nINI * len(Iatom) +
                          i * Numatom[3]).rjust(6),
                      MBcharge[j].rjust(14),
                      line[7].rjust(12),
                      file=fout_itp)
            if j == reactpoint[3][1]:
                print(str(
                    int(line[0]) + len(Atom) + nINI * len(Iatom) +
                    i * Numatom[3]).rjust(6),
                      line[1].rjust(6),
                      str(int(Atom[-1][2]) + nINI + i + 1).rjust(6),
                      'MOB'.rjust(6),
                      'm_p'.rjust(6),
                      str(
                          int(line[5]) + len(Atom) + nINI * len(Iatom) +
                          i * Numatom[3]).rjust(6),
                      MBcharge[j].rjust(14),
                      line[7].rjust(12),
                      file=fout_itp)

    print('\n[ bonds ]', file=fout_itp)
    for i in range(len(bond)):
        if len(bond[i]) > 3:
            print(bond[i][0].rjust(5),
                  bond[i][1].rjust(6),
                  bond[i][2].rjust(6),
                  bond[i][3].rjust(12),
                  bond[i][4].rjust(14),
                  file=fout_itp)
        else:
            print(bond[i][0].rjust(5),
                  bond[i][1].rjust(6),
                  bond[i][2].rjust(6),
                  file=fout_itp)
    if Ibond != []:
        for i in range(nINI):
            for j in range(len(Ibond)):
                if 'mSeminario' in Ibond[j]:
                    lines = Ibond[j]
                    print(str(int(lines[0]) + len(Atom) +
                              i * Numatom[0]).rjust(5),
                          str(int(lines[1]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          lines[2].rjust(6),
                          lines[3].rjust(12),
                          str('{:.2f}'.format(float(lines[4]))).rjust(14),
                          file=fout_itp)
                else:
                    lines = Ibond[j]
                    print(str(int(lines[0]) + len(Atom) +
                              i * Numatom[0]).rjust(5),
                          str(int(lines[1]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          lines[2].rjust(6),
                          file=fout_itp)
    if MBbond != []:
        for i in range(nMB):
            for j in range(len(MBbond)):
                if 'mSeminario' in MBbond[j]:
                    lines = MBbond[j]
                    print(str(
                        int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                        i * Numatom[3]).rjust(5),
                          str(
                              int(lines[1]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          lines[2].rjust(6),
                          lines[3].rjust(12),
                          str('{:.2f}'.format(float(lines[4]))).rjust(14),
                          file=fout_itp)
                else:
                    lines = MBbond[j]
                    print(str(
                        int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                        i * Numatom[3]).rjust(5),
                          str(
                              int(lines[1]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          lines[2].rjust(6),
                          file=fout_itp)

    print('\n[ angles ]', file=fout_itp)
    for i in range(len(angle)):
        if len(angle[i]) > 4:
            print(angle[i][0].rjust(5),
                  angle[i][1].rjust(6),
                  angle[i][2].rjust(6),
                  angle[i][3].rjust(6),
                  angle[i][4].rjust(12),
                  angle[i][5].rjust(14),
                  file=fout_itp)
        else:
            print(angle[i][0].rjust(5),
                  angle[i][1].rjust(6),
                  angle[i][2].rjust(6),
                  angle[i][3].rjust(6),
                  file=fout_itp)
    if Iangle != []:
        for i in range(nINI):
            for j in range(len(Iangle)):
                if 'mSeminario' in Iangle[j]:
                    lines = Iangle[j]
                    print(str(int(lines[0]) + len(Atom) +
                              i * Numatom[0]).rjust(5),
                          str(int(lines[1]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          str(int(lines[2]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          lines[3].rjust(6),
                          lines[4].rjust(12),
                          str('{:.5f}'.format(float(lines[5]))).rjust(14),
                          file=fout_itp)
                else:
                    lines = Iangle[j]
                    print(str(int(lines[0]) + len(Atom) +
                              i * Numatom[0]).rjust(5),
                          str(int(lines[1]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          str(int(lines[2]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          lines[3].rjust(6),
                          file=fout_itp)
    if MBangle != []:
        for i in range(nMB):
            for j in range(len(MBangle)):
                if 'mSeminario' in MBangle[j]:
                    lines = MBangle[j]
                    print(str(
                        int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                        i * Numatom[3]).rjust(5),
                          str(
                              int(lines[1]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          str(
                              int(lines[2]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          lines[3].rjust(6),
                          lines[4].rjust(12),
                          str('{:.5f}'.format(float(lines[5]))).rjust(14),
                          file=fout_itp)
                else:
                    lines = MBangle[j]
                    print(str(
                        int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                        i * Numatom[3]).rjust(5),
                          str(
                              int(lines[1]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          str(
                              int(lines[2]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          lines[3].rjust(6),
                          file=fout_itp)

    print('\n[ dihedrals ]', file=fout_itp)
    for i in range(len(dihedral)):
        if len(dihedral[i]) > 5:
            print(dihedral[i][0].rjust(5),
                  dihedral[i][1].rjust(6),
                  dihedral[i][2].rjust(6),
                  dihedral[i][3].rjust(6),
                  dihedral[i][4].rjust(6),
                  dihedral[i][5].rjust(12),
                  dihedral[i][6].rjust(14),
                  file=fout_itp)
        else:
            print(dihedral[i][0].rjust(5),
                  dihedral[i][1].rjust(6),
                  dihedral[i][2].rjust(6),
                  dihedral[i][3].rjust(6),
                  dihedral[i][4].rjust(6),
                  file=fout_itp)
    if Iproper != []:
        for i in range(nINI):
            for j in range(len(Iproper)):
                if 'mSeminario' in Iproper[j]:
                    lines = Iproper[j]
                    print(str(int(lines[0]) + len(Atom) +
                              i * Numatom[0]).rjust(5),
                          str(int(lines[1]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          str(int(lines[2]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          str(int(lines[3]) + len(Atom) +
                              i * Numatom[0]).rjust(6),
                          lines[4].rjust(6),
                          lines[5].rjust(12),
                          str('{:.5f}'.format(float(lines[6]))).rjust(14),
                          file=fout_itp)
                else:
                    if j != 0:
                        line1 = Iproper[j - 1]
                        line2 = Iproper[j]
                        lines = line2
                        if line1[0] != line2[0] or line1[1] != line2[
                                1] or line1[2] != line2[2] or line1[
                                    3] != line2[3] or line1[4] != line2[4]:
                            print(str(
                                int(lines[0]) + len(Atom) +
                                i * Numatom[0]).rjust(5),
                                  str(
                                      int(lines[1]) + len(Atom) +
                                      i * Numatom[0]).rjust(6),
                                  str(
                                      int(lines[2]) + len(Atom) +
                                      i * Numatom[0]).rjust(6),
                                  str(
                                      int(lines[3]) + len(Atom) +
                                      i * Numatom[0]).rjust(6),
                                  lines[4].rjust(6),
                                  file=fout_itp)
                    else:
                        lines = Iproper[j]
                        print(str(int(lines[0]) + len(Atom) +
                                  i * Numatom[0]).rjust(5),
                              str(int(lines[1]) + len(Atom) +
                                  i * Numatom[0]).rjust(6),
                              str(int(lines[2]) + len(Atom) +
                                  i * Numatom[0]).rjust(6),
                              str(int(lines[3]) + len(Atom) +
                                  i * Numatom[0]).rjust(6),
                              lines[4].rjust(6),
                              file=fout_itp)
    if MBproper != []:
        for i in range(nMB):
            for j in range(len(MBproper)):
                if 'mSeminario' in MBproper[j]:
                    lines = MBproper[j]
                    print(str(
                        int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                        i * Numatom[3]).rjust(5),
                          str(
                              int(lines[1]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          str(
                              int(lines[2]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          str(
                              int(lines[3]) + len(Atom) + nINI * len(Iatom) +
                              i * Numatom[3]).rjust(6),
                          lines[4].rjust(6),
                          lines[5].rjust(12),
                          str('{:.5f}'.format(float(lines[6]))).rjust(14),
                          file=fout_itp)
                else:
                    if j != 0:
                        line1 = MBproper[j - 1]
                        line2 = MBproper[j]
                        lines = line2
                        if line1[0] != line2[0] or line1[1] != line2[
                                1] or line1[2] != line2[2] or line1[
                                    3] != line2[3] or line1[4] != line2[4]:
                            print(str(
                                int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                                i * Numatom[3]).rjust(5),
                                  str(
                                      int(lines[1]) + len(Atom) +
                                      nINI * len(Iatom) +
                                      i * Numatom[3]).rjust(6),
                                  str(
                                      int(lines[2]) + len(Atom) +
                                      nINI * len(Iatom) +
                                      i * Numatom[3]).rjust(6),
                                  str(
                                      int(lines[3]) + len(Atom) +
                                      nINI * len(Iatom) +
                                      i * Numatom[3]).rjust(6),
                                  lines[4].rjust(6),
                                  file=fout_itp)
                    else:
                        lines = MBproper[j]
                        print(str(
                            int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                            i * Numatom[3]).rjust(5),
                              str(
                                  int(lines[1]) + len(Atom) +
                                  nINI * len(Iatom) + i * Numatom[3]).rjust(6),
                              str(
                                  int(lines[2]) + len(Atom) +
                                  nINI * len(Iatom) + i * Numatom[3]).rjust(6),
                              str(
                                  int(lines[3]) + len(Atom) +
                                  nINI * len(Iatom) + i * Numatom[3]).rjust(6),
                              lines[4].rjust(6),
                              file=fout_itp)

    print('\n[ pairs ]', file=fout_itp)
    for i in range(len(dihedral)):
        print(dihedral[i][0].rjust(5),
              dihedral[i][3].rjust(6),
              '1'.rjust(6),
              file=fout_itp)
    if Ipair != []:
        for i in range(nINI):
            for j in range(len(Ipair)):
                lines = Ipair[j]
                print(str(int(lines[0]) + len(Atom) + i * Numatom[0]).rjust(5),
                      str(int(lines[1]) + len(Atom) + i * Numatom[0]).rjust(6),
                      lines[2].rjust(6),
                      file=fout_itp)
    if MBpair != []:
        for i in range(nMB):
            for j in range(len(MBpair)):
                lines = MBpair[j]
                print(str(
                    int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                    i * Numatom[3]).rjust(5),
                      str(
                          int(lines[1]) + len(Atom) + nINI * len(Iatom) +
                          i * Numatom[3]).rjust(6),
                      lines[2].rjust(6),
                      file=fout_itp)

    print('\n[ dihedrals ]', file=fout_itp)
    for i in range(len(improper)):
        print(improper[i][0].rjust(5),
              improper[i][1].rjust(6),
              improper[i][2].rjust(6),
              improper[i][3].rjust(6),
              improper[i][4].rjust(6),
              file=fout_itp)
    if Iimproper != []:
        for i in range(nINI):
            for j in range(len(Iimproper)):
                lines = Iimproper[j]
                print(str(int(lines[0]) + len(Atom) + i * Numatom[0]).rjust(5),
                      str(int(lines[1]) + len(Atom) + i * Numatom[0]).rjust(6),
                      str(int(lines[2]) + len(Atom) + i * Numatom[0]).rjust(6),
                      str(int(lines[3]) + len(Atom) + i * Numatom[0]).rjust(6),
                      lines[4].rjust(6),
                      file=fout_itp)
    if MBimproper != []:
        for i in range(nMB):
            for j in range(len(MBimproper)):
                lines = MBimproper[j]
                print(str(
                    int(lines[0]) + len(Atom) + nINI * len(Iatom) +
                    i * Numatom[3]).rjust(5),
                      str(
                          int(lines[1]) + len(Atom) + nINI * len(Iatom) +
                          i * Numatom[3]).rjust(6),
                      str(
                          int(lines[2]) + len(Atom) + nINI * len(Iatom) +
                          i * Numatom[3]).rjust(6),
                      str(
                          int(lines[3]) + len(Atom) + nINI * len(Iatom) +
                          i * Numatom[3]).rjust(6),
                      lines[4].rjust(6),
                      file=fout_itp)
    fout_itp.close()
    f_sys.close()

    # MB pdb file
    f_MBpdb = open('./topper/MB.pdb', 'w')
    print('REMARK  PISAMD PDB file', file=f_MBpdb)
    for i in range(len(MBpos)):
        lines = MBpos[i].split()
        if i != reactpoint[2][0] and i != reactpoint[2][1]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ', (lines[0] + str(i + 1)).ljust(4, ' '),
                  ' ',
                  'MOB'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_MBpdb)
        if i == reactpoint[2][0]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ',
                  'm_n'.ljust(4, ' '),
                  ' ',
                  'MOB'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_MBpdb)
        if i == reactpoint[2][1]:
            print('ATOM',
                  '  ',
                  str(i + 1).rjust(5, ' '),
                  ' ',
                  'm_p'.ljust(4, ' '),
                  ' ',
                  'MOB'.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  str('{:.3f}'.format(float(lines[1]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[2]))).rjust(8, ' '),
                  str('{:.3f}'.format(float(lines[3]))).rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  lines[0].rjust(2, ' '),
                  sep='',
                  file=f_MBpdb)
    print('TER'.ljust(4, ' '),
          '  ',
          str(len(MBpos) + 1).rjust(5, ' '),
          sep='',
          file=f_MBpdb)
    if MBbond != []:
        for i in range(len(MB_matrix)):
            print('CONECT',
                  str(i + 1).rjust(5, ' '),
                  sep='',
                  end='',
                  file=f_MBpdb)
            for j in range(len(MB_matrix[i])):
                if MB_matrix[i][j] == 1:
                    print(str(j + 1).rjust(5, ' '),
                          sep='',
                          end='',
                          file=f_MBpdb)
            print(file=f_MBpdb)
    print('END', file=f_MBpdb)

    # CTA gro file and transform to pdb
    atom = []
    files = glob.glob('./gromacs/step3/[0-9]*.gro')
    fin_gro = open(files[0], 'r')
    fout_gro = open('./topper/macroCTA.gro', 'w')

    line = fin_gro.readline()
    Natom = fin_gro.readline()
    cout = 0
    while True:
        line = fin_gro.readline()
        molname = line[5:10].strip()
        molid = line[0:5].strip()
        if molid in molid_to_new.keys() and (molname == 'CTA' or molname
                                             == 'INI' or molname == 'MOA'):
            lines = line.split()
            atom.append(lines)
        cout += 1
        if cout == int(Natom):
            line = fin_gro.readline()
            atom.append(line)
            break
    boxsize = atom[-1].split()
    print('System', file=fout_gro)
    print(' ' + str(len(atom) - 1), file=fout_gro)
    for i in range(len(atom) - 1):
        x = float(atom[i][-6])
        y = float(atom[i][-5])
        z = float(atom[i][-4])
        if x > float(boxsize[0]):
            x = round(x - float(boxsize[0]), 3)
        if x < 0:
            x = round(x + float(boxsize[0]), 3)
        if y > float(boxsize[1]):
            y = round(y - float(boxsize[1]), 3)
        if y < 0:
            y = round(y + float(boxsize[1]), 3)
        if z > float(boxsize[2]):
            z = round(z - float(boxsize[2]), 3)
        if z < 0:
            z = round(z + float(boxsize[2]), 3)
        atom[i][-6], atom[i][-5], atom[i][-4] = str(('%.3f' % x)), str(
            ('%.3f' % y)), str(('%.3f' % z))
        atomid = int(Atom[i][0])
        if atomid > 99999:
            atomid = atomid - 100000
        print(Atom[i][2].rjust(5),
              Atom[i][3].ljust(5),
              Atom[i][4].rjust(5),
              str(atomid).rjust(5),
              atom[i][-6].rjust(8),
              atom[i][-5].rjust(8),
              atom[i][-4].rjust(8),
              sep='',
              file=fout_gro)
    print(atom[-1], file=fout_gro)
    fin_gro.close()
    fout_gro.close()
    flow = subprocess.Popen('gmx editconf -f macroCTA.gro -o macroCTA.pdb',
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=True,
                            cwd='./topper')
    flow.communicate()
    if flow.returncode != 0:
        return False
    else:
        return True


def PISA_changeitp(Ilist, Rlist, Clist, MAlist, MBlist, reac):

    fin_itp = open('./gromacs/step6/PISA.itp', 'r')
    atom = []
    breakbond = []
    atomtype = {}  # aid to atype
    reac_molname = [['MOB', 'MOB'], ['MOA', 'MOB'], ['CTA', 'MOB'],
                    ['INI', 'MOB'], ['CTA', 'MOA']]
    reac_list = [['r_p', 'm_n'], ['r_p', 'c_n'], ['i_p', 'm_n'],
                 ['i_p', 'c_n']]

    # atom: atom_id atom_type mol_id mol_name atom_name index charge mass (str)
    while True:
        line = fin_itp.readline()
        if line == '[ atoms ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            lines = line.split()
            atom.append(lines)
            atomtype[lines[0]] = lines[1]

    # adjacent
    adjacent = nx.Graph()
    while True:
        line = fin_itp.readline()
        if line == '[ bonds ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            adjacent.add_edge(int(lines[0]) - 1, int(lines[1]) - 1)
    fin_itp.close()

    # change atom & adjacent
    reactemp = []
    for i in range(len(reac)):
        reac_group = []
        molname_group = []
        reacindex = None
        molnameindex = None
        reac_group.append(reac[i][4])
        reac_group.append(reac[i][5])
        reac_group.sort()
        molname_group.append(atom[int(reac[i][0]) - 1][3])
        molname_group.append(atom[int(reac[i][1]) - 1][3])
        molname_group.sort()
        for j in range(4):
            temp = reac_list[j]
            temp.sort()
            if reac_group == temp:
                reacindex = j
                break
        for j in range(5):
            temp = reac_molname[j]
            temp.sort()
            if molname_group == temp:
                molnameindex = j
                break
        if molnameindex is not None:
            reactemp.append(reac[i])

        mola_id = None
        molb_id = None
        molc_id = None
        mola_range = [0, 0]
        molb_range = [0, 0]
        molc_range = [0, 0]
        # MB: state1-->state3
        if molnameindex == 0 and (reacindex == 0 or reacindex == 2):
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'm_n':
                mola_id = int(reac[i][3])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            if reac[i][4] == 'm_n':
                mola_id = int(reac[i][2])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            for k in range(len(MBlist)):
                if MBlist[k][0] == mola_id:
                    mola_range[0] = MBlist[k][1]
                    mola_range[1] = MBlist[k][2]
                    break
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'm_p':
                    if reacindex == 0:
                        atom[k][4] = 'r_p'
                    if reacindex == 2:
                        atom[k][4] = 'i_p'
                    break
            _PISAchangeMol(atom, mola_range, 'MB', 3)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]

    # MB: state1-->state2 & MA: state3-->state5
        if molnameindex == 2 and reacindex == 0:
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'm_n':
                mola_id = int(reac[i][3])
                molb_id = int(reac[i][2])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            if reac[i][4] == 'm_n':
                mola_id = int(reac[i][2])
                molb_id = int(reac[i][3])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            for k in range(len(MBlist)):
                if MBlist[k][0] == mola_id:
                    mola_range[0] = MBlist[k][1]
                    mola_range[1] = MBlist[k][2]
                    break
            for k in range(len(Rlist)):
                if Rlist[k][0] == molb_id:
                    molb_range[0] = Rlist[k][1]
                    molb_range[1] = Rlist[k][2]
                    break
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'm_p':
                    atom[k][4] = 'r_p'
                    break
            _PISAchangeMol(atom, mola_range, 'MB', 2)
            _PISAchangeMol(atom, molb_range, 'R', 1)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molb_range[0] - 1, molb_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
        if molnameindex == 3 and reacindex == 2:
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'm_n':
                mola_id = int(reac[i][3])
                molb_id = int(reac[i][2])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            if reac[i][4] == 'm_n':
                mola_id = int(reac[i][2])
                molb_id = int(reac[i][3])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            for k in range(len(MBlist)):
                if MBlist[k][0] == mola_id:
                    mola_range[0] = MBlist[k][1]
                    mola_range[1] = MBlist[k][2]
                    break
            for k in range(len(Ilist)):
                if Ilist[k][0] == molb_id:
                    molb_range[0] = Ilist[k][1]
                    molb_range[1] = Ilist[k][2]
                    break
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'm_p':
                    atom[k][4] = 'i_p'
                    break
            _PISAchangeMol(atom, mola_range, 'MB', 2)
            _PISAchangeMol(atom, molb_range, 'I', 1)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molb_range[0] - 1, molb_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
        if molnameindex == 1 and (reacindex == 0 or reacindex == 2):
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'm_n':
                mola_id = int(reac[i][3])
                molb_id = int(reac[i][2])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            if reac[i][4] == 'm_n':
                mola_id = int(reac[i][2])
                molb_id = int(reac[i][3])
                atom[int(reac[i][1]) - 1][4] = 'm_b'
                atom[int(reac[i][0]) - 1][4] = 'm_b'
            for k in range(len(MBlist)):
                if MBlist[k][0] == mola_id:
                    mola_range[0] = MBlist[k][1]
                    mola_range[1] = MBlist[k][2]
                    break
            for k in range(len(MAlist)):
                if MAlist[k][0] == molb_id:
                    molb_range[0] = MAlist[k][1]
                    molb_range[1] = MAlist[k][2]
                    break
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'm_p':
                    if reacindex == 0:
                        atom[k][4] = 'r_p'
                    if reacindex == 2:
                        atom[k][4] = 'i_p'
                    break
            _PISAchangeMol(atom, mola_range, 'MB', 3)
            _PISAchangeMol(atom, molb_range, 'MA', 5)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molb_range[0] - 1, molb_range[1]):
                atomtype[atom[k][0]] = atom[k][1]

    # MB: state3<-->state4
        if (molnameindex == 2 and reacindex == 1) or (molnameindex == 2
                                                      and reacindex == 3):
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'c_n':
                mola_id = int(reac[i][3])
                molb_id = int(reac[i][2])
                if reacindex == 1:
                    atom[int(reac[i][0]) - 1][4] = 'r_b'
                if reacindex == 3:
                    atom[int(reac[i][0]) - 1][4] = 'i_b'
            if reac[i][4] == 'c_n':
                mola_id = int(reac[i][2])
                molb_id = int(reac[i][3])
                if reacindex == 1:
                    atom[int(reac[i][1]) - 1][4] = 'r_b'
                if reacindex == 3:
                    atom[int(reac[i][1]) - 1][4] = 'i_b'
            for k in range(len(Clist)):
                if Clist[k][0] == mola_id:
                    mola_range[0] = Clist[k][1]
                    mola_range[1] = Clist[k][2]
                    break
            for k in range(len(MBlist)):
                if MBlist[k][0] == molb_id:
                    molb_range[0] = MBlist[k][1]
                    molb_range[1] = MBlist[k][2]
                    break
            atoma_id = None
            atomb_id = None
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'c_b':
                    atoma_id = k
            for k in range(len(atom)):
                if adjacent.has_edge(atoma_id, k) and atom[k][4] == 'r_b':
                    atomb_id = k
                    molc_id = int(atom[k][2])
                    atom[k][4] = 'r_p'
                if adjacent.has_edge(atoma_id, k) and atom[k][4] == 'i_b':
                    atomb_id = k
                    molc_id = int(atom[k][2])
                    atom[k][4] = 'i_p'
            adjacent.remove_edge(atoma_id, atomb_id)
            if atoma_id < atomb_id:
                temp = [''] * 2
                temp[0] = str(atoma_id + 1)
                temp[1] = str(atomb_id + 1)
                breakbond.append(temp)
            if atomb_id < atoma_id:
                temp = [''] * 2
                temp[0] = str(atomb_id + 1)
                temp[1] = str(atoma_id + 1)
                breakbond.append(temp)
            _PISAchangeMol(atom, mola_range, 'C', 3)
            _PISAchangeMol(atom, molb_range, 'MB', 4)
            for k in range(len(MAlist)):
                if MAlist[k][0] == molc_id:
                    molc_range[0] = MAlist[k][1]
                    molc_range[1] = MAlist[k][2]
                    break
            if molc_range[0] != 0 and molc_range[1] != 0:
                _PISAchangeMol(atom, molc_range, 'MA', 3)
            else:
                for k in range(len(MBlist)):
                    if MBlist[k][0] == molc_id:
                        molc_range[0] = MBlist[k][1]
                        molc_range[1] = MBlist[k][2]
                        break
                if molc_range[0] != 0 and molc_range[1] != 0:
                    _PISAchangeMol(atom, molc_range, 'MB', 3)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molb_range[0] - 1, molb_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molc_range[0] - 1, molc_range[1]):
                atomtype[atom[k][0]] = atom[k][1]

    # MA: state3<-->state4
        if (molnameindex == 4 and reacindex == 1) or (molnameindex == 4
                                                      and reacindex == 3):
            adjacent.add_edge(int(reac[i][0]) - 1, int(reac[i][1]) - 1)
            if reac[i][5] == 'c_n':
                mola_id = int(reac[i][3])
                molb_id = int(reac[i][2])
                if reacindex == 1:
                    atom[int(reac[i][0]) - 1][4] = 'r_b'
                if reacindex == 3:
                    atom[int(reac[i][0]) - 1][4] = 'i_b'
            if reac[i][4] == 'c_n':
                mola_id = int(reac[i][2])
                molb_id = int(reac[i][3])
                if reacindex == 1:
                    atom[int(reac[i][1]) - 1][4] = 'r_b'
                if reacindex == 3:
                    atom[int(reac[i][1]) - 1][4] = 'i_b'
            for k in range(len(Clist)):
                if Clist[k][0] == mola_id:
                    mola_range[0] = Clist[k][1]
                    mola_range[1] = Clist[k][2]
                    break
            for k in range(len(MAlist)):
                if MAlist[k][0] == molb_id:
                    molb_range[0] = MAlist[k][1]
                    molb_range[1] = MAlist[k][2]
                    break
            atoma_id = None
            atomb_id = None
            for k in range(mola_range[0] - 1, mola_range[1]):
                if atom[k][4] == 'c_b':
                    atoma_id = k
            for k in range(len(atom)):
                if adjacent.has_edge(atoma_id, k) and atom[k][4] == 'r_b':
                    atomb_id = k
                    molc_id = int(atom[k][2])
                    atom[k][4] = 'r_p'
                if adjacent.has_edge(atoma_id, k) and atom[k][4] == 'i_b':
                    atomb_id = k
                    molc_id = int(atom[k][2])
                    atom[k][4] = 'i_p'
            adjacent.remove_edge(atoma_id, atomb_id)
            if atoma_id < atomb_id:
                temp = [''] * 2
                temp[0] = str(atoma_id + 1)
                temp[1] = str(atomb_id + 1)
                breakbond.append(temp)
            if atomb_id < atoma_id:
                temp = [''] * 2
                temp[0] = str(atomb_id + 1)
                temp[1] = str(atoma_id + 1)
                breakbond.append(temp)
            _PISAchangeMol(atom, mola_range, 'C', 2)
            _PISAchangeMol(atom, molb_range, 'MA', 4)
            for k in range(len(MAlist)):
                if MAlist[k][0] == molc_id:
                    molc_range[0] = MAlist[k][1]
                    molc_range[1] = MAlist[k][2]
                    break
            if molc_range[0] != 0 and molc_range[1] != 0:
                _PISAchangeMol(atom, molc_range, 'MA', 3)
            else:
                for k in range(len(MBlist)):
                    if MBlist[k][0] == molc_id:
                        molc_range[0] = MBlist[k][1]
                        molc_range[1] = MBlist[k][2]
                        break
                if molc_range[0] != 0 and molc_range[1] != 0:
                    _PISAchangeMol(atom, molc_range, 'MB', 3)
            for k in range(mola_range[0] - 1, mola_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molb_range[0] - 1, molb_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
            for k in range(molc_range[0] - 1, molc_range[1]):
                atomtype[atom[k][0]] = atom[k][1]
    # change ffbond
    reac = reactemp
    bond, angle, dihedral, improper = _PISAchangeffbond(
        len(atom), adjacent, reac, atomtype, breakbond)

    # write PISA.itp
    fout_itp = open('./gromacs/step6/PISA.itp', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=fout_itp)
    print('\n[ moleculetype ]', file=fout_itp)
    print('PISA     3', file=fout_itp)

    print('\n[ atoms ]', file=fout_itp)
    for i in range(len(atom)):
        print(atom[i][0].rjust(6),
              atom[i][1].rjust(6),
              atom[i][2].rjust(6),
              atom[i][3].rjust(6),
              atom[i][4].rjust(6),
              atom[i][5].rjust(6),
              atom[i][6].rjust(14),
              atom[i][7].rjust(12),
              file=fout_itp)
    print('\n[ bonds ]', file=fout_itp)
    for i in range(len(bond)):
        if len(bond[i]) > 3:
            print(bond[i][0].rjust(5),
                  bond[i][1].rjust(6),
                  bond[i][2].rjust(6),
                  bond[i][3].rjust(12),
                  bond[i][4].rjust(14),
                  file=fout_itp)
        else:
            print(bond[i][0].rjust(5),
                  bond[i][1].rjust(6),
                  bond[i][2].rjust(6),
                  file=fout_itp)
    print('\n[ angles ]', file=fout_itp)
    for i in range(len(angle)):
        if len(angle[i]) > 4:
            print(angle[i][0].rjust(5),
                  angle[i][1].rjust(6),
                  angle[i][2].rjust(6),
                  angle[i][3].rjust(6),
                  angle[i][4].rjust(12),
                  angle[i][5].rjust(14),
                  file=fout_itp)
        else:
            print(angle[i][0].rjust(5),
                  angle[i][1].rjust(6),
                  angle[i][2].rjust(6),
                  angle[i][3].rjust(6),
                  file=fout_itp)
    print('\n[ dihedrals ]', file=fout_itp)
    for i in range(len(dihedral)):
        if len(dihedral[i]) > 5:
            print(dihedral[i][0].rjust(5),
                  dihedral[i][1].rjust(6),
                  dihedral[i][2].rjust(6),
                  dihedral[i][3].rjust(6),
                  dihedral[i][4].rjust(6),
                  dihedral[i][5].rjust(12),
                  dihedral[i][6].rjust(14),
                  file=fout_itp)
        else:
            print(dihedral[i][0].rjust(5),
                  dihedral[i][1].rjust(6),
                  dihedral[i][2].rjust(6),
                  dihedral[i][3].rjust(6),
                  dihedral[i][4].rjust(6),
                  file=fout_itp)
    print('\n[ pairs ]', file=fout_itp)
    for i in range(len(dihedral)):
        print(dihedral[i][0].rjust(5),
              dihedral[i][3].rjust(6),
              '1'.rjust(6),
              file=fout_itp)
    print('\n[ dihedrals ]', file=fout_itp)
    for i in range(len(improper)):
        print(improper[i][0].rjust(5),
              improper[i][1].rjust(6),
              improper[i][2].rjust(6),
              improper[i][3].rjust(6),
              improper[i][4].rjust(6),
              file=fout_itp)
    fout_itp.close()
    f_sys.close()
    return atom


def PISA_changegro(frame, atom):

    system_name = ''
    Natom = ''
    boxsize = ''
    old_atom = []
    fin_gro = open('./gromacs/step6/' + str(frame) + '.gro', 'r')
    system_name = fin_gro.readline()
    Natom = fin_gro.readline()
    for i in range(int(Natom)):
        line = fin_gro.readline()
        temp = [''] * 10
        temp[0] = line[0:5].strip()
        temp[1] = line[5:10].strip()
        temp[2] = line[10:15].strip()
        temp[3] = line[15:20].strip()
        temp[4] = line[20:28].strip()
        temp[5] = line[28:36].strip()
        temp[6] = line[36:44].strip()
        temp[7] = line[44:52].strip()
        temp[8] = line[52:60].strip()
        temp[9] = line[60:68].strip()
        old_atom.append(temp)
    boxsize = fin_gro.readline()
    fin_gro.close()

    for i in range(len(atom)):
        old_atom[i][2] = atom[i][4]

    fout_gro = open('./gromacs/step6/' + str(frame) + '.gro', 'w')
    print(system_name, end='', file=fout_gro)
    print(Natom, end='', file=fout_gro)
    for i in range(len(old_atom)):
        print(old_atom[i][0].rjust(5),
              old_atom[i][1].ljust(5),
              old_atom[i][2].rjust(5),
              old_atom[i][3].rjust(5),
              old_atom[i][4].rjust(8),
              old_atom[i][5].rjust(8),
              old_atom[i][6].rjust(8),
              old_atom[i][7].rjust(8),
              old_atom[i][8].rjust(8),
              old_atom[i][9].rjust(8),
              sep='',
              file=fout_gro)
    print(boxsize, file=fout_gro)
    fout_gro.close()


def PISA_grotodata(frame):

    PISA_atom = []
    PISA_bond = []
    typelist = []
    boxsize = [0] * 3
    Natom = 0
    Nmol = 0
    Nbondtype = 0

    # atomtype to [index, mass]
    atomtype = {}
    fin_itp = open('./topper/ffnonbonded.itp', 'r')
    while True:
        line = fin_itp.readline()
        if line == '[ atomtypes ]\n':
            line = fin_itp.readline()
            break
    index = 1
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            temp = [''] * 2
            lines = line.split()
            temp[0] = str(index)
            temp[1] = lines[2]
            atomtype[lines[0]] = temp
            index += 1
            typelist.append(lines[0])
    fin_itp.close()

    # polymer atom and bond
    fin_itp = open('./gromacs/step6/PISA.itp', 'r')
    while True:
        line = fin_itp.readline()
        if line == '[ atoms ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            temp = [''] * 7  # atomid molid atomtype charge X Y Z
            lines = line.split()
            temp[0] = lines[0]
            temp[1] = lines[2]
            temp[2] = atomtype[lines[1]][0]
            temp[3] = lines[6]
            PISA_atom.append(temp)
    Natom = Natom + len(PISA_atom)
    Nmol = int(PISA_atom[-1][1])
    while True:
        line = fin_itp.readline()
        if line == '[ bonds ]\n':
            break
    while True:
        line = fin_itp.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            temp = [''] * 3  # atomaid atombid bondtype
            lines = line.split()
            temp[0] = lines[0]
            temp[1] = lines[1]
            temp[2] = '1'
            PISA_bond.append(temp)
    fin_itp.close()
    Nbondtype += 1

    # Ion/solvent atom and bond
    mol_list = []
    fin_top = open('./gromacs/step6/system.top', 'r')
    while True:
        line = fin_top.readline()
        if line == '[ molecules ]\n':
            break
    while True:
        line = fin_top.readline()
        if line == '\n' or line == ' \n' or line == '':
            break
        else:
            lines = line.split()
            mol_list.append(lines)
    fin_top.close()
    for i in range(1, len(mol_list)):
        atom_temp = []
        bond_temp = []
        fin_itp = open('./topper/' + mol_list[i][0] + '.itp', 'r')
        while True:
            line = fin_itp.readline()
            if line == '[ atoms ]\n':
                break
        while True:
            line = fin_itp.readline()
            if line == '\n' or line == ' \n' or line == '':
                break
            else:
                lines = line.split()
                atom_temp.append(lines)
        while True:
            if line == '':
                break
            line = fin_itp.readline()
            if line == '[ bonds ]\n':
                break
        while True:
            line = fin_itp.readline()
            if line == '\n' or line == ' \n' or line == '':
                break
            else:
                lines = line.split()
                bond_temp.append(lines)
        for j in range(int(mol_list[i][1])):
            for k in range(len(atom_temp)):
                temp = [''] * 7  # atomid molid atomtype charge X Y Z
                temp[0] = int(atom_temp[k][0]) + Natom + j * len(atom_temp)
                temp[1] = int(atom_temp[k][2]) + Nmol + j
                temp[2] = atomtype[atom_temp[k][1]][0]
                temp[3] = atom_temp[k][6]
                PISA_atom.append(temp)
        if bond_temp != []:
            Nbondtype += 1
            for j in range(int(mol_list[i][1])):
                for k in range(len(bond_temp)):
                    temp = [''] * 3  # atomaid atombid bondtype
                    temp[0] = int(bond_temp[k][0]) + Natom + j * len(atom_temp)
                    temp[1] = int(bond_temp[k][1]) + Natom + j * len(atom_temp)
                    temp[2] = Nbondtype
                    PISA_bond.append(temp)
        Natom = Natom + int(mol_list[i][1]) * len(atom_temp)
        Nmol = Nmol + int(mol_list[i][1])
        fin_top.close()
    # atom position
    fin_gro = open('./gromacs/step6/' + str(frame) + '.gro')
    line = fin_gro.readline()
    line = fin_gro.readline()
    Total_atom = int(line)
    for i in range(Total_atom):
        line = fin_gro.readline()
        lines = line.split()
        PISA_atom[i][4] = lines[-4]
        PISA_atom[i][5] = lines[-5]
        PISA_atom[i][6] = lines[-6]
    line = fin_gro.readline()
    lines = line.split()
    boxsize[0] = lines[0]
    boxsize[1] = lines[1]
    boxsize[2] = lines[2]
    fin_gro.close()

    folder = os.path.exists('./gromacs/step6/traj')
    if not folder:
        os.makedirs('./gromacs/step6/traj')

    fout_data = open('./gromacs/step6/traj/' + str(frame) + '.data', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print('# LAMMPS data file written by PISAMD version ' + version,
          file=fout_data)
    print(len(PISA_atom), 'atoms', sep=' ', file=fout_data)
    print(len(PISA_bond), 'bonds', sep=' ', file=fout_data)
    print(len(typelist), 'atom types', sep=' ', file=fout_data)
    print(Nbondtype, 'bond types', sep=' ', file=fout_data)
    print('0.0', boxsize[0], 'xlo xhi', sep=' ', file=fout_data)
    print('0.0', boxsize[1], 'ylo yhi', sep=' ', file=fout_data)
    print('0.0', boxsize[2], 'zlo zhi', sep=' ', file=fout_data)
    print('\nMasses\n', file=fout_data)
    for i in range(len(typelist)):
        print(atomtype[typelist[i]][0],
              atomtype[typelist[i]][1],
              ' # ',
              typelist[i],
              sep=' ',
              file=fout_data)
    print('\nAtoms  # full\n', file=fout_data)
    for i in range(len(PISA_atom)):
        print(PISA_atom[i][0],
              PISA_atom[i][1],
              PISA_atom[i][2],
              PISA_atom[i][3],
              PISA_atom[i][4],
              PISA_atom[i][5],
              PISA_atom[i][6],
              sep=' ',
              file=fout_data)
    print('\nBonds\n', file=fout_data)
    for i in range(len(PISA_bond)):
        print(str(i + 1),
              PISA_bond[i][2],
              PISA_bond[i][0],
              PISA_bond[i][1],
              sep=' ',
              file=fout_data)
    fout_data.close()
