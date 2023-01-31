import json
import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem


def solvemol():
    Mol = []
    model = ['I.mol', 'R.mol', 'MA.mol', 'MB.mol', 'C.mol']
    for i in range(len(model)):
        mol = Chem.MolFromMolFile('./model/' + model[i])
        mol = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        Mol.append(mol)
    return Mol


def _Numatom():
    Numatom = []
    Mols = solvemol()
    for mol in Mols:
        temp = mol.GetNumAtoms()
        Numatom.append(temp)
    return Numatom


def _loadstr(f_in, begin, end):

    # Load Atom
    index_to_atomtype = {}
    T_type = []
    T_bond = []
    T_angle = []
    T_proper = []

    while True:
        line = f_in.readline()
        if line == '[ atoms ]\n':
            line = f_in.readline()
            break
    while True:
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            atoms = line.split()
            if int(atoms[0]) >= begin and int(atoms[0]) <= end:
                T_type.append(atoms[1])
            index_to_atomtype[atoms[0]] = atoms[1]

    # Load bond
    while True:
        line = f_in.readline()
        if line == '[ bonds ]\n':
            line = f_in.readline()
            break
    while True:
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            temp = line.split(';')
            m = re.search('mSeminario', temp[1])
            if m:
                bonds = temp[0].split()
                if (int(bonds[0]) >= begin and int(bonds[0]) <= end) or (
                        int(bonds[1]) >= begin and int(bonds[1]) <= end):
                    bonds[0] = index_to_atomtype[bonds[0]]
                    bonds[1] = index_to_atomtype[bonds[1]]
                    bonds[-1] = str('{:.2f}'.format(float(bonds[-1])))
                    T_bond.append(bonds)

    # Load angle
    while True:
        line = f_in.readline()
        if line == '[ angles ]\n':
            line = f_in.readline()
            break
    while True:
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            temp = line.split(';')
            m = re.search('mSeminario', temp[1])
            if m:
                angles = temp[0].split()
                if (int(angles[0]) >= begin and int(angles[0]) <= end) or (
                        int(angles[1]) >= begin and int(angles[1]) <= end) or (
                            int(angles[2]) >= begin and int(angles[2]) <= end):
                    angles[0] = index_to_atomtype[angles[0]]
                    angles[1] = index_to_atomtype[angles[1]]
                    angles[2] = index_to_atomtype[angles[2]]
                    angles[-1] = str('{:.5f}'.format(float(angles[-1])))
                    T_angle.append(angles)

    # Load proper
    while True:
        line = f_in.readline()
        if line == '[ dihedrals ] ; propers\n':
            line = f_in.readline()
            break
    while True:
        line = f_in.readline()
        if line == '\n' or line == ' \n':
            break
        else:
            temp = line.split(';')
            if temp[0] != '':
                m = re.search('mSeminario', temp[1])
                if m:
                    propers = temp[0].split()
                    if (int(propers[0]) >= begin and int(propers[0]) <= end
                        ) or (int(propers[1]) >= begin and int(propers[1]) <=
                              end) or (int(propers[2]) >= begin
                                       and int(propers[2]) <= end) or (
                                           int(propers[3]) >= begin
                                           and int(propers[3]) <= end):
                        propers[0] = index_to_atomtype[propers[0]]
                        propers[1] = index_to_atomtype[propers[1]]
                        propers[2] = index_to_atomtype[propers[2]]
                        propers[3] = index_to_atomtype[propers[3]]
                        T_proper.append(propers)

    f_in.close()
    return T_type, T_bond, T_angle, T_proper


def _loadchg(f_in, begin, end):

    T_charge = []
    s_charge = []
    while True:
        line = f_in.readline()
        if line == '':
            break
        charge = line.split()
        T_charge.append(charge[4])
    for i in range(begin, end):
        s_charge.append(T_charge[i])

    f_in.close()
    return s_charge


def _prmfile(type, charge, bond, angle, proper, f_out, x_state):

    print('\n;' + x_state, file=f_out)

    print('[ atoms ]', file=f_out)
    for i in range(len(type)):
        print(type[i], end=' ', file=f_out)
    print(file=f_out)

    print('\n[ charges ]', file=f_out)
    for i in range(len(charge)):
        print('{:.8f}'.format(float(charge[i])), end=' ', file=f_out)
    print(file=f_out)

    print('\n[ bonds ]', file=f_out)
    if bond == []:
        pass
    else:
        for i in range(len(bond)):
            for j in range(len(bond[i])):
                print(bond[i][j], end=' ', file=f_out)
        print(file=f_out)

    print('\n[ angles ]', file=f_out)
    if angle == []:
        pass
    else:
        for i in range(len(angle)):
            for j in range(len(angle[i])):
                print(angle[i][j], end=' ', file=f_out)
            print(file=f_out)

    print('\n[ dihedrals ]', file=f_out)
    if proper == []:
        pass
    else:
        for i in range(len(proper)):
            for j in range(len(proper[i])):
                print(proper[i][j], end=' ', file=f_out)
            print(file=f_out)


def rectpoint(Mol):
    reactlist = []
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


def connectmol(Mol, Relist):

    Numatom = []
    NewMol = []
    for mol in Mol:
        Numatom.append(mol.GetNumAtoms())

    # I-C(0-4)
    Ir = Relist[0][0]
    Cr = Numatom[0] + Relist[4][0]
    combo = Chem.CombineMols(Mol[0], Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Ir, Cr, order=Chem.rdchem.BondType.SINGLE)
    # edcombo.RemoveBond(Numatom[0]+Relist[4][0], Numatom[0]+Relist[4][1])
    # edcombo.AddBond(Numatom[0]+Relist[4][0], Numatom[0]+Relist[4][1], order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Ir).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol1 = edcombo

    # R-C(1-4)
    Rr = Relist[1][0]
    Cr = Numatom[1] + Relist[4][0]
    combo = Chem.CombineMols(Mol[1], Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Rr, Cr, order=Chem.rdchem.BondType.SINGLE)
    # edcombo.RemoveBond(Numatom[0]+Relist[4][0], Numatom[0]+Relist[4][1])
    # edcombo.AddBond(Numatom[0]+Relist[4][0], Numatom[0]+Relist[4][1], order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Rr).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol2 = edcombo

    # I-MA-C(0-2-4)
    Ir = Relist[0][0]
    MAr1 = Numatom[0] + Relist[2][0]
    MAr2 = Numatom[0] + Relist[2][1]
    Cr = Numatom[0] + Numatom[2] + Relist[4][0]
    combo = Chem.CombineMols(Mol[0], Mol[2])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Ir, MAr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr2, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr1, MAr2)
    edcombo.AddBond(MAr1, MAr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Ir).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol3 = edcombo

    # R-MA-C(1-2-4)
    Rr = Relist[1][0]
    MAr1 = Numatom[1] + Relist[2][0]
    MAr2 = Numatom[1] + Relist[2][1]
    Cr = Numatom[1] + Numatom[2] + Relist[4][0]
    combo = Chem.CombineMols(Mol[1], Mol[2])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Rr, MAr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr2, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr1, MAr2)
    edcombo.AddBond(MAr1, MAr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Rr).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol4 = edcombo

    # I-MB-C(0-3-4)
    Ir = Relist[0][0]
    MBr1 = Numatom[0] + Relist[3][0]
    MBr2 = Numatom[0] + Relist[3][1]
    Cr = Numatom[0] + Numatom[3] + Relist[4][0]
    combo = Chem.CombineMols(Mol[0], Mol[3])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Ir, MBr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr2, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr1, MBr2)
    edcombo.AddBond(MBr1, MBr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Ir).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol5 = edcombo

    # R-MB-C(1-3-4)
    Rr = Relist[1][0]
    MBr1 = Numatom[1] + Relist[3][0]
    MBr2 = Numatom[1] + Relist[3][1]
    Cr = Numatom[1] + Numatom[3] + Relist[4][0]
    combo = Chem.CombineMols(Mol[1], Mol[3])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Rr, MBr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr2, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr1, MBr2)
    edcombo.AddBond(MBr1, MBr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Rr).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol6 = edcombo

    # I-MA-MA-C(0-2-2-4)
    Ir = Relist[0][0]
    MAr1 = Numatom[0] + Relist[2][0]
    MAr2 = Numatom[0] + Relist[2][1]
    MAr3 = Numatom[0] + Numatom[2] + Relist[2][0]
    MAr4 = Numatom[0] + Numatom[2] + Relist[2][1]
    Cr = Numatom[0] + Numatom[2] + Numatom[2] + Relist[4][0]
    combo = Chem.CombineMols(Mol[0], Mol[2])
    combo = Chem.CombineMols(combo, Mol[2])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Ir, MAr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr2, MAr3, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr4, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr1, MAr2)
    edcombo.AddBond(MAr1, MAr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr3, MAr4)
    edcombo.AddBond(MAr3, MAr4, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Ir).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol7 = edcombo

    # I-MB-MB-C(0-3-3-4)
    Ir = Relist[0][0]
    MBr1 = Numatom[0] + Relist[3][0]
    MBr2 = Numatom[0] + Relist[3][1]
    MBr3 = Numatom[0] + Numatom[3] + Relist[3][0]
    MBr4 = Numatom[0] + Numatom[3] + Relist[3][1]
    Cr = Numatom[0] + Numatom[3] + Numatom[3] + Relist[4][0]
    combo = Chem.CombineMols(Mol[0], Mol[3])
    combo = Chem.CombineMols(combo, Mol[3])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Ir, MBr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr2, MBr3, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr4, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr1, MBr2)
    edcombo.AddBond(MBr1, MBr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr3, MBr4)
    edcombo.AddBond(MBr3, MBr4, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Ir).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol8 = edcombo

    # R-MA-MA-MA-C(1-2-2-2-4)
    Rr = Relist[1][0]
    MAr1 = Numatom[1] + Relist[2][0]
    MAr2 = Numatom[1] + Relist[2][1]
    MAr3 = Numatom[1] + Numatom[2] + Relist[2][0]
    MAr4 = Numatom[1] + Numatom[2] + Relist[2][1]
    MAr5 = Numatom[1] + Numatom[2] + Numatom[2] + Relist[2][0]
    MAr6 = Numatom[1] + Numatom[2] + Numatom[2] + Relist[2][1]
    Cr = Numatom[1] + Numatom[2] + Numatom[2] + Numatom[2] + Relist[4][0]
    combo = Chem.CombineMols(Mol[1], Mol[2])
    combo = Chem.CombineMols(combo, Mol[2])
    combo = Chem.CombineMols(combo, Mol[2])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Rr, MAr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr2, MAr3, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr4, MAr5, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr6, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr1, MAr2)
    edcombo.AddBond(MAr1, MAr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr3, MAr4)
    edcombo.AddBond(MAr3, MAr4, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr5, MAr6)
    edcombo.AddBond(MAr5, MAr6, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Rr).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol9 = edcombo

    # R-MB-MB-MB-C(1-3-3-3-4)
    Rr = Relist[1][0]
    MBr1 = Numatom[1] + Relist[3][0]
    MBr2 = Numatom[1] + Relist[3][1]
    MBr3 = Numatom[1] + Numatom[3] + Relist[3][0]
    MBr4 = Numatom[1] + Numatom[3] + Relist[3][1]
    MBr5 = Numatom[1] + Numatom[3] + Numatom[3] + Relist[3][0]
    MBr6 = Numatom[1] + Numatom[3] + Numatom[3] + Relist[3][1]
    Cr = Numatom[1] + Numatom[3] + Numatom[3] + Numatom[3] + Relist[4][0]
    combo = Chem.CombineMols(Mol[1], Mol[3])
    combo = Chem.CombineMols(combo, Mol[3])
    combo = Chem.CombineMols(combo, Mol[3])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Rr, MBr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr2, MBr3, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr4, MBr5, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr6, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr1, MBr2)
    edcombo.AddBond(MBr1, MBr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr3, MBr4)
    edcombo.AddBond(MBr3, MBr4, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr5, MBr6)
    edcombo.AddBond(MBr5, MBr6, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Rr).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol10 = edcombo

    # R-MA-MB-MB-C(1-2-3-3-4)
    Rr = Relist[1][0]
    MAr1 = Numatom[1] + Relist[2][0]
    MAr2 = Numatom[1] + Relist[2][1]
    MBr3 = Numatom[1] + Numatom[2] + Relist[3][0]
    MBr4 = Numatom[1] + Numatom[2] + Relist[3][1]
    MBr5 = Numatom[1] + Numatom[2] + Numatom[3] + Relist[3][0]
    MBr6 = Numatom[1] + Numatom[2] + Numatom[3] + Relist[3][1]
    Cr = Numatom[1] + Numatom[2] + Numatom[3] + Numatom[3] + Relist[4][0]
    combo = Chem.CombineMols(Mol[1], Mol[2])
    combo = Chem.CombineMols(combo, Mol[3])
    combo = Chem.CombineMols(combo, Mol[3])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Rr, MAr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr2, MBr3, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr4, MBr5, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr6, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr1, MAr2)
    edcombo.AddBond(MAr1, MAr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr3, MBr4)
    edcombo.AddBond(MBr3, MBr4, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr5, MBr6)
    edcombo.AddBond(MBr5, MBr6, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Rr).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol11 = edcombo

    # R-MA-MA-MB-C(1-2-2-3-4)
    Rr = Relist[1][0]
    MAr1 = Numatom[1] + Relist[2][0]
    MAr2 = Numatom[1] + Relist[2][1]
    MAr3 = Numatom[1] + Numatom[2] + Relist[2][0]
    MAr4 = Numatom[1] + Numatom[2] + Relist[2][1]
    MBr5 = Numatom[1] + Numatom[2] + Numatom[2] + Relist[3][0]
    MBr6 = Numatom[1] + Numatom[2] + Numatom[2] + Relist[3][1]
    Cr = Numatom[1] + Numatom[2] + Numatom[2] + Numatom[3] + Relist[4][0]
    combo = Chem.CombineMols(Mol[1], Mol[2])
    combo = Chem.CombineMols(combo, Mol[2])
    combo = Chem.CombineMols(combo, Mol[3])
    combo = Chem.CombineMols(combo, Mol[4])
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Rr, MAr1, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr2, MAr3, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MAr4, MBr5, order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(MBr6, Cr, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr1, MAr2)
    edcombo.AddBond(MAr1, MAr2, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MAr3, MAr4)
    edcombo.AddBond(MAr3, MAr4, order=Chem.rdchem.BondType.SINGLE)
    edcombo.RemoveBond(MBr5, MBr6)
    edcombo.AddBond(MBr5, MBr6, order=Chem.rdchem.BondType.SINGLE)
    edcombo = edcombo.GetMol()
    edcombo.GetAtomWithIdx(Rr).SetFormalCharge(0)
    edcombo.GetAtomWithIdx(Cr).SetFormalCharge(0)
    mol12 = edcombo

    NewMol.append(mol1)
    NewMol.append(mol2)
    NewMol.append(mol3)
    NewMol.append(mol4)
    NewMol.append(mol5)
    NewMol.append(mol6)
    NewMol.append(mol7)
    NewMol.append(mol8)
    NewMol.append(mol9)
    NewMol.append(mol10)
    NewMol.append(mol11)
    NewMol.append(mol12)

    return NewMol


def GAFFforcefile():

    # GAFF力场
    T_type = []  # The atom type of GAFF
    T_bond = []  # The bond created by GAFF
    T_angle = []  # The angle created by GAFF
    T_proper = []  # The proper created by GAFF
    T_improper = []  # The improper created by GAFF
    folder = os.path.exists('./topper')
    if not folder:
        os.makedirs('./topper')

    filelist = [
        'IC', 'RC', 'IMAC', 'RMAC', 'IMBC', 'RMBC', 'IMA2C', 'IMB2C', 'RMA3C',
        'RMB3C', 'RMAMB2C', 'RMA2MBC', 'MA', 'MB'
    ]
    for filename in filelist:

        f_in = open('./gaussian/' + filename + '.itp', 'r')

        # load T_Atomtype
        while True:
            line = f_in.readline()
            if line == '[ atomtypes ]\n':
                line = f_in.readline()
                break
        while True:
            line = f_in.readline()
            if line == '\n' or line == ' \n':
                break
            else:
                types = line.split()
                T_type.append(types)

        # Load Atom
        index_to_atomtype = {}
        while True:
            line = f_in.readline()
            if line == '[ atoms ]\n':
                line = f_in.readline()
                break
        while True:
            line = f_in.readline()
            if line == '\n' or line == ' \n':
                break
            else:
                atoms = line.split()
                index_to_atomtype[atoms[0]] = atoms[1]

        # Load bond
        while True:
            if line == '':
                break
            line = f_in.readline()
            if line == '[ bonds ]\n':
                line = f_in.readline()
                break
        while True:
            if line == '':
                break
            line = f_in.readline()
            if line == '\n' or line == ' \n':
                break
            else:
                temp = line.split(';')
                m = re.search('mSeminario', temp[1])
                if not m:
                    bonds = temp[0].split()
                    bonds[0] = index_to_atomtype[bonds[0]]
                    bonds[1] = index_to_atomtype[bonds[1]]
                    T_bond.append(bonds)

        # Load angle
        while True:
            if line == '':
                break
            line = f_in.readline()
            if line == '[ angles ]\n':
                line = f_in.readline()
                break
        while True:
            if line == '':
                break
            line = f_in.readline()
            if line == '\n' or line == ' \n':
                break
            else:
                temp = line.split(';')
                m = re.search('mSeminario', temp[1])
                if not m:
                    angles = temp[0].split()
                    angles[0] = index_to_atomtype[angles[0]]
                    angles[1] = index_to_atomtype[angles[1]]
                    angles[2] = index_to_atomtype[angles[2]]
                    T_angle.append(angles)

        # Load proper
        while True:
            if line == '':
                break
            line = f_in.readline()
            if line == '[ dihedrals ] ; propers\n':
                line = f_in.readline()
                break
        while True:
            if line == '':
                break
            line = f_in.readline()
            if line == '\n' or line == ' \n':
                break
            else:
                temp = line.split(';')
                if temp[0] != '':
                    m = re.search('mSeminario', temp[1])
                    if not m:
                        propers = temp[0].split()
                        propers[0] = index_to_atomtype[propers[0]]
                        propers[1] = index_to_atomtype[propers[1]]
                        propers[2] = index_to_atomtype[propers[2]]
                        propers[3] = index_to_atomtype[propers[3]]
                        T_proper.append(propers)

        # Load improper
        while True:
            if line == '':
                break
            line = f_in.readline()
            if line == '[ dihedrals ] ; impropers\n':
                line = f_in.readline()
                break
        while True:
            if line == '':
                break
            line = f_in.readline()
            if line == '':
                break
            else:
                temp = line.split(';')
                if temp[0] != '':
                    m = re.search('mSeminario', temp[1])
                    if not m:
                        impropers = temp[0].split()
                        impropers[0] = index_to_atomtype[impropers[0]]
                        impropers[1] = index_to_atomtype[impropers[1]]
                        impropers[2] = index_to_atomtype[impropers[2]]
                        impropers[3] = index_to_atomtype[impropers[3]]
                        T_improper.append(impropers)
        f_in.close()

    # Remove the same
    T_typesort = list(set([tuple(t) for t in T_type]))
    T_bondsort = list(set([tuple(t) for t in T_bond]))
    T_anglesort = list(set([tuple(t) for t in T_angle]))
    T_propersort = list(set([tuple(t) for t in T_proper]))
    T_propersort.sort(key=lambda ele: (ele[0], ele[1], ele[2], ele[3], ele[7]))
    T_impropersort = list(set([tuple(t) for t in T_improper]))

    f_pair = open('./topper/ffnonbonded.itp', 'w')
    f_bond = open('./topper/ffbonded.itp', 'w')
    f_force = open('./topper/forcefield.itp', 'w')
    f_sys = open('./control/system.json', )  # Read version

    print('#define _FF_GAFF', file=f_force)
    print('\n[ defaults ]', file=f_force)
    print(';nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ', file=f_force)
    print('   1         2         yes       0.5     0.8333', file=f_force)
    print(file=f_force)
    print('#include "ffnonbonded.itp"', file=f_force)
    print('#include "ffbonded.itp"', file=f_force)

    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=f_pair)
    print('\n[ atomtypes ]', file=f_pair)
    print(
        '; name   at.num      mass       charge   ptype     sigma (nm)    epsilon (kJ/mol)',
        file=f_pair)
    for i in range(len(T_typesort)):
        type = T_typesort[i][0]
        atnum = T_typesort[i][1]
        mass = T_typesort[i][2]
        charge = T_typesort[i][3]
        ptype = T_typesort[i][4]
        sigma = str('{:.8f}'.format(float(T_typesort[i][5])))
        epsilon = str('{:.8f}'.format(float(T_typesort[i][6])))
        print(type.rjust(5),
              atnum.rjust(7),
              mass.rjust(13),
              charge.rjust(11),
              ptype.rjust(4),
              sigma.rjust(16),
              epsilon.rjust(16),
              file=f_pair)

    print(';created with PISAMD version ' + version, file=f_bond)
    print(file=f_bond)
    print('[ bondtypes ]', file=f_bond)
    print(';  atom_i  atom_j   functype    r0 (nm)      k (kJ/mol/nm^2)',
          file=f_bond)
    for i in range(len(T_bondsort)):
        atom_i = T_bondsort[i][0]
        atom_j = T_bondsort[i][1]
        functype = T_bondsort[i][2]
        r0 = T_bondsort[i][3]
        k = str('{:.2f}'.format(float(T_bondsort[i][4])))
        print(atom_i.rjust(7),
              atom_j.rjust(7),
              functype.rjust(8),
              r0.rjust(14),
              k.rjust(17),
              file=f_bond)

    print('\n[ angletypes ]', file=f_bond)
    print(
        ';  atom_i  atom_j  atom_k   functype   a0 (Deg.)      k (kJ/mol/rad^2)',
        file=f_bond)
    for i in range(len(T_anglesort)):
        atom_i = T_anglesort[i][0]
        atom_j = T_anglesort[i][1]
        atom_k = T_anglesort[i][2]
        functype = T_anglesort[i][3]
        a0 = T_anglesort[i][4]
        k = str('{:.5f}'.format(float(T_anglesort[i][5])))
        print(atom_i.rjust(7),
              atom_j.rjust(7),
              atom_k.rjust(7),
              functype.rjust(8),
              a0.rjust(14),
              k.rjust(18),
              file=f_bond)

    print('\n[ dihedraltypes ]', file=f_bond)
    print(
        ';  atom_i  atom_j  atom_k  atom_l   functype   phase (Deg.)    kd (kJ/mol)   pn',
        file=f_bond)
    for i in range(len(T_propersort)):
        atom_i = T_propersort[i][0]
        atom_j = T_propersort[i][1]
        atom_k = T_propersort[i][2]
        atom_l = T_propersort[i][3]
        functype = T_propersort[i][4]
        phase = T_propersort[i][5]
        kd = T_propersort[i][6]
        pn = T_propersort[i][7]
        print(atom_i.rjust(7),
              atom_j.rjust(7),
              atom_k.rjust(7),
              atom_l.rjust(7),
              functype.rjust(8),
              phase.rjust(15),
              kd.rjust(16),
              pn.rjust(5),
              file=f_bond)

    print('\n[ dihedraltypes ]', file=f_bond)
    print(
        ';  atom_i  atom_j  atom_k  atom_l   functype   phase (Deg.)    kd (kJ/mol)   pn',
        file=f_bond)
    for i in range(len(T_impropersort)):
        atom_i = T_impropersort[i][0]
        atom_j = T_impropersort[i][1]
        atom_k = T_impropersort[i][2]
        atom_l = T_impropersort[i][3]
        functype = T_impropersort[i][4]
        phase = T_impropersort[i][5]
        kd = T_impropersort[i][6]
        pn = T_impropersort[i][7]
        print(atom_i.rjust(7),
              atom_j.rjust(7),
              atom_k.rjust(7),
              atom_l.rjust(7),
              functype.rjust(8),
              phase.rjust(15),
              kd.rjust(16),
              pn.rjust(5),
              file=f_bond)
    f_bond.close()
    f_pair.close()


def statefile():

    Numatom = _Numatom()

    # I1 & I2
    begin = 0
    end = Numatom[0]
    f_out = open('./topper/I.prm', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=f_out)

    f_itp = open('./gaussian/IMAC.itp', 'r')
    f_chg = open('./gaussian/I_state1.chg', 'r')  # I1
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'I1')

    begin = 0
    end = Numatom[0]
    f_itp = open('./gaussian/IC.itp', 'r')
    f_chg = open('./gaussian/I_state2.chg', 'r')  # I2
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'I2')

    f_out.close()

    # R1 & R2
    begin = 0
    end = Numatom[1]
    f_out = open('./topper/R.prm', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=f_out)

    f_itp = open('./gaussian/RMAC.itp', 'r')
    f_chg = open('./gaussian/R_state1.chg', 'r')  # I1
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'R1')

    begin = 0
    end = Numatom[1]
    f_itp = open('./gaussian/RC.itp', 'r')
    f_chg = open('./gaussian/R_state2.chg', 'r')  # I2
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'R2')

    f_out.close()

    # C1 & C2 & C3
    begin = Numatom[1]
    end = Numatom[1] + Numatom[4]
    f_out = open('./topper/C.prm', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=f_out)

    f_itp = open('./gaussian/RC.itp', 'r')
    f_chg = open('./gaussian/C_state1.chg', 'r')  # C1
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'C1')

    begin = Numatom[1] + Numatom[2]
    end = Numatom[1] + Numatom[2] + Numatom[4]
    f_itp = open('./gaussian/RMAC.itp', 'r')
    f_chg = open('./gaussian/C_state2.chg', 'r')  # C2
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'C2')

    begin = Numatom[1] + Numatom[3]
    end = Numatom[1] + Numatom[3] + Numatom[4]
    f_itp = open('./gaussian/RMBC.itp', 'r')
    f_chg = open('./gaussian/C_state3.chg', 'r')  # C3
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'C3')

    f_out.close()

    # MA1 & MA2 & MA3 & MA4 & MA5
    begin = 0
    end = Numatom[2]
    f_out = open('./topper/MA.prm', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=f_out)

    f_itp = open('./gaussian/MA.itp', 'r')
    f_chg = open('./gaussian/MA_state1.chg', 'r')  # MA1
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MA1')

    begin = Numatom[1]
    end = Numatom[1] + Numatom[2]
    f_itp = open('./gaussian/RMA3C.itp', 'r')
    f_chg = open('./gaussian/MA_state2.chg', 'r')  # MA2
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MA2')

    begin = Numatom[1] + Numatom[2]
    end = Numatom[1] + Numatom[2] + Numatom[2]
    f_itp = open('./gaussian/RMA3C.itp', 'r')
    f_chg = open('./gaussian/MA_state3.chg', 'r')  # MA3
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MA3')

    begin = Numatom[1] + Numatom[2] + Numatom[2]
    end = Numatom[1] + Numatom[2] + Numatom[2] + Numatom[2]
    f_itp = open('./gaussian/RMA3C.itp', 'r')
    f_chg = open('./gaussian/MA_state4.chg', 'r')  # MA4
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MA4')

    begin = Numatom[1] + Numatom[2]
    end = Numatom[1] + Numatom[2] + Numatom[2]
    f_itp = open('./gaussian/RMA2MBC.itp', 'r')
    f_chg = open('./gaussian/MA_state5.chg', 'r')  # MA5
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MA5')

    f_out.close()

    # MB1 & MB2 & MB3 & MB4
    begin = 0
    end = Numatom[3]
    f_out = open('./topper/MB.prm', 'w')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=f_out)

    f_itp = open('./gaussian/MB.itp', 'r')
    f_chg = open('./gaussian/MB_state1.chg', 'r')  # MB1
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MB1')

    begin = Numatom[1] + Numatom[2]
    end = Numatom[1] + Numatom[2] + Numatom[3]
    f_itp = open('./gaussian/RMAMB2C.itp', 'r')
    f_chg = open('./gaussian/MB_state2.chg', 'r')  # MB2
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MB2')

    begin = Numatom[1] + Numatom[3]
    end = Numatom[1] + Numatom[3] + Numatom[3]
    f_itp = open('./gaussian/RMB3C.itp', 'r')
    f_chg = open('./gaussian/MB_state3.chg', 'r')  # MB3
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MB3')

    begin = Numatom[1] + Numatom[2] + Numatom[2]
    end = Numatom[1] + Numatom[2] + Numatom[2] + Numatom[3]
    f_itp = open('./gaussian/RMA2MBC.itp', 'r')
    f_chg = open('./gaussian/MB_state4.chg', 'r')  # MB4
    type, bond, angle, proper = _loadstr(f_itp, begin + 1, end)
    charge = _loadchg(f_chg, begin, end)
    _prmfile(type, charge, bond, angle, proper, f_out, 'MB4')

    f_out.close()


def reactpot_typelist(filename):

    atomtypedict = {}
    reactpotlist = []
    reactlist = ['i_p', 'r_p', 'c_n', 'm_n']

    f_itp = open('./topper/ffnonbonded.itp', 'r')
    while True:
        line = f_itp.readline()
        if line == '[ atomtypes ]\n':
            line = f_itp.readline()
            break
    while True:
        line = f_itp.readline()
        if line == '' or line == '\n':
            break
        lines = line.split()
        atomtypedict[lines[0]] = lines[5]
    f_itp.close()

    f_pisa = open('./gromacs/' + filename + '/PISA.itp', 'r')
    while True:
        line = f_pisa.readline()
        if line == '[ atoms ]\n':
            break
    while True:
        line = f_pisa.readline()
        if line == '' or line == '\n' or line == ' \n':
            break
        lines = line.split()
        if lines[4] in reactlist:
            temp = [0] * 7
            temp[0], temp[1], temp[2], temp[3] = lines[0], lines[2], lines[
                4], atomtypedict[lines[1]]
            reactpotlist.append(temp)
    f_pisa.close()

    return reactpotlist


def reactpot_poslist(filename, frame, reactpotlist):

    reactlist = ['i_p', 'r_p', 'c_n', 'm_n']
    f_gro = open('./gromacs/' + filename + '/' + str(frame) + '.gro', 'r')
    line = f_gro.readline()
    line = f_gro.readline()
    cout = 0
    while True:
        line = f_gro.readline()
        atomname = line[10:15].strip()
        x = line[20:28].strip()
        y = line[28:36].strip()
        z = line[36:44].strip()
        if atomname in reactlist:
            reactpotlist[cout][4], reactpotlist[cout][5], reactpotlist[cout][
                6] = x, y, z
            cout += 1
        if cout == len(reactpotlist):
            break

    f_gro.close()
    return reactpotlist


def CTA_Mollist():
    # Get Natom of models
    Mol = []
    model = ['I.mol', 'R.mol', 'C.mol', 'MA.mol']
    for i in range(len(model)):
        mol = Chem.MolFromMolFile('./model/' + model[i])
        mol = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        Mol.append(mol)
    Numatom = []
    for mol in Mol:
        temp = mol.GetNumAtoms()
        Numatom.append(temp)
    # Get Nmol
    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    boxsize = system['Boxsize']
    polymer = system['PISA']
    x, y, z = boxsize[0], boxsize[1], boxsize[2]
    Pini = 0
    if 'Pini_in_macroCTA' in system:
        Pini = system['Pini_in_macroCTA']
    mCTA, mINI, mMA = polymer[0]['CTA'], polymer[0]['INI'], polymer[0]['MA']
    nINI = int(round(0.6023 * x * y * z * mINI + 0.5))
    if nINI == 0:
        nINI = 1
    nCTA = int(round(nINI * (mCTA / mINI)))
    if nCTA == 0:
        nCTA = 1
    nMA = int(round(nCTA * (mMA / mCTA)))
    if nMA == 0:
        nMA = 1
    nI = int(nCTA * Pini)
    nR = nCTA - nI
    nC = nCTA * 100
    if nR == 0:
        nR = 1
    if nI == 0:
        nI = 1
    f_sys.close()

    Ilist = []
    Rlist = []
    Clist = []
    MAlist = []
    # R & C
    for i in range(nI):
        Itemp = [0] * 3
        Itemp[0] = i + 1
        Itemp[1] = i * Numatom[0] + 1
        Itemp[2] = (i + 1) * Numatom[0]
        Ilist.append(Itemp)
    for i in range(nR):
        Rtemp = [0] * 3
        Rtemp[0] = i + 1 + nI
        Rtemp[1] = nI * Numatom[0] + i * Numatom[1] + 1
        Rtemp[2] = nI * Numatom[0] + (i + 1) * Numatom[1]
        Rlist.append(Rtemp)
    for i in range(nC):
        Ctemp = [0] * 3
        Ctemp[0] = i + 1 + nI + nR
        Ctemp[1] = nI * Numatom[0] + nR * Numatom[1] + i * Numatom[2] + 1
        Ctemp[2] = nI * Numatom[0] + nR * Numatom[1] + (i + 1) * Numatom[2]
        Clist.append(Ctemp)
    for i in range(nMA):
        MAtemp = [0] * 3
        MAtemp[0] = i + 1 + nI + nR + nC
        MAtemp[1] = nI * Numatom[0] + nR * Numatom[1] + nC * Numatom[
            2] + i * Numatom[3] + 1
        MAtemp[2] = nI * Numatom[0] + nR * Numatom[1] + nC * Numatom[2] + (
            i + 1) * Numatom[3]
        MAlist.append(MAtemp)
    return Ilist, Rlist, Clist, MAlist


def PISA_Mollist():
    # Get Natom of models
    Mol = []
    model = ['I.mol', 'R.mol', 'C.mol', 'MA.mol', 'MB.mol']
    for i in range(len(model)):
        mol = Chem.MolFromMolFile('./model/' + model[i])
        mol = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        Mol.append(mol)
    Numatom = []
    for mol in Mol:
        temp = mol.GetNumAtoms()
        Numatom.append(temp)
    # Get Nmol
    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    boxsize = system['Boxsize']
    polymer = system['PISA']
    x, y, z = boxsize[0], boxsize[1], boxsize[2]
    Pini = 0
    if 'Pini_in_macroCTA' in system:
        Pini = system['Pini_in_macroCTA']
    mCTA, mINI, mMB = polymer[0]['CTA'], polymer[0]['INI'], polymer[0]['MB']
    nINI = int(round(0.6023 * x * y * z * mINI + 0.5))
    if nINI == 0:
        nINI = 1
    nCTA = int(round(nINI * (mCTA / mINI)))
    if nCTA == 0:
        nCTA = 1
    nMB = int(round(nCTA * (mMB / mCTA)))
    if nMB == 0:
        nMB = 1
    nI = int(nCTA * Pini)
    nR = nCTA - nI
    if nR == 0:
        nR = 1
    if nI == 0:
        nI = 1
    f_sys.close()
    nMA = 0
    fin_itp = open('./topper/PISA.itp', 'r')
    while True:
        line = fin_itp.readline()
        if line == '[ atoms ]\n':
            break
    while True:
        line = fin_itp.readline()
        lines = line.split()
        if lines[3] == 'MOB':
            nMA = int(lines[2]) - 1 - nI - nR - nCTA - nINI
            break
    fin_itp.close()

    Ilist = []
    Rlist = []
    Clist = []
    MAlist = []
    MBlist = []
    # R & C
    for i in range(nI):
        Itemp = [0] * 3
        Itemp[0] = i + 1
        Itemp[1] = i * Numatom[0] + 1
        Itemp[2] = (i + 1) * Numatom[0]
        Ilist.append(Itemp)
    for i in range(nR):
        Rtemp = [0] * 3
        Rtemp[0] = i + 1 + nI
        Rtemp[1] = nI * Numatom[0] + i * Numatom[1] + 1
        Rtemp[2] = nI * Numatom[0] + (i + 1) * Numatom[1]
        Rlist.append(Rtemp)
    for i in range(nCTA):
        Ctemp = [0] * 3
        Ctemp[0] = i + 1 + nI + nR
        Ctemp[1] = nI * Numatom[0] + nR * Numatom[1] + i * Numatom[2] + 1
        Ctemp[2] = nI * Numatom[0] + nR * Numatom[1] + (i + 1) * Numatom[2]
        Clist.append(Ctemp)
    for i in range(nMA):
        MAtemp = [0] * 3
        MAtemp[0] = i + 1 + nI + nR + nCTA
        MAtemp[1] = nI * Numatom[0] + nR * Numatom[1] + nCTA * Numatom[
            2] + i * Numatom[3] + 1
        MAtemp[2] = nI * Numatom[0] + nR * Numatom[1] + nCTA * Numatom[2] + (
            i + 1) * Numatom[3]
        MAlist.append(MAtemp)
    for i in range(nINI):
        Itemp = [0] * 3
        Itemp[0] = i + 1 + nI + nR + nCTA + nMA
        Itemp[1] = nI * Numatom[0] + nR * Numatom[1] + nCTA * Numatom[
            2] + nMA * Numatom[3] + i * Numatom[0] + 1
        Itemp[2] = nI * Numatom[0] + nR * Numatom[1] + nCTA * Numatom[
            2] + nMA * Numatom[3] + (i + 1) * Numatom[0]
        Ilist.append(Itemp)
    for i in range(nMB):
        MBtemp = [0] * 3
        MBtemp[0] = i + 1 + nI + nR + nCTA + nMA + nINI
        MBtemp[1] = nI * Numatom[0] + nR * Numatom[1] + nCTA * Numatom[
            2] + nMA * Numatom[3] + nINI * Numatom[0] + i * Numatom[4] + 1
        MBtemp[2] = nI * Numatom[0] + nR * Numatom[1] + nCTA * Numatom[
            2] + nMA * Numatom[3] + nINI * Numatom[0] + (i + 1) * Numatom[4]
        MBlist.append(MBtemp)
    return Ilist, Rlist, Clist, MAlist, MBlist
