import json
import os
import re
import subprocess
import numpy as np
from rdkit import Chem
from bin.jointmol import Normalsystem
from bin.force import Normalgaussian

Ionlist = {
    '[Be+2]': ['BE2P', 4, 9.012200, 0.007127189745, 0.4250944],
    '[Al+3]': ['AL3P', 13, 26.982000, 0.039377723342, 0.4460144],
    '[Sc+3]': ['SC3P', 21, 44.956000, 0.195445360786, 0.5899440],
    '[Ti+]': ['TI1P', 22, 47.867000, 0.344011631023, 1.7016328],
    '[Ti+3]': ['TI3P', 22, 47.867000, 0.147693189493, 0.5376440],
    '[V+2]': ['V2P', 23, 50.942000, 0.196656983042, 0.6732056],
    '[V+3]': ['V3P', 23, 50.942000, 0.109437998536, 0.5171424],
    '[Cr+2]': ['CR2P', 24, 51.996000, 0.188727984451, 0.6974728],
    '[Cr+3]': ['CR3P', 24, 51.996000, 0.150223341853, 0.5045904],
    '[Mn+2]': ['MN2P', 25, 54.938000, 0.214474957405, 0.7050040],
    '[Mn+3]': ['MN3P', 25, 54.938000, 0.068492293451, 0.5238368],
    '[Fe+2]': ['FE2P', 26, 55.854000, 0.192451941093, 0.6673480],
    '[Fe+3]': ['FE3P', 26, 55.854000, 0.100903188817, 0.5238368],
    '[Co+2]': ['CO2P', 27, 58.933000, 0.165279530189, 0.6439176],
    '[Co+3]': ['CO3P', 27, 58.933000, 0.045988191830, 0.4983144],
    '[Ni+2]': ['NI2P', 28, 58.693000, 0.139033653953, 0.5999856],
    '[Cu+2]': ['CU2P', 29, 63.546000, 0.129198132105, 0.6296920],
    '[Ga+3]': ['GA3P', 31, 69.723000, 0.040428983829, 0.5045904],
    '[Sr+2]': ['SR2P', 38, 87.620000, 0.330683786199, 0.9572992],
    '[Y+3]': ['Y3P', 39, 88.906000, 0.251839249644, 0.6928704],
    '[Rh+3]': ['RH3P', 45, 102.910000, 0.082194315736, 0.5376440],
    '[Pd+2]': ['PD2P', 46, 106.420000, 0.165190440318, 0.7288528],
    '[Ag+2]': ['AG2P', 47, 107.870000, 0.186910551066, 0.7752952],
    '[In+3]': ['IN3P', 49, 114.820000, 0.154357111905, 0.6163032],
    '[Sn+2]': ['SN2P', 50, 118.710000, 0.295422014935, 0.7861736],
    '[La+3]': ['LA3P', 57, 138.910000, 0.315859231529, 0.7932864],
    '[Ce+3]': ['CE3P', 58, 140.120000, 0.304954631219, 0.7702744],
    '[Pr+3]': ['PR3P', 59, 140.910000, 0.295618012653, 0.7627432],
    '[Nd+3]': ['ND3P', 60, 144.240000, 0.288633366703, 0.7497728],
    '[Pm+3]': ['PM3P', 61, 144.910000, 0.295279471140, 0.7426600],
    '[Sm+2]': ['SM2P', 62, 150.360000, 0.330790694046, 1.0242432],
    '[Sm+3]': ['SM3P', 62, 150.360000, 0.279154204342, 0.7351288],
    '[Eu+2]': ['EU2P', 63, 151.960000, 0.327922000173, 0.9907712],
    '[Eu+3]': ['EU3P', 63, 151.960000, 0.271367749546, 0.7280160],
    '[Gd+3]': ['GD3P', 64, 157.250000, 0.268160514160, 0.7209032],
    '[Tb+3]': ['TB3P', 65, 158.930000, 0.262815121851, 0.7083512],
    '[Dy+3]': ['DY3P', 66, 162.500000, 0.257362821696, 0.6999832],
    '[Ho+3]': ['HO3P', 67, 164.930000, 0.247794569464, 0.6928704],
    '[Er+3]': ['ER3P', 68, 167.280000, 0.243625163463, 0.6845024],
    '[Tm+3]': ['TM3P', 69, 168.930000, 0.240293202257, 0.6778080],
    '[Yb+2]': ['YB2P', 70, 173.050000, 0.286691207498, 0.8593936],
    '[Yb+3]': ['YB3P', 70, 173.050000, 0.230778403947, 0.6711136],
    '[Lu+3]': ['LU3P', 71, 174.970000, 0.240756469590, 0.6644192],
    '[Pt+2]': ['PT2P', 78, 195.080000, 0.142027073646, 0.6824104],
    '[Au+]': ['AU1P', 79, 196.970000, 0.112502690127, 1.4907592],
    '[Au+3]': ['AU3P', 79, 196.970000, 0.065106878322, 0.5564720],
    '[Hg+2]': ['HG2P', 80, 200.590000, 0.210572821020, 0.8593936],
    '[Tl+3]': ['TL3P', 81, 204.380000, 0.154891651136, 0.6815736],
    '[Pb+2]': ['PB2P', 82, 207.200000, 0.313239989298, 1.0016496],
    '[Bi+3]': ['BI3P', 83, 208.980000, 0.243785525232, 0.7740400],
    '[Ra+2]': ['RA2P', 88, 226.030000, 0.372377846208, 1.2631496],
    '[U+3]': ['U3P', 92, 238.030000, 0.303119379860, 0.7932864],
    '[Pu+3]': ['PU3P', 94, 244.060000, 0.297506717936, 0.7702744],
    '[Mg+2]': ['MG', 12, 24.305000, 0.211142996199, 0.0627600],
    '[Ca+2]': ['CAL', 20, 40.080000, 0.243571709540, 0.5020800],
    '[Zn+2]': ['ZN', 30, 65.370000, 0.194215920555, 1.0460000],
    '[Cd+2]': ['CAD', 48, 112.411000, 0.241789912103, 0.5020800],
    '[Ba+2]': ['BAR', 56, 137.327000, 0.336759715457, 0.6276000],
    '[Li+]': ['LIT', 3, 6.941000, 0.254000000000,
              0.013625196],  # for opc3 water
    '[Na+]': ['SOD', 11, 22.989770, 0.293800000000,
              0.126028774],  # for opc3 water
    '[K+]': ['POT', 19, 39.098300, 0.340600000000,
             0.586672238],  # for opc3 water
    '[Rb+]': ['RUB', 37, 85.467800, 0.360000000000,
              0.89173069],  # for opc3 water
    '[Cs+]': ['CES', 55, 132.905450, 0.393000000000,
              1.496323711],  # for opc3 water
    '[Tl+]': ['TL3P', 81, 204.380000, 0.336400000000,
              0.528388699],  # for opc3 water
    '[Cu+]': ['CU1P', 29, 63.546000, 0.240200000000,
              0.004698632],  # for opc3 water
    '[Ag+]': ['AG1P', 47, 107.870000, 0.267000000000,
              0.031871411],  # for opc3 water
    '[F-]': ['FA', 9, 18.998403, 0.363600000000,
             0.953762046],  # for opc3 water
    '[Cl-]': ['CLA', 17, 35.453200, 0.461200000000,
              2.687243956],  # for opc3 water
    '[Br-]': ['BRA', 35, 79.904000, 0.500200000000,
              3.18294026],  # for opc3 water
    '[I-]': ['IA', 53, 126.904473, 0.556000000000,
             3.647436058]  # for opc3 water
}


def _getstr(filename):

    f_in = open('./topper/' + filename + '.itp', 'r')
    type = []
    atom = []
    bond = []
    angle = []
    pair = []
    proper = []
    improper = []

    # type
    while True:
        line = f_in.readline()
        if line == '[ atomtypes ]\n':
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
            type.append(line)

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
            atom.append(line)

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
            bond.append(line)

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
            angle.append(line)

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
            proper.append(line)

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
            pair.append(line)

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
            improper.append(line)

    f_in.close()
    return type, atom, bond, angle, proper, pair, improper


def _getinformation(filename):

    f_in = open('./topper/' + filename + 'step2.chg', 'r')
    info = []
    while True:
        line = f_in.readline()
        if line == '':
            break
        lines = line.split()
        info.append(lines)

    f_in.close()
    return info


def _Formcharge(Mol):

    formcharge = 0
    for atom in Mol.GetAtoms():
        formcharge += atom.GetFormalCharge()

    return formcharge


def _getconnectmatrix(Natom, bondlist):

    matrix = np.zeros((Natom, Natom), dtype=int)
    for i in range(len(bondlist)):
        lines = bondlist[i].split()
        matrix[int(lines[0]) - 1][int(lines[1]) - 1] = 1
        matrix[int(lines[1]) - 1][int(lines[0]) - 1] = 1

    return matrix


def MDfile(filename):

    if len(filename) >= 3:
        Molname = filename[:3]
    else:
        index = 3 - len(filename)
        Molname = 'solvent'[:index] + filename
    Molname = Molname.upper()
    type, atom, bond, angle, proper, pair, improper = _getstr(filename)
    info = _getinformation(filename)
    matrix = _getconnectmatrix(len(atom), bond)

    f_out = open('./topper/' + filename + '.itp', 'w')
    fin_type = open('./topper/ffnonbonded.itp', 'r')
    types = fin_type.readlines()
    fin_type.close()
    fout_type = open('./topper/ffnonbonded.itp', 'a')
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    print(';created with PISAMD version ' + version, file=f_out)

    for j in range(len(type)):
        line = type[j].split()
        str_match = [x for x in types if re.search(line[0], x)]
        if str_match == []:
            print(line[0].rjust(5),
                  line[1].rjust(7),
                  line[2].rjust(13),
                  line[3].rjust(11),
                  line[4].rjust(4),
                  str('{:.8f}'.format(float(line[5]))).rjust(16),
                  str('{:.8f}'.format(float(line[6]))).rjust(16),
                  file=fout_type)
    fout_type.close()

    print('\n[ moleculetype ]', file=f_out)
    print(filename + '     3', file=f_out)

    print('\n[ atoms ]', file=f_out)
    for j in range(len(atom)):
        line = atom[j].split()
        print(line[0].rjust(6),
              line[1].rjust(6),
              line[2].rjust(6),
              Molname.rjust(6), (line[4] + line[0]).rjust(6),
              line[5].rjust(6),
              str('{:.8f}'.format(float(info[j][4]))).rjust(14),
              line[7].rjust(12),
              file=f_out)

    if bond != []:
        print('\n[ bonds ]', file=f_out)
        for j in range(len(bond)):
            lines = bond[j].split()
            print(lines[0].rjust(5),
                  lines[1].rjust(6),
                  lines[2].rjust(6),
                  lines[3].rjust(12),
                  str('{:.2f}'.format(float(lines[4]))).rjust(14),
                  file=f_out)

    if angle != []:
        print('\n[ angles ]', file=f_out)
        for j in range(len(angle)):
            lines = angle[j].split()
            print(lines[0].rjust(5),
                  lines[1].rjust(6),
                  lines[2].rjust(6),
                  lines[3].rjust(6),
                  lines[4].rjust(12),
                  str('{:.5f}'.format(float(lines[5]))).rjust(14),
                  file=f_out)

    if proper != []:
        print('\n[ dihedrals ]', file=f_out)
        for j in range(len(proper)):
            lines = proper[j].split()
            m = re.search('mSeminario', proper[j])
            if m:
                print(lines[0].rjust(5),
                      lines[1].rjust(6),
                      lines[2].rjust(6),
                      lines[3].rjust(6),
                      lines[4].rjust(6),
                      lines[5].rjust(12),
                      str('{:.5f}'.format(float(lines[6]))).rjust(14),
                      file=f_out)
            else:
                print(lines[0].rjust(5),
                      lines[1].rjust(6),
                      lines[2].rjust(6),
                      lines[3].rjust(6),
                      lines[4].rjust(6),
                      lines[5].rjust(12),
                      lines[6].rjust(14),
                      lines[7].rjust(6),
                      file=f_out)
    if pair != []:
        print('\n[ pairs ]', file=f_out)
        for j in range(len(pair)):
            lines = pair[j].split()
            print(lines[0].rjust(5),
                  lines[1].rjust(6),
                  lines[2].rjust(6),
                  file=f_out)

    if improper != []:
        print('\n[ dihedrals ]', file=f_out)
        for j in range(len(improper)):
            lines = improper[j].split()
            print(lines[0].rjust(5),
                  lines[1].rjust(6),
                  lines[2].rjust(6),
                  lines[3].rjust(6),
                  lines[4].rjust(6),
                  lines[5].rjust(12),
                  lines[6].rjust(14),
                  lines[7].rjust(6),
                  file=f_out)

    f_out.close()

    f_pdb = open('./topper/' + filename + '.pdb', 'w')
    print('REMARK  PISAMD PDB file', file=f_pdb)
    for i in range(len(info)):
        lines = info[i]
        print('ATOM',
              '  ',
              str(i + 1).rjust(5, ' '),
              ' ', (lines[0] + str(i + 1)).ljust(4, ' '),
              ' ',
              Molname.ljust(3, ' '),
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
              file=f_pdb)
    print('TER'.ljust(4, ' '),
          '  ',
          str(len(info) + 1).rjust(5, ' '),
          sep='',
          file=f_pdb)
    if bond != []:
        for i in range(len(matrix)):
            print('CONECT',
                  str(i + 1).rjust(5, ' '),
                  sep='',
                  end='',
                  file=f_pdb)
            for j in range(len(matrix[i])):
                if matrix[i][j] == 1:
                    print(str(j + 1).rjust(5, ' '), sep='', end='', file=f_pdb)
            print(file=f_pdb)
    print('END', file=f_pdb)
    f_pdb.close()

    subprocess.Popen('rm -f ' + '*.chg', shell=True, cwd='./topper')


def IonMDfile(filename):

    Ion = Chem.MolFromMolFile('./model/' + filename + '.mol')
    Ionsmi = Chem.MolToSmiles(Ion)
    Naion = Ion.GetNumAtoms()
    Ioncharge = _Formcharge(Ion)
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    f_sys.close()

    if Naion == 1:

        atom = Ion.GetAtoms()
        asym = atom[0].GetSymbol()
        resname = asym.upper() + str(abs(Ioncharge))

        if Ionsmi in Ionlist.keys():
            f_itp = open('./topper/' + filename + '.itp', 'w')
            f_pdb = open('./topper/' + filename + '.pdb', 'w')
            fout_type = open('./topper/ffnonbonded.itp', 'a')

            print(';created with PISAMD version ' + version, file=f_itp)

            print(Ionlist[Ionsmi][0].rjust(5),
                  str(Ionlist[Ionsmi][1]).rjust(7),
                  str('{:.6f}'.format(Ionlist[Ionsmi][2])).rjust(13),
                  '0.000000'.rjust(11),
                  'A'.rjust(4),
                  str('{:.8f}'.format(Ionlist[Ionsmi][3])).rjust(16),
                  str('{:.8f}'.format(Ionlist[Ionsmi][4])).rjust(16),
                  file=fout_type)

            print('\n[ moleculetype ]', file=f_itp)
            print(filename + '     3', file=f_itp)

            print('\n[ atoms ]', file=f_itp)
            print('1'.rjust(6),
                  Ionlist[Ionsmi][0].rjust(6),
                  '1'.rjust(6),
                  resname.rjust(6),
                  asym.rjust(6),
                  '1'.rjust(6),
                  str('{:.8f}'.format(float(Ioncharge))).rjust(14),
                  str('{:.6f}'.format(Ionlist[Ionsmi][2])).rjust(12),
                  file=f_itp)

            print('REMARK  PISAMD PDB file', file=f_pdb)
            print('ATOM',
                  '  ',
                  '1'.rjust(5, ' '),
                  ' ',
                  asym.ljust(4, ' '),
                  ' ',
                  resname.ljust(3, ' '),
                  '1'.rjust(6, ' '),
                  '    ',
                  '1.000'.rjust(8, ' '),
                  '1.000'.rjust(8, ' '),
                  '1.000'.rjust(8, ' '),
                  '1.00'.rjust(6, ' '),
                  '0.00'.rjust(6, ' '),
                  ' '.rjust(10, ' '),
                  asym.rjust(2, ' '),
                  sep='',
                  file=f_pdb)
            print('TER'.ljust(4, ' '),
                  '  ',
                  '2'.rjust(5, ' '),
                  sep='',
                  file=f_pdb)

            f_itp.close()
            f_pdb.close()
            fout_type.close()
            return True

        else:

            step1 = subprocess.Popen('obabel -imol ' + filename + '.mol' +
                                     ' -omol2 -O ' + filename + '.mol2',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True,
                                     cwd='./model')
            step1.wait()

            step2 = subprocess.Popen('mv *.mol2 ../topper',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True,
                                     cwd='./model')

            step2.wait()

            sobtop = open('./sobtop/' + filename + '.sh', 'w')
            print('chmod 777 sobtop',
                  'chmod 777 atomtype',
                  sep='\n',
                  file=sobtop)
            print('./sobtop  <<EOF', file=sobtop)
            print('../topper/', filename + '.mol2', sep='', file=sobtop)
            print(1, 1, 0, 1, sep='\n', file=sobtop)
            print('../topper/', filename + '.top', sep='', file=sobtop)
            print('../topper/', filename + '.itp', sep='', file=sobtop)
            print(0, file=sobtop)
            print('EOF', end='', file=sobtop)
            sobtop.close()

            step3 = subprocess.Popen('sh ' + filename + '.sh',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True,
                                     cwd='./sobtop')
            step3.wait()
            if not os.access('./topper/' + filename + '.itp', os.F_OK):
                print('Interaction caculating is wrong (' + filename + ')')
                return False

            else:

                f_in = open('./topper/' + filename + '.itp', 'r')
                type = []
                atom = []

                # type
                while True:
                    line = f_in.readline()
                    if line == '[ atomtypes ]\n':
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
                        type.append(line)

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
                        atom.append(line)
                f_in.close()

                f_out = open('./topper/' + filename + '.itp', 'w')
                f_pdb = open('./topper/' + filename + '.pdb', 'w')
                fout_type = open('./topper/ffnonbonded.itp', 'a')

                print(';created with PISAMD version ' + version, file=f_out)
                lines = type[0].split()
                print(lines[0].rjust(5),
                      lines[1].rjust(7),
                      lines[2].rjust(13),
                      lines[3].rjust(11),
                      lines[4].rjust(4),
                      str('{:.8f}'.format(float(lines[5]))).rjust(16),
                      str('{:.8f}'.format(float(lines[6]))).rjust(16),
                      file=fout_type)

                print('\n[ moleculetype ]', file=f_out)
                print(filename + '     3', file=f_out)

                print('\n[ atoms ]', file=f_out)

                print('1'.rjust(6),
                      lines[0].rjust(9),
                      '1'.rjust(6),
                      resname.rjust(6),
                      asym.rjust(6),
                      '1'.rjust(6),
                      str('{:.8f}'.format(float(Ioncharge))).rjust(14),
                      str('{:.6f}'.format(float(lines[2]))).rjust(12),
                      file=f_out)

                print('REMARK  PISAMD PDB file', file=f_pdb)
                print('ATOM',
                      '  ',
                      '1'.rjust(5, ' '),
                      ' ',
                      asym.ljust(4, ' '),
                      ' ',
                      resname.ljust(3, ' '),
                      '1'.rjust(6, ' '),
                      '    ',
                      '1.000'.rjust(8, ' '),
                      '1.000'.rjust(8, ' '),
                      '1.000'.rjust(8, ' '),
                      '1.00'.rjust(6, ' '),
                      '0.00'.rjust(6, ' '),
                      ' '.rjust(10, ' '),
                      asym.rjust(2, ' '),
                      sep='',
                      file=f_pdb)
                print('TER'.ljust(4, ' '),
                      '  ',
                      '2'.rjust(5, ' '),
                      sep='',
                      file=f_pdb)

                f_out.close()
                f_pdb.close()
                fout_type.close()

                return True
    else:
        Normalgaussian.gaussianfile(filename, mem=30)
        flow1 = Normalgaussian.gaussiancalculate(filename)
        if not flow1:
            return False

        flow2 = Normalgaussian.itpcalculate(filename)
        if not flow2:
            return False

        flow3 = Normalgaussian.chargecalculate(filename)
        if not flow3:
            return False

        if flow1 and flow2 and flow3:
            Normalsystem.MDfile(filename)
            return True


def waterMDfile(filename):

    water = Chem.MolFromMolFile('./model/' + filename + '.mol')
    watersmi = Chem.MolToSmiles(water)
    f_sys = open('./control/system.json', )  # Read version
    system = json.load(f_sys)
    version = '?'
    if 'Version' in system:
        version = system['Version']
    f_sys.close()

    if watersmi == 'O':
        f_itp = open('./topper/' + filename + '.itp', 'w')
        f_pdb = open('./topper/' + filename + '.pdb', 'w')
        fout_type = open('./topper/ffnonbonded.itp', 'a')

        print(';created with PISAMD version ' + version, file=f_itp)
        print(
            '  OW3       8     15.999400   -0.895170    A       0.31742700       0.68369000',
            file=fout_type)
        print(
            '  HW3       1      1.008000    0.447585    A       0.00000000       0.00000000',
            file=fout_type)

        print('\n[ moleculetype ]', file=f_itp)
        print(filename + '     2 ', file=f_itp)

        print('\n[ atoms ]', file=f_itp)
        print(' 1     OW3    1    WAT     OW     1    -0.895170', file=f_itp)
        print(' 2     HW3    1    WAT     HW1    1     0.447585', file=f_itp)
        print(' 3     HW3    1    WAT     HW2    1     0.447585', file=f_itp)

        print('\n[ bonds ]', file=f_itp)
        print(' 1   2    1     0.097888  502416.0', file=f_itp)
        print(' 1   3    1     0.097888  502416.0', file=f_itp)

        print('\n[ angles ]', file=f_itp)
        print(' 2   1   3    1     109.47  628.02', file=f_itp)

        print('REMARK  PISAMD PDB file', file=f_pdb)
        print('ATOM',
              '  ',
              '1'.rjust(5, ' '),
              ' ',
              'OW'.ljust(4, ' '),
              ' ',
              'WAT'.ljust(3, ' '),
              '1'.rjust(6, ' '),
              '    ',
              '-1.145'.rjust(8, ' '),
              '0.558'.rjust(8, ' '),
              '0.000'.rjust(8, ' '),
              '1.00'.rjust(6, ' '),
              '0.00'.rjust(6, ' '),
              ' '.rjust(10, ' '),
              'O'.rjust(2, ' '),
              sep='',
              file=f_pdb)
        print('ATOM',
              '  ',
              '2'.rjust(5, ' '),
              ' ',
              'HW1'.ljust(4, ' '),
              ' ',
              'WAT'.ljust(3, ' '),
              '1'.rjust(6, ' '),
              '    ',
              '-1.471'.rjust(8, ' '),
              '1.210'.rjust(8, ' '),
              '0.653'.rjust(8, ' '),
              '1.00'.rjust(6, ' '),
              '0.00'.rjust(6, ' '),
              ' '.rjust(10, ' '),
              'H'.rjust(2, ' '),
              sep='',
              file=f_pdb)
        print('ATOM',
              '  ',
              '3'.rjust(5, ' '),
              ' ',
              'HW2'.ljust(4, ' '),
              ' ',
              'WAT'.ljust(3, ' '),
              '1'.rjust(6, ' '),
              '    ',
              '-0.166'.rjust(8, ' '),
              '0.558'.rjust(8, ' '),
              '0.000'.rjust(8, ' '),
              '1.00'.rjust(6, ' '),
              '0.00'.rjust(6, ' '),
              ' '.rjust(10, ' '),
              'H'.rjust(2, ' '),
              sep='',
              file=f_pdb)
        print('TER'.ljust(4, ' '), '  ', '4'.rjust(5, ' '), sep='', file=f_pdb)
        print('CONECT    1    2    3', file=f_pdb)
        print('CONECT    2    1', file=f_pdb)
        print('CONECT    3    1', file=f_pdb)
        print('END', file=f_pdb)
        f_itp.close()
        f_pdb.close()
        fout_type.close()

        return True
    else:
        return False
