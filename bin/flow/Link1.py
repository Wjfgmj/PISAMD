import json
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from bin.jointmol import RAFTcreat
from bin.force import RAFTgaussian


def go():

    # Step One: built molecular models
    if os.access("./gaussian/__Gaussian__", os.F_OK):
        if os.access("./gaussian/__Itp__", os.F_OK):
            if os.access("./gaussian/__State__", os.F_OK):
                print("|", "Link1 Pass".center(40, "-"), "|", flush=True)
                return True
    print("|", "Link1".center(40, "-"), "|", flush=True)
    print("|", "Loading Models".center(40, "-"), "|", flush=True)
    Loading = 0
    cout = 0
    if os.access("./model/C.mol", os.F_OK):
        if os.access("./model/I.mol", os.F_OK):
            if os.access("./model/MA.mol", os.F_OK):
                if os.access("./model/MB.mol", os.F_OK):
                    if os.access("./model/R.mol", os.F_OK):
                        Loading += 1
    if os.access("./model/C.mol", os.R_OK):
        if os.access("./model/I.mol", os.R_OK):
            if os.access("./model/MA.mol", os.R_OK):
                if os.access("./model/MB.mol", os.R_OK):
                    if os.access("./model/R.mol", os.R_OK):
                        Loading += 1
    if os.access("./model/C.mol", os.W_OK):
        if os.access("./model/I.mol", os.W_OK):
            if os.access("./model/MA.mol", os.W_OK):
                if os.access("./model/MB.mol", os.W_OK):
                    if os.access("./model/R.mol", os.W_OK):
                        Loading += 1
    if Loading != 3:
        print('Missing the model files is!\n', flush=True)
        print(
            'The model files: C.mol(CTA[-S-C(Z)=S]), I.mol(Initiator[I*]), MA.mol(Solvophilic monomer), \
        MB.mol(Solvophobic monomer), R.mol(R-end in CTA)\n',
            flush=True)
    if Loading == 3:
        Relist = []  # reaction ponit
        NewMol = []  # molrcular list after connection
        Mol = RAFTcreat.solvemol()  # Model molecules
        Relist = RAFTcreat.rectpoint(Mol)
        f_sys = open('./control/system.json', )  # Read nthread
        system = json.load(f_sys)
        nthread = 4
        if 'CPU' in system:
            nthread = int(system['CPU'])
        else:
            print('The ntheard uses the default value(4cores)!', flush=True)

        for list in Relist:
            if list != []:
                Loading += 1
        if Loading == 8:
            NewMol = RAFTcreat.connectmol(Mol, Relist)

            # Step Two: built .itp files
            Mol2 = Chem.AddHs(Mol[2])
            Mol3 = Chem.AddHs(Mol[3])
            NewMol.append(Mol2)
            NewMol.append(Mol3)
            list = [
                'IC', 'RC', 'IMAC', 'RMAC', 'IMBC', 'RMBC', 'IMA2C', 'IMB2C',
                'RMA3C', 'RMB3C', 'RMAMB2C', 'RMA2MBC', 'MA', 'MB'
            ]

            Gstart = 0
            Istart = 0
            if os.access("./gaussian/IRQ.log", os.F_OK):
                if os.access("./gaussian/IRQ.log", os.R_OK):
                    f_log = open("./gaussian/IRQ.log", 'r')
                    line = f_log.readline()
                    lines = line.split()
                    if lines[0] == 'Gaussian':
                        Gstart = int(lines[3])
                    if lines[0] == 'Itp':
                        Istart = int(lines[3])
                    f_log.close()

            if not os.access("./gaussian/__Gaussian__", os.F_OK):
                print("|",
                      "Gaussian caculating".center(40, "-"),
                      "|",
                      flush=True)
                for i in range(Gstart, len(list)):
                    # structure optimization
                    for j in range(3):
                        if j == 0:
                            AllChem.EmbedMolecule(NewMol[i],
                                                  randomSeed=0xf00d,
                                                  useRandomCoords=True)
                        if j == 1:
                            AllChem.EmbedMolecule(NewMol[i],
                                                  randomSeed=0xf00d,
                                                  useRandomCoords=True)
                            AllChem.MMFFOptimizeMolecule(NewMol[i])
                        if j == 2:
                            AllChem.EmbedMolecule(NewMol[i],
                                                  randomSeed=0xf00d,
                                                  useRandomCoords=True)
                            AllChem.MMFFOptimizeMolecule(NewMol[i],
                                                         maxIters=10000)

                        RAFTgaussian.gaussianfile(NewMol[i], list[i], 30,
                                                  nthread)
                        flow1 = RAFTgaussian.gaussiancalculate(list[i])
                        if flow1:
                            subprocess.Popen('rm -f ' + '*.log',
                                             shell=True,
                                             cwd='./gaussian')
                            subprocess.Popen('rm -f ' + 'fort.7',
                                             shell=True,
                                             cwd='./gaussian')
                            subprocess.Popen('rm -f ' + '*.gjf',
                                             shell=True,
                                             cwd='./gaussian')
                            subprocess.Popen('rm -f ' + list[i] + 'em.chk',
                                             shell=True,
                                             cwd='./gaussian')
                            subprocess.Popen('rm -f ' + list[i] + 'step1.chk',
                                             shell=True,
                                             cwd='./gaussian')
                            break
                        else:
                            subprocess.Popen('rm -f ' + list[i] + 'em.chk',
                                             shell=True,
                                             cwd='./gaussian')
                            subprocess.Popen('rm -f ' + list[i] + 'step1.chk',
                                             shell=True,
                                             cwd='./gaussian')
                            if j == 2:
                                f_log = open("./gaussian/IRQ.log", 'w')
                                print('Gaussian calculate sample ' + str(i) +
                                      ' (' + list[i] + ') is failed',
                                      file=f_log)
                                f_log.close()
                                cout += 1
                                break
            if cout == 0:
                f_G = open("./gaussian/__Gaussian__", 'w')
                f_G.close()

            if not os.access("./gaussian/__Itp__", os.F_OK) and cout == 0:
                print("|",
                      "Force field caculating".center(40, "-"),
                      "|",
                      flush=True)
                for i in range(Istart, len(list)):
                    flow2 = RAFTgaussian.itpcalculate(list[i])
                    if flow2:
                        subprocess.Popen('rm -f ' + '*.top',
                                         shell=True,
                                         cwd='./gaussian')
                        subprocess.Popen('rm -f ' + '*.sh',
                                         shell=True,
                                         cwd='./sobtop')
                        subprocess.Popen('rm -f ' + list[i] + '*step2.chk',
                                         shell=True,
                                         cwd='./gaussian')
                        subprocess.Popen('rm -f ' + '*.mol2',
                                         shell=True,
                                         cwd='./gaussian')
                    else:
                        f_log = open("./gaussian/IRQ.log", 'w')
                        print('Itp calculate sample ' + str(i) + ' (' +
                              list[i] + ') is failed',
                              file=f_log)
                        f_log.close()
                        cout += 1
                        break
            if cout == 0:
                f_I = open("./gaussian/__Itp__", 'w')
                f_I.close()

            if not os.access("./gaussian/__State__", os.F_OK) and cout == 0:
                flow3 = RAFTgaussian.chargecalculate()
                if flow3 is True:
                    subprocess.Popen('rm -f ' + '*.sh',
                                     shell=True,
                                     cwd='./gaussian')
                    subprocess.Popen('rm -f ' + '*.txt',
                                     shell=True,
                                     cwd='./gaussian')
                if flow3 is None:
                    cout += 1
                if flow3 is False:
                    cout += 1
            if cout == 0:
                f_C = open("./gaussian/__State__", 'w')
                f_C.close()

            if cout == 0:
                print("|", "Link1 Success".center(40, "-"), "|", flush=True)
            else:
                print("|", "Link1 Fail".center(40, "-"), "|", flush=True)
        else:
            print("Can't recognize reaction point", flush=True)
            print("|", "Link1 Fail".center(40, "-"), "|", flush=True)
    if Loading == 8 and cout == 0:
        return True
    else:
        return False
