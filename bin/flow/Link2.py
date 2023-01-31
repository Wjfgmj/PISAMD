import os
import json
import subprocess
from bin.jointmol import RAFTcreat, RAFTsystem, Normalsystem
from bin.force import Normalgaussian


def go():

    if not os.access("./topper/__Done__", os.F_OK):
        cout = 0
        print("|", "Link2".center(40, "-"), "|", flush=True)
        print("|", "Building force file".center(40, "-"), "|", flush=True)
        RAFTcreat.GAFFforcefile()
        RAFTcreat.statefile()
        print("|", "Building polymer file".center(40, "-"), "|", flush=True)
        flow = RAFTsystem.RAFTsteponefile()
        if not flow:
            print('Building polymer file is failed!', flush=True)
            cout += 1
        print("|", "Building solvent file".center(40, "-"), "|", flush=True)
        f_sys = open('./control/system.json', )
        system = json.load(f_sys)
        if 'Solvent_CTA' in system:
            Solvent = system['Solvent_CTA']
            if Solvent[0] != []:
                for i in range(len(Solvent[0])):
                    if os.access("./model/" + str(Solvent[0][i]) + ".mol",
                                 os.F_OK):
                        flow0 = Normalsystem.waterMDfile(str(Solvent[0][i]))
                        if not flow0:
                            Normalgaussian.gaussianfile(str(Solvent[0][i]),
                                                        mem=30)
                            flow1 = Normalgaussian.gaussiancalculate(
                                str(Solvent[0][i]))
                            if flow1:
                                subprocess.Popen('rm -f ' + '*.log',
                                                 shell=True,
                                                 cwd='./topper')
                                subprocess.Popen('rm -f ' + 'fort.7',
                                                 shell=True,
                                                 cwd='./topper')
                                subprocess.Popen('rm -f ' + '*.gjf',
                                                 shell=True,
                                                 cwd='./topper')
                                subprocess.Popen('rm -f ' +
                                                 str(Solvent[0][i]) + 'em.chk',
                                                 shell=True,
                                                 cwd='./topper')
                                subprocess.Popen('rm -f ' +
                                                 str(Solvent[0][i]) +
                                                 'step1.chk',
                                                 shell=True,
                                                 cwd='./topper')
                            else:
                                cout += 1
                                break
                            flow2 = Normalgaussian.itpcalculate(
                                str(Solvent[0][i]))
                            if flow2:
                                subprocess.Popen('rm -f ' + '*.top',
                                                 shell=True,
                                                 cwd='./topper')
                                subprocess.Popen('rm -f ' + '*.sh',
                                                 shell=True,
                                                 cwd='./sobtop')
                                subprocess.Popen('rm -f ' +
                                                 str(Solvent[0][i]) +
                                                 '*step2.chk',
                                                 shell=True,
                                                 cwd='./topper')
                                subprocess.Popen('rm -f ' + '*.mol2',
                                                 shell=True,
                                                 cwd='./topper')
                            else:
                                cout += 1
                                break
                            flow3 = Normalgaussian.chargecalculate(
                                str(Solvent[0][i]))
                            if flow3:
                                subprocess.Popen('rm -f ' + '*.sh',
                                                 shell=True,
                                                 cwd='./topper')
                                subprocess.Popen('rm -f ' + '*.txt',
                                                 shell=True,
                                                 cwd='./topper')
                                subprocess.Popen('rm -f ' + '*.fchk',
                                                 shell=True,
                                                 cwd='./topper')
                            else:
                                cout += 1
                                break

                            if cout == 0:
                                Normalsystem.MDfile(str(Solvent[0][i]))
                    else:
                        print('There is not the solvent file in <model> fold!',
                              flush=True)
                        cout += 1
        else:
            print('The solvent is missing!', flush=True)
            cout += 1

        if 'Ion_CTA' in system:
            print("|", "Building ion file".center(40, "-"), "|", flush=True)
            Ion = system['Ion_CTA']
            if Ion[0] != []:
                for i in range(len(Ion[0])):
                    if os.access("./model/" + str(Ion[0][i]) + ".mol",
                                 os.F_OK):
                        flow4 = Normalsystem.IonMDfile(str(Ion[0][i]))
                        if flow4:
                            subprocess.Popen('rm -f ' + '*.log',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + 'fort.7',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + '*.gjf',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + str(Ion[0][i]) +
                                             'em.chk',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + str(Ion[0][i]) +
                                             'step1.chk',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + str(Ion[0][i]) +
                                             'step2.chk',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + str(Ion[0][i]) +
                                             'step2.fchk',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + '*.top',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + '*.sh',
                                             shell=True,
                                             cwd='./sobtop')
                            subprocess.Popen('rm -f ' + '*.sh',
                                             shell=True,
                                             cwd='./topper')
                            subprocess.Popen('rm -f ' + '*.mol2',
                                             shell=True,
                                             cwd='./topper')
                        else:
                            cout += 1
                            break
                    else:
                        print('There is not the ion file in <model> fold!',
                              flush=True)
                        cout += 1

        f_sys.close()
        if cout == 0:
            f_topper = open("./topper/__Done__", 'w')
            print("|", "Link2 Success".center(40, "-"), "|", flush=True)
            f_topper.close()
            fin_itp = open('./topper/ffnonbonded.itp', 'r')
            line = fin_itp.readlines()
            lines = sorted(set(line), key=line.index)
            fin_itp.close()
            fout_itp = open('./topper/ffnonbonded.itp', 'w')
            for i in range(len(lines)):
                print(lines[i], end='', file=fout_itp)
            fout_itp.close()
            return True
        else:
            print("|", "Link2 Fail".center(40, "-"), "|", flush=True)
            return False
    else:
        print("|", "Link2 Pass".center(40, "-"), "|", flush=True)
        return True
