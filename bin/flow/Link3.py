import os
import subprocess
from bin.force import RAFTgromacs


def go():

    step1 = False
    step2 = False
    if os.access("./gromacs/__step1__", os.F_OK):
        if os.access("./gromacs/__step2__", os.F_OK):
            print("|", "Link3 Pass".center(40, "-"), "|", flush=True)
            return True
    print("|", "Link3".center(40, "-"), "|", flush=True)
    print("|", "MD modeling".center(40, "-"), "|", flush=True)
    if not os.access("./gromacs/__step1__", os.F_OK):
        print("|", "Step1".center(40, "-"), "|", flush=True)
        flow1 = RAFTgromacs.step1_1()
        if flow1:
            flow2 = RAFTgromacs.step1_2()
            if flow2 is not False:
                step1 = RAFTgromacs.step1_3(flow2)
                if step1:
                    f_step = open('./gromacs/__step1__', 'w')
                    f_step.close()
                    subprocess.Popen('rm -f *.edr',
                                     shell=True,
                                     cwd='./gromacs/step1')
                    subprocess.Popen('rm -f *.log',
                                     shell=True,
                                     cwd='./gromacs/step1')
                    subprocess.Popen('rm -f *.tpr',
                                     shell=True,
                                     cwd='./gromacs/step1')
                    subprocess.Popen('rm -f *.trr',
                                     shell=True,
                                     cwd='./gromacs/step1')
                    subprocess.Popen('rm -f *.cpt',
                                     shell=True,
                                     cwd='./gromacs/step1')
                    subprocess.Popen('rm -f *.mdp',
                                     shell=True,
                                     cwd='./gromacs/step1')
                    subprocess.Popen('rm -f *.pdb',
                                     shell=True,
                                     cwd='./gromacs/step1')
                    subprocess.Popen('rm -f *.top',
                                     shell=True,
                                     cwd='./gromacs/step1')
                else:
                    print('MD modeling Step1 is failed!', flush=True)
                    print("|", "Link3 Fail".center(40, "-"), "|", flush=True)
                    return False
    else:
        print("|", "Step1 pass".center(40, "-"), "|", flush=True)
        step1 = True

    print("|", "Step2".center(40, "-"), "|", flush=True)
    flow3 = RAFTgromacs.step2_1()
    if flow3 is not False:
        step2 = RAFTgromacs.step2_2(flow3)
        if step2:
            f_step = open('./gromacs/__step2__', 'w')
            f_step.close()
            subprocess.Popen('rm -f *.edr', shell=True, cwd='./gromacs/step2')
            subprocess.Popen('rm -f *.log', shell=True, cwd='./gromacs/step2')
            subprocess.Popen('rm -f *.tpr', shell=True, cwd='./gromacs/step2')
            subprocess.Popen('rm -f *.trr', shell=True, cwd='./gromacs/step2')
            subprocess.Popen('rm -f *.cpt', shell=True, cwd='./gromacs/step2')
            subprocess.Popen('rm -f *.mdp', shell=True, cwd='./gromacs/step2')
            subprocess.Popen('rm -f *.pdb', shell=True, cwd='./gromacs/step2')
            print("|", "Link3 Success".center(40, "-"), "|", flush=True)
            return True
        else:
            print('MD modeling Step2 is failed!', flush=True)
            print("|", "Link3 Fail".center(40, "-"), "|", flush=True)
            return False
