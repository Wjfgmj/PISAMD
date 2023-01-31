import os
import subprocess
from bin.force import RAFTgromacs


def go():

    step1 = False
    step2 = False
    if os.access("./gromacs/__step4__", os.F_OK):
        if os.access("./gromacs/__step5__", os.F_OK):
            print("|", "Link6 Pass".center(40, "-"), "|", flush=True)
            return True
    print("|", "Link6".center(40, "-"), "|", flush=True)
    print("|", "MD modeling".center(40, "-"), "|", flush=True)
    if not os.access("./gromacs/__step4__", os.F_OK):
        print("|", "Step4".center(40, "-"), "|", flush=True)
        flow1 = RAFTgromacs.step4_1()
        if flow1:
            flow2 = RAFTgromacs.step4_2()
            if flow2 is not False:
                step1 = RAFTgromacs.step4_3(flow2)
                if step1:
                    f_step = open('./gromacs/__step4__', 'w')
                    f_step.close()
                    subprocess.Popen('rm -f *.edr',
                                     shell=True,
                                     cwd='./gromacs/step4')
                    subprocess.Popen('rm -f *.log',
                                     shell=True,
                                     cwd='./gromacs/step4')
                    subprocess.Popen('rm -f *.tpr',
                                     shell=True,
                                     cwd='./gromacs/step4')
                    subprocess.Popen('rm -f *.trr',
                                     shell=True,
                                     cwd='./gromacs/step4')
                    subprocess.Popen('rm -f *.cpt',
                                     shell=True,
                                     cwd='./gromacs/step4')
                    subprocess.Popen('rm -f *.mdp',
                                     shell=True,
                                     cwd='./gromacs/step4')
                    subprocess.Popen('rm -f *.pdb',
                                     shell=True,
                                     cwd='./gromacs/step4')
                    subprocess.Popen('rm -f *.top',
                                     shell=True,
                                     cwd='./gromacs/step4')
                else:
                    print('MD modeling Step4 is failed!', flush=True)
                    print("|", "Link6 Fail".center(40, "-"), "|", flush=True)
                    return False
    else:
        print("|", "Step4 pass".center(40, "-"), "|", flush=True)
        step1 = True

    print("|", "Step5".center(40, "-"), "|", flush=True)
    flow3 = RAFTgromacs.step5_1()
    if flow3 is not False:
        step2 = RAFTgromacs.step5_2(flow3)
        if step2:
            f_step = open('./gromacs/__step5__', 'w')
            f_step.close()
            subprocess.Popen('rm -f *.edr', shell=True, cwd='./gromacs/step5')
            subprocess.Popen('rm -f *.log', shell=True, cwd='./gromacs/step5')
            subprocess.Popen('rm -f *.tpr', shell=True, cwd='./gromacs/step5')
            subprocess.Popen('rm -f *.trr', shell=True, cwd='./gromacs/step5')
            subprocess.Popen('rm -f *.cpt', shell=True, cwd='./gromacs/step5')
            subprocess.Popen('rm -f *.mdp', shell=True, cwd='./gromacs/step5')
            subprocess.Popen('rm -f *.pdb', shell=True, cwd='./gromacs/step5')
            print("|", "Link6 Success".center(40, "-"), "|", flush=True)
            return True
        else:
            print('MD modeling Step5 is failed!', flush=True)
            print("|", "Link6 Fail".center(40, "-"), "|", flush=True)
            return False
