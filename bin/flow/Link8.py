import os
import glob
from bin.force import RAFTgromacs


def go():

    if os.access("./gromacs/__step7__", os.F_OK):
        print("|", "Link8 Pass".center(40, "-"), "|", flush=True)
        return True
    print("|", "Link8".center(40, "-"), "|", flush=True)
    print("|", "PISA (equilibrium)".center(40, "-"), "|", flush=True)
    if os.access("./gromacs/step6/", os.F_OK):
        files = glob.glob('./gromacs/step6/[0-9]*.gro')
        if len(files) == 1:
            filename = os.path.splitext(os.path.split(files[0])[-1])[0]
            step = RAFTgromacs.step7_1(filename)
            if step:
                f_step = open('./gromacs/__step7__', 'w')
                f_step.close()
                print("|", "Link8 Success".center(40, "-"), "|", flush=True)
                return True
            else:
                print('The MD simulation is wrong!', flush=True)
                return False
        else:
            print('The PISA process does not end normally!', flush=True)
            return False
