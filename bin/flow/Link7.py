import os
import json
import glob
import subprocess
from bin.force import RAFTgromacs, RAFT
from bin.jointmol import RAFTcreat, RAFTsystem


def go():

    if os.access("./gromacs/__step6__", os.F_OK):
        print("|", "Link7 Pass".center(40, "-"), "|", flush=True)
        return True
    print("|", "Link7".center(40, "-"), "|", flush=True)
    print("|", "PISA (self-assembly)".center(40, "-"), "|", flush=True)
    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Boxsize = system['Boxsize']
    reac = []
    conversion = 0.0
    cout = 0
    frame = 0
    if os.access("./gromacs/step6/", os.F_OK):
        files = glob.glob('./gromacs/step6/[0-9]*.gro')
        if len(files) == 1:
            filename = os.path.splitext(os.path.split(files[0])[-1])[0]
            frame = int(filename)
            print('The MD simulation will start at ' + filename + 'th frame!',
                  flush=True)
        if len(files) >= 2:
            filename = []
            for i in range(len(files)):
                temp = os.path.splitext(os.path.split(files[i])[-1])[0]
                filename.append(int(temp))
            frame = max(filename)
            print('The MD simulation will start at ' + str(frame) +
                  'th frame!',
                  flush=True)
    if frame == 0:
        RAFTgromacs.step6_1()

    Ilist, Rlist, Clist, MAlist, MBlist = RAFTcreat.PISA_Mollist(
    )  # mid, aid_start, aid_end (int)
    reactpotlist = RAFTcreat.reactpot_typelist('step6')
    RAFTsystem.PISA_grotodata(frame)

    if frame != 0:
        fin_gro = open('./gromacs/step6/' + str(frame) + '.gro', 'r')
        line = fin_gro.readline()
        line = fin_gro.readline()
        cout_MOB = 0
        up_molid = '0'
        while True:
            line = fin_gro.readline()
            Molid = line[:5].strip()
            Molname = line[5:10].strip()
            atomname = line[10:15].strip()
            if Molname == 'MOB' and Molid != up_molid:
                cout_MOB += 1
                up_molid = Molid
                if atomname == 'm_b':
                    cout += 1
            if cout_MOB == len(MBlist):
                break
        conversion = round(cout / len(MBlist) * 100, 2)
        fin_gro.close()

    while True:

        md = False
        atom = []
        f_sys = open('./control/system.json', )
        system = json.load(f_sys)
        con_MB = 0.95
        if 'Conversion' in system:
            Conver = system['Conversion']
        else:
            print('The Conversion uses the default value! (C_MA & C_MB=0.95)')
        c_MB = Conver[0].get('MB')
        if c_MB is None:
            print('The Conversion uses the default value! (C_MB=0.95)')
        if c_MB is not None:
            con_MB = c_MB
        f_sys.close()
        print(str(frame) + 'th step is running, ' + str(conversion) +
              '% monomers were converted',
              end='\n',
              flush=True)
        if conversion <= con_MB * 100:
            reactlist = RAFTcreat.reactpot_poslist('step6', frame,
                                                   reactpotlist)
            reac = RAFT.Find(
                reactlist, Boxsize
            )  # atoma_id atomb_id mola_id molb_id atoma_name atomb_name (str)
            real_reac = []
            if reac != []:
                real_reac = []
                if len(reac) > 1:
                    for i in range(len(reac) - 1):
                        real = True
                        for j in range(i + 1, len(reac)):
                            if reac[i][0] == reac[j][0] or \
                               reac[i][1] == reac[j][1] or \
                               reac[i][1] == reac[j][0] or \
                               reac[i][0] == reac[j][1]:
                                real = False
                                break
                        if real:
                            real_reac.append(reac[i])
                    real_reac.append(reac[-1])
                    reac = real_reac
                for i in range(len(reac)):
                    if 'm_n' in reac[i]:
                        cout += 1
                conversion = round(cout / len(MBlist) * 100, 2)
                atom = RAFTsystem.PISA_changeitp(Ilist, Rlist, Clist, MAlist,
                                                 MBlist, reac)
                RAFTsystem.PISA_changegro(frame, atom)
                reactpotlist = RAFTcreat.reactpot_typelist('step6')
            md = RAFTgromacs.step6_2(frame, reac)
            if md:
                RAFTsystem.PISA_grotodata(frame + 1)
                subprocess.Popen('rm -f *.edr',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f *.log',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f *.tpr',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f *.trr',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f *.cpt',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f *.pdb',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f *.1#',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f mdout.mdp',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f em' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f eq' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step6')
                subprocess.Popen('rm -f ' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step6')
                frame += 1
            else:
                print('The MD simulation is wrong!', flush=True)
                print("|", "Link7 Fail".center(40, "-"), "|", flush=True)
                return False
        else:
            f_step = open('./gromacs/__step6__', 'w')
            f_step.close()
            print("|", "Link7 Success".center(40, "-"), "|", flush=True)
            return True
