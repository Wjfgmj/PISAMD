import os
import json
import glob
import subprocess
from bin.force import RAFTgromacs, RAFT
from bin.jointmol import RAFTcreat, RAFTsystem


def go():

    if os.access("./gromacs/__step3__", os.F_OK):
        print("|", "Link4 Pass".center(40, "-"), "|", flush=True)
        return True
    print("|", "Link4".center(40, "-"), "|", flush=True)
    print("|", "PISA (prepare CTA)".center(40, "-"), "|", flush=True)
    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Boxsize = system['Boxsize']
    reac = []
    conversion = 0.0
    bond_CTA = 0
    CTAlist = []
    cout = 0
    frame = 0
    if os.access("./gromacs/step3/", os.F_OK):
        files = glob.glob('./gromacs/step3/[0-9]*.gro')
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
        RAFTgromacs.step3_1()

    Ilist, Rlist, Clist, MAlist = RAFTcreat.CTA_Mollist(
    )  # mid, aid_start, aid_end (int)
    reactpotlist = RAFTcreat.reactpot_typelist('step3')
    RAFTsystem.CTA_grotodata(frame)
    nCTA = len(Rlist) + len(Ilist)

    if frame != 0:
        fin_gro = open('./gromacs/step3/' + str(frame) + '.gro', 'r')
        line = fin_gro.readline()
        line = fin_gro.readline()
        cout_MOA = 0
        up_molid = '0'
        while True:
            line = fin_gro.readline()
            Molid = line[:5].strip()
            Molname = line[5:10].strip()
            atomname = line[10:15].strip()
            if Molname == 'CTA' and atomname == 'c_b':
                bond_CTA += 1
            if Molname == 'MOA' and Molid != up_molid:
                cout_MOA += 1
                up_molid = Molid
                if atomname == 'm_b':
                    cout += 1
            if cout_MOA == len(MAlist):
                break
        conversion = round(cout / len(MAlist) * 100, 2)
        fin_gro.close()

    while True:

        md = False
        atom = []
        f_sys = open('./control/system.json', )
        system = json.load(f_sys)
        con_MA = 0.95
        if 'Conversion' in system:
            Conver = system['Conversion']
        else:
            print('The Conversion uses the default value! (C_MA=0.95)')
        c_MA = Conver[0].get('MA')
        if c_MA is None:
            print('The Conversion uses the default value! (C_MA=0.95)')
        if c_MA is not None:
            con_MA = c_MA
        f_sys.close()

        if conversion <= con_MA * 100:
            print(str(frame) + 'th step is running, ' + str(conversion) +
                  '% monomers were polymerized',
                  end='\n',
                  flush=True)
            reactlist = RAFTcreat.reactpot_poslist('step3', frame,
                                                   reactpotlist)
            reac = RAFT.Find(
                reactlist, Boxsize
            )  # atoma_id atomb_id mola_id molb_id atoma_name atomb_name (str)
            real_reac = []
            for i in range(len(reac)):
                real = True
                if reac[i][5] == 'c_n' or reac[i][4] == 'c_n':
                    real = False
                if real:
                    real_reac.append(reac[i])
            reac = real_reac
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

                conversion = round(cout / len(MAlist) * 100, 2)
                atom = RAFTsystem.CTA_changeitp(Ilist, Rlist, Clist, MAlist,
                                                reac)
                RAFTsystem.CTA_changegro(frame, atom)
                reactpotlist = RAFTcreat.reactpot_typelist('step3')
            md = RAFTgromacs.step3_2(frame, reac)
            if md:
                RAFTsystem.CTA_grotodata(frame + 1)
                subprocess.Popen('rm -f *.edr',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.log',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.tpr',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.trr',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.cpt',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.pdb',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.1#',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f mdout.mdp',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f em' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f eq' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f ' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step3')
                frame += 1
            else:
                print('The MD simulation is wrong!', flush=True)
                print("|", "Link4 Fail".center(40, "-"), "|", flush=True)
                return False
        if conversion > con_MA * 100 and bond_CTA != nCTA:
            Con_CTA = round(bond_CTA / nCTA * 100, 2)
            print(str(frame) + 'th step is running, ' + str(Con_CTA) +
                  '% macroCTA were prepared',
                  end='\n',
                  flush=True)
            reactlist = RAFTcreat.reactpot_poslist('step3', frame,
                                                   reactpotlist)
            reac = RAFT.Find(
                reactlist, Boxsize
            )  # atoma_id atomb_id mola_id molb_id atoma_name atomb_name (str)
            real_reac = []
            for i in range(len(reac)):
                real = False
                if reac[i][5] == 'c_n':
                    if reac[i][3] not in CTAlist:
                        real = True
                        CTAlist.append(reac[i][3])
                if reac[i][4] == 'c_n':
                    if reac[i][2] not in CTAlist:
                        real = True
                        CTAlist.append(reac[i][2])
                if real:
                    real_reac.append(reac[i])
            reac = real_reac
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
                bond_CTA += len(reac)
                atom = RAFTsystem.CTA_changeitp(Ilist, Rlist, Clist, MAlist,
                                                reac)
                RAFTsystem.CTA_changegro(frame, atom)
                reactpotlist = RAFTcreat.reactpot_typelist('step3')
            md = RAFTgromacs.step3_2(frame, reac)
            if md:
                RAFTsystem.CTA_grotodata(frame + 1)
                subprocess.Popen('rm -f *.edr',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.log',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.tpr',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.trr',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.cpt',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.pdb',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f *.1#',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f mdout.mdp',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f em' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f eq' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step3')
                subprocess.Popen('rm -f ' + str(frame) + '.gro',
                                 shell=True,
                                 cwd='./gromacs/step3')
                frame += 1
            else:
                print('The MD simulation is wrong!', flush=True)
                print("|", "Link4 Fail".center(40, "-"), "|", flush=True)
                return False
        if conversion > con_MA * 100 and bond_CTA == nCTA:
            f_step = open('./gromacs/__step3__', 'w')
            f_step.close()
            print("|", "Link4 Success".center(40, "-"), "|", flush=True)
            return True
