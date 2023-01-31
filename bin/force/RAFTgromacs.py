import os
import re
import json
import subprocess


def step1_1():

    cout = 0
    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    if 'PISA' in system:
        PISA = system['PISA']
        if 'CTA' not in PISA[0].keys():
            print('CTA is missing!')
            cout += 1
        if 'INI' not in PISA[0].keys():
            print('INI is missing!')
            cout += 1
        if 'MA' not in PISA[0].keys():
            print('MA is missing!')
            cout += 1
        if 'MB' not in PISA[0].keys():
            print('MB is missing!')
            cout += 1
    else:
        print('PISA system(mol/L) is missing')
        print(
            'Please add: "PISA": [{"CTA": C_CTA(Mol/L), "INI": C_INI(Mol/L), "MA": C_MA(Mol/L), "MB": C_MB(Mol/L)}]'
        )
        cout += 1

    if not os.access("./topper/ffbonded.itp", os.F_OK):
        print('ffbonded.itp is missing!')
        cout += 1
    if not os.access("./topper/ffnonbonded.itp", os.F_OK):
        print('ffnonbonded.itp is missing!')
        cout += 1
    if not os.access("./topper/forcefield.itp", os.F_OK):
        print('forcefield.itp is missing!')
        cout += 1
    if not os.access("./topper/PISA.itp", os.F_OK):
        print('PISA.itp is missing!')
        cout += 1

    if 'Solvent_CTA' in system:
        Solvent = system['Solvent_CTA']
        if len(Solvent) != 2:
            print(
                'Please add: "Solvent_CTA": [["filename_1", "filename_2"], [volume ratio_1, volume ratio_2]]'
            )
            cout += 1
        else:
            for file in Solvent[0]:
                if not os.access("./topper/" + str(file) + '.pdb', os.F_OK):
                    if not os.access("./topper/" + str(file) + '.itp',
                                     os.F_OK):
                        print('The files of solvent (' + file +
                              ') are not existed!')
                        cout += 1
            total = 0.0
            for ratio in Solvent[1]:
                total = total + ratio
            if int(total * 10) / 10 != 1.0:
                print('The sum of the volume ratios is not 1.0!')
                cout += 1
    else:
        print('The Solvent is missing!')
        print(
            'Please add: "Solvent_CTA": [["filename_1", "filename_2"], [volume ratio_1, volume ratio_2]]'
        )
        print('e.g. "Solvent_CTA": [["water", "alcohol"], [0.7, 0.3]]')
        cout += 1

    if 'Ion_CTA' in system:
        Ion = system['Ion_CTA']
        if len(Ion) != 2:
            print(
                'Please add: "Ion_CTA": [["filename_1", "filename_2"], [concentration_1 (Mol/L), concentration_2 (Mol/L)]]'
            )
            print('e.g. "Ion_CTA": [["Na", "Cl"], [0.153, 0.153]]')
            cout += 1
        else:
            for file in Ion[0]:
                if not os.access("./topper/" + str(file) + '.pdb', os.F_OK):
                    if not os.access("./topper/" + str(file) + '.itp',
                                     os.F_OK):
                        print('The files of ion (' + file +
                              ') are not existed!')
                        cout += 1
    f_sys.close()
    if cout == 0:
        return True
    else:
        return False


def step1_2():

    folder = os.path.exists('./gromacs/step1')
    if not folder:
        os.makedirs('./gromacs/step1')

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Solvent = system['Solvent_CTA']

    size = [30] * len(Solvent[0])
    for i in range(len(Solvent[0])):
        flow1 = subprocess.Popen('cp ./topper/' + str(Solvent[0][i]) +
                                 '.pdb ./gromacs/step1',
                                 shell=True)
        flow1.wait()
        f_packmol = open('./gromacs/step1/' + str(Solvent[0][i]) + '.inp', 'w')
        print('tolerance 2.4\n',
              'filetype pdb',
              'output  Sys_' + str(Solvent[0][i]) + '.pdb',
              'add_box_sides 1.2',
              '\nstructure ' + str(Solvent[0][i]) + '.pdb',
              '    number 500',
              '    inside box 0 0 0' + ' ' + str(size[i]) + ' ' +
              str(size[i]) + ' ' + str(size[i]),
              'end structure',
              sep='\n',
              file=f_packmol)
        f_packmol.close()

    cout = 0
    for i in range(len(Solvent[0])):
        pack = False
        while not pack:
            flow1 = subprocess.Popen('packmol < ' + str(Solvent[0][i]) +
                                     '.inp',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step1')
            try:
                (msg, errs) = flow1.communicate(timeout=120)
            except subprocess.TimeoutExpired:
                flow1.kill()
                flow1.terminate()
                f_packmol = open(
                    './gromacs/step1/' + str(Solvent[0][i]) + '.inp', 'w')
                size[i] = size[i] * 1.5
                print('tolerance 2.4\n',
                      'filetype pdb',
                      'output  Sys_' + str(Solvent[0][i]) + '.pdb',
                      'add_box_sides 1.2',
                      '\nstructure ' + str(Solvent[0][i]) + '.pdb',
                      '    number 500',
                      '    inside box 0 0 0' + ' ' + str(size[i]) + ' ' +
                      str(size[i]) + ' ' + str(size[i]),
                      'end structure',
                      sep='\n',
                      file=f_packmol)
                f_packmol.close()
            if flow1.returncode == 0:
                pack = True
                cout = 0
                subprocess.Popen('rm -f ' + str(Solvent[0][i]) + '.inp',
                                 shell=True,
                                 cwd='./gromacs/step1')
                subprocess.Popen('rm -f ' + str(Solvent[0][i]) + '.pdb',
                                 shell=True,
                                 cwd='./gromacs/step1')
            else:
                cout += 1
            if cout > 10:
                break

    f_sys.close()
    if cout > 10:
        return False
    else:
        return size


def step1_3(size):

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Solvent = system['Solvent_CTA']
    cout = 0
    CPU = 4
    GPU = 0
    GPU_id = ''
    if 'CPU' in system:
        CPU = system['CPU']
    else:
        print('The ntheard uses the default value! (CPU=4)')
    if 'GPU' in system:
        GPU = system['GPU']
    else:
        print('The number of GPU uses the default value! (GPU=0)')
    if GPU != 0:
        for id in range(int(GPU)):
            GPU_id = GPU_id + str(id)

    for i in range(len(Solvent[0])):
        f_em = open('./gromacs/step1/em' + str(i) + '.mdp', 'w')
        print('integrator              = steep',
              'emtol                   = 1000.0',
              'nsteps                  = 10000',
              'nstlist                 = 10',
              'cutoff-scheme           = Verlet',
              'rlist                   = 1.2',
              'vdwtype                 = Cut-off',
              'vdw-modifier            = Force-switch',
              'rvdw_switch             = 1.0',
              'rvdw                    = 1.2',
              'coulombtype             = pme',
              'rcoulomb                = 1.2',
              'constraints             = h-bonds',
              'constraint_algorithm    = LINCS',
              sep='\n',
              file=f_em)

        f_top = open('./gromacs/step1/' + str(Solvent[0][i]) + '.top', 'w')
        print('#include "../../topper/forcefield.itp"',
              '#include "../../topper/' + str(Solvent[0][i]) + '.itp”',
              '\n[ system ]',
              str(Solvent[0][i]),
              '\n[ molecules ]',
              str(Solvent[0][i]) + '     500',
              sep='\n',
              file=f_top)

        f_em.close()
        f_top.close()

        box = str('{:.5f}'.format(size[i] / 10))
        flow1 = subprocess.Popen('gmx editconf -f Sys_' + str(Solvent[0][i]) +
                                 '.pdb -o ' + str(Solvent[0][i]) +
                                 '.gro -box' + ' ' + box + ' ' + box + ' ' +
                                 box,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step1')
        flow1.wait()
        flow2 = subprocess.Popen('gmx grompp -f em' + str(i) + '.mdp -c ' +
                                 str(Solvent[0][i]) + '.gro -p ' +
                                 str(Solvent[0][i]) + '.top -o em_' +
                                 str(Solvent[0][i]) + '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step1')
        flow2.wait()
        if GPU == 0:
            flow3 = subprocess.Popen('gmx mdrun -v -deffnm em_' +
                                     str(Solvent[0][i]) + ' -ntmpi 1 -ntomp ' +
                                     str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step1')
            flow3.communicate()
            if flow3.returncode != 0:
                cout += 1
        else:
            flow3 = subprocess.Popen(
                'gmx mdrun -v -deffnm em_' + str(Solvent[0][i]) +
                ' -pin on -ntmpi 1 -ntomp ' + str(CPU) + ' -nb gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step1')
            flow3.communicate()
            if flow3.returncode != 0:
                cout += 1

    T = 300
    P = 1.0
    if 'Temperature' in system:
        T = system['Temperature']
    else:
        print('The Temperature uses the default value! (Temperature=300)')
    if 'Pressure' in system:
        P = system['Pressure']
    else:
        print('The Pressure used the default value! (Pressure=1.0)')

    for i in range(len(Solvent[0])):

        f_eq = open('./gromacs/step1/eq' + str(i) + '.mdp', 'w')
        print('integrator              = md',
              'dt                      = 0.001',
              'nsteps                  = 1000000',
              'nstlog                  = 0',
              'nstcalcenergy           = 0',
              'nstenergy               = 0',
              'cutoff-scheme           = Verlet',
              'nstlist                 = 20',
              'rlist                   = 1.2',
              'coulombtype             = pme',
              'rcoulomb                = 1.2',
              'vdwtype                 = Cut-off',
              'vdw-modifier            = Force-switch',
              'rvdw_switch             = 1.0',
              'rvdw                    = 1.2',
              'tcoupl                  = berendsen',
              'tc_grps                 = SYSTEM',
              'tau_t                   = 1.0',
              'ref_t                   = ' + str(T),
              'pcoupl                  = berendsen',
              'pcoupltype              = isotropic',
              'tau_p                   = 5.0',
              'compressibility         = 4.5e-5',
              'ref_p                   = ' + str(P),
              'constraints             = h-bonds',
              'constraint_algorithm    = LINCS',
              'continuation            = no',
              'nstcomm                 = 100',
              'comm_mode               = Linear',
              'comm_grps               = SYSTEM',
              'gen-vel                 = yes',
              'gen-temp                = ' + str(T),
              'gen-seed                = -1',
              sep='\n',
              file=f_eq)
        f_eq.close()

        flow4 = subprocess.Popen('gmx grompp -f eq' + str(i) + '.mdp -c em_' +
                                 str(Solvent[0][i]) + '.gro -p ' +
                                 str(Solvent[0][i]) + '.top -o eq_' +
                                 str(Solvent[0][i]) + '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step1')
        flow4.wait()
        if GPU == 0:
            flow5 = subprocess.Popen('gmx mdrun -v -deffnm eq_' +
                                     str(Solvent[0][i]) + ' -ntmpi 1 -ntomp ' +
                                     str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step1')
            flow5.communicate()
            if flow5.returncode != 0:
                cout += 1
        else:
            flow5 = subprocess.Popen(
                'gmx mdrun -v -deffnm eq_' + str(Solvent[0][i]) +
                ' -pin on -ntmpi 1 -ntomp ' + str(CPU) + ' -nb gpu -pme gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step1')
            flow5.communicate()
            if flow5.returncode != 0:
                cout += 1

    for i in range(len(Solvent[0])):
        f_md = open('./gromacs/step1/md' + str(i) + '.mdp', 'w')
        print('integrator              = md',
              'dt                      = 0.001',
              'nsteps                  = 4000000',
              'nstlog                  = 0',
              'nstcalcenergy           = 0',
              'nstenergy               = 0',
              'cutoff-scheme           = Verlet',
              'nstlist                 = 20',
              'rlist                   = 1.2',
              'coulombtype             = pme',
              'rcoulomb                = 1.2',
              'vdwtype                 = Cut-off',
              'vdw-modifier            = Force-switch',
              'rvdw_switch             = 1.0',
              'rvdw                    = 1.2',
              'tcoupl                  = Nose-Hoover',
              'tc_grps                 = SYSTEM',
              'tau_t                   = 1.0',
              'ref_t                   = ' + str(T),
              'pcoupl                  = Parrinello-Rahman',
              'pcoupltype              = isotropic',
              'tau_p                   = 5.0',
              'compressibility         = 4.5e-5',
              'ref_p                   = ' + str(P),
              'constraints             = h-bonds',
              'constraint_algorithm    = LINCS',
              'continuation            = no',
              'nstcomm                 = 100',
              'comm_mode               = Linear',
              'comm_grps               = SYSTEM',
              sep='\n',
              file=f_md)
        f_md.close()

        flow6 = subprocess.Popen('gmx grompp -f md' + str(i) + '.mdp -c eq_' +
                                 str(Solvent[0][i]) + '.gro -p ' +
                                 str(Solvent[0][i]) + '.top -o md_' +
                                 str(Solvent[0][i]) + '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step1')
        flow6.wait()
        if GPU == 0:
            flow7 = subprocess.Popen('gmx mdrun -v -deffnm md_' +
                                     str(Solvent[0][i]) + ' -ntmpi 1 -ntomp ' +
                                     str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step1')
            flow7.communicate()
            if flow7.returncode != 0:
                cout += 1
        else:
            flow7 = subprocess.Popen(
                'gmx mdrun -v -deffnm md_' + str(Solvent[0][i]) +
                ' -pin on -ntmpi 1 -ntomp ' + str(CPU) + ' -nb gpu -pme gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step1')
            flow7.communicate()
            if flow7.returncode != 0:
                cout += 1

    if cout != 0:
        return False
    else:
        return True


def step2_1():

    folder = os.path.exists('./gromacs/step2')
    if not folder:
        os.makedirs('./gromacs/step2')

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Solvent = system['Solvent_CTA']
    boxsize = system['Boxsize']
    polymer = system['PISA']
    Pini = 0
    if 'Pini_in_macroCTA' in system:
        Pini = system['Pini_in_macroCTA']

    Ion = []

    for i in range(len(Solvent[0])):
        flow1 = subprocess.Popen('cp ./gromacs/step1/md_' +
                                 str(Solvent[0][i]) + '.gro ./gromacs/step2',
                                 shell=True)
        flow1.wait()

    flow2 = subprocess.Popen('cp ./topper/I.pdb ./gromacs/step2', shell=True)
    flow2.wait()
    flow3 = subprocess.Popen('cp ./topper/R.pdb ./gromacs/step2', shell=True)
    flow3.wait()
    flow4 = subprocess.Popen('cp ./topper/C.pdb ./gromacs/step2', shell=True)
    flow4.wait()
    flow5 = subprocess.Popen('cp ./topper/MA.pdb ./gromacs/step2', shell=True)
    flow5.wait()

    x, y, z = boxsize[0], boxsize[1], boxsize[2]
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
    f_packmol = open('./gromacs/step2/system.inp', 'w')
    print('tolerance 2.0\n',
          'filetype pdb',
          'output  system.pdb',
          'add_box_sides 1.2',
          '\nstructure I.pdb',
          '    number ' + str(nI),
          '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' + str(y * 10) +
          ' ' + str(z * 10),
          'end structure',
          '\nstructure R.pdb',
          '    number ' + str(nR),
          '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' + str(y * 10) +
          ' ' + str(z * 10),
          'end structure',
          '\nstructure C.pdb',
          '    number ' + str(nC),
          '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' + str(y * 10) +
          ' ' + str(z * 10),
          'end structure',
          '\nstructure MA.pdb',
          '    number ' + str(nMA),
          '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' + str(y * 10) +
          ' ' + str(z * 10),
          'end structure',
          sep='\n',
          file=f_packmol)

    if 'Ion_CTA' in system:
        Ion = system['Ion_CTA']
        for i in range(len(Ion[0])):
            flow2 = subprocess.Popen('cp ./topper/' + str(Ion[0][i]) +
                                     '.pdb ./gromacs/step2',
                                     shell=True)
            flow2.wait()
            nIon = int(0.6023 * x * y * z * Ion[1][i] + 0.5)
            print('\nstructure ' + str(Ion[0][i]) + '.pdb',
                  '    number ' + str(nIon),
                  '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' +
                  str(y * 10) + ' ' + str(z * 10),
                  'end structure',
                  sep='\n',
                  file=f_packmol)

    f_packmol.close()

    flow6 = subprocess.Popen('packmol < system.inp',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True,
                             cwd='./gromacs/step2')
    try:
        (msg, errs) = flow6.communicate(timeout=180)
    except subprocess.TimeoutExpired:
        flow6.kill()
        flow6.terminate()
    if flow6.returncode != 0:
        print('Packmol is failed, please check the .inp file!')
        return False
    else:
        subprocess.Popen('rm -f system.inp', shell=True, cwd='./gromacs/step2')

    flow7 = subprocess.Popen('gmx editconf -f system.pdb -o system.gro -box' +
                             ' ' + str(x) + ' ' + str(y) + ' ' + str(z),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True,
                             cwd='./gromacs/step2')
    flow7.wait()

    N_sol = [1] * len(Solvent[0])
    for i in range(len(Solvent[0])):
        Nms = 0
        if i == 0:
            flow8 = subprocess.Popen(
                'gmx solvate -cp system.gro -cs md_' + str(Solvent[0][i]) +
                '.gro -o add_' + str(i) + '.gro -box' + ' ' + str(x) + ' ' +
                str(y) + ' ' + str(z),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step2')
            stdout, _ = flow8.communicate()
            out = stdout.decode('utf-8')
            Num = re.search(r'Number of solvent molecules:.*[0-9]', out)
            Nms = int(Num.group().split(':')[1])
            if flow8.returncode != 0:
                print('Add solvent' + str(Solvent[0][i]) + ' is failed!')
                return False
        if i > 0:
            flow8 = subprocess.Popen(
                'gmx solvate -cp add_' + str(i - 1) + '.gro -cs md_' +
                str(Solvent[0][i]) + '.gro -o add_' + str(i) + '.gro -box' +
                ' ' + str(x) + ' ' + str(y) + ' ' + str(z),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step2')
            stdout, _ = flow8.communicate()
            out = stdout.decode('utf-8')
            Num = re.search(r'Number of solvent molecules:.*[0-9]', out)
            Nms = int(Num.group().split(':')[1])
            if flow8.returncode != 0:
                print('Add solvent' + str(Solvent[0][i]) + ' is failed!')
                return False
            if i == len(Solvent[0]) - 1:
                N_sol[i] = Nms
        des = 0
        f_gro = open('./gromacs/step1/' + str(Solvent[0][i]) + '.gro', 'r')
        line = f_gro.readline()
        line = f_gro.readline()
        Natom = int(int(line) / 500)
        f_gro.close()
        if i == 0:
            des = int((1.0 - Solvent[1][i]) * Nms)
            if des == 0 and len(N_sol) > 1:
                des = 1
            N_sol[i] = Nms - des
            fr_sol = open('./gromacs/step2/add_' + str(i) + '.gro', 'r')
            lines = fr_sol.readlines()
            grofile = [
                lines[j] for j in range(len(lines))
                if (j < (len(lines) - 1 - des * Natom) or j == len(lines) - 1)
            ]
            grofile[1] = str(int(grofile[1]) - des * Natom) + '\n'
            fw_sol = open('./gromacs/step2/add_' + str(i) + '.gro', 'w')
            for k in range(len(grofile)):
                print(grofile[k], end='', file=fw_sol)
            fr_sol.close()
            fw_sol.close()
        if i > 0 and i != len(Solvent[0]) - 1:
            pde = 0
            for j in range(i + 1):
                pde = pde + Solvent[1][j]
            des = int((1 - pde) / Solvent[1][i] * Nms)
            if des == 0:
                des = 1
            N_sol[i] = Nms - des
            fr_sol = open('./gromacs/step2/add_' + str(i) + '.gro', 'r')
            lines = fr_sol.readlines()
            grofile = [
                lines[j] for j in range(len(lines))
                if (j < (len(lines) - 1 - des * Natom) or j == len(lines) - 1)
            ]
            fw_sol = open('./gromacs/step2/add_' + str(i) + '.gro', 'w')
            grofile[1] = str(int(grofile[1]) - des * Natom) + '\n'
            for k in range(len(grofile)):
                print(grofile[k], end='', file=fw_sol)
            fr_sol.close()
            fw_sol.close()

    return N_sol


def step2_2(Nsol):

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Solvent = system['Solvent_CTA']
    boxsize = system['Boxsize']
    x, y, z = boxsize[0], boxsize[1], boxsize[2]
    Ion = []
    cout = 0
    CPU = 4
    GPU = 0
    GPU_id = ''
    if 'CPU' in system:
        CPU = system['CPU']
    else:
        print('The ntheard uses the default value! (CPU=4)')
    if 'GPU' in system:
        GPU = system['GPU']
    else:
        print('The number of GPU uses the default value! (GPU=0)')
    if GPU != 0:
        for id in range(int(GPU)):
            GPU_id = GPU_id + str(id)

    f_em = open('./gromacs/step2/em.mdp', 'w')
    print('integrator              = steep',
          'emtol                   = 1000.0',
          'nsteps                  = 10000',
          'nstlist                 = 10',
          'cutoff-scheme           = Verlet',
          'rlist                   = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          sep='\n',
          file=f_em)

    f_top = open('./gromacs/step2/system.top', 'w')
    print('#include "../../topper/forcefield.itp"',
          '#include "../../topper/PISA.itp"',
          sep='\n',
          file=f_top)
    if 'Ion_CTA' in system:
        Ion = system['Ion_CTA']
        for i in range(len(Ion[0])):
            print('#include "../../topper/' + str(Ion[0][i]) + '.itp"',
                  file=f_top)
    for i in range(len(Solvent[0])):
        print('#include "../../topper/' + str(Solvent[0][i]) + '.itp"',
              file=f_top)

    print('\n[ system ]',
          'System',
          '\n[ molecules ]',
          'PISA     1',
          sep='\n',
          file=f_top)

    if 'Ion_CTA' in system:
        nIon = [0] * len(Ion[0])
        for i in range(len(Ion[0])):
            nIon = int(0.6023 * x * y * z * Ion[1][i] + 0.5)
            print(str(Ion[0][i]) + '     ' + str(nIon), file=f_top)

    for i in range(len(Solvent[0])):
        print(str(Solvent[0][i]) + '     ' + str(Nsol[i]), file=f_top)

    f_em.close()
    f_top.close()

    flow2 = subprocess.Popen('gmx grompp -f em.mdp -c add_' +
                             str(len(Solvent[0]) - 1) +
                             '.gro -p system.top -o em.tpr -maxwarn -1',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True,
                             cwd='./gromacs/step2')
    flow2.communicate()

    if GPU != 0:
        flow3 = subprocess.Popen(
            'gmx mdrun -v -deffnm em -pin on -ntmpi 1 -ntomp ' + str(CPU) +
            ' -nb gpu',
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True,
            cwd='./gromacs/step2')
        flow3.communicate()
        if flow3.returncode != 0:
            cout += 1
    else:
        flow3 = subprocess.Popen('gmx mdrun -v -deffnm em -ntmpi 1 -ntomp ' +
                                 str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step2')
        flow3.communicate()
        if flow3.returncode != 0:
            cout += 1

    T = 300
    P = 1.0
    Time = 0
    if 'Temperature' in system:
        T = system['Temperature']
    else:
        print('The Temperature uses the default value! (Temperature=300)')
    if 'Pressure' in system:
        P = system['Pressure']
    else:
        print('The Pressure used the default value! (Pressure=1.0)')

    if len(Solvent[0]) == 1:
        Time = 5000000
    if len(Solvent[0]) > 1:
        Time = 10000000

    f_eq = open('./gromacs/step2/eq.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.001',
          'nsteps                  = ' + str(Time),
          'nstlog                  = 0',
          'nstcalcenergy           = 0',
          'nstenergy               = 0',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = berendsen',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = 323',
          'pcoupl                  = berendsen',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          'gen-vel                 = yes',
          'gen-temp                = ' + str(T),
          'gen-seed                = -1',
          sep='\n',
          file=f_eq)
    f_eq.close()

    flow4 = subprocess.Popen(
        'gmx grompp -f eq.mdp -c em.gro -p system.top -o eq.tpr -maxwarn -1',
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd='./gromacs/step2')
    flow4.wait()

    if GPU != 0:
        if GPU > 1:
            flow5 = subprocess.Popen(
                'gmx mdrun -v -deffnm eq -pin on -ntmpi 4 -ntomp ' +
                str(int(CPU / 4)) + ' -pme gpu -npme 1 -nb gpu -gpu_id ' +
                GPU_id,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step2')
            flow5.communicate()
            if flow5.returncode != 0:
                cout += 1
        if GPU == 1:
            flow5 = subprocess.Popen(
                'gmx mdrun -v -deffnm eq -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                ' -pme gpu -nb gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step2')
            flow5.communicate()
            if flow5.returncode != 0:
                cout += 1
    else:
        flow5 = subprocess.Popen('gmx mdrun -v -deffnm eq -ntmpi 1 -ntomp ' +
                                 str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step2')
        flow5.communicate()
        if flow5.returncode != 0:
            cout += 1

    f_md = open('./gromacs/step2/md.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.001',
          'nsteps                  = ' + str(Time),
          'nstlog                  = 0',
          'nstcalcenergy           = 0',
          'nstenergy               = 0',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = Nose-Hoover',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = ' + str(T),
          'pcoupl                  = Parrinello-Rahman',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          sep='\n',
          file=f_md)
    f_md.close()

    flow6 = subprocess.Popen(
        'gmx grompp -f md.mdp -c eq.gro -p system.top -o md.tpr -maxwarn -1',
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd='./gromacs/step2')
    flow6.wait()

    if GPU != 0:
        if GPU > 1:
            flow7 = subprocess.Popen(
                'gmx mdrun -v -deffnm md -pin on -ntmpi 4 -ntomp ' +
                str(int(CPU / 4)) + ' -pme gpu -npme 1 -nb gpu -gpu_id ' +
                GPU_id,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step2')
            flow7.communicate()
            if flow7.returncode != 0:
                cout += 1
        if GPU == 1:
            flow7 = subprocess.Popen(
                'gmx mdrun -v -deffnm md -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                ' -pme gpu -nb gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step2')
            flow7.communicate()
            if flow7.returncode != 0:
                cout += 1
    else:
        flow7 = subprocess.Popen('gmx mdrun -v -deffnm md -ntmpi 1 -ntomp ' +
                                 str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step2')
        flow7.communicate()
        if flow7.returncode != 0:
            cout += 1

    if cout != 0:
        return False
    else:
        return True


def step3_1():

    folder = os.path.exists('./gromacs/step3')
    if not folder:
        os.makedirs('./gromacs/step3')

    flow1 = subprocess.Popen('cp ./system.top ../step3/system.top',
                             shell=True,
                             cwd='./gromacs/step2')
    flow1.wait()
    flow2 = subprocess.Popen('cp ./md.gro ../step3/0.gro',
                             shell=True,
                             cwd='./gromacs/step2')
    flow2.wait()
    flow3 = subprocess.Popen('cp ./PISA.itp ../gromacs/step3/PISA.itp',
                             shell=True,
                             cwd='./topper')
    flow3.wait()

    fin_top = open("./gromacs/step3/system.top", 'r')
    lines = fin_top.readlines()
    fin_top.close()
    fout_top = open("./gromacs/step3/system.top", 'w')
    for i in range(len(lines)):
        if i == 1:
            print('#include "./PISA.itp"', file=fout_top)
        else:
            print(lines[i], end='', file=fout_top)
    fout_top.close()

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    T = 300
    P = 1.0
    Time = 100000
    if 'Temperature' in system:
        T = system['Temperature']
    else:
        print('The Temperature uses the default value! (Temperature=300)')
    if 'Pressure' in system:
        P = system['Pressure']
    else:
        print('The Pressure used the default value! (Pressure=1.0)')
    if 'Polymerization rate' in system:
        TotalPr = system['Polymerization rate']
        if TotalPr != []:
            if 'MA' in TotalPr[0]:
                Time = 1000000 * TotalPr[0]['MA']
            else:
                print(
                    'The Polymerization rate of MA uses the default value! (Polymerization rate=0.1)'
                )
        else:
            print(
                'The Polymerization rate uses the default value! (Polymerization rate=0.1)'
            )
    else:
        print(
            'The Polymerization rate uses the default value! (Polymerization rate=0.1)'
        )
    f_sys.close()

    f_md = open('./gromacs/step3/raft.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.001',
          'nsteps                  = ' + str(int(Time)),
          'nstlog                  = 0',
          'nstcalcenergy           = 0',
          'nstenergy               = 0',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = berendsen',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = ' + str(T),
          'pcoupl                  = berendsen',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          'gen-vel                 = yes',
          'gen-temp                = ' + str(T),
          'gen-seed                = -1',
          'cos_acceleration        = 0.1',
          sep='\n',
          file=f_md)
    f_md.close()

    f_eq = open('./gromacs/step3/eq.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.0001',
          'nsteps                  = 10000',
          'nstlog                  = 0',
          'nstcalcenergy           = 0',
          'nstenergy               = 0',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = berendsen',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = ' + str(T),
          'pcoupl                  = berendsen',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          'gen-vel                 = yes',
          'gen-temp                = ' + str(T),
          'gen-seed                = -1',
          sep='\n',
          file=f_eq)
    f_eq.close()

    f_em = open('./gromacs/step3/em.mdp', 'w')
    print('integrator              = steep',
          'emtol                   = 1000.0',
          'nsteps                  = 1000',
          'nstlist                 = 10',
          'cutoff-scheme           = Verlet',
          'rlist                   = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          sep='\n',
          file=f_em)
    f_em.close()


def step3_2(frame, reac):

    cout = 0
    CPU = 4
    GPU = 0
    GPU_id = ''
    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    if 'CPU' in system:
        CPU = system['CPU']
    else:
        print('The ntheard uses the default value! (CPU=4)')
    if 'GPU' in system:
        GPU = system['GPU']
    else:
        print('The number of GPU uses the default value! (GPU=0)')
    if GPU != 0:
        for id in range(int(GPU)):
            GPU_id = GPU_id + str(id)
    f_sys.close()

    flow = subprocess.Popen('gmx grompp -f em.mdp -c ' + str(frame) +
                            '.gro -p system.top -o em' + str(frame) +
                            '.tpr -maxwarn -1',
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=True,
                            cwd='./gromacs/step3')
    flow.communicate()

    if GPU != 0:
        flow0 = subprocess.Popen('gmx mdrun -v -deffnm em' + str(frame) +
                                 ' -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                                 ' -nb gpu',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step3')
        flow0.communicate()
        if flow0.returncode != 0:
            cout += 1
    else:
        flow0 = subprocess.Popen('gmx mdrun -v -deffnm em' + str(frame) +
                                 ' -ntmpi 1 -ntomp ' + str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step3')
        flow0.communicate()
        if flow0.returncode != 0:
            cout += 1
    if reac != []:
        flow1 = subprocess.Popen('gmx grompp -f eq.mdp -c em' + str(frame) +
                                 '.gro -p system.top -o eq' + str(frame) +
                                 '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step3')
        flow1.wait()

        if GPU != 0:
            if GPU > 1:
                flow2 = subprocess.Popen(
                    'gmx mdrun -v -deffnm eq' + str(frame) +
                    ' -pin on -ntmpi 4 -ntomp ' + str(int(CPU / 4)) +
                    ' -pme gpu -npme 1 -nb gpu -gpu_id ' + GPU_id,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step3')
                flow2.communicate()
                if flow2.returncode != 0:
                    cout += 1
            if GPU == 1:
                flow2 = subprocess.Popen(
                    'gmx mdrun -v -deffnm eq' + str(frame) +
                    ' -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                    ' -pme gpu -nb gpu',
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step3')
                flow2.communicate()
                if flow2.returncode != 0:
                    cout += 1
        else:
            flow2 = subprocess.Popen('gmx mdrun -v -deffnm eq' + str(frame) +
                                     ' -ntmpi 1 -ntomp ' + str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step3')
            flow2.communicate()
            if flow2.returncode != 0:
                cout += 1

        flow3 = subprocess.Popen('gmx grompp -f raft.mdp -c eq' + str(frame) +
                                 '.gro -p system.top -o ' + str(frame + 1) +
                                 '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step3')
        flow3.wait()

        if GPU != 0:
            if GPU > 1:
                flow4 = subprocess.Popen(
                    'gmx mdrun -v -deffnm ' + str(frame + 1) +
                    ' -pin on -ntmpi 4 -ntomp ' + str(int(CPU / 4)) +
                    ' -pme gpu -npme 1 -nb gpu -gpu_id ' + GPU_id,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step3')
                flow4.communicate()
                if flow4.returncode != 0:
                    cout += 1
            if GPU == 1:
                flow4 = subprocess.Popen(
                    'gmx mdrun -v -deffnm ' + str(frame + 1) +
                    ' -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                    ' -pme gpu -nb gpu',
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step3')
                flow4.communicate()
                if flow4.returncode != 0:
                    cout += 1
        else:
            flow4 = subprocess.Popen('gmx mdrun -v -deffnm ' + str(frame + 1) +
                                     ' -ntmpi 1 -ntomp ' + str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step3')
            flow4.communicate()
            if flow4.returncode != 0:
                cout += 1
    else:
        flow3 = subprocess.Popen('gmx grompp -f raft.mdp -c em' + str(frame) +
                                 '.gro -p system.top -o ' + str(frame + 1) +
                                 '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step3')
        flow3.wait()

        if GPU != 0:
            if GPU > 1:
                flow4 = subprocess.Popen(
                    'gmx mdrun -v -deffnm ' + str(frame + 1) +
                    ' -pin on -ntmpi 4 -ntomp ' + str(int(CPU / 4)) +
                    ' -pme gpu -npme 1 -nb gpu -gpu_id ' + GPU_id,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step3')
                flow4.communicate()
                if flow4.returncode != 0:
                    cout += 1
            if GPU == 1:
                flow4 = subprocess.Popen(
                    'gmx mdrun -v -deffnm ' + str(frame + 1) +
                    ' -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                    ' -pme gpu -nb gpu',
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step3')
                flow4.communicate()
                if flow4.returncode != 0:
                    cout += 1
        else:
            flow4 = subprocess.Popen('gmx mdrun -v -deffnm ' + str(frame + 1) +
                                     ' -ntmpi 1 -ntomp ' + str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step3')
            flow4.communicate()
            if flow4.returncode != 0:
                cout += 1
    if cout != 0:
        return False
    else:
        return True


def step4_1():

    cout = 0
    f_sys = open('./control/system.json', )
    system = json.load(f_sys)

    if 'Solvent_PISA' in system:
        Solvent = system['Solvent_PISA']
        if len(Solvent) != 2:
            print(
                'Please add: "Solvent_PISA": [["filename_1", "filename_2"], [volume ratio_1, volume ratio_2]]'
            )
            print('e.g. "Solvent_PISA": [["water", "alcohol"], [0.7, 0.3]]')
            cout += 1
        else:
            for file in Solvent[0]:
                if not os.access("./topper/" + str(file) + '.pdb', os.F_OK):
                    if not os.access("./topper/" + str(file) + '.itp',
                                     os.F_OK):
                        print('The files of solvent (' + file +
                              ') are not existed!')
                        cout += 1
            total = 0.0
            for ratio in Solvent[1]:
                total = total + ratio
            if int(total * 10) / 10 != 1.0:
                print('The sum of the volume ratios is not 1.0!')
                cout += 1
    else:
        print('The Solvent is missing!')
        print(
            'Please add: "Solvent_PISA": [["filename_1", "filename_2"], [volume ratio_1, volume ratio_2]]'
        )
        print('e.g. "Solvent_PISA": [["water", "alcohol"], [0.7, 0.3]]')
        cout += 1

    if 'Ion_PISA' in system:
        Ion = system['Ion_PISA']
        if len(Ion) != 2:
            print(
                'Please add: "Ion_PISA": [["filename_1", "filename_2"], [concentration_1 (Mol/L), concentration_2 (Mol/L)]]'
            )
            print('e.g. "Ion_PISA": [["Na", "Cl"], [0.153, 0.153]]')
            cout += 1
        else:
            for file in Ion[0]:
                if not os.access("./topper/" + str(file) + '.pdb', os.F_OK):
                    if not os.access("./topper/" + str(file) + '.itp',
                                     os.F_OK):
                        print('The files of ion (' + file +
                              ') are not existed!')
                        cout += 1
    f_sys.close()
    if cout == 0:
        return True
    else:
        return False


def step4_2():

    folder = os.path.exists('./gromacs/step4')
    if not folder:
        os.makedirs('./gromacs/step4')

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Solvent = system['Solvent_PISA']
    origin_Solvent = system['Solvent_CTA']
    for sol in Solvent[0][:]:
        if sol in origin_Solvent[0]:
            Solvent[0].remove(sol)
            flow = subprocess.Popen('cp ./gromacs/step1/md_' + str(sol) +
                                    '.gro ./gromacs/step4',
                                    shell=True)
            flow.wait()
    size = [30] * len(Solvent[0])
    for i in range(len(Solvent[0])):
        flow1 = subprocess.Popen('cp ./topper/' + str(Solvent[0][i]) +
                                 '.pdb ./gromacs/step4',
                                 shell=True)
        flow1.wait()
        f_packmol = open('./gromacs/step4/' + str(Solvent[0][i]) + '.inp', 'w')
        print('tolerance 2.4\n',
              'filetype pdb',
              'output  Sys_' + str(Solvent[0][i]) + '.pdb',
              'add_box_sides 1.2',
              '\nstructure ' + str(Solvent[0][i]) + '.pdb',
              '    number 500',
              '    inside box 0 0 0' + ' ' + str(size[i]) + ' ' +
              str(size[i]) + ' ' + str(size[i]),
              'end structure',
              sep='\n',
              file=f_packmol)
        f_packmol.close()

    cout = 0
    for i in range(len(Solvent[0])):
        pack = False
        while not pack:
            flow1 = subprocess.Popen('packmol < ' + str(Solvent[0][i]) +
                                     '.inp',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step4')
            try:
                (msg, errs) = flow1.communicate(timeout=120)
            except subprocess.TimeoutExpired:
                flow1.kill()
                flow1.terminate()
                f_packmol = open(
                    './gromacs/step4/' + str(Solvent[0][i]) + '.inp', 'w')
                size[i] = size[i] * 1.5
                print('tolerance 2.4\n',
                      'filetype pdb',
                      'output  Sys_' + str(Solvent[0][i]) + '.pdb',
                      'add_box_sides 1.2',
                      '\nstructure ' + str(Solvent[0][i]) + '.pdb',
                      '    number 500',
                      '    inside box 0 0 0' + ' ' + str(size[i]) + ' ' +
                      str(size[i]) + ' ' + str(size[i]),
                      'end structure',
                      sep='\n',
                      file=f_packmol)
                f_packmol.close()
            if flow1.returncode == 0:
                pack = True
                cout = 0
                subprocess.Popen('rm -f ' + str(Solvent[0][i]) + '.inp',
                                 shell=True,
                                 cwd='./gromacs/step4')
                subprocess.Popen('rm -f ' + str(Solvent[0][i]) + '.pdb',
                                 shell=True,
                                 cwd='./gromacs/step4')
            else:
                cout += 1
            if cout > 10:
                break

    f_sys.close()
    if cout > 10:
        return False
    else:
        return size


def step4_3(size):

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Solvent = system['Solvent_PISA']
    origin_Solvent = system['Solvent_CTA']
    for sol in Solvent[0][:]:
        if sol in origin_Solvent[0]:
            Solvent[0].remove(sol)
    cout = 0
    CPU = 4
    GPU = 0
    GPU_id = ''
    if 'CPU' in system:
        CPU = system['CPU']
    else:
        print('The ntheard uses the default value! (CPU=4)')
    if 'GPU' in system:
        GPU = system['GPU']
    else:
        print('The number of GPU uses the default value! (GPU=0)')
    if GPU != 0:
        for id in range(int(GPU)):
            GPU_id = GPU_id + str(id)

    for i in range(len(Solvent[0])):
        f_em = open('./gromacs/step4/em' + str(i) + '.mdp', 'w')
        print('integrator              = steep',
              'emtol                   = 1000.0',
              'nsteps                  = 10000',
              'nstlist                 = 10',
              'cutoff-scheme           = Verlet',
              'rlist                   = 1.2',
              'vdwtype                 = Cut-off',
              'vdw-modifier            = Force-switch',
              'rvdw_switch             = 1.0',
              'rvdw                    = 1.2',
              'coulombtype             = pme',
              'rcoulomb                = 1.2',
              'constraints             = h-bonds',
              'constraint_algorithm    = LINCS',
              sep='\n',
              file=f_em)

        f_top = open('./gromacs/step4/' + str(Solvent[0][i]) + '.top', 'w')
        print('#include "../../topper/forcefield.itp"',
              '#include "../../topper/' + str(Solvent[0][i]) + '.itp',
              '\n[ system ]',
              str(Solvent[0][i]),
              '\n[ molecules ]',
              str(Solvent[0][i]) + '     500',
              sep='\n',
              file=f_top)

        f_em.close()
        f_top.close()

        box = str('{:.5f}'.format(size[i] / 10))
        flow1 = subprocess.Popen('gmx editconf -f Sys_' + str(Solvent[0][i]) +
                                 '.pdb -o ' + str(Solvent[0][i]) +
                                 '.gro -box' + ' ' + box + ' ' + box + ' ' +
                                 box,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step4')
        flow1.wait()
        flow2 = subprocess.Popen('gmx grompp -f em' + str(i) + '.mdp -c ' +
                                 str(Solvent[0][i]) + '.gro -p ' +
                                 str(Solvent[0][i]) + '.top -o em_' +
                                 str(Solvent[0][i]) + '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step4')
        flow2.wait()
        if GPU == 0:
            flow3 = subprocess.Popen('gmx mdrun -v -deffnm em_' +
                                     str(Solvent[0][i]) + ' -ntmpi 1 -ntomp ' +
                                     str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step4')
            flow3.communicate()
            if flow3.returncode != 0:
                cout += 1
        else:
            flow3 = subprocess.Popen(
                'gmx mdrun -v -deffnm em_' + str(Solvent[0][i]) +
                ' -pin on -ntmpi 1 -ntomp ' + str(CPU) + ' -nb gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step4')
            flow3.communicate()
            if flow3.returncode != 0:
                cout += 1

    T = 300
    P = 1.0
    if 'Temperature' in system:
        T = system['Temperature']
    else:
        print('The Temperature uses the default value! (Temperature=300)')
    if 'Pressure' in system:
        P = system['Pressure']
    else:
        print('The Pressure used the default value! (Pressure=1.0)')

    for i in range(len(Solvent[0])):

        f_eq = open('./gromacs/step4/eq' + str(i) + '.mdp', 'w')
        print('integrator              = md',
              'dt                      = 0.001',
              'nsteps                  = 1000000',
              'nstlog                  = 0',
              'nstcalcenergy           = 0',
              'nstenergy               = 0',
              'cutoff-scheme           = Verlet',
              'nstlist                 = 20',
              'rlist                   = 1.2',
              'coulombtype             = pme',
              'rcoulomb                = 1.2',
              'vdwtype                 = Cut-off',
              'vdw-modifier            = Force-switch',
              'rvdw_switch             = 1.0',
              'rvdw                    = 1.2',
              'tcoupl                  = berendsen',
              'tc_grps                 = SYSTEM',
              'tau_t                   = 1.0',
              'ref_t                   = ' + str(T),
              'pcoupl                  = berendsen',
              'pcoupltype              = isotropic',
              'tau_p                   = 5.0',
              'compressibility         = 4.5e-5',
              'ref_p                   = ' + str(P),
              'constraints             = h-bonds',
              'constraint_algorithm    = LINCS',
              'continuation            = no',
              'nstcomm                 = 100',
              'comm_mode               = Linear',
              'comm_grps               = SYSTEM',
              'gen-vel                 = yes',
              'gen-temp                = ' + str(T),
              'gen-seed                = -1',
              sep='\n',
              file=f_eq)
        f_eq.close()

        flow4 = subprocess.Popen('gmx grompp -f eq' + str(i) + '.mdp -c em_' +
                                 str(Solvent[0][i]) + '.gro -p ' +
                                 str(Solvent[0][i]) + '.top -o eq_' +
                                 str(Solvent[0][i]) + '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step4')
        flow4.wait()
        if GPU == 0:
            flow5 = subprocess.Popen('gmx mdrun -v -deffnm eq_' +
                                     str(Solvent[0][i]) + ' -ntmpi 1 -ntomp ' +
                                     str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step4')
            flow5.communicate()
            if flow5.returncode != 0:
                cout += 1
        else:
            flow5 = subprocess.Popen(
                'gmx mdrun -v -deffnm eq_' + str(Solvent[0][i]) +
                ' -pin on -ntmpi 1 -ntomp ' + str(CPU) + ' -nb gpu -pme gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step4')
            flow5.communicate()
            if flow5.returncode != 0:
                cout += 1

    for i in range(len(Solvent[0])):
        f_md = open('./gromacs/step4/md' + str(i) + '.mdp', 'w')
        print('integrator              = md',
              'dt                      = 0.001',
              'nsteps                  = 4000000',
              'nstlog                  = 0',
              'nstcalcenergy           = 0',
              'nstenergy               = 0',
              'cutoff-scheme           = Verlet',
              'nstlist                 = 20',
              'rlist                   = 1.2',
              'coulombtype             = pme',
              'rcoulomb                = 1.2',
              'vdwtype                 = Cut-off',
              'vdw-modifier            = Force-switch',
              'rvdw_switch             = 1.0',
              'rvdw                    = 1.2',
              'tcoupl                  = Nose-Hoover',
              'tc_grps                 = SYSTEM',
              'tau_t                   = 1.0',
              'ref_t                   = ' + str(T),
              'pcoupl                  = Parrinello-Rahman',
              'pcoupltype              = isotropic',
              'tau_p                   = 5.0',
              'compressibility         = 4.5e-5',
              'ref_p                   = ' + str(P),
              'constraints             = h-bonds',
              'constraint_algorithm    = LINCS',
              'continuation            = no',
              'nstcomm                 = 100',
              'comm_mode               = Linear',
              'comm_grps               = SYSTEM',
              sep='\n',
              file=f_md)
        f_md.close()

        flow6 = subprocess.Popen('gmx grompp -f md' + str(i) + '.mdp -c eq_' +
                                 str(Solvent[0][i]) + '.gro -p ' +
                                 str(Solvent[0][i]) + '.top -o md_' +
                                 str(Solvent[0][i]) + '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step4')
        flow6.wait()
        if GPU == 0:
            flow7 = subprocess.Popen('gmx mdrun -v -deffnm md_' +
                                     str(Solvent[0][i]) + ' -ntmpi 1 -ntomp ' +
                                     str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step4')
            flow7.communicate()
            if flow7.returncode != 0:
                cout += 1
        else:
            flow7 = subprocess.Popen(
                'gmx mdrun -v -deffnm md_' + str(Solvent[0][i]) +
                ' -pin on -ntmpi 1 -ntomp ' + str(CPU) + ' -nb gpu -pme gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step4')
            flow7.communicate()
            if flow7.returncode != 0:
                cout += 1

    if cout != 0:
        return False
    else:
        return True


def step5_1():

    folder = os.path.exists('./gromacs/step5')
    if not folder:
        os.makedirs('./gromacs/step5')

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Solvent = system['Solvent_PISA']
    fin_pdb = open('./topper/macroCTA.pdb', 'r')
    line = fin_pdb.readline()
    line = fin_pdb.readline()
    line = fin_pdb.readline()
    lines = line.split()
    boxsize = [0] * 3
    boxsize[0], boxsize[1], boxsize[2] = float(lines[1]) / 10, float(
        lines[2]) / 10, float(lines[3]) / 10
    polymer = system['PISA']
    Ion = []
    fin_pdb.close()

    for i in range(len(Solvent[0])):
        flow1 = subprocess.Popen('cp ./gromacs/step4/md_' +
                                 str(Solvent[0][i]) + '.gro ./gromacs/step5',
                                 shell=True)
        flow1.wait()

    flow2 = subprocess.Popen('cp ./topper/macroCTA.pdb ./gromacs/step5',
                             shell=True)
    flow2.wait()
    flow3 = subprocess.Popen('cp ./topper/I.pdb ./gromacs/step5', shell=True)
    flow3.wait()
    flow4 = subprocess.Popen('cp ./topper/MB.pdb ./gromacs/step5', shell=True)
    flow4.wait()

    x, y, z = boxsize[0], boxsize[1], boxsize[2]
    mCTA, mINI, mMB = polymer[0]['CTA'], polymer[0]['INI'], polymer[0]['MB']

    nINI = int(
        round(0.6023 * system['Boxsize'][0] * system['Boxsize'][1] *
              system['Boxsize'][2] * mINI + 0.5))
    if nINI == 0:
        nINI = 1
    nCTA = int(round(nINI * (mCTA / mINI)))
    if nCTA == 0:
        nCTA = 1
    nMB = int(round(nCTA * (mMB / mCTA)))
    if nMB == 0:
        nMB = 1

    f_packmol = open('./gromacs/step5/system.inp', 'w')
    print('tolerance 2.0\n',
          'filetype pdb',
          'output  system.pdb',
          '\nstructure macroCTA.pdb',
          '    number 1',
          '    fixed 0. 0. 0. 0. 0. 0.',
          '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' + str(y * 10) +
          ' ' + str(z * 10),
          'end structure',
          '\nstructure I.pdb',
          '    number ' + str(nINI),
          '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' + str(y * 10) +
          ' ' + str(z * 10),
          'end structure',
          '\nstructure MB.pdb',
          '    number ' + str(nMB),
          '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' + str(y * 10) +
          ' ' + str(z * 10),
          'end structure',
          sep='\n',
          file=f_packmol)

    if 'Ion_PISA' in system:
        Ion = system['Ion_PISA']
        for i in range(len(Ion[0])):
            flow2 = subprocess.Popen('cp ./topper/' + str(Ion[0][i]) +
                                     '.pdb ./gromacs/step5',
                                     shell=True)
            flow2.wait()
            nIon = int(0.6023 * system['Boxsize'][0] * system['Boxsize'][1] *
                       system['Boxsize'][2] * Ion[1][i] + 0.5)
            print('\nstructure ' + str(Ion[0][i]) + '.pdb',
                  '    number ' + str(nIon),
                  '    inside box 0 0 0' + ' ' + str(x * 10) + ' ' +
                  str(y * 10) + ' ' + str(z * 10),
                  'end structure',
                  sep='\n',
                  file=f_packmol)

    f_packmol.close()

    flow6 = subprocess.Popen('packmol < system.inp',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True,
                             cwd='./gromacs/step5')
    try:
        (msg, errs) = flow6.communicate(timeout=180)
    except subprocess.TimeoutExpired:
        flow6.kill()
        flow6.terminate()
    if flow6.returncode != 0:
        print('Packmol is failed, please check the .inp file!')
        return False
    else:
        subprocess.Popen('rm -f system.inp', shell=True, cwd='./gromacs/step5')

    flow7 = subprocess.Popen('gmx editconf -f system.pdb -o system.gro -box' +
                             ' ' + str(x) + ' ' + str(y) + ' ' + str(z),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True,
                             cwd='./gromacs/step5')
    flow7.wait()

    N_sol = [1] * len(Solvent[0])
    for i in range(len(Solvent[0])):
        Nms = 0
        if i == 0:
            flow8 = subprocess.Popen(
                'gmx solvate -cp system.gro -cs md_' + str(Solvent[0][i]) +
                '.gro -o add_' + str(i) + '.gro -box' + ' ' + str(x) + ' ' +
                str(y) + ' ' + str(z),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step5')
            stdout, _ = flow8.communicate()
            out = stdout.decode('utf-8')
            Num = re.search(r'Number of solvent molecules:.*[0-9]', out)
            Nms = int(Num.group().split(':')[1])
            if flow8.returncode != 0:
                print('Add solvent' + str(Solvent[0][i]) + ' is failed!')
                return False
        if i > 0:
            flow8 = subprocess.Popen(
                'gmx solvate -cp add_' + str(i - 1) + '.gro -cs md_' +
                str(Solvent[0][i]) + '.gro -o add_' + str(i) + '.gro -box' +
                ' ' + str(x) + ' ' + str(y) + ' ' + str(z),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step5')
            stdout, _ = flow8.communicate()
            out = stdout.decode('utf-8')
            Num = re.search(r'Number of solvent molecules:.*[0-9]', out)
            Nms = int(Num.group().split(':')[1])
            if flow8.returncode != 0:
                print('Add solvent' + str(Solvent[0][i]) + ' is failed!')
                return False
            if i == len(Solvent[0]) - 1:
                N_sol[i] = Nms
        des = 0
        f_gro = open('./gromacs/step4/md_' + str(Solvent[0][i]) + '.gro', 'r')
        line = f_gro.readline()
        line = f_gro.readline()
        Natom = int(int(line) / 500)
        f_gro.close()
        if i == 0:
            des = int((1.0 - Solvent[1][i]) * Nms)
            if des == 0 and len(N_sol) > 1:
                des = 1
            N_sol[i] = Nms - des
            fr_sol = open('./gromacs/step5/add_' + str(i) + '.gro', 'r')
            lines = fr_sol.readlines()
            grofile = [
                lines[j] for j in range(len(lines))
                if (j < (len(lines) - 1 - des * Natom) or j == len(lines) - 1)
            ]
            grofile[1] = str(int(grofile[1]) - des * Natom) + '\n'
            fw_sol = open('./gromacs/step5/add_' + str(i) + '.gro', 'w')
            for k in range(len(grofile)):
                print(grofile[k], end='', file=fw_sol)
            fr_sol.close()
            fw_sol.close()
        if i > 0 and i != len(Solvent[0]) - 1:
            pde = 0
            for j in range(i + 1):
                pde = pde + Solvent[1][j]
            des = int((1 - pde) / Solvent[1][i] * Nms)
            if des == 0:
                des = 1
            N_sol[i] = Nms - des
            fr_sol = open('./gromacs/step5/add_' + str(i) + '.gro', 'r')
            lines = fr_sol.readlines()
            grofile = [
                lines[j] for j in range(len(lines))
                if (j < (len(lines) - 1 - des * Natom) or j == len(lines) - 1)
            ]
            fw_sol = open('./gromacs/step5/add_' + str(i) + '.gro', 'w')
            grofile[1] = str(int(grofile[1]) - des * Natom) + '\n'
            for k in range(len(grofile)):
                print(grofile[k], end='', file=fw_sol)
            fr_sol.close()
            fw_sol.close()

    return N_sol


def step5_2(Nsol):

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    Solvent = system['Solvent_PISA']
    boxsize = system['Boxsize']
    x, y, z = boxsize[0], boxsize[1], boxsize[2]
    Ion = []
    cout = 0
    CPU = 4
    GPU = 0
    GPU_id = ''
    if 'CPU' in system:
        CPU = system['CPU']
    else:
        print('The ntheard uses the default value! (CPU=4)')
    if 'GPU' in system:
        GPU = system['GPU']
    else:
        print('The number of GPU uses the default value! (GPU=0)')
    if GPU != 0:
        for id in range(int(GPU)):
            GPU_id = GPU_id + str(id)

    f_em = open('./gromacs/step5/em.mdp', 'w')
    print('integrator              = steep',
          'emtol                   = 1000.0',
          'nsteps                  = 10000',
          'nstlist                 = 10',
          'cutoff-scheme           = Verlet',
          'rlist                   = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          sep='\n',
          file=f_em)

    f_top = open('./gromacs/step5/system.top', 'w')
    print('#include "../../topper/forcefield.itp"',
          '#include "../../topper/PISA.itp"',
          sep='\n',
          file=f_top)
    if 'Ion_PISA' in system:
        Ion = system['Ion_PISA']
        for i in range(len(Ion[0])):
            print('#include "../../topper/' + str(Ion[0][i]) + '.itp',
                  file=f_top)
    for i in range(len(Solvent[0])):
        print('#include "../../topper/' + str(Solvent[0][i]) + '.itp',
              file=f_top)

    print('\n[ system ]',
          'System',
          '\n[ molecules ]',
          'PISA     1',
          sep='\n',
          file=f_top)

    if 'Ion_PISA' in system:
        nIon = [0] * len(Ion[0])
        for i in range(len(Ion[0])):
            nIon = int(0.6023 * x * y * z * Ion[1][i] + 0.5)
            print(str(Ion[0][i]) + '     ' + str(nIon), file=f_top)

    for i in range(len(Solvent[0])):
        print(str(Solvent[0][i]) + '     ' + str(Nsol[i]), file=f_top)

    f_em.close()
    f_top.close()

    flow2 = subprocess.Popen('gmx grompp -f em.mdp -c add_' +
                             str(len(Solvent[0]) - 1) +
                             '.gro -p system.top -o em.tpr -maxwarn -1',
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True,
                             cwd='./gromacs/step5')
    flow2.communicate()

    if GPU != 0:
        flow3 = subprocess.Popen(
            'gmx mdrun -v -deffnm em -pin on -ntmpi 1 -ntomp ' + str(CPU) +
            ' -nb gpu',
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True,
            cwd='./gromacs/step5')
        flow3.communicate()
        if flow3.returncode != 0:
            cout += 1
    else:
        flow3 = subprocess.Popen('gmx mdrun -v -deffnm em -ntmpi 1 -ntomp ' +
                                 str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step5')
        flow3.communicate()
        if flow3.returncode != 0:
            cout += 1

    T = 300
    P = 1.0
    Time = 0
    if 'Temperature' in system:
        T = system['Temperature']
    else:
        print('The Temperature uses the default value! (Temperature=300)')
    if 'Pressure' in system:
        P = system['Pressure']
    else:
        print('The Pressure used the default value! (Pressure=1.0)')

    if len(Solvent[0]) == 1:
        Time = 5000000
    if len(Solvent[0]) > 1:
        Time = 10000000

    f_eq = open('./gromacs/step5/eq.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.001',
          'nsteps                  = ' + str(Time),
          'nstlog                  = 0',
          'nstcalcenergy           = 0',
          'nstenergy               = 0',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = berendsen',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = 323',
          'pcoupl                  = berendsen',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          'gen-vel                 = yes',
          'gen-temp                = 323',
          'gen-seed                = -1',
          sep='\n',
          file=f_eq)
    f_eq.close()

    flow4 = subprocess.Popen(
        'gmx grompp -f eq.mdp -c em.gro -p system.top -o eq.tpr -maxwarn -1',
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd='./gromacs/step5')
    flow4.wait()

    if GPU != 0:
        if GPU > 1:
            flow5 = subprocess.Popen(
                'gmx mdrun -v -deffnm eq -pin on -ntmpi 4 -ntomp ' +
                str(int(CPU / 4)) + ' -pme gpu -npme 1 -nb gpu -gpu_id ' +
                GPU_id,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step5')
            flow5.communicate()
            if flow5.returncode != 0:
                cout += 1
        if GPU == 1:
            flow5 = subprocess.Popen(
                'gmx mdrun -v -deffnm eq -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                ' -pme gpu -nb gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step5')
            flow5.communicate()
            if flow5.returncode != 0:
                cout += 1
    else:
        flow5 = subprocess.Popen('gmx mdrun -v -deffnm eq -ntmpi 1 -ntomp ' +
                                 str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step5')
        flow5.communicate()
        if flow5.returncode != 0:
            cout += 1

    f_md = open('./gromacs/step5/md.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.001',
          'nsteps                  = ' + str(Time),
          'nstlog                  = 0',
          'nstcalcenergy           = 0',
          'nstenergy               = 0',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = Nose-Hoover',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = ' + str(T),
          'pcoupl                  = Parrinello-Rahman',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          sep='\n',
          file=f_md)
    f_md.close()

    flow6 = subprocess.Popen(
        'gmx grompp -f md.mdp -c eq.gro -p system.top -o md.tpr -maxwarn -1',
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd='./gromacs/step5')
    flow6.wait()

    if GPU != 0:
        if GPU > 1:
            flow7 = subprocess.Popen(
                'gmx mdrun -v -deffnm md -pin on -ntmpi 4 -ntomp ' +
                str(int(CPU / 4)) + ' -pme gpu -npme 1 -nb gpu -gpu_id ' +
                GPU_id,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step5')
            flow7.communicate()
            if flow7.returncode != 0:
                cout += 1
        if GPU == 1:
            flow7 = subprocess.Popen(
                'gmx mdrun -v -deffnm md -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                ' -pme gpu -nb gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step5')
            flow7.communicate()
            if flow7.returncode != 0:
                cout += 1
    else:
        flow7 = subprocess.Popen('gmx mdrun -v -deffnm md -ntmpi 1 -ntomp ' +
                                 str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step5')
        flow7.communicate()
        if flow7.returncode != 0:
            cout += 1

    if cout != 0:
        return False
    else:
        return True


def step6_1():

    folder = os.path.exists('./gromacs/step6')
    if not folder:
        os.makedirs('./gromacs/step6')

    flow1 = subprocess.Popen('cp ./system.top ../step6/system.top',
                             shell=True,
                             cwd='./gromacs/step5')
    flow1.wait()
    flow2 = subprocess.Popen('cp ./md.gro ../step6/0.gro',
                             shell=True,
                             cwd='./gromacs/step5')
    flow2.wait()
    flow3 = subprocess.Popen('cp ./PISA.itp ../gromacs/step6/PISA.itp',
                             shell=True,
                             cwd='./topper')
    flow3.wait()

    fin_top = open("./gromacs/step6/system.top", 'r')
    lines = fin_top.readlines()
    fin_top.close()
    fout_top = open("./gromacs/step6/system.top", 'w')
    for i in range(len(lines)):
        if i == 1:
            print('#include "./PISA.itp"', file=fout_top)
        else:
            print(lines[i], end='', file=fout_top)
    fout_top.close()

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    T = 300
    P = 1.0
    Time = 100000
    if 'Temperature' in system:
        T = system['Temperature']
    else:
        print('The Temperature uses the default value! (Temperature=300)')
    if 'Pressure' in system:
        P = system['Pressure']
    else:
        print('The Pressure used the default value! (Pressure=1.0)')
    if 'Polymerization rate' in system:
        TotalPr = system['Polymerization rate']
        if TotalPr != []:
            if 'MB' in TotalPr[0]:
                Time = 1000000 * TotalPr[0]['MB']
            else:
                print(
                    'The Polymerization rate of MB uses the default value! (Polymerization rate=0.1)'
                )
        else:
            print(
                'The Polymerization rate uses the default value! (Polymerization rate=0.1)'
            )
    else:
        print(
            'The Polymerization rate uses the default value! (Polymerization rate=0.1)'
        )
    f_sys.close()

    f_md = open('./gromacs/step6/raft.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.001',
          'nsteps                  = ' + str(int(Time)),
          'nstlog                  = 0',
          'nstcalcenergy           = 0',
          'nstenergy               = 0',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = berendsen',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = ' + str(T),
          'pcoupl                  = berendsen',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          'gen-vel                 = yes',
          'gen-temp                = ' + str(T),
          'gen-seed                = -1',
          sep='\n',
          file=f_md)
    f_md.close()

    f_eq = open('./gromacs/step6/eq.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.0001',
          'nsteps                  = 10000',
          'nstlog                  = 0',
          'nstcalcenergy           = 0',
          'nstenergy               = 0',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = berendsen',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = ' + str(T),
          'pcoupl                  = berendsen',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          'gen-vel                 = yes',
          'gen-temp                = ' + str(T),
          'gen-seed                = -1',
          sep='\n',
          file=f_eq)
    f_eq.close()

    f_em = open('./gromacs/step6/em.mdp', 'w')
    print('integrator              = steep',
          'emtol                   = 1000.0',
          'nsteps                  = 1000',
          'nstlist                 = 10',
          'cutoff-scheme           = Verlet',
          'rlist                   = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          sep='\n',
          file=f_em)
    f_em.close()


def step6_2(frame, reac):

    cout = 0
    CPU = 4
    GPU = 0
    GPU_id = ''
    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    if 'CPU' in system:
        CPU = system['CPU']
    else:
        print('The ntheard uses the default value! (CPU=4)')
    if 'GPU' in system:
        GPU = system['GPU']
    else:
        print('The number of GPU uses the default value! (GPU=0)')
    if GPU != 0:
        for id in range(int(GPU)):
            GPU_id = GPU_id + str(id)
    f_sys.close()

    flow = subprocess.Popen('gmx grompp -f em.mdp -c ' + str(frame) +
                            '.gro -p system.top -o em' + str(frame) +
                            '.tpr -maxwarn -1',
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=True,
                            cwd='./gromacs/step6')
    flow.communicate()

    if GPU != 0:
        flow0 = subprocess.Popen('gmx mdrun -v -deffnm em' + str(frame) +
                                 ' -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                                 ' -nb gpu',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step6')
        flow0.communicate()
        if flow0.returncode != 0:
            cout += 1
    else:
        flow0 = subprocess.Popen('gmx mdrun -v -deffnm em' + str(frame) +
                                 ' -ntmpi 1 -ntomp ' + str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step6')
        flow0.communicate()
        if flow0.returncode != 0:
            cout += 1
    if reac != []:
        flow1 = subprocess.Popen('gmx grompp -f eq.mdp -c em' + str(frame) +
                                 '.gro -p system.top -o eq' + str(frame) +
                                 '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step6')
        flow1.wait()

        if GPU != 0:
            if GPU > 1:
                flow2 = subprocess.Popen(
                    'gmx mdrun -v -deffnm eq' + str(frame) +
                    ' -pin on -ntmpi 4 -ntomp ' + str(int(CPU / 4)) +
                    ' -pme gpu -npme 1 -nb gpu -gpu_id ' + GPU_id,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step6')
                flow2.communicate()
                if flow2.returncode != 0:
                    cout += 1
            if GPU == 1:
                flow2 = subprocess.Popen(
                    'gmx mdrun -v -deffnm eq' + str(frame) +
                    ' -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                    ' -pme gpu -nb gpu',
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step6')
                flow2.communicate()
                if flow2.returncode != 0:
                    cout += 1
        else:
            flow2 = subprocess.Popen('gmx mdrun -v -deffnm eq' + str(frame) +
                                     ' -ntmpi 1 -ntomp ' + str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step6')
            flow2.communicate()
            if flow2.returncode != 0:
                cout += 1

        flow3 = subprocess.Popen('gmx grompp -f raft.mdp -c eq' + str(frame) +
                                 '.gro -p system.top -o ' + str(frame + 1) +
                                 '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step6')
        flow3.wait()

        if GPU != 0:
            if GPU > 1:
                flow4 = subprocess.Popen(
                    'gmx mdrun -v -deffnm ' + str(frame + 1) +
                    ' -pin on -ntmpi 4 -ntomp ' + str(int(CPU / 4)) +
                    ' -pme gpu -npme 1 -nb gpu -gpu_id ' + GPU_id,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step6')
                flow4.communicate()
                if flow4.returncode != 0:
                    cout += 1
            if GPU == 1:
                flow4 = subprocess.Popen(
                    'gmx mdrun -v -deffnm ' + str(frame + 1) +
                    ' -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                    ' -pme gpu -nb gpu',
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step6')
                flow4.communicate()
                if flow4.returncode != 0:
                    cout += 1
        else:
            flow4 = subprocess.Popen('gmx mdrun -v -deffnm ' + str(frame + 1) +
                                     ' -ntmpi 1 -ntomp ' + str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step6')
            flow4.communicate()
            if flow4.returncode != 0:
                cout += 1
    else:
        flow3 = subprocess.Popen('gmx grompp -f raft.mdp -c em' + str(frame) +
                                 '.gro -p system.top -o ' + str(frame + 1) +
                                 '.tpr -maxwarn -1',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step6')
        flow3.wait()

        if GPU != 0:
            if GPU > 1:
                flow4 = subprocess.Popen(
                    'gmx mdrun -v -deffnm ' + str(frame + 1) +
                    ' -pin on -ntmpi 4 -ntomp ' + str(int(CPU / 4)) +
                    ' -pme gpu -npme 1 -nb gpu -gpu_id ' + GPU_id,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step6')
                flow4.communicate()
                if flow4.returncode != 0:
                    cout += 1
            if GPU == 1:
                flow4 = subprocess.Popen(
                    'gmx mdrun -v -deffnm ' + str(frame + 1) +
                    ' -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                    ' -pme gpu -nb gpu',
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    cwd='./gromacs/step6')
                flow4.communicate()
                if flow4.returncode != 0:
                    cout += 1
        else:
            flow4 = subprocess.Popen('gmx mdrun -v -deffnm ' + str(frame + 1) +
                                     ' -ntmpi 1 -ntomp ' + str(CPU) +
                                     ' -pin on -pinstride 1 -pinoffset 0',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True,
                                     cwd='./gromacs/step6')
            flow4.communicate()
            if flow4.returncode != 0:
                cout += 1
    if cout != 0:
        return False
    else:
        return True


def step7_1(filename):

    folder = os.path.exists('./gromacs/step7')
    if not folder:
        os.makedirs('./gromacs/step7')

    flow1 = subprocess.Popen('cp ./' + filename + '.gro ../step7/pisa.gro',
                             shell=True,
                             cwd='./gromacs/step6')
    flow1.wait()
    flow1 = subprocess.Popen('cp ./PISA.itp ../step7/PISA.itp',
                             shell=True,
                             cwd='./gromacs/step6')
    flow1.wait()
    flow1 = subprocess.Popen('cp ./system.top ../step7/system.top',
                             shell=True,
                             cwd='./gromacs/step6')
    flow1.wait()

    f_sys = open('./control/system.json', )
    system = json.load(f_sys)
    cout = 0
    CPU = 4
    GPU = 0
    GPU_id = ''
    if 'CPU' in system:
        CPU = system['CPU']
    else:
        print('The ntheard uses the default value! (CPU=4)')
    if 'GPU' in system:
        GPU = system['GPU']
    else:
        print('The number of GPU uses the default value! (GPU=0)')
    if GPU != 0:
        for id in range(int(GPU)):
            GPU_id = GPU_id + str(id)
    T = 300
    P = 1.0
    Balance_Time = 10000000
    if 'Temperature' in system:
        T = system['Temperature']
    else:
        print('The Temperature uses the default value! (Temperature=300)')
    if 'Pressure' in system:
        P = system['Pressure']
    else:
        print('The Pressure used the default value! (Pressure=1.0)')
    if 'Balance time' in system:
        Balance_Time = system['Balance time'] * 1000000
    else:
        print('The Balance time uses the default value! (Balance time=10)')
    f_sys.close()

    f_em = open('./gromacs/step7/em.mdp', 'w')
    print('integrator              = steep',
          'emtol                   = 1000.0',
          'nsteps                  = 10000',
          'nstlist                 = 10',
          'cutoff-scheme           = Verlet',
          'rlist                   = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          sep='\n',
          file=f_em)
    f_em.close()

    f_md = open('./gromacs/step7/md.mdp', 'w')
    print('integrator              = md',
          'dt                      = 0.001',
          'nsteps                  = ' + str(Balance_Time),
          'nstlog                  = 100000',
          'nstxout                 = 100000',
          'nstvout                 = 100000',
          'nstfout                 = 100000',
          'nstxout-compressed      = 100000',
          'nstcalcenergy           = 10000',
          'nstenergy               = 10000',
          'cutoff-scheme           = Verlet',
          'nstlist                 = 20',
          'rlist                   = 1.2',
          'coulombtype             = pme',
          'rcoulomb                = 1.2',
          'vdwtype                 = Cut-off',
          'vdw-modifier            = Force-switch',
          'rvdw_switch             = 1.0',
          'rvdw                    = 1.2',
          'tcoupl                  = berendsen',
          'tc_grps                 = SYSTEM',
          'tau_t                   = 1.0',
          'ref_t                   = ' + str(T),
          'pcoupl                  = berendsen',
          'pcoupltype              = isotropic',
          'tau_p                   = 5.0',
          'compressibility         = 4.5e-5',
          'ref_p                   = ' + str(P),
          'constraints             = h-bonds',
          'constraint_algorithm    = LINCS',
          'continuation            = no',
          'nstcomm                 = 100',
          'comm_mode               = Linear',
          'comm_grps               = SYSTEM',
          sep='\n',
          file=f_md)
    f_md.close()

    flow = subprocess.Popen(
        'gmx grompp -f em.mdp -c pisa.gro -p system.top -o em.tpr -maxwarn -1',
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd='./gromacs/step7')
    flow.communicate()

    if GPU != 0:
        flow0 = subprocess.Popen(
            'gmx mdrun -v -deffnm em -pin on -ntmpi 1 -ntomp ' + str(CPU) +
            ' -nb gpu',
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True,
            cwd='./gromacs/step7')
        flow0.communicate()
        if flow0.returncode != 0:
            cout += 1
    else:
        flow0 = subprocess.Popen('gmx mdrun -v -deffnm em -ntmpi 1 -ntomp ' +
                                 str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step7')
        flow0.communicate()
        if flow0.returncode != 0:
            cout += 1

    flow2 = subprocess.Popen(
        'gmx grompp -f md.mdp -c em.gro -p system.top -o md.tpr -maxwarn -1',
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd='./gromacs/step7')
    flow2.wait()

    if GPU != 0:
        if GPU > 1:
            flow3 = subprocess.Popen(
                'gmx mdrun -v -deffnm md -pin on -ntmpi 4 -ntomp ' +
                str(int(CPU / 4)) + ' -pme gpu -npme 1 -nb gpu -gpu_id ' +
                GPU_id,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step7')
            flow3.communicate()
            if flow3.returncode != 0:
                cout += 1
        if GPU == 1:
            flow3 = subprocess.Popen(
                'gmx mdrun -v -deffnm md -pin on -ntmpi 1 -ntomp ' + str(CPU) +
                ' -pme gpu -nb gpu',
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd='./gromacs/step7')
            flow3.communicate()
            if flow3.returncode != 0:
                cout += 1
    else:
        flow3 = subprocess.Popen('gmx mdrun -v -deffnm md -ntmpi 1 -ntomp ' +
                                 str(CPU) +
                                 ' -pin on -pinstride 1 -pinoffset 0',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True,
                                 cwd='./gromacs/step7')
        flow3.communicate()
        if flow3.returncode != 0:
            cout += 1

    if cout != 0:
        return False
    else:
        return True
