import getopt
import sys
import bin.flow.Link1 as L1
import bin.flow.Link2 as L2
import bin.flow.Link3 as L3
import bin.flow.Link4 as L4
import bin.flow.Link5 as L5
import bin.flow.Link6 as L6
import bin.flow.Link7 as L7
import bin.flow.Link8 as L8

if __name__ == '__main__':

    model = None
    argv = sys.argv[1:]

    opts, args = getopt.getopt(argv, "m:")
    for opt, arg in opts:
        if opt == '-m':
            model = arg

    if model == '1':
        l1 = L1.go()
        if l1:
            l2 = L2.go()
            if l2:
                l3 = L3.go()
                if l3:
                    l4 = L4.go()
                    if l4:
                        l5 = L5.go()
                        if l5:
                            l6 = L6.go()
                            if l6:
                                l7 = L7.go()
                                if l7:
                                    l8 = L8.go()

    if model == '2':
        pass  # TODO: not finish

    if opts == []:
        sys.stdout = sys.stderr
        print("Please select model! (-m 1 or 2)")
        print('1: The PISA (RAFT)')
        print('2: The PISA (ROP)')
        sys.exit(2)
