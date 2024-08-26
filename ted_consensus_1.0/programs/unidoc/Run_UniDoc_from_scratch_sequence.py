import math,sys,os,gc,re,time
import argparse

parser=argparse.ArgumentParser(
    add_help=False,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='Unified Domain Cutter for sequence-based domain parsing',
    epilog="""An example:\npython Run_UniDoc_from_scratch_sequence.py -f seq.ss2 -s seq.horiz -n seq.npz \n """
    )
arghelp=parser.add_argument_group('help information')
arghelp.add_argument("-h", "--help",action="help",help="show this message and exit")

argrequired=parser.add_argument_group('mandatory arguments')
argrequired.add_argument("-n",dest='NPZ', required=True, type=str,help="the npz file")
argrequired.add_argument("-f",dest='SEQ', required=True, type=str, help="the *.ss2 file from PSIPRED")
argrequired.add_argument("-s",dest='HOR', required=True, type=str, help="the *.horiz file from PSIPRED")

argoptional=parser.add_argument_group('optional arguments')
argoptional.add_argument("-o",type=str,dest="OUT",default=0,help="the ouput file")

args=parser.parse_args()

npzfile=args.NPZ
ss2=args.SEQ
horiz = args.HOR
out = args.OUT

bindir = './bin'


def logo():
    print("""\
``````_```_```````_`____`````````````````
`````|`|`|`|_`__`(_)``_`\``___```___`````
`````|`|`|`|`'_`\|`|`|`|`|/`_`\`/`__|````
`````|`|_|`|`|`|`|`|`|_|`|`(_)`|`(__`````
``````\___/|_|`|_|_|____/`\___/`\___|````
`````````````````````````````````````````""")
    print("""\
*****************************************
***** UniDoc: Unified Domain Cutter *****
***(for sequence-based domain parsing)***
*****************************************""")


def caculate_ss(ss2,horiz,bindir):
    # convert the psipred format 
    assert os.path.exists(ss2)
    assert os.path.exists(horiz)
    binpath2 = os.path.join(bindir,'format_ss.pl')
    return os.system('%s %s %s'%(binpath2,ss2,horiz))


def npz_to_dist(npzfile,bindir,name):
    binpath = os.path.join(bindir,'npz2dist.py')
    assert os.path.exists(npzfile)
    return os.system('python %s -n %s -o %s.dist.mat'%(binpath,npzfile,name))


def parse_domain(name,bindir,out):
    binpath = os.path.join(bindir,'UniDoc_sequence')
    if out:
        return os.system('%s %s.dist.mat seq.dat > %s'%(binpath,name,out))
    else:
        return os.system('%s %s.dist.mat seq.dat'%(binpath,name))



def main():
    run_code = 0
    logo()
    print("reading input files...")
    bname = os.path.basename(npzfile)
    name = ".".join(bname.split('.')[:-1])

    if run_code == 0:
        # step 1: caculate the secondary structure with PSIPRED
        print("step 1: caculate the secondary structure with PSIPRED")
        run_code += caculate_ss(ss2,horiz,bindir)

    if run_code == 0:
        # step 2: NPZ file to the distance matrix 
        print("step 2: convert NPZ file to the distance matrix ")
        run_code += npz_to_dist(npzfile,bindir,name)

    if run_code == 0:
        # step 3: parse domain 
        print("step 3: parse domain ")
        run_code += parse_domain(name,bindir,out)



if __name__ == "__main__":
    main()