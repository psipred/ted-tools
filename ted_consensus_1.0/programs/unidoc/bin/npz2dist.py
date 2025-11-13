#!/public/home/sry/opt/miniconda3/envs/bio/bin/python
import sys
import numpy as np
import heapq,os
import argparse

parser=argparse.ArgumentParser(
    add_help=False,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='NPZ file to distance matrix',
    epilog="""An example:\npython npz2dist.py -n seq.npz -o the output file\n """
    )
arghelp=parser.add_argument_group('help information')
arghelp.add_argument("-h", "--help",action="help",help="show this message and exit")

argrequired=parser.add_argument_group('mandatory arguments')
argrequired.add_argument("-n",dest='NPZ', required=True, type=str,help="npz file")
argrequired.add_argument("-o",dest='OUT', required=True, type=str, help="the output file")

args=parser.parse_args()

npz = args.NPZ
out = args.OUT # output path




labels = ['dist','omega','theta','phi']
nbins = [36,24,24,12]


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    if(len(values)==1):return values[0],0

    average = np.average(values, weights=weights)
    variance=np.array(values).std()
    return average, variance


def dist(npz):
    dat=npz['dist']
    nres=int(dat.shape[0])
    pcut=0.05
    bin_p=0.01
    mat=np.zeros((nres, nres))
    for i in range(0, nres):
        for j in range(0, nres):
            if(j == i):
                mat[i][i]=0
                continue
            if(j<i):
                mat[i][j]=mat[j][i]
                continue

            #check probability
            Praw = dat[i][j]#37 bins
            first_bin=2
            first_d=2.25 #4-->3.75, 5-->4.25
            weight=0

            for P in Praw[first_bin:]:
            #from 5th - last
                if(P>bin_p): weight += P
            if(weight < pcut):
                mat[i][j]=20
                continue
            Pnorm = [P for P in Praw[first_bin:]]
            probs=[]
            dists=[]
            Plargest_index = list(map(Pnorm.index,heapq.nlargest(5,Pnorm)))
            e_dis=8;
            e_std=0;
            sum = 0
            for k in Plargest_index:
                d = first_d + k*0.5
                if Pnorm[k] > bin_p:
                    sum += Pnorm[k]
                    probs.append(Pnorm[k])
                    dists.append(d)
                if sum > 0.70:
                    break
            prob = [P/sum for P in probs]
            e_dis, e_std=weighted_avg_and_std(dists, prob)#加权距离
            mat[i][j] = e_dis
    return (mat)
# #def dist
#Distance
def main():
    npz_file = np.load(npz)
    bname = os.path.basename(npz) # base name of npz file
    name = bname.split('.')[0]
    img=dist(npz_file)
    np.savetxt(out,img,fmt='%f',delimiter='\t') 


if __name__ == "__main__":
    main()


