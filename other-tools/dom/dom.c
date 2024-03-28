/*

cc -O dom.c -o dom util/wt/util.o util/wt/sort.o util/wt/geom.o util/aa/pdbprot.o -lm

*/ 

#include <alloca.h>
#include <string.h>
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"
#include "util/aa/incl/pdbprot.h"

#define NCYCLES 200
#define NALLOC 2000
#define NACID 30

#define MINDOM 40
#define MAXDOM 400

typedef struct {
	Vec	v;
	float	d;
	Vec	cos;
} Tri;

typedef struct {
	int	a, b;
	float	c;
} Pairs;

typedef struct {
	char	*res;
	float	*acc;
	int	*dom;
	int	*rid;
	Vec	*ca, *cb;
	int	len;
} Seq;

double	supermac();
double  drand48();

int	spread, nruns, subdom, sheet;

float	compare();
float	betacut();

#define IT 5

main(argc,argv)
int argc; char *argv[];
{
Tri	**a[1];
Seq	seq[1];
char	file[255];
Pdbentry_ *prot, *prot2;
int	i, j, len, nn;
int	*doms[4], domn[4];
float	**mat, **net, *dom;
float	pct[4][4];
char	tag[4][10];
long    rseed = 234;
        srand48(rseed);
	prot = get_pdb(argv[1],1,1);
	if (argc > 2) {
		sscanf(argv[2],"%d", &spread);
	} else { spread = 15; }
	if (argc > 3) {
		sscanf(argv[3],"%d", &nruns);
	} else { nruns = 1; }
	if (argc > 4) {
		sscanf(argv[4],"%d", &subdom);
	} else { subdom = 0; }
	Pi(spread) Pi(nruns) Pi(subdom) NL
/*						nruns	       -3  -2  -1   0   1   2
						structures	N   S   N   S  NS  NS
						beta-bias	-   -   B   B   -   B
*/
	Ps(prot->Compound) NL
	protin(prot,seq,0,a,1.0,0);
	len = seq->len;
	for (i=1; i<len; i++) { float d = vdif(seq->ca[i],seq->ca[i+1]);
		if ( d > 8.0 ) printf("*NB* chain BREAK at %d/%d\n", i,i+1);
	}
	if (subdom && len<250) {
		printf("Only domains of 250 and over are considered for reparsing\n");
		exit(0);
	}
	nn = len+2;
	mat = (float**)malloc(sizeof(float*)*nn); TEST(mat)
	for (i=0; i<nn; i++) { mat[i] = (float*)malloc(sizeof(float)*nn); TEST(mat[i]) }
	net = (float**)malloc(sizeof(float*)*nn); TEST(net)
	for (i=0; i<nn; i++) { net[i] = (float*)malloc(sizeof(float)*nn); TEST(net[i]) }
	dom = (float*)malloc(sizeof(float)*nn); TEST(dom)
	for (i=0; i<4; i++) { doms[i] = (int*)malloc(sizeof(int)*nn); TEST(doms[i]) }
	sheet = 0;
	if (nruns < 1) { float	cuts;
	        for (i=1; i<=len; i++) vcopy(seq->ca[i], seq->cb+i);
		if (nruns == 0 || nruns == -2) {
        		flatten(seq->cb,len,5);
        		spread += 3;
		}
		if (nruns > -2) beta(net,mat,seq->ca,dom,len);
        		else for (i=1; i<=len; i++) dom[i] = (float)i;
		nruns = 0;
        	Pi(spread) NL
        	for (i=1; i<=IT && define(0,i,seq,seq->cb,dom,net,mat,len, domn); i++);
		for (i=1; i<=len; i++) doms[0][i] = seq->dom[i];
		if (sheet) {
			cuts = betacut (net,doms[0],len);
			printf ("Error in beta split =%7.2f\n", cuts);
		}
		exit(0);
	}
	printf("\n\n===============  normal  ===============\n\n");
	Pi(spread) NL
	if (nruns == 2) beta(net,mat,seq->ca,dom,len);
        	else for (i=1; i<=len; i++) dom[i] = (float)i;
	for (i=1; i<=IT && define(0,i,seq,seq->ca,dom,net,mat,len, domn+0); i++);
	for (i=1; i<=len; i++) doms[0][i] = seq->dom[i];
	printf("\n\n=============== smoothed ===============\n\n");
	for (i=1; i<=len; i++) vcopy(seq->ca[i], seq->cb+i);
	flatten(seq->cb,len,5);
	spread += 3;
	Pi(spread) NL
	if (nruns == 2) beta(net,mat,seq->ca,dom,len);
        	else for (i=1; i<=len; i++) dom[i] = (float)i;
	for (i=1; i<=IT && define(2,i,seq,seq->cb,dom,net,mat,len, domn+2); i++);
	for (i=1; i<=len; i++) doms[2][i] = seq->dom[i];
	NLL
	compare(net, doms[0],doms[2], domn[0],domn[2], len);
}

float betacut (net,doms,n) float **net; int *doms, n;
{
int     i, j, in[200];
float   cut = 0.0;
	for (i=1; i<199; i++) in[i] = 0;
        for (i=1; i<n; i++) in[doms[i]]++;
        for (i=1; i<n; i++) {
                for (j=i+1; j<=n; j++) {
                        if (net[i][j] < 0.0) continue;
                        if (doms[i] == doms[j]) continue;
			if (in[doms[j]] <= MINDOM) continue;
		if (in[doms[i]] <= MINDOM) continue;
                        cut += net[i][j];
                }
        }
	return cut;
}

float	compare (net, doms1, doms2, domn1,domn2, len)
float	**net;
int	*doms1, *doms2, domn1, domn2;
{
int	i, j, match, best, maxdd, **dd;
float	pct, bad;
	maxdd = 0;
	for (i=1; i<=len; i++) {
		if (doms1[i] > maxdd) maxdd = doms1[i];
		if (doms2[i] > maxdd) maxdd = doms2[i];
	}
	maxdd++;
	dd = (int**)alloca(sizeof(int*)*maxdd);
	for (i=0; i<maxdd; i++) dd[i] = (int*)alloca(sizeof(int)*maxdd);
	for (i=0; i<maxdd; i++) for (j=0; j<maxdd; j++) dd[i][j] = 0;
	for (i=1; i<=len; i++) {
		dd[doms1[i]][doms2[i]]++;
	}
	match = 0;
	while (best = getbest(dd,maxdd)) match += best;
	pct = 100.0*(float)match/(float)len;
	if (sheet) { float cuts1, cuts2;
		bad = sqrt((float)len);  /* *2.0 times 2 if 1 added to net */
		cuts1 = betacut(net,doms1,len);
		printf ("Error in beta split 1 =%7.2f\n", cuts1);
		cuts2 = betacut(net,doms2,len);
		printf ("Error in beta split 2 =%7.2f\n", cuts2);
		printf("(error limit =%f7.3)\n", bad);
		if (cuts1 > bad || cuts2 > bad ) pct = -pct;
	}
	printf("%d on %d domains, %3d pct. agree\n", domn1, domn2, (int)pct);
}

getbest (dd,n) int **dd, n;
{
int	i, j, row, col, max = 0;
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (dd[i][j] > max) { max = dd[i][j]; row = j; col = i; }
		}
	}
	for (i=0; i<n; i++) dd[col][i] = -1;
	for (i=0; i<n; i++) dd[i][row] = -1;
	return max;
}

define (is, it, seq, cas, dom, net, mat, len, doms)
int is, it;
Seq *seq; Vec *cas; float *dom, **net, **mat; int len, *doms;
{
int	i, j, k, n, p[NALLOC], in[200], redo, last, mark, ndom;
float	min, max, f;
char	file[255];
FILE	*out;
float	*tmp;
	printf("\n\n**** CYCLE %d ****\n\n", it);
	tmp = (float*)alloca(sizeof(float)*(len+2)); TEST(tmp)
	domain(mat,net,cas,dom,len);
	for (i=1; i<=len; i++) dom[i] = -dom[i];
	dom[0] = 999.9;
	sort(0,dom,0,p,len+1,1);
	seq->dom[p[0]] = n = 0;
	for (i=1; i<199; i++) in[i] = 0;
	for (i=1; i<=len; i++)
	{ float	d1 = dom[p[i-1]],
		d2 = dom[p[i]];
		if (d1-d2 > 1.5) { n++; in[n] = 0; }
		seq->dom[p[i]] = n;
		in[n]++;
	}
	ndom = n;
	f = (float)(n-1);
	for (i=1; i<=len; i++) seq->acc[i] = (float)(seq->dom[i]-1)/f;
	NL
	redo = 0;
	if (n==1) {
		printf("assign   1 domain\n");
		*doms = 1;
		return 0;
	}
	printf("%d domains\n", ndom);
	*doms = 0;
	for (i=1; i<=n; i++) {
		printf("   %5d = %d res.\n", i, in[i]);
		if (in[i] < MINDOM) redo = 1;
			else	(*doms)++;
	}
	if (*doms == 0) { printf ("NO domains\n"); return 0; }
	last = 0;
	if (!redo) last = 1;
	if (it==IT) last = 1;
	mark = 0;
	if (last && !is) mark = 1;
	if (mark && !nruns) {
		printf("assign   %d domains\n", *doms);
		if (subdom) sprintf(file,"dom%d-0.out", subdom);
		       else sprintf(file,"dom0.out");
		out = fopen(file,"w");
		putpdb(seq->ca,seq,out,0);
		fclose(out);
		j = 0;
		for (i=1; i<=n; i++) {
			if (in[i] < MINDOM) continue;
			printf("assign   %5d = %d res.\n", i, in[i]);
			j++;
			if (subdom) sprintf(file,"dom%d-%d.out", subdom, j);
			       else sprintf(file,"dom%d.out", j);
			out = fopen(file,"w");
			putpdb(seq->ca,seq,out,i);
			fclose(out);
		}
	} NL
	n = 1;
        if (mark && !nruns) printf("assign   ");
	printf("segment 1 in domain %d = 1 (%c%d)", seq->dom[1], seq->res[1], seq->rid[1]);
	for (i=2; i<=len; i++) {
		if (seq->dom[i-1] == seq->dom[i]) continue;
		n++;
		printf(" --- %d (%c%d) \n", i-1, seq->res[i-1], seq->rid[i-1]);
        	if (mark && !nruns) printf("assign   ");
		printf("segment %d in domain %d = %d (%c%d)", n, seq->dom[i], i, seq->res[i], seq->rid[i]);
	}
	printf(" --- %d (%c%d)\n", len, seq->res[len], seq->rid[len]);
        if (mark && !nruns) printf("assign   ");
	printf("%d segments\n", n);
	for (i=1; i<=len; i++) {
		if (in[seq->dom[i]] < MINDOM) dom[i] = -1.0;
			else dom[i] = 100.0*(float)(seq->dom[i]);
	}
	return redo;
}

domain (mat,net,ca,ave,n)
float **mat, **net; Vec *ca; float *ave; int n;
{
int	i, j, k;
	if (sheet) {
		printf("setting sheet\n");
		ising(net,ave,n,n,0,1);
	}
	printf("defining      ");
	for (i=1; i<n; i++) for (j=i+1; j<=n; j++) mat[i][j] = mat[j][i] = vdif(ca[i],ca[j]);
	for (i=1; i<=n; i++) {
		mat[i][i] = 0.0;
		for (j=1; j<=n; j++)
		{ float d = mat[i][j];
			if (i==j) continue;
			mat[i][j] = (float)spread/d;
			if (d>(float)spread) mat[i][j] = -1.0;
				else	mat[i][i] += 1.0;
			if (ave[i]<0.0 && ave[j]<0.0) mat[i][j] = -1.0;
		}
	}
	ising(mat,ave,n,n,0,0);
	smooth(ave,n,21);
}

beta (net,mat,ca,ave,n) float **net, **mat; Vec *ca; float *ave; int n;
{
int	i, j, k, m, nn;
float	cut = 7.5;
	NLL
        for (i=1; i<=n; i++) ave[i] = (float)i;
	for (i=1; i<n; i++) for (j=i+1; j<=n; j++) mat[i][j] = mat[j][i] = vdif(ca[i],ca[j]);
	for (i=0; i<n+2; i++) for (j=0; j<n+2; j++) net[i][j] = 0.0;
	for (i=2; i<n; i++)
	{ int	p[3], q[3], r, s[3];
	  float dmin;
		p[0]=p[1]=p[2]=q[0]=q[1]=q[2] = 0;
		s[0]=i-1, s[1]=i, s[2]=i+1;
		dmin = 999.9;
		for (j=2; j<n; j++)
		{ float	a = mat[i][j];
			if (abs(i-j)<4) continue;
			if (a > cut) continue;
			if (a > dmin) continue;
			p[1]=j; dmin=a;
		}
		dmin = 999.9;
		for (j=2; j<n; j++)
		{ float	a = mat[i][j];
			if (abs(i-j)<4) continue;
			if (a > cut) continue;
			if (a > dmin) continue;
			if (abs(j-p[1])<6) continue;
			q[1]=j; dmin=a;
		}
		if (!p[1] || !q[1]) continue;
		if (fmin(mat[s[0]][p[1]-1],mat[s[2]][p[1]+1])
		  < fmin(mat[s[2]][p[1]-1],mat[s[0]][p[1]+1])) {
			p[0]=p[1]-1; p[2]=p[1]+1;
		} else {
			p[0]=p[1]+1; p[2]=p[1]-1;
		}
		if (fmin(mat[s[0]][q[1]-1],mat[s[2]][q[1]+1])
		  < fmin(mat[s[2]][q[1]-1],mat[s[0]][q[1]+1])) {
			q[0]=q[1]-1; q[2]=q[1]+1;
		} else {
			q[0]=q[1]+1; q[2]=q[1]-1;
		}
		if (mat[s[0]][p[0]] > cut) continue;
		if (mat[s[2]][p[2]] > cut) continue;
		if (mat[s[0]][q[0]] > cut) continue;
		if (mat[s[2]][q[2]] > cut) continue;
		if (abs(s[2]-p[2])<4) continue;
		if (abs(s[0]-p[0])<4) continue;
		if (abs(s[2]-p[0])<4) continue;
		if (abs(s[0]-p[2])<4) continue;
		if (abs(s[2]-q[2])<4) continue;
		if (abs(s[0]-q[0])<4) continue;
		if (abs(s[2]-q[0])<4) continue;
		if (abs(s[0]-q[2])<4) continue;
		if (abs(q[2]-p[2])<4) continue;
		if (abs(q[0]-p[0])<4) continue;
		if (abs(q[2]-p[0])<4) continue;
		if (abs(q[0]-p[2])<4) continue;
		printf("sheet: %4d %4d %4d\n", p[1],s[1],q[1]);
		/* link across chain */
		net[s[0]][q[0]] += 1.0; net[q[0]][s[0]] += 1.0;
		net[s[0]][p[0]] += 1.0; net[p[0]][s[0]] += 1.0;
		net[s[1]][q[1]] += 1.0; net[q[1]][s[1]] += 1.0;
		net[s[1]][p[1]] += 1.0; net[p[1]][s[1]] += 1.0;
		net[s[2]][q[2]] += 1.0; net[q[2]][s[2]] += 1.0;
		net[s[2]][p[2]] += 1.0; net[p[2]][s[2]] += 1.0;
		/* link along chain */
		net[s[1]][s[0]] += 1.0; net[s[0]][s[1]] += 1.0;
		net[s[1]][s[2]] += 1.0; net[s[2]][s[1]] += 1.0;
		net[p[1]][p[0]] += 1.0; net[p[0]][p[1]] += 1.0;
		net[p[1]][p[2]] += 1.0; net[p[2]][p[1]] += 1.0;
		net[q[1]][q[0]] += 1.0; net[q[0]][q[1]] += 1.0;
		net[q[1]][q[2]] += 1.0; net[q[2]][q[1]] += 1.0;
		/* extra links along chain 
		net[s[0]][s[0]-1] += 1.0; net[s[0]-1][s[0]] += 1.0;
		net[s[2]][s[2]+1] += 1.0; net[s[2]+1][s[2]] += 1.0;
		*/
	}
	NL
	sheet = 0;
	for (i=1; i<=n; i++) {
		for (j=1; j<=n; j++) {
			if (net[i][j]>0.5) {
				net[i][i] += 1.0;
/*
				net[j][i] += 1.0;
				net[i][j] += 1.0;
*/
				sheet = 1;
			} else {
				if (i != j) net[i][j] = -1.0;
			}
		}
	}
}

ising (mat, ave, n, cycles, limit, smooth)
float **mat, *ave;
int	n, cycles, limit, smooth;
{
int	i, j, k;
float	*tmp, rms, sum, shift, step, mean, wt;
int	in, jumps, **list;
float	*old, *new;
	old = (float*)alloca(sizeof(float)*(n+1)); TEST(old)
	new = (float*)alloca(sizeof(float)*(n+1)); TEST(new)
	for (j=1; j<=n; j++) old[j] = new[j] = ave[j];
	list = (int**)alloca(sizeof(int*)*(n+1)); TEST(list)
	for (i=1; i<=n; i++) {
		list[i] = (int*)alloca(sizeof(int)*(int)(mat[i][i]+1.5)); TEST(list[i])
		k = 1;
		for (j=1; j<=n; j++) {
			if (i==j) continue;
			if (mat[i][j] < 0.0) continue;
			list[i][k] = j;
			k++;
		}
		list[i][0] = k;
	}
	shift = 1.0;
	step = 2.0/(float)n;
	for (k=1; k<cycles; k++) { int ii;
		for (j=1; j<=n; j++) {
			in = 0;
			mean = sum = wt = 0.0;
			for (ii=1; ii<list[j][0]; ii++) {
				i = list[j][ii];
                                if (limit && abs(i-j) > limit ) continue;
				if (old[i]>old[j]) sum += mat[i][j];
				if (old[i]<old[j]) sum -= mat[i][j];
				if (old[i] > 0.0) {
					mean += old[i]*mat[i][j];
					wt += mat[i][j];
					in++;
				}
			}
			if (sum>0.0) new[j] += shift;
			if (sum<0.0) new[j] -= shift;
			if (in && (old[j]<0.0 || smooth)) new[j] = mean/wt;
		}
		rms = 0.0;
		jumps = 0;
		for (j=1; j<=n; j++) { float a;
			if (j>1 && abs(new[j]-new[j-1]) > 2) jumps++;
			a = ave[j];
			ave[j] = 0.5*(new[j]+old[j]);
			a -= ave[j];
			rms += a*a;
			old[j] = new[j];
		}
		rms = sqrt(rms/(float)n);
		/* printf("%d %5d %f\n", k, jumps, rms); */
		if (rms<0.000001 || k>n/2) shift -= step; 
		if (shift < 0.0) break;
	}
	printf("%d %5d %f\n", k, jumps, rms);
}

smooth (dat,n,win) float *dat; int n, win;
{
float	*old, *new;
float	rms;
int	i, j, k;
int	p[50], jumps, w = win/2;
	/* printf("smoothing\n"); */
	old = (float*)alloca(sizeof(float)*(n+win*2+1)); TEST(old)
	new = (float*)alloca(sizeof(float)*(n+win*2+1)); TEST(new)
	for (i=1; i<=n; i++) old[i+win] = dat[i];
	for (i=0; i<=win; i++) old[i] = new[i] = dat[win-i+1];
	for (i=0; i<=win; i++) old[i+n+win+1] = new[i+n+win+1] = dat[n-i];
	for (k=1; k<win; k++) {
		rms = 0.0;
		jumps = 0;
		for (i=w; i<=n+win+w; i++) {
			sort(0,old+i-w,0,p,win,1);
			new[i] = old[i+p[w]-w];
		}
		dat[0] = new[win];
		for (i=1; i<=n; i++)
		{ float a = old[i+win]-new[i+win];
			rms += a*a;
			dat[i] = new[i+win];
			if (i && fabs(dat[i]-dat[i-1]) > 2.0) jumps++;
		}
		rms = sqrt(rms/(float)n);
		/* printf("%d %5d %f\n", k, jumps, rms); */
		if (rms<0.000001) break;
		for (i=0; i<n+win*2; i++) old[i] = new[i];
	}
}

flatten (a,n,times) Vec *a; int n, times; 
{
int	i, j;
Vec	p1, p2, p3, q;
	for (i=0; i<times; i++) {
		vcopy(a[1],&p1); vcopy(a[2],&p2); vcopy(a[3],&p3);
		for (j=2; j<n; j++) {
			vinit(&q);
			vsum(p1,&q); vsum(p2,&q); vsum(p3,&q);
			vdiv(&q,3.0);
			vcopy(a[j],&p1); vcopy(a[j+1],&p2); vcopy(a[j+2],&p3);
			vcopy(q,a+j);
		}
	}
}

protin (prot,seq,id,m,z,flip)
Pdbentry_ *prot;
Seq *seq; Tri ***m; float z; int flip, id;
{
FILE	*pdb;
int	i, j, len;
Tri	**mat;
Vec	cog;
	len = copyca(prot->Chains,seq,flip,z);
	vinit(&cog);
	for (i=1; i<=len; i++) vsum(seq->ca[i],&cog);
	vdiv(&cog,(float)len);
	for (i=1; i<=len; i++) vsub(seq->ca[i],cog,seq->ca+i);
        mat = (Tri**)malloc(sizeof(Tri*)*(len+2));
        for (i=0; i<=len+1; i++) {
                mat[i] = (Tri*)malloc(sizeof(Tri)*(len+2));
        }
	add_cb(seq);
	set_vect(seq->ca,seq->cb,mat,len);
	set_cbcb(seq->ca,seq->cb,mat,len);
	*m = mat;
	return len;
}

add_cb (seq)
Seq	*seq;
{       int     i;
        for (i=1; i<=seq->len; i++)
        { Vec   n, c, b;
          float d, bond = 3.0;
                vsub(seq->ca[i],seq->ca[i-1],&n);
                vnorm(&n);
                vsub(seq->ca[i],seq->ca[i+1],&c);
                vnorm(&c);
                vadd(n,c,&b);
                vnorm(&b);
                vmul(&b,bond);
                vadd(seq->ca[i],b,&seq->cb[i]);
        }
}

set_vect (a,b,m,l) Vec *a, *b; Tri **m; int l; {
int     i, j;
Mat     frame;
        for (i=1; i<=l; i++) {
                setframe(a[i-1],a[i],a[i+1],&frame);
                for (j=1; j<=l; j++) { Vec s, t;
                        m[i][j].d = vdif(a[i],a[j]);
                        vinit(&(m[i][j].v));
                        if (i==j) continue;
                        vsub(a[j],a[i],&s);
                        VmulM(&frame,s,&(m[i][j].v));
                }
        }
}

set_cbcb (a,b,m,l) Vec *a, *b; Tri **m; int l; {
int	i, j;
	for (i=1; i<=l; i++)
	{ Vec ai, bi, ci;
		vsub(a[i+1],a[i-1],&ai);
		vnorm(&ai);
		vsub(b[i],a[i],&bi);
		vnorm(&bi);
		vprod(ai,bi,&ci);
		for (j=1; j<=l; j++)
		{ Vec aj, bj, cj;
			vsub(a[j+1],a[j-1],&aj);
			vnorm(&aj);
			vsub(b[j],a[j],&bj);
			vnorm(&bj);
			vprod(aj,bj,&cj);
			m[i][j].cos.x = vdot(ai,aj);
			m[i][j].cos.y = vdot(bi,bj);
			m[i][j].cos.z = vdot(ci,cj);
		}
	}
}

extend (res,i,j,k,new)
Vec	*res;
int	i, j, k, new;
{
	Vec	m, v;
	vave(res[j],res[k],&m);
	vsub(m,res[i],&v);
	vadd(m,v,&res[new]);
}
 
copyca (pdb,s,flip,z)
Chain_  *pdb;
Seq	*s;
int	flip;
float	z;
{	int	i, n;
	char	*seq;
	Vec	*ca, *cb;
	float	*acc;
	int	*dom;
	int	*rid;
	n = pdb->Aano;
	seq = (char*)malloc(sizeof(char)*(n+3));
	acc = (float*)malloc(sizeof(float)*(n+3));
	dom = (int*)malloc(sizeof(float)*(n+3));
	rid = (int*)malloc(sizeof(float)*(n+3));
	ca = (Vec*)malloc(sizeof(Vec)*(n+3));
	cb = (Vec*)malloc(sizeof(Vec)*(n+3));
	for (i=0; i<n; i++) {
		ca[i+1].x = pdb->Atoms[i].X;
		ca[i+1].y = pdb->Atoms[i].Y;
		ca[i+1].z = pdb->Atoms[i].Z;
		acc[i+1] = pdb->Atoms[i].Bfact;
		rid[i+1] = pdb->Atoms[i].Resno;
		seq[i+1] = pdb->Atoms[i].Aa;
		if (seq[i+1]<'A' || seq[i+1]>'Z') {
			printf("*NB* funny aa = %c\n", seq[i+1]);
			seq[i+1] = 'X';
		}
	}
	seq[0] = 'n';
        extend(ca,3,2,1,0);    
        extend(ca,n-2,n-1,n,n+1);
	seq[n+1] = 'c';
	seq[n+2] = 0;
	for (i=0; i<=n+1; i++) ca[i].z *= z;
	if (flip) flipseq(ca,seq,acc,n);
	s->res = seq;
	s->acc = acc;
	s->dom = dom;
	s->rid = rid;
	s->ca = ca;
	s->cb = cb;
	s->len = n;
	return n;
}

flipseq (ca,seq,acc,n) Vec *ca; char *seq; float *acc; int n;
{
int	i;
	for (i=0; i<=n/2; i++)
	{ Vec r; char c; float a;
	  int j = n+1-i;
		r = ca[i]; ca[i] = ca[j]; ca[j] = r;
		c = seq[i]; seq[i] = seq[j]; seq[j] = c;
		a = acc[i]; acc[i] = acc[j]; acc[j] = a;
	}
}
 
getca (res,pdb)
Vec    *res;
FILE	*pdb;
{	int	i = 1;
	char	line[225], junk[30];
        while(!feof(pdb)) {
		read_line(pdb,line);
		if (!strncmp(line,"TER",3)) break;
		if (strncmp(line,"ATOM",4)) continue;
		if (strncmp(line+13,"CA ",3)) continue;
		sscanf(line,"%30c %f%f%f",
                       	junk, &res[i].x, &res[i].y, &res[i].z);
		i++;
	}
	i--;
        extend(res,3,2,1,0);    
        extend(res,i-2,i-1,i,i+1);
	return i;
}
 
putpdb (atom,seq,out,id)
Vec	*atom;
Seq     *seq;
FILE    *out;
int	id;
{
int     ii, i, j, n = 0;
int     len = seq->len;
int	Ndom, Cdom, offset, insert = seq->rid[len]+100;
int	start, end;
char    aa3[80], aaa[4];
        strcpy(aa3,
        "ALAASXCYSASPGLUPHEGLYHISILEACELYSLEUMETASNPCAPROGLNARGSERTHRUNKVALTRPXXXTYRGLX");
	if (subdom) offset = 100+insert; else offset = insert;
	Ndom = Cdom = 0;
	start = 0;
	end = len;
        for (i=len; i>0; i--) if (id && id==seq->dom[i] && seq->rid[i] < insert) { end = i; break; }
        for (i=1; i<=end; i++) {
                strncpy(aaa,aa3+3*(seq->res[i]-'A'),3); aaa[3] = 0;
                if (id && id!=seq->dom[i]) {
			Cdom = i;
			continue;
		}
		if (!start && seq->rid[i] >= insert) continue;
		start = 1;
		if (Ndom && Cdom) { Vec	cent;
			vinit(&cent);
			for (j=Ndom+1; j<Cdom; j++) vsum(atom[j], &cent);
			vdiv(&cent,(float)(Cdom-Ndom-1));
			loop(out,id,&n,atom[Ndom],atom[Cdom],cent,offset);
			Ndom = Cdom = 0;
		}
		n++;
                fprintf(out,"ATOM%7d  CA  %s %c%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n",
                        n, aaa, 'A'+id, seq->rid[i], atom[i].x, atom[i].y, atom[i].z, seq->acc[i]);
		if (id && seq->dom[i] != seq->dom[i+1]) Ndom = i+1;
        }
        fprintf(out,"TER\n");
}

loop (out, id, n, oldN, oldC, cent, offset)
FILE *out;
int	id, *n;
Vec	oldN, oldC, cent;
int	offset;
{
int	i,j;
Vec	newN, newC;
	fprintf(out,"ATOM%7d  CA  UNK %c%4d     %7.3f %7.3f %7.3f 99.99  0.00\n",
                        ++(*n), 'A'+id, (*n)+offset, oldN.x, oldN.y, oldN.z);
	if (vdif(oldN,oldC) > 5.0) { Vec ave;
		vsub(cent,oldN,&newN);	vsub(cent,oldC,&newC);
		vnorm(&newN);		vnorm(&newC);
		vmul(&newN,3.8);	vmul(&newC,3.8);
		vadd(oldN,newN,&newN);	vadd(oldC,newC,&newC);
		vave(newN,newC,&ave);
		if (vdif(newN,newC) < 3.0) {
			fprintf(out,"ATOM%7d  CA  UNK %c%4d     %7.3f %7.3f %7.3f 99.99  0.00\n",
                        	++(*n), 'A'+id, (*n)+offset, ave.x, ave.y, ave.z);
		} else {
			vave(ave,cent,&cent);
			loop(out,id,n,newN,newC,cent,offset);
		}
	}
	fprintf(out,"ATOM%7d  CA  UNK %c%4d     %7.3f %7.3f %7.3f 99.99  0.00\n",
                        ++(*n), 'A'+id, (*n)+offset, oldC.x, oldC.y, oldC.z);
}

setframe (a, b, c, frame)
    Vec a, b, c;
    Mat *frame;
{
    int    i;
    Vec    x, y, z ;
	vsub(c,a,&x);
	vave(c,a,&c);
	vsub(c,b,&y);
	vprod(y,x,&z);
	vprod(z,x,&y);
	vnorm(&x);
	vnorm(&y);
	vnorm(&z);
	VtoM(x,y,z,frame);
}

