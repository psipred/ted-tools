#include <stdio.h>
#include <math.h>

#define min(i,j) ((i)<(j) ? (i) : (j))
#define max(i,j) ((i)>(j) ? (i) : (j))

unsigned int pack(i,j,n)	int	i, j, n;
/*	references top half (+diag) of a 2D array (N*N) as a vector	*/
{	unsigned int	ip, jp;
	if(n <= 0) { printf("Bad dimension in IPACK\n"); return -1; }
	if(i <= 0) { printf("Bad I in IPACK: i=%d\n",i); return -2; }
	if(j <= 0) { printf("Bad J in IPACK: j=%d\n",j); return -3; }
	ip = max(i,j);
	jp = min(i,j)-1;
	return (unsigned int)(ip-jp+n*jp-(jp*jp-jp)/2);
}

unpack(ip,jp,n,id)	int *ip, *jp, n; unsigned int id;
/*	references top half (+diag) of a 2D array (N*N) from a vector	*/
{	unsigned int	i, j;
	double 	b = n+0.5,
		a = b*b-2*id;
	if(n>65535) printf("Error in UNPACK: n=%d is too big\n",n);
	if(a<=0.0) printf("Error in UNPACK: n=%d, id=%d\n", n,id);
	j = (unsigned int)(b-sqrt(a)-0.00001);
	i = j+id-n*j+(j*j-j)/2;
	j++;
	if (j<=0 || j>n ) printf("Bad J in UNPACK: j=%d (n=%d id=%d)\n",j,n,id);
	if (i<=0 || i>n ) printf("Bad I in UNPACK: i=%d (n=%d id=%d)\n",i,n,id);
	*ip = min(i,j);
	*jp = max(i,j);
}

read_line(file,string) FILE *file; char *string;
{	char	c;
	int	i=0;
	*string = 0;
	while(c=getc(file)) {
printf("%d >%c<\n", c,c);
		if (feof(file)) return -i-1;
		if (c=='\n') return i;
		string[i] = c;
		string[++i] = 0;
	}
}
next_line(file) FILE *file;
{	char	c;
	while(c=getc(file)) {
		if (feof(file)) return 0;
		if (c=='\n') return 1;
	}
}
