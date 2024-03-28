#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <strings.h>
#include <stdlib.h>
#include <alloca.h>

#define LOWER(c)	((c)<'a'?((c)-'A'+'a'):(c))
#define UPPER(c)	((c)<'a'?(c):((c)-'a'+'A'))
#define DIST3v(v)	sqrt( (v).x*(v).x + (v).y*(v).y + (v).z*(v).z )
#define DDIST3v(v)	(v).x*(v).x + (v).y*(v).y + (v).z*(v).z
#define DIST3(x,y,z)	sqrt((x)*(x)+(y)*(y)+(z)*(z))
#define PRINTi(i)     printf(" i = %d",i);
#define Pi(i) printf(" " #i " = %d",i);
#define PRINTr(r)     printf(" r = %f",r);
#define Pr(r) printf(" " #r " = %f",r);
#define PRINTc(x)     printf(" x = %c",x);
#define Pc(x) printf(" " #x " = %c",x);
#define PRINTs(x)     printf(" x = %s",x);
#define Ps(x) printf(" " #x " = %s",x);
#define	PRINTv(v)	printf(" v = %f %f %f", v.x, v.y, v.z); 
#define	Pv(v)	printf(" " #v " = %f %f %f", v.x, v.y, v.z); 
#define	PRINTe(v)	printf(" v = %f %f %f  %f", v.x, v.y, v.z, v.e); 
#define DO(i,n)		for ( i = 1; i <= n; i++)
#define SP		printf(" ");
#define SPP		printf("  ");
#define SPPP		printf("   ");
#define NL		printf("\n");
#define NLL		printf("\n\n");
#define NLLL		printf("\n\n\n");
#define TEST(obj) 	if (obj == NULL) { printf("malloc fail for " #obj "\n"); exit(1); }
#define REST(obj) 	if (obj == NULL) { printf("realloc fail for " #obj "\n"); exit(1); }

# define YES 1
# define NO  0
# define SET   1
# define UNSET 0
# define LIVE 1
# define DEAD 0
# define TRUE  1
# define FALSE 0
# define LEFT -1
# define RIGHT 1
# define MAX 2147483647
# define MIN -MAX
# define BIG 1000000000
# define WEE -BIG
# define BYTE 255

//float fmin();
//float fmax();
