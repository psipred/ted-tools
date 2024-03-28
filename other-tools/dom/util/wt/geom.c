#include "incl/util.h"
#include "incl/geom.h"

vinit (c)		Vec	*c;
{
	c->x = 0.0;
	c->y = 0.0;
	c->z = 0.0;
}

vcopy (b, c)		Vec	b, *c;
{
	c->x = b.x;
	c->y = b.y;
	c->z = b.z;
}

vunit (b, c)		Vec	b, *c;
{	float	d = DIST3v((b));
	c->x = b.x/d;
	c->y = b.y/d;
	c->z = b.z/d;
}

vnorm (c)		Vec	*c;
{	float	d = DIST3v((*c));
	c->x = c->x/d;
	c->y = c->y/d;
	c->z = c->z/d;
}

vave (a, b, c)	Vec	a, b, *c;
{
	c->x = 0.5 * (a.x + b.x);
	c->y = 0.5 * (a.y + b.y);
	c->z = 0.5 * (a.z + b.z);
}	 

vsum (a, c)	Vec	a, *c;
{
	c->x += a.x;
	c->y += a.y;
	c->z += a.z;
}	 

vadd (a, b, c)	Vec	a, b, *c;
{
	c->x = a.x + b.x;
	c->y = a.y + b.y;
	c->z = a.z + b.z;
}	 

vsub (a, b, c)	Vec	a, b, *c;
{
	c->x = a.x - b.x;
	c->y = a.y - b.y;
	c->z = a.z - b.z;
}	 

vmul (c, s)		Vec	*c;
			float	s;
{
	c->x = (c->x)*s;
	c->y = (c->y)*s;
	c->z = (c->z)*s;
}

vdiv (c, s)		Vec	*c;
			float	s;
{
	c->x = (c->x)/s;
	c->y = (c->y)/s;
	c->z = (c->z)/s;
}

float	vdif (a, b)	Vec	a, b;
{	Vec	c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vddif (a, b)	Vec	a, b;
{	Vec	c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return (c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vdad (a, b)	Vec	a, b;
{	Vec	c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vddad (a, b)	Vec	a, b;
{	Vec	c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return (c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vsqr (a)	Vec	a;
{
	return (a.x*a.x + a.y*a.y + a.z*a.z);
}	 

float	vmod (a)	Vec	a;
{
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}	 

float	vdot (a, b)	Vec	a, b;
{
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}	 

vprod (a, b, c)	Vec	a, b, *c;
{
      c->x = a.y*b.z - b.y*a.z;
      c->y = a.z*b.x - b.z*a.x;
      c->z = a.x*b.y - b.x*a.y;
}

float	vtri (a, b, c)	Vec	a, b, c;
{	Vec	d;
	vprod(a,b,&d);
	return vdot(c,d);
}	 

float	pdotp (a, b, c, d)	Vec	a, b, c, d;
{	Vec	x, y;
	vsub(b,a,&x);
	vsub(d,c,&y);
	return vdot(x,y);
}	 

float	pvol (a, b, c, d)	Vec	a, b, c, d;
{	Vec	x, y, z;
	vsub(b,a,&x);
	vsub(d,c,&y);
	vsub(c,b,&z);
	return vtri(x,y,z);
}	 

float	phand (a, b, c, d)	Vec	a, b, c, d;
{	Vec	x, y, z;
	vsub(b,a,&x);
	vsub(d,c,&y);
	vsub(c,b,&z);
	vnorm(&z);
	return vtri(x,y,z);
}	 

VtoM (a,b,c,M) Vec a,b,c; Mat *M; {
        M->A.x = a.x; M->A.y = a.y; M->A.z = a.z;
        M->B.x = b.x; M->B.y = b.y; M->B.z = b.z;
        M->C.x = c.x; M->C.y = c.y; M->C.z = c.z;
}

Mprint (m) Mat *m; {
Mat M = *m;
        printf("\n");
	printf("%9.4f %9.4f %9.4f\n", M.A.x,M.A.y,M.A.z);
	printf("%9.4f %9.4f %9.4f\n", M.B.x,M.B.y,M.B.z);
	printf("%9.4f %9.4f %9.4f\n", M.C.x,M.C.y,M.C.z);
        printf("\n");
}

MmulM (p,q,R)  Mat *p, *q, *R; {
Mat P = *p, Q = *q;
Vec a = P.A, b = P.B, c = P.C;
Vec A = Q.A, B = Q.B, C = Q.C;

	R->A.x = a.x*A.x + b.x*A.y + c.x*A.z;
        R->B.x = a.x*B.x + b.x*B.y + c.x*B.z;
        R->C.x = a.x*C.x + b.x*C.y + c.x*C.z;

	R->A.y = a.y*A.x + b.y*A.y + c.y*A.z;
        R->B.y = a.y*B.x + b.y*B.y + c.y*B.z;
        R->C.y = a.y*C.x + b.y*C.y + c.y*C.z;

	R->A.z = a.z*A.x + b.z*A.y + c.z*A.z;
        R->B.z = a.z*B.x + b.z*B.y + c.z*B.z;
        R->C.z = a.z*C.x + b.z*C.y + c.z*C.z;
}

Minv (m,W,d) Mat *m, *W; float d;
{
Mat M = *m;
Vec a = M.A, b = M.B, c = M.C;

	W->A.x =  (b.y*c.z - b.z*c.y)/d;
        W->A.y = -(a.y*c.z - a.z*c.y)/d;
	W->A.z =  (a.y*b.z - a.z*b.y)/d;

	W->B.x = -(b.x*c.z - b.z*c.x)/d;
        W->B.y =  (a.x*c.z - a.z*c.x)/d;
        W->B.z = -(a.x*b.z - a.z*b.x)/d;

	W->C.x =  (b.x*c.y - b.y*c.x)/d;
        W->C.y = -(a.x*c.y - a.y*c.x)/d;
        W->C.z =  (a.x*b.y - a.y*b.x)/d;
}

float Mdet (m) Mat *m; {
Mat M = *m;
Vec a = M.A, b = M.B, c = M.C;
	return a.x * (b.y*c.z - b.z*c.y)
             - b.x * (a.y*c.z - a.z*c.y)
             + c.x * (a.y*b.z - a.z*b.y);
}

MmulV (m,d,e) Mat *m; Vec d, *e; {
Mat M = *m;
Vec a = M.A, b = M.B, c = M.C;
        e->x = d.x*a.x + d.y*b.x + d.z*c.x;
        e->y = d.x*a.y + d.y*b.y + d.z*c.y;
        e->z = d.x*a.z + d.y*b.z + d.z*c.z;
}

VmulM (m,d,e) Mat *m; Vec d, *e; {
Mat M = *m;
Vec a = M.A, b = M.B, c = M.C;
	e->x = d.x*a.x + d.y*a.y + d.z*a.z;
	e->y = d.x*b.x + d.y*b.y + d.z*b.z;
	e->z = d.x*c.x + d.y*c.y + d.z*c.z;
}

line2tri (a,b,c,d,e) Vec a,b,c,d,e; {
/* TRUE if line segment d-e cuts triangle a-b-c */
Vec x, y, z;
Mat M[1], W[1];
float det;
	vsub(b,a,&x);
	vsub(c,a,&y);
	vsub(d,e,&z);
	VtoM(x,y,z,M);
	det = Mdet(M);
	if (fabs(det) < 0.000001) return 0;
	Minv(M,W,det);
	vsub(d,a,&d);
	MmulV(W,d,&e);
	if (e.x < 0.0) return 0;
	if (e.y < 0.0) return 0;
	if (e.z < 0.0) return 0;
	if (e.z > 1.0) return 0;
	if (e.x + e.y > 1.0) return 0;
	return 1;
}

line2line (a,b,c,d,s) Vec a,b,c,d; float s; {
/* true if line segments a-b and d-e are closer than s */
Vec x, y, z;
Mat M[1], W[1];
float det;
Vec e;
        vsub(b,a,&x);
        vsub(d,c,&y);
        vprod(x,y,&z);
        vnorm(&z);
        VtoM(x,y,z,M);
        det = Mdet(M);
        Minv(M,W,det);
        vsub(d,a,&d);
        MmulV(W,d,&e);
        if (e.x < 0.0) return 0;
        if (e.y < 0.0) return 0;
        if (e.x > 1.0) return 0;
        if (e.y > 1.0) return 0;
        if (e.z >  s ) return 0;
        if (e.z < -s ) return 0;
        return 1;
}

float lineOline (a,b,c,d,box) Vec a,b,c,d, *box; {
/* return overlap length for line segments a-b and c-d */
Vec x, y, z;
Mat M[1], W[1];
float lap, det, aa, ga, gb, hc, hd, r, s, t[4];
Vec e,f,g,h, pox[4];
int i, key[4], ley[4];
	vcopy(a,box+0); vcopy(b,box+1); vcopy(c,box+2); vcopy(d,box+3);
	vsub(b,a,&x); vnorm(&x);
	vsub(d,c,&y); vnorm(&y);
	vprod(x,y,&z); vnorm(&z);
	if (vdot(x,y) < 0.0001) { 
		vcopy(x,&e); vmul(&e,0.0001);
		vsub(c,e,&c); vadd(d,e,&d);
		vsub(d,c,&y); vnorm(&y);
	 }
	VtoM(x,y,z,M);
	det = Mdet(M);
	Minv(M,W,det);
	vsub(d,a,&e);
	MmulV(W,e,&f);
	vmul(&x,f.x);
	vadd(a,x,&g); /* g = top of mut.perp. line */
	vmul(&y,f.y);
	vsub(d,y,&h); /* h = bot of mut.perp. line */
	ga = vdif(g,a); gb = vdif(g,b);
	hc = vdif(h,c); hd = vdif(h,d);
	vsub(a,g,&a); vsub(b,g,&b);
	vsub(c,h,&c); vsub(d,h,&d);
	aa = vsqr(a);
	r = vdot(a,c)/aa;
	s = vdot(a,d)/aa;
	if (vdot(a,b)<0.0) gb = -gb;
	if (r<0.0) hc = -hc;
	if (s<0.0) hd = -hd;
	t[0] = ga; t[1] = gb; t[2] = hc; t[3] = hd; 
	sort(0,t,0,key,4,1);
	for (i=0; i<4; i++) if (key[i]<2) ley[i] = 0; else ley[i] = 1;
	if (ley[0]==ley[1]) return 0.0;
	lap = fabs(t[key[1]]-t[key[2]]);
	if (box)
	{ int ke1 = key[1], ke2 = key[2],
	      le1 = ley[1], le2 = ley[2];
	  Vec tmp1, tmp2, tmp;
		if (ga > gb) vsub(box[0],box[1],&x);
			else vsub(box[1],box[0],&x);
		if (hc > hd) vsub(box[2],box[3],&y);
			else vsub(box[3],box[2],&y);
		vnorm(&x); vnorm(&y);
		vcopy(box[ke1],&tmp1);
		vcopy(box[ke2],&tmp2);
		vcopy(tmp1,box+0);
		if (le1==le2) { /* contained */
			vcopy(tmp2,box+1);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+3);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+3); }
		} else { /* staggered */
			vcopy(tmp2,box+3);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+1);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+1); }
		}
		if (le1==0) {
			 vcopy(y,&tmp); vmul(&tmp,t[ke1]); vadd(h,tmp,box+2);
		} else { vcopy(x,&tmp); vmul(&tmp,t[ke1]); vadd(g,tmp,box+2); }
	}
	/* restore ab cd line order in box (but box lines run parallel */
	if (ley[1]) {
		pox[0] = box[2]; pox[1] = box[3];
		pox[2] = box[0]; pox[3] = box[1];
		for (i=0; i<4; i++) box[i] = pox[i];
	}
	return lap;
}

line2dot (a,b,c,s) Vec a,b,c; float s; {
/* TRUE if c lies over line segment a-b closer than s */
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	if (d < 0.0 || d > 1.0) return 0;
	vmul(&p,d);
	vsub(q,p,&q);
	d = vsqr(q);
	if (d>s*s) return 0;
	return 1;
}

float dotOline (a,b,c,e) Vec a,b,c, *e; {
/* returns distance of c to an extended line a-b and image of c on a-b in e */
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	vmul(&p,d);
	vsub(q,p,&q);
	d = vsqr(q);
	vadd(a,p,&p);
	if (e) { e->x = p.x; e->y = p.y; e->z = p.z; }
	return sqrt(d);
}

float torsion (a,b,c,d)  Vec a, b, c, d; {
/* returns the torsion angle down b-c */
Vec x, y, z, p, q, r;
float vol, cos, sin, tor, cxy, cyz;
	vsub(a,b,&x);
	vsub(b,c,&y);
	vsub(c,d,&z);
	vprod(x,y,&p);
	vnorm(&p);
	vprod(y,z,&q);
	vnorm(&q);
	cos = vdot(p,q);
	vprod(p,q,&r);
	vol = vtri(p,q,y);
	sin = vmod(r);
	if (vol < 0.0) sin = -sin;
	tor = angle1pi(sin,cos);
	return tor;
}

float	angle1pi (s,c) float s, c;
{
	if (s>=0.0 && c>=0.0) {
		if (s<0.5) return asin(s);
		      else return acos(c);
	}
	if (s>=0.0 && c<=0.0) {
		if (s<0.5) return PI - asin(s);
		      else return PI - acos(-c);
	}
	if (s<=0.0 && c<=0.0) {
		if (s>-.5) return asin(-s) - PI;
		      else return acos(-c) - PI;
	}
	if (s<=0.0 && c>=0.0) {
		if (s>-.5) return -asin(-s);
		      else return -acos(c);
	}
	printf("angle1pi(s,c) out of range: s = %f, c = %f\n", s,c);
	exit(1);
}

float	angle2pi (s,c) float s, c;
{
	if (s>=0.0 && c>=0.0) {
		if (s<0.5) return asin(s);
		      else return acos(c);
	}
	if (s>=0.0 && c<=0.0) {
		if (s<0.5) return PI - asin(s);
		      else return PI - acos(-c);
	}
	if (s<=0.0 && c<=0.0) {
		if (s>-.5) return PI + asin(-s);
		      else return PI + acos(-c);
	}
	if (s<=0.0 && c>=0.0) {
		if (s>-.5) return twoPI - asin(-s);
		      else return twoPI - acos(c);
	}
	printf("angle2pi(s,c) out of range: s = %f, c = %f\n", s,c);
	exit(1);
}

float angdif (a,b) float a, b;
{
        if (a<0.0 && b>0.0)
        { float d = b-a;
                if (d<PI) return d;
                     else return twoPI-d;
        }
        if (a>0.0 && b<0.0)
        { float d = a-b;
                if (d<PI) return d;
                     else return twoPI-d;
        }
        return fabs(a-b);
}
