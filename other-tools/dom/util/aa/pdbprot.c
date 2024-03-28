/* ==== FUNCTIONS pdbprot.c ==== */

/* 
#define sqrtf sqrt
#define expf exp
*/

/* Protein Data Bank I/O routines. Replaces "pdb.h" */

/* ANSI C, IRIX 5.3, 21. June 1996. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <unistd.h> /* SEEK_SET on Suns */

/* ---- INCLUDE FILES ---- */

#include "incl/pdbprot.h"	/* own header */

/* ---- DEFINITIONS ---- */

#define LINELEN 2000  /* line buffer length */

/* ---- ALLOCATION MACROS ---- */

/* The following "general alloc" macros are just shorthand for the
 * tedious C memory allocation syntax.
 */
#define G_MALLOC(PTR, TYPE) (PTR)=(TYPE*) malloc(sizeof(TYPE))
#define G_CALLOC(PTR, TYPE, SIZE) (PTR)=(TYPE*) calloc((SIZE), sizeof(TYPE))
#define G_REALLOC(PTR, TYPE, SIZE) \
    (PTR)=(TYPE*) realloc((PTR), (SIZE)*sizeof(TYPE))
#define GFREE(PTR) if ((PTR)!=NULL) { free(PTR); (PTR)=NULL; }

/* ---- FILE-SCOPE GLOBALS ---- */

static long Lastpos=0L;	/* filepos where last good line was found */

/* ---- PROTOTYPES ---- */

static char *atom_seq(const Atom_ Atoms[], int Atomno, int Aano);
static int ascend_res(const void *Sec1, const void *Sec2);
static int get_record(FILE *Pdb, char *Line, const char *Label);
static void print_card(FILE *Pdb, char *Card, const char *Code, int Linecnt);
static int secstr_outorder(const void *Ssb1, const void *Ssb2);
static int ssbond_outorder(const void *Ssb1, const void *Ssb2);

/* ==== FUNCTIONS ==== */

/* aa_code31: converts 3-letter all-uppercase amino acid abbreviations
 * into 1-char codes. Returns 'X' if it cannot translate Aa3.
 * In addition to the 20 standard AA-s, the following translations
 * are accepted:
 * ALB -> 'A' (beta-alanine)
 * CYH, CSH, CSS, CYX -> 'C' (cysteine and cystine)
 * HYP, PR0, PRZ -> 'P' (hydroxyproline and various Pro synonyms)
 * ILU -> 'I' (isoleucine)
 * TRY-> 'W' (tryptophan)
 * UNK -> 'U' (unknown)
 * any other code is translated to X.
 */
char aa_code31(const char *Aa3)
{
    static char *Names3="ALAALBARGASNASPASXCYSCYHCSHCSSCYXGLNGLUGLXGLYHISILEILULEULYSMETPHEPROPR0PRZHYPSERTHRTRPTRYTYRVALUNKXXX";
    static char *Letters1="AARNDBCCCCCQEZGHIILKMFPPPPSTWWYVUX";
    char *Code;
    
    if (strlen(Aa3)<3) return('X'); /* 1,2-letter codes not accepted */
    Code=strstr(Names3, Aa3);
    if (Code==NULL) return('X');    /* not found */
    else return(Letters1[(Code-Names3)/3]);
}
/* END of aa_code31 */   

/* aa_code13: converts the 1-letter codes to 3-letter standards.
 * Almost the inverse of aa_code31(). Non-standard translations:
 * J->JJJ, O->OOO, U->UNK, X->XXX
 */
char *aa_code13(char Aa1)
{
    static char *Aacode[26]={"ALA","ASX","CYS","ASP","GLU","PHE",
	    "GLY","HIS","ILE","JJJ","LYS","LEU",
	    "MET","ASN","OOO","PRO","GLN","ARG",
	    "SER","THR","UNK","VAL","TRP","XXX",
	    "TYR","GLX"};
    
    if (Aa1<'A' || Aa1>'Z') return(Aacode[23]);	/* 'X' for out-of-range */
    else return(Aacode[Aa1-'A']);
}
/* END of aa_code13 */

/* ---- PDB INPUT ---- */

/* get_pdb: reads a PDB file called Pdbfn and puts the result
 * in a Pdbentry_ which is created inside and returned.
 * If something went wrong then NULL is returned.
***
 * If Ca=1 then only the C-alpha atoms are read in without alternate positions;
 * If Ca=2 then only the C-alpha atoms are read in including alternate positions;
***
 * if Strict!=0 then only the 20 standard amino acids are read.
 */
#define CHUNK 256
Pdbentry_ *get_pdb(const char *Pdbfn, int Ca,int Strict)
{
    Pdbentry_ *Entry=NULL;	    /* the complete PDB entry */
    FILE *Pdb;
    Secstr_ *Secs=NULL;  /* temp array for sec. structure elements */
    Ssbond_ *Ssbs=NULL;	    /* temp array for S-S bonds */
    Chain_ *Chains=NULL;    /* temp array for chains */
    Atom_ *Atoms=NULL;    /* temp array for atom records */
    Hbond_ *Hbonds=NULL;	    /* non-redundant D->A storage for H-bonds */
    char Line[LINELEN], Endmdline[LINELEN];      /* read buffers */
    char *Cptr=NULL;	/* aux. char pointer */
    char Aa3[4];       /* 3-letter AA code */
    /* alternative pos, res ID char like in "27A", prev chain ID */
    char Oldres,Oldrid,Oldchain, Creatchn,Endf,Protein,Calpha,Thisca, Model, Ter;       
    int i,j, Sno, Atno, Hbno, Ssbno, Partner, Chainno,Aano,Oldaano, Lastatno=0, 
	Lowlim, Uplim, At, Don, Acc;
    long Oldlastpos=0L, Modelpos=0L,  Terpos=0L;

    /* open PDB file */
    if (NULL==(Pdb=fopen(Pdbfn, "r")))
    {
	fprintf(stderr, "? get_pdb(): Cannot open %s\n", Pdbfn);
	return(NULL);
    }
    
    /* go to the beginning of file */
    rewind(Pdb); Lastpos=0L;
    
    /* allocate the entry record */
    G_MALLOC(Entry, Pdbentry_);
    if (Entry==NULL)
    {
	fputs("? get_pdb(): Out of memory\n", stderr);
	return(NULL);
    }
    
    /* do some zeroing */
    Entry->Chains=NULL; Entry->Hbonds=NULL; Entry->Ssbs=NULL;
    Entry->Chainno=Entry->Hbno=Entry->Ssbno=0;
    
    /* get the HEADER card */
    if (get_record(Pdb, Line, "HEADER"))
    {
	strncpy(Entry->Header, Line+10, 40); Entry->Header[40]='\0';
	strncpy(Entry->Date, Line+50, 9); Entry->Date[9]='\0';
	strncpy(Entry->Pdbcode, Line+62, 4); Entry->Pdbcode[4]='\0';
    }
    else Entry->Header[0]=Entry->Date[0]=Entry->Pdbcode[0]='\0';   /* missing HEADER */
    
    /* get the COMPND cards */
    Entry->Compound=NULL;
    for (i=0; get_record(Pdb, Line, "COMPND"); i++)
    {
	if (i==0) { /* just read store ist line WRT */
	G_REALLOC(Entry->Compound, char, 60*(i+1)+1);
	sscanf(Line+10, "%60c", Entry->Compound+60*i);
	Entry->Compound[60*(i+1)]='\0';	/* not appended when "%nc" is scanned... */
    } }
    if (Entry->Compound==NULL)	    /* no COMPND cards */
	G_CALLOC(Entry->Compound, char, 1); /* empty string */
    
    /* get the SOURCE cards */
    Entry->Source=NULL;
    for (i=0; get_record(Pdb, Line, "SOURCE"); i++)
    {
	G_REALLOC(Entry->Source, char, 60*(i+1)+1);
	sscanf(Line+10, "%60c", Entry->Source+60*i);
	Entry->Source[60*(i+1)]='\0';
    }
    if (Entry->Source==NULL)	    /* no SOURCE cards */
	G_CALLOC(Entry->Source, char, 1); /* empty string */
    
    /* get the experimental technique */
    Entry->Expdta[0]='\0';
    if (!get_record(Pdb, Line, "EXPDTA"))
	strcpy(Entry->Expdta, "X-RAY DIFFRACTION"); /* default technique */
    else
	sscanf(Line+10, "%60c", Entry->Expdta);
    
    /* get resolution from REMARK 2, store -1.0 for "NOT APPLICABLE" */
    Entry->Resol=-1.0;
    while (get_record(Pdb, Line, "REMARK"))
    {
	if (strncmp(Line, "REMARK   2", 10)) continue;	/* not second remark */
	Cptr=strstr(Line, "RESOLUTION.");
	if (Cptr==NULL) continue;
	if (NULL!=strstr(Line, "NOT APPLICABLE"))
	{ Entry->Resol=-1.0; break; }
	if (sscanf(Cptr+11, "%f", &(Entry->Resol)))   /* found resolution */
	    break;
    }
    
    /* get secondary structure records en masse */
    Sno=0;
    while (get_record(Pdb, Line, "HELIX"))	/* read HELIX entries */
    {
	G_REALLOC(Secs, Secstr_, Sno+1);
	Secs[Sno].Sectype=HELIX;
	sscanf(Line+7, "%d %s", &(Secs[Sno].No), Secs[Sno].Id);
	Secs[Sno].Chid=Line[19];
	sscanf(Line+15, "%3s", Aa3);
	Secs[Sno].Begaa=aa_code31(Aa3);
	sscanf(Line+21, "%d%c", &(Secs[Sno].Beg), &(Secs[Sno].Begrid));
	sscanf(Line+27, "%3s", Aa3);
	Secs[Sno].Endaa=aa_code31(Aa3);
	sscanf(Line+33, "%d%c", &(Secs[Sno].End), &(Secs[Sno].Endrid));
	sscanf(Line+38, "%d", &(Secs[Sno].Type));
	Sno++;
    }
    while (get_record(Pdb, Line, "SHEET"))	/* read SHEET entries */
    {
	G_REALLOC(Secs, Secstr_, Sno+1);
	Secs[Sno].Sectype=SHEET;
	sscanf(Line+7, "%d %s%d", &(Secs[Sno].No), Secs[Sno].Id, &(Secs[Sno].Strandno));
	Secs[Sno].Chid=Line[21];
	sscanf(Line+17, "%3s", Aa3);
	Secs[Sno].Begaa=aa_code31(Aa3);
	sscanf(Line+22, "%d%c", &(Secs[Sno].Beg), &(Secs[Sno].Begrid));
	sscanf(Line+28, "%3s", Aa3);
	Secs[Sno].Endaa=aa_code31(Aa3);
	sscanf(Line+33, "%d%c", &(Secs[Sno].End), &(Secs[Sno].Endrid));
	sscanf(Line+38, "%d", &(Secs[Sno].Type));
	
	/* read registration if this is not the first strand */
	if (Secs[Sno].Type)
	{
	    /* quite a number of PDB entries (e.g. 1ACX) contain no
	     * information for the beta-strand registrations. .Other
	     * will be set to 0 to indicate this sad condition and
	     * a warning printed to stderr
	     */
	    Secs[Sno].Other=0;
	    
	    /* atom name, residue name, pos and insert ID for the
	     * registration point on THIS strand. For atom names, 
	     * see the hack comment at the ATOM record read below
	     */
	    if (Line[41]==' ') sscanf(Line+42, "%3s", Secs[Sno].Thisat);
	    else sscanf(Line+41, "%4s", Secs[Sno].Thisat);
	    sscanf(Line+45, "%3s",Aa3);
	    Secs[Sno].Thisaa=aa_code31(Aa3);
	    sscanf(Line+50, "%d%c", &(Secs[Sno].This), &(Secs[Sno].Thisrid));

	    /* same + chain ID for the OTHER strand */
	    if (Line[56]==' ') sscanf(Line+57, "%3s", Secs[Sno].Otherat);
	    else sscanf(Line+56, "%4s", Secs[Sno].Otherat);
	    sscanf(Line+60, "%3s", Aa3);
	    Secs[Sno].Otheraa=aa_code31(Aa3);
	    sscanf(Line+64, "%c%d%c", &(Secs[Sno].Otherchid), 
		    &(Secs[Sno].Other), &(Secs[Sno].Otherid));
	}
	else	/* first strand, zero the This/Other fields */
	{
	    Secs[Sno].Thisaa=Secs[Sno].Thisrid=Secs[Sno].Thisat[0]='\0';
	    Secs[Sno].Otheraa=Secs[Sno].Otherid=Secs[Sno].Otherat[0]='\0';
	    Secs[Sno].This=Secs[Sno].Other=0;
	    Secs[Sno].Otherchid='\0';
	}
	Sno++;
    }
    while (get_record(Pdb, Line, "TURN"))	/* read TURN entries */
    {
	G_REALLOC(Secs, Secstr_, Sno+1);
	Secs[Sno].Sectype=TURN;
	sscanf(Line+7, "%d %s", &(Secs[Sno].No), Secs[Sno].Id);
	sscanf(Line+15, "%3s", Aa3);
	Secs[Sno].Chid=Line[19];
	Secs[Sno].Begaa=aa_code31(Aa3);
	sscanf(Line+20, "%4d%c", &(Secs[Sno].Beg), &(Secs[Sno].Begrid));
	sscanf(Line+26, "%3s", Aa3);
	Secs[Sno].Endaa=aa_code31(Aa3);
	sscanf(Line+31, "%4d%c", &(Secs[Sno].End), &(Secs[Sno].Endrid));
	Secs[Sno].Type=0;
	Sno++;
    }
    
    /* sort secondary structure records in ascending residue order */
    if (Sno)
	qsort(Secs, Sno, sizeof(Secstr_), ascend_res);
    
    /* read S-S bond records in bulk */
    for (Ssbno=0; get_record(Pdb, Line, "SSBOND"); Ssbno++)
    {
	G_REALLOC(Ssbs, Ssbond_, Ssbno+1);
	sscanf(Line+7, "%3d", &(Ssbs[Ssbno].No));
	sscanf(Line+15, "%c%4d%c", 
	    &(Ssbs[Ssbno].Ch1), &(Ssbs[Ssbno].Pos1), &(Ssbs[Ssbno].Rid1));
	sscanf(Line+29, "%c%4d%c", 
	    &(Ssbs[Ssbno].Ch2), &(Ssbs[Ssbno].Pos2), &(Ssbs[Ssbno].Rid2));
    }
    
    /* read chains and atom coordinates: since TER is sometimes missing,
     * the chain boundaries are also detected when the chain ID char
     * changes between two subsequent records. Chain boundaries
     * may also be indicated by ENDMDL cards (NMR structures).
     */
    Chainno=0; i=0; Model=get_record(Pdb, Line, "MODEL");
    Oldchain='\0';
    
    Oldlastpos=Lastpos;	    /* find next TER if any */
    if (Ter=get_record(Pdb, Line, "TER"))   /* = intended */
    { Terpos=Lastpos; Lastpos=Oldlastpos; }
    
    if (Model)
    {	    /* find next ENDMDL */
	Oldlastpos=Lastpos;
	if(get_record(Pdb, Line, "ENDMDL"))
	{ Modelpos=Lastpos; Lastpos=Oldlastpos; }
    }
    
    do  /* read ATOM records */
    {
	/* find first line beginning with "ATOM" */
	Endf=!get_record(Pdb, Line, "ATOM");
	
	/* start new chain if the chain indicator has changed
	 * or a TER or ENDMDL card was passed
	 */
	if (!Endf) 
	{
	    Creatchn=(Oldchain!=Line[21]);
	    if (Ter && Lastpos>Terpos)	    /* TER has been passed */
	    {
		/* get next TER if any */
		Oldlastpos=Lastpos; Lastpos=Terpos;
		if (get_record(Pdb, Endmdline, "TER"))
		    Terpos=Lastpos;
		else Ter=0;	/* don't check any more */
		Lastpos=Oldlastpos;
		Creatchn=1;
	    }
	    
	    if (Model && Lastpos>Modelpos)
	    {	    /* find next ENDMDL */
		Oldlastpos=Lastpos; Lastpos=Modelpos;
		if(get_record(Pdb, Endmdline, "ENDMDL"))
		    Modelpos=Lastpos;
		else Model=0;	/* no more checks */
		Lastpos=Oldlastpos;
		Creatchn=1;
	    }
	}

	/* save prev chain before exiting or starting a new one */
	if (Oldchain && (Endf || Creatchn))
	{
	    Chains[Chainno].Atomno=i; /* save no. of atoms */
	    Chains[Chainno].Aano=Aano; /* no. of residues */
	    Chains[Chainno].Chid=Oldchain; /* chain ID */
	    if (!Protein) Chains[Chainno].Type='X';  /* not protein */
	    else if (Calpha) Chains[Chainno].Type='A'; /* C-alpha only */
	    else Chains[Chainno].Type='P';  /* proper protein */

	    Chains[Chainno].Atoms=       /* data pointer */
	    (Atom_ *) realloc(Atoms,i*sizeof(Atom_));
	    Chains[Chainno].Seq=atom_seq(Chains[Chainno].Atoms, i, Aano);
	    
	    /* H-bond and secstr info will come later */
	    Chains[Chainno].Hbonds=NULL; Chains[Chainno].Secs=NULL;
	    Chains[Chainno].Ssbs=NULL;
	    Chains[Chainno].Hbno=Chains[Chainno].Secsno=
		Chains[Chainno].Ssbno=0;
	    Chainno++;
	}

	if (Endf) break; /* exit */

	if (Creatchn)  /* start new chain */
	{
	    G_REALLOC(Chains, Chain_, Chainno+1);
	    G_CALLOC(Atoms, Atom_, CHUNK);
	    i=Aano=0; Oldres=Oldrid=' '; Oldchain=Line[21]; Oldaano=-9999;
	    Protein=0; Calpha=1;
	}

	/* make room for Atoms[] for further growth */
	if (i && !(i % CHUNK))
	    G_REALLOC(Atoms, Atom_, (i/CHUNK+1)*CHUNK);
	
	/* start processing line */
	sscanf(Line,"ATOM%d", &(Atoms[i].Atno));
	
	/* an ugly hack: in the PDB format, atom names are built of
	 * 4 chars: 
	 * 1-2=chem.symbol,  RIGHT justified;
	 * 3=remoteness (alpha, optional)
	 * 4=branch designator (numeric, optional)
	 * When the 1st pos (Line[12]) is blank, we read in a 3-char
	 * string only in order not to pick up the (optional) Alt character. 
	 * Otherwise a 4-char string is read. My .Id strings are
	 * always left-justified internally
	 */
	if (Line[12]==' ') sscanf(Line+13, "%3s", Atoms[i].Id);
	else sscanf(Line+12, "%4s", Atoms[i].Id);
	
	Atoms[i].Alt=Line[16]; 

	/* skip non-C-alphas when the Ca parameter is set */
	Thisca=(!strcmp(Atoms[i].Id,"CA"));	/* this is C-alpha */
	if ((Ca>0) && !Thisca) continue;
	if ((Atoms[i].Alt != ' ') && (Atoms[i].Alt != 'A')) {
		if (Ca==1 || Ca==0) continue; 
	}
	sscanf(Line+17,"%3s",Aa3);	/* get AA name */
	Atoms[i].Aa=aa_code31(Aa3);  /* convert to 1-letter code */
	/* accepted as protein if at least one AA is not 'X' */
	/* else skip non-standard AA atoms if Strict is on */
	if (Atoms[i].Aa!='X') Protein=1;
	else if (Strict) continue;

	/* not C-alpha if at least one atom is not "CA" */
	if (!Thisca) Calpha=0;

	/* read residue number */
	sscanf(Line+22,"%4d",&(Atoms[i].Resno));

	/* read coordinates, occupancy and B-factor */
	sscanf(Line+27,"%f%f%f%f%f",
	    &(Atoms[i].X),&(Atoms[i].Y),&(Atoms[i].Z),
	    &(Atoms[i].Occu), &(Atoms[i].Bfact));

	/* some residues are not numbered consecutively, like in
	1TIM (1,2,4,5,...) or there are more residues with the
	same number but a different residue ID like 27,27A,27B...
	The following if() is intended to take care of these
	anomalies. Note that in 1ABB, the A-chain has the
	following "order": Met-679, Lys-910, Phe-681....
	This is a nightmare. */
	Atoms[i].Rid=Line[26];
	if (Atoms[i].Resno!=Oldaano || Oldrid!=Atoms[i].Rid || Oldres!=Atoms[i].Aa)
	{                            
	    Oldaano=Atoms[i].Resno;
	    Oldrid=Atoms[i].Rid;  
	    Oldres=Atoms[i].Aa;
	    Aano++;
	}
	Lastatno=Atoms[i].Atno;  /* last "good" atom seen so far */
	i++; /* next atom: not inc'd if Strict skips */
    }
    while (1);     /* end of do..while : break if Endf */
    
    /* read CONECT records for H-bonds */
    Hbno=0;
    while (get_record(Pdb, Line, "CONECT"))
    {
	/* the atom no. at positions 6..10 (Atno) can have acceptors
	 * at positions 31..35 and 36..40, and donors at positions
	 * 46..50, 51..55. The rest of the CONECT record is ignored.
	 * H-bonds that connect
	 * atoms outside the 1..Lastatno range will be ignored.
	 */
	sscanf(Line+6, "%d", &Atno);
	if (Atno>Lastatno) continue;	/* outside chain atom range */
	
	At=Atno;
	Line[41]=Line[56]='\0';
	
	/* read 1 acceptor and 1 donor only */
	if(1==sscanf(Line+31, "%d", &Partner) && Partner<=Lastatno)
	    Acc=Partner; else Acc=0;
	if(1==sscanf(Line+46, "%d", &Partner) && Partner<=Lastatno)
	    Don=Partner; else Don=0;
	/* store the corresponding pairs */
	if (Don>0)
	{
	    for (j=0; j<Hbno; j++)
		if (Don==Hbonds[j].Don && At==Hbonds[j].Acc) break;
	    if (!Hbno || j>=Hbno)
	    {
		G_REALLOC(Hbonds, Hbond_, Hbno+1);
		Hbonds[Hbno].Don=Don; Hbonds[Hbno].Acc=At;
		Hbno++;
	    }
	}
	if (Acc>0)
	{
	    for (j=0; j<Hbno; j++)
		if (At==Hbonds[j].Don && Acc==Hbonds[j].Acc) break;
	    if (!Hbno || j>=Hbno)
	    {
		G_REALLOC(Hbonds, Hbond_, Hbno+1);
		Hbonds[Hbno].Don=At; Hbonds[Hbno].Acc=Acc;
		Hbno++;
	    }
	}
    }
    fclose(Pdb);
    
    /* separate secondary structure records by chains */
    for (i=0; i<Sno; i++)
    {
	/* locate the chain where the i-th secstr record belongs:
	 * a given secstr may belong to separate chains as in
	 * NMR structural models e.g. 3CI2
	 */
	for (j=0; j<Chainno; j++)
	    if (Chains[j].Chid==Secs[i].Chid) /* found, make room,store */
	    {
		G_REALLOC(Chains[j].Secs, Secstr_, Chains[j].Secsno+1);
		Chains[j].Secs[Chains[j].Secsno++]=Secs[i];
	    }
    }
    GFREE(Secs);
    
    /* separate S-S bonds by chains and build interchain S-S bonds */
    for (i=0; i<Ssbno; i++)
    {
	if (Ssbs[i].Ch1!=Ssbs[i].Ch2)	/* interchain */
	{
	    G_REALLOC(Entry->Ssbs, Ssbond_, Entry->Ssbno+1);
	    Entry->Ssbs[Entry->Ssbno++]=Ssbs[i];
	}
	else	    /* intrachain */
	{
	    /* locate the chain(s) with the same ID */
	    for (j=0; j<Chainno; j++)
		if (Chains[j].Chid==Ssbs[i].Ch1)
		{
		    G_REALLOC(Chains[j].Ssbs, Ssbond_, Chains[j].Ssbno+1);
		    Chains[j].Ssbs[Chains[j].Ssbno++]=Ssbs[i];
		}
	}
    }
    GFREE(Ssbs);
    
    /* separate H-bonds by chains and build interchain H-bonds */
    for (i=0; i<Hbno; i++)
    {
	/* locate chain: check if the H-donors (or acceptors)
	 * fall within the range of Atno indices in the j-th
	 * chain. If not, then the bond is interchain
	 */
	Don=Hbonds[i].Don; Acc=Hbonds[i].Acc;
	for (j=0; j<Chainno; j++)
	{
	    Lowlim=Chains[j].Atoms[0].Atno;
	    Uplim=Chains[j].Atoms[Chains[j].Atomno-1].Atno;
	    if (Don>=Lowlim && Don<=Uplim &&
		Acc>=Lowlim && Acc<=Uplim)
		break;
	}
	if (j<Chainno)	    /* intra-chain */
	{
	    /* make room */
	    G_REALLOC(Chains[j].Hbonds, Hbond_, Chains[j].Hbno+1);
	    Chains[j].Hbonds[Chains[j].Hbno++]=Hbonds[i];
	}
	else		    /* inter-chain */
	{
	    /* make room */
	    G_REALLOC(Entry->Hbonds, Hbond_, Entry->Hbno+1);
	    Entry->Hbonds[Entry->Hbno++]=Hbonds[i];
	}
    }
    GFREE(Hbonds);
    
    /* cleanup and final data transfer */
    Entry->Chains=Chains; Entry->Chainno=Chainno;
    return(Entry);
}
/* END of get_pdb */

/* get_record: returns the first record (in Line) from Pdb which
 * begins with the string Label. Line is assumed to be longer than
 * LINELEN. Stores the file position where the last matching
 * line was found in Lastpos (global).
 * Return value: 1 if a matching record was found, 0 if no more
 * records were found. In this latter case, Line=="".
 */
static int get_record(FILE *Pdb, char *Line, const char *Label)
{
    int Lablen=strlen(Label);
    
    if (ftell(Pdb)!=Lastpos)
	fseek(Pdb, Lastpos, SEEK_SET);    /* go where last was found */
    Line[0]='\0';	/* zero line */
    while (NULL!=fgets(Line, LINELEN, Pdb))
	if (!strncmp(Line, Label, Lablen))
	{ int i, l = strlen(Line);
	    if (l<72) for (i=l-1; i<72; i++) Line[i] = ' ';
	    Lastpos=ftell(Pdb); 	/* found */
	    Line[72]='\0';	/* chop off chars 73-80 "0XXX1234" */
	    return(1);
	}
    return(0);	    /* not found */
}
/* END of get_record */

/* ascend_res: for sorting Secstr_ records in increasing order.
 * Alphabetical ordering for chain identifiers, and increasing
 * residue number (beginning) order within the same chain.
 */
static int ascend_res(const void *Sec1, const void *Sec2)
{
    register Secstr_ *S1=(Secstr_ *)Sec1, *S2=(Secstr_ *)Sec2;
    
    if (S1->Chid<S2->Chid) return(-1);
    else if (S1->Chid==S2->Chid) 
    {
	if (S1->Beg<S2->Beg) return(-1);
	else if (S1->Beg>S2->Beg) return(1);
	else return(0);
    }
    else return(0);
}
/* END of ascend_res */

/* atom_seq: given a PDB chain atom array together with its
 * length Atomno and the number of amino acids Aano, the 
 * sequence (in 1-letter code) is returned in a \0-terminated string.
 */
static char *atom_seq(const Atom_ Atoms[], int Atomno, int Aano)
{
  register int i,j,Resno;
  char *Seq;
  char Rid;
  
  Seq= (char *) calloc(Aano+1,sizeof(char));
 
  for (j=0,i=0; i<Atomno; j++)
  {
    Seq[j]=Atoms[i].Aa;
    Rid=Atoms[i].Rid; Resno=Atoms[i++].Resno;
    while (i<Atomno && Resno==Atoms[i].Resno && Rid==Atoms[i].Rid) i++;
  }
  return(Seq);
}
/* END of atom_seq */

/* free_pdb: frees up memory allocated to Entry and its arrays.
 * Does not set Entry to NULL. (Destructor in C++)
 */
void free_pdb(Pdbentry_ *Entry)
{
    register int i;
    
    if (Entry==NULL) return;
    
    for (i=0; i<Entry->Chainno; i++)
    {
	if (Entry->Chains[i].Aano) free(Entry->Chains[i].Seq);
	if (Entry->Chains[i].Atomno) free(Entry->Chains[i].Atoms);
	if (Entry->Chains[i].Secsno) free(Entry->Chains[i].Secs);
	if (Entry->Chains[i].Hbno) free(Entry->Chains[i].Hbonds);
	if (Entry->Chains[i].Ssbno) free(Entry->Chains[i].Ssbs);
    }
    if (Entry->Chainno) free(Entry->Chains);
    if (Entry->Ssbno) free(Entry->Ssbs);
    if (Entry->Hbno) free(Entry->Hbonds);
    free(Entry->Compound); free(Entry->Source);
    free(Entry);
}
/* END of free_pdb */

/* atom_dist: returns the distance between two Atom_ entries. */
float atom_dist(const Atom_ *At1, const Atom_ *At2)
{
    float Dx,D;

    D=0.0;
    Dx=At1->X-At2->X;
    D+=Dx*Dx;
    Dx=At1->Y-At2->Y;
    D+=Dx*Dx;
    Dx=At1->Z-At2->Z;
    D+=Dx*Dx;
    return(sqrtf(D));
}
/* END of atom_dist */

/* ---- OUTPUT ---- */

/* put_pdb: writes the atomic coordinate information in Entry
 * to the file Pdbfn in PDB format. Currently supported records are
 * HEADER, COMPND, SOURCE, REMARK, SEQRES, HELIX, SHEET, TURN, 
 * SSBOND, ATOM, CONECT. Remno remarks are supplied in a char array.
 */
void put_pdb(const char *Pdbfn, Pdbentry_ *Entry, 
	char *Remarks[], int Remno)
{
    /* 3-atom H-bond record type for CONECT output */
    typedef struct
    {
	int Don, At, Acc;   /* Don or Acc may be -1 */
    }
    Hb3_ ;
    
    /* index for the counters in the MASTER record */
    typedef enum {REM, FTN, HET, HLX, SHT, TRN, 
		SIT, TRS, ATM, TER, CON, SQR} Masteridx_ ;
    
    FILE *Pdb;		/* output file */
    char Card[LINELEN];	/* output line str buffer */
    Chain_ *Chains;	/* temp ptrs to Entry's (sub)fields */
    Atom_ *Atoms;
    Hbond_ *Hbonds=NULL;    /* merged list of Don->Acc H-bonds */
    Hb3_ *Hb3s=NULL;	/* list of Don->At->Acc output H-bonds */
    Secstr_ *Secs=NULL;	/* merged list of secondary structures */
    Ssbond_ *Ssbs=NULL;	/* merged list of S-S bonds */
    register int i, j, k, se, ss, hb;
    int Chno, Atomno, Aano, Hbno, Hb3no, Ssbno, Secsno, Don, Acc, Linecnt, Atomcnt;
    char *Seq, *Code, *Cmp, Chid, Dside, Aside;
    int Mc[12]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};	/* record counter */
    Masteridx_ Mx;

    /* check data quickly and open file */
    if (Entry==NULL || !Entry->Chainno) return;    /* possibly empty */
    if (NULL==(Pdb=fopen(Pdbfn, "w")))
    {
	fprintf(stderr, "? put_pdb: Cannot open %s\n", Pdbfn);
	return;
    }
    
    Linecnt=1;
    /* compose header */
    sprintf(Card, "HEADER    %-40s%9s   %4s",
	Entry->Header, Entry->Date, Code=Entry->Pdbcode);
    print_card(Pdb, Card, Code, Linecnt++);

    /* compound description (may span several lines) */
    for (i=1, Cmp=Entry->Compound, k=strlen(Cmp); k>0; i++, Cmp+=60, k-=60)
    {
	if (i==1) sprintf(Card, "COMPND    %-60.60s", Cmp);
	else sprintf(Card, "COMPND  %2d%-60.60s", i, Cmp);
	print_card(Pdb, Card, Code, Linecnt++);
    }

    /* source description (may span several lines) */
    for (i=1, Cmp=Entry->Source, k=strlen(Cmp); k>0; i++, Cmp+=60, k-=60)
    {
	if (i==1) sprintf(Card, "SOURCE    %-60.60s", Cmp);
	else sprintf(Card, "SOURCE  %2d%-60.60s", i, Cmp);
	print_card(Pdb, Card, Code, Linecnt++);
    }

    /* experimental technique */
    sprintf(Card, "EXPDTA    %-60.60s", Entry->Expdta);
    print_card(Pdb, Card, Code, Linecnt++);
    
    /* insert a warning and print resolution */
    sprintf(Card, "REMARK   1 NOT A GENUINE PDB ENTRY!");
    print_card(Pdb, Card, Code, Linecnt++);
    if (Entry->Resol<0)
	strcpy(Card, "REMARK   2 RESOLUTION. NOT APPLICABLE.");
    else
	sprintf(Card, "REMARK   2 RESOLUTION. %.2f ANGSTROMS.", Entry->Resol);
    print_card(Pdb, Card, Code, Linecnt++);
    Mc[REM]=2;
    
    /* add user-supplied remarks starting at 4 because 3 is reserved for refinement */
    for (j=0; j<Remno; j++, Mc[REM]++)
    {
	if (strlen(Remarks[j])>=62) continue;	/* skip longs */
	sprintf(Card,"REMARK   4 %s",Remarks[j]);
	print_card(Pdb, Card, Code, Linecnt++);
    }

    /* list sequence info */
    Chains=Entry->Chains;
    for (Chno=0; Chno<Entry->Chainno; Chno++)
    {
	Seq=Chains[Chno].Seq;
	Aano=Chains[Chno].Aano; Chid=Chains[Chno].Chid;
	for (i=0; 13*i<Aano; Linecnt++,i++, Mc[SQR]++)
	{
	    sprintf(Card, "SEQRES%4d %c %4d ", i+1, Chid, Aano);
	    for (j=0; j<13 && 13*i+j<Aano; j++)
		sprintf(Card+18+4*j, "%4s", aa_code13(Seq[13*i+j]));
	    print_card(Pdb, Card, Code, Linecnt);
	}
    }

    /* the secstr, S-S bond and H-bond records are distributed
     * among the chains and have to be pooled for output. The
     * records are counted here, temp space is allocated, and
     * after copying the records are sorted
     */
    Secsno=0; Ssbno=Entry->Ssbno; Hbno=Entry->Hbno;
    for (Chno=0; Chno<Entry->Chainno; Chno++)
    {
	Secsno+=Chains[Chno].Secsno;
	Ssbno+=Chains[Chno].Ssbno;
	Hbno+=Chains[Chno].Hbno;
    }
    G_CALLOC(Secs, Secstr_, Secsno);
    G_CALLOC(Ssbs, Ssbond_, Ssbno);
    G_CALLOC(Hbonds, Hbond_, Hbno);
    for (se=ss=hb=0, Chno=0; Chno<Entry->Chainno; Chno++)
    {
	memcpy(Secs+se, Chains[Chno].Secs, 
	    sizeof(Secstr_)*Chains[Chno].Secsno);
	se+=Chains[Chno].Secsno;
	memcpy(Ssbs+ss, Chains[Chno].Ssbs, 
	    sizeof(Ssbond_)*Chains[Chno].Ssbno);
	ss+=Chains[Chno].Ssbno;
	memcpy(Hbonds+hb, Chains[Chno].Hbonds, 
	    sizeof(Hbond_)*Chains[Chno].Hbno);
	hb+=Chains[Chno].Hbno;
    }
    memcpy(Ssbs+ss, Entry->Ssbs, sizeof(Ssbond_)*Entry->Ssbno);
    memcpy(Hbonds+hb, Entry->Hbonds, sizeof(Hbond_)*Entry->Hbno);
    qsort(Secs, Secsno, sizeof(Secstr_), secstr_outorder);
    qsort(Ssbs, Ssbno, sizeof(Ssbond_), ssbond_outorder);
    
    /* list sec. struct. info */
    for (i=0; i<Secsno; i++)
	if (Secs[i].Sectype==HELIX)	    /* helices */
	{
	    sprintf(Card, "HELIX  %3d %3s %3s %c %4d%c %3s %c %4d%c%2d", 
		Secs[i].No, Secs[i].Id, aa_code13(Secs[i].Begaa), 
		Secs[i].Chid, Secs[i].Beg, Secs[i].Begrid, 
		aa_code13(Secs[i].Endaa), Secs[i].Chid, Secs[i].End,
		Secs[i].Endrid, Secs[i].Type);
	    print_card(Pdb, Card, Code, Linecnt++); Mc[HLX]++;
	}
	else if(Secs[i].Sectype==SHEET)	/* sheets */
	{
	    sprintf(Card, "SHEET  %3d %3s%2d %3s %c%4d%c %3s %c%4d%c%2d ", 
		Secs[i].No, Secs[i].Id, Secs[i].Strandno,
		aa_code13(Secs[i].Begaa), Secs[i].Chid,
		Secs[i].Beg, Secs[i].Begrid, 
		aa_code13(Secs[i].Endaa), Secs[i].Chid, Secs[i].End,
		Secs[i].Endrid, Secs[i].Type);
		
	    /* print registration if present */
	    if (Secs[i].Type && Secs[i].Other)
	    {
		sprintf(Card+41, (strlen(Secs[i].Thisat)==4)? "%4s": " % -3.3s", 
		    Secs[i].Thisat);
		sprintf(Card+45, "%3s %c%4d%c ",  
		    aa_code13(Secs[i].Thisaa),
		    Secs[i].Chid, Secs[i].This, Secs[i].Thisrid);
		sprintf(Card+56, (strlen(Secs[i].Otherat)==4)? "%4s": " % -3.3s", 
		    Secs[i].Otherat);
		sprintf(Card+60, "%3s %c%4d%c",  
		    aa_code13(Secs[i].Otheraa),
		    Secs[i].Otherchid, Secs[i].Other, Secs[i].Otherid);
	    }
	    print_card(Pdb, Card, Code, Linecnt++); Mc[SHT]++;
	}
	else	/* turns */
	{
	    sprintf(Card, "TURN   %3d %3s %3s %c%4d%c %3s %c%4d%c", 
		Secs[i].No, Secs[i].Id, aa_code13(Secs[i].Begaa), 
		Secs[i].Chid, Secs[i].Beg, Secs[i].Begrid, 
		aa_code13(Secs[i].Endaa), Secs[i].Chid, Secs[i].End,
		Secs[i].Endrid);
	    print_card(Pdb, Card, Code, Linecnt++); Mc[TRN]++;
	}
    GFREE(Secs);
    
    /* list disulfide bridges */
    for (i=0; i<Ssbno; i++)
    {
	sprintf(Card, "SSBOND %3d CYS %c %4d%c   CYS %c %4d%c", 
	    Ssbs[i].No, 
	    Ssbs[i].Ch1, Ssbs[i].Pos1, Ssbs[i].Rid1, 
	    Ssbs[i].Ch2, Ssbs[i].Pos2, Ssbs[i].Rid2);
	print_card(Pdb, Card, Code, Linecnt++);
    }
    GFREE(Ssbs);
    
    /* list atomic coordinates for each chain (ATOM, TER) */
    Atomcnt=1;
    for (Chno=0; Chno<Entry->Chainno; Chno++)
    {
	Aano=Chains[Chno].Aano; Chid=Chains[Chno].Chid;
	Atoms=Chains[Chno].Atoms; Atomno=Chains[Chno].Atomno;
	for (i=0; i<Atomno; Atomcnt++,i++)
	{
	    sprintf(Card, "ATOM  %5d ", Atoms[i].Atno);
	    sprintf(Card+12, (strlen(Atoms[i].Id)==4)? "%4s": " %-3.3s", 
		Atoms[i].Id);
	    sprintf(Card+16, "%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f",
		Atoms[i].Alt, aa_code13(Atoms[i].Aa),
		Chid,Atoms[i].Resno,Atoms[i].Rid, 
		Atoms[i].X,Atoms[i].Y,Atoms[i].Z,
		Atoms[i].Occu,Atoms[i].Bfact);
	    print_card(Pdb, Card, Code, Linecnt++); Mc[ATM]++;
	}
	sprintf(Card, "TER   %5d      %3s %c%4d",
	    Atomcnt++,aa_code13(Atoms[Atomno-1].Aa),
	    Chid,Atoms[Atomno-1].Resno); Mc[TER]++;
	print_card(Pdb, Card, Code, Linecnt++);
    }
    
    /* construct the 3-member Don->At->Acc Hb3_ records from the general
     * Don->Acc Hbond_ records. This is to comply with the
     * rather awkward CONECT format which lists atoms rather
     * than bonds and is redundant
     */
    for (Hb3no=i=0; i<Hbno; i++)
    {
	Dside=Aside=1;
	Don=Hbonds[i].Don; Acc=Hbonds[i].Acc;
	for (j=0; j<Hb3no && (Dside || Aside); j++)
	{
	    if (Dside && Acc==Hb3s[j].At && !Hb3s[j].Don)
	    {	/* join to j-th on the donor side */
		Hb3s[j].Don=Don; Dside=0;
	    }
	    else Dside=Dside && !(Don==Hb3s[j].Don && Acc==Hb3s[j].At);
	    if (Aside && Don==Hb3s[j].At && !Hb3s[j].Acc)
	    {	/* join to j-th on the acceptor side */
		Hb3s[j].Acc=Acc; Aside=0;
	    }
	    else Aside=Aside && !(Don==Hb3s[j].At && Acc==Hb3s[j].Acc);
	}
	
	/* store if unpaired still */
	if (Dside)
	{
	    G_REALLOC(Hb3s, Hb3_, Hb3no+1);
	    Hb3s[Hb3no].Don=Don; Hb3s[Hb3no].At=Acc; 
	    Hb3s[Hb3no++].Acc=0; 
	}
	if (Aside)
	{
	    G_REALLOC(Hb3s, Hb3_, Hb3no+1);
	    Hb3s[Hb3no].Don=0; Hb3s[Hb3no].At=Don; 
	    Hb3s[Hb3no++].Acc=Acc; 
	}
    }
    GFREE(Hbonds);

    /* list H-bonds from Hb3s[] */
    for (i=0; i<Hb3no; i++)
    {
	memset(Card, 0, LINELEN);
	/* atom no and blank covalent bonds (not std but...) */
	sprintf(Card, "CONECT %4d", Hb3s[i].At+1);
	memset(Card+11, ' ', 20);
	if (!Hb3s[i].Acc)
	    memset(Card+31, ' ', 15);
	else sprintf(Card+31, "%5d          ", Hb3s[i].Acc);
	if (!Hb3s[i].Don)
	    memset(Card+46, ' ', 15);
	else sprintf(Card+46, "%5d", Hb3s[i].Don);
	print_card(Pdb, Card, Code, Linecnt++); Mc[CON]++;
    }
    GFREE(Hb3s);
    
    /* finish listing with MASTER and END */
    strcpy(Card, "MASTER    ");
    for (Mx=REM; Mx<=SQR; Mx++)
	sprintf(Card+10+5*Mx, "%5d", Mc[Mx]);
    print_card(Pdb, Card, Code, Linecnt++);
    strcpy(Card, "END");
    print_card(Pdb, Card, Code, Linecnt);
    fclose(Pdb);
}
/* END of put_pdb */

/* print_card: completes a PDB output line ("card" in FORTRAN terminology).
 * The incomplete line is in Card, to which padding spaces, the Code
 * and the line count Linecnt are appended to make it LINELEN(==80)
 * chars long. Then Card is then printed to Pdb followed by a newline.
 */
static void print_card(FILE *Pdb, char *Card, const char *Code, int Linecnt)
{
    register int i;
    
    i=strlen(Card);
    memset(Card+i, ' ', 72-i);	/* pad with spaces */
    sprintf(Card+72, "%4s%4d", Code, Linecnt);
    fprintf(Pdb, "%s\n", Card);
}
/* END of print_card */

/* ---- ENTRY SORTING FOR OUTPUT ---- */

/* secstr_outorder: sorts the secstruct elements into PDB output order:
 * helix-sheet-turn, and same types are sorted
 * separately according to the order number (No).
 */
static int secstr_outorder(const void *Sec1, const void *Sec2)
{
    register Secstr_ *S1=(Secstr_ *)Sec1, *S2=(Secstr_ *)Sec2;
    register int Ord;
    
    Ord=(int)(S1->Sectype)-(int)(S2->Sectype);
    if (Ord) return(Ord);	/* helix-sheet-turn */
    else if (S1->Sectype==SHEET)    /* both are sheets */
    {
	Ord=strcmp(S1->Id, S2->Id); /* sort by alphabetical ID order */
	if (Ord) return(Ord);
	else return(S1->No-S2->No); /* same sheet, different strands */
    }
    else return(S1->No-S2->No);	    /* both helices or turns */
}
/* END of secstr_outorder */

/* ssbond_outorder: for sorting the S-S bond records according
 * their No field.
 */
static int ssbond_outorder(const void *Ssb1, const void *Ssb2)
{
    register Ssbond_ *S1=(Ssbond_ *)Ssb1, *S2=(Ssbond_ *)Ssb2;
    return(S1->No-S2->No);
}
/* END of ssbond_outorder */

#undef LINELEN
#undef G_MALLOC
#undef G_CALLOC
#undef G_REALLOC
#undef GFREE

/* ==== END OF FUNCTIONS pdbprot.c ==== */

