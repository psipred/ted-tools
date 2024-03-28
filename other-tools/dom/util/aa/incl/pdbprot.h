#ifndef PDBPROT_H
#define PDBPROT_H

/* ==== HEADER pdbprot.h ==== */

/* Protein Data Bank I/O routines. Replaces "pdb.h" */

/* ANSI C, IRIX 4.0.5, 9. Nov. 1995. Andris Aszodi */

/* ---- GLOBAL DEFINITIONS ---- */

#define ALLATOMS 0     /* read all atoms */
#define CALPHA 1       /* read C-alpha backbone only */
#define RELAXED 0      /* read non-standard AA-s ('X') */
#define STRICT 1       /* read the 20 standard AA-s + B,Z only */

/* ---- GLOBAL PDB TYPES ---- */

typedef char Str4_[5];  /* 4-letter words... */

typedef struct    /* entry for an atom */
{
  int Atno;	/* atom serial number as in the PDB entry */
  Str4_ Id;  /* atom type like CA or OD2: 4 chars max. + \0 */
  char Alt;	/* alternate location indicator */
  char Aa;     /* 1-letter AA code or ' ' for OXT */
  int Resno;  /* residue no. within a chain */
  char Rid;   /* insertion code like in "27A" */
  float X,Y,Z; /* the coordinates */
  float Occu, Bfact; /* occupancy and temperature factor */
} Atom_ ;

typedef struct	    /* entry for a H-bond */
{
    int Don, Acc;   /* atom nos of donor and acceptor */
} Hbond_ ;

typedef enum {HELIX, SHEET, TURN} Sectype_ ;
typedef struct	    /* entry for a secondary structure element */
{
    Sectype_  Sectype;	/* HELIX, SHEET, TURN */
    int No;	    /* serial number of structure */
    Str4_ Id;	    /* identifier string */
    int Beg, End;   /* residue nos (1..N) */
    char Chid;	    /* chain ID, for begin only */
    char Begaa, Endaa;	/* amino acid 1-letter codes */
    char Begrid, Endrid;    /* insertion codes */
    int Type;	    /* 1..7 for helices, -1..1 for sheets */
	    /* these entries are for sheets only */
    int Strandno;	/* number of strands in this sheet */ 
    Str4_ Thisat, Otherat;  /* atom names for beta registration */
    char Thisaa, Otheraa;   /* amino acid codes for beta regs */
    int This, Other;	/* registration AA positions */
    char Thisrid, Otherid, Otherchid;	/* ins codes and chain ID for the OTHER */
} Secstr_ ;

typedef struct	    /* entry for a disulfide bridge */
{
    int No;
    int Pos1, Pos2; /* CYS seq. positions */
    char Ch1, Ch2, Rid1, Rid2;	/* chain ID-s and insertion codes */
} Ssbond_ ;

typedef struct   /* entry for a chain */
{
  Atom_ *Atoms; /* array of atoms in this chain */
  int Atomno;   /* no. of atoms in chain */
  Secstr_ *Secs;    /* secondary structure info */
  int Secsno;	/* length of Secs[] */
  Hbond_ *Hbonds;   /* array of H-bonds */
  int Hbno;	/* length of Hbonds[] */
  Ssbond_ *Ssbs;    /* array of S-S bonds */
  int Ssbno;	/* length of Ssbs[] */
  int Aano;     /* no. of amino acids in chain */
  char Chid;    /* chain ID (for multichains) */
  char Type;    /* chain type: 'P'-rotein, 'A'-lpha or 'X' */
  char *Seq;	/* chain sequence in 1-letter code, \0-termin. */
} Chain_ ;

typedef struct	/* entry for a PDB record */
{
    char Header[41];	/* header info: function description */
    char Date[10];	/* date of deposition: "20-AUG-64" format */
    Str4_ Pdbcode;	/* standard PDB code "0XXX" */
    char *Compound;	/* compound descriptor string */
    char *Source;	/* species, organ, tissue, mutant */
    char Expdta[61];	/* experimental technique: default "X-RAY DIFFRACTION" */
    float Resol;	/* resolution; <0 means "NOT APPLICABLE" */
    Chain_ *Chains;	/* array of chain records */
    int Chainno;	/* no. of chains */
    Hbond_ *Hbonds;	/* array of interchain H-bonds */
    int Hbno;		/* no. of interchain H-bonds */
    Ssbond_ *Ssbs;	/* array of interchain S-S bonds */
    int Ssbno;	/* no, of interchain S-S bonds */
} Pdbentry_ ;

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

char aa_code31(const char *Aa3);
char *aa_code13(char Aa1);
Pdbentry_ *get_pdb(const char *Pdbfn, int Ca,int Strict);
float atom_dist(const Atom_ *At1, const Atom_ *At2);
void put_pdb(const char *Pdbfn, Pdbentry_ *Entry, 
	char *Remarks[], int Remno);
void free_pdb(Pdbentry_ *Entry);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER pdbprot.h ==== */
#endif	/* PDBPROT_H */
