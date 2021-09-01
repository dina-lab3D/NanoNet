
//
// PULCHRA
// Protein Chain Restoration Algorithm
//
// Version 3.06
// December 2007
// Contact: Piotr Rotkiewicz, piotr -at- pirx -dot- com
//
// to compile:
// cc -O3 pulchra.c pulchra_data.c -lm -o pulchra
//
// to use:
// ./pulchra file.pdb
//
// to display available options:
// ./pulchra
//

#define COMPILE_BB
#define COMPILE_ROT

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/timeb.h>

#define uchar unsigned char
#define uint unsigned int
#define real double

#include "pulchra_common.h"

#define PULCHRA_VERSION 3.06
#define MAX_BUF_SIZE 1000

#define FILE_SUCCESS     0
#define FILE_NOT_FOUND  -1
#define FILE_WARNING    -2

#define FATAL_MAE -1

#define FLAG_BACKBONE  1
#define FLAG_CALPHA    2
#define FLAG_SIDECHAIN 4
#define FLAG_SCM       8
#define FLAG_INITIAL  16

#define FLAG_PROTEIN  1
#define FLAG_DNA      2
#define FLAG_RNA      4
#define FLAG_CHYDRO   8

#define RADDEG 180./M_PI
#define DEGRAD M_PI/180.

int _VERBOSE = 0;
int _BB_REARRANGE = 1;
int _BB_OPTIMIZE = 0;
int _CA_OPTIMIZE = 1;
int _CA_RANDOM = 0;
int _CA_ITER = 100;
int _CA_TRAJECTORY = 0;
int _CISPRO = 0;
int _CHIRAL = 1;
int _CENTER_CHAIN = 0;
int _REBUILD_BB = 1;
int _REBUILD_SC = 1;
int _REBUILD_H = 0;
int _PDB_SG = 0;
int _TIME_SEED = 0;
int _XVOLUME = 1;
int _XVOL_ITER = 3;
int _PRESERVE = 1;
real _CA_START_DIST = 3.0;
real _CA_XVOL_DIST = 3.5;
real _SG_XVOL_DIST = 1.6;

#define CALC_C_ALPHA
#define CALC_C_ALPHA_ANGLES
#define CALC_C_ALPHA_START
#define CALC_C_ALPHA_XVOL

real CA_K=10.0;
real CA_ANGLE_K=20.0;
real CA_START_K=0.01;
real CA_XVOL_K=10.00;

#define CA_DIST 3.8
#define CA_DIST_TOL 0.1
#define CA_DIST_CISPRO 2.9
#define CA_DIST_CISPRO_TOL 0.1
#define E_EPS 1e-10

#ifndef bool
#define bool int
#endif

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

int **RBINS = NULL;
real **X_COORDS = NULL;
real **C_ALPHA = NULL;

// grid resolution (used for fast clash detection)
#define GRID_RES 6.0

int chain_length = 0;

char AA_NAMES[21][4] =
  { "GLY", "ALA", "SER", "CYS", "VAL",
    "THR", "ILE", "PRO", "MET", "ASP",
    "ASN", "LEU", "LYS", "GLU", "GLN",
    "ARG", "HIS", "PHE", "TYR", "TRP",
    "UNK" };

char SHORT_AA_NAMES[22] = { "GASCVTIPMDNLKEQRHFYWX" };

int AA_NUMS[256];

int nheavy[20] = { 0, 1, 2, 2, 3, 3, 4, 3, 4, 4, 4, 4, 5, 5, 5, 7, 6, 7, 8, 10};

char *backbone_atoms[4] = { "N  ", "CA ", "C  ", "O  " };

char *heavy_atoms[200]= {
/* GLY */  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ALA */ "CB ", NULL,   NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* SER */ "CB ", "OG ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* CYS */ "CB ", "SG ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* VAL */ "CB ", "CG1", "CG2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* THR */ "CB ", "OG1", "CG2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ILE */ "CB ", "CG1", "CG2", "CD1",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* PRO */ "CB ", "CG ", "CD ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* MET */ "CB ", "CG ", "SD ", "CE ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ASP */ "CB ", "CG ", "OD1", "OD2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ASN */ "CB ", "CG ", "OD1", "ND2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* LEU */ "CB ", "CG ", "CD1", "CD2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* LYS */ "CB ", "CG ", "CD ", "CE ", "NZ ",  NULL,  NULL,  NULL,  NULL,  NULL,
/* GLU */ "CB ", "CG ", "CD ", "OE1", "OE2",  NULL,  NULL,  NULL,  NULL,  NULL,
/* GLN */ "CB ", "CG ", "CD ", "OE1", "NE2",  NULL,  NULL,  NULL,  NULL,  NULL,
/* ARG */ "CB ", "CG ", "CD ", "NE ", "CZ ", "NH1", "NH2",  NULL,  NULL,  NULL,
/* HIS */ "CB ", "CG ", "ND1", "CD2", "CE1", "NE2",  NULL,  NULL,  NULL,  NULL,
/* PHE */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ",  NULL,  NULL,  NULL,
/* TYR */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ", "OH ",  NULL,  NULL,
/* TRP */ "CB ", "CG ", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"};

/* reads full-atom pdb file */

struct _res_type;

typedef struct _atom_type {
  struct _atom_type *next;
  real x, y, z;
  char *name;
  int num, locnum;
  int flag;
  char cispro;
  int gx, gy, gz;
  struct _res_type *res;
  struct _atom_type *prev;
} atom_type;

typedef struct _res_type {
  struct _res_type *next;
  atom_type *atoms;
  int num, locnum, natoms;
  int type;
  char pdbsg;
  char protein;
  char *name;
  char chain;
  real sgx, sgy, sgz;
  real cmx, cmy, cmz;
  struct _res_type *prev;
} res_type;

typedef struct _mol_type {
  struct _mol_type *next;
  res_type *residua;
  int nres;
  unsigned char *r14;
  char *name;
  uchar *seq;
  char **contacts;
  real **cutoffs;
  struct _mol_type *prev;
} mol_type;

#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

mol_type *chain = NULL;


real rnd(void)
{
  return 0.001*(real)(rand()%1000);
}

// atom/res/mol manipulation functions

atom_type *new_atom(void)
{
  atom_type *tmpatom;

    tmpatom = (atom_type*) calloc(sizeof(atom_type),1);
    if (tmpatom) {
      tmpatom->x=tmpatom->y=tmpatom->z=0.;
      tmpatom->name=NULL;
      tmpatom->num=tmpatom->locnum=tmpatom->flag=0;
      tmpatom->next=tmpatom->prev=NULL;
    }

  return tmpatom;
}

res_type* new_res(void)
{
  res_type *tmpres;

    tmpres = (res_type*) calloc(sizeof(res_type),1);
    if (tmpres) {
      tmpres->num=0;
      tmpres->name=NULL;
      tmpres->atoms=NULL;
      tmpres->chain=' ';
      tmpres->next=tmpres->prev=NULL;
    }

  return tmpres;
}

mol_type *new_mol(void)
{
  mol_type *tmpmol;

    tmpmol = (mol_type*) calloc(sizeof(mol_type),1);
    if (tmpmol) {
      tmpmol->name=NULL;
      tmpmol->residua=NULL;
      tmpmol->next=tmpmol->prev=NULL;
    }

  return tmpmol;
}

void add_atom(atom_type* atomlist, atom_type* newatom)
{
  atom_type *tmpatom;

    if (!atomlist)
      atomlist=newatom;
    else {
      tmpatom=atomlist->next;
      atomlist->next=newatom;
      newatom->prev=atomlist;
      newatom->next=tmpatom;
      if (tmpatom) tmpatom->prev=newatom;
    }
}

void add_res(res_type* reslist, res_type* newres)
{
  res_type *tmpres;

    if (!reslist)
      reslist=newres;
    else {
      tmpres=reslist->next;
      reslist->next=newres;
      newres->prev=reslist;
      newres->next=tmpres;
      if (tmpres) tmpres->prev=newres;
    }
}

void add_mol(mol_type* mollist, mol_type* newmol)
{
  mol_type *tmpmol;

    if (!mollist)
      mollist=newmol;
    else {
      tmpmol=mollist->next;
      mollist->next=newmol;
      newmol->prev=mollist;
      newmol->next=tmpmol;
      if (tmpmol) tmpmol->prev=newmol;
    }
}

void delete_atom(atom_type* atom)
{
  atom_type *tmpatom;

    if (atom->prev) atom->prev->next=atom->next;
    if (atom->next) atom->next->prev=atom->prev;
    if (atom->name) free(atom->name);
    free(atom);
    atom=NULL;
}

void delete_res(res_type* res)
{
  res_type *tmpres;
  atom_type *tmpatom;

    if (res->prev) res->prev->next=res->next;
    if (res->next) res->next->prev=res->prev;
    if (res->name) free(res->name);
    if (res->atoms) {
      while (res->atoms) {
        tmpatom = res->atoms->next;
        delete_atom(res->atoms);
        res->atoms=tmpatom;
      }
    }
    free(res);
    res=NULL;
}

void delete_mol(mol_type* mol)
{
  mol_type *tmpmol;
  res_type *tmpres;
  int i;

    if (mol->prev) mol->prev->next=mol->next;
    if (mol->next) mol->next->prev=mol->prev;
    if (mol->name) free(mol->name);
    if (mol->residua) {
      while (mol->residua) {
        tmpres = mol->residua->next;
        delete_res(mol->residua);
        mol->residua=tmpres;
      }
    }
    if (mol->contacts) {
      for (i=0; i<mol->nres; i++) free(mol->contacts[i]);
      free(mol->contacts);
    }
    if (mol->cutoffs) {
      for (i=0; i<mol->nres; i++) free(mol->cutoffs[i]);
      free(mol->cutoffs);
    }
    free(mol);
    mol=NULL;
}


atom_type* get_last_atom(atom_type* atom)
{
    while (atom->next) atom=atom->next;

  return atom;
}

res_type* get_last_res(res_type* res)
{
    while (res->next) res=res->next;

  return res;
}

mol_type *get_last_mol(mol_type* mol)
{
    while (mol->next) mol=mol->next;

  return mol;
}

// single-aa from 3-letter code
char setseq(char* aaname)
{
  int i;

    for (i=0; i<21; i++) {
      if ((aaname[0]==AA_NAMES[i][0]) &&
          (aaname[1]==AA_NAMES[i][1]) &&
          (aaname[2]==AA_NAMES[i][2]))
         break;
    }
    
    if (i==21) {
      if (!strcmp(aaname, "GLX"))
        return 'E';
      if (!strcmp(aaname, "ASX"))
        return 'D';
      if (!strcmp(aaname, "HID"))
        return 'H';
      if (!strcmp(aaname, "MSE"))
        return 'M';
      if (!strcmp(aaname, "SEP"))
        return 'S';
      if (!strcmp(aaname, "TPO"))
        return 'T';
      if (!strcmp(aaname, "PTR"))
        return 'Y';
      i--;
    }

  return SHORT_AA_NAMES[i];
}

// side chain - side chain orientation
int orient(res_type *res1, res_type *res2)
{
  real x1, y1, z1;
  real x2, y2, z2;
  real cax, cay, caz;
  real len, vect, angle;
  atom_type *atom;

    if (!res1 || !res2) return 0;

    atom=res1->atoms;
    cax=cay=caz=0.;
    while (atom) {
      if (!strncmp(atom->name,"CA",2)) {
        cax=atom->x; cay=atom->y; caz=atom->z;
      }
      atom=atom->next;
    }
    x1=res1->sgx-cax; y1=res1->sgy-cay; z1=res1->sgz-caz;
    if (x1==0. && y1==0. && z1==0.) x1+=1.0;

    atom=res2->atoms;
    cax=cay=caz=0.;
    while (atom) {
      if (!strncmp(atom->name,"CA",2)) {
        cax=atom->x; cay=atom->y; caz=atom->z;
      }
      atom=atom->next;
    }
    x2=res2->sgx-cax; y2=res2->sgy-cay; z2=res2->sgz-caz;
    if (x2==0. && y2==0. && z2==0.) x2+=1.0;

    vect = x1*x2+y1*y2+z1*z2;
    len = sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2+z2*z2);
    if (len) vect /= len;

    angle=RADDEG*acos(vect);

    if (angle>120.) return 1; /*anti*/
    if (angle>60.) return 2;  /*mid*/

  return 3; /*par*/
}

int res_contact(res_type *res1, res_type *res2) {
  atom_type *atoms1, *atoms2;
  real dx, dy, dz;

    atoms1 = res1->atoms;
    while (atoms1) {
      atoms2 = res2->atoms;
      while (atoms2) {
        dx=atoms1->x-atoms2->x;
        dy=atoms1->y-atoms2->y;
        dz=atoms1->z-atoms2->z;
        if ((atoms1->flag & FLAG_SIDECHAIN) && (atoms2->flag & FLAG_SIDECHAIN) && (dx*dx+dy*dy+dz*dz<20.25)) {
          return 1;
        }
        atoms2=atoms2->next;
      }
      atoms1=atoms1->next;
    }

  return 0;
}


int read_pdb_file(char* filename, mol_type* molecules, char *realname)
{
  FILE *inp;
  char buffer[1000];
  char atmname[10];
  char resname[10];
  char version;
  int prevresnum, resnum, atmnum, locatmnum, num, locnum=0, i, j;
  atom_type *prevatom1, *prevatom2, *prevatom3, *prevatom4;
  int sgnum, cc, nres, ok, natom;
  real sgx, sgy, sgz;
  res_type *res, *test1, *test2;
  atom_type *atom;
  real x, y, z;
  real dist;
  unsigned char bin;
  int warn=0;
  real cutoff;

    if (_VERBOSE) printf("Reading input file %s...\n", filename);

    inp = fopen(filename, "r");
    if (!inp) {
      if (_VERBOSE) printf("ERROR: can't open %s !!!\n", filename);
      return FILE_NOT_FOUND;
    }

    molecules->nres=0;
    molecules->name=(char*)calloc(strlen(realname)+1,1);
    strcpy(molecules->name, realname);

    atmname[3]=0;
    resname[3]=0;
    prevresnum=-666;
    locatmnum=0;
    sgnum=0;
    sgx=sgy=sgz=0.;
    res=NULL;
    while (!feof(inp)) {
      if (fgets(buffer, 1000, inp)!=buffer) break;
      if (!strncmp(buffer, "END", 3) || !strncmp(buffer, "TER", 3)) break; // end of file; only single molecule read
      if (!strncmp(buffer, "ATOM", 4) || !strncmp(buffer, "HETATM", 6)) {
        if (buffer[16]!=' ' && buffer[16]!='A') continue;
        sscanf(&buffer[22], "%d", &resnum);
        strncpy(resname, &buffer[17], 3);
        strncpy(atmname, &buffer[13], 3);
        if (resnum==prevresnum && !strncmp(atmname, "N ", 2)) {
          if (_VERBOSE) printf("WARNING: fault in numeration at residuum %s[%d]\n", resname, resnum);
          warn=1;
        }
        if (atmname[0]=='H') continue;
        if (resnum!=prevresnum || !strncmp(atmname, "N ", 2)) {
          prevresnum=resnum;
          if (res)
            if (sgnum) {
              res->sgx=sgx/(real)sgnum;
              res->sgy=sgy/(real)sgnum;
              res->sgz=sgz/(real)sgnum;
            } else {
              res->sgx=res->sgy=res->sgz=0.;
            }
          locatmnum=0;
          version=' ';
          res = new_res();
          sgnum=0;
          sgx=sgy=sgz=0.;
          molecules->nres++;
          res->name = calloc(strlen(resname)+1, 1);
          res->type = AA_NUMS[setseq(resname)];
          res->locnum=locnum++;
          res->num = resnum;
          res->natoms=0;
          res->chain = buffer[21];
          strcpy(res->name, resname);
          if (molecules->residua) {
            add_res(get_last_res(molecules->residua), res);
          } else {
            molecules->residua = res;
          }
        }
        atom = new_atom();
        atom->res = res;
        atom->flag |= FLAG_INITIAL;
        res->natoms++;
        locatmnum++;
        sscanf(&buffer[7], "%d", &atmnum);
        sscanf(&buffer[30], "%lf", &x);
        sscanf(&buffer[38], "%lf", &y);
        sscanf(&buffer[46], "%lf", &z);
        version = buffer[16];
        atom->name = calloc(strlen(atmname)+1,1);
        strcpy(atom->name, atmname);
        atom->x=x; atom->y=y; atom->z=z;
        atom->num = atmnum;
        atom->locnum = locatmnum;
        if ((atmname[0]=='S' && atmname[1]=='C')||(atmname[0]=='C' && atmname[1]=='M')) {
          res->cmx = x;
          res->cmy = y;
          res->cmz = z;
          res->pdbsg=1;
          if (res->type<20) {
            res->protein=1;
          }
        } else
        if (!( ((atmname[0]=='C' || atmname[0]=='N' || atmname[0]=='O') && atmname[1]==' ') ||
               (atmname[0]=='H') ||
               (atmname[0]=='C' && atmname[1]=='A') ||
               (atmname[0]=='O' && atmname[1]=='X' && atmname[2]=='T') ) ) {
          sgx+=x;
          sgy+=y;
          sgz+=z;
          sgnum++;
          atom->flag |= FLAG_SIDECHAIN;
        } else
          atom->flag |= FLAG_BACKBONE;
        if (atmname[0]=='C' && atmname[1]=='A') {
          atom->flag |= FLAG_BACKBONE;
          if (res->type<20) {
            res->protein=1;
          }
          if (!res->pdbsg) {
            res->cmx = x;
            res->cmy = y;
            res->cmz = z;
          }
        }
        if (atmname[0]=='C' && atmname[1]=='M') {
          atom->flag |= FLAG_SCM;
        }
        if (atmname[0]=='S' && atmname[1]=='C') {
          atom->flag |= FLAG_SCM;
        }
        if (res->atoms) {
          add_atom(get_last_atom(res->atoms), atom);
        } else {
          res->atoms = atom;
        }
      }
    }

    if (res)
      if (sgnum) {
        res->sgx=sgx/(real)sgnum;
        res->sgy=sgy/(real)sgnum;
        res->sgz=sgz/(real)sgnum;
      } else {
        res->sgx=res->sgy=res->sgz=0.;
      }

    fclose(inp);

    molecules->seq = (uchar*)calloc(sizeof(uchar)*molecules->nres+1,1);
    res=molecules->residua; i=0;
    while (res) {
      molecules->seq[i++]=(uchar)AA_NUMS[(int)setseq(res->name)];
      res=res->next;
    }

  if (!warn) return FILE_SUCCESS; else return FILE_WARNING;
}

// energy calculation for C-alpha optimizer
real calc_ca_energy(atom_type **c_alpha, real **new_c_alpha, real **init_c_alpha, real **gradient, real alpha, real *ene, bool calc_gradient)
{
  int i, j;
  real dx, dy, dz;
  real dist, ddist, ddist2;
  real new_e_pot;
  real theta0, tdif, th, aa, bb, ab;
  real ff0, ff2, dth, m0, m2, grad, f0[3], f2[3];
  real adiff[3], bdiff[3];
  real deriv, theta, dtheta, len1, len2, cos_theta, sin_theta;
  real dx1, dy1, dz1;
  real dx2, dy2, dz2;
  real dx3, dy3, dz3;
  real vx1, vy1, vz1;
  real vx2, vy2, vz2;
  real vx3, vy3, vz3;

  real r12x, r12y, r12z;
  real r32x, r32y, r32z;
  real d12, d32, d12inv, d32inv, c1, c2, diff;
  real f1x, f1y, f1z;
  real f2x, f2y, f2z;
  real f3x, f3y, f3z;

        for (i=0; i<chain_length; i++) {
          new_c_alpha[i][0]=c_alpha[i]->x+alpha*gradient[i][0];
          new_c_alpha[i][1]=c_alpha[i]->y+alpha*gradient[i][1];
          new_c_alpha[i][2]=c_alpha[i]->z+alpha*gradient[i][2];
        }

        new_e_pot = 0.0;

        ene[0]=ene[1]=ene[2]=ene[3]=0.0;

        for (i=0; i<chain_length; i++) {
#ifdef CALC_C_ALPHA_START
          dx=new_c_alpha[i][0]-init_c_alpha[i][0];
          dy=new_c_alpha[i][1]-init_c_alpha[i][1];
          dz=new_c_alpha[i][2]-init_c_alpha[i][2];
          dist=sqrt(dx*dx+dy*dy+dz*dz);
          ddist = -dist;
          if (dist>_CA_START_DIST) {
            ddist2=dist*dist;
            new_e_pot+=CA_START_K*ddist2;
            ene[1] += CA_START_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_START_K)/dist;
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
            }
          }

#endif


#ifdef CALC_C_ALPHA
          if (i>0) {
            dx=new_c_alpha[i][0]-new_c_alpha[i-1][0];
            dy=new_c_alpha[i][1]-new_c_alpha[i-1][1];
            dz=new_c_alpha[i][2]-new_c_alpha[i-1][2];
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (c_alpha[i]->cispro) {
              ddist=CA_DIST_CISPRO-dist;
//              if (fabs(ddist)<CA_DIST_CISPRO_TOL) ddist=0.0;
            } else {
              ddist=CA_DIST-dist;
//              if (fabs(ddist)<CA_DIST_TOL) ddist=0.0;
            }
            ddist2=ddist*ddist;
            new_e_pot+=CA_K*ddist2;
            ene[0] += CA_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_K)/dist;
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
              gradient[i-1][0]+=grad*dx;
              gradient[i-1][1]+=grad*dy;
              gradient[i-1][2]+=grad*dz;
            }
          }
#endif

#ifdef CALC_C_ALPHA_XVOL
          for (j=0;j<i;j++) {
            if (abs(i-j)>2) {
              dx=new_c_alpha[i][0]-new_c_alpha[j][0];
              dy=new_c_alpha[i][1]-new_c_alpha[j][1];
              dz=new_c_alpha[i][2]-new_c_alpha[j][2];
              dist=sqrt(dx*dx+dy*dy+dz*dz);
              ddist = dist-_CA_XVOL_DIST;
              if (dist<_CA_XVOL_DIST) {
                ddist2 = dist*dist;
                new_e_pot+=CA_XVOL_K*ddist2;
                ene[3] += CA_XVOL_K*ddist2;
                if (calc_gradient) {
                  grad = ddist*(8.0*CA_XVOL_K)/dist;
                  gradient[i][0]-=grad*dx;
                  gradient[i][1]-=grad*dy;
                  gradient[i][2]-=grad*dz;
                  gradient[j][0]+=grad*dx;
                  gradient[j][1]+=grad*dy;
                  gradient[j][2]+=grad*dz;
                }
              }
            }
          }
#endif

#ifdef CALC_C_ALPHA_ANGLES

        if (i>0 && i<chain_length-1) {
          r12x=new_c_alpha[i-1][0]-new_c_alpha[i][0];
          r12y=new_c_alpha[i-1][1]-new_c_alpha[i][1];
          r12z=new_c_alpha[i-1][2]-new_c_alpha[i][2];
          r32x=new_c_alpha[i+1][0]-new_c_alpha[i][0];
          r32y=new_c_alpha[i+1][1]-new_c_alpha[i][1];
          r32z=new_c_alpha[i+1][2]-new_c_alpha[i][2];
          d12 = sqrt(r12x*r12x+r12y*r12y+r12z*r12z);
          d32 = sqrt(r32x*r32x+r32y*r32y+r32z*r32z);
          cos_theta = (r12x*r32x+r12y*r32y+r12z*r32z)/(d12*d32);
          if (cos_theta>1.0)
            cos_theta = 1.0;
          else
          if (cos_theta<-1.0)
            cos_theta = -1.0;
          sin_theta = sqrt(1.0-cos_theta*cos_theta);
          theta = acos(cos_theta);

          if (RADDEG*theta<80.)
            diff = theta-80.*DEGRAD;
          else
          if (RADDEG*theta>150.)
            diff = theta-150.*DEGRAD;
          else
            diff = 0.0;

          new_e_pot += CA_ANGLE_K*diff*diff;
          ene[2] += CA_ANGLE_K*diff*diff;
          if (calc_gradient) {
            d12inv = 1.0/d12;
            d32inv = 1.0/d32;
            diff *= (-2.0*CA_ANGLE_K)/sin_theta;
            c1 = diff*d12inv;
            c2 = diff*d32inv;
            f1x = c1*(r12x*(d12inv*cos_theta)-r32x*d32inv);
            f1y = c1*(r12y*(d12inv*cos_theta)-r32y*d32inv);
            f1z = c1*(r12z*(d12inv*cos_theta)-r32z*d32inv);
            f3x = c2*(r32x*(d32inv*cos_theta)-r12x*d12inv);
            f3y = c2*(r32y*(d32inv*cos_theta)-r12y*d12inv);
            f3z = c2*(r32z*(d32inv*cos_theta)-r12z*d12inv);
            f2x = -f1x-f3x;
            f2y = -f1y-f3y;
            f2z = -f1z-f3z;
            gradient[i-1][0]+=f1x;
            gradient[i-1][1]+=f1y;
            gradient[i-1][2]+=f1z;
            gradient[i][0]+=f2x;
            gradient[i][1]+=f2y;
            gradient[i][2]+=f2z;
            gradient[i+1][0]+=f3x;
            gradient[i+1][1]+=f3y;
            gradient[i+1][2]+=f3z;
          }
        }

#endif

        }

//printf("ene[3] = %f\n", ene[3]);

  return new_e_pot;
}

/*
 *  Steepest gradient optimization using v=k*(r0-r)^2
 *  k = CA_K, r0 = CA_DIST
 */
void ca_optimize(char *tname, char *iname)
{
  char buf[1000];
  int i, j, hx, my_iter;
  real dx, dy, dz, dd, dist, dist2, dist3, ddist, ddist2;
  real e_pot, new_e_pot, grad, alpha, e_pot1, e_pot2, e_pot3;
  real adiff[3], bdiff[3];
  real ff0, ff2, aa, ab, bb, th, tdif, dth, m0, m2;
  real theta0, deg_th, maxgrad, sum;
  real f0[3], f2[3];
  real x, y, z;
  int numsteps, numsteps2, msteps;
  int *sec;
  real **new_c_alpha, **gradient, **init_c_alpha, last_alpha, tmp, last_good_alpha, d_alpha, last_e_pot;
  atom_type *atom, **c_alpha;
  res_type *res;
  FILE *inp, *out;
  int mnum, init, ok;
  real alpha1, alpha2, alpha3, a0;
  real ene1, ene2, ene3, e0;
  real energies[4];
  real w1, w2, w3, eps;
  real gnorm, last_gnorm;
  int mode, fcnt;


    if (_CA_TRAJECTORY) {
      out = fopen(tname,"w");
      if (out) fclose(out);
    }

    if (_VERBOSE) printf("Alpha carbons optimization...\n");
      
    new_c_alpha = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    init_c_alpha = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    for (i=0;i<=chain_length;i++) {
      new_c_alpha[i] = (real*)calloc(sizeof(real)*3,1);
      init_c_alpha[i] = (real*)calloc(sizeof(real)*3,1);
    }
    gradient = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    for (i=0;i<=chain_length;i++) {
      gradient[i] = (real*)calloc(sizeof(real)*3,1);
    }

    c_alpha = (atom_type**)calloc(sizeof(atom_type*)*(chain_length+1),1);

    i = 0;
    res = chain->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A') {
          if (i<chain_length) {
            c_alpha[i] = atom;
            i++;
            break;
          } else {
            if (_VERBOSE) printf("WARNING: number of C-alpha atoms exceeds the chain length!\n");
            break;
          }
        }
        atom = atom->next;
      }
      res = res->next;
    }

    if (i<chain_length) chain_length = i;
      
    for (i=0; i<chain_length; i++) {
      init_c_alpha[i][0] = c_alpha[i]->x;
      init_c_alpha[i][1] = c_alpha[i]->y;
      init_c_alpha[i][2] = c_alpha[i]->z;
    }

    if (_CISPRO) {
      for (i=1; i<chain_length; i++) {
        dx = c_alpha[i]->x-c_alpha[i-1]->x;
        dy = c_alpha[i]->y-c_alpha[i-1]->y;
        dz = c_alpha[i]->z-c_alpha[i-1]->z;
        dd = sqrt(dx*dx+dy*dy+dz*dz);
        if ((setseq(c_alpha[i]->res->name)=='P') && (dd>CA_DIST_CISPRO-5*CA_DIST_CISPRO_TOL) && (dd<CA_DIST_CISPRO+5*CA_DIST_CISPRO_TOL)) {
          if (_VERBOSE) printf("Probable cis-proline found at postion %d\n", c_alpha[i]->res->num);
          c_alpha[i]->cispro = 1;
        }
      }
    }

    if (_CA_RANDOM) {
      if (_VERBOSE) printf("Generating random C-alpha coordinates...\n");
      c_alpha[0]->x = 0.0;
      c_alpha[0]->y = 0.0;
      c_alpha[0]->z = 0.0;
      for (i=1;i<chain_length;i++) {
        dx = 0.01*(100-rand()%200);
        dy = 0.01*(100-rand()%200);
        dz = 0.01*(100-rand()%200);
        dd = 3.8/sqrt(dx*dx+dy*dy+dz*dz);
        dx *= dd;
        dy *= dd;
        dz *= dd;
        c_alpha[i]->x = c_alpha[i-1]->x+dx;
        c_alpha[i]->y = c_alpha[i-1]->y+dy;
        c_alpha[i]->z = c_alpha[i-1]->z+dz;
      }
    }

    if (iname) {
      inp = fopen(iname,"r");
      if (inp) {
        if (_VERBOSE) printf("Reading initial structure %s...\n", iname);
        i = 0;
        while (!feof(inp)) {
          if (fgets(buf,1000,inp)==buf && buf[13]=='C' && buf[14]=='A') {
            if (i<chain_length) {
              if (sscanf(&buf[30],"%lf%lf%lf",&x,&y,&z)==3) {
                c_alpha[i]->x = x;
                c_alpha[i]->y = y;
                c_alpha[i]->z = z;
                i++;
              }
            } else {
              if (_VERBOSE) printf("WARNING: number of ini-file C-alpha atoms exceeds the chain length!\n");
              break;
            }
          }
        }
        fclose(inp);
      } else
        if (_VERBOSE) printf("WARNING: can't read initial corrdinates %s\n", iname);
    }

    mnum = 1;
    mode = 0;
    init = 0;
    numsteps=numsteps2=0;
    last_alpha = 0.0;


    if (_VERBOSE) printf("Optimizing alpha carbons...\n");

    eps = 0.5;

    fcnt=0;

    last_gnorm = 1000.;

    do {
      last_e_pot = e_pot;

      if (_CA_TRAJECTORY) {
        out = fopen(tname,"a");
        if (out) {
          fprintf(out,"MODEL  %d\n",mnum++);
          for (i=0; i<chain_length; i++) {
            fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                    i+1, "CA ", c_alpha[i]->res->name, ' ', c_alpha[i]->res->num,
  				          c_alpha[i]->x, c_alpha[i]->y, c_alpha[i]->z);

          }
          fprintf(out,"ENDMDL\n");
          fclose(out);
        }
      }

// calculate gradients

      e_pot=e_pot1=e_pot2=e_pot3=0.;

      for (i=0; i<chain_length; i++)
        gradient[i][0]=gradient[i][1]=gradient[i][2]=0.;

      e_pot = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, 0.0, energies, true);

      if (_VERBOSE && !init) {
        printf("Initial energy: bond=%.5lf angle=%.5f restraints=%.5f xvol=%.5f total=%.5f\n", energies[0], energies[2], energies[1], energies[3], e_pot);
      }

      if (!init) init=1;

// LINE SEARCH

      alpha1 = -1.0;
      alpha2 = 0.0;
      alpha3 = 1.0;

      ene1 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
      ene2 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha2, energies, false);
      ene3 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);

      msteps = 0;
      while (ene2>MIN(ene1,ene3) && msteps<_CA_ITER) {
        msteps++;
        alpha1 *= 2.0;
        alpha3 *= 2.0;
        ene1 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
        ene3 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);
      }

      msteps = 0;
      do {
        if (alpha3-alpha2>alpha2-alpha1) {
          a0 = 0.5*(alpha2+alpha3);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0-1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0+1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          if (e0<ene2) {
            alpha1 = alpha2;
            alpha2 = a0;
            ene1 = ene2;
            ene2 = e0;
          } else {
            alpha3 = a0;
            ene3 = e0;
          }
        } else {
          a0 = 0.5*(alpha1+alpha2);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0-1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0+1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          if (e0<ene2) {
            alpha3 = alpha2;
            alpha2 = a0;
            ene3 = ene2;
            ene2 = e0;
          } else {
            alpha1 = a0;
            ene1 = e0;
          }
        }
        msteps++;
      } while (alpha3-alpha1>1e-6 && msteps<20);

      last_alpha = alpha2;
      e_pot = ene2;

      for (i=0; i<chain_length; i++) {
        c_alpha[i]->x=c_alpha[i]->x+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][0];
        c_alpha[i]->y=c_alpha[i]->y+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][1];
        c_alpha[i]->z=c_alpha[i]->z+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][2];
      }

      e_pot = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, 0.0, energies, false);

      eps *= 0.75;
      if (eps<1e-3) eps=0.0;

      numsteps++;

      gnorm = 0.0;

	  for (i=0; i<chain_length; i++) {
        gnorm += gradient[i][0]*gradient[i][0] + gradient[i][1]*gradient[i][1] + gradient[i][2]*gradient[i][2];
      }

      gnorm = sqrt(gnorm/(double)chain_length);

      if (last_gnorm-gnorm<1e-3) fcnt++;

      last_gnorm = gnorm;

    } while ( (fcnt<3) &&  (gnorm>0.01) && (numsteps<_CA_ITER));


     if (_VERBOSE) {
        for (i=0; i<chain_length; i++) {

#ifdef CALC_C_ALPHA
          if (i>0) {
            dx=c_alpha[i]->x-c_alpha[i-1]->x;
            dy=c_alpha[i]->y-c_alpha[i-1]->y;
            dz=c_alpha[i]->z-c_alpha[i-1]->z;
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (c_alpha[i]->cispro) {
              ddist=CA_DIST_CISPRO-dist;
              if (fabs(ddist)<CA_DIST_CISPRO_TOL) ddist=0.0;
            } else {
              ddist=CA_DIST-dist;
              if (fabs(ddist)<CA_DIST_TOL) ddist=0.0;
            }
            ddist2=ddist*ddist;
       	    if (fabs(ddist)>=CA_DIST_TOL) printf("WARNING: distance %d = %.3lf A\n", i, dist);
          }
#endif
        }

        for (i=0; i<chain_length; i++) {
#ifdef CALC_C_ALPHA_ANGLES
          if (i>0 && i<chain_length-1) {
            aa=ab=bb=0.0;
            adiff[0]=c_alpha[i-1]->x-c_alpha[i]->x;
            bdiff[0]=c_alpha[i+1]->x-c_alpha[i]->x;
            aa+=adiff[0]*adiff[0];
            ab+=adiff[0]*bdiff[0];
            bb+=bdiff[0]*bdiff[0];
            adiff[1]=c_alpha[i-1]->y-c_alpha[i]->y;
            bdiff[1]=c_alpha[i+1]->y-c_alpha[i]->y;
            aa+=adiff[1]*adiff[1];
            ab+=adiff[1]*bdiff[1];
            bb+=bdiff[1]*bdiff[1];
            adiff[2]=c_alpha[i-1]->z-c_alpha[i]->z;
            bdiff[2]=c_alpha[i+1]->z-c_alpha[i]->z;
            aa+=adiff[2]*adiff[2];
            ab+=adiff[2]*bdiff[2];
            bb+=bdiff[2]*bdiff[2];

            th=ab/sqrt(aa*bb);
            if (th<-1.0) th=-1.0;
            if (th>1.0) th=1.0;
            th=acos(th);
            deg_th=RADDEG*th;
            if (deg_th>150.) theta0=DEGRAD*150.; else
            if (deg_th<75.) theta0=DEGRAD*75.; else
            theta0=th;
       	    if (fabs(deg_th-RADDEG*theta0)>1.0) printf("WARNING: angle %d = %.3lf degrees\n", i, deg_th);
          }
#endif
        }
      }

    if (_VERBOSE) printf("Optimization done after %d step(s).\nFinal energy: bond=%.5lf angle=%.5f restraints=%.5f xvol=%.5f total=%.5f\n", numsteps, energies[0], energies[2], energies[1], energies[3], e_pot);

    if (_CA_TRAJECTORY) {
      out = fopen(tname,"a");
      if (out) {
        fprintf(out,"END\n");
      }
    }

    for (i=0;i<chain_length+1;i++) {
      free(init_c_alpha[i]);
      free(new_c_alpha[i]);
      free(gradient[i]);
    }
    free(new_c_alpha);
    free(gradient);
    free(c_alpha);
    free(init_c_alpha);
}

void center_chain(mol_type *mol)
{
  real cx, cy, cz;
  int natom;
  res_type *res;
  atom_type *atom;

    cx = cy = cz = 0.0;
    natom = 0;

    res = mol->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        cx += atom->x;
        cy += atom->y;
        cz += atom->z;
        natom++;
  			atom=atom->next;
  		}
      res = res->next;
    }

    cx /= (real)natom;
    cy /= (real)natom;
    cz /= (real)natom;

    if (_VERBOSE) printf("Molecule center: %8.3f %8.3f %8.3f -> 0.000 0.000 0.000\n", cx, cy, cz);

    res = mol->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (!(_PRESERVE && (atom->flag & FLAG_INITIAL))) {        
          atom->x -= cx;
          atom->y -= cy;
          atom->z -= cz;
        }
        natom++;
  		atom=atom->next;
  	  }
      res = res->next;
    }

}

void write_pdb(char *name, mol_type *mol)
{
  FILE *out;
  res_type *res;
  atom_type *atom, *oxt;
  int anum;

    oxt = NULL;
    out = fopen(name,"w");
    if (!out) {
      if (_VERBOSE) printf("Can't open output file!\n");
      return;
    }
    fprintf(out,"REMARK 999 REBUILT BY PULCHRA V.%.2f\n", PULCHRA_VERSION);
    anum=1;
    res = mol->residua;
    while (res) {
      if (res->protein) {
        if (!_BB_REARRANGE) {
          atom = res->atoms;
          while (atom) {
            if (!(atom->name[0]=='D' && atom->name[1]=='U') &&
                !(atom->name[0]=='S' && atom->name[1]=='C') &&
                !(atom->name[0]=='C' && atom->name[1]=='M') &&
                !(atom->name[0]=='O' && atom->name[1]=='X') &&
                !(atom->name[0]=='H' && !_REBUILD_H))
              fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                            anum++, atom->name, res->name, ' ', res->num,
    	    				          atom->x, atom->y, atom->z);
            if (atom->name[0]=='O' && atom->name[1]=='X') oxt = atom;      	  
      	    atom=atom->next;
	  }
      	} else {
          atom = res->atoms;
          while (atom) {
            if (!(atom->name[0]=='D' && atom->name[1]=='U') &&
                !(atom->name[0]=='S' && atom->name[1]=='C') &&
                !(atom->name[0]=='C' && atom->name[1]==' ') &&
                !(atom->name[0]=='O' && atom->name[1]==' ') &&
                !(atom->name[0]=='C' && atom->name[1]=='M') &&
                !(atom->name[0]=='O' && atom->name[1]=='X') &&
                !(atom->name[0]=='H' && !_REBUILD_H))
              fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                            anum++, atom->name, res->name, ' ', res->num,
    	    				          atom->x, atom->y, atom->z);
            if (atom->name[0]=='O' && atom->name[1]=='X') oxt = atom;      	  
      	    atom=atom->next;
      	  }
          atom = res->atoms;
          while (atom) {
            if (((atom->name[0]=='C' && atom->name[1]==' ') ||
                (atom->name[0]=='O' && atom->name[1]==' ')) &&
               !(atom->name[0]=='H' && !_REBUILD_H))
              fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                            anum++, atom->name, res->name, ' ', res->num,
    	    				          atom->x, atom->y, atom->z);
      			atom=atom->next;
      		}
      	}
      }
      if (!res->next && oxt) {
        atom = oxt;
        fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                     anum++, atom->name, res->name, ' ', res->num, atom->x, atom->y, atom->z);
      }
      res = res->next;
    }
    fprintf(out,"TER\nEND\n");
    fclose(out);
}


void write_pdb_sg(char *name, mol_type *mol)
{
  FILE *out;
  res_type *res;
  atom_type *atom;
  int anum;

    out = fopen(name,"w");
    if (!out) {
      if (_VERBOSE) printf("Can't open output file!\n");
      return;
    }
    fprintf(out,"REMARK 999 REBUILT BY PULCHRA V.%.2f\n", PULCHRA_VERSION);
    anum=1;
    res = mol->residua;
    while (res) {
      if (res->protein) {
        atom = res->atoms;
        while (atom) {
          if ((atom->name[0]=='C' && atom->name[1]=='A'))
            fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                          anum++, atom->name, res->name, ' ', res->num,
    	  				          atom->x, atom->y, atom->z);
    			atom=atom->next;
    		}
        fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                      anum++, "CM ", res->name, ' ', res->num,
    				          res->cmx, res->cmy, res->cmz);
    	}
      res = res->next;
    }
    fprintf(out,"TER\nEND\n");
    fclose(out);
}

// distance
real calc_distance(real x1, real y1, real z1,
							 		  real x2, real y2, real z2)
{
  real dx,dy,dz;
  real dist2;

    dx = (x1) - (x2);
    dy = (y1) - (y2);
    dz = (z1) - (z2);
    if (dx || dy || dz ) {
      dist2 = dx*dx+dy*dy+dz*dz;
      return (sqrt(dist2));
    } else
      return 0.0;
}

// r14 chiral distance
real calc_r14(real x1, real y1, real z1,
							 real x2, real y2, real z2,
							 real x3, real y3, real z3,
							 real x4, real y4, real z4)
{
  real r, dx, dy, dz;
  real vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3;
  real hand;

    dx = x4-x1;
    dy = y4-y1;
    dz = z4-z1;

    r = sqrt(dx*dx+dy*dy+dz*dz);

    vx1=x2-x1;
    vy1=y2-y1;
    vz1=z2-z1;
    vx2=x3-x2;
    vy2=y3-y2;
    vz2=z3-z2;
    vx3=x4-x3;
    vy3=y4-y3;
    vz3=z4-z3;

    hand = (vy1*vz2-vy2*vz1)*vx3+
           (vz1*vx2-vz2*vx1)*vy3+
           (vx1*vy2-vx2*vy1)*vz3;

    if (hand<0) r=-r;

  return r;
}

// superimposition of two sets for coordinates + optional transformation of tpoints

real superimpose2(real **coords1, real **coords2, int npoints, real **tpoints, int ntpoints)
{
  real mat_s[3][3], mat_a[3][3], mat_b[3][3], mat_g[3][3];
  real mat_u[3][3], tmp_mat[3][3];
  real val, d, alpha, beta, gamma, x, y, z;
  real cx1, cy1, cz1, cx2, cy2, cz2, tmpx, tmpy, tmpz;
  int i, j, k, n;

    cx1=cy1=cz1=cx2=cy2=cz2=0.;

    for (i=0; i<npoints; i++) {
      cx1+=coords1[i][0];
      cy1+=coords1[i][1];
      cz1+=coords1[i][2];
      cx2+=coords2[i][0];
      cy2+=coords2[i][1];
      cz2+=coords2[i][2];
    }

    cx1/=(real)npoints;
    cy1/=(real)npoints;
    cz1/=(real)npoints;

    cx2/=(real)npoints;
    cy2/=(real)npoints;
    cz2/=(real)npoints;

    for (i=0; i<npoints; i++) {
      coords1[i][0]-=cx1;
      coords1[i][1]-=cy1;
      coords1[i][2]-=cz1;
      coords2[i][0]-=cx2;
      coords2[i][1]-=cy2;
      coords2[i][2]-=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]-=cx2;
      tpoints[i][1]-=cy2;
      tpoints[i][2]-=cz2;
    }

    for (i=0; i<3; i++)
      for (j=0; j<3; j++) {
        if (i==j)
          mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=1.0;
        else
          mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=0.0;
        mat_u[i][j]=0.;
      }

    for (n=0; n<npoints; n++) {
      mat_u[0][0]+=coords1[n][0]*coords2[n][0];
      mat_u[0][1]+=coords1[n][0]*coords2[n][1];
      mat_u[0][2]+=coords1[n][0]*coords2[n][2];
      mat_u[1][0]+=coords1[n][1]*coords2[n][0];
      mat_u[1][1]+=coords1[n][1]*coords2[n][1];
      mat_u[1][2]+=coords1[n][1]*coords2[n][2];
      mat_u[2][0]+=coords1[n][2]*coords2[n][0];
      mat_u[2][1]+=coords1[n][2]*coords2[n][1];
      mat_u[2][2]+=coords1[n][2]*coords2[n][2];
    }

    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        tmp_mat[i][j]=0.;

    do {
      d=mat_u[2][1]-mat_u[1][2];
      if (d==0) alpha=0; else alpha=atan(d/(mat_u[1][1]+mat_u[2][2]));
      if (cos(alpha)*(mat_u[1][1]+mat_u[2][2])+sin(alpha)*(mat_u[2][1]-mat_u[1][2])<0.0)       alpha+=M_PI;
      mat_a[1][1]=mat_a[2][2]=cos(alpha);
      mat_a[2][1]=sin(alpha);
      mat_a[1][2]=-mat_a[2][1];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_a[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_a[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[0][2]-mat_u[2][0];
      if (d==0) beta=0; else beta=atan(d/(mat_u[0][0]+mat_u[2][2]));
      if (cos(beta)*(mat_u[0][0]+mat_u[2][2])+sin(beta)*(mat_u[0][2]-mat_u[2][0])<0.0) beta+=M_PI;
      mat_b[0][0]=mat_b[2][2]=cos(beta);
      mat_b[0][2]=sin(beta);
      mat_b[2][0]=-mat_b[0][2];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_b[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_b[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[1][0]-mat_u[0][1];
      if (d==0) gamma=0; else gamma=atan(d/(mat_u[0][0]+mat_u[1][1]));
      if (cos(gamma)*(mat_u[0][0]+mat_u[1][1])+sin(gamma)*(mat_u[1][0]-mat_u[0][1])<0.0)
        gamma+=M_PI;
      mat_g[0][0]=mat_g[1][1]=cos(gamma);
      mat_g[1][0]=sin(gamma);
      mat_g[0][1]=-mat_g[1][0];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_g[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_g[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      val=fabs(alpha)+fabs(beta)+fabs(gamma);
    } while (val>0.001);

    val=0.;
    for (i=0; i<npoints; i++) {
      x=coords2[i][0];
      y=coords2[i][1];
      z=coords2[i][2];
      tmpx=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tmpy=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tmpz=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
      x=coords1[i][0]-tmpx;
      y=coords1[i][1]-tmpy;
      z=coords1[i][2]-tmpz;
      val+=x*x+y*y+z*z;
    }

    for (i=0; i<ntpoints; i++) {
      x=tpoints[i][0];
      y=tpoints[i][1];
      z=tpoints[i][2];
      tpoints[i][0]=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tpoints[i][1]=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tpoints[i][2]=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
    }

    for (i=0; i<npoints; i++) {
      coords1[i][0]+=cx1;
      coords1[i][1]+=cy1;
      coords1[i][2]+=cz1;
      coords2[i][0]+=cx2;
      coords2[i][1]+=cy2;
      coords2[i][2]+=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]+=cx1;
      tpoints[i][1]+=cy1;
      tpoints[i][2]+=cz1;
    }

  return sqrt(val/(real)npoints);
}


atom_type *find_atom(res_type *res, char *aname)
{
  atom_type *atom;

    atom = res->atoms;
    while (atom) {
      if (atom->name[0]==aname[0] && atom->name[1]==aname[1] && atom->name[2]==aname[2]) {
        return atom;
        break;
      }
      atom = atom->next;
    }

  return NULL;
}


void add_replace(res_type *res, char *aname, real x, real y, real z, int flags)
{
  atom_type *atom, *newatom;

    atom = res->atoms;
    while (atom) {
      if (atom->name[0]==aname[0] && atom->name[1]==aname[1] && atom->name[2]==aname[2]) {
        if (!(_PRESERVE && (atom->flag & FLAG_INITIAL))) { atom->x = x; atom->y = y; atom->z = z; }
        atom->flag |= flags;
        break;
      }
      atom = atom->next;
    }

    if (!atom) {
      newatom = (atom_type*)calloc(sizeof(atom_type),1);
      newatom->x = x;
      newatom->y = y;
      newatom->z = z;
      newatom->flag |= flags;
      newatom->res = res;
      newatom->name = (char*)calloc(4,1);
      strcpy(newatom->name,aname);

      atom = res->atoms;
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A')
          break;
        atom = atom->next;
      }
      if (aname[0]=='N' && aname[1]==' ') {
        newatom->next = res->atoms;
        res->atoms = newatom;
      } else {
        while (atom->next) atom=atom->next;
        atom->next = newatom;
      }
    }
}



void prepare_rbins(void)
{
  int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real r13_1, r13_2, r14;
  real **cacoords, **tmpcoords, **tmpstat;
  res_type *res, *prevres;
  atom_type *atom;

  if (!RBINS) {
    RBINS = (int**)calloc(sizeof(int*)*(chain_length+1),1);
    for (i=0;i<chain_length+1;i++)
      RBINS[i] = (int*)calloc(sizeof(int)*3,1);

    X_COORDS = (real**)calloc(sizeof(real*)*(chain_length+10),1);
    for (i=0;i<chain_length+10;i++)
      X_COORDS[i] = (real*)calloc(sizeof(real)*3,1);

    i = 5;

    res = chain->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A') {
          X_COORDS[i][0] = atom->x;
          X_COORDS[i][1] = atom->y;
          X_COORDS[i][2] = atom->z;
          i++;
        }
        atom = atom->next;
      }
      res = res->next;
    }

    C_ALPHA = &X_COORDS[5];

    cacoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpcoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpstat = (real**)calloc(sizeof(real*)*(8),1);
    for (i=0;i<8;i++) {
      cacoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpcoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpstat[i] = (real*)calloc(sizeof(real)*3,1);;
    }

    // rebuild ends...

    for (i=0,j=0;i<5;i++,j++)
      for (k=0;k<3;k++)
        tmpcoords[j][k] = C_ALPHA[i][k];
    for (i=2,j=0;i<5;i++,j++)
      for (k=0;k<3;k++)
        cacoords[j][k] = C_ALPHA[i][k];
    for (i=0,j=0;i<3;i++,j++)
      for (k=0;k<3;k++)
        tmpstat[j][k] = C_ALPHA[i][k];

    superimpose2(tmpstat,cacoords,3,tmpcoords,5);

    for (i=-2,j=0;i<0;i++,j++)
      for (k=0;k<3;k++)
        C_ALPHA[i][k] = tmpcoords[j][k];

    for (i=chain_length-5,j=0;i<chain_length;i++,j++)
      for (k=0;k<3;k++)
        tmpcoords[j][k] = C_ALPHA[i][k];
    for (i=chain_length-5,j=0;i<chain_length-2;i++,j++)
      for (k=0;k<3;k++)
        cacoords[j][k] = C_ALPHA[i][k];
    for (i=chain_length-3,j=0;i<chain_length;i++,j++)
      for (k=0;k<3;k++)
        tmpstat[j][k] = C_ALPHA[i][k];

    superimpose2(tmpstat,cacoords,3,tmpcoords,5);

    for (i=chain_length-3,j=0;i<chain_length;i++,j++)
      for (k=0;k<3;k++)
        C_ALPHA[i+3][k] = tmpcoords[j+3][k];

    for (i=0;i<chain_length+1;i++) {
    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];

    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];

    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];

    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];

    	r13_1 = calc_distance(x1, y1, z1, x3, y3, z3);
    	r13_2 = calc_distance(x2, y2, z2, x4, y4, z4);
    	r14 = calc_r14(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

      bin13_1 = (int)((r13_1-4.6)/0.3);
      bin13_2 = (int)((r13_2-4.6)/0.3);
      bin14 = (int)((r14+11.)/0.3);

      if (bin13_1<0) bin13_1=0;
      if (bin13_2<0) bin13_2=0;
      if (bin14<0) bin14=0;
      if (bin13_1>9) bin13_1=9;
      if (bin13_2>9) bin13_2=9;
      if (bin14>73) bin14=73;

      RBINS[i][0] = bin13_1;
      RBINS[i][1] = bin13_2;
      RBINS[i][2] = bin14;
    }
  }
}

#ifdef COMPILE_BB

void rebuild_backbone(void)
{

  res_type *res, *prevres;
  atom_type *atom;
  real **cacoords, **tmpcoords, **tmpstat;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real besthit, hit;
  int bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real rmsd, total, maxrms;
  FILE *debug, *out;

    if (_VERBOSE) printf("Rebuilding backbone...\n");

    prepare_rbins();

    cacoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpcoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpstat = (real**)calloc(sizeof(real*)*(8),1);
    for (i=0;i<8;i++) {
      cacoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpcoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpstat[i] = (real*)calloc(sizeof(real)*3,1);;
    }


    prevres = NULL;
    res = chain->residua;


    total = maxrms = 0.0;

    for (i=0;i<chain_length+1;i++) {
    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];

    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];

    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];

    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];

    	cacoords[0][0] = x1;
    	cacoords[0][1] = y1;
    	cacoords[0][2] = z1;

    	cacoords[1][0] = x2;
     	cacoords[1][1] = y2;
     	cacoords[1][2] = z2;

     	cacoords[2][0] = x3;
     	cacoords[2][1] = y3;
     	cacoords[2][2] = z3;

     	cacoords[3][0] = x4;
     	cacoords[3][1] = y4;
     	cacoords[3][2] = z4;

      bin13_1 = RBINS[i][0];
      bin13_2 = RBINS[i][1];
      bin14 = RBINS[i][2];

      pro = 0;

      if (prevres && !strncmp(prevres->name,"PRO",3)) {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat_pro[j].bins[0]-bin13_1)+fabs(nco_stat_pro[j].bins[1]-bin13_2)+0.2*fabs(nco_stat_pro[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat_pro[j].bins[0]>=0 && hit>1e-3);
        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  }
     	  }
        for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  }
     	  }
      } else {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat[j].bins[0]-bin13_1)+fabs(nco_stat[j].bins[1]-bin13_2)+0.2*fabs(nco_stat[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat[j].bins[0]>=0 && hit>1e-3);
        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat[bestpos].data[j][k];
       	  }
     	  }
        for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat[bestpos].data[j][k];
       	  }
     	  }
      }

     	rmsd=superimpose2(cacoords, tmpstat, 4, tmpcoords, 8);

     	total += rmsd;
     	if (rmsd>maxrms) maxrms=rmsd;

// add-or-replace

      if (prevres) {
        add_replace(prevres, "C  ", tmpcoords[4][0], tmpcoords[4][1], tmpcoords[4][2], FLAG_BACKBONE);
        add_replace(prevres, "O  ", tmpcoords[5][0], tmpcoords[5][1], tmpcoords[5][2], FLAG_BACKBONE);
      }

      if (res) {
        add_replace(res, "N  ", tmpcoords[6][0], tmpcoords[6][1], tmpcoords[6][2], FLAG_BACKBONE);
      } else { // terminal oxygen instead of nitrogen
        add_replace(prevres, "OXT", tmpcoords[6][0], tmpcoords[6][1], tmpcoords[6][2], FLAG_BACKBONE);
      }

      prevres = res;
      if (res)
        res = res->next;
    }

    if (_VERBOSE) printf("Backbone rebuilding deviation: average = %.3f, max = %.3f\n", total/(real)chain_length, maxrms);
}

#endif


#ifdef COMPILE_ROT

typedef struct _rot_struct {
  int r13_1, r13_2, r14;
  int nc;
  real ***coords;
  struct _rot_struct *next;
} rot_struct;

rot_struct *rotamers[20];

/* this is obsolete in a standalone version of PULCHRA */
int read_rotamers(void)
{
  FILE *inp;
  char buf[1000];
  char dum[100];
  int aa, i, j, k, l, n;
  rot_struct *new_rot, *last_rot;
  real x, y, z;

    if (_VERBOSE) printf("Reading rotamer library...\n");

    inp = fopen("NEWROT","r");
    last_rot=NULL;
    while (!feof(inp)) {
      if (fgets(buf,1000,inp)==buf) {
        if (buf[0]=='A') {
          sscanf(buf,"%s %d", dum, &aa);
          if (last_rot) last_rot->next = NULL;
          last_rot = NULL;
          if (fgets(buf,1000,inp)!=buf) break;
        }
//        printf("aa: %d\n", aa);
        if (aa==20) break;
        sscanf(buf,"%d %d %d %s %d", &i, &j, &k, dum, &l);
        new_rot = (rot_struct*)calloc(sizeof(rot_struct),1);
//        printf("%d %d %d nc: %d\n", i, j, k, l);
        new_rot->r13_1 = i;
        new_rot->r13_2 = j;
        new_rot->r14 = k;
        new_rot->nc = l;
        new_rot->next = NULL;
        new_rot->coords = (real***)calloc(sizeof(real**)*l,1);
        for (i=0;i<l;i++) {
          new_rot->coords[i]=(real**)calloc(sizeof(real*)*(nheavy[aa]+1),1);
          for (j=0;j<(nheavy[aa]+1);j++) {
            new_rot->coords[i][j]=(real*)calloc(sizeof(real)*3,1);
          }
        }
        for (i=0;i<l;i++) {
          fgets(buf,1000,inp);
          for (j=0;j<(nheavy[aa]+1);j++) {
            fgets(buf,1000,inp);
            sscanf(buf,"%lf%lf%lf",&x, &y, &z);
            new_rot->coords[i][j][0]=x;
            new_rot->coords[i][j][1]=y;
            new_rot->coords[i][j][2]=z;
          }
          if (last_rot) {
            last_rot->next = new_rot;
          } else {
            rotamers[aa] = new_rot;
          }
          last_rot = new_rot;
        }
      }
    }
    fclose(inp);
}


void cross(real *v1, real *v2, real *v3)
{
  v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
  v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

void norm(real *v)
{
  real d;

    d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
}


int check_xvol(res_type *res)
{
  res_type *res2;
  atom_type *atom1, *atom2;
  real dx, dy, dz, dd;

    res2 = chain->residua;

    while (res2) {
      atom2 = res2->atoms;
      if (res!=res2) {
        while (atom2) {
          atom1 = res->atoms;
          while (atom1) {
            if (atom1->flag & FLAG_SIDECHAIN) {
              dx = atom1->x-atom2->x;
              dy = atom1->y-atom2->y;
              dz = atom1->z-atom2->z;
              dd = dx*dx+dy*dy+dz*dz;
              if (dd<(1.7*1.7)) {
                return 1;
              }
            }
            atom1=atom1->next;
          }
          atom2=atom2->next;
        }
      }
      res2=res2->next;
    }

  return 0;
}


real ***SORTED_ROTAMERS;

void rebuild_sidechains(void)
{
  FILE *out;
  res_type *res, *prevres, *testres;
  atom_type *atom, *atom1, *atom2;
  real **cacoords, **tmpcoords, **tmpstat;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real x5, y5, z5;
  real r14, r13_1, r13_2;
  real dx, dy, dz, dd;
  real hit, besthit;
  int exvol, bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14;
  real rmsd, total;
  real v1[3], v2a[3], v2b[3], v2[3], v3[3];
  int nsc, nca;
  real cax, cay, caz;
  real **lsys, **vv, **sc;
  char scn[12][4];
  rot_struct *rot;
  int ok, last_a, last_b, last_c, last_d, jpos;
  int jx, jy, jz, jxi, jyi, jzi, b13_1, b13_2, b14, jm;
  int crot, bestrot, minexvol, totexvol, rtried, pos, cpos;
  real cmx, cmy, cmz, ddx, ddy, ddz, ddd, bestdd;
  real sort_rot[100][2];

    if (_VERBOSE) printf("Rebuilding side chains...\n");

    prepare_rbins();

    lsys = (real**)calloc(sizeof(real*)*3,1);
    vv = (real**)calloc(sizeof(real*)*3,1);
    sc = (real**)calloc(sizeof(real*)*12,1);
    for (i=0;i<12;i++)
      sc[i] = (real*)calloc(sizeof(real)*3,1);
    for (i=0;i<3;i++) {
      lsys[i] = (real*)calloc(sizeof(real)*3,1);
      vv[i] = (real*)calloc(sizeof(real)*3,1);
    }

    SORTED_ROTAMERS = (real***)calloc(sizeof(real**)*(chain_length+1),1);
    for (i=0;i<chain_length+1;i++) {
      SORTED_ROTAMERS[i] = (real**)calloc(sizeof(real*)*10,1);
      for (j=0;j<10;j++) {
        SORTED_ROTAMERS[i][j] = (real*)calloc(sizeof(real)*2,1);
      }
    }

    prevres = NULL;
    res = chain->residua;
    totexvol = 0;

    for (i=0;i<chain_length;i++) {
      if (!strncmp(res->name,"GLY",3) || !res->protein) {
        if (res->next) res = res->next;
        continue;
      }

  		x1 = C_ALPHA[i-2][0];
  		y1 = C_ALPHA[i-2][1];
  		z1 = C_ALPHA[i-2][2];
  		x2 = C_ALPHA[i-1][0];
  		y2 = C_ALPHA[i-1][1];
  		z2 = C_ALPHA[i-1][2];
  		x3 = C_ALPHA[i][0];
  		y3 = C_ALPHA[i][1];
  		z3 = C_ALPHA[i][2];
  		x4 = C_ALPHA[i+1][0];
  		y4 = C_ALPHA[i+1][1];
  		z4 = C_ALPHA[i+1][2];

      bin13_1 = RBINS[i][0];
      bin13_2 = RBINS[i][1];
      bin14 = RBINS[i][2];

      v1[0] = x4-x2;
      v1[1] = y4-y2;
      v1[2] = z4-z2;

      v2a[0] = x4-x3;
      v2a[1] = y4-y3;
      v2a[2] = z4-z3;

      v2b[0] = x3-x2;
      v2b[1] = y3-y2;
      v2b[2] = z3-z2;

      cross(v2a, v2b, v2);
      cross(v1, v2, v3);

      norm(v1);
      norm(v2);
      norm(v3);

// gather 10 closest rotamer conformations...

      for (j=0;j<10;j++)
        SORTED_ROTAMERS[i][j][0] = 500.;

      j = 0;
      besthit = 1000.;
      bestpos = 0;
      do {
        if (rot_stat_idx[j][0]==res->type) {
          hit = fabs(rot_stat_idx[j][1]-bin13_1)+fabs(rot_stat_idx[j][2]-bin13_2)+0.2*fabs(rot_stat_idx[j][3]-bin14);
          if (hit<SORTED_ROTAMERS[i][9][0]) {
            k = 9;
            while (k>=0 && hit<SORTED_ROTAMERS[i][k][0]) {
              k--;
            }
            k++;
            // k = hit
            for (l=9;l>k;l--) {
              SORTED_ROTAMERS[i][l][0]=SORTED_ROTAMERS[i][l-1][0];
              SORTED_ROTAMERS[i][l][1]=SORTED_ROTAMERS[i][l-1][1];
            }
            SORTED_ROTAMERS[i][k][0]=hit;
            SORTED_ROTAMERS[i][k][1]=j;
          }
        }
        j++;
      } while (rot_stat_idx[j][0]>=0);


      besthit = SORTED_ROTAMERS[i][0][0];
      bestpos = SORTED_ROTAMERS[i][0][1];


// new rebuild...

      pos = rot_stat_idx[bestpos][5];
      nsc = nheavy[res->type]+1;

      if (_PDB_SG) { // more than one rotamer - check SC
        bestdd = 100.; crot = 0;
        for (l=0;l<2;l++) { // check two closest conformations
          cpos = SORTED_ROTAMERS[i][l][1];
          for (m=0;m<rot_stat_idx[cpos][4];m++) {
            for (j=0;j<3;j++) {
              vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
              for (k=0;k<3;k++) {
                if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
              }
            }
            pos = rot_stat_idx[cpos][5]+nsc*m;
            for (j=0;j<nsc;j++) {
              for (k=0;k<3;k++) {
                sc[j][k] = rot_stat_coords[pos+j][k];
              }
            }
            superimpose2(vv,lsys,3,sc,nsc);
            for (j=0;j<nsc;j++) {
              sc[j][0] += x3;
              sc[j][1] += y3;
              sc[j][2] += z3;
            }
            cmx = 0.; cmy = 0.; cmz = 0.;
            for (j=0;j<nsc;j++) {
              cmx += sc[j][0];
              cmy += sc[j][1];
              cmz += sc[j][2];
            }
            cmx /= (real) nsc;
            cmy /= (real) nsc;
            cmz /= (real) nsc;
            ddx = res->cmx-cmx;
            ddy = res->cmy-cmy;
            ddz = res->cmz-cmz;
            ddx *= ddx;
            ddy *= ddy;
            ddz *= ddz;
            ddd = ddx+ddy+ddz;
            if (ddd<bestdd) {
              bestdd = ddd;
              crot = pos; // closest rotamer position
            }
          }
        }
        pos = crot;
      } // PDB_SG

      for (j=0;j<3;j++) {
        vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
        for (k=0;k<3;k++) {
          if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
        }
      }

      for (j=0;j<nsc;j++) {
        for (k=0;k<3;k++) {
          sc[j][k] = rot_stat_coords[pos+j][k];
        }
      }

      superimpose2(vv,lsys,3,sc,nsc);

      for (j=0;j<nsc;j++) {
        sc[j][0] += x3;
        sc[j][1] += y3;
        sc[j][2] += z3;
      }

      for (j=1;j<nsc;j++) {
        add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
      }

      if (res->next) res = res->next;

    } // i++, next res

    for (i=0;i<12;i++)
      free(sc[i]);
    for (i=0;i<3;i++) {
      free(lsys[i]);
      free(vv[i]);
    }
    free(sc);
	free(lsys);
	free(vv);
}


typedef struct _atom_list {
  atom_type *atom;
  struct _atom_list *next;
} atom_list;

// finds steric conflicts

int get_conflicts(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid)
{
  atom_list *llist;
  atom_type *atom, *atom2;
  int i, j, k, x, y, z;
  int ii, jj, kk, con, iter, maxcon, merged;
  real dx, dy, dz, dd;

    con = 0;
    atom = res->atoms;
    while (atom) {
      i = atom->gx;
      j = atom->gy;
      k = atom->gz;
      for (ii=i-2;ii<=i+2;ii++)
        for (jj=j-2;jj<=j+2;jj++)
          for (kk=k-2;kk<=k+2;kk++) {
            if (ii>=0 && ii<xgrid && jj>=0 && jj<ygrid && kk>=0 && kk<zgrid) {
              llist = grid[ii][jj][kk];
              while (llist) {
                atom2 = llist->atom;
                if (atom && atom2 && res && atom2->res) {
                  merged=0;
                  if (res==atom2->res) { // self-xvol
                    if (atom->flag & FLAG_SIDECHAIN && atom2->flag & FLAG_SIDECHAIN) merged=1;
                    if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='A' && atom2->name[0]=='C' && atom2->name[1]=='B') merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='B' && atom2->name[0]=='C' && atom2->name[1]=='A') merged=1;
                    if (res->name[0]=='P') {
                      if (atom->name[0]=='C' && atom->name[1]=='D' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                      if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]=='D') merged=1;
                    }

                    if (!merged) {
//                      printf("merged: %s[%d] %s-%s %d %d\n", res->name,res->num,atom->name,atom2->name,atom->flag,atom2->flag);
                    }
                  } else
                  if (res->next==atom2->res || res==atom2->res->next) {
                    if (atom->name[0]=='C' && atom->name[1]==' ' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                    if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]==' ') merged=1;
                  }
                  if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1; // for now
                  if (atom->flag & FLAG_SCM || atom2->flag & FLAG_SCM) merged=1; // for now
                  if (!merged) {
                    dx = atom->x-atom2->x;
                    dx*=dx;
                    dy = atom->y-atom2->y;
                    dy*=dy;
                    dz = atom->z-atom2->z;
                    dz*=dz;
                    dd = dx+dy+dz;
                    if (dd<_SG_XVOL_DIST*_SG_XVOL_DIST) {
                      con++;
                    }
                  }
                }
                llist = llist->next;
              }
            }
          }
      atom = atom->next;
    }

  return con;
}

int display_conflicts(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid)
{
  atom_list *llist;
  atom_type *atom, *atom2;
  int i, j, k, x, y, z;
  int ii, jj, kk, con, iter, maxcon, merged;
  real dx, dy, dz, dd;

    con = 0;
    atom = res->atoms;
    while (atom) {
      i = atom->gx;
      j = atom->gy;
      k = atom->gz;
      for (ii=i-2;ii<=i+2;ii++)
        for (jj=j-2;jj<=j+2;jj++)
          for (kk=k-2;kk<=k+2;kk++) {
            if (ii>=0 && ii<xgrid && jj>=0 && jj<ygrid && kk>=0 && kk<zgrid) {
              llist = grid[ii][jj][kk];
              while (llist) {
                atom2 = llist->atom;
                if (atom && atom2 && res && atom2->res) {
                  merged=0;
                  if (res==atom2->res) { // self-xvol
                    if (atom->flag & FLAG_SIDECHAIN && atom2->flag & FLAG_SIDECHAIN) merged=1;
                    if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='A' && atom2->name[0]=='C' && atom2->name[1]=='B') merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='B' && atom2->name[0]=='C' && atom2->name[1]=='A') merged=1;
                    if (res->name[0]=='P') {
                      if (atom->name[0]=='C' && atom->name[1]=='D' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                      if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]=='D') merged=1;
                    }
                    if (!merged) {
//                      printf("merged: %s[%d] %s-%s %d %d\n", res->name,res->num,atom->name,atom2->name,atom->flag,atom2->flag);
                    }
                  } else
                  if (res->next==atom2->res || res==atom2->res->next) {
                    if (atom->name[0]=='C' && atom->name[1]==' ' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                    if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]==' ') merged=1;
                  }
                  if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1; // for now
                  if (atom->flag & FLAG_SCM || atom2->flag & FLAG_SCM) merged=1; // for now
                  if (!merged) {
                    dx = atom->x-atom2->x;
                    dx*=dx;
                    dy = atom->y-atom2->y;
                    dy*=dy;
                    dz = atom->z-atom2->z;
                    dz*=dz;
                    dd = dx+dy+dz;
                    if (dd<1.6*1.6) {
                      printf("STERIC CONFLICT: %s[%d]%s-%s[%d]%s\n", atom->res->name,atom->res->num,atom->name,atom2->res->name,atom2->res->num,atom2->name);
                      con++;
                    }
                  }
                }
                llist = llist->next;
              }
            }
          }
      atom = atom->next;
    }

  return con;
}


void allocate_grid(atom_list *****grid_, int *xgrid_, int *ygrid_, int *zgrid_)
{
  static int xgrid, ygrid, zgrid;
  static atom_list ****grid = NULL;
  atom_list *llist, *alist;
  real min[3], max[3];
  res_type *res, *worst;
  atom_type *atom, *atom2;
  int i, j, x, y, z;

    if (!grid && chain->residua && chain->residua->atoms) {
 	    res = chain->residua;
      min[0]=max[0]=res->atoms->x;
      min[1]=max[1]=res->atoms->y;
      min[2]=max[2]=res->atoms->z;
	    while (res) {
	      atom = res->atoms;
	      while (atom) {
	        if (atom->x<min[0]) min[0]=atom->x;
	        if (atom->y<min[1]) min[1]=atom->y;
	        if (atom->z<min[2]) min[2]=atom->z;
	        if (atom->x>max[0]) max[0]=atom->x;
	        if (atom->y>max[1]) max[1]=atom->y;
	        if (atom->z>max[2]) max[2]=atom->z;
	        atom = atom->next;
	      }
	      res = res->next;
	    }

	    xgrid = (max[0]-min[0])/GRID_RES;
	    ygrid = (max[1]-min[1])/GRID_RES;
	    zgrid = (max[2]-min[2])/GRID_RES;

	    if (_VERBOSE) printf("Allocating grid (%d %d %d)...\n", xgrid, ygrid, zgrid);

	   grid = (atom_list****)calloc(sizeof(atom_list***)*(xgrid+1),1);
	   for (i=0;i<xgrid+1;i++) {
	     grid[i] = (atom_list***)calloc(sizeof(atom_list**)*(ygrid+1),1);
	     for (j=0;j<ygrid+1;j++) {
	       grid[i][j] = (atom_list**)calloc(sizeof(atom_list*)*(zgrid+1),1);
	     }
	   }

	   res = chain->residua;
	   while (res) {
	     atom = res->atoms;
	     while (atom) {
	       x = xgrid*(atom->x-min[0])/(max[0]-min[0]);
	       y = ygrid*(atom->y-min[1])/(max[1]-min[1]);
	       z = zgrid*(atom->z-min[2])/(max[2]-min[2]);
	       alist = (atom_list*)calloc(sizeof(atom_list),1);
	       alist->atom = atom;
	       atom->gx = x;
	       atom->gy = y;
	       atom->gz = z;
	       if (grid[x][y][z]!=NULL) {
	         llist = grid[x][y][z];
	         while (llist->next) llist=llist->next;
	         llist->next = alist;
	       } else {
	         grid[x][y][z]=alist;
	       }
	       atom = atom->next;
	     }
	     res = res->next;
	   }
	} else {
	   if (_VERBOSE) printf("Grid already allocated (%d %d %d)\n", xgrid, ygrid, zgrid);
	}

	*grid_ = grid;
	*xgrid_ = xgrid;
	*ygrid_ = ygrid;
	*zgrid_ = zgrid;
}

void optimize_exvol(void)
{
  real min[3], max[3];
  res_type *res, *worst;
  atom_type *atom, *atom2;
  int xgrid, ygrid, zgrid;
  atom_list ****grid, *llist, *alist;
  int i, j, k, l, m, x, y, z;
  int ii, jj, kk, con, iter, maxcon, totcon;
  int cpos, bestpos, pos, con0;
  real v1[3], v2a[3], v2b[3], v2[3], v3[3];
  int nsc, nca;
  real cax, cay, caz;
  real **lsys, **vv, **sc;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;

    min[0]=1e5;
    min[1]=1e5;
    min[2]=1e5;
    max[0]=-1e5;
    max[1]=-1e5;
    max[2]=-1e5;

    lsys = (real**)calloc(sizeof(real*)*3,1);
    vv = (real**)calloc(sizeof(real*)*3,1);
    sc = (real**)calloc(sizeof(real*)*12,1);
    for (i=0;i<12;i++)
      sc[i] = (real*)calloc(sizeof(real)*3,1);
    for (i=0;i<3;i++) {
      lsys[i] = (real*)calloc(sizeof(real)*3,1);
      vv[i] = (real*)calloc(sizeof(real)*3,1);
    }

  	allocate_grid(&grid, &xgrid, &ygrid, &zgrid);

    if (_VERBOSE) printf("Finding excluded volume conflicts...\n", xgrid, ygrid, zgrid);

	iter = 0;

	do {
//printf("ITER: %d\n", iter);

      maxcon = 0;
      totcon=0;

      res = chain->residua;
      while (res) {
        if (res->protein) {
          con = get_conflicts(res, grid, xgrid, ygrid, zgrid);
          if (con>0) {
            totcon+=con;
            if (con>maxcon) {
              maxcon = con;
              worst = res;
            }
          }
        }
        res = res->next;
      }

      if (_VERBOSE && iter==0) {
        printf("Total number of conflicts: %d\n", totcon);
      }

      if (totcon==0) break;

      if (_VERBOSE && iter==0) {
        printf("Maximum number of conflicts: %s[%d] : %d\n", worst->name, worst->num, maxcon);
      }

      totcon=0;

      if (maxcon>0) {

// try to fix...

    res = chain->residua;
    for (i=0;i<chain_length;i++) {
      if (!strncmp(res->name,"GLY",3) || !res->protein) {
        if (res->next) res = res->next;
        continue;
      }

      nsc = nheavy[res->type]+1;

    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];
    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];
    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];
    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];

      v1[0] = x4-x2;
      v1[1] = y4-y2;
      v1[2] = z4-z2;

      v2a[0] = x4-x3;
      v2a[1] = y4-y3;
      v2a[2] = z4-z3;

      v2b[0] = x3-x2;
      v2b[1] = y3-y2;
      v2b[2] = z3-z2;

      cross(v2a, v2b, v2);
      cross(v1, v2, v3);

      norm(v1);
      norm(v2);
      norm(v3);

      con = get_conflicts(res, grid, xgrid, ygrid, zgrid);

      if (con>0) {

        bestpos=0;
        con0 = 100;
        for (l=0;l<10;l++) { // check two closest conformations
          cpos = SORTED_ROTAMERS[i][l][1];
          for (m=0;m<rot_stat_idx[cpos][4];m++) {
            for (j=0;j<3;j++) {
              vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
              for (k=0;k<3;k++) {
                if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
              }
            }
            pos = rot_stat_idx[cpos][5]+nsc*m;
            for (j=0;j<nsc;j++) {
              for (k=0;k<3;k++) {
                sc[j][k] = rot_stat_coords[pos+j][k];
              }
            }
            superimpose2(vv,lsys,3,sc,nsc);
            for (j=0;j<nsc;j++) {
              sc[j][0] += x3;
              sc[j][1] += y3;
              sc[j][2] += z3;
            }
            for (j=1;j<nsc;j++) {
              add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
            }
            con = get_conflicts(res, grid, xgrid, ygrid, zgrid);
//printf("test: %d\n", con);

            if (con<con0) {
              con0 = con;
              bestpos = pos;
            }
            if (con==0) break;
          }
          if (con==0) break;
        }

		totcon += con0;

        for (j=0;j<3;j++) {
          vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
          for (k=0;k<3;k++) {
            if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
          }
        }
        pos = bestpos;
        for (j=0;j<nsc;j++) {
          for (k=0;k<3;k++) {
            sc[j][k] = rot_stat_coords[pos+j][k];
          }
        }
        superimpose2(vv,lsys,3,sc,nsc);
        for (j=0;j<nsc;j++) {
          sc[j][0] += x3;
          sc[j][1] += y3;
          sc[j][2] += z3;
        }
        for (j=1;j<nsc;j++) {
          add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
        }
      }



      	res=res->next;

    	} // i
	}

	iter++;

	} while (iter<_XVOL_ITER);


    if (_VERBOSE) {
      if (totcon>0)
        printf("WARNING: %d steric conflict(s) are still there.\n", totcon);
      else
        printf("All steric conflicts removed.\n");
    }

    for (i=0;i<12;i++)
      free(sc[i]);
    for (i=0;i<3;i++) {
      free(lsys[i]);
      free(vv[i]);
    }
    free(sc);
	free(lsys);
	free(vv);


}


void vcross(real ax,real ay,real az,real bx,real by,real bz,real *cx,real *cy,real *cz)
{
    *cx = ay * bz - by * az;
    *cy = az * bx - bz * ax;
    *cz = ax * by - bx * ay;
}

real vdot(real ax,real ay,real az,real bx,real by,real bz)
{
    return ax*bx+ay*by+az*bz;
}

real calc_torsion(atom_type *a1, atom_type *a2, atom_type *a3, atom_type *a4)
{
  real v12x, v12y, v12z;
  real v43x, v43y, v43z;
  real zx, zy, zz;
  real px, py, pz;
  real xx, xy, xz;
  real yx, yy, yz;
  real u, v, angle;

    v12x = a1->x-a2->x;
    v12y = a1->y-a2->y;
    v12z = a1->z-a2->z;

    v43x = a4->x-a3->x;
    v43y = a4->y-a3->y;
    v43z = a4->z-a3->z;

    zx = a2->x-a3->x;
    zy = a2->y-a3->y;
    zz = a2->z-a3->z;

    vcross(zx,zy,zz,v12x,v12y,v12z,&px,&py,&pz);
    vcross(zx,zy,zz,v43x,v43y,v43z,&xx,&xy,&xz);
    vcross(zx,zy,zz,xx,xy,xz,&yx,&yy,&yz);

    u = vdot(xx,xy,xz,xx,xy,xz);
    v = vdot(yx,yy,yz,yx,yy,yz);

    angle = 360.;

    if (u<0. || v<0.) return angle;

    u = vdot(px,py,pz,xx,xy,xz) / sqrt(u);
    v = vdot(px,py,pz,yx,yy,yz) / sqrt(v);

    if (u != 0.0 || v != 0.0) angle = atan2(v, u) * RADDEG;


  return angle;

}


// Ca-N-C-Cb angle should be close to 34 deg
// check and fix 

int chirality_check(void)
{
  int i;
  atom_type *a_ca, *a_n, *a_c, *a_cb;
  atom_type *atom;
  res_type *res;
  real angle;
  real nx, ny, nz;
  real px, py, pz;
  real qx, qy, qz;
  real rx, ry, rz;
  real xx, xy, xz;
  real yx, yy, yz;
  real dd, costheta, sintheta;

    if (_VERBOSE) printf("Checking chirality...\n");
    res = chain->residua;
    while (res) {
      a_ca = a_n = a_c = a_cb = NULL;
      a_ca = find_atom(res,"CA ");
      a_n = find_atom(res,"N  ");
      a_c = find_atom(res,"C  ");
      a_cb = find_atom(res,"CB ");
      if (a_ca && a_n && a_c && a_cb) {
        angle = calc_torsion(a_ca, a_n, a_c, a_cb);
        if (angle<0.) {
          if (_VERBOSE) printf("WARNING: D-aa detected at %s %3d : %5.2f", res->name, res->num, angle);
          xx = a_ca->x-a_n->x;
          xy = a_ca->y-a_n->y;
          xz = a_ca->z-a_n->z;
          yx = a_c->x-a_ca->x;
          yy = a_c->y-a_ca->y;
          yz = a_c->z-a_ca->z;
          vcross(xx,xy,xz,yx,yy,yz,&nx,&ny,&nz);
          dd = sqrt(nx*nx+ny*ny+nz*nz);
          nx /= dd;
          ny /= dd;
          nz /= dd;
          // nx, ny, nz = reflection plane normal
          rx = xx-yx;
          ry = xy-yy;
          rz = xz-yz;
          dd = sqrt(rx*rx+ry*ry+rz*rz);
          rx /= dd;
          ry /= dd;
          rz /= dd;
          costheta = -1.;
          sintheta = 0.;
          atom = res->atoms;
          while (atom) {
            if (atom->flag & FLAG_SIDECHAIN) {
              px = atom->x-a_ca->x;
              py = atom->y-a_ca->y;
              pz = atom->z-a_ca->z;
              qx = qy = qz = 0.;
              qx += (costheta + (1 - costheta) * rx * rx) * px;
              qx += ((1 - costheta) * rx * ry - rz * sintheta) * py;
              qx += ((1 - costheta) * rx * rz + ry * sintheta) * pz;
              qy += ((1 - costheta) * rx * ry + rz * sintheta) * px;
              qy += (costheta + (1 - costheta) * ry * ry) * py;
              qy += ((1 - costheta) * ry * rz - rx * sintheta) * pz;
              qz += ((1 - costheta) * rx * rz - ry * sintheta) * px;
              qz += ((1 - costheta) * ry * rz + rx * sintheta) * py;
              qz += (costheta + (1 - costheta) * rz * rz) * pz;
              qx += a_ca->x;
              qy += a_ca->y;
              qz += a_ca->z;
              atom->x = qx;
              atom->y = qy;
              atom->z = qz;
            }
            atom = atom->next;
          }
          angle = calc_torsion(a_ca, a_n, a_c, a_cb);
          if (_VERBOSE) printf(", fixed : %5.2f\n", angle);
        }
      }
      res = res->next;
    }
}


#endif

// DSSP energy of petide-peptide HB

real hb_energy(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid)
{
	atom_type *atom, *c_atom1, *o_atom1, *n_atom1, *c_atom2, *o_atom2, *n_atom2, *tmp_atom;
	atom_type h_atom;
	int i, j, k, ii, jj, kk;
  atom_list *llist, *alist;
  real dx, dy, dz, dist, min_dist1, min_dist2;
  real hx1, hy1, hz1, dd;
  real dno, dnc, dho, dhc;
  real ene, Q;

    ene = 1e3;

    if (!res || !res->prev) return ene;
            
    Q = -27888.0; // DSSP h-bond energy constant

		c_atom1 = o_atom1 = n_atom1 = NULL;

  	atom = res->prev->atoms;
  	while (atom) {
			if (atom->name[0]=='C' && atom->name[1]==' ') c_atom1 = atom;
			if (atom->name[0]=='O' && atom->name[1]==' ') o_atom1 = atom;
			atom = atom->next;
		}

  	atom = res->atoms;
  	while (atom) {
			if (atom->name[0]=='N' && atom->name[1]==' ') { n_atom1 = atom; break; }
			atom = atom->next;
		}

// first bond

    min_dist2 = 1e10;
    o_atom2 = c_atom2 = NULL;
		if (n_atom1) {
			i = n_atom1->gx;
			j = n_atom1->gy;
			k = n_atom1->gz;
			for (ii=i-1;ii<=i+1;ii++) {
				for (jj=j-1;jj<=j+1;jj++) {
					for (kk=k-1;kk<=k+1;kk++) {
						if (ii>=0 && ii<xgrid && jj>=0 && jj<ygrid && kk>=0 && kk<=zgrid) {
							llist = grid[ii][jj][kk];
							while (llist) {
							  if (llist->atom->name[0]=='O' && llist->atom->name[1]==' ' && abs(llist->atom->res->locnum-n_atom1->res->locnum)>2)  {
							    tmp_atom = llist->atom;
							    dx = n_atom1->x-tmp_atom->x;
							    dy = n_atom1->y-tmp_atom->y;
							    dz = n_atom1->z-tmp_atom->z;
							    dist = dx*dx+dy*dy+dz*dz;
							    if (dist<min_dist2 && dist<25.0) {
							      o_atom2=tmp_atom;
							      min_dist2 = dist;
   						    } 
							  } 
							  llist = llist->next;
		          }					
						}
					}
				}
			}
		}

    if (o_atom2) {
      atom = o_atom2->res->atoms;
    	while (atom) {
  			if (atom->name[0]=='C' && atom->name[1]==' ') { c_atom2 = atom; break; }
  			atom = atom->next;
  		}                      
      if (c_atom2) {    
    		hx1 = o_atom1->x-c_atom1->x;
    		hy1 = o_atom1->y-c_atom1->y;
    		hz1 = o_atom1->z-c_atom1->z;
    		dd = -1.081f/sqrt(hx1*hx1+hy1*hy1+hz1*hz1);
    		hx1 *= dd;
    		hy1 *= dd;
    		hz1 *= dd;
    		
    		hx1 += n_atom1->x;
    		hy1 += n_atom1->y;
    		hz1 += n_atom1->z;
        
        add_replace(n_atom1->res, "H  ", hx1, hy1, hz1, FLAG_BACKBONE);

  // dno
        dx = n_atom1->x-o_atom2->x;
        dy = n_atom1->y-o_atom2->y;
        dz = n_atom1->z-o_atom2->z;
        dno = sqrt(dx*dx+dy*dy+dz*dz);
  
  // dnc
        dx = n_atom1->x-c_atom2->x;
        dy = n_atom1->y-c_atom2->y;
        dz = n_atom1->z-c_atom2->z;
        dnc = sqrt(dx*dx+dy*dy+dz*dz);
  
  // dho
        dx = hx1-o_atom2->x;
        dy = hy1-o_atom2->y;
        dz = hz1-o_atom2->z;
        dho = sqrt(dx*dx+dy*dy+dz*dz);
  
  // dhc
        dx = hx1-c_atom2->x;
        dy = hy1-c_atom2->y;
        dz = hz1-c_atom2->z;
        dhc = sqrt(dx*dx+dy*dy+dz*dz);
        if (dho<0.01F || dhc<0.01F || dnc<0.01F || dno<0.01F) {
          ene = -10.0;
        } else {
          ene = 0.001*(Q/dho - Q/dhc + Q/dnc - Q/dno);
        }       
      }
    }

  return ene;
}

// rotates a point around a vector
void rot_point_vector(real *x, real *y, real *z, real u, real v, real w, real angle)
{
  real ux, uy, uz, vx, vy, vz, wx, wy, wz, sa, ca;
  
    sa = sinf(10.0*M_PI*angle/180.0);
    ca = cosf(10.0*M_PI*angle/180.0);
    
    ux = u**x;
    uy = u**y;
    uz = u**z;
    vx = v**x;
    vy = v**y;
    vz = v**z;
    wx = w**x;
    wy = w**y;
    wz = w**z;

    *x = u*(ux+vy+wz)+(*x*(v*v+w*w)-u*(vy+wz))*ca+(-wy+vz)*sa;
    *y = v*(ux+vy+wz)+(*y*(u*u+w*w)-v*(ux+wz))*ca+( wx-uz)*sa;
    *z = w*(ux+vy+wz)+(*z*(u*u+v*v)-w*(ux+vy))*ca+(-vx+uy)*sa;
}


// rotates a peptide plate

void rot_peptide(res_type *res, real angle)
{
	atom_type *atom, *c_atom, *o_atom, *n_atom, *ca_atom1, *ca_atom2;
  real u, v, w, x, y, z, dd;

    if (!res || !res->prev) return;
      
    c_atom = o_atom = n_atom = ca_atom1 = ca_atom2 = NULL;
    
  	atom = res->prev->atoms;
  	while (atom) {
			if (atom->name[0]=='C' && atom->name[1]=='A') ca_atom1 = atom;
			if (atom->name[0]=='C' && atom->name[1]==' ') c_atom = atom;
			if (atom->name[0]=='O' && atom->name[1]==' ') o_atom = atom;
			atom = atom->next;
		}
		
    atom = res->atoms;
  	while (atom) {
			if (atom->name[0]=='C' && atom->name[1]=='A') ca_atom2 = atom;
			if (atom->name[0]=='N' && atom->name[1]==' ') n_atom = atom;
			atom = atom->next;
		}
		
    if (c_atom && o_atom && n_atom && ca_atom1 && ca_atom2) {
      u = ca_atom2->x-ca_atom1->x;
      v = ca_atom2->y-ca_atom1->y;
      w = ca_atom2->z-ca_atom1->z;
      dd = 1.0f/sqrt(u*u+v*v+w*w); 
      u*=dd; v*=dd; w*=dd; // normalize ca-ca vector
      x = n_atom->x-ca_atom1->x;
      y = n_atom->y-ca_atom1->y;
      z = n_atom->z-ca_atom1->z;
      rot_point_vector(&x, &y, &z, u, v, w, angle);
      n_atom->x = x+ca_atom1->x;
      n_atom->y = y+ca_atom1->y;
      n_atom->z = z+ca_atom1->z;
      x = c_atom->x-ca_atom1->x;
      y = c_atom->y-ca_atom1->y;
      z = c_atom->z-ca_atom1->z;
      rot_point_vector(&x, &y, &z, u, v, w, angle);
      c_atom->x = x+ca_atom1->x;
      c_atom->y = y+ca_atom1->y;
      c_atom->z = z+ca_atom1->z;
      x = o_atom->x-ca_atom1->x;
      y = o_atom->y-ca_atom1->y;
      z = o_atom->z-ca_atom1->z;
      rot_point_vector(&x, &y, &z, u, v, w, angle);
      o_atom->x = x+ca_atom1->x;
      o_atom->y = y+ca_atom1->y;
      o_atom->z = z+ca_atom1->z;
    }  		
  
}

void optimize_backbone(mol_type *chain)
{
  int xgrid, ygrid, zgrid;
  atom_list ****grid;
  atom_type *atom;
	res_type *res;
  real ene, min_ene, tot1, tot2;
  int i, k, best;
FILE *out;
   
  	if (_VERBOSE) printf("Optimizing backbone...\n", xgrid, ygrid, zgrid);

		allocate_grid(&grid, &xgrid, &ygrid, &zgrid);

    tot1 = tot2 = 0.0;
    
		res = chain->residua;
		while (res) {
			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
			if (ene<-0.5) tot1 += ene;
			res = res->next;
		}

		res = chain->residua;
		while (res) {
		  if (res->type!=7) {
  			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
  			if (ene<1.0) { // try to optimize
  			  min_ene = ene;
  			  rot_peptide(res, -1.1);
  			  best = 0;
  			  for (i=-10;i<10;i++) {
  			    rot_peptide(res, 0.1);
      			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
  			    if (ene<min_ene) {
  			      best = i;
  			      min_ene = ene;
  			    }			  
  			  }	
  			  rot_peptide(res,-0.9);		
    			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
  			  if (min_ene<ene) {
  			    rot_peptide(res,0.1*best);
      			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
  			  } 
  			}			
  		}
			res = res->next;
		}

		res = chain->residua;
		while (res) {
			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
			if (ene<-0.5) tot2 += ene;
			res = res->next;
		}

    if (_VERBOSE) printf("Backbone HB energy: before %g, after: %g, difference: %g\n", tot1, tot2, tot2-tot1);

}

int main(int argc, char **argv)
{
  int i, j, next;
  char buf[100];
  char *name=NULL, *ini_name=NULL;
  char *ptr, out_name[1000];
  real f;
  mol_type *mol;
  struct timeb time0, time1;

    for (i=1; i<argc; i++) {
      if (argv[i][0]=='-') {
        next=0;
        for (j=1; j<(int)strlen(argv[i]); j++) {
          switch(argv[i][j]) {
            case 'v': _VERBOSE=1; break;
            case 'c': _CA_OPTIMIZE=0; break;
            case 'e': _BB_REARRANGE=1; break;
            case 'r': _CA_RANDOM=1; break;
            case 'z': _CHIRAL=0; break;
            case 't': _CA_TRAJECTORY=1; break;
            case 'n': _CENTER_CHAIN=1; break;
            case 'b': _REBUILD_BB=0; break;
            case 's': _REBUILD_SC=0; break;
            case 'i': ini_name = argv[++i]; next=1; break;
            case 'g': _PDB_SG=1; break;
            case 'x': _TIME_SEED=1; break;
            case 'o': _XVOLUME=0; break;
            case 'h': _REBUILD_H=0; break;
            case 'q': _BB_OPTIMIZE=1; break;
            case 'p': _CISPRO=1; break;
            case 'f': _PRESERVE=1; break;
            case 'u':
              if (sscanf(argv[++i],"%lf",&f)==1) {
                _CA_START_DIST = f;
              }
              next=1;
            break;
            default: {
              printf("Unknown option: %c\n", argv[i][j]);
              return -1;
            }
          }
          if (next) break;
        }
      } else {
        if (!name) name=argv[i];
      }
    }

    if (!name) {
      printf("PULCHRA Protein Chain Restoration Algorithm version %4.2f\n", PULCHRA_VERSION);
      printf("Usage: %s [options] <pdb_file>\n", argv[0]);
      printf("The program default input is a PDB file.\n");
      printf("Output file <pdb_file.rebuild.pdb> will be created as a result.\n");
      printf("Valid options are:\n\n");
      printf("  -v : verbose output (default: off)\n");
      printf("  -n : center chain (default: off)\n");
      printf("  -x : time-seed random number generator (default: off)\n");
      printf("  -g : use PDBSG as an input format (CA=C-alpha, SC or CM=side chain c.m.)\n\n");
      printf("  -c : skip C-alpha positions optimization (default: on)\n");
      printf("  -p : detect cis-prolins (default: off)\n");
      printf("  -r : start from a random chain (default: off)\n");
      printf("  -i pdbfile : read the initial C-alpha coordinates from a PDB file\n");
      printf("  -t : save chain optimization trajectory to file <pdb_file.pdb.trajectory>\n");
      printf("  -u value : maximum shift from the restraint coordinates (default: 0.5A)\n\n");
      printf("  -e : rearrange backbone atoms (C, O are output after side chain) (default: off)\n");
      printf("  -f : preserve initial coordinates (default: off, implies '-c' on and '-n' off)\n");

#ifdef COMPILE_BB
      printf("  -b : skip backbone reconstruction (default: on)\n");
      printf("  -q : optimize backbone hydrogen bonds pattern (default: off)\n");
      printf("  -h : outputs hydrogen atoms (default: off)\n");
#endif

#ifdef COMPILE_ROT
      printf("  -s : skip side chains reconstruction (default: on)\n");
      printf("  -o : don't attempt to fix excluded volume conflicts (default: on)\n");
      printf("  -z : don't check amino acid chirality (default: on)\n");
#endif
	  printf("\n");
      return -1;
    }

    for (i=0; i<255; i++) /* prepare hash table*/
      AA_NUMS[i] = 20; /* dummy aa code */
    for (i=0; i<20; i++)
      AA_NUMS[(int)SHORT_AA_NAMES[i]] = i;

    setbuf(stdout,0);

    if (_TIME_SEED) srand(time(NULL)); else srand(1234);

    if (_VERBOSE) printf("PULCHRA Protein Chain Restoration Algorithm version %4.2f\n", PULCHRA_VERSION);

    ftime(&time0);

    chain = new_mol();

    if (read_pdb_file(name,chain,"chain")==FILE_NOT_FOUND) {
      if (_VERBOSE) printf("Can't read the input file!\n");
      return -1;
    }

    if (_VERBOSE) printf("%d residua read.\n", chain->nres);

    if (_PRESERVE) printf("Initial coordinates will be preserved.\n");

    chain_length = chain->nres;

    if (_CA_OPTIMIZE && !_PRESERVE) {
      snprintf(out_name,1000,"%s.tra",name);
      ca_optimize(out_name, ini_name);
    }

#ifdef COMPILE_BB
    if (_REBUILD_BB) {
      rebuild_backbone();
      if (_BB_OPTIMIZE) {
        optimize_backbone(chain);
      }
    }
#endif

#ifdef COMPILE_ROT
    if (_REBUILD_SC) {
      rebuild_sidechains();
      if (_XVOLUME)
        optimize_exvol();
      if (_CHIRAL)
        chirality_check();
    }
#endif

    if (_CENTER_CHAIN) {
      center_chain(chain);
    }


    if (_BB_REARRANGE) {
      if (_VERBOSE) printf("Rearranging backbone atoms...\n", out_name);
    }

    ptr = strstr(name,".pdb");
    if (ptr) ptr[0]=0;
      
    snprintf(out_name,1000,"%s.rebuilt.pdb",name);
    
    if (_VERBOSE) printf("Writing output file %s...\n", out_name);
    write_pdb(out_name, chain);

    ftime(&time1);

    if (_VERBOSE) printf("Done. Reconstruction finished in %.3f s.\n", (real)0.001*(1000.*(time1.time-time0.time)+(time1.millitm-time0.millitm)));

  return 0;
}

