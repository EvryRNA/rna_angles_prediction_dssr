#ifndef _X3DNA_H
#define _X3DNA_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#define NR_END 1
#define FREE_ARG char*
#define NO_MATCH -1L
#define DUMMY -1L
#define TRUE 1L
#define FALSE 0L
#define BUF32 32
#define BUF512 512
#define BUF1K 1024
#define BUF2K 2048
#define BUFBIG 8192
#define UNUSED_PARAMETER(x) (void)(x)
#define NELEMS(x) ((sizeof (x))/(sizeof ((x)[0])))
enum slre_option { SLRE_CASE_SENSITIVE = 0, SLRE_CASE_INSENSITIVE = 1 };
enum slre_capture { SLRE_STRING, SLRE_INT, SLRE_FLOAT };
typedef struct {
    long min_base_hb;
    double hb_lower;
    double hb_dist1;
    double hb_dist2;
    char hb_atoms[BUF512];
    long hb_idx[BUF512];
    char alt_list[BUF512];
    double max_dorg;
    double min_dorg;
    double max_dv;
    double min_dv;
    double max_plane_angle;
    double min_plane_angle;
    double max_dNN;
    double min_dNN;
    double helix_break;
    double std_curved;
    double water_dist;
    double water_dlow;
    char water_atoms[BUF512];
    long water_idx[BUF512];
    double o3p_dist;
} miscPars;
typedef struct {
    long DEBUG;
    long VERBOSE;
    long NUM_ELE;
    long CHAIN_CASE;
    long ALL_MODEL;
    long ATTACH_RESIDUE;
    long THREE_LETTER_NTS;
    long PDBV3;
    long ORIGINAL_COORDINATE;
    long OCCUPANCY;
    long HEADER;
    long mmcif;
    double NT_CUTOFF;
    char X3DNA_VER[BUF512];
    char X3DNA_HOMEDIR[BUF512];
    char CHAIN_MARKERS[BUF512];
    char REBUILD_CHAIN_IDS[BUF512];
    char *PROGNAME;
    char **ATOM_NAMES;
    long NUM_SATOM;
    char **ATOMLIST;
    long NUM_SBASE;
    char **BASELIST;
    char **AtomName0;
    char **ResName0;
    long Name0;
    long label_RC8_YC6;
    miscPars misc_pars;
} struct_Gvars;
extern struct_Gvars Gvars;
#define DEBUG_LEVEL 6
#define SEPC '\t'
#define WSPACES " ,\t\n"
#define SKIPS "#\0"
#define UNKATM "XX"
#define DUMSTR "XXXXXX"
#define SFACTOR    2L
#define PI 3.141592653589793
#define XEPS 1.0e-7
#define XBIG 1.0e+18
#define XBIG_CUTOFF 1.0e+16
#define MFACTOR 10000.0
#define NMISC 34
#define NATOMCOL 11
#define NBASECOL 7
#define PS_DFTSIZE 500
#define PS_BOUND 10
#define FIG_DFTSIZE 8333
#define FIG_BOUND 166
#define PAR_FILE "misc_3dna.par"
#define BASE_FILE "baselist.dat"
#define ATOM_FILE "atomlist.dat"
#define HELP3DNA "help3dna.dat"
#define REF_FILE "ref_frames.dat"
#define MREF_FILE "mref_frames.dat"
#define POC_FILE "poc_haxis.r3d"
#define MUL_FILE "multiplets.pdb"
#define ALLP_FILE "allpairs.pdb"
#define BESTP_FILE "bestpairs.pdb"
#define STACK_FILE "stacking.pdb"
#define HSTACK_FILE "hstacking.pdb"
#define HLXREG_FILE "hel_regions.pdb"
#define BPORDER_FILE "bp_order.dat"
#define COLCHN_FILE "col_chains.scr"
#define COLHLX_FILE "col_helices.scr"
#define AUX_FILE "auxiliary.par"
#define BPSTEP_FILE "bp_step.par"
#define HLXSTEP_FILE "bp_helical.par"
#define SEVEN_FILE "cf_7methods.par"
#define HB_FILE "hbonds_info.dat"
#define LKG_FILE "bonds_lkg.dat"
#define SNUM_FILE "serial_num.pdb"
#define ATOMALC_FILE "atom_lkg.alc"
#define BBLKALC_FILE "bblk_lkg.alc"
#define TMP_FILE "tmp_file"
#define MULBP_FILE "mulbp.inp"
#define ROTMAT_FILE "rotmat.dat"
#define VIEW1_FILE "pmiview1"
#define VIEW2_FILE "pmiview2"
#define VIEW3_FILE "pmiview3"
#define NP 101L
#define BOND_UPPER_LIMIT 2.5
#define HTWIST0 0.05
#define BOND_FACTOR 1.15
#define NBOND_FNUM  2.0
#define NON_WC_IDX  6
#define AXIS_LENGTH 3.5
#define NUM_BASE_ATOMS BUF512
#define NUM_RESIDUE_ATOMS BUF512
#define NUM_DINUCLEOTIDE_ATOMS BUFBIG
#define EMPTY_NUMBER -9999.99
#define EMPTY_CRITERION -9999
#define MAXBASE 30000
#define NELE 12
#define O3P_UPPER 2.5
#define RTNNUM 37
#define PSTNUM 29
#define OLCRT 1.2
#define MBASES 50
#define MAXCH 100
#define END_STACK_XANG 125.0
#define MAXCLEN 52
#define RA_LIST " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
#define WC_LIST "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
#define CX_LIST "ACGITUX"
#define CB_LIST "ACGITU"
#define NT_LIST "  A", "  C", "  G", "  I", "  T", "  U", \
                "ADE", "CYT", "GUA", "INO", "THY", "URA", \
                " +A", " +C", " +G", " +I", " +T", " +U"
#define WATER_LIST "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP"
#define NUM_SAA 20
#define NUM_ATM 12
#define AA_LIST "ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU",   \
                "ASP", "GLU", "LYS", "ARG",   \
                "SER", "THR", "TYR", "HIS", "CYS", "ASN", "GLN", "TRP",   \
                "GLY"
#define OVERLAP 0.01
#define PBLKALC_FILE "pblk_lkg.alc"
#define SNAP_PEP_PDB "snap_pep.pdb"
#define SNAP_PEP_ALC "snap_pep.alc"
#define TRSP_RMS 0.25
#define DUMCHAR '@'
#define LBIG 1000000
#define DNA_BASE "ACGT"
#define SNAP_AAA "snap_aaa.pdb"
#define ABLKALC_FILE "ablk_lkg.alc"
#define SNAP_NTS "snap_nts.pdb"
#define SNAP_OPTS "snap_options"
#define WITH_BASE 1
#define WITH_BKBN 2
#define WITH_BOTH 3
#define PDBX "PDBx:"
#define O2_STACK "o2_stack.dat"
#include "x3dna_fncs.h"
#endif
