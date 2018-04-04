// ============ GENALG.H ================

#ifndef GENALG_H
#define GENALG_H

// ================= Include =================

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "pberr.h"
#include "pbmath.h"
#include "gset.h"

// ================= Define ==================

#define GENALG_NBENTITIES 100
#define GENALG_NBELITES 20
#define GENALG_INBREEDINGTHRESHOLD 0.1

// ------------- GenAlgAdn

// ================= Data structure ===================

typedef struct GenAlg GenAlg;

typedef struct GenAlgAdn {
  // ID
  int _id;
  // Age
  int _age;
  // Adn for floating point value
  VecFloat* _adnF;
  // Delta Adn during mutation
  VecFloat* _deltaAdnF;
  // Adn for integer point value
  VecShort* _adnI;
} GenAlgAdn;

// ================ Functions declaration ====================

// Create a new GenAlgAdn with ID 'id', 'lengthAdnF' and 'lengthAdnI'
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlgAdn* GenAlgAdnCreate(int id, int lengthAdnF, 
  int lengthAdnI);

// Free memory used by the GenAlgAdn 'that'
void GenAlgAdnFree(GenAlgAdn** that);

// Return the adn for floating point values of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
VecFloat* GAAdnAdnF(GenAlgAdn* that);

// Return the delta of adn for floating point values of the 
// GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
VecFloat* GAAdnDeltaAdnF(GenAlgAdn* that);

// Return the adn for integer values of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
VecShort* GAAdnAdnI(GenAlgAdn* that);

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga'
void GAAdnInit(GenAlgAdn* that, GenAlg* ga);

// Get the 'iGene'-th gene of the adn for floating point values of the
// GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
float GAAdnGetGeneF(GenAlgAdn* that, int iGene);

// Get the delta of the 'iGene'-th gene of the adn for floating point 
// values of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
float GAAdnGetDeltaGeneF(GenAlgAdn* that, int iGene);

// Get the 'iGene'-th gene of the adn for int values of the
// GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
int GAAdnGetGeneI(GenAlgAdn* that, int iGene);

// Set the 'iGene'-th gene of the adn for floating point values of the
// GenAlgAdn 'that' to 'gene'
#if BUILDMODE != 0
inline
#endif
void GAAdnSetGeneF(GenAlgAdn* that, int iGene, float gene);

// Set the delta of the 'iGene'-th gene of the adn for floating point 
// values of the GenAlgAdn 'that' to 'delta'
#if BUILDMODE != 0
inline
#endif
void GAAdnSetDeltaGeneF(GenAlgAdn* that, int iGene, float delta);

// Set the 'iGene'-th gene of the adn for int values of the
// GenAlgAdn 'that'to 'gene'
#if BUILDMODE != 0
inline
#endif
void GAAdnSetGeneI(GenAlgAdn* that, int iGene, short gene);

// Get the id of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
int GAAdnGetId(GenAlgAdn* that);

// Get the age of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
int GAAdnGetAge(GenAlgAdn* that);

// Print the information about the GenAlgAdn 'that' on the 
// stream 'stream'
void GAAdnPrintln(GenAlgAdn* that, FILE* stream);

// ------------- GenAlg

// ================= Define ===================

#define GABestAdnF(that) GAAdnAdnF(GAAdn(that, 0))
#define GABestAdnI(that) GAAdnAdnI(GAAdn(that, 0))

// ================= Data structure ===================

typedef struct GenAlg {
  // GSet of GenAlgAdn, sortval == score so the head of the set is the 
  // worst adn and the tail of the set is the best
  GSet* _adns;
  // Current epoch
  int _curEpoch;
  // Nb elite entities in population
  int _nbElites;
  // Id of the next new GenAlgAdn
  int _nextId;
  // Length of adn for floating point value
  int _lengthAdnF;
  // Length of adn for integer value
  int _lengthAdnI;
  // Bounds (min, max) for floating point values adn
  VecFloat2D* _boundsF;
  // Bounds (min, max) for integer values adn
  VecShort2D* _boundsI;
} GenAlg;

// ================ Functions declaration ====================

// Create a new GenAlg with 'nbEntities', 'nbElites', 'lengthAdnF' 
// and 'lengthAdnI'
// 'nbEntities' must greater than 2
// 'nbElites' must greater than 1
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlg* GenAlgCreate(int nbEntities, int nbElites, int lengthAdnF, 
  int lengthAdnI);

// Free memory used by the GenAlg 'that'
void GenAlgFree(GenAlg** that);

// Return the GSet of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
GSet* GAAdns(GenAlg* that);

// Return the nb of entities of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetNbAdns(GenAlg* that);

// Return the nb of elites of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetNbElites(GenAlg* that);

// Return the current epoch of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetCurEpoch(GenAlg* that);

// Set the nb of entities of the GenAlg 'that' to 'nb'
// 'nb' must be greater than 1, if 'nb' is lower than the current nb 
// of elite the number of elite is set to 'nb' - 1
void GASetNbEntities(GenAlg* that, int nb);

// Set the nb of elites of the GenAlg 'that' to 'nb'
// 'nb' must be greater than 0, if 'nb' is greater or equal to the 
// current nb of entities the number of entities is set to 'nb' + 1
void GASetNbElites(GenAlg* that, int nb);

// Get the length of adn for floating point value
#if BUILDMODE != 0
inline
#endif
int GAGetLengthAdnFloat(GenAlg* that);

// Get the length of adn for integer value
#if BUILDMODE != 0
inline
#endif
int GAGetLengthAdnInt(GenAlg* that);

// Get the bounds for the 'iGene'-th gene of adn for floating point 
// values
#if BUILDMODE != 0
inline
#endif
VecFloat2D* GABoundsAdnFloat(GenAlg* that, int iGene);

// Get the bounds for the 'iGene'-th gene of adn for integer values
#if BUILDMODE != 0
inline
#endif
VecShort2D* GABoundsAdnInt(GenAlg* that, int iGene);

// Set the bounds for the 'iGene'-th gene of adn for floating point 
// values to a copy of 'bounds'
#if BUILDMODE != 0
inline
#endif
void GASetBoundsAdnFloat(GenAlg* that, int iGene, VecFloat2D* bounds);

// Set the bounds for the 'iGene'-th gene of adn for integer values
// to a copy of 'bounds'
#if BUILDMODE != 0
inline
#endif
void GASetBoundsAdnInt(GenAlg* that, int iGene, VecShort2D* bounds);

// Get the GenAlgAdn of the GenAlg 'that' currently at rank 'iRank'
#if BUILDMODE != 0
inline
#endif
GenAlgAdn* GAAdn(GenAlg* that, int iRank);

// Init the GenAlg 'that'
// Must be called after the bounds have been set
// The random generator must have been initialised before calling this
// function
void GAInit(GenAlg* that);

// Step an epoch for the GenAlg 'that' with the current ranking of
// GenAlgAdn
void GAStep(GenAlg* that); 

// Print the information about the GenAlg 'that' on the stream 'stream'
void GAPrintln(GenAlg* that, FILE* stream);

// Get the level of inbreeding of curent entities of the GenAlg 'that'
// The return value is in [0.0, 1.0]
// 0.0 means all the elite entities have exactly the same adns 
float GAGetInbreeding(GenAlg* that);

// Load the GenAlg 'that' from the stream 'stream'
// If the GenAlg is already allocated, it is freed before loading
// Return true in case of success, else false
bool GALoad(GenAlg** that, FILE* stream);

// Save the GenAlg 'that' to the stream 'stream'
// Return true in case of success, else false
bool GASave(GenAlg* that, FILE* stream);

// Set the value of the GenAlgAdn 'adn' of the GenAlg 'that' to 'val'
#if BUILDMODE != 0
inline
#endif
void GASetAdnValue(GenAlg* that, GenAlgAdn* adn, float val);

// ================= Polymorphism ==================

// ================ Inliner ====================

#if BUILDMODE != 0
#include "genalg-inline.c"
#endif


#endif
