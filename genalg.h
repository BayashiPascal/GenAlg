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
#include "elorank.h"

// ================= Define ==================

#define GENALG_RUNPEREPOCH 100
#define GENALG_NBENTITIES 100
#define GENALG_NBELITES 20
#define GENALG_INBREEDINGTHRESHOLD 0.1

// ------------- GenAlgEntity

// ================= Data structure ===================

typedef struct GenAlg GenAlg;

typedef struct GenAlgEntity {
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
} GenAlgEntity;

// ================ Functions declaration ====================

// Create a new GenAlgEntity with ID 'id', 'lengthAdnF' and 'lengthAdnI'
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlgEntity* GenAlgEntityCreate(int id, int lengthAdnF, 
  int lengthAdnI);

// Free memory used by the GenAlgEntity 'that'
void GenAlgEntityFree(GenAlgEntity** that);

// Return the adn for floating point values of the GenAlgEntity 'that'
#if BUILDMODE != 0
inline
#endif
VecFloat* GAEntAdnF(GenAlgEntity* that);

// Return the delta of adn for floating point values of the 
// GenAlgEntity 'that'
#if BUILDMODE != 0
inline
#endif
VecFloat* GAEntDeltaAdnF(GenAlgEntity* that);

// Return the adn for integer values of the GenAlgEntity 'that'
#if BUILDMODE != 0
inline
#endif
VecShort* GAEntAdnI(GenAlgEntity* that);

// Initialise randomly the genes of the GenAlgEntity 'that' of the 
// GenAlg 'ga'
void GAEntInit(GenAlgEntity* that, GenAlg* ga);

// Get the 'iGene'-th gene of the adn for floating point values of the
// GenAlgEntity 'that'
#if BUILDMODE != 0
inline
#endif
float GAEntGetGeneF(GenAlgEntity* that, int iGene);

// Get the delta of the 'iGene'-th gene of the adn for floating point 
// values of the GenAlgEntity 'that'
#if BUILDMODE != 0
inline
#endif
float GAEntGetDeltaGeneF(GenAlgEntity* that, int iGene);

// Get the 'iGene'-th gene of the adn for int values of the
// GenAlgEntity 'that'
#if BUILDMODE != 0
inline
#endif
int GAEntGetGeneI(GenAlgEntity* that, int iGene);

// Set the 'iGene'-th gene of the adn for floating point values of the
// GenAlgEntity 'that' to 'gene'
#if BUILDMODE != 0
inline
#endif
void GAEntSetGeneF(GenAlgEntity* that, int iGene, float gene);

// Set the delta of the 'iGene'-th gene of the adn for floating point 
// values of the GenAlgEntity 'that' to 'delta'
#if BUILDMODE != 0
inline
#endif
void GAEntSetDeltaGeneF(GenAlgEntity* that, int iGene, float delta);

// Set the 'iGene'-th gene of the adn for int values of the
// GenAlgEntity 'that'to 'gene'
#if BUILDMODE != 0
inline
#endif
void GAEntSetGeneI(GenAlgEntity* that, int iGene, short gene);

// Get the id of the GenAlgEntity 'that'
#if BUILDMODE != 0
inline
#endif
int GAEntGetId(GenAlgEntity* that);

// Get the age of the GenAlgEntity 'that'
#if BUILDMODE != 0
inline
#endif
int GAEntGetAge(GenAlgEntity* that);

// Print the information about the GenAlgEntity 'that' on the 
// stream 'stream'
void GAEntPrintln(GenAlgEntity* that, FILE* stream);

// ------------- GenAlg

// ================= Data structure ===================

typedef struct GenAlg {
  // ELORank of GenAlgEntity
  ELORank* _elo;
  // Nb runs per epoch
  int _runsPerEpoch;
  // Current run
  int _curRun;
  // Current epoch
  int _curEpoch;
  // Nb entities in population
  int _nbEntities;
  // Nb elite entities in population
  int _nbElites;
  // Id of the next new GenAlgEntity
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

// Return the ELORank of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
ELORank* GAEloRank(GenAlg* that);

// Return the nb of runs per epoch of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetRunsPerEpoch(GenAlg* that);

// Return the nb of entities of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetNbEntities(GenAlg* that);

// Return the nb of elites of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetNbElites(GenAlg* that);

// Return the current run of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetCurRun(GenAlg* that);

// Return the current epoch of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetCurEpoch(GenAlg* that);

// Set the nb of runs per epoch of the GenAlg 'that' to 'runs'
// 'runs' must be greater than 0
#if BUILDMODE != 0
inline
#endif
void GASetRunsPerEpoch(GenAlg* that, int runs);

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

// Get the GenAlgEntity of the GenAlg 'that' currently at rank 'iRank'
#if BUILDMODE != 0
inline
#endif
GenAlgEntity* GAEntity(GenAlg* that, int iRank);

// Init the GenAlg 'that'
// Must be called after the bounds have been set
// The random generator must have been initialised before calling this
// function
void GAInit(GenAlg* that);

// Step a run for the GenAlg 'that' with ranking of GenAlgEntity given
// in the GSet of GenAlgEntity 'rank' (from best to worst, ie _sortVal
// from greater to lower) 
void GAStepRun(GenAlg* that, GSet* rank);

// Step an epoch for the GenAlg 'that' with the current ranking of
// GenAlgEntity
void GAStepEpoch(GenAlg* that); 

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

// ================= Polymorphism ==================

// ================ Inliner ====================

#if BUILDMODE != 0
#include "genalg-inline.c"
#endif


#endif
