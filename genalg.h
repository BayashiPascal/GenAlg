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
#include "respublish.h"

// ================= Define ==================

#define GABestAdnF(that) GAAdnAdnF(GABestAdn(that))
#define GABestAdnI(that) GAAdnAdnI(GABestAdn(that))

#define GENALG_NBENTITIES 100
#define GENALG_NBELITES 20

#define GENALG_TXTOMETER_NBADNDISPLAYED 40
#define GENALG_TXTOMETER_LINE1 "Epoch #xxxxxx  KTEvent #xxxxxx  \n"
#define GENALG_TXTOMETER_FORMAT1 "Epoch #%06lu  KTEvent #%06lu\n"
#define GENALG_TXTOMETER_LINE2 "Id        Age      Val\n"
#define GENALG_TXTOMETER_LINE3 "xxxxxxxx  xxxxxx  +xxxxxx.xxxxxx\n"
#define GENALG_TXTOMETER_FORMAT3 "%08lu  %06lu  %+06.6f\n"
#define GENALG_TXTOMETER_LINE4 "--------------------------------\n"
#define GENALG_TXTOMETER_LINE5 "Diversity +xxxxxx.xxxxxx        \n"
#define GENALG_TXTOMETER_FORMAT5 "Diversity %+06.6f  \n"
#define GENALG_TXTOMETER_LINE6 "Size pool xxxxxx                \n"
#define GENALG_TXTOMETER_FORMAT6 "Size pool %06d  \n"

// ------------- GenAlgAdn

// ================= Data structure ===================

typedef struct GenAlg GenAlg;

typedef struct GenAlgAdn {
  // ID
  unsigned long _id;
  // Age
  unsigned long _age;
  // Adn for floating point value
  VecFloat* _adnF;
  // Delta Adn during mutation for floating point value
  VecFloat* _deltaAdnF;
  // Adn for integer value
  VecLong* _adnI;
  // Value
  float _val;
  // Mutability of adn for floating point value
  VecFloat* _mutabilityF;
  // Mutability of adn for integer value
  VecFloat* _mutabilityI;
} GenAlgAdn;

// ================ Functions declaration ====================

// Create a new GenAlgAdn with ID 'id', 'lengthAdnF' and 'lengthAdnI'
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlgAdn* GenAlgAdnCreate(const unsigned long id, const long lengthAdnF, 
  const long lengthAdnI);

// Free memory used by the GenAlgAdn 'that'
void GenAlgAdnFree(GenAlgAdn** that);

// Return the adn for floating point values of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
const VecFloat* GAAdnAdnF(const GenAlgAdn* const that);

// Return the delta of adn for floating point values of the 
// GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
const VecFloat* GAAdnDeltaAdnF(const GenAlgAdn* const that);

// Return the adn for integer values of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
VecLong* GAAdnAdnI(const GenAlgAdn* const that);

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga' according to the type of the GenAlg
void GAAdnInit(const GenAlgAdn* const that, const GenAlg* ga);

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga', version used to calculate the parameters of a NeuraNet
void GAAdnInitNeuraNet(const GenAlgAdn* const that, const GenAlg* ga);

// Get the 'iGene'-th gene of the adn for floating point values of the
// GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
float GAAdnGetGeneF(const GenAlgAdn* const that, const long iGene);

// Get the delta of the 'iGene'-th gene of the adn for floating point 
// values of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
float GAAdnGetDeltaGeneF(const GenAlgAdn* const that, const long iGene);

// Get the 'iGene'-th gene of the adn for int values of the
// GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
int GAAdnGetGeneI(const GenAlgAdn* const that, const long iGene);

// Set the 'iGene'-th gene of the adn for floating point values of the
// GenAlgAdn 'that' to 'gene'
#if BUILDMODE != 0
inline
#endif
void GAAdnSetGeneF(GenAlgAdn* const that, const long iGene, 
  const float gene);

// Set the delta of the 'iGene'-th gene of the adn for floating point 
// values of the GenAlgAdn 'that' to 'delta'
#if BUILDMODE != 0
inline
#endif
void GAAdnSetDeltaGeneF(GenAlgAdn* const that, const long iGene, 
  const float delta);

// Set the 'iGene'-th gene of the adn for int values of the
// GenAlgAdn 'that'to 'gene'
#if BUILDMODE != 0
inline
#endif
void GAAdnSetGeneI(GenAlgAdn* const that, const long iGene, 
  const long gene);

// Get the id of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
unsigned long GAAdnGetId(const GenAlgAdn* const that);

// Get the age of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
unsigned long GAAdnGetAge(const GenAlgAdn* const that);

// Get the value of the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
float GAAdnGetVal(const GenAlgAdn* const that);

// Print the information about the GenAlgAdn 'that' on the 
// stream 'stream'
void GAAdnPrintln(const GenAlgAdn* const that, FILE* const stream);

// Return true if the GenAlgAdn 'that' is new, i.e. is age equals 1
// Return false
#if BUILDMODE != 0
inline
#endif
bool GAAdnIsNew(const GenAlgAdn* const that);

// Copy the GenAlgAdn 'tho' into the GenAlgAdn 'that'
#if BUILDMODE != 0
inline
#endif
void GAAdnCopy(GenAlgAdn* const that, const GenAlgAdn* const tho);

// Set the mutability vectors for the GenAlgAdn 'that' to 'mutability'
#if BUILDMODE != 0
inline
#endif
void GAAdnSetMutabilityInt(GenAlgAdn* const that, 
  const VecFloat* const mutability);
#if BUILDMODE != 0
inline
#endif
void GAAdnSetMutabilityFloat(GenAlgAdn* const that, 
  const VecFloat* const mutability);

// ------------- GenAlg

// ================= Data structure ===================

typedef enum GenAlgType {
  genAlgTypeDefault,
  genAlgTypeNeuraNet,
  genAlgTypeNeuraNetConv
} GenAlgType;

// Data used when GenAlg is applied to a NeuraNet
typedef struct GANeuraNet {
  // Nb of input, hidden and output of the NeuraNet
  int _nbIn;
  int _nbHid;
  int _nbOut;
  long _nbBaseConv;
  long _nbBaseCellConv;
  long _nbLink;
} GANeuraNet;

typedef struct GenAlg {
  // GSet of GenAlgAdn, sortval == score so the head of the set is the 
  // worst adn and the tail of the set is the best
  GSet* _adns;
  // Copy of the best adn
  GenAlgAdn* _bestAdn;
  // Type of the GenAlg
  GenAlgType _type;
  // Current epoch
  unsigned long _curEpoch;
  // Nb elite entities in population
  int _nbElites;
  // Id of the next new GenAlgAdn
  unsigned long _nextId;
  // Length of adn for floating point value
  const long _lengthAdnF;
  // Length of adn for integer value
  const long _lengthAdnI;
  // Bounds (min, max) for floating point values adn
  VecFloat2D* _boundsF;
  // Bounds (min, max) for integer values adn
  VecLong2D* _boundsI;
  // Norm of the range value for adns (optimization for diversity
  // calculation)
  float _normRangeFloat;
  float _normRangeInt;
  // Data used if the GenAlg is applied to a NeuraNet
  GANeuraNet _NNdata;
  // Number of ktevent
  unsigned long _nbKTEvent;
  // Flag to remember if we display info via a TextOMeter 
  // about the population
  bool _flagTextOMeter;
  // DiversityThreshold
  float _diversityThreshold;
  // TextOMeter to display information about the population
  // If the TextOMEter is used, its ocntent is refreshed at each call 
  // of the function GAStep();
  TextOMeter* _textOMeter;
  // Nb min of adns
  int _nbMinAdn;
  // Nb max of adns
  int _nbMaxAdn;
} GenAlg;

// ================ Functions declaration ====================

// Create a new GenAlg with 'nbEntities', 'nbElites', 'lengthAdnF' 
// and 'lengthAdnI'
// 'nbEntities' must greater than 2
// 'nbElites' must greater than 1
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlg* GenAlgCreate(const int nbEntities, const int nbElites, 
  const long lengthAdnF, const long lengthAdnI);

// Free memory used by the GenAlg 'that'
void GenAlgFree(GenAlg** that);

// Get the type of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
GenAlgType GAGetType(const GenAlg* const that);

// Set the type of the GenAlg 'that' to genAlgTypeNeuraNet, the GenAlg
// will be used with a NeuraNet having 'nbIn' inputs, 'nbHid' hidden 
// values and 'nbOut' outputs
#if BUILDMODE != 0
inline
#endif
void GASetTypeNeuraNet(GenAlg* const that, const int nbIn, 
  const int nbHid, const int nbOut);

// Set the type of the GenAlg 'that' to genAlgTypeNeuraNetConv, 
// the GenAlg will be used with a NeuraNet having 'nbIn' inputs, 
// 'nbHid' hidden values, 'nbOut' outputs, 'nbBaseConv' bases function,
// 'nbLink' links dedicated to the convolution and 'nbBaseCellConv' bases function per cell of convolution
#if BUILDMODE != 0
inline
#endif
void GASetTypeNeuraNetConv(GenAlg* const that, const int nbIn, 
  const int nbHid, const int nbOut, const long nbBaseConv,
  const long nbBaseCellConv, const long nbLink);

// Return the GSet of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
GSet* GAAdns(const GenAlg* const that);

// Return the max nb of adns of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetNbMaxAdn(const GenAlg* const that);

// Return the min nb of adns of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetNbMinAdn(const GenAlg* const that);

// Set the min nb of adns of the GenAlg 'that' to 'nb'
#if BUILDMODE != 0
inline
#endif
void GASetNbMaxAdn(GenAlg* const that, const int nb);

// Set the min nb of adns of the GenAlg 'that' to 'nb'
#if BUILDMODE != 0
inline
#endif
void GASetNbMinAdn(GenAlg* const that, const int nb);

// Return the nb of entities of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetNbAdns(const GenAlg* const that);

// Get the diversity threshold of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
float GAGetDiversityThreshold(const GenAlg* const that);

// Set the diversity threshold of the GenAlg 'that' to 'threshold'
#if BUILDMODE != 0
inline
#endif
void GASetDiversityThreshold(GenAlg* const that, const float threshold);

// Return the nb of elites of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
int GAGetNbElites(const GenAlg* const that);

// Return the current epoch of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
unsigned long GAGetCurEpoch(const GenAlg* const that);

// Return the number of KTEvent of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
unsigned long GAGetNbKTEvent(const GenAlg* const that);

// Set the nb of entities of the GenAlg 'that' to 'nb'
// 'nb' must be greater than 1, if 'nb' is lower than the current nb 
// of elite the number of elite is set to 'nb' - 1
void GASetNbEntities(GenAlg* const that, const int nb);

// Set the nb of elites of the GenAlg 'that' to 'nb'
// 'nb' must be greater than 0, if 'nb' is greater or equal to the 
// current nb of entities the number of entities is set to 'nb' + 1
void GASetNbElites(GenAlg* const that, const int nb);

// Get the length of adn for floating point value
#if BUILDMODE != 0
inline
#endif
long GAGetLengthAdnFloat(const GenAlg* const that);

// Get the length of adn for integer value
#if BUILDMODE != 0
inline
#endif
long GAGetLengthAdnInt(const GenAlg* const that);

// Get the bounds for the 'iGene'-th gene of adn for floating point 
// values
#if BUILDMODE != 0
inline
#endif
const VecFloat2D* GABoundsAdnFloat(const GenAlg* const that, 
  const long iGene);

// Get the bounds for the 'iGene'-th gene of adn for integer values
#if BUILDMODE != 0
inline
#endif
const VecLong2D* GABoundsAdnInt(const GenAlg* const that, 
  const long iGene);

// Set the bounds for the 'iGene'-th gene of adn for floating point 
// values to a copy of 'bounds'
#if BUILDMODE != 0
inline
#endif
void GASetBoundsAdnFloat(GenAlg* const that, const long iGene, 
  const VecFloat2D* const bounds);

// Set the bounds for the 'iGene'-th gene of adn for integer values
// to a copy of 'bounds'
#if BUILDMODE != 0
inline
#endif
void GASetBoundsAdnInt(GenAlg* const that, const long iGene, 
  const VecLong2D* bounds);

// Get the GenAlgAdn of the GenAlg 'that' currently at rank 'iRank'
#if BUILDMODE != 0
inline
#endif
GenAlgAdn* GAAdn(const GenAlg* const that, const int iRank);

// Init the GenAlg 'that'
// Must be called after the bounds have been set
// The random generator must have been initialised before calling this
// function
void GAInit(GenAlg* const that);

// Step an epoch for the GenAlg 'that' with the current ranking of
// GenAlgAdn
void GAStep(GenAlg* const that); 

// Print the information about the GenAlg 'that' on the stream 'stream'
void GAPrintln(const GenAlg* const that, FILE* const stream);

// Print a summary about the elite entities of the GenAlg 'that'
// on the stream 'stream'
void GAEliteSummaryPrintln(const GenAlg* const that, 
  FILE* const stream);

// Get the diversity of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
float GAGetDiversity(const GenAlg* const that);

// Function which return the JSON encoding of 'that' 
JSONNode* GAEncodeAsJSON(const GenAlg* const that);

// Function which decode from JSON encoding 'json' to 'that'
bool GADecodeAsJSON(GenAlg** that, const JSONNode* const json);

// Load the GenAlg 'that' from the stream 'stream'
// If the GenAlg is already allocated, it is freed before loading
// Return true in case of success, else false
bool GALoad(GenAlg** that, FILE* const stream);

// Save the GenAlg 'that' to the stream 'stream'
// If 'compact' equals true it saves in compact form, else it saves in 
// readable form
// Return true in case of success, else false
bool GASave(const GenAlg* const that, FILE* const stream, 
  const bool compact);

// Set the value of the GenAlgAdn 'adn' of the GenAlg 'that' to 'val'
#if BUILDMODE != 0
inline
#endif
void GASetAdnValue(GenAlg* const that, GenAlgAdn* const adn, 
  const float val);

// Update the norm of the range value for adans of the GenAlg 'that'
void GAUpdateNormRange(GenAlg* const that);

// Reset the GenAlg 'that'
// Randomize all the gene except those of the first adn
void GAKTEvent(GenAlg* const that);

// Return the best adn of the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
const GenAlgAdn* GABestAdn(const GenAlg* const that);

// Return the flag memorizing if the TextOMeter is displayed for
// the GenAlg 'that'
#if BUILDMODE != 0
inline
#endif
bool GAIsTextOMeterActive(const GenAlg* const that);

// Set the flag memorizing if the TextOMeter is displayed for
// the GenAlg 'that' to 'flag'
void GASetTextOMeterFlag(GenAlg* const that, bool flag);

// ================= Polymorphism ==================

// ================ Inliner ====================

#if BUILDMODE != 0
#include "genalg-inline.c"
#endif


#endif
