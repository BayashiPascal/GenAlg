// ============ GENALG.C ================

// ================= Include =================

#include "genalg.h"
#if BUILDMODE == 0
#include "genalg-inline.c"
#endif

// ------------- GenAlgAdn

// ================ Functions declaration ====================

// Get the diversity value of 'adnA' against 'adnB'
// The diversity is equal to 
float GAAdnGetDiversity(const GenAlgAdn* const adnA, 
  const GenAlgAdn* const adnB, const GenAlg* const ga);

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga'
void GAAdnInitDefault(const GenAlgAdn* const that, const GenAlg* ga);

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga', version used to calculate the parameters of a NeuraNet
void GAAdnInitNeuraNet(const GenAlgAdn* const that, const GenAlg* ga);

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga', version used to calculate the parameters of a NeuraNet
// with convolution
void GAAdnInitNeuraNetConv(const GenAlgAdn* const that, 
  const GenAlg* const ga);
    
// ================ Functions implementation ====================

// Create a new GenAlgAdn with ID 'id', 'lengthAdnF' and 'lengthAdnI'
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlgAdn* GenAlgAdnCreate(const unsigned long id, 
  const long lengthAdnF, const long lengthAdnI) {
#if BUILDMODE == 0
  if (lengthAdnF < 0) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'lengthAdnF' is invalid (%ld>=0)", 
      lengthAdnF);
    PBErrCatch(GenAlgErr);
  }
  if (lengthAdnI < 0) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'lengthAdnI' is invalid (%ld>=0)", 
      lengthAdnI);
    PBErrCatch(GenAlgErr);
  }
#endif
  // Allocate memory
  GenAlgAdn* that = PBErrMalloc(GenAlgErr, sizeof(GenAlgAdn));
  // Set the properties
  that->_age = 1;
  that->_id = id;
  that->_val = 0.0;
  if (lengthAdnF > 0) {
    that->_adnF = VecFloatCreate(lengthAdnF);
    that->_deltaAdnF = VecFloatCreate(lengthAdnF);
    that->_mutabilityF = VecFloatCreate(lengthAdnF);
  } else {
    that->_adnF = NULL;
    that->_deltaAdnF = NULL;
    that->_mutabilityF = NULL;
  }
  if (lengthAdnI > 0) {
    that->_adnI = VecLongCreate(lengthAdnI);
    that->_mutabilityI = VecFloatCreate(lengthAdnI);
  } else {
    that->_adnI = NULL;
    that->_mutabilityI = NULL;
  }
  // Return the new GenAlgAdn
  return that;
}

// Free memory used by the GenAlgAdn 'that'
void GenAlgAdnFree(GenAlgAdn** that) {
  // Check the argument
  if (that == NULL || *that == NULL) return;
  // Free memory
  if ((*that)->_adnF != NULL)
    VecFree(&((*that)->_adnF));
  if ((*that)->_deltaAdnF != NULL)
    VecFree(&((*that)->_deltaAdnF));
  if ((*that)->_adnI != NULL)
    VecFree(&((*that)->_adnI));
  if ((*that)->_mutabilityF != NULL)
    VecFree(&((*that)->_mutabilityF));
  if ((*that)->_mutabilityI != NULL)
    VecFree(&((*that)->_mutabilityI));
  free(*that);
  // Set the pointer to null
  *that = NULL;
}

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga' according to the type of GenAlg
void GAAdnInit(const GenAlgAdn* const that, const GenAlg* const ga) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  switch (GAGetType(ga)) {
    case genAlgTypeNeuraNet:
      GAAdnInitNeuraNet(that, ga);
      break;
    case genAlgTypeNeuraNetConv:
      GAAdnInitNeuraNetConv(that, ga);
      break;
    case genAlgTypeDefault:
    default:
      GAAdnInitDefault(that, ga);
  }
}

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga'
void GAAdnInitDefault(const GenAlgAdn* const that, 
  const GenAlg* const ga) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // For each floating point value gene
  for (long iGene = GAGetLengthAdnFloat(ga); iGene--;) {
    float min = VecGet(GABoundsAdnFloat(ga, iGene), 0);
    float max = VecGet(GABoundsAdnFloat(ga, iGene), 1);
    float val = min + (max - min) * rnd();
    VecSet(that->_adnF, iGene, val);
    VecSet(that->_mutabilityF, iGene, 1.0);
  }
  // For each integer value gene
  for (long iGene = GAGetLengthAdnInt(ga); iGene--;) {
    long min = VecGet(GABoundsAdnInt(ga, iGene), 0);
    long max = VecGet(GABoundsAdnInt(ga, iGene), 1);
    long val = (long)round((float)min + (float)(max - min) * rnd());
    VecSet(that->_adnI, iGene, val);
    VecSet(that->_mutabilityI, iGene, 1.0);
  }
}

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga', version used to calculate the parameters of a NeuraNet
// with convolution
void GAAdnInitNeuraNetConv(const GenAlgAdn* const that, 
  const GenAlg* const ga) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // For each floating point value gene
  for (long iGene = GAGetLengthAdnFloat(ga); iGene--;) {
    float min = VecGet(GABoundsAdnFloat(ga, iGene), 0);
    float max = VecGet(GABoundsAdnFloat(ga, iGene), 1);
    float val = min + (max - min) * rnd();
    VecSet(that->_adnF, iGene, val);
    VecSet(that->_mutabilityF, iGene, 1.0);
  }
}

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga', version used to calculate the parameters of a NeuraNet
void GAAdnInitNeuraNet(const GenAlgAdn* const that, const GenAlg* ga) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Init the base functions randomly
  // For each floating point value gene
  for (long iGene = GAGetLengthAdnFloat(ga); iGene--;) {
    float min = VecGet(GABoundsAdnFloat(ga, iGene), 0);
    float max = VecGet(GABoundsAdnFloat(ga, iGene), 1);
    float val = min + (max - min) * rnd();
    VecSet(that->_adnF, iGene, val);
    VecSet(that->_mutabilityF, iGene, 1.0);
  }
  // Init the links by ensuring there is at least one link reaching 
  // each output and use inputs as start of the initial links
  // For each integer value gene
  int shiftOut = ga->_NNdata._nbIn + ga->_NNdata._nbHid;
  for (long iGene = GAGetLengthAdnInt(ga); iGene--;) {
    VecSet(that->_adnI, iGene, -1);
    VecSet(that->_mutabilityI, iGene, 1.0);
  }
  for (int iOut = 0; iOut < ga->_NNdata._nbOut; ++iOut) {
    // The base function is randomly choosen but can't be an 
    // inactive link
    long min = 0;
    long max = VecGet(GABoundsAdnInt(ga, iOut * 3), 1);
    long val = (long)round((float)min + (float)(max - min) * rnd());
    VecSet(that->_adnI, iOut * 3, val);
    // The start of the link is randomly choosen amongst inputs
    min = 0;
    max = ga->_NNdata._nbIn - 1;
    val = (long)round((float)min + (float)(max - min) * rnd());
    VecSet(that->_adnI, iOut * 3 + 1, val);
    // The end of the link is choosen sequencially amongst outputs
    VecSet(that->_adnI, iOut * 3 + 2, iOut + shiftOut);
  }
}

// Print the information about the GenAlgAdn 'that' on the 
// stream 'stream'
void GAAdnPrintln(const GenAlgAdn* const that, FILE* const stream) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (stream == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'stream' is null");
    PBErrCatch(GenAlgErr);
  }
#endif  
  fprintf(stream, "id:%lu age:%lu", GAAdnGetId(that), GAAdnGetAge(that));
  fprintf(stream, "\n");
  fprintf(stream, "  adnF:");
  VecFloatPrint(GAAdnAdnF(that), stream,6);
  fprintf(stream, "\n");
  fprintf(stream, "  deltaAdnF:");
  VecFloatPrint(GAAdnDeltaAdnF(that), stream,6);
  fprintf(stream, "\n");
  fprintf(stream, "  adnI:");
  VecPrint(GAAdnAdnI(that), stream);
  fprintf(stream, "\n");
}

// ------------- GenAlg

// ================ Functions declaration ====================

// Select the rank of two parents for the SRM algorithm
// Return the ranks in 'parents', with parents[0] <= parents[1]
void GASelectParents(const GenAlg* const that, int* const parents);

// Set the genes of the entity at rank 'iChild' as a 50/50 mix of the 
// genes of entities at ranks 'parents[0]' and 'parents[1]'
void GAReproduction(GenAlg* const that, const int* const parents, 
  const int iChild);

// Set the genes of the entity at rank 'iChild' as a 50/50 mix of the 
// genes of entities at ranks 'parents[0]' and 'parents[1]'
void GAReproductionDefault(GenAlg* const that, 
  const int* const parents, const int iChild);

// Set the genes of the adn at rank 'iChild' as a 50/50 mix of the 
// genes of adns at ranks 'parents[0]' and 'parents[1]'
// This version is optimised to calculate the parameters of a NeuraNet
// with convolution by inheriting whole bases from parents
void GAReproductionNeuraNetConv(GenAlg* const that, 
  const int* const parents, const int iChild);
  
// Set the genes of the entity at rank 'iChild' as a 50/50 mix of the 
// genes of entities at ranks 'parents[0]' and 'parents[1]'
// This version is optimised to calculate the parameters of a NeuraNet
// by inheriting whole bases and links from parents
void GAReproductionNeuraNet(GenAlg* const that, 
  const int* const parents, const int iChild);

// Router toward the appropriate Mute function according to the type 
// of GenAlg
void GAMute(GenAlg* const that, const int* const parents, 
  const int iChild);
  
// Mute the genes of the entity at rank 'iChild'
void GAMuteDefault(GenAlg* const that, const int* const parents, 
  const int iChild);

// Mute the genes of the entity at rank 'iChild'
// This version is optimised to calculate the parameters of a NeuraNet
// by ensuring coherence in links: outputs have at least one link
// and there is no dead link
void GAMuteNeuraNet(GenAlg* const that, const int* const parents, 
  const int iChild);

// Mute the genes of the entity at rank 'iChild'
// This version is optimised to calculate the parameters of a NeuraNet
// with convolution by muting bases function per cell
void GAMuteNeuraNetConv(GenAlg* const that, const int* const parents, 
  const int iChild);
  
// ================ Functions implementation ====================

// Create a new GenAlg with 'nbEntities', 'nbElites', 'lengthAdnF' 
// and 'lengthAdnI'
// 'nbEntities' must greater than 2
// 'nbElites' must greater than 1
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlg* GenAlgCreate(const int nbEntities, const int nbElites, 
  const long lengthAdnF, const long lengthAdnI) {
  // Allocate memory
  GenAlg* that = PBErrMalloc(GenAlgErr, sizeof(GenAlg));
  // Set the properties
  that->_type = genAlgTypeDefault;
  that->_adns = GSetCreate();
  that->_curEpoch = 0;
  that->_nbKTEvent = 0;
  that->_bestAdn = GenAlgAdnCreate(0, lengthAdnF, lengthAdnI);
  *(long*)&(that->_lengthAdnF) = lengthAdnF;
  *(long*)&(that->_lengthAdnI) = lengthAdnI;
  if (lengthAdnF > 0) {
    that->_boundsF = 
      PBErrMalloc(GenAlgErr, sizeof(VecFloat2D) * lengthAdnF);
    for (long iGene = lengthAdnF; iGene--;)
      that->_boundsF[iGene] = VecFloatCreateStatic2D();
  } else
    that->_boundsF = NULL;
  if (lengthAdnI > 0) {
    that->_boundsI = 
      PBErrMalloc(GenAlgErr, sizeof(VecLong2D) * lengthAdnI);
    for (long iGene = lengthAdnI; iGene--;)
      that->_boundsI[iGene] = VecLongCreateStatic2D();
  } else
    that->_boundsI = NULL;
  that->_normRangeFloat = 1.0;
  that->_normRangeInt = 1.0;
  that->_nbElites = 0;
  that->_nextId = 0;
  GASetNbEntities(that, nbEntities);
  GASetNbElites(that, nbElites);
  // Return the new GenAlg
  return that;
}

// Free memory used by the GenAlg 'that'
void GenAlgFree(GenAlg** that) {
  // Check the argument
  if (that == NULL || *that == NULL) return;
  // Free memory
  GSetIterForward iter = GSetIterForwardCreateStatic(GAAdns(*that));
  do {
    GenAlgAdn* gaEnt = GSetIterGet(&iter);
    GenAlgAdnFree(&gaEnt);
  } while (GSetIterStep(&iter));
  GSetFree(&((*that)->_adns));
  if ((*that)->_boundsF != NULL)
    free((*that)->_boundsF);
  if ((*that)->_boundsI != NULL)
    free((*that)->_boundsI);
  GenAlgAdnFree(&((*that)->_bestAdn));
  free(*that);
  // Set the pointer to null
  *that = NULL;
}

// Set the nb of entities of the GenAlg 'that' to 'nb'
// 'nb' must be greater than 1, if 'nb' is lower than the current nb 
// of elite the number of elite is set to 'nb' - 1
void GASetNbEntities(GenAlg* const that, const int nb) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (nb <= 1) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'nb' is invalid (%d>1)", nb);
    PBErrCatch(GenAlgErr);
  }
#endif
  while (GSetNbElem(GAAdns(that)) > nb) {
    GenAlgAdn* gaEnt = GSetPop(GAAdns(that));
    GenAlgAdnFree(&gaEnt);
  }
  while (GSetNbElem(GAAdns(that)) < nb) {
    GenAlgAdn* ent = GenAlgAdnCreate(that->_nextId++,
      GAGetLengthAdnFloat(that), GAGetLengthAdnInt(that));
    GSetPush(GAAdns(that), ent);
  }
  if (GAGetNbElites(that) >= nb)
    GASetNbElites(that, nb - 1);
}

// Set the nb of elites of the GenAlg 'that' to 'nb'
// 'nb' must be greater than 0, if 'nb' is greater or equal to the 
// current nb of entities the number of entities is set to 'nb' + 1
void GASetNbElites(GenAlg* const that, const int nb) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (nb <= 1) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'nb' is invalid (%d>1)", nb);
    PBErrCatch(GenAlgErr);
  }
#endif
  if (GAGetNbAdns(that) <= nb)
    GASetNbEntities(that, nb + 1);
  that->_nbElites = nb;
}

// Init the GenAlg 'that'
// Must be called after the bounds have been set
// The random generator must have been initialised before calling this
// function
void GAInit(GenAlg* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // For each adn
  GSetIterForward iter = GSetIterForwardCreateStatic(GAAdns(that));
  do {
    // Get the adn
    GenAlgAdn* adn = GSetIterGet(&iter);
    // Initialise randomly the genes of the adn
    GAAdnInit(adn, that);
  } while (GSetIterStep(&iter));
  GAAdnCopy(that->_bestAdn, GAAdn(that, 0));
}

// Reset the GenAlg 'that'
// Randomize all the gene except those of the best adn
void GAKTEvent(GenAlg* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  ++(that->_nbKTEvent);
  GenAlgAdn* adn = GAAdn(that, 0);
  unsigned long int age = adn->_age;
  GAAdnCopy(adn, GABestAdn(that));
  adn->_age = age;
  int parents[2] = {0};
  GAMute(that, parents, 0);
  adn->_age = 1;
  for (int iEnt = 1; iEnt < GAGetNbAdns(that); 
    ++iEnt) {
    GenAlgAdn* adn = GAAdn(that, iEnt);
    GAAdnInit(adn, that);
    adn->_age = 1;
    adn->_id = (that->_nextId)++;
  }
}

// Step an epoch for the GenAlg 'that' with the current ranking of
// GenAlgAdn
void GAStep(GenAlg* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Selection, Reproduction, Mutation
  // Ensure the set of adns is sorted
  GSetSort(GAAdns(that));
  // Update the best adn if necessary
  if (that->_curEpoch == 1 || 
    GAAdnGetVal(GAAdn(that, 0)) > GAAdnGetVal(GABestAdn(that))) {
    GAAdnCopy(that->_bestAdn, GAAdn(that, 0));
  }
  // Declare a variable to memorize the parents
  int parents[2];
  // Get the diversity level
  float diversity = GAGetDiversity(that);
  // Correct the diversity level with the age of the best adn
  //diversity *= 
    //1.0 - fsquare((float)(GAAdnGetAge(GAAdn(that, 0))) / 1000.0);
  // If the diversity level is too low
  if (diversity < PBMATH_EPSILON || GAAdnGetAge(GAAdn(that, 0)) > 200) {
    // Renew diversity by applying a KT event (in memory of 
    // chickens' grand pa and grand ma)
    GAKTEvent(that);
  // Else, the diversity level is ok
  } else {
    // For each adn which is an elite
    for (int iAdn = 0; iAdn < GAGetNbElites(that); ++iAdn) {
      // Increment age
      ++(GAAdn(that, iAdn)->_age);
    }
    // For each adn which is not an elite
    for (int iAdn = GAGetNbElites(that); iAdn < GAGetNbAdns(that); 
      ++iAdn) {
      // Select two parents for this adn
      GASelectParents(that, parents);
      // Set the genes of the adn as a 50/50 mix of parents' genes
      GAReproduction(that, parents, iAdn);
      // Mute the genes of the adn
      GAMute(that, parents, iAdn);
    }
  }
  // Increment the number of epochs
  ++(that->_curEpoch);
}

// Select the rank of two parents for the SRM algorithm
// Return the ranks in 'parents', with parents[0] <= parents[1]
void GASelectParents(const GenAlg* const that, int* const parents) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Declare a variable to memorize the parents' rank
  int p[2];
  do {
    for (int i = 2; i--;)
      // p[i] below may be equal to the rank of the highest non elite 
      // adn, but it's not a problem so leave it and let's call that 
      // the Hawking radiation of this function in memory of this great 
      // man.
      p[i] = (int)floor(rnd() * (float)GAGetNbElites(that));
  } while (p[0] == p[1]);
  // Memorize the sorted parents' rank
  if (p[0] < p[1]) {
    parents[0] = p[0];
    parents[1] = p[1];
  } else {
    parents[0] = p[1];
    parents[1] = p[0];
  }
}

// Set the genes of the adn at rank 'iChild' as a 50/50 mix of the 
// genes of adns at ranks 'parents[0]' and 'parents[1]'
void GAReproduction(GenAlg* const that, 
  const int* const parents, const int iChild) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
  if (iChild < 0 || iChild >= GAGetNbAdns(that)) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, GAGetNbAdns(that));
    PBErrCatch(GenAlgErr);
  }
#endif
  switch (GAGetType(that)) {
    case genAlgTypeNeuraNet:
      GAReproductionNeuraNet(that, parents, iChild);
      break;
    case genAlgTypeNeuraNetConv:
      GAReproductionNeuraNetConv(that, parents, iChild);
      break;
    case genAlgTypeDefault:
    default:
      GAReproductionDefault(that, parents, iChild);
  }
}

// Set the genes of the adn at rank 'iChild' as a 50/50 mix of the 
// genes of adns at ranks 'parents[0]' and 'parents[1]'
// This version is optimised to calculate the parameters of a NeuraNet
// by inheriting whole bases and links from parents
void GAReproductionNeuraNet(GenAlg* const that, 
  const int* const parents, const int iChild) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
  if (iChild < 0 || iChild >= GAGetNbAdns(that)) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, GAGetNbAdns(that));
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the parents and child
  GenAlgAdn* parentA = GAAdn(that, parents[0]);
  GenAlgAdn* parentB = GAAdn(that, parents[1]);
  GenAlgAdn* child = GAAdn(that, iChild);
  // For each gene of the adn for floating point value
  for (long iGene = 0; iGene < GAGetLengthAdnFloat(that); iGene += 3) {
    // Get the gene from one parent or the other with equal 
    // probabililty
    if (rnd() < 0.5) {
      for (long jGene = 3; jGene--;) {
        VecSet(child->_adnF, iGene + jGene, 
          VecGet(parentA->_adnF, iGene + jGene));
        VecSet(child->_deltaAdnF, iGene + jGene, 
          VecGet(parentA->_deltaAdnF, iGene + jGene));
      }
    } else {
      for (long jGene = 3; jGene--;) {
        VecSet(child->_adnF, iGene + jGene, 
          VecGet(parentB->_adnF, iGene + jGene));
        VecSet(child->_deltaAdnF, iGene + jGene, 
          VecGet(parentB->_deltaAdnF, iGene + jGene));
      }
    }
  }
  // For each gene of the adn for int value
  for (long iGene = 0; iGene < GAGetLengthAdnInt(that); iGene += 3) {
    // Get the gene from one parent or the other with equal probabililty
    if (rnd() < 0.5) {
      for (long jGene = 3; jGene--;)
        VecSet(child->_adnI, iGene + jGene, 
          VecGet(parentA->_adnI, iGene + jGene));
    } else { 
      for (long jGene = 3; jGene--;)
        VecSet(child->_adnI, iGene + jGene, 
          VecGet(parentB->_adnI, iGene + jGene));
    }
  }
  // Reset the age of the child
  child->_age = 1;
  // Set the id of the child
  child->_id = (that->_nextId)++;
}

// Set the genes of the adn at rank 'iChild' as a 50/50 mix of the 
// genes of adns at ranks 'parents[0]' and 'parents[1]'
// This version is optimised to calculate the parameters of a NeuraNet
// with convolution by inheriting whole bases from parents
void GAReproductionNeuraNetConv(GenAlg* const that, 
  const int* const parents, const int iChild) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
  if (iChild < 0 || iChild >= GAGetNbAdns(that)) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, GAGetNbAdns(that));
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the parents and child
  GenAlgAdn* parentA = GAAdn(that, parents[0]);
  GenAlgAdn* parentB = GAAdn(that, parents[1]);
  GenAlgAdn* child = GAAdn(that, iChild);
  // For each gene of the adn for floating point value of convolution
  // base functions
  for (long iGene = 0; 
    iGene < that->_NNdata._nbBaseConv * 3; 
    iGene += that->_NNdata._nbBaseCellConv * 3) {
    // Get the gene from one parent or the other with equal probabililty
    if (rnd() < 0.5) {
      for (long jGene = that->_NNdata._nbBaseCellConv * 3;
        jGene--;) {
        VecSet(child->_adnF, iGene + jGene, 
          VecGet(parentA->_adnF, iGene + jGene));
        VecSet(child->_deltaAdnF, iGene + jGene, 
          VecGet(parentA->_deltaAdnF, iGene + jGene));
      }
    } else {
      for (long jGene = that->_NNdata._nbBaseCellConv * 3;
        jGene--;) {
        VecSet(child->_adnF, iGene + jGene, 
          VecGet(parentB->_adnF, iGene + jGene));
        VecSet(child->_deltaAdnF, iGene + jGene, 
          VecGet(parentB->_deltaAdnF, iGene + jGene));
      }
    }
  }
  // For each gene of the adn for floating point value of convolution
  // base functions
  for (long iGene = that->_NNdata._nbBaseConv * 3; 
    iGene < GAGetLengthAdnFloat(that); iGene += 3) {
    // Get the gene from one parent or the other with equal probabililty
    if (rnd() < 0.5) {
      for (long jGene = 3; --jGene;) {
        VecSet(child->_adnF, iGene + jGene, 
          VecGet(parentA->_adnF, iGene + jGene));
        VecSet(child->_deltaAdnF, iGene + jGene, 
          VecGet(parentA->_deltaAdnF, iGene + jGene));
      }
    } else {
      for (long jGene = 3; --jGene;) {
        VecSet(child->_adnF, iGene + jGene, 
          VecGet(parentB->_adnF, iGene + jGene));
        VecSet(child->_deltaAdnF, iGene + jGene, 
          VecGet(parentB->_deltaAdnF, iGene + jGene));
      }
    }
  }
  // Reset the age of the child
  child->_age = 1;
  // Set the id of the child
  child->_id = (that->_nextId)++;
}

// Set the genes of the adn at rank 'iChild' as a 50/50 mix of the 
// genes of adns at ranks 'parents[0]' and 'parents[1]'
void GAReproductionDefault(GenAlg* const that, 
  const int* const parents, const int iChild) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
  if (iChild < 0 || iChild >= GAGetNbAdns(that)) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, GAGetNbAdns(that));
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the parents and child
  GenAlgAdn* parentA = GAAdn(that, parents[0]);
  GenAlgAdn* parentB = GAAdn(that, parents[1]);
  GenAlgAdn* child = GAAdn(that, iChild);
  // For each gene of the adn for floating point value
  for (long iGene = GAGetLengthAdnFloat(that); iGene--;) {
    // Get the gene from one parent or the other with equal probabililty
    if (rnd() < 0.5) {
      VecSet(child->_adnF, iGene, VecGet(parentA->_adnF, iGene));
      VecSet(child->_deltaAdnF, iGene, 
        VecGet(parentA->_deltaAdnF, iGene));
    } else {
      VecSet(child->_adnF, iGene, VecGet(parentB->_adnF, iGene));
      VecSet(child->_deltaAdnF, iGene, 
        VecGet(parentB->_deltaAdnF, iGene));
    }
  }
  // For each gene of the adn for int value
  for (long iGene = GAGetLengthAdnInt(that); iGene--;) {
    // Get the gene from one parent or the other with equal probabililty
    if (rnd() < 0.5)
      VecSet(child->_adnI, iGene, VecGet(parentA->_adnI, iGene));
    else
      VecSet(child->_adnI, iGene, VecGet(parentB->_adnI, iGene));
  }
  // Reset the age of the child
  child->_age = 1;
  // Set the id of the child
  child->_id = (that->_nextId)++;
}

// Router toward the appropriate Mute function according to the type 
// of GenAlg
void GAMute(GenAlg* const that, const int* const parents, 
  const int iChild) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
  if (iChild < 0 || iChild >= GAGetNbAdns(that)) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, GAGetNbAdns(that));
    PBErrCatch(GenAlgErr);
  }
#endif
  switch (GAGetType(that)) {
    case genAlgTypeNeuraNet:
      GAMuteNeuraNet(that, parents, iChild);
      break;
    case genAlgTypeNeuraNetConv:
      GAMuteNeuraNetConv(that, parents, iChild);
      break;
    case genAlgTypeDefault:
    default:
      GAMuteDefault(that, parents, iChild);
  }
}

// Mute the genes of the entity at rank 'iChild'
// This version is optimised to calculate the parameters of a NeuraNet
// by ensuring coherence in links: outputs have at least one link
// and there is no dead link
void GAMuteNeuraNet(GenAlg* const that, const int* const parents, 
  const int iChild) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
  if (iChild < 0 || iChild >= GAGetNbAdns(that)) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, GAGetNbAdns(that));
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the first parent and child
  GenAlgAdn* parentA = GAAdn(that, parents[0]);
  GenAlgAdn* child = GAAdn(that, iChild);
  // Get the proba and amplitude of mutation
  float probMute = sqrt(((float)iChild) / ((float)GAGetNbAdns(that)));
  float amp = 1.0 - sqrt(1.0 / (float)(parentA->_age + 1));
  probMute /= (float)(GAGetLengthAdnInt(that));
  probMute += (float)(parentA->_age) / 10000.0;
  // Ensure the proba is not null
  if (probMute < PBMATH_EPSILON)
    probMute = PBMATH_EPSILON;
  // Declare a variable to memorize if there has been mutation
  bool hasMuted = false;
  // Declare a variable to memorize the used values amongst input and 
  // hidden
  long nbMaxUsedVal = that->_NNdata._nbIn + that->_NNdata._nbHid;
  char* isUsed = PBErrMalloc(GenAlgErr, sizeof(char) * nbMaxUsedVal);
  // Loop until there has been at least one mutation
  do {
    // Reset the used values
    memset(isUsed, 0, sizeof(char) * nbMaxUsedVal);
    memset(isUsed, 1, sizeof(char) * that->_NNdata._nbIn);
    // For each gene of the adn for int value (links definitions)
    for (long iGene = 0; iGene < GAGetLengthAdnInt(that); iGene += 3) {
      // If the link mutes
      if (rnd() < probMute) {
        hasMuted= true;
        // If this link is currently inactivated
        if (GAAdnGetGeneI(child, iGene) == -1) {
          // Base function
          long iBase = (int)round((float)iGene / 3.0);
          GAAdnSetGeneI(child, iGene, iBase);
          // Input
          long min = 
            VecGet(GABoundsAdnInt(that, iGene + 1), 0);
          long max = 
            VecGet(GABoundsAdnInt(that, iGene + 1), 1);
          long val = min;
          // Ensure the input is a used value
          do {
            val = (long)round((float)min + 
              (float)(max - min) * rnd());
          } while (isUsed[val] == 0);
          GAAdnSetGeneI(child, iGene + 1, val);
          // Output
          min = MAX(val, VecGet(GABoundsAdnInt(that, iGene + 2), 0));
          max = VecGet(GABoundsAdnInt(that, iGene + 2), 1);
          val = (long)round((float)min + (float)(max - min) * rnd());
          GAAdnSetGeneI(child, iGene + 2, val);
          if (val < nbMaxUsedVal)
            isUsed[val] = 1;
        // Else, this link is currently activated
        } else {
          // Choose between inactivation or mutation
          if (rnd() < 0.5) {
            // Inactivate the link
            GAAdnSetGeneI(child, iGene, -1);
          } else {
            // Input
            long min = 
              VecGet(GABoundsAdnInt(that, iGene + 1), 0);
            long max = 
              VecGet(GABoundsAdnInt(that, iGene + 1), 1);
            long val = min;
            // Ensure the input is a used value
            do {
              val = (long)round((float)min + 
                (float)(max - min) * rnd());
            } while (isUsed[val] == 0);
            GAAdnSetGeneI(child, iGene + 1, val);
            // Output
            min = MAX(val, VecGet(GABoundsAdnInt(that, iGene + 2), 0));
            max = VecGet(GABoundsAdnInt(that, iGene + 2), 1);
            val = (long)round((float)min + (float)(max - min) * rnd());
            GAAdnSetGeneI(child, iGene + 2, val);
            if (val < nbMaxUsedVal)
              isUsed[val] = 1;
          }
        }
      }
      // Get the index of the base function
      long baseFun = GAAdnGetGeneI(child, iGene);
      // If the link is active
      if (baseFun != -1) {
        // If the associated base function mutes
        if (rnd() < probMute) {
          long baseFunGene = baseFun * 3;
          for (long jGene = 3; jGene--;) {
            // Get the bounds
            const VecFloat2D* const bounds = 
              GABoundsAdnFloat(that, baseFunGene + jGene);
            // Declare a variable to memorize the previous value 
            // of the gene
            float prevVal = GAAdnGetGeneF(child, baseFunGene + jGene);
            // Apply the mutation
            GAAdnSetGeneF(child, baseFunGene + jGene, 
              GAAdnGetGeneF(child, baseFunGene + jGene) + 
              (VecGet(bounds, 1) - VecGet(bounds, 0)) * amp * 
              //VecGet(parentA->_mutabilityF, baseFun * 3) * 
              (rnd() - 0.5) + 
              GAAdnGetDeltaGeneF(child, baseFunGene + jGene));
            // Keep the gene value in bounds
            while (GAAdnGetGeneF(child, baseFunGene + jGene) < 
              VecGet(bounds, 0) ||
              GAAdnGetGeneF(child, baseFunGene + jGene) > 
              VecGet(bounds, 1)) {
              if (GAAdnGetGeneF(child, baseFunGene + jGene) > 
                VecGet(bounds, 1))
                GAAdnSetGeneF(child, baseFunGene + jGene, 
                  2.0 * VecGet(bounds, 1) - 
                  GAAdnGetGeneF(child, baseFunGene + jGene));
              else if (GAAdnGetGeneF(child, baseFunGene + jGene) < 
                VecGet(bounds, 0))
                GAAdnSetGeneF(child, baseFunGene + jGene, 
                  2.0 * VecGet(bounds, 0) - 
                  GAAdnGetGeneF(child, baseFunGene + jGene));
            }
            // Update the deltaAdn
            GAAdnSetDeltaGeneF(child, baseFunGene + jGene, 
              GAAdnGetGeneF(child, baseFunGene + jGene) - prevVal);
          }
        }
      }
    }
  } while (hasMuted == false);
  free(isUsed);
}

// Mute the genes of the entity at rank 'iChild'
void GAMuteDefault(GenAlg* const that, const int* const parents, 
  const int iChild) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
  if (iChild < 0 || iChild >= GAGetNbAdns(that)) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, GAGetNbAdns(that));
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the first parent and child
  GenAlgAdn* parentA = GAAdn(that, parents[0]);
  GenAlgAdn* child = GAAdn(that, iChild);
  // Get the proba amplitude of mutation
  float probMute = sqrt(((float)iChild) / ((float)GAGetNbAdns(that)));
  float amp = 1.0 - sqrt(1.0 / (float)(parentA->_age));
  probMute /= (float)(MAX(GAGetLengthAdnInt(that), 
    GAGetLengthAdnFloat(that)));
  probMute += (float)(parentA->_age) / 10000.0;
  if (probMute < PBMATH_EPSILON)
    probMute = PBMATH_EPSILON;
  bool hasMuted = false;
  do {
    // For each gene of the adn for floating point value
    for (long iGene = GAGetLengthAdnFloat(that); iGene--;) {
      // If this gene mutes
      if (rnd() < probMute) {
        hasMuted = true;
        // Get the bounds
        const VecFloat2D* const bounds = GABoundsAdnFloat(that, iGene);
        // Declare a variable to memorize the previous value of the gene
        float prevVal = GAAdnGetGeneF(child, iGene);
        // Apply the mutation
        GAAdnSetGeneF(child, iGene, GAAdnGetGeneF(child, iGene) + 
          (VecGet(bounds, 1) - VecGet(bounds, 0)) * amp * 
          (rnd() - 0.5) + GAAdnGetDeltaGeneF(child, iGene));
        // Keep the gene value in bounds
        while (GAAdnGetGeneF(child, iGene) < VecGet(bounds, 0) ||
          GAAdnGetGeneF(child, iGene) > VecGet(bounds, 1)) {
          if (GAAdnGetGeneF(child, iGene) > VecGet(bounds, 1))
            GAAdnSetGeneF(child, iGene, 
              2.0 * VecGet(bounds, 1) - GAAdnGetGeneF(child, iGene));
          else if (GAAdnGetGeneF(child, iGene) < VecGet(bounds, 0))
            GAAdnSetGeneF(child, iGene, 
              2.0 * VecGet(bounds, 0) - GAAdnGetGeneF(child, iGene));
        }
        // Update the deltaAdn
        GAAdnSetDeltaGeneF(child, iGene, 
          GAAdnGetGeneF(child, iGene) - prevVal);
      }
    }
    // For each gene of the adn for int value
    for (long iGene = GAGetLengthAdnInt(that); iGene--;) {
      // If this gene mutes
      if (rnd() < probMute) {
        hasMuted = true;
        // Get the bounds
        const VecLong2D* const boundsI = GABoundsAdnInt(that, iGene);
        VecFloat2D bounds = VecLongToFloat2D(boundsI);
        // Apply the mutation (as it is int value, ensure the amplitude
        // is big enough to have an effect
        float ampI = MIN(2.0, 
          (float)(VecGet(&bounds, 1) - VecGet(&bounds, 0)) * amp);
        GAAdnSetGeneI(child, iGene, GAAdnGetGeneI(child, iGene) +
          (long)round(ampI * (rnd() - 0.5)));
        // Keep the gene value in bounds
        while (GAAdnGetGeneI(child, iGene) < VecGet(&bounds, 0) ||
          GAAdnGetGeneI(child, iGene) > VecGet(&bounds, 1)) {
          if (GAAdnGetGeneI(child, iGene) > VecGet(&bounds, 1))
            GAAdnSetGeneI(child, iGene, 
              2 * VecGet(&bounds, 1) - GAAdnGetGeneI(child, iGene));
          else if (GAAdnGetGeneI(child, iGene) < VecGet(&bounds, 0))
            GAAdnSetGeneI(child, iGene, 
              2 * VecGet(&bounds, 0) - GAAdnGetGeneI(child, iGene));
        }
      }
    }
  } while (hasMuted == false);
}

// Mute the genes of the entity at rank 'iChild'
// This version is optimised to calculate the parameters of a NeuraNet
// with convolution by muting bases function per cell
void GAMuteNeuraNetConv(GenAlg* const that, const int* const parents, 
  const int iChild) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (parents == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'parents' is null");
    PBErrCatch(GenAlgErr);
  }
  if (iChild < 0 || iChild >= GAGetNbAdns(that)) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, GAGetNbAdns(that));
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the first parent and child
  GenAlgAdn* parentA = GAAdn(that, parents[0]);
  GenAlgAdn* child = GAAdn(that, iChild);
  // Get the proba amplitude of mutation
  float probMute = sqrt(((float)iChild) / ((float)GAGetNbAdns(that)));
  float amp = 1.0 - sqrt(1.0 / (float)(parentA->_age));
  probMute /= (float)(that->_NNdata._nbLink);
  probMute += (float)(parentA->_age) / 10000.0;
  if (probMute < PBMATH_EPSILON)
    probMute = PBMATH_EPSILON;
  bool hasMuted = false;
  int nbTry = 0;
  do {
    // For each gene of the adn for floating point value
    for (long iGene = GAGetLengthAdnFloat(that); iGene--;) {
      // If this gene mutes
      if (rnd() < probMute * VecGet(parentA->_mutabilityF, iGene)) {
        hasMuted = true;
        // Get the bounds
        const VecFloat2D* const bounds = GABoundsAdnFloat(that, iGene);
        // Declare a variable to memorize the previous value of the gene
        float prevVal = GAAdnGetGeneF(child, iGene);
        // Apply the mutation
        GAAdnSetGeneF(child, iGene, GAAdnGetGeneF(child, iGene) + 
          (VecGet(bounds, 1) - VecGet(bounds, 0)) * amp * 
          (rnd() - 0.5) + GAAdnGetDeltaGeneF(child, iGene));
        // Keep the gene value in bounds
        while (GAAdnGetGeneF(child, iGene) < VecGet(bounds, 0) ||
          GAAdnGetGeneF(child, iGene) > VecGet(bounds, 1)) {
          if (GAAdnGetGeneF(child, iGene) > VecGet(bounds, 1))
            GAAdnSetGeneF(child, iGene, 
              2.0 * VecGet(bounds, 1) - GAAdnGetGeneF(child, iGene));
          else if (GAAdnGetGeneF(child, iGene) < VecGet(bounds, 0))
            GAAdnSetGeneF(child, iGene, 
              2.0 * VecGet(bounds, 0) - GAAdnGetGeneF(child, iGene));
        }
        // Update the deltaAdn
        GAAdnSetDeltaGeneF(child, iGene, 
          GAAdnGetGeneF(child, iGene) - prevVal);
      }
    }
    ++nbTry;
  } while (hasMuted == false && nbTry < 10);
}

// Print the information about the GenAlg 'that' on the stream 'stream'
void GAPrintln(const GenAlg* const that, FILE* const stream) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (stream == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'stream' is null");
    PBErrCatch(GenAlgErr);
  }
#endif  
  fprintf(stream, "epoch:%lu\n", GAGetCurEpoch(that));
  fprintf(stream, "%d entities, %d elites\n", GAGetNbAdns(that), 
    GAGetNbElites(that));
  GSetIterBackward iter = GSetIterBackwardCreateStatic(GAAdns(that));
  int iEnt = 0;
  do {
    GenAlgAdn* ent = GSetIterGet(&iter);
    fprintf(stream, "#%d value:%f ", iEnt,
      GSetIterGetElem(&iter)->_sortVal);
    if (iEnt < GAGetNbElites(that))
      fprintf(stream, "elite ");
    GAAdnPrintln(ent, stream);
    ++iEnt;
  } while (GSetIterStep(&iter));
}

// Print a summary about the elite entities of the GenAlg 'that'
// on the stream 'stream'
void GAEliteSummaryPrintln(const GenAlg* const that, 
  FILE* const stream) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (stream == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'stream' is null");
    PBErrCatch(GenAlgErr);
  }
#endif  
  GSetIterBackward iter = GSetIterBackwardCreateStatic(GAAdns(that));
  int iEnt = 0;
  GenAlgAdn* leader = GSetIterGet(&iter);
  fprintf(stream, "(age,val,div) ");
  do {
    GenAlgAdn* ent = GSetIterGet(&iter);
    fprintf(stream, "(%lu,%.3f,%.3f) ", GAAdnGetAge(ent), 
      GSetIterGetElem(&iter)->_sortVal, 
      GAAdnGetDiversity(ent, leader, that));
    ++iEnt;
  } while (GSetIterStep(&iter) && iEnt < GAGetNbElites(that));
  fprintf(stream, "\n");
}

// Update the norm of the range value for adans of the GenAlg 'that'
void GAUpdateNormRange(GenAlg* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // If there are float adn
  if (GAGetLengthAdnFloat(that) > 0) {
    // Declare a vector to memorize the ranges in float gene values
    VecFloat* range = VecFloatCreate(GAGetLengthAdnFloat(that)); 
    // Calculate the ranges in gene values
    for (long iGene = GAGetLengthAdnFloat(that); iGene--;)
      VecSet(range, iGene, 
        VecGet(GABoundsAdnFloat(that, iGene), 1) - 
        VecGet(GABoundsAdnFloat(that, iGene), 0));
    // Calculate the norm of the range
    that->_normRangeFloat = VecNorm(range);
    // Free memory
    VecFree(&range);
  }

  // If there are int adn
  if (GAGetLengthAdnInt(that) > 0) {
    // Declare a vector to memorize the ranges in int gene values
    VecFloat* range = VecFloatCreate(GAGetLengthAdnInt(that)); 
    // Calculate the ranges in gene values
    for (long iGene = GAGetLengthAdnInt(that); iGene--;)
      VecSet(range, iGene, 
        VecGet(GABoundsAdnInt(that, iGene), 1) - 
        VecGet(GABoundsAdnInt(that, iGene), 0));
    // Calculate the norm of the range
    that->_normRangeInt = VecNorm(range);
    // Free memory
    VecFree(&range);
  }
}


// Get the diversity value of 'adnA' against 'adnB'
// The diversity is equal to 
float GAAdnGetDiversity(const GenAlgAdn* const adnA, 
  const GenAlgAdn* const adnB, const GenAlg* const ga) {
#if BUILDMODE == 0
  if (adnA == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'adnA' is null");
    PBErrCatch(GenAlgErr);
  }
  if (adnB == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'adnB' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Declare a variable to memorize the result
  float diversity = 0.0;
  // If there are adn for floating point values
  if (GAAdnAdnF(adnA) != NULL && GAAdnAdnF(adnB) != NULL) {
    // Get the difference in adn with the first entity
    VecFloat* diff = 
      VecGetOp(GAAdnAdnF(adnA), 1.0, GAAdnAdnF(adnB), -1.0);
    // Calculate the diversity
    diversity += VecNorm(diff) / ga->_normRangeFloat;
    // Free memory
    VecFree(&diff);
  }
  // If there are adn for int values
  if (GAAdnAdnI(adnA) != NULL && GAAdnAdnI(adnB) != NULL) {
    // Get the difference in adn with the first entity
    VecLong* diffI = 
      VecGetOp(GAAdnAdnI(adnA), 1, GAAdnAdnI(adnB), -1);
    VecFloat* diff = VecLongToFloat(diffI);
    // Calculate the diversity
    diversity += VecNorm(diff) / ga->_normRangeInt;
    // Free memory
    VecFree(&diffI);
    VecFree(&diff);
  }
  // Correct diversity if there was both float and int adns
  if (GAAdnAdnF(adnA) != NULL && GAAdnAdnF(adnB) != NULL && 
    GAAdnAdnI(adnA) != NULL && GAAdnAdnI(adnB) != NULL)
    diversity /= 2.0;
  // Return the result
  return diversity;
}

// Function which return the JSON encoding of 'that' 
JSONNode* GAAdnEncodeAsJSON(const GenAlgAdn* const that, 
  const float elo) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBMathErr->_type = PBErrTypeNullPointer;
    sprintf(PBMathErr->_msg, "'that' is null");
    PBErrCatch(PBMathErr);
  }
#endif
  // Create the JSON structure
  JSONNode* json = JSONCreate();
  // Declare a buffer to convert value into string
  char val[100];
  // Encode the id
  sprintf(val, "%lu", that->_id);
  JSONAddProp(json, "_id", val);
  // Encode the age
  sprintf(val, "%lu", that->_age);
  JSONAddProp(json, "_age", val);
  // Encode the elo
  sprintf(val, "%f", elo);
  JSONAddProp(json, "_elo", val);
  // Encode the value
  sprintf(val, "%f", that->_val);
  JSONAddProp(json, "_val", val);
  // Encode the genes
  if (that->_adnF != NULL) {
    JSONAddProp(json, "_adnF", VecEncodeAsJSON(that->_adnF));
    JSONAddProp(json, "_deltaAdnF", VecEncodeAsJSON(that->_deltaAdnF));
  }
  if (that->_adnI != NULL)
    JSONAddProp(json, "_adnI", VecEncodeAsJSON(that->_adnI));
  // Return the created JSON 
  return json;
}

// Function which return the JSON encoding of 'that' 
JSONNode* GAEncodeAsJSON(const GenAlg* const that) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBMathErr->_type = PBErrTypeNullPointer;
    sprintf(PBMathErr->_msg, "'that' is null");
    PBErrCatch(PBMathErr);
  }
#endif
  // Create the JSON structure
  JSONNode* json = JSONCreate();
  // Declare a buffer to convert value into string
  char val[100];
  // Encode the type
  sprintf(val, "%d", GAGetType(that));
  JSONAddProp(json, "_type", val);
  switch (GAGetType(that)) {
    case genAlgTypeNeuraNet:
      sprintf(val, "%d", that->_NNdata._nbIn);
      JSONAddProp(json, "NN_nbIn", val);
      sprintf(val, "%d", that->_NNdata._nbHid);
      JSONAddProp(json, "NN_nbHid", val);
      sprintf(val, "%d", that->_NNdata._nbOut);
      JSONAddProp(json, "NN_nbOut", val);
      break;
    case genAlgTypeNeuraNetConv:
      sprintf(val, "%d", that->_NNdata._nbIn);
      JSONAddProp(json, "NN_nbIn", val);
      sprintf(val, "%d", that->_NNdata._nbHid);
      JSONAddProp(json, "NN_nbHid", val);
      sprintf(val, "%d", that->_NNdata._nbOut);
      JSONAddProp(json, "NN_nbOut", val);
      sprintf(val, "%ld", that->_NNdata._nbBaseConv);
      JSONAddProp(json, "NN_nbBaseConv", val);
      sprintf(val, "%ld", that->_NNdata._nbBaseCellConv);
      JSONAddProp(json, "NN_nbBaseCellConv", val);
      sprintf(val, "%ld", that->_NNdata._nbLink);
      JSONAddProp(json, "NN_nbLink", val);
      break;
    default:
      break;
  }
  // Encode the nb adns
  sprintf(val, "%d", GAGetNbAdns(that));
  JSONAddProp(json, "_nbAdns", val);
  // Encode the nb elites
  sprintf(val, "%d", GAGetNbElites(that));
  JSONAddProp(json, "_nbElites", val);
  // Encode the length adn float
  sprintf(val, "%ld", GAGetLengthAdnFloat(that));
  JSONAddProp(json, "_lengthAdnF", val);
  // Encode the length adn int
  sprintf(val, "%ld", GAGetLengthAdnInt(that));
  JSONAddProp(json, "_lengthAdnI", val);
  // Encode the epoch
  sprintf(val, "%lu", GAGetCurEpoch(that));
  JSONAddProp(json, "_curEpoch", val);
  // Encode the next id
  sprintf(val, "%lu", that->_nextId);
  JSONAddProp(json, "_nextId", val);
  // Encode the bounds
  JSONArrayStruct setBoundFloat = JSONArrayStructCreateStatic();
  if (GAGetLengthAdnFloat(that) > 0) {
    for (long iBound = 0; iBound < GAGetLengthAdnFloat(that); ++iBound)
      JSONArrayStructAdd(&setBoundFloat, 
        VecEncodeAsJSON((VecFloat*)GABoundsAdnFloat(that, iBound)));
    JSONAddProp(json, "_boundFloat", &setBoundFloat);
  }
  JSONArrayStruct setBoundInt = JSONArrayStructCreateStatic();
  if (GAGetLengthAdnInt(that) > 0) {
    for (long iBound = 0; iBound < GAGetLengthAdnInt(that); ++iBound)
      JSONArrayStructAdd(&setBoundInt, 
        VecEncodeAsJSON((VecLong*)GABoundsAdnInt(that, iBound)));
    JSONAddProp(json, "_boundInt", &setBoundInt);
  }
  // Save the adns
  JSONArrayStruct setAdn = JSONArrayStructCreateStatic();
  for (int iEnt = 0; iEnt < GAGetNbAdns(that); ++iEnt) {
    GenAlgAdn* ent = GSetElemData(GSetElement(GAAdns(that), iEnt));
    float sortVal = GSetElemGetSortVal(GSetElement(GAAdns(that), iEnt));
    JSONArrayStructAdd(&setAdn, GAAdnEncodeAsJSON(ent, sortVal));
  }
  JSONAddProp(json, "_adns", &setAdn);
  // Save the best adn
  JSONAddProp(json, "_bestAdn", 
    GAAdnEncodeAsJSON(GABestAdn(that), 0.0));
  // Free memory
  JSONArrayStructFlush(&setBoundFloat);
  JSONArrayStructFlush(&setBoundInt);
  JSONArrayStructFlush(&setAdn);
  // Return the created JSON 
  return json;
}

// Function which decode from JSON encoding 'json' to 'that'
bool GAAdnDecodeAsJSON(GenAlgAdn** that, const JSONNode* const json) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBMathErr->_type = PBErrTypeNullPointer;
    sprintf(PBMathErr->_msg, "'that' is null");
    PBErrCatch(PBMathErr);
  }
  if (json == NULL) {
    PBMathErr->_type = PBErrTypeNullPointer;
    sprintf(PBMathErr->_msg, "'json' is null");
    PBErrCatch(PBMathErr);
  }
#endif
  // If 'that' is already allocated
  if (*that != NULL)
    // Free memory
    GenAlgAdnFree(that);
  // Get the id from the JSON
  JSONNode* prop = JSONProperty(json, "_id");
  if (prop == NULL) {
    return false;
  }
  unsigned long id = strtoul(JSONLabel(JSONValue(prop, 0)), NULL, 10);
  // Get the lengthAdnF from the JSON
  long lengthAdnF = 0;
  prop = JSONProperty(json, "_adnF");
  if (prop != NULL) {
    JSONNode* subprop = JSONProperty(prop, "_dim");
    lengthAdnF = atol(JSONLabel(JSONValue(subprop, 0)));
  }
  // Get the lengthAdnI from the JSON
  long lengthAdnI = 0;
  prop = JSONProperty(json, "_adnI");
  if (prop != NULL) {
    JSONNode* subprop = JSONProperty(prop, "_dim");
    lengthAdnI = atol(JSONLabel(JSONValue(subprop, 0)));
  }
  // Allocate memory
  *that = GenAlgAdnCreate(id, lengthAdnF, lengthAdnI);
  // Get the age from the JSON
  prop = JSONProperty(json, "_age");
  if (prop == NULL) {
    return false;
  }
  (*that)->_age = strtoul(JSONLabel(JSONValue(prop, 0)), NULL, 10);
  // Get the adnF from the JSON
  prop = JSONProperty(json, "_adnF");
  if (prop != NULL) {
    if (!VecDecodeAsJSON(&((*that)->_adnF), prop)) {
      return false;
    }
    prop = JSONProperty(json, "_deltaAdnF");
    if (prop == NULL) {
      return false;
    }
    if (!VecDecodeAsJSON(&((*that)->_deltaAdnF), prop)) {
      return false;
    }
  }
  // Get the adnI from the JSON
  prop = JSONProperty(json, "_adnI");
  if (prop != NULL)
    if (!VecDecodeAsJSON(&((*that)->_adnI), prop)) {
      return false;
    }
  // Get the value
  prop = JSONProperty(json, "_val");
  if (prop == NULL) {
    return false;
  }
  (*that)->_val = atof(JSONLabel(JSONValue(prop, 0)));
  // Return the success code
  return true;
}

// Function which decode from JSON encoding 'json' to 'that'
bool GADecodeAsJSON(GenAlg** that, const JSONNode* const json) {
#if BUILDMODE == 0
  if (that == NULL) {
    PBMathErr->_type = PBErrTypeNullPointer;
    sprintf(PBMathErr->_msg, "'that' is null");
    PBErrCatch(PBMathErr);
  }
  if (json == NULL) {
    PBMathErr->_type = PBErrTypeNullPointer;
    sprintf(PBMathErr->_msg, "'json' is null");
    PBErrCatch(PBMathErr);
  }
#endif
  // If 'that' is already allocated
  if (*that != NULL)
    // Free memory
    GenAlgFree(that);
  // Decode the nb adns
  JSONNode* prop = JSONProperty(json, "_nbAdns");
  if (prop == NULL) {
    return false;
  }
  int nbAdns = atoi(JSONLabel(JSONValue(prop, 0)));
  // Decode the nb elites
  prop = JSONProperty(json, "_nbElites");
  if (prop == NULL) {
    return false;
  }
  int nbElites = atoi(JSONLabel(JSONValue(prop, 0)));
  // Decode the length adn float
  prop = JSONProperty(json, "_lengthAdnF");
  if (prop == NULL) {
    return false;
  }
  long lengthAdnF = atol(JSONLabel(JSONValue(prop, 0)));
  // Decode the length adn int
  prop = JSONProperty(json, "_lengthAdnI");
  if (prop == NULL) {
    return false;
  }
  long lengthAdnI = atol(JSONLabel(JSONValue(prop, 0)));
  // Allocate memory
  *that = GenAlgCreate(nbAdns, nbElites, lengthAdnF, lengthAdnI);
  // Decode the type
  prop = JSONProperty(json, "_type");
  if (prop == NULL) {
    return false;
  }
  int type = atoi(JSONLabel(JSONValue(prop, 0)));
  int nbIn = 0;
  int nbOut = 0;
  int nbHid = 0;
  switch (type) {
    case genAlgTypeNeuraNet:
      prop = JSONProperty(json, "NN_nbIn");
      if (prop == NULL) {
        return false;
      }
      nbIn = atoi(JSONLabel(JSONValue(prop, 0)));
      prop = JSONProperty(json, "NN_nbOut");
      if (prop == NULL) {
        return false;
      }
      nbOut = atoi(JSONLabel(JSONValue(prop, 0)));
      prop = JSONProperty(json, "NN_nbHid");
      if (prop == NULL) {
        return false;
      }
      nbHid = atoi(JSONLabel(JSONValue(prop, 0)));
      prop = JSONProperty(json, "NN_nbBaseConv");
      if (prop == NULL) {
        return false;
      }
      GASetTypeNeuraNet(*that, nbIn, nbHid, nbOut);
      break;
    case genAlgTypeNeuraNetConv:
      prop = JSONProperty(json, "NN_nbIn");
      if (prop == NULL) {
        return false;
      }
      nbIn = atoi(JSONLabel(JSONValue(prop, 0)));
      prop = JSONProperty(json, "NN_nbOut");
      if (prop == NULL) {
        return false;
      }
      nbOut = atoi(JSONLabel(JSONValue(prop, 0)));
      prop = JSONProperty(json, "NN_nbHid");
      if (prop == NULL) {
        return false;
      }
      nbHid = atoi(JSONLabel(JSONValue(prop, 0)));
      prop = JSONProperty(json, "NN_nbBaseConv");
      if (prop == NULL) {
        return false;
      }
      long nbBaseConv = atol(JSONLabel(JSONValue(prop, 0)));
      prop = JSONProperty(json, "NN_nbBaseCellConv");
      if (prop == NULL) {
        return false;
      }
      long nbBaseCellConv = atol(JSONLabel(JSONValue(prop, 0)));
      prop = JSONProperty(json, "NN_nbLink");
      if (prop == NULL) {
        return false;
      }
      long nbLink = atol(JSONLabel(JSONValue(prop, 0)));
      GASetTypeNeuraNetConv(*that, nbIn, nbHid, nbOut, nbBaseConv, 
        nbBaseCellConv, nbLink);
      break;
    default:
      break;
  }
  // Decode the epoch
  prop = JSONProperty(json, "_curEpoch");
  if (prop == NULL) {
    return false;
  }
  (*that)->_curEpoch = 
    strtoul(JSONLabel(JSONValue(prop, 0)), NULL, 10);
  // Decode the next id
  prop = JSONProperty(json, "_nextId");
  if (prop == NULL) {
    return false;
  }
  (*that)->_nextId = strtoul(JSONLabel(JSONValue(prop, 0)), NULL, 10);
  // Decode the bounds
  prop = JSONProperty(json, "_boundFloat");
  if (prop != NULL) {
    if (JSONGetNbValue(prop) != GAGetLengthAdnFloat(*that))
      return false;
    for (long iBound = 0; iBound < GAGetLengthAdnFloat(*that); ++iBound) {
      JSONNode* val = JSONValue(prop, iBound);
      VecFloat2D* b = NULL;
      if (!VecDecodeAsJSON((VecFloat**)&b, val)) {
        return false;
      }
      GASetBoundsAdnFloat(*that, iBound, b);
      VecFree((VecFloat**)&b);
    }
  }
  prop = JSONProperty(json, "_boundInt");
  if (prop != NULL) {
    if (JSONGetNbValue(prop) != GAGetLengthAdnInt(*that))
      return false;
    for (long iBound = 0; iBound < GAGetLengthAdnInt(*that); ++iBound) {
      JSONNode* val = JSONValue(prop, iBound);
      VecLong2D* b = NULL;
      if (!VecDecodeAsJSON((VecLong**)&b, val)) {
        return false;
      }
      GASetBoundsAdnInt(*that, iBound, b);
      VecFree((VecLong**)&b);
    }
  }
  // Upadte the norm of the range values
  GAUpdateNormRange(*that);
  // Decode the adns
  prop = JSONProperty(json, "_adns");
  if (prop == NULL) {
    return false;
  }
  if (JSONGetNbValue(prop) != GAGetNbAdns(*that))
    return false;
  for (int iEnt = 0; iEnt < GAGetNbAdns(*that); ++iEnt) {
    JSONNode* val = JSONValue(prop, iEnt);
    if (!GAAdnDecodeAsJSON(
      (GenAlgAdn**)&(GSetElement(GAAdns(*that), iEnt)->_data), val)) {
      return false;
    }
  }
  // Decode the best adn
  prop = JSONProperty(json, "_bestAdn");
  if (prop == NULL) {
    return false;
  }
  if (!GAAdnDecodeAsJSON((GenAlgAdn**)&((*that)->_bestAdn), prop)) {
    return false;
  }
  
  // Return the success code
  return true;
}

// Load the GenAlg 'that' from the stream 'stream'
// If the GenAlg is already allocated, it is freed before loading
// Return true in case of success, else false
bool GALoad(GenAlg** that, FILE* const stream) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (stream == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'stream' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Declare a json to load the encoded data
  JSONNode* json = JSONCreate();
  // Load the whole encoded data
  if (!JSONLoad(json, stream)) {
    return false;
  }
  // Decode the data from the JSON
  if (!GADecodeAsJSON(that, json)) {
    return false;
  }
  // Free the memory used by the JSON
  JSONFree(&json);
  // Return the success code
  return true;
}

// Save the GenAlg 'that' to the stream 'stream'
// If 'compact' equals true it saves in compact form, else it saves in 
// readable form
// Return true in case of success, else false
bool GASave(const GenAlg* const that, FILE* const stream, 
  const bool compact) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (stream == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'stream' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the JSON encoding
  JSONNode* json = GAEncodeAsJSON(that);
  // Save the JSON
  if (!JSONSave(json, stream, compact)) {
    return false;
  }
  // Free memory
  JSONFree(&json);
  // Return success code
  return true;
}

