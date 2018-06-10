// ============ GENALG.C ================

// ================= Include =================

#include "genalg.h"
#if BUILDMODE == 0
#include "genalg-inline.c"
#endif

// ------------- GenAlgAdn

// ================ Functions declaration ====================

// ================ Functions implementation ====================

// Create a new GenAlgAdn with ID 'id', 'lengthAdnF' and 'lengthAdnI'
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlgAdn* GenAlgAdnCreate(int id, int lengthAdnF, 
  int lengthAdnI) {
#if BUILDMODE == 0
  if (lengthAdnF < 0) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'lengthAdnF' is invalid (%d>=0)", 
      lengthAdnF);
    PBErrCatch(GenAlgErr);
  }
  if (lengthAdnI < 0) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'lengthAdnI' is invalid (%d>=0)", 
      lengthAdnI);
    PBErrCatch(GenAlgErr);
  }
#endif
  // Allocate memory
  GenAlgAdn* that = PBErrMalloc(GenAlgErr, sizeof(GenAlgAdn));
  // Set the properties
  that->_age = 1;
  that->_id = id;
  if (lengthAdnF > 0) {
    that->_adnF = VecFloatCreate(lengthAdnF);
    that->_deltaAdnF = VecFloatCreate(lengthAdnF);
  } else {
    that->_adnF = NULL;
    that->_deltaAdnF = NULL;
  }
  if (lengthAdnI > 0)
    that->_adnI = VecShortCreate(lengthAdnI);
  else
    that->_adnI = NULL;
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
  free(*that);
  // Set the pointer to null
  *that = NULL;
}

// Initialise randomly the genes of the GenAlgAdn 'that' of the 
// GenAlg 'ga'
void GAAdnInit(GenAlgAdn* that, GenAlg* ga) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // For each floating point value gene
  for (int iGene = GAGetLengthAdnFloat(ga); iGene--;) {
    float min = VecGet(GABoundsAdnFloat(ga, iGene), 0);
    float max = VecGet(GABoundsAdnFloat(ga, iGene), 1);
    float val = min + (max - min) * rnd();
    VecSet(that->_adnF, iGene, val);
  }
  // For each integer value gene
  for (int iGene = GAGetLengthAdnInt(ga); iGene--;) {
    short min = VecGet(GABoundsAdnInt(ga, iGene), 0);
    short max = VecGet(GABoundsAdnInt(ga, iGene), 1);
    short val = (short)round((float)min + (float)(max - min) * rnd());
    VecSet(that->_adnI, iGene, val);
  }
}

// Print the information about the GenAlgAdn 'that' on the 
// stream 'stream'
void GAAdnPrintln(GenAlgAdn* that, FILE* stream) {
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
void GASelectParents(GenAlg* that, int* parents);

// Set the genes of the entity at rank 'iChild' as a 50/50 mix of the 
// genes of entities at ranks 'parents[0]' and 'parents[1]'
void GAReproduction(GenAlg* that, int* parents, int iChild);

// Mute the genes of the entity at rank 'iChild'
// The probability of mutation for one gene is equal to 
// 'rankChild'/'that'->_nbEntities
// The amplitude of the mutation
// is equal to (max-min).(gauss(0.0, 1.0)+deltaAdn).ln('parents[0]'.age)
void GAMute(GenAlg* that, int* parents, int iChild);

// Reset the GenAlg 'that'
// Randomize all the gene except those of the first adn
void GAKTEvent(GenAlg* that);

// ================ Functions implementation ====================

// Create a new GenAlg with 'nbEntities', 'nbElites', 'lengthAdnF' 
// and 'lengthAdnI'
// 'nbEntities' must greater than 2
// 'nbElites' must greater than 1
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlg* GenAlgCreate(int nbEntities, int nbElites, int lengthAdnF, 
  int lengthAdnI) {
  // Allocate memory
  GenAlg* that = PBErrMalloc(GenAlgErr, sizeof(GenAlg));
  // Set the properties
  that->_adns = GSetCreate();
  that->_curEpoch = 0;
  that->_lengthAdnF = lengthAdnF;
  that->_lengthAdnI = lengthAdnI;
  if (lengthAdnF > 0) {
    that->_boundsF = 
      PBErrMalloc(GenAlgErr, sizeof(VecFloat2D) * lengthAdnF);
    for (int iGene = lengthAdnF; iGene--;)
      that->_boundsF[iGene] = VecFloatCreateStatic2D();
  } else
    that->_boundsF = NULL;
  if (lengthAdnI > 0) {
    that->_boundsI = 
      PBErrMalloc(GenAlgErr, sizeof(VecShort2D) * lengthAdnI);
    for (int iGene = lengthAdnI; iGene--;)
      that->_boundsI[iGene] = VecShortCreateStatic2D();
  } else
    that->_boundsI = NULL;
  that->_nbElites = 0;
  that->_nextId = 0;
  that->_diversityThreshold = GENALG_DIVERSITYTHRESHOLD;
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
  free(*that);
  // Set the pointer to null
  *that = NULL;
}

// Set the nb of entities of the GenAlg 'that' to 'nb'
// 'nb' must be greater than 1, if 'nb' is lower than the current nb 
// of elite the number of elite is set to 'nb' - 1
void GASetNbEntities(GenAlg* that, int nb) {
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
void GASetNbElites(GenAlg* that, int nb) {
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
void GAInit(GenAlg* that) {
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
}

// Reset the GenAlg 'that'
// Randomize all the gene except those of the best adn
void GAKTEvent(GenAlg* that) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // For each adn except the best one
  GSetIterBackward iter = GSetIterBackwardCreateStatic(GAAdns(that));
  GSetIterStep(&iter);
  do {
    // Get the adn
    GenAlgAdn* adn = GSetIterGet(&iter);
    // Initialise randomly the genes of the adn
    GAAdnInit(adn, that);
    // Reset the age of the child
    adn->_age = 1;
    // Set the id of the child
    adn->_id = (that->_nextId)++;
  } while (GSetIterStep(&iter));
}

// Step an epoch for the GenAlg 'that' with the current ranking of
// GenAlgAdn
void GAStep(GenAlg* that) {
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
  // Declare a variable to memorize the parents
  int parents[2];
  // Get the diversity level
  float diversity = GAGetDiversity(that);
  // If the diversity level is too low
  if (diversity < GENALG_DIVERSITYTHRESHOLD) {
    // Break the diversity by applying KT event (in memory of 
    // chickens' grand pa and grand ma)
    GAKTEvent(that);
  // Else, the diversity level is ok
  } else {
    // For each adn which is an elite
    for (int iAdn = 0; iAdn < GAGetNbElites(that); ++iAdn) {
      // Increment age
      (GAAdn(that, iAdn)->_age)++;
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
void GASelectParents(GenAlg* that, int* parents) {
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
  for (int i = 2; i--;)
    // p[i] below may be equal to the rank of the highest non elite 
    // adn, but it's not a problem so leave it and let's call that 
    // the Hawking radiation of this function in memory of this great 
    // man.
    p[i] = (int)floor(rnd() * (float)GAGetNbElites(that));
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
void GAReproduction(GenAlg* that, int* parents, int iChild) {
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
  for (int iGene = GAGetLengthAdnFloat(that); iGene--;) {
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
  for (int iGene = GAGetLengthAdnInt(that); iGene--;) {
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

// Mute the genes of the entity at rank 'iChild'
// The probability of mutation for one gene is equal to 
// 'rankChild'/'that'->_nbEntities
// The amplitude of the mutation
// is equal to (max-min).(gauss(0.0, 1.0)+deltaAdn).ln('parents[0]'.age)
void GAMute(GenAlg* that, int* parents, int iChild) {
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
  float probMute = ((float)iChild) / ((float)GAGetNbAdns(that));
  float amp = 1.0 - 1.0 / sqrt((float)(parentA->_age + 1));
  // For each gene of the adn for floating point value
  for (int iGene = GAGetLengthAdnFloat(that); iGene--;) {
    // If this gene mutes
    if (rnd() < probMute) {
      // Get the bounds
      VecFloat2D* bounds = GABoundsAdnFloat(that, iGene);
      // Declare a variable to memorize the previous value of the gene
      float prevVal = GAAdnGetGeneF(child, iGene);
      // Apply the mutation
      GAAdnSetGeneF(child, iGene, GAAdnGetGeneF(child, iGene) + 
        (VecGet(bounds, 1) - VecGet(bounds, 0)) * amp * 
        (rnd() - 0.5 + GAAdnGetDeltaGeneF(child, iGene)));
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
  for (int iGene = GAGetLengthAdnInt(that); iGene--;) {
    // If this gene mutes
    if (rnd() < probMute) {
      // Get the bounds
      VecShort2D* boundsI = GABoundsAdnInt(that, iGene);
      VecFloat2D bounds = VecShortToFloat2D(boundsI);
      // Apply the mutation (as it is int value, ensure the amplitude
      // is big enough to have an effect
      float ampI = MIN(2.0, 
        (float)(VecGet(&bounds, 1) - VecGet(&bounds, 0)) * amp);
      GAAdnSetGeneI(child, iGene, GAAdnGetGeneI(child, iGene) +
        (short)round(ampI * (rnd() - 0.5)));
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
}

// Print the information about the GenAlg 'that' on the stream 'stream'
void GAPrintln(GenAlg* that, FILE* stream) {
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

// Get the level of diversity of curent entities of the GenAlg 'that'
// The return value is in [0.0, 1.0]
// 0.0 means all the elite entities have exactly the same adns 
float GAGetDiversity(GenAlg* that) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Declare a variable to memorize the result
  float diversity = 0.0;
  // Declare a variable for calculation
  int nb = 1;
  // If there are adn for floating point values
  if (GAGetLengthAdnFloat(that) > 0) {
    // Declare a vector to memorize the ranges in gene values
    VecFloat* range = VecFloatCreate(GAGetLengthAdnFloat(that)); 
    // Calculate the ranges in gene values
    for (int iGene = GAGetLengthAdnFloat(that); iGene--;)
      VecSet(range, iGene, 
        VecGet(GABoundsAdnFloat(that, iGene), 1) - 
        VecGet(GABoundsAdnFloat(that, iGene), 0));
    // Calculate the norm of the range
    float normRange = VecNorm(range);
    // For each elite entity except the first one
    for (int iEnt = 1; iEnt < GAGetNbElites(that); ++iEnt) {
      // Get the difference in adn with the first entity
      VecFloat* diff = VecGetOp(GAAdnAdnF(GAAdn(that, iEnt)), 1.0, 
        GAAdnAdnF(GAAdn(that, 0)), -1.0);
      // Calculate the diversity
      diversity += VecNorm(diff) / 
        (normRange * (float)GAAdnGetAge(GAAdn(that, iEnt)));
      // Free memory
      VecFree(&diff);
    }
    // Calculate the diversity
    nb += GAGetNbElites(that);
    // Free memory
    VecFree(&range);
  }
  // If there are adn for floating point values
  if (GAGetLengthAdnInt(that) > 0) {
    // Declare a vector to memorize the ranges in gene values
    VecFloat* range = VecFloatCreate(GAGetLengthAdnInt(that));
    // Calculate the ranges in gene values
    for (int iGene = GAGetLengthAdnInt(that); iGene--;)
      VecSet(range, iGene, 
        (float)(VecGet(GABoundsAdnInt(that, iGene), 1) - 
        VecGet(GABoundsAdnInt(that, iGene), 0)));
    // Calculate the norm of the range
    float normRange = VecNorm(range);
    // For each elite entity except the first one
    for (int iEnt = 1; iEnt < GAGetNbElites(that); ++iEnt) {
      // Get the difference in adn with the first entity
      VecShort* diff = VecGetOp(GAAdnAdnI(GAAdn(that, iEnt)), 1, 
        GAAdnAdnI(GAAdn(that, 0)), -1);
      VecFloat* diffF = VecShortToFloat(diff);
      // Calculate the diversity
      diversity += VecNorm(diffF) / 
        (normRange * (float)GAAdnGetAge(GAAdn(that, iEnt)));
      // Free memory
      VecFree(&diffF);
      VecFree(&diff);
    }
    // Calculate the diversity
    nb += GAGetNbElites(that);
    // Free memory
    VecFree(&range);
  }
  // Calculate the diversity
  diversity /= (float)nb;
  // Return the result
  return diversity;
}

// Function which return the JSON encoding of 'that' 
JSONNode* GAAdnEncodeAsJSON(GenAlgAdn* that, float elo) {
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
JSONNode* GAEncodeAsJSON(GenAlg* that) {
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
  // Encode the diversity threshold
  sprintf(val, "%f", GAGetDiversityThreshold(that));
  JSONAddProp(json, "_diversityThreshold", val);
  // Encode the nb adns
  sprintf(val, "%d", GAGetNbAdns(that));
  JSONAddProp(json, "_nbAdns", val);
  // Encode the nb elites
  sprintf(val, "%d", GAGetNbElites(that));
  JSONAddProp(json, "_nbElites", val);
  // Encode the length adn float
  sprintf(val, "%d", GAGetLengthAdnFloat(that));
  JSONAddProp(json, "_lengthAdnF", val);
  // Encode the length adn int
  sprintf(val, "%d", GAGetLengthAdnInt(that));
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
    for (int iBound = 0; iBound < GAGetLengthAdnFloat(that); ++iBound)
      JSONArrayStructAdd(&setBoundFloat, 
        VecEncodeAsJSON((VecFloat*)GABoundsAdnFloat(that, iBound)));
    JSONAddProp(json, "_boundFloat", &setBoundFloat);
  }
  JSONArrayStruct setBoundInt = JSONArrayStructCreateStatic();
  if (GAGetLengthAdnInt(that) > 0) {
    for (int iBound = 0; iBound < GAGetLengthAdnInt(that); ++iBound)
      JSONArrayStructAdd(&setBoundInt, 
        VecEncodeAsJSON((VecShort*)GABoundsAdnInt(that, iBound)));
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
  // Free memory
  JSONArrayStructFlush(&setBoundFloat);
  JSONArrayStructFlush(&setBoundInt);
  JSONArrayStructFlush(&setAdn);
  // Return the created JSON 
  return json;
}

// Function which decode from JSON encoding 'json' to 'that'
bool GAAdnDecodeAsJSON(GenAlgAdn** that, JSONNode* json) {
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
  int id = strtoul(JSONLabel(JSONValue(prop, 0)), NULL, 10);
  // Get the lengthAdnF from the JSON
  int lengthAdnF = 0;
  prop = JSONProperty(json, "_adnF");
  if (prop != NULL) {
    JSONNode* subprop = JSONProperty(prop, "_dim");
    lengthAdnF = atoi(JSONLabel(JSONValue(subprop, 0)));
  }
  // Get the lengthAdnF from the JSON
  int lengthAdnI = 0;
  prop = JSONProperty(json, "_adnI");
  if (prop != NULL) {
    JSONNode* subprop = JSONProperty(prop, "_dim");
    lengthAdnI = atoi(JSONLabel(JSONValue(subprop, 0)));
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
    if (JSONGetNbValue(prop) != lengthAdnF)
      return false;
    if (!VecDecodeAsJSON(&((*that)->_adnF), prop))
      return false;
    prop = JSONProperty(json, "_deltaAdnF");
    if (prop == NULL)
      return false;
    if (!VecDecodeAsJSON(&((*that)->_deltaAdnF), prop))
      return false;
  }
  // Get the adnI from the JSON
  prop = JSONProperty(json, "_adnI");
  if (prop != NULL)
    if (JSONGetNbValue(prop) != lengthAdnI)
      return false;
    if (!VecDecodeAsJSON(&((*that)->_adnI), prop))
      return false;
  // Return the success code
  return true;
}

// Function which decode from JSON encoding 'json' to 'that'
bool GADecodeAsJSON(GenAlg** that, JSONNode* json) {
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
  int lengthAdnF = atoi(JSONLabel(JSONValue(prop, 0)));
  // Decode the length adn int
  prop = JSONProperty(json, "_lengthAdnI");
  if (prop == NULL) {
    return false;
  }
  int lengthAdnI = atoi(JSONLabel(JSONValue(prop, 0)));
  // Allocate memory
  *that = GenAlgCreate(nbAdns, nbElites, lengthAdnF, lengthAdnI);
  // Decode the diversity threshold
  prop = JSONProperty(json, "_diversityThreshold");
  if (prop == NULL) {
    return false;
  }
  (*that)->_diversityThreshold = atof(JSONLabel(JSONValue(prop, 0)));
  // Decode the epoch
  prop = JSONProperty(json, "_curEpoch");
  if (prop == NULL) {
    return false;
  }
  (*that)->_curEpoch = strtoul(JSONLabel(JSONValue(prop, 0)), NULL, 10);
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
    for (int iBound = 0; iBound < GAGetLengthAdnFloat(*that); ++iBound) {
      JSONNode* val = JSONValue(prop, iBound);
      VecFloat2D* b = NULL;
      if (!VecDecodeAsJSON((VecFloat**)&b, val))
        return false;
      GASetBoundsAdnFloat(*that, iBound, b);
      VecFree((VecFloat**)&b);
    }
  }
  prop = JSONProperty(json, "_boundInt");
  if (prop != NULL) {
    if (JSONGetNbValue(prop) != GAGetLengthAdnInt(*that))
      return false;
    for (int iBound = 0; iBound < GAGetLengthAdnInt(*that); ++iBound) {
      JSONNode* val = JSONValue(prop, iBound);
      VecShort2D* b = NULL;
      if (!VecDecodeAsJSON((VecShort**)&b, val))
        return false;
      GASetBoundsAdnInt(*that, iBound, b);
      VecFree((VecShort**)&b);
    }
  }
  // Decode the adns
  prop = JSONProperty(json, "_adns");
  if (prop == NULL) {
    return false;
  }
  if (JSONGetNbValue(prop) != GAGetNbAdns(*that))
    return false;
  for (int iEnt = 0; iEnt < GAGetNbAdns(*that); ++iEnt) {
    JSONNode* val = JSONValue(prop, iEnt);
    GenAlgAdn* data = GSetElemData(GSetElement(GAAdns(*that), iEnt));
    if (!GAAdnDecodeAsJSON(&data, val))
      return false;
  }
  // Return the success code
  return true;
}

// Load the GenAlg 'that' from the stream 'stream'
// If the GenAlg is already allocated, it is freed before loading
// Return true in case of success, else false
bool GALoad(GenAlg** that, FILE* stream) {
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
bool GASave(GenAlg* that, FILE* stream, bool compact) {
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

