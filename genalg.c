// ============ GENALG.C ================

// ================= Include =================

#include "genalg.h"
#if BUILDMODE == 0
#include "genalg-inline.c"
#endif

// ------------- GenAlgEntity

// ================ Functions declaration ====================

// Select the rank of two parents for the SRM algorithm
// Return the ranks in 'parents', with parents[0] <= parents[1]
void GASelectParents(GenAlg* that, int* parents);

// Set the genes of the entity at rank 'iChild' as a 50/50 mix of the 
// genes of entities at ranks 'parents[0]' and 'parents[1]'
void GAReproduction(GenAlg* that, int* parents, int iChild);

// Mute the genes of the entity at rank 'iChild' 
// The probability of mutation for one gene is equal to 
// 'iChild'/'that'->_nbEntities and the amplitude of the mutation
// is equal to (max-min).(gauss(0.0, 1.0)+deltaAdn).ln('parents[0]'.age)
void GAMute(GenAlg* that, int* parents, int iChild);

// ================ Functions implementation ====================

// Create a new GenAlgEntity with ID 'id', 'lengthAdnF' and 'lengthAdnI'
// 'lengthAdnF' and 'lengthAdnI' must be greater than or equal to 0
GenAlgEntity* GenAlgEntityCreate(int id, int lengthAdnF, 
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
  GenAlgEntity* that = PBErrMalloc(GenAlgErr, sizeof(GenAlgEntity));
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
  // Return the new GenAlgEntity
  return that;
}

// Free memory used by the GenAlgEntity 'that'
void GenAlgEntityFree(GenAlgEntity** that) {
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

// Initialise randomly the genes of the GenAlgEntity 'that' of the 
// GenAlg 'ga'
void GAEntInit(GenAlgEntity* that, GenAlg* ga) {
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

// Print the information about the GenAlgEntity 'that' on the 
// stream 'stream'
void GAEntPrintln(GenAlgEntity* that, FILE* stream) {
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
  fprintf(stream, "id:%d age:%d", GAEntGetId(that), GAEntGetAge(that));
  fprintf(stream, "\n");
  fprintf(stream, "  adnF:");
  VecPrint(GAEntAdnF(that), stream);
  fprintf(stream, "\n");
  fprintf(stream, "  deltaAdnF:");
  VecPrint(GAEntDeltaAdnF(that), stream);
  fprintf(stream, "\n");
  fprintf(stream, "  adnI:");
  VecPrint(GAEntAdnI(that), stream);
  fprintf(stream, "\n");
}

// ------------- GenAlg

// ================ Functions declaration ====================

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
  that->_elo = ELORankCreate();
  that->_runsPerEpoch = GENALG_RUNPEREPOCH;
  that->_curRun = 0;
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
  that->_nbEntities = 0;
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
  for (int iEnt = (*that)->_nbEntities; iEnt--;) {
    GenAlgEntity* gaEnt = GAEntity(*that, iEnt);
    GenAlgEntityFree(&gaEnt);
  }
  ELORankFree(&((*that)->_elo));
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
  while (that->_nbEntities > nb) {
    ELOEntity* ent = ELORankGetRanked(GAEloRank(that), 
      that->_nbEntities - 1);
    GenAlgEntity* gaEnt = (GenAlgEntity*)(ent->_data);
    ELORankRemove(GAEloRank(that), ent->_data);
    GenAlgEntityFree(&gaEnt);
    that->_nbEntities = ELORankGetNb(GAEloRank(that));
  }
  while (that->_nbEntities < nb) {
    GenAlgEntity* ent = GenAlgEntityCreate(that->_nextId++,
      GAGetLengthAdnFloat(that), GAGetLengthAdnInt(that));
    ELORankAdd(GAEloRank(that), ent);
    that->_nbEntities = ELORankGetNb(GAEloRank(that));
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
  if (GAGetNbEntities(that) <= nb)
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
  // For each entity
  GSetIterForward iter = 
    GSetIterForwardCreateStatic(&(GAEloRank(that)->_set));
  do {
    // Get the entity
    GenAlgEntity* ent = ((ELOEntity*)GSetIterGet(&iter))->_data;
    // Initialise randomly the genes of the entity
    GAEntInit(ent, that);
  } while (GSetIterStep(&iter));
}

// Step a run for the GenAlg 'that' with ranking of GenAlgEntity given
// in the GSet of GenAlgEntity 'rank'  (from best to worst, ie _sortVal
// from greater to lower)
void GAStepRun(GenAlg* that, GSet* rank) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
  if (rank == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'rank' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Update the ELORank
  ELORankUpdate(GAEloRank(that), rank);
  // Increment the number of runs
  ++(that->_curRun);
  // Increment the age of the entities of this run
  GSetIterForward iter = GSetIterForwardCreateStatic(rank);
  do {
    GenAlgEntity* ent = GSetIterGet(&iter);
    ++(ent->_age);
  } while (GSetIterStep(&iter));
  // If we have reached the end of the current epoch
  if (that->_curRun >= that->_runsPerEpoch)
    // Step the epoch
    GAStepEpoch(that);
}

// Step an epoch for the GenAlg 'that' with the current ranking of
// GenAlgEntity
void GAStepEpoch(GenAlg* that) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Selection, Reproduction, Mutation
  // Declare a GSet of GenAlgEntity to memorize the childs
  GSet* childs = GSetCreate();
  // Declare a variable to memorize the parents
  int parents[2];
  // Get the inbreeding level
  float inbreeding = GAGetInbreeding(that);
  // If the inbreeding level is too high
  // Break the inbreeding by applying mutation to elites except the 
  // best one
  if (inbreeding < GENALG_INBREEDINGTHRESHOLD) {
    // For each entity which is not an elite
    for (int iEnt = 1; iEnt < GAGetNbElites(that); ++iEnt) {
      // Reproduce with itself and mute the genes of the entity
      parents[1] = parents[0] = iEnt;
      GAReproduction(that, parents, iEnt);
      GAMute(that, parents, iEnt);
      // Add this entity to the set of childs too
      GSetAppend(childs, GAEntity(that, iEnt));
    }
  }
  // For each entity which is not an elite
  for (int iEnt = GAGetNbElites(that); iEnt < GAGetNbEntities(that); 
    ++iEnt) {
    // Select two parents for this entity
    GASelectParents(that, parents);
    // Set the genes of the entity as a 50/50 mix of parents' genes
    GAReproduction(that, parents, iEnt);
    // Mute the genes of the entity
    GAMute(that, parents, iEnt);
    // Add the child to the set of childs
    GSetAppend(childs, GAEntity(that, iEnt));
  }
  // Remove and re-add the childs from/to the ELOrank to reset their
  // position and ELO
  GSetIterForward iter = GSetIterForwardCreateStatic(childs);
  do {
    ELORankRemove(GAEloRank(that), GSetIterGet(&iter));
    ELORankAdd(GAEloRank(that), GSetIterGet(&iter));
  } while (GSetIterStep(&iter));
  // Increment the number of epochs
  ++(that->_curEpoch);
  // Reset the current run
  that->_curRun = 0;
  // Free memory
  GSetFree(&childs);
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
    // entity, but it's not a problem so leave it and let's call that 
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

// Set the genes of the entity at rank 'iChild' as a 50/50 mix of the 
// genes of entities at ranks 'parents[0]' and 'parents[1]'
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
  if (iChild < 0 || iChild >= that->_nbEntities) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, that->_nbEntities);
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the parents and child
  GenAlgEntity* parentA = GAEntity(that, parents[0]);
  GenAlgEntity* parentB = GAEntity(that, parents[1]);
  GenAlgEntity* child = GAEntity(that, iChild);
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
// 'iChild'/'that'->_nbEntities and the amplitude of the mutation
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
  if (iChild < 0 || iChild >= that->_nbEntities) {
    GenAlgErr->_type = PBErrTypeInvalidArg;
    sprintf(GenAlgErr->_msg, "'child' is invalid (0<=%d<%d)",
      iChild, that->_nbEntities);
    PBErrCatch(GenAlgErr);
  }
#endif
  // Get the first parent and child
  GenAlgEntity* parentA = GAEntity(that, parents[0]);
  GenAlgEntity* child = GAEntity(that, iChild);
  // Get the probability of mutation
  float probMute = ((float)iChild) / ((float)GAGetNbEntities(that));
  // Get the amplitude of mutation
  float amp = 1.0 - 1.0 / sqrt((float)(parentA->_age + 1));
  // For each gene of the adn for floating point value
  for (int iGene = GAGetLengthAdnFloat(that); iGene--;) {
    // If this gene mutes
    if (rnd() < probMute) {
      // Get the bounds
      VecFloat2D* bounds = GABoundsAdnFloat(that, iGene);
      // Declare a variable to memorize the previous value of the gene
      float prevVal = GAEntGetGeneF(child, iGene);
      // Apply the mutation
      GAEntSetGeneF(child, iGene, GAEntGetGeneF(child, iGene) + 
        (VecGet(bounds, 1) - VecGet(bounds, 0)) * amp * 
        (rnd() - 0.5 + GAEntGetDeltaGeneF(child, iGene)));
      // Keep the gene value in bounds
      while (GAEntGetGeneF(child, iGene) < VecGet(bounds, 0) ||
        GAEntGetGeneF(child, iGene) > VecGet(bounds, 1)) {
        if (GAEntGetGeneF(child, iGene) > VecGet(bounds, 1))
          GAEntSetGeneF(child, iGene, 
            2.0 * VecGet(bounds, 1) - GAEntGetGeneF(child, iGene));
        else if (GAEntGetGeneF(child, iGene) < VecGet(bounds, 0))
          GAEntSetGeneF(child, iGene, 
            2.0 * VecGet(bounds, 0) - GAEntGetGeneF(child, iGene));
      }
      // Update the deltaAdn
      GAEntSetDeltaGeneF(child, iGene, 
        GAEntGetGeneF(child, iGene) - prevVal);
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
      GAEntSetGeneI(child, iGene, GAEntGetGeneI(child, iGene) +
        (short)round(ampI * (rnd() - 0.5)));
      // Keep the gene value in bounds
      while (GAEntGetGeneI(child, iGene) < VecGet(&bounds, 0) ||
        GAEntGetGeneI(child, iGene) > VecGet(&bounds, 1)) {
        if (GAEntGetGeneI(child, iGene) > VecGet(&bounds, 1))
          GAEntSetGeneI(child, iGene, 
            2 * VecGet(&bounds, 1) - GAEntGetGeneI(child, iGene));
        else if (GAEntGetGeneI(child, iGene) < VecGet(&bounds, 0))
          GAEntSetGeneI(child, iGene, 
            2 * VecGet(&bounds, 0) - GAEntGetGeneI(child, iGene));
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
  fprintf(stream, "epoch:%d - run:%d\n", 
    GAGetCurEpoch(that), GAGetCurRun(that));
  fprintf(stream, "%d entities, %d elites\n", GAGetNbEntities(that), 
    GAGetNbElites(that));
  for (int iEnt = 0; iEnt < GAGetNbEntities(that); ++iEnt) {
    GenAlgEntity* ent = GAEntity(that, iEnt);
    fprintf(stream, "#%d elo:%f ", iEnt,
      ELORankGetELO(GAEloRank(that), ent));
    if (iEnt < GAGetNbElites(that))
      fprintf(stream, "elite ");
    GAEntPrintln(ent, stream);
  }
}

// Get the level of inbreeding of curent entities of the GenAlg 'that'
// The return value is in [0.0, 1.0]
// 0.0 means all the elite entities have exactly the same adns 
float GAGetInbreeding(GenAlg* that) {
#if BUILDMODE == 0
  if (that == NULL) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "'that' is null");
    PBErrCatch(GenAlgErr);
  }
#endif
  // Declare a variable to memorize the result
  float inbreeding = 0.0;
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
      VecFloat* diff = VecGetOp(GAEntAdnF(GAEntity(that, iEnt)), 1.0, 
        GAEntAdnF(GAEntity(that, 0)), -1.0);
      // Calculate the inbreeding
      inbreeding += VecNorm(diff) / normRange;
      // Free memory
      VecFree(&diff);
    }
    // Calculate the inbreeding
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
      VecShort* diff = VecGetOp(GAEntAdnI(GAEntity(that, iEnt)), 1, 
        GAEntAdnI(GAEntity(that, 0)), -1);
      VecFloat* diffF = VecShortToFloat(diff);
      // Calculate the inbreeding
      inbreeding += VecNorm(diffF) / normRange;
      // Free memory
      VecFree(&diffF);
      VecFree(&diff);
    }
    // Calculate the inbreeding
    nb += GAGetNbElites(that);
    // Free memory
    VecFree(&range);
  }
  // Calculate the inbreeding
  inbreeding /= (float)nb;
  // Return the result
  return inbreeding;
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
  // If 'that' is already allocated
  if (*that != NULL) {
    // Free memory
    GenAlgFree(that);
  }
  // Load the number of entity and elite, and the length of adn
  int nbEnt, nbElite, lenAdnF, lenAdnI;
  int ret = fscanf(stream, "%d %d %d %d", &nbEnt, &nbElite, 
    &lenAdnF, &lenAdnI);
  // If we couldn't fscanf
  if (ret == EOF)
    return false;
  // Check the data
  if (nbEnt < 3 || nbElite < 2 || lenAdnF < 0 || lenAdnI < 0)
    return false;
  // Allocate memory
  *that = GenAlgCreate(nbEnt, nbElite, lenAdnF, lenAdnI);
  // Load the run, epoch, nbRunPerEpoch, nextId
  ret = fscanf(stream, "%d %d %d %d", &((*that)->_curRun), 
    &((*that)->_curEpoch), &((*that)->_runsPerEpoch),
    &((*that)->_nextId));
  // If we couldn't fscanf
  if (ret == EOF)
    return false;
  // Load the bounds
  for (int iBound = 0; iBound < lenAdnF; ++iBound) {
    VecFloat* b = NULL;
    if (VecLoad(&b, stream) == false)
      return false;
    VecCopy(GABoundsAdnFloat(*that, iBound), b);
    VecFree(&b);
  }
  for (int iBound = 0; iBound < lenAdnI; ++iBound) {
    VecShort* b = NULL;
    if (VecLoad(&b, stream) == false)
      return false;
    VecCopy(GABoundsAdnInt(*that, iBound), b);
    VecFree(&b);
  }
  // Load the ELO rank
  for (int iEnt = 0; iEnt < nbEnt; ++iEnt) {
    GSetElem* setElem = GSetGetElem(&(GAEloRank(*that)->_set), iEnt);
    ELOEntity* eloEnt = (ELOEntity*)(setElem->_data);
    GenAlgEntity* ent = (GenAlgEntity*)(eloEnt->_data);
    // Load the id, age and elo
    int id, age;
    float elo;
    int ret = fscanf(stream, "%d %d %f", &id, &age, &elo);
    // If we couldn't fscanf
    if (ret == EOF)
      return false;
    // Set the id and elo
    ent->_id = id;
    ent->_age = age;
    setElem->_sortVal = elo;
    // Load the genes
    if (lenAdnF > 0) {
      if (VecLoad(&(ent->_adnF), stream) == false)
        return false;
      if (VecLoad(&(ent->_deltaAdnF), stream) == false)
        return false;
    }
    if (lenAdnI > 0)
      if (VecLoad(&(ent->_adnI), stream) == false)
        return false;
  }
  // Return success code
  return true;
}

// Save the GenAlg 'that' to the stream 'stream'
// Return true in case of success, else false
bool GASave(GenAlg* that, FILE* stream) {
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
  // Save the number of entity and elite, and the length of adn
  int ret = fprintf(stream, "%d %d %d %d\n", GAGetNbEntities(that),
    GAGetNbElites(that), GAGetLengthAdnFloat(that), 
    GAGetLengthAdnInt(that));
  // If we couldn't fprintf
  if (ret < 0)
    return false;
  // Save the run, epoch, nbRunPerEpoch, nextId
  ret = fprintf(stream, "%d %d %d %d\n", GAGetCurRun(that),
    GAGetCurEpoch(that), GAGetRunsPerEpoch(that), that->_nextId);
  // If we couldn't fprintf
  if (ret < 0)
    return false;
  // Save the bounds
  for (int iBound = 0; iBound < GAGetLengthAdnFloat(that); ++iBound)
    if (VecSave(GABoundsAdnFloat(that, iBound), stream) == false)
      return false;
  for (int iBound = 0; iBound < GAGetLengthAdnInt(that); ++iBound)
    if (VecSave(GABoundsAdnInt(that, iBound), stream) == false)
      return false;
  // Save the ELO rank
  for (int iEnt = 0; iEnt < GAGetNbEntities(that); ++iEnt) {
    GSetElem* setElem = GSetGetElem(&(GAEloRank(that)->_set), iEnt);
    ELOEntity* eloEnt = (ELOEntity*)(setElem->_data);
    GenAlgEntity* ent = (GenAlgEntity*)(eloEnt->_data);
    // Save the id, age and elo
    int ret = fprintf(stream, "%d %d %f\n", ent->_id, ent->_age, 
      setElem->_sortVal);
    // If we couldn't fprintf
    if (ret < 0)
      return false;
    // Save the genes
    if (GAGetLengthAdnFloat(that) > 0) {
      if (VecSave(ent->_adnF, stream) == false)
        return false;
      if (VecSave(ent->_deltaAdnF, stream) == false)
        return false;
    }
    if (GAGetLengthAdnInt(that) > 0)
      if (VecSave(ent->_adnI, stream) == false)
        return false;
  }
  // Return success code
  return true;
}

