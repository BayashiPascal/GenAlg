#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include "genalg.h"

#define RANDOMSEED 2

void UnitTestGenAlgEntityCreateFree() {
  int id = 1;
  int lengthAdnF = 2;
  int lengthAdnI = 3;
  GenAlgEntity* ent = GenAlgEntityCreate(id, lengthAdnF, lengthAdnI);
  if (ent->_age != 1 ||
    ent->_id != id ||
    VecGetDim(ent->_adnF) != lengthAdnF ||
    VecGetDim(ent->_deltaAdnF) != lengthAdnF ||
    VecGetDim(ent->_adnI) != lengthAdnI) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GenAlgEntityCreate failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgEntityFree(&ent);
  if (ent != NULL) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GenAlgEntityFree failed");
    PBErrCatch(GenAlgErr);
  }
  printf("UnitTestGenAlgEntityCreateFree OK\n");
}

void UnitTestGenAlgEntityGetSet() {
  int id = 1;
  int lengthAdnF = 2;
  int lengthAdnI = 3;
  GenAlgEntity* ent = GenAlgEntityCreate(id, lengthAdnF, lengthAdnI);
  if (GAEntAdnF(ent) != ent->_adnF) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntAdnF failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAEntDeltaAdnF(ent) != ent->_deltaAdnF) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntDeltaAdnF failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAEntAdnI(ent) != ent->_adnI) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntAdnI failed");
    PBErrCatch(GenAlgErr);
  }
  GAEntSetGeneF(ent, 0, 1.0);
  if (ISEQUALF(VecGet(ent->_adnF, 0), 1.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntSetGeneF failed");
    PBErrCatch(GenAlgErr);
  }
  if (ISEQUALF(GAEntGetGeneF(ent, 0), 1.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntGetGeneF failed");
    PBErrCatch(GenAlgErr);
  }
  GAEntSetDeltaGeneF(ent, 0, 2.0);
  if (ISEQUALF(VecGet(ent->_deltaAdnF, 0), 2.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntSetDeltaGeneF failed");
    PBErrCatch(GenAlgErr);
  }
  if (ISEQUALF(GAEntGetDeltaGeneF(ent, 0), 2.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntGetDeltaGeneF failed");
    PBErrCatch(GenAlgErr);
  }
  GAEntSetGeneI(ent, 0, 3);
  if (VecGet(ent->_adnI, 0) != 3) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntSetGeneI failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAEntGetGeneI(ent, 0) != 3) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntGetGeneI failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAEntGetAge(ent) != 1) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntGetAge failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAEntGetId(ent) != id) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntGetId failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgEntityFree(&ent);
  printf("UnitTestGenAlgEntityGetSet OK\n");
}

void UnitTestGenAlgEntityInit() {
  srandom(5);
  int id = 1;
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlgEntity* ent = GenAlgEntityCreate(id, lengthAdnF, lengthAdnI);
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES, 
    lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecShort2D boundsI = VecShortCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  VecCopy(GABoundsAdnFloat(ga, 0), &boundsF);
  VecCopy(GABoundsAdnFloat(ga, 1), &boundsF);
  VecCopy(GABoundsAdnInt(ga, 0), &boundsI);
  VecCopy(GABoundsAdnInt(ga, 1), &boundsI);
  GAEntInit(ent, ga);
  if (ISEQUALF(VecGet(ent->_adnF, 0), -0.907064) == false ||
    ISEQUALF(VecGet(ent->_adnF, 1), -0.450509) == false ||
    VecGet(ent->_adnI, 0) != 2 ||
    VecGet(ent->_adnI, 1) != 10) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEntInit failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  GenAlgEntityFree(&ent);
  printf("UnitTestGenAlgEntityInit OK\n");
}

void UnitTestGenAlgEntity() {
  UnitTestGenAlgEntityCreateFree();
  UnitTestGenAlgEntityGetSet();
  UnitTestGenAlgEntityInit();
  printf("UnitTestGenAlgEntity OK\n");
}

void UnitTestGenAlgCreateFree() {
  int lengthAdnF = 2;
  int lengthAdnI = 3;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES,
    lengthAdnF, lengthAdnI);
  if (ga->_runsPerEpoch != GENALG_RUNPEREPOCH ||
    ga->_curRun != 0 ||
    ga->_curEpoch != 0 ||
    ga->_nextId != GENALG_NBENTITIES ||
    ga->_nbEntities != GENALG_NBENTITIES ||
    ga->_nbElites != GENALG_NBELITES ||
    ga->_lengthAdnF != lengthAdnF ||
    ga->_lengthAdnI != lengthAdnI ||
    ELORankGetNb(GAEloRank(ga)) != GENALG_NBENTITIES) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GenAlgCreate failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  if (ga != NULL) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GenAlgFree failed");
    PBErrCatch(GenAlgErr);
  }
  printf("UnitTestGenAlgCreateFree OK\n");
}

void UnitTestGenAlgGetSet() {
  int lengthAdnF = 2;
  int lengthAdnI = 3;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES,
    lengthAdnF, lengthAdnI);
  if (GAEloRank(ga) != ga->_elo) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEloRank failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetRunsPerEpoch(ga) != GENALG_RUNPEREPOCH) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetRunsPerEpoch failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNbEntities(ga) != GENALG_NBENTITIES) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNbEntities failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNbElites(ga) != GENALG_NBELITES) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNbElites failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetCurRun(ga) != 0) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetCurRun failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetCurEpoch(ga) != 0) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetCurEpoch failed");
    PBErrCatch(GenAlgErr);
  }
  GASetRunsPerEpoch(ga, 10);
  if (GAGetRunsPerEpoch(ga) != 10) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetRunsPerEpoch failed");
    PBErrCatch(GenAlgErr);
  }
  GASetNbEntities(ga, 10);
  if (GAGetNbEntities(ga) != 10 ||
    GAGetNbElites(ga) != 9 ||
    ELORankGetNb(GAEloRank(ga)) != 10) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetNbEntities failed");
    PBErrCatch(GenAlgErr);
  }
  GASetNbElites(ga, 20);
  if (GAGetNbEntities(ga) != 21 ||
    GAGetNbElites(ga) != 20 ||
    ELORankGetNb(GAEloRank(ga)) != 21) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetNbElites failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetLengthAdnFloat(ga) != lengthAdnF) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetLengthAdnFloat failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetLengthAdnInt(ga) != lengthAdnI) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetLengthAdnInt failed");
    PBErrCatch(GenAlgErr);
  }
  if (GABoundsAdnFloat(ga, 1) != ga->_boundsF + 1) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GABoundsAdnFloat failed");
    PBErrCatch(GenAlgErr);
  }
  if (GABoundsAdnInt(ga, 1) != ga->_boundsI + 1) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GABoundsAdnInt failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  printf("UnitTestGenAlgGetSet OK\n");
}

void UnitTestGenAlgInit() {
  srandom(5);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES,
    lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecShort2D boundsI = VecShortCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  VecCopy(GABoundsAdnFloat(ga, 0), &boundsF);
  VecCopy(GABoundsAdnFloat(ga, 1), &boundsF);
  VecCopy(GABoundsAdnInt(ga, 0), &boundsI);
  VecCopy(GABoundsAdnInt(ga, 1), &boundsI);
  GAInit(ga);
  GenAlgEntity* ent = 
    (GenAlgEntity*)(
    ((ELOEntity*)(GAEloRank(ga)->_set._head->_data))->_data);
  if (ISEQUALF(VecGet(ent->_adnF, 0), -0.907064) == false ||
    ISEQUALF(VecGet(ent->_adnF, 1), -0.450509) == false ||
    VecGet(ent->_adnI, 0) != 2 ||
    VecGet(ent->_adnI, 1) != 10) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAInit failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  printf("UnitTestGenAlgInit OK\n");
}

void UnitTestGenAlgPrint() {
  srandom(5);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(3, 2, lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecShort2D boundsI = VecShortCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  VecCopy(GABoundsAdnFloat(ga, 0), &boundsF);
  VecCopy(GABoundsAdnFloat(ga, 1), &boundsF);
  VecCopy(GABoundsAdnInt(ga, 0), &boundsI);
  VecCopy(GABoundsAdnInt(ga, 1), &boundsI);
  GAInit(ga);
  GAPrintln(ga, stdout);
  GenAlgFree(&ga);
  printf("UnitTestGenAlgInit OK\n");
}

void UnitTestGenAlgStepRun() {
  srandom(5);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES,
    lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecShort2D boundsI = VecShortCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  VecCopy(GABoundsAdnFloat(ga, 0), &boundsF);
  VecCopy(GABoundsAdnFloat(ga, 1), &boundsF);
  VecCopy(GABoundsAdnInt(ga, 0), &boundsI);
  VecCopy(GABoundsAdnInt(ga, 1), &boundsI);
  GASetNbElites(ga, 2);
  GASetNbEntities(ga, 3);
  GAInit(ga);
  GSet* rank = GSetCreate();
  for (int i = 3; i--;)
    GSetAddSort(rank, GAEntity(ga, i), 3.0 - (float)i);
  GAStepRun(ga, rank);
  if (GAGetCurRun(ga) != 1) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAStepRun failed");
    PBErrCatch(GenAlgErr);
  }
  for (int i = 3; i--;) {
    if (ELORankGetRank(ga->_elo, GSetGet(rank, i)) != 2 - i ||
      ISEQUALF(ELORankGetELO(ga->_elo, GSetGet(rank, i)), 
        92.0 + (float)i * 8.0) == false ||
      GAEntity(ga, i)->_age != 2) {
      GenAlgErr->_type = PBErrTypeUnitTestFailed;
      sprintf(GenAlgErr->_msg, "GAStepRun failed");
      PBErrCatch(GenAlgErr);
    }
  }
  GSetFree(&rank);
  GenAlgFree(&ga);
  printf("UnitTestGenAlgStepRun OK\n");
}

void UnitTestGenAlgGetInbreeding() {
  srandom(5);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES,
    lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecShort2D boundsI = VecShortCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  VecCopy(GABoundsAdnFloat(ga, 0), &boundsF);
  VecCopy(GABoundsAdnFloat(ga, 1), &boundsF);
  VecCopy(GABoundsAdnInt(ga, 0), &boundsI);
  VecCopy(GABoundsAdnInt(ga, 1), &boundsI);
  GASetNbElites(ga, 2);
  GASetNbEntities(ga, 3);
  GAInit(ga);
  if (ISEQUALF(GAGetInbreeding(ga), 0.182041) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetInbreeding failed");
    PBErrCatch(GenAlgErr);
  }
  VecCopy(GAEntity(ga, 1)->_adnF, GAEntity(ga, 0)->_adnF);
  VecCopy(GAEntity(ga, 1)->_adnI, GAEntity(ga, 0)->_adnI);
  if (ISEQUALF(GAGetInbreeding(ga), 0.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetInbreeding failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  printf("UnitTestGenAlgGetInbreeding OK\n");
}

void UnitTestGenAlgStepEpoch() {
  srandom(5);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(3, 2, lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecShort2D boundsI = VecShortCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  VecCopy(GABoundsAdnFloat(ga, 0), &boundsF);
  VecCopy(GABoundsAdnFloat(ga, 1), &boundsF);
  VecCopy(GABoundsAdnInt(ga, 0), &boundsI);
  VecCopy(GABoundsAdnInt(ga, 1), &boundsI);
  GAInit(ga);
  GSet* rank = GSetCreate();
  for (int i = 3; i--;)
    GSetAddSort(rank, GAEntity(ga, i), 3.0 - (float)i);
  GAStepRun(ga, rank);
  printf("Before StepEpoch:\n");
  GAPrintln(ga, stdout);
  GenAlgEntity* child = GAEntity(ga, 2);
  GAStepEpoch(ga);
  printf("After StepEpoch:\n");
  GAPrintln(ga, stdout);
  if (ga->_nextId != 4 || GAEntGetId(child) != 3 || 
    GAEntGetAge(child) != 1 ||
    ISEQUALF(GAEntGetGeneF(child, 0), 0.755265) == false ||
    ISEQUALF(GAEntGetGeneF(child, 1), -0.209552) == false ||
    ISEQUALF(GAEntGetDeltaGeneF(child, 0), -0.032739) == false ||
    ISEQUALF(GAEntGetDeltaGeneF(child, 1), -0.206048) == false ||
    GAEntGetGeneI(child, 0) != 4 ||
    GAEntGetGeneI(child, 1) != 1 ||
    GAEntity(ga, 1) != child ||
    GAEntGetAge(GAEntity(ga, 0)) != 2 ||
    GAEntGetAge(GAEntity(ga, 2)) != 2 ||
    GAEntGetId(GAEntity(ga, 0)) != 2 ||
    GAEntGetId(GAEntity(ga, 2)) != 1 ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), child), 100.0) == false ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), GAEntity(ga, 0)), 
      108.0) == false ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), GAEntity(ga, 2)), 
      100.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAStepEpoch failed");
    PBErrCatch(GenAlgErr);
  }
  VecCopy(GAEntity(ga, 1)->_adnF, GAEntity(ga, 0)->_adnF);
  VecCopy(GAEntity(ga, 1)->_adnI, GAEntity(ga, 0)->_adnI);
  GAStepEpoch(ga);
  printf("After StepEpoch with interbreeding:\n");
  GAPrintln(ga, stdout);
  if (ga->_nextId != 6 || GAEntGetId(child) != 4 || 
    GAEntGetAge(child) != 1 ||
    ISEQUALF(GAEntGetGeneF(child, 0), 0.788004) == false ||
    ISEQUALF(GAEntGetGeneF(child, 1), -0.003504) == false ||
    ISEQUALF(GAEntGetDeltaGeneF(child, 0), -0.032739) == false ||
    ISEQUALF(GAEntGetDeltaGeneF(child, 1), -0.206048) == false ||
    GAEntGetGeneI(child, 0) != 3 ||
    GAEntGetGeneI(child, 1) != 1 ||
    GAEntity(ga, 2) != child ||
    GAEntGetAge(GAEntity(ga, 0)) != 2 ||
    GAEntGetAge(GAEntity(ga, 1)) != 1 ||
    GAEntGetId(GAEntity(ga, 0)) != 2 ||
    GAEntGetId(GAEntity(ga, 1)) != 5 ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), child), 100.0) == false ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), GAEntity(ga, 0)), 
      108.0) == false ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), GAEntity(ga, 1)), 
      100.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAStepEpoch failed");
    PBErrCatch(GenAlgErr);
  }
  GSetFree(&rank);
  GenAlgFree(&ga);
  printf("UnitTestGenAlgStepEpoch OK\n");
}

void UnitTestGenAlgLoadSave() {
  srandom(5);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(3, 2, lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecShort2D boundsI = VecShortCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  VecCopy(GABoundsAdnFloat(ga, 0), &boundsF);
  VecCopy(GABoundsAdnFloat(ga, 1), &boundsF);
  VecCopy(GABoundsAdnInt(ga, 0), &boundsI);
  VecCopy(GABoundsAdnInt(ga, 1), &boundsI);
  GAInit(ga);
  GAStepEpoch(ga);
  GSet* rank = GSetCreate();
  for (int i = 3; i--;)
    GSetAddSort(rank, GAEntity(ga, i), 3.0 - (float)i);
  GAStepRun(ga, rank);
  FILE* stream = fopen("./UnitTestGenAlgLoadSave.txt", "w");
  if (GASave(ga, stream) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASave failed");
    PBErrCatch(GenAlgErr);
  }
  fclose(stream);
  stream = fopen("./UnitTestGenAlgLoadSave.txt", "r");
  GenAlg* gaLoad = NULL;
  if (GALoad(&gaLoad, stream) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GALoad failed");
    PBErrCatch(GenAlgErr);
  }
  fclose(stream);
  if (ga->_nextId != gaLoad->_nextId||
    ga->_curEpoch != gaLoad->_curEpoch ||
    ga->_curRun != gaLoad->_curRun ||
    ga->_runsPerEpoch != gaLoad->_runsPerEpoch ||
    ga->_nbEntities != gaLoad->_nbEntities ||
    ga->_nbElites != gaLoad->_nbElites ||
    ga->_lengthAdnF != gaLoad->_lengthAdnF ||
    ga->_lengthAdnI != gaLoad->_lengthAdnI ||
    VecIsEqual(ga->_boundsF, gaLoad->_boundsF) == false ||
    VecIsEqual(ga->_boundsF + 1, gaLoad->_boundsF + 1) == false ||
    VecIsEqual(ga->_boundsI, gaLoad->_boundsI) == false ||
    VecIsEqual(ga->_boundsI + 1, gaLoad->_boundsI + 1) == false ||
    GAEntGetId(GAEntity(ga, 0)) != GAEntGetId(GAEntity(gaLoad, 0)) ||
    GAEntGetId(GAEntity(ga, 1)) != GAEntGetId(GAEntity(gaLoad, 1)) ||
    GAEntGetId(GAEntity(ga, 2)) != GAEntGetId(GAEntity(gaLoad, 2)) ||
    GAEntGetAge(GAEntity(ga, 0)) != GAEntGetAge(GAEntity(gaLoad, 0)) ||
    GAEntGetAge(GAEntity(ga, 1)) != GAEntGetAge(GAEntity(gaLoad, 1)) ||
    GAEntGetAge(GAEntity(ga, 2)) != GAEntGetAge(GAEntity(gaLoad, 2)) ||
    VecIsEqual(GAEntity(ga, 0)->_adnF, 
      GAEntity(gaLoad, 0)->_adnF) == false ||
    VecIsEqual(GAEntity(ga, 0)->_deltaAdnF, 
      GAEntity(gaLoad, 0)->_deltaAdnF) == false ||
    VecIsEqual(GAEntity(ga, 0)->_adnI, 
      GAEntity(gaLoad, 0)->_adnI) == false ||
    VecIsEqual(GAEntity(ga, 1)->_adnF, 
      GAEntity(gaLoad, 1)->_adnF) == false ||
    VecIsEqual(GAEntity(ga, 1)->_deltaAdnF, 
      GAEntity(gaLoad, 1)->_deltaAdnF) == false ||
    VecIsEqual(GAEntity(ga, 1)->_adnI, 
      GAEntity(gaLoad, 1)->_adnI) == false ||
    VecIsEqual(GAEntity(ga, 2)->_adnF, 
      GAEntity(gaLoad, 2)->_adnF) == false ||
    VecIsEqual(GAEntity(ga, 2)->_deltaAdnF, 
      GAEntity(gaLoad, 2)->_deltaAdnF) == false ||
    VecIsEqual(GAEntity(ga, 2)->_adnI, 
      GAEntity(gaLoad, 2)->_adnI) == false ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), GAEntity(ga, 0)), 
      ELORankGetELO(GAEloRank(gaLoad), GAEntity(gaLoad, 0))) == false ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), GAEntity(ga, 1)), 
      ELORankGetELO(GAEloRank(gaLoad), GAEntity(gaLoad, 1))) == false ||
    ISEQUALF(ELORankGetELO(GAEloRank(ga), GAEntity(ga, 2)), 
      ELORankGetELO(GAEloRank(gaLoad), GAEntity(gaLoad, 2))) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "UnitTestGenAlgLoadSave failed");
    PBErrCatch(GenAlgErr);
  }
  GSetFree(&rank);
  GenAlgFree(&ga);
  GenAlgFree(&gaLoad);
  printf("UnitTestGenAlgLoadSave OK\n");
}

float ftarget(float x) {
  return -0.5 * fastpow(x, 3) + 0.314 * fastpow(x, 2) - 0.7777 * x + 0.1;
}

float evaluate(VecFloat* adnF, VecShort* adnI) {
  float delta = 0.02;
  int nb = (int)round(4.0 / delta);
  float res = 0.0;
  float x = -2.0;
  for (int i = 0; i < nb; ++i, x += delta) {
    float y = 0.0;
    for (int j = 4; j--;)
      y += VecGet(adnF, j) * fastpow(x, VecGet(adnI, j));
    res += fabs(ftarget(x) - y);
  }
  return res / (float)nb;
}

void UnitTestGenAlgTest() {
  srandom(5);
  int lengthAdnF = 4;
  int lengthAdnI = lengthAdnF;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES, 
    lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecShort2D boundsI = VecShortCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 0); VecSet(&boundsI, 1, 4);
  for (int i = lengthAdnF; i--;) {
    VecCopy(GABoundsAdnFloat(ga, i), &boundsF);
    VecCopy(GABoundsAdnInt(ga, i), &boundsI);
  }
  GAInit(ga);
  GASetRunsPerEpoch(ga, 1);
  GSet* rank = GSetCreate();
//float best = 1.0;
//int step = 0;
  do {
    for (int iEnt = GAGetNbEntities(ga); iEnt--;)
      GSetAddSort(rank, GAEntity(ga, iEnt), 
        -1.0 * evaluate(GAEntAdnF(GAEntity(ga, iEnt)), 
        GAEntAdnI(GAEntity(ga, iEnt))));
    GAStepRun(ga, rank);
    GSetFlush(rank);
/*float ev = evaluate(GAEntAdnF(GAEntity(ga, 0)), 
GAEntAdnI(GAEntity(ga, 0)));
if (step == 10){
  printf("%d %f %f\n",GAGetCurEpoch(ga), ev, GAGetInbreeding(ga));
  step = 0;
} else step++;*/
/*if (best - ev > PBMATH_EPSILON) {
  best = ev;
  printf("%d %f ", GAGetCurEpoch(ga), best);
  VecFloatPrint(GAEntAdnF(GAEntity(ga, 0)), stdout, 6);
  printf(" ");
  VecPrint(GAEntAdnI(GAEntity(ga, 0)), stdout);
  printf("\n");
}*/
  } while (GAGetCurEpoch(ga) < 20000 || 
    evaluate(GAEntAdnF(GAEntity(ga, 0)), 
      GAEntAdnI(GAEntity(ga, 0))) < PBMATH_EPSILON);
  printf("target: -0.5*x^3 + 0.314*x^2 - 0.7777*x + 0.1\n");
  printf("approx: \n");
  GAEntPrintln(GAEntity(ga, 0), stdout);
  printf("error: %f\n", evaluate(GAEntAdnF(GAEntity(ga, 0)), 
    GAEntAdnI(GAEntity(ga, 0))));
  GSetFree(&rank);
  GenAlgFree(&ga);
  printf("UnitTestGenAlgTest OK\n");
}

void UnitTestGenAlg() {
  UnitTestGenAlgCreateFree();
  UnitTestGenAlgGetSet();
  UnitTestGenAlgInit();
  UnitTestGenAlgPrint();
  UnitTestGenAlgStepRun();
  UnitTestGenAlgGetInbreeding();
  UnitTestGenAlgStepEpoch();
  UnitTestGenAlgLoadSave();
  UnitTestGenAlgTest();
  printf("UnitTestGenAlg OK\n");
}

void UnitTestAll() {
  UnitTestGenAlgEntity();
  UnitTestGenAlg();
  printf("UnitTestAll OK\n");
}

int main() {
  UnitTestAll();
  // Return success code
  return 0;
}

