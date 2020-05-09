#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include "genalg.h"

#define RANDOMSEED 2

void UnitTestGenAlgAdnCreateFree() {
  unsigned long int id = 1;
  int lengthAdnF = 2;
  int lengthAdnI = 3;
  GenAlgAdn* ent = GenAlgAdnCreate(id, lengthAdnF, lengthAdnI);
  if (ent->_age != 1 ||
    ent->_id != id ||
    VecGetDim(ent->_adnF) != lengthAdnF ||
    VecGetDim(ent->_deltaAdnF) != lengthAdnF ||
    VecGetDim(ent->_adnI) != lengthAdnI) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GenAlgAdnCreate failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgAdnFree(&ent);
  if (ent != NULL) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GenAlgAdnFree failed");
    PBErrCatch(GenAlgErr);
  }
  printf("UnitTestGenAlgAdnCreateFree OK\n");
}

void UnitTestGenAlgAdnGetSet() {
  unsigned long int id = 1;
  int lengthAdnF = 2;
  int lengthAdnI = 3;
  GenAlgAdn* ent = GenAlgAdnCreate(id, lengthAdnF, lengthAdnI);
  if (GAAdnAdnF(ent) != ent->_adnF) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnAdnF failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAAdnDeltaAdnF(ent) != ent->_deltaAdnF) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnDeltaAdnF failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAAdnAdnI(ent) != ent->_adnI) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnAdnI failed");
    PBErrCatch(GenAlgErr);
  }
  GAAdnSetGeneF(ent, 0, 1.0);
  if (ISEQUALF(VecGet(ent->_adnF, 0), 1.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnSetGeneF failed");
    PBErrCatch(GenAlgErr);
  }
  if (ISEQUALF(GAAdnGetGeneF(ent, 0), 1.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnGetGeneF failed");
    PBErrCatch(GenAlgErr);
  }
  GAAdnSetDeltaGeneF(ent, 0, 2.0);
  if (ISEQUALF(VecGet(ent->_deltaAdnF, 0), 2.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnSetDeltaGeneF failed");
    PBErrCatch(GenAlgErr);
  }
  if (ISEQUALF(GAAdnGetDeltaGeneF(ent, 0), 2.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnGetDeltaGeneF failed");
    PBErrCatch(GenAlgErr);
  }
  GAAdnSetGeneI(ent, 0, 3);
  if (VecGet(ent->_adnI, 0) != 3) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnSetGeneI failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAAdnGetGeneI(ent, 0) != 3) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnGetGeneI failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAAdnGetAge(ent) != 1) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnGetAge failed");
    PBErrCatch(GenAlgErr);
  }
  ent->_val = 2.0;
  if (ISEQUALF(GAAdnGetVal(ent), 2.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnGetVal failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAAdnGetId(ent) != id) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnGetId failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAAdnIsNew(ent) != true) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnIsNew failed");
    PBErrCatch(GenAlgErr);
  }
  ent->_age = 2;
  if (GAAdnIsNew(ent) != false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnIsNew failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgAdnFree(&ent);
  printf("UnitTestGenAlgAdnGetSet OK\n");
}

void UnitTestGenAlgAdnInit() {
  srandom(5);
  unsigned long int id = 1;
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlgAdn* ent = GenAlgAdnCreate(id, lengthAdnF, lengthAdnI);
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES, 
    lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecLong2D boundsI = VecLongCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  GASetBoundsAdnFloat(ga, 0, &boundsF);
  GASetBoundsAdnFloat(ga, 1, &boundsF);
  GASetBoundsAdnInt(ga, 0, &boundsI);
  GASetBoundsAdnInt(ga, 1, &boundsI);
  GAAdnInit(ent, ga);
  if (ISEQUALF(VecGet(ent->_adnF, 0), -0.907064) == false ||
    ISEQUALF(VecGet(ent->_adnF, 1), -0.450509) == false ||
    VecGet(ent->_adnI, 0) != 2 ||
    VecGet(ent->_adnI, 1) != 10) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAAdnInit failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  GenAlgAdnFree(&ent);
  printf("UnitTestGenAlgAdnInit OK\n");
}

void UnitTestGenAlgAdn() {
  UnitTestGenAlgAdnCreateFree();
  UnitTestGenAlgAdnGetSet();
  UnitTestGenAlgAdnInit();
  printf("UnitTestGenAlgAdn OK\n");
}

void UnitTestGenAlgCreateFree() {
  int lengthAdnF = 2;
  int lengthAdnI = 3;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES,
    lengthAdnF, lengthAdnI);
  if (ga->_type != genAlgTypeDefault ||
    ga->_curEpoch != 0 ||
    ga->_nbKTEvent != 0 ||
    ga->_nextId != GENALG_NBENTITIES ||
    ga->_nbElites != GENALG_NBELITES ||
    ga->_lengthAdnF != lengthAdnF ||
    ga->_lengthAdnI != lengthAdnI ||
    ga->_flagTextOMeter != false ||
    ga->_nbMinAdn != GENALG_NBENTITIES ||
    ga->_nbMaxAdn != GENALG_NBENTITIES ||
    ISEQUALF(ga->_diversityThreshold, PBMATH_EPSILON) != true ||
    ga->_textOMeter != NULL ||
    GSetNbElem(GAAdns(ga)) != GENALG_NBENTITIES) {
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
  if (GAGetType(ga) != ga->_type) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetType failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAAdns(ga) != ga->_adns) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAEloRank failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNbAdns(ga) != GENALG_NBENTITIES) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNbAdns failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNbElites(ga) != GENALG_NBELITES) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNbElites failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetCurEpoch(ga) != 0) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetCurEpoch failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNbKTEvent(ga) != 0) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNbKTEvent failed");
    PBErrCatch(GenAlgErr);
  }
  if (ISEQUALF(GAGetDiversityThreshold(ga), 
    ga->_diversityThreshold) != true) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetDiversityThreshold failed");
    PBErrCatch(GenAlgErr);
  }
  GASetDiversityThreshold(ga, 2.0);
  if (ISEQUALF(GAGetDiversityThreshold(ga), 2.0) != true) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetDiversityThrehsold failed");
    PBErrCatch(GenAlgErr);
  }
  GASetNbEntities(ga, 10);
  if (GAGetNbAdns(ga) != 10 ||
    GAGetNbElites(ga) != 9 ||
    GSetNbElem(GAAdns(ga)) != 10) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetNbEntities failed");
    PBErrCatch(GenAlgErr);
  }
  GASetNbElites(ga, 20);
  if (GAGetNbAdns(ga) != 21 ||
    GAGetNbElites(ga) != 20 ||
    GSetNbElem(GAAdns(ga)) != 21) {
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
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  GASetBoundsAdnFloat(ga, 1, &boundsF);
  if (VecIsEqual(GABoundsAdnFloat(ga, 1), &boundsF) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetBoundsAdnFloat failed");
    PBErrCatch(GenAlgErr);
  }
  VecLong2D boundsS = VecLongCreateStatic2D();
  VecSet(&boundsS, 0, -1); VecSet(&boundsS, 1, 1);
  GASetBoundsAdnInt(ga, 1, &boundsS);
  if (VecIsEqual(GABoundsAdnInt(ga, 1), &boundsS) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetBoundsAdnInt failed");
    PBErrCatch(GenAlgErr);
  }
  if (GABoundsAdnInt(ga, 1) != ga->_boundsI + 1) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GABoundsAdnInt failed");
    PBErrCatch(GenAlgErr);
  }
  GASetAdnValue(ga, GAAdn(ga, 0), 1.0);
  if (ISEQUALF(GAAdn(ga, 0)->_val, 1.0) == false ||
    ISEQUALF(ga->_adns->_tail->_sortVal, 1.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetAdnValue failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNbMaxAdn(ga) != ga->_nbMaxAdn) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNbMaxAdn failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNbMinAdn(ga) != ga->_nbMinAdn) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNbMinAdn failed");
    PBErrCatch(GenAlgErr);
  }
  GASetNbMaxAdn(ga, 100);
  if (GAGetNbMaxAdn(ga) != 100) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetNbMaxAdn failed");
    PBErrCatch(GenAlgErr);
  }
  GASetNbMinAdn(ga, 100);
  if (GAGetNbMinAdn(ga) != 100) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetNbMinAdn failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES, 3, 3);
  GASetTypeNeuraNet(ga, 1, 2, 3);
  if (GAGetType(ga) != genAlgTypeNeuraNet ||
    ga->_NNdata._nbIn != 1 ||
    ga->_NNdata._nbHid != 2 ||
    ga->_NNdata._nbOut != 3) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetTypeNeuraNet failed");
    PBErrCatch(GenAlgErr);
  }
  GASetNeuraNetLinkMutability(ga, true);
  if (ga->_NNdata._flagMutableLink != true) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetNeuraNetLinkMutability failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNeuraNetLinkMutability(ga) != true) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNeuraNetLinkMutability failed");
    PBErrCatch(GenAlgErr);
  }
  GASetNeuraNetLinkMutability(ga, false);
  if (ga->_NNdata._flagMutableLink != false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GASetNeuraNetLinkMutability failed");
    PBErrCatch(GenAlgErr);
  }
  if (GAGetNeuraNetLinkMutability(ga) != false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetNeuraNetLinkMutability failed");
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
  VecLong2D boundsI = VecLongCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  GASetBoundsAdnFloat(ga, 0, &boundsF);
  GASetBoundsAdnFloat(ga, 1, &boundsF);
  GASetBoundsAdnInt(ga, 0, &boundsI);
  GASetBoundsAdnInt(ga, 1, &boundsI);
  GAInit(ga);
  GenAlgAdn* ent = (GenAlgAdn*)(GAAdns(ga)->_head->_data);
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
  VecLong2D boundsI = VecLongCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  GASetBoundsAdnFloat(ga, 0, &boundsF);
  GASetBoundsAdnFloat(ga, 1, &boundsF);
  GASetBoundsAdnInt(ga, 0, &boundsI);
  GASetBoundsAdnInt(ga, 1, &boundsI);
  GAInit(ga);
  GAPrintln(ga, stdout);
  GAEliteSummaryPrintln(ga, stdout);
  GenAlgFree(&ga);
  printf("UnitTestGenAlgInit OK\n");
}

void UnitTestGenAlgGetDiversity() {
  srandom(5);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES,
    lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecLong2D boundsI = VecLongCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  GASetBoundsAdnFloat(ga, 0, &boundsF);
  GASetBoundsAdnFloat(ga, 1, &boundsF);
  GASetBoundsAdnInt(ga, 0, &boundsI);
  GASetBoundsAdnInt(ga, 1, &boundsI);
  GASetNbElites(ga, 2);
  GASetNbEntities(ga, 3);
  GAInit(ga);
  if (ISEQUALF(GAGetDiversity(ga), 0.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetDiversity failed");
    PBErrCatch(GenAlgErr);
  }
  VecCopy(GAAdn(ga, 1)->_adnF, GAAdn(ga, 0)->_adnF);
  VecCopy(GAAdn(ga, 1)->_adnI, GAAdn(ga, 0)->_adnI);
  if (ISEQUALF(GAGetDiversity(ga), 0.0) == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAGetDiversity failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  printf("UnitTestGenAlgGetDiversity OK\n");
}

void UnitTestGenAlgStep() {
  srandom(2);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(3, 2, lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecLong2D boundsI = VecLongCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  GASetBoundsAdnFloat(ga, 0, &boundsF);
  GASetBoundsAdnFloat(ga, 1, &boundsF);
  GASetBoundsAdnInt(ga, 0, &boundsI);
  GASetBoundsAdnInt(ga, 1, &boundsI);
  GAInit(ga);
  for (int i = 3; i--;)
    GASetAdnValue(ga, GAAdn(ga, i), 3.0 - (float)i);
  printf("Before Step:\n");
  GAPrintln(ga, stdout);
  GenAlgAdn* child = GAAdn(ga, 2);
  GAStep(ga);
  printf("After Step:\n");
  GAPrintln(ga, stdout);
  if (ga->_nextId != 4 || GAAdnGetId(child) != 3 || 
    GAAdnGetAge(child) != 1 ||
    ISEQUALF(GAAdnGetGeneF(child, 0), 0.285933) == false ||
    ISEQUALF(GAAdnGetGeneF(child, 1), 0.287805) == false ||
    ISEQUALF(GAAdnGetDeltaGeneF(child, 0), 0.0) == false ||
    ISEQUALF(GAAdnGetDeltaGeneF(child, 1), 0.112841) == false ||
    GAAdnGetGeneI(child, 0) != 4 ||
    GAAdnGetGeneI(child, 1) != 10 ||
    GAAdn(ga, 2) != child ||
    GAAdnGetAge(GAAdn(ga, 0)) != 2 ||
    GAAdnGetAge(GAAdn(ga, 1)) != 2 ||
    GAAdnGetId(GAAdn(ga, 0)) != 0 ||
    GAAdnGetId(GAAdn(ga, 1)) != 1) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAStep failed");
    PBErrCatch(GenAlgErr);
  }
  GenAlgFree(&ga);
  printf("UnitTestGenAlgStep OK\n");
}

void UnitTestGenAlgLoadSave() {
  srandom(5);
  int lengthAdnF = 2;
  int lengthAdnI = 2;
  GenAlg* ga = GenAlgCreate(3, 2, lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecLong2D boundsI = VecLongCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 1); VecSet(&boundsI, 1, 10);
  GASetBoundsAdnFloat(ga, 0, &boundsF);
  GASetBoundsAdnFloat(ga, 1, &boundsF);
  GASetBoundsAdnInt(ga, 0, &boundsI);
  GASetBoundsAdnInt(ga, 1, &boundsI);
  GAInit(ga);
  GAStep(ga);
  GSet* rank = GSetCreate();
  for (int i = 3; i--;)
    GSetAddSort(rank, GAAdn(ga, i), 3.0 - (float)i);
  FILE* stream = fopen("./UnitTestGenAlgLoadSave.txt", "w");
  if (GASave(ga, stream, false) == false) {
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
  if (ga->_nextId != gaLoad->_nextId ||
    ga->_curEpoch != gaLoad->_curEpoch ||
    ga->_nbElites != gaLoad->_nbElites ||
    ga->_type != genAlgTypeDefault ||
    ga->_lengthAdnF != gaLoad->_lengthAdnF ||
    ga->_lengthAdnI != gaLoad->_lengthAdnI ||
    VecIsEqual(ga->_boundsF, gaLoad->_boundsF) == false ||
    VecIsEqual(ga->_boundsF + 1, gaLoad->_boundsF + 1) == false ||
    VecIsEqual(ga->_boundsI, gaLoad->_boundsI) == false ||
    VecIsEqual(ga->_boundsI + 1, gaLoad->_boundsI + 1) == false ||
    GAAdnGetId(GAAdn(ga, 0)) != GAAdnGetId(GAAdn(gaLoad, 0)) ||
    GAAdnGetId(GAAdn(ga, 1)) != GAAdnGetId(GAAdn(gaLoad, 1)) ||
    GAAdnGetId(GAAdn(ga, 2)) != GAAdnGetId(GAAdn(gaLoad, 2)) ||
    GAAdnGetAge(GAAdn(ga, 0)) != GAAdnGetAge(GAAdn(gaLoad, 0)) ||
    GAAdnGetAge(GAAdn(ga, 1)) != GAAdnGetAge(GAAdn(gaLoad, 1)) ||
    GAAdnGetAge(GAAdn(ga, 2)) != GAAdnGetAge(GAAdn(gaLoad, 2)) ||
    VecIsEqual(GAAdn(ga, 0)->_adnF, 
      GAAdn(gaLoad, 0)->_adnF) == false ||
    VecIsEqual(GAAdn(ga, 0)->_deltaAdnF, 
      GAAdn(gaLoad, 0)->_deltaAdnF) == false ||
    VecIsEqual(GAAdn(ga, 0)->_adnI, 
      GAAdn(gaLoad, 0)->_adnI) == false ||
    VecIsEqual(GAAdn(ga, 1)->_adnF, 
      GAAdn(gaLoad, 1)->_adnF) == false ||
    VecIsEqual(GAAdn(ga, 1)->_deltaAdnF, 
      GAAdn(gaLoad, 1)->_deltaAdnF) == false ||
    VecIsEqual(GAAdn(ga, 1)->_adnI, 
      GAAdn(gaLoad, 1)->_adnI) == false ||
    VecIsEqual(GAAdn(ga, 2)->_adnF, 
      GAAdn(gaLoad, 2)->_adnF) == false ||
    VecIsEqual(GAAdn(ga, 2)->_deltaAdnF, 
      GAAdn(gaLoad, 2)->_deltaAdnF) == false ||
    VecIsEqual(GAAdn(ga, 2)->_adnI, 
      GAAdn(gaLoad, 2)->_adnI) == false) {
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

float evaluate(const VecFloat* adnF, const VecLong* adnI) {
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
  srandom(0);
  int lengthAdnF = 4;
  int lengthAdnI = lengthAdnF;
  GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES, 
    lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecLong2D boundsI = VecLongCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 0); VecSet(&boundsI, 1, 4);
  for (int i = lengthAdnF; i--;) {
    GASetBoundsAdnFloat(ga, i, &boundsF);
    GASetBoundsAdnInt(ga, i, &boundsI);
  }
  GASetFlagHistory(ga, true);
  GASetHistoryPath(ga, "./history.json");
  GAInit(ga);
  GASetTextOMeterFlag(ga, true);
  GASetNbMinAdn(ga, GENALG_NBELITES * 2);
  GASetNbMaxAdn(ga, GENALG_NBENTITIES);
  float best = 1.0;
  unsigned long nbMaxEpoch = 2000;
  do {
    for (int iEnt = GAGetNbAdns(ga); iEnt--;)
      if (GAAdnIsNew(GAAdn(ga, iEnt))) {
        GASetAdnValue(ga, GAAdn(ga, iEnt), 
          -1.0 * evaluate(GAAdnAdnF(GAAdn(ga, iEnt)), 
          GAAdnAdnI(GAAdn(ga, iEnt))));
      }
    GAStep(ga);
    // Slow down the process to have time to read the TextOMeter
    unsigned int microseconds = 10000;
    usleep(microseconds);
    //sleep(1);
    // Display info if there is improvment
    float ev = evaluate(GABestAdnF(ga), GABestAdnI(ga));
    if (best - ev > PBMATH_EPSILON) {
      best = ev;
      printf("%lu %f ", GAGetCurEpoch(ga), best);
      VecFloatPrint(GABestAdnF(ga), stdout, 6);
      printf(" ");
      VecPrint(GABestAdnI(ga), stdout);
      printf("\n");
    }
  } while (GAGetCurEpoch(ga) < nbMaxEpoch && best > PBMATH_EPSILON);
  // Save the history
  bool ret = GASaveHistory(ga);
  if (ret == false) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "Couldn't save the history");
    PBErrCatch(GenAlgErr);
  }
  printf("target: -0.5*x^3 + 0.314*x^2 - 0.7777*x + 0.1\n");
  printf("approx: \n");
  GAAdnPrintln(GABestAdn(ga), stdout);
  printf("error: %f\n", evaluate(GABestAdnF(ga), GABestAdnI(ga)));
  GenAlgFree(&ga);
  printf("UnitTestGenAlgTest OK\n");
}

void UnitTestGenAlgPerf() {
  int nbRun = 10;
  unsigned long int nbMaxEpoch = 2000;
  float maxEv = 0.0;
  float bestEv = 0.0;
  float sumEv = 0.0;
  float avgEv = 0.0;
  for (int iRun = 0; iRun < nbRun; ++iRun) {
    srandom(time(NULL));
    int lengthAdnF = 4;
    int lengthAdnI = lengthAdnF;
    GenAlg* ga = GenAlgCreate(GENALG_NBENTITIES, GENALG_NBELITES, 
      lengthAdnF, lengthAdnI);
    VecFloat2D boundsF = VecFloatCreateStatic2D();
    VecLong2D boundsI = VecLongCreateStatic2D();
    VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
    VecSet(&boundsI, 0, 0); VecSet(&boundsI, 1, 4);
    for (int i = lengthAdnF; i--;) {
      GASetBoundsAdnFloat(ga, i, &boundsF);
      GASetBoundsAdnInt(ga, i, &boundsI);
    }
    GAInit(ga);
    GASetNbMinAdn(ga, GENALG_NBELITES * 2);
    GASetNbMaxAdn(ga, GENALG_NBENTITIES);
    float ev = 0.0;
    do {
      for (int iEnt = GAGetNbAdns(ga); iEnt--;)
        if (GAAdnIsNew(GAAdn(ga, iEnt)))
          GASetAdnValue(ga, GAAdn(ga, iEnt), 
            -1.0 * evaluate(GAAdnAdnF(GAAdn(ga, iEnt)), 
            GAAdnAdnI(GAAdn(ga, iEnt))));
      GAStep(ga);
      ev = evaluate(GABestAdnF(ga), GABestAdnI(ga));
    } while (GAGetCurEpoch(ga) < nbMaxEpoch && ev > PBMATH_EPSILON);
    sumEv += ev;
    if (iRun == 0 || bestEv > ev)
      bestEv = ev;
    if (iRun == 0 || maxEv < ev)
      maxEv = ev;
    avgEv = sumEv / (float)(iRun + 1);
    printf("best: %f, worst: %f, avg: %f, ktevent: %lu\n", 
      bestEv, maxEv, avgEv, ga->_nbKTEvent);
    GenAlgFree(&ga);
  }
  avgEv = sumEv / (float)nbRun;
  printf("in %d runs, %lu epochs, best: %f, worst: %f, avg: %f\n", 
    nbRun, nbMaxEpoch, bestEv, maxEv, avgEv);
  printf("UnitTestGenAlgPerf OK\n");
}

void UnitTestGenAlgHistory() {
  srandom(0);
  int lengthAdnF = 4;
  int lengthAdnI = lengthAdnF;
  GenAlg* ga = GenAlgCreate(8, 4, lengthAdnF, lengthAdnI);
  VecFloat2D boundsF = VecFloatCreateStatic2D();
  VecLong2D boundsI = VecLongCreateStatic2D();
  VecSet(&boundsF, 0, -1.0); VecSet(&boundsF, 1, 1.0);
  VecSet(&boundsI, 0, 0); VecSet(&boundsI, 1, 4);
  for (int i = lengthAdnF; i--;) {
    GASetBoundsAdnFloat(ga, i, &boundsF);
    GASetBoundsAdnInt(ga, i, &boundsI);
  }
  GASetFlagHistory(ga, true);
  GASetHistoryPath(ga, "./UnitTestGenAlgHistory.json");
  GAInit(ga);
  GASetNbMinAdn(ga, 8);
  GASetNbMaxAdn(ga, 16);
  float best = 1.0;
  do {
    for (int iEnt = GAGetNbAdns(ga); iEnt--;)
      if (GAAdnIsNew(GAAdn(ga, iEnt)))
        GASetAdnValue(ga, GAAdn(ga, iEnt), 
          -1.0 * evaluate(GAAdnAdnF(GAAdn(ga, iEnt)), 
          GAAdnAdnI(GAAdn(ga, iEnt))));
    GAStep(ga);
    // Display info if there is improvment
    float ev = evaluate(GABestAdnF(ga), GABestAdnI(ga));
    if (best - ev > PBMATH_EPSILON) {
      best = ev;
    }
  } while (GAGetCurEpoch(ga) < 10 && best > PBMATH_EPSILON);
  // Save the history
  bool ret = GASaveHistory(ga);
  if (ret == false) {
    GenAlgErr->_type = PBErrTypeNullPointer;
    sprintf(GenAlgErr->_msg, "Couldn't save the history");
    PBErrCatch(GenAlgErr);
  }

  GAHistory history = GAHistoryCreateStatic();
  FILE* stream = fopen(GAGetHistoryPath(ga), "r");
  ret = GAHistoryLoad(&history, stream);
  if (ret == false) {
    GenAlgErr->_type = PBErrTypeUnitTestFailed;
    sprintf(GenAlgErr->_msg, "GAHistoryLoad failed");
    PBErrCatch(GenAlgErr);
  }
  fclose(stream);
  GSetIterForward iterA = 
    GSetIterForwardCreateStatic(&(ga->_history._genealogy));
  GSetIterForward iterB = 
    GSetIterForwardCreateStatic(&(history._genealogy));
  do {
    GAHistoryBirth* birthA = GSetIterGet(&iterA);
    GAHistoryBirth* birthB = GSetIterGet(&iterB);
    if (birthA->_epoch != birthB->_epoch &&
      birthA->_idParents[0] != birthB->_idParents[0] &&
      birthA->_idParents[1] != birthB->_idParents[1] &&
      birthA->_idChild != birthB->_idChild) {
      GenAlgErr->_type = PBErrTypeUnitTestFailed;
      sprintf(GenAlgErr->_msg, "GAHistoryLoad/Save failed");
      PBErrCatch(GenAlgErr);
    }
    
  } while (GSetIterStep(&iterA) && GSetIterStep(&iterB));
  GAHistoryFree(&history);
  GenAlgFree(&ga);
  printf("UnitTestGenAlgHistory OK\n");
}

void UnitTestGenAlg() {
  UnitTestGenAlgCreateFree();
  UnitTestGenAlgGetSet();
  UnitTestGenAlgInit();
  UnitTestGenAlgPrint();
  UnitTestGenAlgGetDiversity();
  UnitTestGenAlgStep();
  UnitTestGenAlgLoadSave();
  UnitTestGenAlgTest();
  UnitTestGenAlgPerf();
  UnitTestGenAlgHistory();
  printf("UnitTestGenAlg OK\n");
}

void UnitTestAll() {
  UnitTestGenAlgAdn();
  UnitTestGenAlg();
  printf("UnitTestAll OK\n");
}

int main() {
  UnitTestAll();
  // Return success code
  return 0;
}

