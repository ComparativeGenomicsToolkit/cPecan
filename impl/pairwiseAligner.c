/*
 * pairwiseAligner.c
 *
 *  Created on: 1 Mar 2012
 *      Author: benedictpaten
 */
#include "safesort.h"

//This is being included to make popen work!
#ifndef _XOPEN_SOURCE
# define _XOPEN_SOURCE 500
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include "bioioC.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "pairwiseAlignment.h"

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#endif

///////////////////////////////////
///////////////////////////////////
//Diagonal
//
//Structure for working with x-y diagonal of dp matrix
///////////////////////////////////
///////////////////////////////////

const char *PAIRWISE_ALIGNMENT_EXCEPTION_ID = "PAIRWISE_ALIGNMENT_EXCEPTION";

Diagonal diagonal_construct(int64_t xay, int64_t xmyL, int64_t xmyR) {
    if ((xay + xmyL) % 2 != 0 || (xay + xmyR) % 2 != 0 || xmyL > xmyR) {
        stThrowNew(PAIRWISE_ALIGNMENT_EXCEPTION_ID,
                   "Attempt to create diagonal with invalid coordinates: xay %" PRIi64 " xmyL %" PRIi64 " xmyR %" PRIi64 "",
                xay, xmyL, xmyR);
    }
    Diagonal diagonal;
    diagonal.xay = xay;
    diagonal.xmyL = xmyL;
    diagonal.xmyR = xmyR;
    assert(xmyL <= xmyR);
    assert(xay >= 0);
    return diagonal;
}

inline int64_t diagonal_getXay(Diagonal diagonal) {
    return diagonal.xay;
}

inline int64_t diagonal_getMinXmy(Diagonal diagonal) {
    return diagonal.xmyL;
}

inline int64_t diagonal_getMaxXmy(Diagonal diagonal) {
    return diagonal.xmyR;
}

inline int64_t diagonal_getWidth(Diagonal diagonal) {
    return (diagonal.xmyR - diagonal.xmyL) / 2 + 1;
}

inline int64_t diagonal_getXCoordinate(int64_t xay, int64_t xmy) {
    assert((xay + xmy) % 2 == 0);
    return (xay + xmy) / 2;
}

inline int64_t diagonal_equals(Diagonal diagonal1, Diagonal diagonal2) {
    return diagonal1.xay == diagonal2.xay && diagonal1.xmyL == diagonal2.xmyL && diagonal1.xmyR == diagonal2.xmyR;
}

inline int64_t diagonal_getYCoordinate(int64_t xay, int64_t xmy) {
    assert((xay - xmy) % 2 == 0);
    return (xay - xmy) / 2;
}

inline char *diagonal_getString(Diagonal diagonal) {
    return stString_print("Diagonal, xay: %" PRIi64 " xmyL %" PRIi64 ", xmyR: %" PRIi64 "", diagonal_getXay(diagonal),
            diagonal_getMinXmy(diagonal), diagonal_getMaxXmy(diagonal));
}

///////////////////////////////////
///////////////////////////////////
//Band Iterator
//
//Iterator for walking along x+y diagonals in banded fashion
//(using a set of anchor constraints)
///////////////////////////////////
///////////////////////////////////

struct _band {
    Diagonal *diagonals;
    int64_t lXalY;
};

static int64_t band_avoidOffByOne(int64_t xay, int64_t xmy) {
    return (xay + xmy) % 2 == 0 ? xmy : xmy + 1;
}

static void band_setCurrentDiagonalP(int64_t *xmy, int64_t i, int64_t j, int64_t k) {
    if (i < j) {
        *xmy += (int64_t) (2 * ((int64_t) j - i) * (int64_t) k);
    }
}

static Diagonal band_setCurrentDiagonal(int64_t xay, int64_t xL, int64_t yL, int64_t xU, int64_t yU) {
    int64_t xmyL = xL - yL;
    int64_t xmyR = xU - yU;

    assert(xay >= xL + yU);
    assert(xay <= xU + yL);

    //Avoid in-between undefined x,y coordinate positions when intersecting xay and xmy.
    xmyL = band_avoidOffByOne(xay, xmyL);
    xmyR = band_avoidOffByOne(xay, xmyR);

    //Bound the xmy coordinates by the xL, yL and xU, yU band boundaries
    band_setCurrentDiagonalP(&xmyL, diagonal_getXCoordinate(xay, xmyL), xL, 1);
    band_setCurrentDiagonalP(&xmyL, yL, diagonal_getYCoordinate(xay, xmyL), 1);
    band_setCurrentDiagonalP(&xmyR, xU, diagonal_getXCoordinate(xay, xmyR), -1);
    band_setCurrentDiagonalP(&xmyR, diagonal_getYCoordinate(xay, xmyR), yU, -1);

    return diagonal_construct(xay, xmyL, xmyR);
}

static int64_t band_boundCoordinate(int64_t z, int64_t lZ) {
    return z < 0 ? 0 : (z > lZ ? lZ : z);
}

Band *band_constructDynamic(stList *anchorPairs, int64_t lX, int64_t lY) {
    //Prerequisities
    assert(lX >= 0);
    assert(lY >= 0);

    Band *band = st_malloc(sizeof(Band));
    band->diagonals = st_malloc(sizeof(Diagonal) * (lX + lY + 1));
    band->lXalY = lX + lY;

    //Now initialise the diagonals
    int64_t anchorPairIndex = 0;
    int64_t xay = 0;
    int64_t pxay = 0, pxmy = 0;
    int64_t nxay = 0, nxmy = 0;
    int64_t xL = 0, yL = 0, xU = 0, yU = 0, expansion = 0;

    while (xay <= band->lXalY) {
        band->diagonals[xay] = band_setCurrentDiagonal(xay, xL, yL, xU, yU);
        if (nxay == xay++) {
            //The previous diagonals become the next
            pxay = nxay;
            pxmy = nxmy;

            int64_t x = lX, y = lY;
            if (anchorPairIndex < stList_length(anchorPairs)) {
                stIntTuple *anchorPair = stList_get(anchorPairs, anchorPairIndex++);
                x = stIntTuple_get(anchorPair, 0) + 1; //Plus ones, because matrix coordinates are +1 the sequence ones
                y = stIntTuple_get(anchorPair, 1) + 1;
                expansion = stIntTuple_get(anchorPair, 2);

                //Check the anchor pairs
                assert(x > diagonal_getXCoordinate(pxay, pxmy));
                assert(y > diagonal_getYCoordinate(pxay, pxmy));
                assert(x <= lX);
                assert(y <= lY);
                assert(x > 0);
                assert(y > 0);
                assert(expansion >= 0);
                assert(expansion % 2 == 0);
            }

            nxay = x + y;
            nxmy = x - y;

            //Now call to set the lower and upper x,y coordinates
            xL = band_boundCoordinate(diagonal_getXCoordinate(pxay, pxmy - expansion), lX);
            yL = band_boundCoordinate(diagonal_getYCoordinate(nxay, nxmy - expansion), lY);
            xU = band_boundCoordinate(diagonal_getXCoordinate(nxay, nxmy + expansion), lX);
            yU = band_boundCoordinate(diagonal_getYCoordinate(pxay, pxmy + expansion), lY);
        }
    }

    return band;
}

Band *band_construct(stList *anchorPairs, int64_t lX, int64_t lY, int64_t expansion) {
    //Prerequisities
    assert(lX >= 0);
    assert(lY >= 0);
    assert(expansion % 2 == 0);

    Band *band = st_malloc(sizeof(Band));
    band->diagonals = st_malloc(sizeof(Diagonal) * (lX + lY + 1));
    band->lXalY = lX + lY;

    //Now initialise the diagonals
    int64_t anchorPairIndex = 0;
    int64_t xay = 0;
    int64_t pxay = 0, pxmy = 0;
    int64_t nxay = 0, nxmy = 0;
    int64_t xL = 0, yL = 0, xU = 0, yU = 0;

    while (xay <= band->lXalY) {
        band->diagonals[xay] = band_setCurrentDiagonal(xay, xL, yL, xU, yU);
        if (nxay == xay++) {
            //The previous diagonals become the next
            pxay = nxay;
            pxmy = nxmy;

            int64_t x = lX, y = lY;
            if (anchorPairIndex < stList_length(anchorPairs)) {
                stIntTuple *anchorPair = stList_get(anchorPairs, anchorPairIndex++);
                x = stIntTuple_get(anchorPair, 0) + 1; //Plus ones, because matrix coordinates are +1 the sequence ones
                y = stIntTuple_get(anchorPair, 1) + 1;

                //Check the anchor pairs
                assert(x > diagonal_getXCoordinate(pxay, pxmy));
                assert(y > diagonal_getYCoordinate(pxay, pxmy));
                assert(x <= lX);
                assert(y <= lY);
                assert(x > 0);
                assert(y > 0);
            }

            nxay = x + y;
            nxmy = x - y;

            //Now call to set the lower and upper x,y coordinates
            xL = band_boundCoordinate(diagonal_getXCoordinate(pxay, pxmy - expansion), lX);
            yL = band_boundCoordinate(diagonal_getYCoordinate(nxay, nxmy - expansion), lY);
            xU = band_boundCoordinate(diagonal_getXCoordinate(nxay, nxmy + expansion), lX);
            yU = band_boundCoordinate(diagonal_getYCoordinate(pxay, pxmy + expansion), lY);
        }
    }

    return band;
}

void band_destruct(Band *band) {
    free(band->diagonals);
    free(band);
}

struct _bandIterator {
    Band *band;
    int64_t index;
};

BandIterator *bandIterator_construct(Band *band) {
    BandIterator *bandIterator = st_malloc(sizeof(BandIterator));
    bandIterator->band = band;
    bandIterator->index = 0;
    return bandIterator;
}

BandIterator *bandIterator_clone(BandIterator *bandIterator) {
    BandIterator *bandIterator2 = st_malloc(sizeof(BandIterator));
    memcpy(bandIterator2, bandIterator, sizeof(BandIterator));
    return bandIterator2;
}

void bandIterator_destruct(BandIterator *bandIterator) {
    free(bandIterator);
}

Diagonal bandIterator_getNext(BandIterator *bandIterator) {
    Diagonal diagonal = bandIterator->band->diagonals[
            bandIterator->index > bandIterator->band->lXalY ? bandIterator->band->lXalY : bandIterator->index];
    if (bandIterator->index <= bandIterator->band->lXalY) {
        bandIterator->index++;
    }
    return diagonal;
}

Diagonal bandIterator_getPrevious(BandIterator *bandIterator) {
    if (bandIterator->index > 0) {
        bandIterator->index--;
    }
    return bandIterator->band->diagonals[bandIterator->index];
}

///////////////////////////////////
///////////////////////////////////
//Log Add functions
//
//Interpolation function for doing log add
///////////////////////////////////
///////////////////////////////////

#define logUnderflowThreshold 7.5
#define posteriorMatchThreshold 0.01

static inline double lookup(double x) {
    //return log (exp (x) + 1);
    assert(x >= 0.00f);
    assert(x <= logUnderflowThreshold);
    if (x <= 1.00f)
        return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
    if (x <= 2.50f)
        return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
    if (x <= 4.50f)
        return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
    return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

double logAdd(double x, double y) {
    if (x < y)
        return (x == LOG_ZERO || y - x >= logUnderflowThreshold) ? y : lookup(y - x) + x;
    return (y == LOG_ZERO || x - y >= logUnderflowThreshold) ? x : lookup(x - y) + y;
}

///////////////////////////////////
///////////////////////////////////
//Symbols
//
//Emissions probs/functions to convert to symbol sequence
///////////////////////////////////
///////////////////////////////////

Symbol symbol_convertCharToSymbol(char i) {
    switch (i) {
        case 'A':
        case 'a':
            return a;
        case 'C':
        case 'c':
            return c;
        case 'G':
        case 'g':
            return g;
        case 'T':
        case 't':
            return t;
        default:
            return n;
    }
}

char symbol_convertSymbolToChar(Symbol i) {
    switch (i) {
        case a:
            return 'A';
        case c:
            return 'C';
        case g:
            return 'G';
        case t:
            return 'T';
        default:
            return 'N';
    }
}

Symbol *symbol_convertStringToSymbols(const char *s, int64_t sL) {
    assert(sL >= 0);
    assert(strlen(s) == sL);
    Symbol *cS = st_malloc(sL * sizeof(Symbol));
    for (int64_t i = 0; i < sL; i++) {
        cS[i] = symbol_convertCharToSymbol(s[i]);
    }
    return cS;
}

SymbolString symbolString_construct(const char *sequence, int64_t length) {
    SymbolString symbolString;
    symbolString.sequence = symbol_convertStringToSymbols(sequence, length);
    symbolString.length = length;
    return symbolString;
}

void symbolString_destruct(SymbolString s) {
    free(s.sequence);
}

///////////////////////////////////
///////////////////////////////////
//Cell calculations
//
//A cell is a set of states associated with an x, y coordinate.
//These functions do the forward/backward calculations for the pairwise
//alignment model.
///////////////////////////////////
///////////////////////////////////

static inline void doTransitionForward(double *fromCells, double *toCells, int64_t from, int64_t to, double eP,
                                       double tP, void *extraArgs) {
    toCells[to] = logAdd(toCells[to], fromCells[from] + (eP + tP));
}

void cell_calculateForward(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
                           void *extraArgs) {
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY, doTransitionForward, extraArgs);
}

static inline void doTransitionBackward(double *fromCells, double *toCells, int64_t from, int64_t to, double eP,
                                        double tP, void *extraArgs) {
    fromCells[from] = logAdd(fromCells[from], toCells[to] + (eP + tP));
}

void cell_calculateBackward(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
                            void *extraArgs) {
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY, doTransitionBackward, extraArgs);
}

double cell_dotProduct(double *cell1, double *cell2, int64_t stateNumber) {
    double totalProb = cell1[0] + cell2[0];
    for (int64_t i = 1; i < stateNumber; i++) {
        totalProb = logAdd(totalProb, cell1[i] + cell2[i]);
    }
    return totalProb;
}

double cell_dotProduct2(double *cell, StateMachine *sM, double (*getStateValue)(StateMachine *, int64_t)) {
    double totalProb = cell[0] + getStateValue(sM, 0);
    for (int64_t i = 1; i < sM->stateNumber; i++) {
        totalProb = logAdd(totalProb, cell[i] + getStateValue(sM, i));
    }
    return totalProb;
}

static inline void updateExpectations(double *fromCells, double *toCells, int64_t from, int64_t to, double eP,
                                      double tP, void *extraArgs) {
    //void *extraArgs2[2] = { &totalProbability, hmmExpectations };
    double totalProbability = *((double *) ((void **) extraArgs)[0]);
    Hmm *hmmExpectations = ((void **) extraArgs)[1];
    Symbol x = *((Symbol *)((void **) extraArgs)[2]);
    Symbol y = *((Symbol *)((void **) extraArgs)[3]);
    //Calculate posterior probability of the transition/emission pair
    double p = exp(fromCells[from] + toCells[to] + (eP + tP) - totalProbability);
    //Add in the expectation of the transition
    hmm_addToTransitionExpectation(hmmExpectations, from, to, p);
    if(x < SYMBOL_NUMBER_NO_N && y < SYMBOL_NUMBER_NO_N) { //Ignore gaps involving Ns.
        hmm_addToEmissionsExpectation(hmmExpectations, to, x, y, p);
    }
}

static void cell_calculateExpectation(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
                                      void *extraArgs) {
    void *extraArgs2[4] = { ((void **)extraArgs)[0], ((void **)extraArgs)[1], &cX, &cY };
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY, updateExpectations, extraArgs2);
}

///////////////////////////////////
///////////////////////////////////
//DpDiagonal
//
//Structure for storing a x-y diagonal of the dp matrix
///////////////////////////////////
///////////////////////////////////

struct _dpDiagonal {
    Diagonal diagonal;
    int64_t stateNumber;
    double *cells;
};

DpDiagonal *dpDiagonal_construct(Diagonal diagonal, int64_t stateNumber) {
    DpDiagonal *dpDiagonal = st_malloc(sizeof(DpDiagonal));
    dpDiagonal->diagonal = diagonal;
    dpDiagonal->stateNumber = stateNumber;
    assert(diagonal_getWidth(diagonal) >= 0);
    dpDiagonal->cells = st_malloc(sizeof(double) * stateNumber * (int64_t) diagonal_getWidth(diagonal));
    return dpDiagonal;
}

DpDiagonal *dpDiagonal_clone(DpDiagonal *diagonal) {
    DpDiagonal *diagonal2 = dpDiagonal_construct(diagonal->diagonal, diagonal->stateNumber);
    memcpy(diagonal2->cells, diagonal->cells, sizeof(double) * diagonal_getWidth(diagonal->diagonal) * diagonal->stateNumber);
    return diagonal2;
}

bool dpDiagonal_equals(DpDiagonal *diagonal1, DpDiagonal *diagonal2) {
    if (!diagonal_equals(diagonal1->diagonal, diagonal2->diagonal)) {
        return 0;
    }
    if(diagonal1->stateNumber != diagonal2->stateNumber) {
        return 0;
    }
    for (int64_t i = 0; i < diagonal_getWidth(diagonal1->diagonal) * diagonal1->stateNumber; i++) {
        if (diagonal1->cells[i] != diagonal2->cells[i]) {
            return 0;
        }
    }
    return 1;
}

void dpDiagonal_destruct(DpDiagonal *dpDiagonal) {
    free(dpDiagonal->cells);
    free(dpDiagonal);
}

double *dpDiagonal_getCell(DpDiagonal *dpDiagonal, int64_t xmy) {
    if (xmy < dpDiagonal->diagonal.xmyL || xmy > dpDiagonal->diagonal.xmyR) {
        return NULL;
    }
    assert((diagonal_getXay(dpDiagonal->diagonal) + xmy) % 2 == 0);
    return &dpDiagonal->cells[((xmy - dpDiagonal->diagonal.xmyL) / 2) * dpDiagonal->stateNumber];
}

void dpDiagonal_zeroValues(DpDiagonal *diagonal) {
    for (int64_t i = 0; i < diagonal_getWidth(diagonal->diagonal) * diagonal->stateNumber; i++) {
        diagonal->cells[i] = LOG_ZERO;
    }
}

void dpDiagonal_initialiseValues(DpDiagonal *diagonal, StateMachine *sM, double (*getStateValue)(StateMachine *, int64_t)) {
    for (int64_t i = diagonal_getMinXmy(diagonal->diagonal); i <= diagonal_getMaxXmy(diagonal->diagonal); i += 2) {
        double *cell = dpDiagonal_getCell(diagonal, i);
        assert(cell != NULL);
        for (int64_t j = 0; j < diagonal->stateNumber; j++) {
            cell[j] = getStateValue(sM, j);
        }
    }
}

double dpDiagonal_dotProduct(DpDiagonal *diagonal1, DpDiagonal *diagonal2) {
    double totalProbability = LOG_ZERO;
    Diagonal diagonal = diagonal1->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        totalProbability = logAdd(totalProbability,
                                  cell_dotProduct(dpDiagonal_getCell(diagonal1, xmy), dpDiagonal_getCell(diagonal2, xmy), diagonal1->stateNumber));
        xmy += 2;
    }
    return totalProbability;
}

///////////////////////////////////
///////////////////////////////////
//DpMatrix
//
//Structure for storing dp-matrix
///////////////////////////////////
///////////////////////////////////

struct _dpMatrix {
    DpDiagonal **diagonals;
    int64_t diagonalNumber;
    int64_t activeDiagonals;
    int64_t stateNumber;
};

DpMatrix *dpMatrix_construct(int64_t diagonalNumber, int64_t stateNumber) {
    assert(diagonalNumber >= 0);
    DpMatrix *dpMatrix = st_malloc(sizeof(DpMatrix));
    dpMatrix->diagonalNumber = diagonalNumber;
    dpMatrix->diagonals = st_calloc(dpMatrix->diagonalNumber + 1, sizeof(DpDiagonal *));
    dpMatrix->activeDiagonals = 0;
    dpMatrix->stateNumber = stateNumber;
    return dpMatrix;
}

void dpMatrix_destruct(DpMatrix *dpMatrix) {
    assert(dpMatrix->activeDiagonals == 0);
    free(dpMatrix->diagonals);
    free(dpMatrix);
}

DpDiagonal *dpMatrix_getDiagonal(DpMatrix *dpMatrix, int64_t xay) {
    if (xay < 0 || xay > dpMatrix->diagonalNumber) {
        return NULL;
    }
    return dpMatrix->diagonals[xay];
}

int64_t dpMatrix_getActiveDiagonalNumber(DpMatrix *dpMatrix) {
    return dpMatrix->activeDiagonals;
}

DpDiagonal *dpMatrix_createDiagonal(DpMatrix *dpMatrix, Diagonal diagonal) {
    assert(diagonal.xay >= 0);
    assert(diagonal.xay <= dpMatrix->diagonalNumber);
    assert(dpMatrix_getDiagonal(dpMatrix, diagonal.xay) == NULL);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal, dpMatrix->stateNumber);
    dpMatrix->diagonals[diagonal_getXay(diagonal)] = dpDiagonal;
    dpMatrix->activeDiagonals++;
    return dpDiagonal;
}

void dpMatrix_deleteDiagonal(DpMatrix *dpMatrix, int64_t xay) {
    assert(xay >= 0);
    assert(xay <= dpMatrix->diagonalNumber);
    if (dpMatrix->diagonals[xay] != NULL) {
        dpMatrix->activeDiagonals--;
        assert(dpMatrix->activeDiagonals >= 0);
        dpDiagonal_destruct(dpMatrix->diagonals[xay]);
        dpMatrix->diagonals[xay] = NULL;
    }
}

///////////////////////////////////
///////////////////////////////////
//Diagonal DP Calculations
//
//Functions which do forward/backward/posterior calculations
//between diagonal rows of a dp-matrix
///////////////////////////////////
///////////////////////////////////

static Symbol getXCharacter(const SymbolString sX, int64_t xay, int64_t xmy) {
    int64_t x = diagonal_getXCoordinate(xay, xmy);
    assert(x >= 0 && x <= sX.length);
    return x > 0 ? sX.sequence[x - 1] : n;
}

static Symbol getYCharacter(const SymbolString sY, int64_t xay, int64_t xmy) {
    int64_t y = diagonal_getYCoordinate(xay, xmy);
    assert(y >= 0 && y <= sY.length);
    return y > 0 ? sY.sequence[y - 1] : n;
}

static void diagonalCalculation(StateMachine *sM, DpDiagonal *dpDiagonal, DpDiagonal *dpDiagonalM1, DpDiagonal *dpDiagonalM2,
                                const SymbolString sX, const SymbolString sY,
                                void (*cellCalculation)(StateMachine *, double *, double *, double *, double *, Symbol, Symbol, void *), void *extraArgs) {
    Diagonal diagonal = dpDiagonal->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        Symbol x = getXCharacter(sX, diagonal_getXay(diagonal), xmy);
        Symbol y = getYCharacter(sY, diagonal_getXay(diagonal), xmy);
        double *current = dpDiagonal_getCell(dpDiagonal, xmy);
        double *lower = dpDiagonalM1 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM1, xmy - 1);
        double *middle = dpDiagonalM2 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM2, xmy);
        double *upper = dpDiagonalM1 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM1, xmy + 1);
        cellCalculation(sM, current, lower, middle, upper, x, y, extraArgs);
        xmy += 2;
    }
}

void diagonalCalculationForward(StateMachine *sM, int64_t xay, DpMatrix *dpMatrix, const SymbolString sX, const SymbolString sY) {
    diagonalCalculation(sM, dpMatrix_getDiagonal(dpMatrix, xay), dpMatrix_getDiagonal(dpMatrix, xay - 1),
                        dpMatrix_getDiagonal(dpMatrix, xay - 2), sX, sY, cell_calculateForward, NULL);
}

void diagonalCalculationBackward(StateMachine *sM, int64_t xay, DpMatrix *dpMatrix, const SymbolString sX, const SymbolString sY) {
    diagonalCalculation(sM, dpMatrix_getDiagonal(dpMatrix, xay), dpMatrix_getDiagonal(dpMatrix, xay - 1),
                        dpMatrix_getDiagonal(dpMatrix, xay - 2), sX, sY, cell_calculateBackward, NULL);
}

double diagonalCalculationTotalProbability(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
                                           const SymbolString sX, const SymbolString sY) {
    //Get the forward and backward diagonals
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    double totalProbability = dpDiagonal_dotProduct(forwardDiagonal, backDiagonal);
    //Now calculate the contribution of matches through xay.
    forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay - 1);
    backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay + 1);
    if (backDiagonal != NULL && forwardDiagonal != NULL) {
        DpDiagonal *matchDiagonal = dpDiagonal_clone(backDiagonal);
        dpDiagonal_zeroValues(matchDiagonal);
        diagonalCalculation(sM, matchDiagonal, NULL, forwardDiagonal, sX, sY, cell_calculateForward, NULL);
        totalProbability = logAdd(totalProbability, dpDiagonal_dotProduct(matchDiagonal, backDiagonal));
        dpDiagonal_destruct(matchDiagonal);
    }
    return totalProbability;
}

void addPosteriorProb(int64_t x, int64_t y, double posteriorProbability, stList *posteriorProbs, PairwiseAlignmentParameters *p) {
    if (posteriorProbability >= p->threshold) {
        if (posteriorProbability > 1.0) {
            posteriorProbability = 1.0;
        }
        posteriorProbability = floor(posteriorProbability * PAIR_ALIGNMENT_PROB_1);

        stList_append(posteriorProbs, stIntTuple_construct3((int64_t) posteriorProbability, x - 1, y - 1));
    }
}

void diagonalCalculationPosteriorMatchProbs(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
                                            const SymbolString sX, const SymbolString sY, double totalProbability, PairwiseAlignmentParameters *p,
                                            void *extraArgs) {
    assert(p->threshold >= 0.0);
    assert(p->threshold <= 1.0);
    stList *alignedPairs = ((void **) extraArgs)[0];
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    Diagonal diagonal = forwardDiagonal->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);
    //Walk over the cells computing the posteriors
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        int64_t x = diagonal_getXCoordinate(diagonal_getXay(diagonal), xmy);
        int64_t y = diagonal_getYCoordinate(diagonal_getXay(diagonal), xmy);
        if (x > 0 && y > 0) {
            double *cellForward = dpDiagonal_getCell(forwardDiagonal, xmy);
            double *cellBackward = dpDiagonal_getCell(backDiagonal, xmy);
            double posteriorProbability = exp(
                    (cellForward[sM->matchState] + cellBackward[sM->matchState]) - totalProbability);
            addPosteriorProb(x, y, posteriorProbability, alignedPairs, p);
        }
        xmy += 2;
    }
}

void diagonalCalculationPosteriorProbs(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
                                       const SymbolString sX, const SymbolString sY, double totalProbability, PairwiseAlignmentParameters *p,
                                       void *extraArgs) {
    assert(p->threshold >= 0.0);
    assert(p->threshold <= 1.0);

    stList *alignedPairs = ((void **) extraArgs)[0];
    stList *gapXPairs = ((void **) extraArgs)[2];
    stList *gapYPairs = ((void **) extraArgs)[4];

    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    Diagonal diagonal = forwardDiagonal->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);
    //Walk over the cells computing the posteriors
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        int64_t x = diagonal_getXCoordinate(diagonal_getXay(diagonal), xmy);
        int64_t y = diagonal_getYCoordinate(diagonal_getXay(diagonal), xmy);

        double *cellForward = dpDiagonal_getCell(forwardDiagonal, xmy);
        double *cellBackward = dpDiagonal_getCell(backDiagonal, xmy);
        if (x > 0 && y > 0) {
            // Posterior match prob
            double posteriorProbability = exp(
                    (cellForward[sM->matchState] + cellBackward[sM->matchState]) - totalProbability);
            addPosteriorProb(x, y, posteriorProbability, alignedPairs, p);
        }

        if(x > 0) {
            double posteriorProbability = exp(
                    (cellForward[sM->gapXState] + cellBackward[sM->gapXState]) - totalProbability);
            addPosteriorProb(x, y, posteriorProbability, gapXPairs, p);
        }

        if(y > 0) {
            double posteriorProbability = exp(
                    (cellForward[sM->gapYState] + cellBackward[sM->gapYState]) - totalProbability);
            addPosteriorProb(x, y, posteriorProbability, gapYPairs, p);
        }

        xmy += 2;
    }
}

static void diagonalCalculationExpectations(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
                                            const SymbolString sX, const SymbolString sY, double totalProbability, PairwiseAlignmentParameters *p,
                                            void *extraArgs) {
    /*
     * Updates the expectations of the transitions/emissions for the given diagonal.
     */
    Hmm *hmmExpectations = extraArgs;
    void *extraArgs2[2] = { &totalProbability, hmmExpectations };
    hmmExpectations->likelihood += totalProbability; //We do this once per diagonal, which is a hack, rather than for the whole matrix. The correction factor is approximately 1/number of diagonals.
    diagonalCalculation(sM, dpMatrix_getDiagonal(backwardDpMatrix, xay), dpMatrix_getDiagonal(forwardDpMatrix, xay - 1),
                        dpMatrix_getDiagonal(forwardDpMatrix, xay - 2), sX, sY, cell_calculateExpectation, extraArgs2);
}

///////////////////////////////////
///////////////////////////////////
//Fast banded forward-backward
//
//The log space recursion below goes through three levels of function pointer per cell (the
//diagonal, the cell and the transition) and an approximate log-add per transition.  For the
//posterior probabilities -- which is all bar asks for -- the same banded, chunked recursion is
//computed here in probability space instead: each diagonal is rescaled so its largest value is
//one and the scale is kept as a log, the state machine is flattened into arrays once per call,
//and each state of a diagonal is a plain loop over its cells.  The band, the traceback points
//and which diagonals get posteriors in which chunk are exactly those of the log space version;
//only the arithmetic differs, and this is the more accurate of the two, since logAdd is a
//piecewise cubic fit that also drops any term under e^-7.5 of the larger.
///////////////////////////////////
///////////////////////////////////

#if defined(__SSE2__)
#include <xmmintrin.h>
#endif
#if defined(__AVX2__)
#include <immintrin.h>
#endif

#define PAIR_HMM_MAX_STATES 5

enum { PAIR_HMM_MATCH = 0, PAIR_HMM_GAP_X = 1, PAIR_HMM_GAP_Y = 2 };

typedef struct _pairHmmTerms { // the nonzero transitions into (or out of) a state
    int64_t number;
    int64_t state[PAIR_HMM_MAX_STATES];
    double t[PAIR_HMM_MAX_STATES];
} PairHmmTerms;

typedef struct _pairHmm {
    int64_t stateNumber;
    int64_t type[PAIR_HMM_MAX_STATES]; // what a state emits: a match, a gap in y (x only) or a gap in x (y only)
    double t[PAIR_HMM_MAX_STATES][PAIR_HMM_MAX_STATES]; // transition probabilities, [from][to]
    double eMatch[SYMBOL_NUMBER * SYMBOL_NUMBER], eGapX[SYMBOL_NUMBER], eGapY[SYMBOL_NUMBER];
    double start[PAIR_HMM_MAX_STATES], raggedStart[PAIR_HMM_MAX_STATES];
    double end[PAIR_HMM_MAX_STATES], raggedEnd[PAIR_HMM_MAX_STATES];
    PairHmmTerms into[PAIR_HMM_MAX_STATES], outOf[PAIR_HMM_MAX_STATES]; // t, sparse, by destination and by source
    bool ok; // false if the state machine does not fit this form, in which case the log space code is used
} PairHmm;

typedef struct {
    PairHmm *hmm;
    double *lower, *middle, *upper;
    int64_t x, y;
} PairHmmExtraction;

static void pairHmm_setOnce(double *slot, double value, PairHmm *hmm) {
    if (*slot < 0.0) {
        *slot = value;
    } else if (*slot != value) {
        hmm->ok = 0;
    }
}

/*
 * A doTransition that records the transition instead of computing with it.  Which neighbour
 * the transition reads from says what the destination state emits.
 */
static void pairHmm_recordTransition(double *fromCells, double *toCells, int64_t from, int64_t to, double eP, double tP,
                                     void *extraArgs) {
    PairHmmExtraction *e = extraArgs;
    PairHmm *hmm = e->hmm;
    if (from < 0 || from >= hmm->stateNumber || to < 0 || to >= hmm->stateNumber) {
        hmm->ok = 0;
        return;
    }
    int64_t type = fromCells == e->lower ? PAIR_HMM_GAP_X : (fromCells == e->middle ? PAIR_HMM_MATCH : PAIR_HMM_GAP_Y);
    if (hmm->type[to] == -1) {
        hmm->type[to] = type;
    } else if (hmm->type[to] != type) { // the recursion needs each state entered from one neighbour only
        hmm->ok = 0;
    }
    pairHmm_setOnce(&hmm->t[from][to], exp(tP), hmm); // must not depend on the symbols
    double *emission = type == PAIR_HMM_MATCH ? &hmm->eMatch[e->x * SYMBOL_NUMBER + e->y] :
                       (type == PAIR_HMM_GAP_X ? &hmm->eGapX[e->x] : &hmm->eGapY[e->y]);
    pairHmm_setOnce(emission, exp(eP), hmm); // must not depend on the state transitioned from
}

/*
 * Flattens the state machine by running its cell calculation once for every pair of symbols
 * and recording the transitions, rather than by reading its fields, so it works for any state
 * machine whose states each emit in one way and whose transitions do not depend on the symbols.
 */
static void pairHmm_construct(PairHmm *hmm, StateMachine *sM) {
    hmm->ok = sM->stateNumber > 0 && sM->stateNumber <= PAIR_HMM_MAX_STATES;
    if (!hmm->ok) {
        return;
    }
    hmm->stateNumber = sM->stateNumber;
    for (int64_t s = 0; s < PAIR_HMM_MAX_STATES; s++) {
        hmm->type[s] = -1;
        for (int64_t s2 = 0; s2 < PAIR_HMM_MAX_STATES; s2++) {
            hmm->t[s][s2] = -1.0;
        }
    }
    for (int64_t i = 0; i < SYMBOL_NUMBER * SYMBOL_NUMBER; i++) {
        hmm->eMatch[i] = -1.0;
    }
    for (int64_t i = 0; i < SYMBOL_NUMBER; i++) {
        hmm->eGapX[i] = -1.0;
        hmm->eGapY[i] = -1.0;
    }
    double current[PAIR_HMM_MAX_STATES], lower[PAIR_HMM_MAX_STATES], middle[PAIR_HMM_MAX_STATES], upper[PAIR_HMM_MAX_STATES];
    PairHmmExtraction e = { hmm, lower, middle, upper, 0, 0 };
    for (e.x = 0; e.x < SYMBOL_NUMBER; e.x++) {
        for (e.y = 0; e.y < SYMBOL_NUMBER; e.y++) {
            sM->cellCalculate(sM, current, lower, middle, upper, (Symbol) e.x, (Symbol) e.y, pairHmm_recordTransition, &e);
        }
    }
    for (int64_t s = 0; s < hmm->stateNumber; s++) {
        if (hmm->type[s] == -1) {
            hmm->ok = 0;
        }
        for (int64_t s2 = 0; s2 < hmm->stateNumber; s2++) {
            if (hmm->t[s][s2] < 0.0) {
                hmm->t[s][s2] = 0.0; // no such transition
            }
        }
        hmm->start[s] = exp(sM->startStateProb(sM, s));
        hmm->raggedStart[s] = exp(sM->raggedStartStateProb(sM, s));
        hmm->end[s] = exp(sM->endStateProb(sM, s));
        hmm->raggedEnd[s] = exp(sM->raggedEndStateProb(sM, s));
    }
    for (int64_t i = 0; i < SYMBOL_NUMBER * SYMBOL_NUMBER; i++) {
        hmm->ok = hmm->ok && hmm->eMatch[i] >= 0.0;
    }
    for (int64_t i = 0; i < SYMBOL_NUMBER; i++) {
        hmm->ok = hmm->ok && hmm->eGapX[i] >= 0.0 && hmm->eGapY[i] >= 0.0;
    }
    for (int64_t s = 0; s < hmm->stateNumber; s++) {
        hmm->into[s].number = 0;
        hmm->outOf[s].number = 0;
    }
    for (int64_t from = 0; from < hmm->stateNumber; from++) { // in increasing order of state, both ways
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            double t = hmm->t[from][to];
            if (t > 0.0) {
                PairHmmTerms *i = &hmm->into[to], *o = &hmm->outOf[from];
                i->state[i->number] = from;
                i->t[i->number++] = t;
                o->state[o->number] = to;
                o->t[o->number++] = t;
            }
        }
    }}

/*
 * A diagonal of the scaled dp matrix.  Cells are stored state-major, so the cell for state s
 * at xmy = xmyL + 2i is cells[s * width + i], and the true value is that times exp(logScale).
 */
enum { FAST_STORED, FAST_RING, FAST_BLOCK }; // where a forward diagonal's cells are, see FastForward

typedef struct _fastDiagonal {
    int64_t xmyL, width;
    int64_t where, offset; // which buffer holds the cells, and where in it (for the ring, which of its three)
    double logScale;
} FastDiagonal;

/*
 * The dp between one traceback and the next.  A chunk computes forward values for diagonals
 * (start, end], walks the backward values from end down to tracedBackTo, takes posteriors up
 * to tracedBackFrom, and hands the forward values from tracedBackFrom on to the next chunk.
 *
 * A large chunk does not keep all of its forward diagonals: only the ones it hands on and a
 * pair before every blockSize'th, from which the walk back recomputes each block of them when
 * it gets there.  That is a second forward pass over the chunk, but it holds a few hundred
 * diagonals rather than thousands, and what the walk reads is still in cache.  Keeping them
 * all, an unanchored 3000x3000 matrix holds 360MB per thread and runs at the speed of memory.
 */
typedef struct _fastChunk {
    int64_t start, end, tracedBackTo, tracedBackFrom;
    int64_t blockSize; // 0 if every forward diagonal is kept
} FastChunk;

#define FAST_CHECKPOINT_BYTES (16 * 1024 * 1024) // chunks whose forward values would take more are recomputed in blocks

/*
 * Whether the chunk keeps diagonal xay's forward values, rather than computing them only to feed
 * the next two.  Blocks start at start + 1 + j * blockSize.
 */
static bool fastChunk_stores(FastChunk *c, int64_t xay) {
    if (c->blockSize == 0 || xay >= c->tracedBackFrom || xay <= c->start) {
        return 1;
    }
    for (int64_t blockStart = xay + 1; blockStart <= xay + 2; blockStart++) { // the pair before a block
        if (blockStart > c->start + 1 && blockStart < c->tracedBackFrom && (blockStart - c->start - 1) % c->blockSize == 0) {
            return 1;
        }
    }
    return 0;
}

/*
 * The forward diagonals the current chunk holds, and the buffers they are held in.
 */
typedef struct _fastForward {
    FastDiagonal *diagonals; // diagonals[xay - base], for every diagonal of the chunk however it is held
    int64_t base, length, capacity;
    double *cells; // FAST_STORED: kept diagonals, contiguous and in order
    int64_t cellsLength, cellsCapacity;
    double *ring[3]; // FAST_RING: diagonals computed only to feed the next two, round robin
    double *block; // FAST_BLOCK: the block the walk back last recomputed
    int64_t blockStart, blockCellsLength, blockCellsCapacity;
} FastForward;

static FastDiagonal *fastForward_get(FastForward *f, int64_t xay) {
    return xay >= f->base && xay < f->base + f->length ? &f->diagonals[xay - f->base] : NULL;
}

static double *fastForward_cells(FastForward *f, FastDiagonal *d) {
    return d->where == FAST_STORED ? f->cells + d->offset : (d->where == FAST_RING ? f->ring[d->offset] : f->block + d->offset);
}

static FastDiagonal *fastForward_add(FastForward *f, Diagonal diagonal, int64_t stateNumber, bool store) {
    assert(f->length < f->capacity); // all the buffers are sized up front, from the band
    int64_t xay = diagonal_getXay(diagonal);
    assert(xay == f->base + f->length);
    FastDiagonal *d = &f->diagonals[f->length++];
    d->xmyL = diagonal_getMinXmy(diagonal);
    d->width = diagonal_getWidth(diagonal);
    d->logScale = 0.0;
    if (store) {
        d->where = FAST_STORED;
        d->offset = f->cellsLength;
        f->cellsLength += d->width * stateNumber;
        assert(f->cellsLength <= f->cellsCapacity);
    } else {
        d->where = FAST_RING;
        d->offset = xay % 3;
    }
    return d;
}

/*
 * Plans the chunks.  Where the tracebacks fall depends only on the band, so the chunks, and so
 * the buffers they need, are known before any dp is done: nothing is grown as the dp goes.
 */
static FastChunk *fastForward_construct(FastForward *f, Band *band, int64_t diagonalNumber, int64_t stateNumber,
                                        PairwiseAlignmentParameters *p, int64_t *chunkNumber) {
    int64_t chunkCapacity = 16;
    FastChunk *chunks = st_malloc(sizeof(FastChunk) * chunkCapacity);
    *chunkNumber = 0;
    int64_t maxDiagonals = 1, maxCells = 1, maxBlockCells = 1, maxWidth = 1;
    int64_t tracedBackTo = 0, previousEnd = 0;
    for (int64_t xay = 0; xay <= diagonalNumber; xay++) {
        int64_t width = diagonal_getWidth(band->diagonals[xay]);
        maxWidth = width > maxWidth ? width : maxWidth;
        bool atEnd = xay == diagonalNumber;
        bool tracebackPoint = xay >= tracedBackTo + p->minDiagsBetweenTraceBack && width <= p->diagonalExpansion * 2 + 1;
        if (xay == 0 || (!atEnd && !tracebackPoint)) {
            continue;
        }
        FastChunk c = { previousEnd, xay, tracedBackTo, xay - (atEnd ? 0 : p->traceBackDiagonals + 1), 0 };
        int64_t cells = 0;
        for (int64_t d = c.tracedBackTo; d <= c.end; d++) {
            cells += diagonal_getWidth(band->diagonals[d]);
        }
        if (cells * stateNumber * (int64_t) sizeof(double) > FAST_CHECKPOINT_BYTES) {
            c.blockSize = (int64_t) sqrt((double) (c.tracedBackFrom - c.start));
            c.blockSize = c.blockSize < 8 ? 8 : c.blockSize;
            cells = 0;
            int64_t blockCells = 0;
            for (int64_t d = c.tracedBackTo; d <= c.end; d++) {
                int64_t w = diagonal_getWidth(band->diagonals[d]);
                cells += fastChunk_stores(&c, d) ? w : 0;
                if (d > c.start && d < c.tracedBackFrom) {
                    blockCells = (d - c.start - 1) % c.blockSize == 0 ? w : blockCells + w;
                    maxBlockCells = blockCells > maxBlockCells ? blockCells : maxBlockCells;
                }
            }
        }
        maxCells = cells > maxCells ? cells : maxCells;
        maxDiagonals = c.end - c.tracedBackTo + 1 > maxDiagonals ? c.end - c.tracedBackTo + 1 : maxDiagonals;
        if (*chunkNumber == chunkCapacity) {
            chunkCapacity *= 2;
            chunks = st_realloc(chunks, sizeof(FastChunk) * chunkCapacity);
        }
        chunks[(*chunkNumber)++] = c;
        tracedBackTo = c.tracedBackFrom;
        previousEnd = xay;
    }
    f->capacity = maxDiagonals + 1;
    f->diagonals = st_malloc(sizeof(FastDiagonal) * f->capacity);
    f->cellsCapacity = maxCells * stateNumber;
    f->cells = st_malloc(sizeof(double) * f->cellsCapacity);
    for (int64_t i = 0; i < 3; i++) {
        f->ring[i] = st_malloc(sizeof(double) * maxWidth * stateNumber);
    }
    f->blockCellsCapacity = maxBlockCells * stateNumber;
    f->block = st_malloc(sizeof(double) * f->blockCellsCapacity);
    f->blockStart = -1;
    f->blockCellsLength = 0;
    f->base = 0;
    f->length = 0;
    f->cellsLength = 0;
    return chunks;
}

static void fastForward_destruct(FastForward *f) {
    free(f->diagonals);
    free(f->cells);
    for (int64_t i = 0; i < 3; i++) {
        free(f->ring[i]);
    }
    free(f->block);
}

/*
 * Drops the diagonals before xay, moving the rest -- all kept -- to the front of the buffer.
 */
static void fastForward_keepFrom(FastForward *f, int64_t xay) {
    assert(xay >= f->base && xay < f->base + f->length);
    int64_t first = xay - f->base;
    assert(f->diagonals[first].where == FAST_STORED);
    int64_t cellStart = f->diagonals[first].offset;
    memmove(f->cells, f->cells + cellStart, sizeof(double) * (f->cellsLength - cellStart));
    f->cellsLength -= cellStart;
    memmove(f->diagonals, f->diagonals + first, sizeof(FastDiagonal) * (f->length - first));
    f->length -= first;
    f->base = xay;
    for (int64_t i = 0; i < f->length; i++) {
        assert(f->diagonals[i].where == FAST_STORED);
        f->diagonals[i].offset -= cellStart;
    }
    f->blockStart = -1;
}

/*
 * Where cell i of a diagonal finds its neighbour in another diagonal, whose xmy is the cell's
 * plus xmyDelta: the neighbour is cell i + shift, and exists for i in [lo, hi).
 */
static void fastDiagonal_neighbours(FastDiagonal *cur, FastDiagonal *other, int64_t xmyDelta,
                                    int64_t *shift, int64_t *lo, int64_t *hi) {
    if (other == NULL) {
        *shift = 0;
        *lo = 0;
        *hi = 0;
        return;
    }
    assert((cur->xmyL + xmyDelta - other->xmyL) % 2 == 0);
    *shift = (cur->xmyL + xmyDelta - other->xmyL) / 2;
    *lo = *shift < 0 ? -*shift : 0;
    *hi = other->width - *shift < cur->width ? other->width - *shift : cur->width;
    if (*hi < *lo) {
        *hi = *lo;
    }
}

/*
 * The values of a diagonal are kept relative to a scale, exp(logScale), chosen so they neither
 * underflow nor overflow.  Rescaling every diagonal cost a pass over it for nothing: probabilities
 * shrink by a bounded factor per diagonal, so the scale only needs moving every few dozen.  Moves
 * it, if the diagonal's largest value has drifted out of range, and returns the log of the factor
 * the cells were divided by.
 */
#define FAST_SCALE_LOW 1e-100
#define FAST_SCALE_HIGH 1e100

static double fastDiagonal_rescale(double *restrict cells, int64_t n, double max) {
    if (max <= 0.0 || (max >= FAST_SCALE_LOW && max <= FAST_SCALE_HIGH)) {
        return 0.0;
    }
    double r = 1.0 / max;
    for (int64_t i = 0; i < n; i++) {
        cells[i] *= r;
    }
    return log(max);
}

/*
 * out[i] = em[i] * sum_k t[k] * in[k][i] for i in [0, n), or without the em factor if em is NULL,
 * returning the largest out[i].  One loop per number of terms, so each case vectorises and the
 * transitions a state lacks cost nothing.  The terms are summed in order, as the separate passes
 * this replaced did.
 */
#define FAST_SUM_LOOP(EXPR) \
    _Pragma("omp simd reduction(max:max)") \
    for (int64_t i = 0; i < n; i++) { \
        double v = (EXPR); \
        out[i] = v; \
        max = v > max ? v : max; \
    }

static double fastSum(double *restrict out, const double *restrict em, const double *const *in, const double *t,
                      int64_t terms, int64_t n) {
    double max = 0.0;
    const double *restrict i0 = terms > 0 ? in[0] : NULL, *restrict i1 = terms > 1 ? in[1] : NULL;
    const double *restrict i2 = terms > 2 ? in[2] : NULL, *restrict i3 = terms > 3 ? in[3] : NULL;
    const double *restrict i4 = terms > 4 ? in[4] : NULL;
    double t0 = terms > 0 ? t[0] : 0, t1 = terms > 1 ? t[1] : 0, t2 = terms > 2 ? t[2] : 0;
    double t3 = terms > 3 ? t[3] : 0, t4 = terms > 4 ? t[4] : 0;
    if (em != NULL) {
        switch (terms) {
            case 0: FAST_SUM_LOOP(0.0); break;
            case 1: FAST_SUM_LOOP(em[i] * (t0 * i0[i])); break;
            case 2: FAST_SUM_LOOP(em[i] * (t0 * i0[i] + t1 * i1[i])); break;
            case 3: FAST_SUM_LOOP(em[i] * (t0 * i0[i] + t1 * i1[i] + t2 * i2[i])); break;
            case 4: FAST_SUM_LOOP(em[i] * (t0 * i0[i] + t1 * i1[i] + t2 * i2[i] + t3 * i3[i])); break;
            default: FAST_SUM_LOOP(em[i] * (t0 * i0[i] + t1 * i1[i] + t2 * i2[i] + t3 * i3[i] + t4 * i4[i]));
        }
    } else {
        switch (terms) {
            case 0: FAST_SUM_LOOP(0.0); break;
            case 1: FAST_SUM_LOOP(t0 * i0[i]); break;
            case 2: FAST_SUM_LOOP(t0 * i0[i] + t1 * i1[i]); break;
            case 3: FAST_SUM_LOOP(t0 * i0[i] + t1 * i1[i] + t2 * i2[i]); break;
            case 4: FAST_SUM_LOOP(t0 * i0[i] + t1 * i1[i] + t2 * i2[i] + t3 * i3[i]); break;
            default: FAST_SUM_LOOP(t0 * i0[i] + t1 * i1[i] + t2 * i2[i] + t3 * i3[i] + t4 * i4[i]);
        }
    }
    return max;
}

static void fastZero(double *restrict out, int64_t n) {
    for (int64_t i = 0; i < n; i++) {
        out[i] = 0.0;
    }
}

/*
 * The two sequences as bytes, with Y reversed, so that walking along a diagonal -- x up, y
 * down -- reads both forwards, and the emission lookups below vectorise.
 */
typedef struct _fastSymbols {
    uint8_t *x, *yReversed;
    int64_t lX, lY;
} FastSymbols;

static FastSymbols fastSymbols_construct(const SymbolString sX, const SymbolString sY) {
    FastSymbols f;
    f.lX = sX.length;
    f.lY = sY.length;
    f.x = st_malloc(sX.length + 1);
    f.yReversed = st_malloc(sY.length + 1);
    for (int64_t i = 0; i < sX.length; i++) {
        f.x[i] = (uint8_t) sX.sequence[i];
    }
    for (int64_t i = 0; i < sY.length; i++) {
        f.yReversed[i] = (uint8_t) sY.sequence[sY.length - 1 - i];
    }
    return f;
}

static void fastSymbols_destruct(FastSymbols f) {
    free(f.x);
    free(f.yReversed);
}

/*
 * The emission of cells [lo, hi) of a diagonal into a state of the given type, into em[lo, hi).
 * x0/y0 are the coordinates of cell 0, and cell i is at (x0 + i, y0 - i).  For the forward
 * recursion a cell emits its own symbols (the ones before it, at x-1 and y-1); for the backward
 * it is the successor's (at x and y) -- pass x0 and y0 one larger in that case.
 */
static void fastDiagonal_emissions(PairHmm *hmm, int64_t type, FastSymbols *symbols,
                                   int64_t x0, int64_t y0, int64_t lo, int64_t hi, double scale, double *restrict em) {
    double *restrict e = em + lo;
    int64_t n = hi - lo, i = 0;
    // Only the sequences the type emits are read: at the edge of the matrix the other one's index
    // is off its end.  cX[i] is X[x0 - 1 + lo + i] and cY[i] is Y[y0 - 1 - lo - i].
    const uint8_t *restrict cX = type != PAIR_HMM_GAP_Y ? symbols->x + (x0 - 1 + lo) : NULL;
    const uint8_t *restrict cY = type != PAIR_HMM_GAP_X ? symbols->yReversed + (symbols->lY - y0 + lo) : NULL;
#if defined(__AVX2__)
    // Table lookups the compiler will not vectorise on its own
    __m256d vScale = _mm256_set1_pd(scale);
    __m128i vSymbols = _mm_set1_epi32(SYMBOL_NUMBER);
    for (; i + 4 <= n; i += 4) {
        int32_t bytes;
        __m256d v;
        if (type == PAIR_HMM_GAP_X) {
            memcpy(&bytes, cX + i, sizeof(int32_t));
            v = _mm256_i32gather_pd(hmm->eGapX, _mm_cvtepu8_epi32(_mm_cvtsi32_si128(bytes)), 8);
        } else if (type == PAIR_HMM_GAP_Y) {
            memcpy(&bytes, cY + i, sizeof(int32_t));
            v = _mm256_i32gather_pd(hmm->eGapY, _mm_cvtepu8_epi32(_mm_cvtsi32_si128(bytes)), 8);
        } else {
            memcpy(&bytes, cX + i, sizeof(int32_t));
            __m128i vx = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(bytes));
            memcpy(&bytes, cY + i, sizeof(int32_t));
            __m128i vy = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(bytes));
            v = _mm256_mul_pd(_mm256_i32gather_pd(hmm->eMatch, _mm_add_epi32(_mm_mullo_epi32(vx, vSymbols), vy), 8), vScale);
        }
        _mm256_storeu_pd(e + i, v);
    }
#endif
    if (type == PAIR_HMM_GAP_X) {
        for (; i < n; i++) {
            e[i] = hmm->eGapX[cX[i]];
        }
    } else if (type == PAIR_HMM_GAP_Y) {
        for (; i < n; i++) {
            e[i] = hmm->eGapY[cY[i]];
        }
    } else {
        for (; i < n; i++) {
            e[i] = hmm->eMatch[cX[i] * SYMBOL_NUMBER + cY[i]] * scale;
        }
    }
}

/*
 * Forward values for diagonal cur, from the diagonals one (m1) and two (m2) before it.
 */
static void fastForwardDiagonal(PairHmm *hmm, int64_t xay, FastDiagonal *cur, double *restrict curCells,
                                FastDiagonal *m1, const double *m1Cells, FastDiagonal *m2, const double *m2Cells,
                                FastSymbols *symbols, double *restrict em) {
    int64_t w = cur->width;
    int64_t x0 = (xay + cur->xmyL) / 2, y0 = (xay - cur->xmyL) / 2;
    // Matches read from two diagonals back, whose scale can differ from the previous diagonal's,
    // which this one's is taken relative to.  Fold the difference into the match emissions.
    double matchScale = m2 != NULL ? (m2->logScale == m1->logScale ? 1.0 : exp(m2->logScale - m1->logScale)) : 0.0;
    double max = 0.0;
    for (int64_t type = 0; type < 3; type++) {
        FastDiagonal *src = type == PAIR_HMM_MATCH ? m2 : m1;
        const double *srcCells = type == PAIR_HMM_MATCH ? m2Cells : m1Cells;
        int64_t shift, lo, hi;
        fastDiagonal_neighbours(cur, src, type == PAIR_HMM_GAP_X ? -1 : (type == PAIR_HMM_GAP_Y ? 1 : 0), &shift, &lo, &hi);
        if (hi > lo) {
            fastDiagonal_emissions(hmm, type, symbols, x0, y0, lo, hi, matchScale, em);
        }
        for (int64_t s = 0; s < hmm->stateNumber; s++) {
            if (hmm->type[s] != type) {
                continue;
            }
            double *restrict out = curCells + s * w;
            fastZero(out, lo);
            fastZero(out + hi, w - hi);
            if (hi > lo) {
                const double *in[PAIR_HMM_MAX_STATES];
                for (int64_t k = 0; k < hmm->into[s].number; k++) {
                    in[k] = srcCells + hmm->into[s].state[k] * src->width + shift + lo;
                }
                double m = fastSum(out + lo, em + lo, in, hmm->into[s].t, hmm->into[s].number, hi - lo);
                max = m > max ? m : max;
            }
        }
    }
    cur->logScale = m1->logScale + fastDiagonal_rescale(curCells, w * hmm->stateNumber, max);
}

/*
 * Backward values for diagonal cur, from the diagonals one (p1) and two (p2) after it.  w is
 * scratch for stateNumber * width values.
 */
static void fastBackwardDiagonal(PairHmm *hmm, int64_t xay, FastDiagonal *cur, double *restrict curCells,
                                 FastDiagonal *p1, const double *p1Cells, FastDiagonal *p2, const double *p2Cells,
                                 FastSymbols *symbols, double *restrict w, double *restrict em) {
    int64_t width = cur->width;
    int64_t x0 = (xay + cur->xmyL) / 2, y0 = (xay - cur->xmyL) / 2;
    double matchScale = p2 != NULL ? (p2->logScale == p1->logScale ? 1.0 : exp(p2->logScale - p1->logScale)) : 0.0;
    // w[s] = what entering state s from this cell is worth: the successor's emission and backward value
    for (int64_t type = 0; type < 3; type++) {
        FastDiagonal *succ = type == PAIR_HMM_MATCH ? p2 : p1;
        const double *succCells = type == PAIR_HMM_MATCH ? p2Cells : p1Cells;
        int64_t shift, lo, hi;
        fastDiagonal_neighbours(cur, succ, type == PAIR_HMM_GAP_X ? 1 : (type == PAIR_HMM_GAP_Y ? -1 : 0), &shift, &lo, &hi);
        if (hi > lo) {
            fastDiagonal_emissions(hmm, type, symbols, x0 + 1, y0 + 1, lo, hi, matchScale, em);
        }
        for (int64_t s = 0; s < hmm->stateNumber; s++) {
            if (hmm->type[s] != type) {
                continue;
            }
            double *restrict ws = w + s * width;
            fastZero(ws, lo);
            fastZero(ws + hi, width - hi);
            if (hi > lo) {
                const double *restrict in = succCells + s * succ->width + shift + lo;
                double *restrict out = ws + lo;
                const double *restrict e = em + lo;
                for (int64_t i = 0; i < hi - lo; i++) {
                    out[i] = e[i] * in[i];
                }
            }
        }
    }
    double max = 0.0;
    for (int64_t f = 0; f < hmm->stateNumber; f++) {
        const double *in[PAIR_HMM_MAX_STATES];
        for (int64_t k = 0; k < hmm->outOf[f].number; k++) {
            in[k] = w + hmm->outOf[f].state[k] * width;
        }
        double m = fastSum(curCells + f * width, NULL, in, hmm->outOf[f].t, hmm->outOf[f].number, width);
        max = m > max ? m : max;
    }
    cur->logScale = p1->logScale + fastDiagonal_rescale(curCells, width * hmm->stateNumber, max);
}

/*
 * Posteriors of the states the caller wants lists for, over one diagonal.
 */
static void fastPosteriors(StateMachine *sM, int64_t xay, FastDiagonal *fd, const double *fCells,
                           FastDiagonal *bd, const double *bCells, double logTotal, PairwiseAlignmentParameters *p,
                           stList *alignedPairs, stList *gapXPairs, stList *gapYPairs) {
    assert(fd->xmyL == bd->xmyL && fd->width == bd->width);
    int64_t w = fd->width;
    int64_t x0 = (xay + fd->xmyL) / 2, y0 = (xay - fd->xmyL) / 2;
    double logNorm = fd->logScale + bd->logScale - logTotal;
    bool direct = logNorm < 700.0; // else exp(logNorm) would overflow, and it is done cell by cell in log space
    double norm = direct ? exp(logNorm) : 0.0;
    // Most cells fall short of the threshold, so test them against it in their own scale -- one
    // multiply and one compare -- a hair low, so addPosteriorProb makes the exact decision.
    double fbThreshold = direct ? p->threshold * exp(-logNorm) * (1.0 - 1e-9) : 0.0;
    bool everyCell = p->threshold <= 0.0; // then the log space path reports every cell, zeros included
    int64_t states[3] = { sM->matchState, sM->gapXState, sM->gapYState };
    stList *lists[3] = { alignedPairs, gapXPairs, gapYPairs };
    for (int64_t k = 0; k < 3; k++) {
        if (lists[k] == NULL) {
            continue;
        }
        // Matches and gaps in y need x > 0, i.e. i >= 1 - x0; matches and gaps in x need y > 0, i.e. i < y0
        int64_t lo = k != 2 && 1 - x0 > 0 ? 1 - x0 : 0;
        int64_t hi = k != 1 && y0 < w ? y0 : w;
        const double *f = fCells + states[k] * w, *b = bCells + states[k] * w;
        for (int64_t i = lo; i < hi; i++) {
            double fb = f[i] * b[i];
            if (everyCell || (fb >= fbThreshold && fb > 0.0)) {
                addPosteriorProb(x0 + i, y0 - i, direct ? fb * norm : (fb > 0.0 ? exp(log(fb) + logNorm) : 0.0), lists[k], p);
            }
        }
    }
}

/*
 * The fast equivalent of getPosteriorProbsWithBanding with diagonalCalculationPosteriorMatchProbs
 * (gapXPairs and gapYPairs NULL) or diagonalCalculationPosteriorProbs.  Returns false, having done
 * nothing, if the state machine is not one it can use.
 */
/*
 * The machine's flattened form, made the first time it is asked for and kept on the machine:
 * flattening runs the cell calculation for every pair of symbols and takes hundreds of exps, which
 * cost more than the dp itself for the thousands of few-base alignments a bar run makes.  Threads
 * share machines, so a thread that loses the race to store its copy frees it and uses the winner's.
 */
static PairHmm *pairHmm_get(StateMachine *sM) {
    PairHmm *hmm = __atomic_load_n((PairHmm **) &sM->flat, __ATOMIC_ACQUIRE);
    if (hmm == NULL) {
        PairHmm *made = st_malloc(sizeof(PairHmm));
        pairHmm_construct(made, sM);
        PairHmm *expected = NULL;
        if (__atomic_compare_exchange_n((PairHmm **) &sM->flat, &expected, made, 0, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
            hmm = made;
        } else {
            free(made);
            hmm = expected;
        }
    }
    return hmm;
}

static bool getPosteriorProbsWithBandingFast(StateMachine *sM, stList *anchorPairs, const SymbolString sX, const SymbolString sY,
                                             PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd,
                                             stList *alignedPairs, stList *gapXPairs, stList *gapYPairs) {
    PairHmm hmm = *pairHmm_get(sM);
    if (!hmm.ok) {
        return 0;
    }
    int64_t diagonalNumber = sX.length + sY.length;
    if (diagonalNumber == 0) {
        return 1;
    }
    int64_t ns = hmm.stateNumber;

#if defined(__SSE2__)
    // Cells far off the alignment decay towards zero; flush them rather than crawl through denormals
    unsigned int oldCsr = _mm_getcsr();
    _mm_setcsr(oldCsr | 0x8040); // flush-to-zero and denormals-are-zero
#endif

    Band *band = p->dynamicAnchorExpansion ? band_constructDynamic(anchorPairs, sX.length, sY.length) :
                 band_construct(anchorPairs, sX.length, sY.length, p->diagonalExpansion);
    FastForward fwd;
    int64_t chunkNumber;
    FastChunk *chunks = fastForward_construct(&fwd, band, diagonalNumber, ns, p, &chunkNumber);
    FastSymbols symbols = fastSymbols_construct(sX, sY);

    // scratch, sized to the widest diagonal
    int64_t maxWidth = 1;
    for (int64_t d = 0; d <= diagonalNumber; d++) {
        int64_t w = diagonal_getWidth(band->diagonals[d]);
        maxWidth = w > maxWidth ? w : maxWidth;
    }
    double *em = st_malloc(sizeof(double) * maxWidth);
    double *scratch = st_malloc(sizeof(double) * maxWidth * ns);
    FastDiagonal back[3]; // the backward diagonals xay, xay + 1 and xay + 2, round robin
    double *backCells[3];
    for (int64_t i = 0; i < 3; i++) {
        backCells[i] = st_malloc(sizeof(double) * maxWidth * ns);
    }

    FastDiagonal *d0 = fastForward_add(&fwd, band->diagonals[0], ns, 1);
    for (int64_t s = 0; s < ns; s++) {
        double v = alignmentHasRaggedLeftEnd ? hmm.raggedStart[s] : hmm.start[s];
        for (int64_t i = 0; i < d0->width; i++) {
            fwd.cells[d0->offset + s * d0->width + i] = v;
        }
    }

    for (int64_t c = 0; c < chunkNumber; c++) {
        FastChunk *chunk = &chunks[c];
        for (int64_t xay = chunk->start + 1; xay <= chunk->end; xay++) {
            FastDiagonal *cur = fastForward_add(&fwd, band->diagonals[xay], ns, fastChunk_stores(chunk, xay));
            FastDiagonal *m1 = fastForward_get(&fwd, xay - 1), *m2 = fastForward_get(&fwd, xay - 2);
            fastForwardDiagonal(&hmm, xay, cur, fastForward_cells(&fwd, cur), m1, fastForward_cells(&fwd, m1),
                                m2, m2 != NULL ? fastForward_cells(&fwd, m2) : NULL, &symbols, em);
        }

        // Treat every cell of the chunk's last diagonal as an end point, and walk back to where the
        // last traceback stopped.  Posteriors are taken only up to traceBackDiagonals before the
        // end, where that pretence has worn off; the rest are redone by the next traceback.
        bool atEnd = chunk->end == diagonalNumber;
        const double *endProbs = atEnd && alignmentHasRaggedRightEnd ? hmm.raggedEnd : hmm.end;
        FastDiagonal *last = fastForward_get(&fwd, chunk->end);
        const double *lastCells = fastForward_cells(&fwd, last);
        double total = 0.0;
        for (int64_t s = 0; s < ns; s++) {
            for (int64_t i = 0; i < last->width; i++) {
                total += lastCells[s * last->width + i] * endProbs[s];
            }
        }
        double logTotal = last->logScale + log(total);
        bool havePosteriors = total > 0.0;

        for (int64_t xay2 = chunk->end; xay2 > chunk->tracedBackTo; xay2--) {
            FastDiagonal *b = &back[xay2 % 3];
            double *bCells = backCells[xay2 % 3];
            FastDiagonal *f = fastForward_get(&fwd, xay2);
            b->xmyL = f->xmyL;
            b->width = f->width;
            if (xay2 == chunk->end) {
                for (int64_t s = 0; s < ns; s++) {
                    for (int64_t i = 0; i < b->width; i++) {
                        bCells[s * b->width + i] = endProbs[s];
                    }
                }
                b->logScale = 0.0;
            } else {
                FastDiagonal *p1 = &back[(xay2 + 1) % 3], *p2 = xay2 + 2 <= chunk->end ? &back[(xay2 + 2) % 3] : NULL;
                fastBackwardDiagonal(&hmm, xay2, b, bCells, p1, backCells[(xay2 + 1) % 3],
                                     p2, p2 != NULL ? backCells[(xay2 + 2) % 3] : NULL, &symbols, scratch, em);
            }
            if (xay2 > chunk->tracedBackFrom || !havePosteriors) {
                continue;
            }
            if (f->where != FAST_STORED && !(f->where == FAST_BLOCK && xay2 >= fwd.blockStart)) {
                // Recompute the block this diagonal is in, from the kept pair before it.  The same
                // code on the same inputs, so the same values, scales and all, as the first time.
                int64_t blockStart = chunk->start + 1 + (xay2 - chunk->start - 1) / chunk->blockSize * chunk->blockSize;
                fwd.blockStart = blockStart;
                fwd.blockCellsLength = 0;
                for (int64_t xay3 = blockStart; xay3 <= xay2; xay3++) {
                    FastDiagonal *cur = fastForward_get(&fwd, xay3);
                    cur->where = FAST_BLOCK;
                    cur->offset = fwd.blockCellsLength;
                    fwd.blockCellsLength += cur->width * ns;
                    assert(fwd.blockCellsLength <= fwd.blockCellsCapacity);
                    FastDiagonal *m1 = fastForward_get(&fwd, xay3 - 1), *m2 = fastForward_get(&fwd, xay3 - 2);
                    fastForwardDiagonal(&hmm, xay3, cur, fastForward_cells(&fwd, cur), m1, fastForward_cells(&fwd, m1),
                                        m2, m2 != NULL ? fastForward_cells(&fwd, m2) : NULL, &symbols, em);
                }
            }
            fastPosteriors(sM, xay2, f, fastForward_cells(&fwd, f), b, bCells, logTotal, p, alignedPairs, gapXPairs, gapYPairs);
        }
        if (!atEnd) {
            fastForward_keepFrom(&fwd, chunk->tracedBackFrom);
        }
    }

    free(chunks);
    fastForward_destruct(&fwd);
    fastSymbols_destruct(symbols);
    free(em);
    free(scratch);
    for (int64_t i = 0; i < 3; i++) {
        free(backCells[i]);
    }
    band_destruct(band);
#if defined(__SSE2__)
    _mm_setcsr(oldCsr);
#endif
    return 1;
}

///////////////////////////////////
///////////////////////////////////
//Banded alignment routine to calculate posterior match probs
//
//
///////////////////////////////////
///////////////////////////////////

void getPosteriorProbsWithBanding(StateMachine *sM, stList *anchorPairs, const SymbolString sX, const SymbolString sY,
                                  PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd,
                                  void (*diagonalPosteriorProbFn)(StateMachine *, int64_t, DpMatrix *, DpMatrix *, const SymbolString, const SymbolString, double,
                                                                  PairwiseAlignmentParameters *, void *), void *extraArgs) {
    //Prerequisites
    assert(p->traceBackDiagonals >= 1);
    assert(p->diagonalExpansion >= 0);
    assert(p->diagonalExpansion % 2 == 0);
    assert(p->minDiagsBetweenTraceBack >= 2);
    assert(p->traceBackDiagonals + 1 < p->minDiagsBetweenTraceBack);

    int64_t diagonalNumber = sX.length + sY.length;
    if (diagonalNumber == 0) { //Deal with trivial case
        return;
    }

    //The posterior probabilities have a fast path.  The expectations for EM training do not.
    if (diagonalPosteriorProbFn == diagonalCalculationPosteriorMatchProbs &&
        getPosteriorProbsWithBandingFast(sM, anchorPairs, sX, sY, p, alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd,
                                         ((void **) extraArgs)[0], NULL, NULL)) {
        return;
    }
    if (diagonalPosteriorProbFn == diagonalCalculationPosteriorProbs &&
        getPosteriorProbsWithBandingFast(sM, anchorPairs, sX, sY, p, alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd,
                                         ((void **) extraArgs)[0], ((void **) extraArgs)[2], ((void **) extraArgs)[4])) {
        return;
    }

    //Primitives for the forward matrix recursion
    Band *band = p->dynamicAnchorExpansion ? band_constructDynamic(anchorPairs, sX.length, sY.length) : band_construct(anchorPairs, sX.length, sY.length, p->diagonalExpansion);
    BandIterator *forwardBandIterator = bandIterator_construct(band);
    DpMatrix *forwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
    dpDiagonal_initialiseValues(dpMatrix_createDiagonal(forwardDpMatrix, bandIterator_getNext(forwardBandIterator)), sM,
                                alignmentHasRaggedLeftEnd ? sM->raggedStartStateProb : sM->startStateProb); //Initialise forward matrix.

    //Backward matrix.
    DpMatrix *backwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);

    int64_t tracedBackTo = 0;
    int64_t totalPosteriorCalculations = 0;
    while (1) { //Loop that moves through the matrix forward
        Diagonal diagonal = bandIterator_getNext(forwardBandIterator);

        //Forward calculation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(forwardDpMatrix, diagonal));
        diagonalCalculationForward(sM, diagonal_getXay(diagonal), forwardDpMatrix, sX, sY);

        bool atEnd = diagonal_getXay(diagonal) == diagonalNumber; //Condition true at the end of the matrix
        bool tracebackPoint = diagonal_getXay(diagonal) >= tracedBackTo + p->minDiagsBetweenTraceBack
                              && diagonal_getWidth(diagonal) <= p->diagonalExpansion * 2 + 1; //Condition true when we want to do an intermediate traceback.

        //Traceback
        if (atEnd || tracebackPoint) {
            //Initialise the last row (until now) of the backward matrix to represent an end point
            dpDiagonal_initialiseValues(dpMatrix_createDiagonal(backwardDpMatrix, diagonal), sM,
                                        (atEnd && alignmentHasRaggedRightEnd) ? sM->raggedEndStateProb : sM->endStateProb);
            if (diagonal_getXay(diagonal) > tracedBackTo + 1) { //This is a diagonal between the place we trace back to and where we trace back from
                DpDiagonal *j = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal) - 1);
                assert(j != NULL);
                dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, j->diagonal));
            }

            //Do walk back
            BandIterator *backwardBandIterator = bandIterator_clone(forwardBandIterator);
            Diagonal diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            assert(diagonal_getXay(diagonal2) == diagonal_getXay(diagonal));
            int64_t tracedBackFrom = diagonal_getXay(diagonal) - (atEnd ? 0 : p->traceBackDiagonals + 1);
            double totalProbability = LOG_ZERO;
            int64_t totalPosteriorCalculationsThisTraceback = 0;
            while (diagonal_getXay(diagonal2) > tracedBackTo) {
                //Create the earlier diagonal
                if (diagonal_getXay(diagonal2) > tracedBackTo + 2) {
                    DpDiagonal *j = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2) - 2);
                    assert(j != NULL);
                    dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, j->diagonal));
                }
                if (diagonal_getXay(diagonal2) > tracedBackTo + 1) {
                    diagonalCalculationBackward(sM, diagonal_getXay(diagonal2), backwardDpMatrix, sX, sY);
                }
                if (diagonal_getXay(diagonal2) <= tracedBackFrom) {
                    assert(dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)) != NULL);
                    assert(dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)-1) != NULL);
                    assert(dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)) != NULL);
                    if (diagonal_getXay(diagonal2) != diagonalNumber) {
                        assert(dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)+1) != NULL);
                    }
                    if (totalPosteriorCalculationsThisTraceback++ % 10 == 0) {
                        double newTotalProbability = diagonalCalculationTotalProbability(sM, diagonal_getXay(diagonal2),
                                                                                         forwardDpMatrix, backwardDpMatrix, sX, sY);
                        if (totalPosteriorCalculationsThisTraceback != 1) {
                            assert(totalProbability + 1.0 > newTotalProbability);
                            assert(newTotalProbability + 1.0 > newTotalProbability);
                        }
                        totalProbability = newTotalProbability;
                    }

                    diagonalPosteriorProbFn(sM, diagonal_getXay(diagonal2), forwardDpMatrix, backwardDpMatrix, sX, sY,
                                            totalProbability, p, extraArgs);

                    if (diagonal_getXay(diagonal2) < tracedBackFrom || atEnd) {
                        dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)); //Delete forward diagonal after last access in posterior calculation
                    }
                }
                if (diagonal_getXay(diagonal2) + 1 <= diagonalNumber) {
                    dpMatrix_deleteDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2) + 1); //Delete backward diagonal after last access in backward calculation
                }
                diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            }
            tracedBackTo = tracedBackFrom;
            bandIterator_destruct(backwardBandIterator);
            dpMatrix_deleteDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2) + 1);
            dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2));
            //Check memory state.
            assert(dpMatrix_getActiveDiagonalNumber(backwardDpMatrix) == 0);
            totalPosteriorCalculations += totalPosteriorCalculationsThisTraceback;
            if (!atEnd) {
                assert(dpMatrix_getActiveDiagonalNumber(forwardDpMatrix) == p->traceBackDiagonals + 2);
            }
        }

        if (atEnd) {
            break;
        }
    }
    assert(totalPosteriorCalculations == diagonalNumber);
    assert(tracedBackTo == diagonalNumber);
    assert(dpMatrix_getActiveDiagonalNumber(backwardDpMatrix) == 0);
    assert(dpMatrix_getActiveDiagonalNumber(forwardDpMatrix) == 0);
    //Cleanup
    dpMatrix_destruct(forwardDpMatrix);
    dpMatrix_destruct(backwardDpMatrix);
    bandIterator_destruct(forwardBandIterator);
    band_destruct(band);
}

double getForwardProbWithBanding(StateMachine *sM, stList *anchorPairs, const SymbolString sX, const SymbolString sY,
                                 PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    //Prerequisites
    assert(p->traceBackDiagonals >= 1);
    assert(p->diagonalExpansion >= 0);
    assert(p->diagonalExpansion % 2 == 0);
    assert(p->minDiagsBetweenTraceBack >= 2);
    assert(p->traceBackDiagonals + 1 < p->minDiagsBetweenTraceBack);

    int64_t diagonalNumber = sX.length + sY.length;
    if (diagonalNumber == 0) { //Deal with trivial case
        return LOG_ONE;
    }

    //Primitives for the forward matrix recursion
    Band *band = band_construct(anchorPairs, sX.length, sY.length, p->diagonalExpansion);
    BandIterator *forwardBandIterator = bandIterator_construct(band);
    DpMatrix *forwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
    dpDiagonal_initialiseValues(dpMatrix_createDiagonal(forwardDpMatrix, bandIterator_getNext(forwardBandIterator)), sM,
                                alignmentHasRaggedLeftEnd ? sM->raggedStartStateProb : sM->startStateProb); //Initialise forward matrix.

    double totalLogProbability = LOG_ZERO;

    while (1) { //Loop that moves through the matrix forward
        Diagonal diagonal = bandIterator_getNext(forwardBandIterator);

        //Forward calculation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(forwardDpMatrix, diagonal));
        diagonalCalculationForward(sM, diagonal_getXay(diagonal), forwardDpMatrix, sX, sY);

        bool atEnd = diagonal_getXay(diagonal) == diagonalNumber; //Condition true at the end of the matrix
        if (atEnd) {
            //Backward matrix.
            DpMatrix *backwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
            dpDiagonal_initialiseValues(dpMatrix_createDiagonal(backwardDpMatrix, diagonal), sM,
                                        alignmentHasRaggedRightEnd ? sM->raggedEndStateProb : sM->endStateProb);
            totalLogProbability = diagonalCalculationTotalProbability(sM, diagonalNumber,
                                                                      forwardDpMatrix, backwardDpMatrix, sX, sY);
            dpMatrix_deleteDiagonal(backwardDpMatrix, diagonalNumber);
            dpMatrix_destruct(backwardDpMatrix);
            break;
        }
    }
    //Cleanup
    for (int64_t i=0; i<=diagonalNumber; i++) {
        dpMatrix_deleteDiagonal(forwardDpMatrix, i);
    }
    dpMatrix_destruct(forwardDpMatrix);
    bandIterator_destruct(forwardBandIterator);
    band_destruct(band);

    return totalLogProbability;
}

/*
 * Computes for the forward log probability of aligning the two sequences
 */
double computeForwardProbability(char *seqX, char *seqY, stList *anchorPairs, PairwiseAlignmentParameters *p, StateMachine *sM,
                                 bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {

    SymbolString sX = symbolString_construct(seqX, strlen(seqX));
    SymbolString sY = symbolString_construct(seqY, strlen(seqY));

    double totalLogProb = getForwardProbWithBanding(sM, anchorPairs, sX, sY,
                                                    p, alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd);

    symbolString_destruct(sX);
    symbolString_destruct(sY);

    return totalLogProb;
}

///////////////////////////////////
///////////////////////////////////
//Blast anchoring functions
//
//Use lastz to get sets of anchors
///////////////////////////////////
///////////////////////////////////

static char *makeUpperCase(const char *s, int64_t l) {
    char *s2 = stString_copy(s);
    for (int64_t i = 0; i < l; i++) {
        s2[i] = toupper(s[i]);
    }
    return s2;
}

static void writeSequenceToFile(char *file, const char *name, const char *sequence) {
    FILE *fileHandle = fopen(file, "w");
    fastaWrite((char *) sequence, (char *) name, fileHandle);
    fclose(fileHandle);
}

stList *convertPairwiseForwardStrandAlignmentToAnchorPairs(struct PairwiseAlignment *pA, int64_t trim, int64_t expansion) {
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct); //the list to put the output in
    int64_t j = pA->start1;
    int64_t k = pA->start2;
    assert(pA->strand1);
    assert(pA->strand2);
    for (int64_t i = 0; i < pA->operationList->length; i++) {
        struct AlignmentOperation *op = pA->operationList->list[i];
        if (op->opType == PAIRWISE_MATCH) {
            for (int64_t l = trim; l < op->length - trim; l++) {
                stList_append(alignedPairs, stIntTuple_construct3(j + l, k + l, expansion));
            }
        }
        if (op->opType != PAIRWISE_INDEL_Y) {
            j += op->length;
        }
        if (op->opType != PAIRWISE_INDEL_X) {
            k += op->length;
        }
    }

    assert(j == pA->end1);
    assert(k == pA->end2);
    return alignedPairs;
}

stList *getBlastPairs(const char *sX, const char *sY, int64_t lX, int64_t lY, PairwiseAlignmentParameters *p, bool repeatMask) {
    /*
     * Uses lastz to compute a bunch of monotonically increasing pairs such that for any pair of consecutive pairs in the list
     * (x1, y1) (x2, y2) in the set of aligned pairs x1 appears before x2 in X and y1 appears before y2 in Y.
     */
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct); //the list to put the output in

    if (lX == 0 || lY == 0) {
        return alignedPairs;
    }

    if (!repeatMask) { // Optionally remove repeat masking
        sX = makeUpperCase(sX, lX);
        sY = makeUpperCase(sY, lY);
    }

    // Get temporary files - notably these functions are hopefully thread safe now
    char *tempFile1 = getTempFile();
    char *tempFile2 = lY > 1000 ? getTempFile() : NULL;

    // Write the sequences to be aligned to the temporary files and construct the lastz command using the temporary files
    writeSequenceToFile(tempFile1, "a", sX);
    char *command;
    if (lY > 1000) {
        writeSequenceToFile(tempFile2, "b", sY);
        command = stString_print("cPecanLastz --hspthresh=800 --chain --strand=plus --gapped --format=cigar --ambiguous=iupac,100,100 %s %s",
                                 tempFile1, tempFile2);
    } else {
        command = stString_print("echo '>b\n%s\n' | cPecanLastz --hspthresh=800 --chain --strand=plus --gapped --format=cigar --ambiguous=iupac,100,100 %s",
                                 sY, tempFile1);
    }

#if defined(_OPENMP) && defined(PECAN_LOCK_POPEN)
    omp_set_lock(&(p->lastzLock));
#endif
    // Call lastz with popen
    FILE *fileHandle = popen(command, "r");
#if defined(_OPENMP) && defined(PECAN_LOCK_POPEN)
    omp_unset_lock(&(p->lastzLock));
#endif

    // Check popen worked
    if (fileHandle == NULL) {
        st_errnoAbort("Problems with lastz pipe");
    }

    //Read from alignmemnts from the stream
    struct PairwiseAlignment *pA;
    while ((pA = cigarRead(fileHandle)) != NULL) {
        assert(strcmp(pA->contig1, "a") == 0);
        assert(strcmp(pA->contig2, "b") == 0);
        stList *alignedPairsForCigar = convertPairwiseForwardStrandAlignmentToAnchorPairs(pA, p->constraintDiagonalTrim, p->diagonalExpansion);
        stList_appendAll(alignedPairs, alignedPairsForCigar);
        stList_setDestructor(alignedPairsForCigar, NULL);
        stList_destruct(alignedPairsForCigar);
        destructPairwiseAlignment(pA);
    }

#if defined(_OPENMP) && defined(PECAN_LOCK_POPEN)
    omp_set_lock(&(p->lastzLock));
#endif
    // Close the stream
    int64_t status = pclose(fileHandle);
#if defined(_OPENMP) && defined(PECAN_LOCK_POPEN)
    omp_unset_lock(&(p->lastzLock));
#endif

    // Check we closed the process properly
    if (status != 0) {
        st_errnoAbort("pclose failed when getting rid of lastz pipe with value %" PRIi64 " and command %s", status,
                command);
    }
    free(command);

    //Remove temporary files
    status = remove(tempFile1);
    if (status != 0) {
        st_errnoAbort("Failed to remove first temporary file when running lastz");
    }
    free(tempFile1);
    if (tempFile2 != NULL) {
        status = remove(tempFile2);
        if (status != 0) {
            st_errnoAbort("Failed to remove second temporary file when running lastz");
        }
    }
    free(tempFile2);

    //Ensure the coordinates are increasing
    //stList_sort(alignedPairs, sortByXPlusYCoordinate);

    // If we removed the repeat masking, clean up the temporary sequences
    if (!repeatMask) {
        free((char *) sX);
        free((char *) sY);
    }

    // Convert to an alignment
    stList_sort(alignedPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
    stList *filteredAnchorPairs = filterToRemoveOverlap(alignedPairs);
    stList_destruct(alignedPairs);
    return filteredAnchorPairs;
}

static void convertBlastPairs(stList *alignedPairs2, int64_t offsetX, int64_t offsetY) {
    /*
     * Convert the coordinates of the computed pairs.
     */
    for (int64_t k = 0; k < stList_length(alignedPairs2); k++) {
        stIntTuple *i = stList_get(alignedPairs2, k);
        assert(stIntTuple_length(i) == 3);
        stList_set(alignedPairs2, k,
                   stIntTuple_construct3(stIntTuple_get(i, 0) + offsetX, stIntTuple_get(i, 1) + offsetY, stIntTuple_get(i, 2)));
        stIntTuple_destruct(i);
    }
}

stList *filterToRemoveOverlap(stList *sortedOverlappingPairs) {
    stList *nonOverlappingPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    //Traverse backwards
    stSortedSet *set = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn, NULL);
    int64_t pX = INT64_MAX, pY = INT64_MAX;
    for (int64_t i = stList_length(sortedOverlappingPairs) - 1; i >= 0; i--) {
        stIntTuple *pair = stList_get(sortedOverlappingPairs, i);
        int64_t x = stIntTuple_get(pair, 0);
        int64_t y = stIntTuple_get(pair, 1);
        if (x < pX && y < pY) {
            stSortedSet_insert(set, pair);
        }
        pX = x < pX ? x : pX;
        pY = y < pY ? y : pY;
    }

    //Traverse forwards to final set of pairs
    pX = INT64_MIN;
    pY = INT64_MIN;
#ifndef NDEBUG
    int64_t pY2 = INT64_MIN;
#endif
    for (int64_t i = 0; i < stList_length(sortedOverlappingPairs); i++) {
        stIntTuple *pair = stList_get(sortedOverlappingPairs, i);
        int64_t x = stIntTuple_get(pair, 0);
        int64_t y = stIntTuple_get(pair, 1);
        if (x > pX && y > pY && stSortedSet_search(set, pair) != NULL) {
            stList_append(nonOverlappingPairs, stIntTuple_construct3(x, y, stIntTuple_get(pair, 2)));
        }
#ifndef NDEBUG
        //Check things are sorted in the input
        assert(x >= pX);
        if (x == pX) {
            assert(y >= pY2);
        }
        pY2 = y;
#endif
        pX = x > pX ? x : pX;
        pY = y > pY ? y : pY;
    }
    stSortedSet_destruct(set);

    return nonOverlappingPairs;
}

static void getBlastPairsForPairwiseAlignmentParametersP(const char *sX, const char *sY, int64_t pX, int64_t pY,
                                                         int64_t x, int64_t y, PairwiseAlignmentParameters *p, stList *combinedAnchorPairs) {
    int64_t lX2 = x - pX;
    assert(lX2 >= 0);
    int64_t lY2 = y - pY;
    assert(lY2 >= 0);
    int64_t matrixSize = (int64_t) lX2 * lY2;
    if (matrixSize > p->anchorMatrixBiggerThanThis) {
        char *sX2 = stString_getSubString(sX, pX, lX2);
        char *sY2 = stString_getSubString(sY, pY, lY2);
        stList *bottomLevelAnchorPairs = getBlastPairs(sX2, sY2, lX2, lY2, p, matrixSize > p->repeatMaskMatrixBiggerThanThis);
        st_logDebug("Got %" PRIi64 " bottom level anchor pairs\n", stList_length(bottomLevelAnchorPairs));
        convertBlastPairs(bottomLevelAnchorPairs, pX, pY);
        free(sX2);
        free(sY2);
        stList_appendAll(combinedAnchorPairs, bottomLevelAnchorPairs);
        stList_setDestructor(bottomLevelAnchorPairs, NULL);
        stList_destruct(bottomLevelAnchorPairs);
    }
}

stList *getBlastPairsForPairwiseAlignmentParameters(const char *sX, const char *sY, const int64_t lX, const int64_t lY,
                                                    PairwiseAlignmentParameters *p) {
    //Anchor pairs
    stList *topLevelAnchorPairs = getBlastPairs(sX, sY, lX, lY, p, 1);
    st_logDebug("Got %" PRIi64 " top level anchor pairs \n", stList_length(topLevelAnchorPairs));

    int64_t pX = 0;
    int64_t pY = 0;
    stList *combinedAnchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(topLevelAnchorPairs); i++) {
        stIntTuple *anchorPair = stList_get(topLevelAnchorPairs, i);
        int64_t x = stIntTuple_get(anchorPair, 0);
        int64_t y = stIntTuple_get(anchorPair, 1);
        assert(x >= 0 && x < lX);
        assert(y >= 0 && y < lY);
        assert(x >= pX);
        assert(y >= pY);
        getBlastPairsForPairwiseAlignmentParametersP(sX, sY, pX, pY, x, y, p, combinedAnchorPairs);
        stList_append(combinedAnchorPairs, anchorPair);
        pX = x + 1;
        pY = y + 1;
    }
    getBlastPairsForPairwiseAlignmentParametersP(sX, sY, pX, pY, lX, lY, p, combinedAnchorPairs);
    stList_setDestructor(topLevelAnchorPairs, NULL);
    stList_destruct(topLevelAnchorPairs);
    st_logDebug("Got %" PRIi64 " combined anchor pairs\n", stList_length(combinedAnchorPairs));
    return combinedAnchorPairs;
}

stList *getAnchorPairsForPairwiseAlignmentParameters(const char *sX, const char *sY, const int64_t lX, const int64_t lY,
                                                    PairwiseAlignmentParameters *p) {
    if ((int64_t) lX * lY <= p->anchorMatrixBiggerThanThis) {
        return stList_construct();
    }
    if (p->anchorMethod == PAIRWISE_ANCHOR_SEED) {
        return getSeedAnchors(sX, sY, lX, lY, p);
    }
    if(p->useMumAnchors) {
        return getAlignedMums(sX, sY, lX, lY, p, 0, 0);
    }
    return getBlastPairsForPairwiseAlignmentParameters(sX, sY, lX, lY, p);
}

///////////////////////////////////
///////////////////////////////////
//Split large gap functions
//
//Functions to split up alignment around gaps in the anchors that are too large.
///////////////////////////////////
///////////////////////////////////

static bool getSplitPointsP(int64_t *x1, int64_t *y1, int64_t x2, int64_t y2, int64_t x3, int64_t y3,
                            stList *splitPoints, int64_t splitMatrixBiggerThanThis, bool skipBlock) {
    /*
     * x2/y2 are the previous anchor point, x3/y3 are the next anchor point. Gaps greater than (x3-x2)*(y3-y2) are split up.
     */
    int64_t lX2 = x3 - x2;
    int64_t lY2 = y3 - y2;
    int64_t matrixSize = lX2 * lY2;
    if (matrixSize > splitMatrixBiggerThanThis) {
        st_logDebug("Split point found at x1: %" PRIi64 " x2: %" PRIi64 " y1: %" PRIi64 " y2: %" PRIi64 "\n", x2, x3,
                y2, y3);
        int64_t maxSequenceLength = sqrt(splitMatrixBiggerThanThis);
        int64_t hX = lX2 / 2 > maxSequenceLength ? maxSequenceLength : lX2 / 2;
        int64_t hY = lY2 / 2 > maxSequenceLength ? maxSequenceLength : lY2 / 2;
        if(!skipBlock) {
            stList_append(splitPoints, stIntTuple_construct4(*x1, *y1, x2 + hX, y2 + hY));
        }
        *x1 = x3 - hX;
        *y1 = y3 - hY;
        return 1;
    }
    return 0;
}

stList *getSplitPoints(stList *anchorPairs, int64_t lX, int64_t lY, int64_t splitMatrixBiggerThanThis,
                       bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    int64_t x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    assert(lX >= 0);
    assert(lY >= 0);
    stList *splitPoints = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *anchorPair = stList_get(anchorPairs, i);
        int64_t x3 = stIntTuple_get(anchorPair, 0), y3 = stIntTuple_get(anchorPair, 1);
        getSplitPointsP(&x1, &y1, x2, y2, x3, y3, splitPoints, splitMatrixBiggerThanThis, alignmentHasRaggedLeftEnd && i == 0);
        assert(x3 >= x2);
        assert(y3 >= y2);
        assert(x3 < lX);
        assert(y3 < lY);
        x2 = x3 + 1;
        y2 = y3 + 1;
    }
    if(!getSplitPointsP(&x1, &y1, x2, y2, lX, lY, splitPoints, splitMatrixBiggerThanThis,
                        alignmentHasRaggedLeftEnd && stList_length(anchorPairs) == 0) || !alignmentHasRaggedRightEnd) {
        stList_append(splitPoints, stIntTuple_construct4(x1, y1, lX, lY));
    }

    if (stList_length(splitPoints) > 1) {
        st_logDebug("For sequences of length %" PRIi64 " and %" PRIi64 " we got %" PRIi64 " splits\n", lX, lY,
                stList_length(splitPoints));
    }
    return splitPoints;
}

static void convertAlignedPairs(stList *alignedPairs2, int64_t offsetX, int64_t offsetY) {
    /*
     * Convert the coordinates of the computed pairs.
     */
    for (int64_t k = 0; k < stList_length(alignedPairs2); k++) {
        stIntTuple *i = stList_get(alignedPairs2, k);
        assert(stIntTuple_length(i) == 3);
        stList_set(alignedPairs2, k,
                   stIntTuple_construct3(stIntTuple_get(i, 0), stIntTuple_get(i, 1) + offsetX,
                                         stIntTuple_get(i, 2) + offsetY));
        stIntTuple_destruct(i);
    }
}

void getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(StateMachine *sM, stList *anchorPairs, const char *sX, const char *sY,
                                                                int64_t lX, int64_t lY, PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd,
                                                                bool alignmentHasRaggedRightEnd,
                                                                void (*diagonalPosteriorProbFn)(StateMachine *, int64_t, DpMatrix *, DpMatrix *, const SymbolString, const SymbolString, double,
                                                                                                PairwiseAlignmentParameters *, void *), void (*coordinateCorrectionFn)(), void *extraArgs) {
    stList *splitPoints = getSplitPoints(anchorPairs, lX, lY, p->splitMatrixBiggerThanThis, alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd);
    int64_t j = 0;
    //Now to the actual alignments
    for (int64_t i = 0; i < stList_length(splitPoints); i++) {
        stIntTuple *subRegion = stList_get(splitPoints, i);
        int64_t x1 = stIntTuple_get(subRegion, 0);
        int64_t y1 = stIntTuple_get(subRegion, 1);
        int64_t x2 = stIntTuple_get(subRegion, 2);
        int64_t y2 = stIntTuple_get(subRegion, 3);

        //Sub sequences
        char *sX2 = stString_getSubString(sX, x1, x2 - x1);
        char *sY2 = stString_getSubString(sY, y1, y2 - y1);
        SymbolString sX3 = symbolString_construct(sX2, x2 - x1);
        SymbolString sY3 = symbolString_construct(sY2, y2 - y1);

        //List of anchor pairs
        stList *subListOfAnchorPoints = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        while (j < stList_length(anchorPairs)) {
            stIntTuple *anchorPair = stList_get(anchorPairs, j);
            int64_t x = stIntTuple_get(anchorPair, 0);
            int64_t y = stIntTuple_get(anchorPair, 1);
            assert(x + y >= x1 + y1);
            if (x + y >= x2 + y2) {
                break;
            }
            assert(x >= x1 && x < x2);
            assert(y >= y1 && y < y2);
            stList_append(subListOfAnchorPoints, stIntTuple_construct3(x - x1, y - y1, stIntTuple_get(anchorPair, 2)));
            j++;
        }

        //Make the alignments
        getPosteriorProbsWithBanding(sM, subListOfAnchorPoints, sX3, sY3, p, (alignmentHasRaggedLeftEnd || i > 0),
                                     (alignmentHasRaggedRightEnd || i < stList_length(splitPoints) - 1), diagonalPosteriorProbFn, extraArgs);
        if (coordinateCorrectionFn != NULL) {
            coordinateCorrectionFn(x1, y1, extraArgs);
        }

        //Clean up
        stList_destruct(subListOfAnchorPoints);
        free(sX2);
        free(sY2);
        symbolString_destruct(sX3);
        symbolString_destruct(sY3);
    }
    assert(j == stList_length(anchorPairs));
    stList_destruct(splitPoints);
}

///////////////////////////////////
///////////////////////////////////
//Core public functions
///////////////////////////////////
///////////////////////////////////

PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters_construct() {
    PairwiseAlignmentParameters *p = st_malloc(sizeof(PairwiseAlignmentParameters));
    p->threshold = 0.01;
    p->minDiagsBetweenTraceBack = 1000;
    p->traceBackDiagonals = 40;
    p->diagonalExpansion = 20;
    p->constraintDiagonalTrim = 14;
    p->anchorMatrixBiggerThanThis = 500 * 500;
    p->repeatMaskMatrixBiggerThanThis = 500 * 500;
    p->splitMatrixBiggerThanThis = (int64_t) 3000 * 3000;
    p->alignAmbiguityCharacters = 0;
    p->gapGamma = 0.5;
    p->dynamicAnchorExpansion = 0;
#if defined(_OPENMP) && defined(PECAN_LOCK_POPEN)
    omp_init_lock(&(p->lastzLock));
#endif
    p->useMumAnchors = 1;
    p->recursiveMums = 1;
    p->k = 50;
    p->u = 1;
    p->anchorMethod = PAIRWISE_ANCHOR_LEGACY;
    p->seedHspThreshold = 2000;
    p->seedHspThresholdMin = 1000;
    p->seedRecursionDepth = 2;
    p->seedXDrop = 910;

    return p;
}

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentParameters *p) {
#if defined(_OPENMP) && defined(PECAN_LOCK_POPEN)
    omp_destroy_lock(&(p->lastzLock));
#endif
    free(p);
}

PairwiseAlignmentParameters *pairwiseAlignmentParameters_jsonParse(char *buf, size_t r) {
    // Setup parser
    jsmntok_t *tokens;
    char *js;
    int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

    PairwiseAlignmentParameters *params = pairwiseAlignmentBandingParameters_construct();

    for(int64_t tokenIndex=1; tokenIndex < tokenNumber; tokenIndex++) {
        jsmntok_t key = tokens[tokenIndex];
        char *keyString = stJson_token_tostr(js, &key);

        if (strcmp(keyString, "threshold") == 0) {
            params->threshold = stJson_parseFloat(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "minDiagsBetweenTraceBack") == 0) {
            params->minDiagsBetweenTraceBack = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "traceBackDiagonals") == 0) {
            params->traceBackDiagonals = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "diagonalExpansion") == 0) {
            params->diagonalExpansion = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "constraintDiagonalTrim") == 0) {
            params->constraintDiagonalTrim = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "anchorMatrixBiggerThanThis") == 0) {
            params->anchorMatrixBiggerThanThis = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "repeatMaskMatrixBiggerThanThis") == 0) {
            params->repeatMaskMatrixBiggerThanThis = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "splitMatrixBiggerThanThis") == 0) {
            params->splitMatrixBiggerThanThis = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "alignAmbiguityCharacters") == 0) {
            params->alignAmbiguityCharacters = stJson_parseBool(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "gapGamma") == 0) {
            params->gapGamma = stJson_parseFloat(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "dynamicAnchorExpansion") == 0) {
            params->dynamicAnchorExpansion = stJson_parseBool(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "anchorMethod") == 0) {
            params->anchorMethod = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "seedHspThreshold") == 0) {
            params->seedHspThreshold = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "seedHspThresholdMin") == 0) {
            params->seedHspThresholdMin = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "seedRecursionDepth") == 0) {
            params->seedRecursionDepth = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "seedXDrop") == 0) {
            params->seedXDrop = stJson_parseInt(js, tokens, ++tokenIndex);
        }
        else {
            st_errAbort("ERROR: Unrecognised key in pairwise alignment parameters json: %s\n", keyString);
        }
    }

    // Cleanup
    free(js);
    free(tokens);

    return params;
}

static void alignedPairCoordinateCorrectionFn(int64_t offsetX, int64_t offsetY, void *extraArgs) {
    stList *subListOfAlignedPairs = ((void **) extraArgs)[0];
    stList *alignedPairs = ((void **) extraArgs)[1];
    convertAlignedPairs(subListOfAlignedPairs, offsetX, offsetY); //Shift back the aligned pairs to the appropriate coordinates
    while (stList_length(subListOfAlignedPairs) > 0) {
        stList_append(alignedPairs, stList_pop(subListOfAlignedPairs));
    }
}

static void pairCoordinateCorrectionFn(int64_t offsetX, int64_t offsetY, void *extraArgs) {
    for(int64_t i=0; i<6; i+=2) {
        stList *subListOfPairs = ((void **) extraArgs)[i];
        stList *pairs = ((void **) extraArgs)[i+1];
        convertAlignedPairs(subListOfPairs, offsetX, offsetY); //Shift back the  pairs to the appropriate coordinates
        while (stList_length(subListOfPairs) > 0) {
            stList_append(pairs, stList_pop(subListOfPairs));
        }
    }
}

stList *getAlignedPairsUsingAnchors(StateMachine *sM, const char *sX, const char *sY, stList *anchorPairs, PairwiseAlignmentParameters *p,
                                    bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    const int64_t lX = strlen(sX);
    const int64_t lY = strlen(sY);

    //This list of pairs to be returned. Not in any order, but points must be unique
    stList *subListOfAlignedPairs = stList_construct();
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[2] = { subListOfAlignedPairs, alignedPairs };

    getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(sM, anchorPairs, sX, sY, lX, lY, p,
                                                               alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd, diagonalCalculationPosteriorMatchProbs,
                                                               alignedPairCoordinateCorrectionFn, extraArgs);

    assert(stList_length(subListOfAlignedPairs) == 0);
    stList_destruct(subListOfAlignedPairs);

    return alignedPairs;
}

void getAlignedPairsWithIndelsUsingAnchors(StateMachine *sM, const char *sX, const char *sY, stList *anchorPairs,
                                           PairwiseAlignmentParameters *p, stList **alignedPairs, stList **gapXPairs, stList **gapYPairs,
                                           bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    const int64_t lX = strlen(sX);
    const int64_t lY = strlen(sY);

    stList *subListOfAlignedPairs = stList_construct();
    *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    stList *subListOfGapXPairs = stList_construct();
    *gapXPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    stList *subListOfGapYPairs = stList_construct();
    *gapYPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    void *extraArgs[6] = { subListOfAlignedPairs, *alignedPairs,
                           subListOfGapXPairs, *gapXPairs, subListOfGapYPairs, *gapYPairs };

    getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(sM, anchorPairs, sX, sY, lX, lY, p,
                                                               alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd, diagonalCalculationPosteriorProbs,
                                                               pairCoordinateCorrectionFn, extraArgs);

    assert(stList_length(subListOfAlignedPairs) == 0);
    stList_destruct(subListOfAlignedPairs);
    assert(stList_length(subListOfGapXPairs) == 0);
    stList_destruct(subListOfGapXPairs);
    assert(stList_length(subListOfGapYPairs) == 0);
    stList_destruct(subListOfGapYPairs);
}

stList *getAlignedPairs(StateMachine *sM, const char *sX, const char *sY, PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd,
                        bool alignmentHasRaggedRightEnd) {
    stList *anchorPairs = getAnchorPairsForPairwiseAlignmentParameters(sX, sY, strlen(sX), strlen(sY), p);
    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, sX, sY, anchorPairs, p, alignmentHasRaggedLeftEnd,
                                                       alignmentHasRaggedRightEnd);
    stList_destruct(anchorPairs);
    return alignedPairs;
}

void getAlignedPairsWithIndels(StateMachine *sM, const char *sX, const char *sY, PairwiseAlignmentParameters *p,
                               stList **alignedPairs, stList **gapXPairs, stList **gapYPairs,
                               bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    stList *anchorPairs = getAnchorPairsForPairwiseAlignmentParameters(sX, sY, strlen(sX), strlen(sY), p);
    getAlignedPairsWithIndelsUsingAnchors(sM, sX, sY, anchorPairs, p,
                                          alignedPairs, gapXPairs, gapYPairs,
                                          alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd);
    stList_destruct(anchorPairs);
}

void getExpectationsUsingAnchors(StateMachine *sM, Hmm *hmmExpectations, const char *sX, const char *sY, stList *anchorPairs,
                                 PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(sM, anchorPairs, sX, sY, strlen(sX), strlen(sY), p,
                                                               alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd, diagonalCalculationExpectations, NULL,
                                                               hmmExpectations);
}

void getExpectations(StateMachine *sM, Hmm *hmmExpectations, const char *sX, const char *sY, PairwiseAlignmentParameters *p,
                     bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    stList *anchorPairs = getAnchorPairsForPairwiseAlignmentParameters(sX, sY, strlen(sX), strlen(sY), p);
    getExpectationsUsingAnchors(sM, hmmExpectations, sX, sY, anchorPairs, p, alignmentHasRaggedLeftEnd,
                                alignmentHasRaggedRightEnd);
    stList_destruct(anchorPairs);
}

/*
 * Functions for adjusting weights to account for probability of alignment to a gap.
 */

int64_t *getIndelProbabilities(stList *alignedPairs, int64_t seqLength, bool xIfTrueElseY) {
    int64_t *indelProbs = st_malloc(seqLength * sizeof(int64_t));
    for(int64_t i=0; i<seqLength; i++) {
        indelProbs[i] = PAIR_ALIGNMENT_PROB_1;
    }
    for(int64_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *j = stList_get(alignedPairs, i);
        indelProbs[stIntTuple_get(j, xIfTrueElseY ? 1 : 2)] -= stIntTuple_get(j, 0);
    }
    for(int64_t i=0; i<seqLength; i++) {
        if(indelProbs[i] < 0) {
            indelProbs[i] = 0;
        }
    }
    return indelProbs;
}

stList *reweightAlignedPairs(stList *alignedPairs,
                             int64_t *indelProbsX, int64_t *indelProbsY, double gapGamma) {
    stList *reweightedAlignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int64_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(aPair, 1);
        int64_t y = stIntTuple_get(aPair, 2);
        int64_t updatedWeight = stIntTuple_get(aPair, 0) - gapGamma * (indelProbsX[x] + indelProbsY[y]);
        stList_append(reweightedAlignedPairs, stIntTuple_construct3(updatedWeight, x, y));
    }
    stList_destruct(alignedPairs);
    return reweightedAlignedPairs;
}

stList *reweightAlignedPairs2(stList *alignedPairs, int64_t seqLengthX, int64_t seqLengthY, double gapGamma) {
    if(gapGamma <= 0.0) {
        return alignedPairs;
    }
    int64_t *indelProbsX = getIndelProbabilities(alignedPairs, seqLengthX, 1);
    int64_t *indelProbsY = getIndelProbabilities(alignedPairs, seqLengthY, 0);
    alignedPairs = reweightAlignedPairs(alignedPairs, indelProbsX, indelProbsY, gapGamma);
    free(indelProbsX);
    free(indelProbsY);
    return alignedPairs;
}

int64_t getNumberOfMatchingAlignedPairs(char *subSeqX, char *subSeqY, stList *alignedPairs) {
    int64_t matches = 0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(aPair, 1), y = stIntTuple_get(aPair, 2);
        matches += toupper(subSeqX[x]) == toupper(subSeqY[y]) && toupper(subSeqX[x]) != 'N';
    }
    return matches;
}

double scoreByIdentity(char *subSeqX, char *subSeqY, int64_t lX, int64_t lY, stList *alignedPairs) {
    int64_t matches = getNumberOfMatchingAlignedPairs(subSeqX, subSeqY, alignedPairs);
    return 100.0 * ((lX + lY) == 0 ? 0 : (2.0 * matches) / (lX + lY));
}

double scoreByIdentityIgnoringGaps(char *subSeqX, char *subSeqY, stList *alignedPairs) {
    int64_t matches = getNumberOfMatchingAlignedPairs(subSeqX, subSeqY, alignedPairs);
    return 100.0 * matches / (double) stList_length(alignedPairs);
}

static double totalScore(stList *alignedPairs) {
    double score = 0.0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        score += stIntTuple_get(aPair, 0);
    }
    return score;
}

double scoreByPosteriorProbability(int64_t lX, int64_t lY, stList *alignedPairs) {
    return 100.0 * ((lX + lY) == 0 ? 0 : (2.0 * totalScore(alignedPairs)) / ((lX + lY) * PAIR_ALIGNMENT_PROB_1));
}

double scoreByPosteriorProbabilityIgnoringGaps(stList *alignedPairs) {
    return 100.0 * totalScore(alignedPairs) / ((double) stList_length(alignedPairs) * PAIR_ALIGNMENT_PROB_1);
}

/*
 * Functions for pairwise alignment creation.
 */

static int64_t *getCumulativeGapProbs(stList *gapPairs, int64_t seqLength, bool seqXNotSeqY) {
    int64_t *gapCumulativeProbs = st_calloc(seqLength, sizeof(int64_t));

    // Work out the per-position gap probability
    for(int64_t i=0; i<stList_length(gapPairs); i++) {
        stIntTuple *gapPair = stList_get(gapPairs, i);
        assert(stIntTuple_get(gapPair, seqXNotSeqY ? 1 : 2) >= 0);
        assert(stIntTuple_get(gapPair, seqXNotSeqY ? 1 : 2) < seqLength);
        gapCumulativeProbs[stIntTuple_get(gapPair, seqXNotSeqY ? 1 : 2)] += stIntTuple_get(gapPair, 0);
    }

    // Make cumulative
    for(int64_t i=1; i<seqLength; i++) {
        gapCumulativeProbs[i] += gapCumulativeProbs[i-1];
    }

    return gapCumulativeProbs;
}

static int64_t getIndelProb(int64_t *gapCumulativeProbs, int64_t start, int64_t length) {
    assert(start >= 0);
    assert(length >= 0);
    return length == 0 ? 0 : (gapCumulativeProbs[start + length - 1] - (start > 0 ? gapCumulativeProbs[start-1] : 0));
}

stList *getMaximalExpectedAccuracyPairwiseAlignment(stList *alignedPairs,
                                                    stList *gapXPairs, stList *gapYPairs,
                                                    int64_t seqXLength, int64_t seqYLength, double *alignmentScore, PairwiseAlignmentParameters *p) {

    int64_t totalPairs = stList_length(alignedPairs); // Total number of aligned pairs

    double *scores = st_calloc(totalPairs+1, sizeof(double)); // MEA alignment score for each aligned pair
    int64_t *backPointers = st_calloc(totalPairs+1, sizeof(int64_t)); // Trace back pointers
    bool *isHighScore = st_calloc(totalPairs+1, sizeof(bool)); // Records if the score for a given aligned pair at index i is larger than any
    // score for an aligned pair at any index less than i

    // Calculate gap array cumulative probs
    int64_t *gapYCumulativeProbs = getCumulativeGapProbs(gapYPairs, seqYLength, FALSE);
    int64_t *gapXCumulativeProbs = getCumulativeGapProbs(gapXPairs, seqXLength, TRUE);

    // Iterate through the aligned pairs in order of increasing sequence coordinate

    double maxScore = 0; // Max score seen so far

    for(int64_t i=0; i<totalPairs+1; i++) {

        int64_t matchProb, x, y;

        if(i == totalPairs) { // Add final aligned pair at the end of the sequences to trace back the final alignment
            matchProb = 0; x = seqXLength; y = seqYLength;
        }
        else {
            stIntTuple *aPair = stList_get(alignedPairs, i);
            matchProb = stIntTuple_get(aPair, 0); x = stIntTuple_get(aPair, 1); y = stIntTuple_get(aPair, 2);
        }

        // The MEA alignment score of the pair with no preceding alignment pair
        double score = matchProb +
                       (getIndelProb(gapXCumulativeProbs, 0, x) + getIndelProb(gapYCumulativeProbs, 0, y)) * p->gapGamma;
        int64_t backPointer = -1;

        // Walk back through previous aligned pairs
        for(int64_t j=i-1; j>= 0; j--) {
            stIntTuple *pPair = stList_get(alignedPairs, j);
            int64_t x2 = stIntTuple_get(pPair, 1), y2 = stIntTuple_get(pPair, 2);

            // If the previous pair, pPair, and aPair can form an alignment
            if(x2 < x && y2 < y) {

                // Calc score of MEA alignment including pPair
                int64_t s = matchProb + scores[j] +
                            (getIndelProb(gapXCumulativeProbs, x2+1, x-x2-1) +
                             getIndelProb(gapYCumulativeProbs, y2+1, y-y2-1)) * p->gapGamma;

                // If score s is highest keep it
                if(s > score) {
                    score = s;
                    backPointer = j;
                }

                // If the score of pPair is a high score then can not increase score by exploring further back
                // pointers
                if(isHighScore[j]) {
                    break;
                }
            }
        }

        // Store the best alignment for aPair
        backPointers[i] = backPointer;
        scores[i] = score;

        // If the score of the alignment ending at aPair is higher than any we've seen to date
        double s = score + ((x < seqXLength ? getIndelProb(gapXCumulativeProbs, x+1, seqXLength-x-1) : 0) +
                            (y < seqYLength ? getIndelProb(gapYCumulativeProbs, y+1, seqYLength-y-1) : 0)) * p->gapGamma;
        if(s >= maxScore) {
            maxScore = s; // Record the max score
            isHighScore[i] = 1; // Record the fact that the score represents a max seen so far.
        }
    }

    // Trace back to build the MEA alignment in reverse
    stList *filteredAlignment = stList_construct3(0, (void(*)(void *))stIntTuple_destruct);
    int64_t i = backPointers[totalPairs];
    while(i >= 0) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        stList_append(filteredAlignment, stIntTuple_construct3(stIntTuple_get(aPair, 0),
                                                               stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2)));
        i = backPointers[i];
    }
    stList_reverse(filteredAlignment); // Flip the order

    // Cleanup
    free(scores);
    free(backPointers);
    free(isHighScore);
    free(gapXCumulativeProbs);
    free(gapYCumulativeProbs);

    *alignmentScore = maxScore;
    return filteredAlignment;
}

stList *leftShiftAlignment(stList *alignedPairs, char *seqX, char *seqY) {
    int64_t seqXLength = strlen(seqX), seqYLength = strlen(seqY);

    stList *leftShiftedAlignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

    int64_t x = seqXLength, y = seqYLength;
    for(int64_t i=stList_length(alignedPairs)-1; i>=0; i--) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int64_t x2 = stIntTuple_get(alignedPair, 1), y2 = stIntTuple_get(alignedPair, 2);

        while((x - x2 > 1 || y - y2 > 1) && toupper(seqX[x-1]) == toupper(seqY[y-1])) { // Insert in seqX or seqY and shift possible
            stList_append(leftShiftedAlignedPairs, stIntTuple_construct3(stIntTuple_get(alignedPair, 0), x-1, y-1)); // Hacks the score by borrowing from the current aligned pair being considered
            x--; y--;

            if(x2 == x || y2 == y) { // We've shifted over an existing aligned pair
                break;
            }
        }
        if(x2 < x && y2 < y) {
            stList_append(leftShiftedAlignedPairs, stIntTuple_construct3(stIntTuple_get(alignedPair, 0), x2, y2));
            x = x2;
            y = y2;
        }
    }

    // Deal with boundary at beginning of alignment
    while(((x > 0) && (y > 0)) && (toupper(seqX[x-1]) == toupper(seqY[y-1]))) {
        int64_t score = stList_length(alignedPairs) > 0 ? stIntTuple_get(stList_get(alignedPairs, 0), 0) : 1;
        stList_append(leftShiftedAlignedPairs, stIntTuple_construct3(score, x-1, y-1));
        x--; y--;
    }

    // Reverse, because built backwards
    stList_reverse(leftShiftedAlignedPairs);

    return leftShiftedAlignedPairs;
}

/*
 * Convenience function that aligns two sequences return a left-shift MEA alignment
 */
stList *getShiftedMEAAlignment(char *seqX, char *seqY, stList *anchorAlignment, PairwiseAlignmentParameters *p, StateMachine *sM,
                               bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd, double *alignmentScore) {
    // Generate the posterior alignment probabilities
    stList *alignedPairs, *gapXPairs, *gapYPairs;
    getAlignedPairsWithIndelsUsingAnchors(sM, seqX, seqY, anchorAlignment,
                                          p, &alignedPairs, &gapXPairs, &gapYPairs,
                                          alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd);

    // Get the MEA alignment
    stList *alignment = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
                                                                    strlen(seqX), strlen(seqY),
                                                                    alignmentScore, p);

    // Left shift the alignment
    stList *leftShiftedAlignment = leftShiftAlignment(alignment, seqX, seqY);

    // Cleanup
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    stList_destruct(alignedPairs);
    stList_destruct(alignment);

    return leftShiftedAlignment;
}

///////////////////////////////////
///////////////////////////////////
// Following functions used maximal unique matches (MUMs) to generate alignment anchors.
///////////////////////////////////
///////////////////////////////////

struct sortKmersArgs {
    const char *sequence;
    int64_t k;
};

int cmpKmers(const char *k1, const char *k2, int64_t k, int64_t *matchLength) {
    for(int64_t i=0; i<k; i++) {
        if(tolower(k1[i]) < tolower(k2[i])) {
            *matchLength = i;
            return -1;
        }
        if(tolower(k1[i]) > tolower(k2[i])) {
            *matchLength = i;
            return 1;
        }
    }
    *matchLength = k;
    return 0;
}

/* converts pointer to pointer into pointer to element */
static int sortKmersCmpFn(const void *a, const void *b, void* args) {
    struct sortKmersArgs *sargs = (struct sortKmersArgs *)args;
    int64_t m; // dummy variable used to store match length
    return cmpKmers(&(sargs->sequence[*(int64_t *)a]), &(sargs->sequence[*(int64_t *)b]), sargs->k, &m);
}

/*
 * Simple/crappy method to build a suffix array like object.
 */
static int64_t *getSortedKmers(const char *sequence, const int64_t length, const int64_t k) {
    if(length - k + 1 <= 0) {
        return NULL;
    }
    int64_t *sortedKmers = st_malloc(sizeof(int64_t) * (length - k + 1));

    // Initialize the array with the start indices of the kmers in the sequence
    for(int64_t i=0; i<length-k+1; i++) {
        sortedKmers[i] = i;
    }

    // Now sort the kmers by kmer content
    struct sortKmersArgs args = {sequence, k};
    safesort(sortedKmers, length-k+1, sizeof(int64_t), sortKmersCmpFn, &args);

    return sortedKmers;
}

static int64_t getMatchLength(int64_t *sortedTargetSuffixes, int64_t targetLength, const char *targetSequence,
                              const char *searchSequence, int64_t k, int64_t m) {
    if(m < 0 || m >= targetLength) {
        return 0;
    }
    cmpKmers(searchSequence, &targetSequence[sortedTargetSuffixes[m]], k, &m);
    return m;
}

/*
 * Returns the index of the unique longest match in the target sequence, if it exists. If there is no
 * unique longest match then returns -1.
 */
int64_t getLongestUniqueMatch(int64_t *sortedTargetSuffixes, int64_t targetLength, const char *targetSequence,
                              const char *searchSequence, int64_t k, int64_t u, int64_t *matchLength) {
    // First use binary search to find the index of the longest match in the target
    int64_t l=0, h=targetLength, p=-2; // interval (l, h) that item can be in, l is inclusive, h is exclusive
    *matchLength = 0; // Set the match length to be 0, initially
    while(l < h) {
        int64_t m = (l + h) / 2, n; // Mid point
        int64_t i = cmpKmers(searchSequence, &targetSequence[sortedTargetSuffixes[m]], k, &n);
        if(n > *matchLength) { // Update the longest match seen so far
            *matchLength = n;
            p = m;
        }
        if(i < 0) { // Item must occur before m in the list
            h = m;
        }
        else if(i > 0) { // Item must occur after m in the list
            l = m+1;
        } else { // else item at index i equals i, so break
            break;
        }
    }

    // Now check the match is longer than the surrounding matches
    return (*matchLength > u + getMatchLength(sortedTargetSuffixes, targetLength, targetSequence, searchSequence, k, p - 1) &&
            *matchLength > u + getMatchLength(sortedTargetSuffixes, targetLength, targetSequence, searchSequence, k, p + 1)) ? p : -1;
}

/*
 * Struct to represent an aligned mum.
 */
typedef struct _mum Mum;
struct _mum {
    int64_t x, y, length, score, refCount;
    Mum *pMum; // The previous mum in the alignment
};

static Mum *mum_construct(int64_t x, int64_t y, int64_t length) {
    Mum *mum = st_calloc(1, sizeof(Mum));
    mum->x = x;
    mum->y = y;
    mum->length = length;
    mum->refCount = 1;
    return mum;
}

void mum_destruct(Mum *mum) {
    // A loop, not recursion: a chain of mums along a megabase is tens of thousands long
    while (mum != NULL) {
        assert(mum->refCount > 0);
        if (--mum->refCount > 0) {
            return;
        }
        Mum *pMum = mum->pMum;
        free(mum);
        mum = pMum;
    }
}

/*
 * Compare mums by y start-coordinate.
 */
static int mum_sweep_cmp(const void *a, const void *b) {
    int64_t y1 = ((Mum *)a)->y + ((Mum *)a)->length, y2 = ((Mum *)b)->y + ((Mum *)b)->length;
    return y1 > y2 ? 1 : (y1 < y2 ? -1 : 0);
}

/*
 * Adds new mums the sweep line, destroys the mumsToAdd list and returns an updated mumsToAdd list in the
 * process.
 */
stList *updateSweepLine(stSortedSet *sweepLine, stList *mumsToAdd, int64_t x) {
    stList *mumsToAdd2 = stList_construct();
    for(int64_t i=0; i<stList_length(mumsToAdd); i++) {
        Mum *mum = stList_get(mumsToAdd, i);

        // Skip the mum if it doesn't end at x
        if(mum->x + mum->length != x) {
            stList_append(mumsToAdd2, mum);
            continue;
        }

        // Check if it is lower scoring than an existing alignment that precedes it on the sweep line
        Mum *mum2 = stSortedSet_searchLessThanOrEqual(sweepLine, mum);
        if(mum2 != NULL && mum2->score >= mum->score) {
            // In this case we don't need to include it, so we delete it and continue
            mum_destruct(mum);
            continue;
        }

        // Remove any points that are eclipsed by the new mum
        while(1) {
            Mum *mum2 = stSortedSet_searchGreaterThanOrEqual(sweepLine, mum);
            if(mum2 == NULL || mum2->score > mum->score) { // If there is no further mum, or its score is larger
                // than the mum being added
                break;
            }
            mum_destruct(stSortedSet_remove(sweepLine, mum2));
        }

        // Now we can add the new mum to the sweep line
        stSortedSet_insert(sweepLine, mum);
    }
    stList_destruct(mumsToAdd);
    return mumsToAdd2;
}

/*
 * Chains the mum with the highest scoring preceding mum
 * Does not currently account for gap penalties.
 */
static void chainMum(stSortedSet *sweepLine, Mum *mum) {
    Mum m; m.length = 0; m.y = mum->y;
    mum->pMum = stSortedSet_searchLessThan(sweepLine, &m);
    if(mum->pMum == NULL) {
        mum->score = mum->length;
        return;
    }
    assert(mum->x - mum->pMum->x - mum->pMum->length >= 0);
    assert(mum->y - mum->pMum->y - mum->pMum->length >= 0);
    mum->pMum->refCount++;
    mum->score = mum->pMum->score + mum->length;
}

static void getAlignedMums2(const char *sX, const char *sY, int64_t lX, int64_t lY, PairwiseAlignmentParameters *p,
                            int64_t oX, int64_t oY, bool recursive, stList *alignedPairs);

/*
 * Traces back the highest scoring alignment of mums to create set of aligned pairs.
 */
static void tracebackMums(stSortedSet *sweepLine,
                          const char *sX, const char *sY, int64_t lX, int64_t lY,
                          PairwiseAlignmentParameters *p,
                          int64_t oX, int64_t oY, bool recursive, stList *alignedPairs) {
    Mum *mum = stSortedSet_getLast(sweepLine);
    while(mum != NULL) { // traceback from last-to-first mum
        if(recursive) {
            int64_t i = mum->x + mum->length;
            int64_t j = mum->y + mum->length;
            assert(lX - i >= 0 && lY - j >= 0);
            if((lX - i) * (lY - j) > p->anchorMatrixBiggerThanThis) {
                getAlignedMums2(&sX[i], &sY[j], lX-i, lY-j, p, oX + i, oY + j, 0, alignedPairs);
            }
        }
        for(int64_t i=mum->length-1; i>=0; i--) { // go backwards, so that when we reverse order pairs will be increasing
            stList_append(alignedPairs, stIntTuple_construct3(oX + mum->x+i, oY + mum->y+i, p->diagonalExpansion));
        }
        lX = mum->x; lY = mum->y;
        mum = mum->pMum;
    }
    if(recursive) {
        if (lX * lY > p->anchorMatrixBiggerThanThis) { // Deal with the gap between the start and the first mum
            getAlignedMums2(sX, sY, lX, lY, p, oX, oY, 0, alignedPairs);
        }
        stList_reverse(alignedPairs); // reverse to get in first-to-last order, but only at the top-level, so we don't
        // do it twice
    }
}

static void getAlignedMums2(const char *sX, const char *sY, int64_t lX, int64_t lY, PairwiseAlignmentParameters *p,
                            int64_t oX, int64_t oY, bool recursive, stList *alignedPairs) {
    // Get the kmers in each sequence sorted lexicographically
    int64_t *sortedKmersY = getSortedKmers(sY, lY, p->k);

    stSortedSet *sweepLine = stSortedSet_construct3(mum_sweep_cmp, (void (*)(void *))mum_destruct);
    stList *mumsToAdd = stList_construct();

    int64_t pDiag=-1, pXEnd=-1;
    for(int64_t x=0; x<lX-p->k+1; x++) {
        // Add any mums that end at x to the sweep line
        mumsToAdd = updateSweepLine(sweepLine, mumsToAdd, x);

        // Get the next longest unique match
        int64_t matchLength;
        int64_t j = getLongestUniqueMatch(sortedKmersY, lY-p->k+1, sY, &sX[x], p->k, p->u, &matchLength);

        // If is a longest unique match
        if(j >= 0 && j < lY-p->k+1) {
            int64_t y = sortedKmersY[j]; // The y coordinate of the start of the match
            if(pXEnd < x || pDiag != x - y) { // If not part of the previous match then is a MUM
                // Make a mum
                Mum *mum = mum_construct(x, y, matchLength);

                // Add the mum to it's best alignment
                chainMum(sweepLine, mum);

                //st_uglyf("Got mum, x: %i, y: %i, length: %i, pMum: %i\n", (int)mum->x, (int)mum->y,
                //         (int)mum->length, (int)mum->pMum);

                // Queue the mum up to add the sweep line when we get to its end x-coordinate
                stList_append(mumsToAdd, mum);

                // Update the coordinates of the last found mum
                pDiag = x - y;
                pXEnd = x + matchLength;
            }
        }
    }

    for(int64_t x=lX-p->k+1; x<lX; x++) {
        // Finish adding the mums to the sweep line
        mumsToAdd = updateSweepLine(sweepLine, mumsToAdd, x);
    }

    // Get highest scoring chain
    tracebackMums(sweepLine, sX, sY, lX, lY, p, oX, oY, recursive, alignedPairs);

    // cleanup
    free(sortedKmersY);
    assert(stList_length(mumsToAdd) == 0);
    stList_destruct(mumsToAdd);
    stSortedSet_destruct(sweepLine);
}

stList *getAlignedMums(const char *sX, const char *sY, int64_t lX, int64_t lY, PairwiseAlignmentParameters *p,
                       int64_t oX, int64_t oY) {
    stList *alignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    getAlignedMums2(sX, sY, lX, lY, p, oX, oY, p->recursiveMums, alignedPairs);
    return alignedPairs;
}


///////////////////////////////////
///////////////////////////////////
//Anchoring by spaced seeds, ungapped extension and chaining
//
//What the lastz call above does -- 12of19 spaced seeds allowing one transition, HOXD70 scored
//ungapped extension, and the heaviest colinear chain of the resulting HSPs -- but in process.
//lastz costs a process, two temporary files and ~50 ms of startup for every pairwise alignment
//it anchors.  The MUM anchors, which were the way around that, accept "unique" matches as short
//as a dozen bases, and at the divergence of a deep branch most of those are random: simulated
//at 0.3 substitutions per site per branch, half the MUM anchors were off the true alignment,
//and the banded dp, which cannot leave the band they define, followed them.
///////////////////////////////////
///////////////////////////////////

#define SEED_SPAN 19
#define SEED_WEIGHT 12
#define SEED_MAX_WORD_COUNT 64 // Y positions a seed word may have before it is taken to be a repeat and skipped

static const int64_t seedCare[SEED_WEIGHT] = { 0, 1, 2, 4, 7, 8, 11, 13, 15, 16, 17, 18 }; // 1110100110010101111

static const int64_t seedHoxd70[5][5] = { // A C G T, then N against anything
        {   91, -114,  -31, -123, -100 },
        { -114,  100, -125,  -31, -100 },
        {  -31, -125,  100, -114, -100 },
        { -123,  -31, -114,   91, -100 },
        { -100, -100, -100, -100, -100 } };

typedef struct _seedSequence {
    uint8_t *codes; // 0-3 for ACGT, 4 for anything else
    int32_t *masked; // masked[i] is the number of lower case bases before i
} SeedSequence;

static SeedSequence seedSequence_construct(const char *s, int64_t l) {
    SeedSequence seq;
    seq.codes = st_malloc(l + 1);
    seq.masked = st_malloc(sizeof(int32_t) * (l + 1));
    seq.masked[0] = 0;
    for (int64_t i = 0; i < l; i++) {
        switch (s[i]) {
            case 'A': case 'a': seq.codes[i] = 0; break;
            case 'C': case 'c': seq.codes[i] = 1; break;
            case 'G': case 'g': seq.codes[i] = 2; break;
            case 'T': case 't': seq.codes[i] = 3; break;
            default: seq.codes[i] = 4;
        }
        seq.masked[i + 1] = seq.masked[i] + (islower((unsigned char) s[i]) ? 1 : 0);
    }
    return seq;
}

static void seedSequence_destruct(SeedSequence seq) {
    free(seq.codes);
    free(seq.masked);
}

/*
 * The seed word at i, over the seed's care positions.  False if the span has an N, or, when
 * repeat masking, a lower case base.
 */
static inline bool seed_word(const SeedSequence *seq, int64_t i, bool repeatMask, uint32_t *word) {
    if (repeatMask && seq->masked[i + SEED_SPAN] != seq->masked[i]) {
        return 0;
    }
    uint32_t w = 0;
    for (int64_t j = 0; j < SEED_WEIGHT; j++) {
        uint8_t c = seq->codes[i + seedCare[j]];
        if (c > 3) {
            return 0;
        }
        w |= ((uint32_t) c) << (2 * j);
    }
    *word = w;
    return 1;
}

typedef struct _seedEntry {
    uint32_t word;
    int32_t pos;
} SeedEntry;

static int seedEntry_cmp(const void *a, const void *b) {
    const SeedEntry *e = a, *f = b;
    if (e->word != f->word) {
        return e->word < f->word ? -1 : 1;
    }
    return e->pos < f->pos ? -1 : (e->pos > f->pos ? 1 : 0);
}

typedef struct _seedHsp {
    int64_t x, y, length, score;
} SeedHsp;

/*
 * Extends a seed hit at (x, y) both ways along its diagonal, stopping each way once the score
 * falls xDrop below the best seen, and returns the best scoring segment.
 */
static SeedHsp seed_extend(const uint8_t *cX, int64_t x0, int64_t x1, const uint8_t *cY, int64_t y0, int64_t y1,
                           int64_t x, int64_t y, int64_t xDrop) {
    int64_t s = 0, best = 0, right = 0;
    for (int64_t k = 0; x + k < x1 && y + k < y1; k++) {
        s += seedHoxd70[cX[x + k]][cY[y + k]];
        if (s > best) {
            best = s;
            right = k + 1;
        } else if (s < best - xDrop) {
            break;
        }
    }
    int64_t s2 = 0, best2 = 0, left = 0;
    for (int64_t k = 1; x - k >= x0 && y - k >= y0; k++) {
        s2 += seedHoxd70[cX[x - k]][cY[y - k]];
        if (s2 > best2) {
            best2 = s2;
            left = k;
        } else if (s2 < best2 - xDrop) {
            break;
        }
    }
    SeedHsp hsp = { x - left, y - left, left + right, best + best2 };
    return hsp;
}

/*
 * The HSPs scoring at least threshold between X[x0, x1) and Y[y0, y1).
 */
static SeedHsp *seed_getHsps(SeedSequence *sX, int64_t x0, int64_t x1, SeedSequence *sY, int64_t y0, int64_t y1,
                             bool repeatMask, int64_t threshold, int64_t xDrop, int64_t *hspNumber) {
    *hspNumber = 0;
    int64_t lX = x1 - x0, lY = y1 - y0;
    if (lX < SEED_SPAN || lY < SEED_SPAN) {
        return NULL;
    }
    // Index the words of Y, sorted, with a table of where each bucket of leading word bits starts
    SeedEntry *entries = st_malloc(sizeof(SeedEntry) * (lY - SEED_SPAN + 1));
    int64_t entryNumber = 0;
    for (int64_t y = y0; y <= y1 - SEED_SPAN; y++) {
        uint32_t w;
        if (seed_word(sY, y, repeatMask, &w)) {
            entries[entryNumber].word = w;
            entries[entryNumber++].pos = (int32_t) y;
        }
    }
    qsort(entries, entryNumber, sizeof(SeedEntry), seedEntry_cmp);
    int64_t bucketBits = 4;
    while (bucketBits < 20 && ((int64_t) 1 << bucketBits) < entryNumber) {
        bucketBits++;
    }
    int64_t shift = 2 * SEED_WEIGHT - bucketBits;
    int64_t *bucketStart = st_malloc(sizeof(int64_t) * (((int64_t) 1 << bucketBits) + 1));
    for (int64_t b = 0, e = 0; b <= ((int64_t) 1 << bucketBits); b++) {
        while (e < entryNumber && (int64_t) (entries[e].word >> shift) < b) {
            e++;
        }
        bucketStart[b] = e;
    }

    // Walk X, looking up each word and the words one transition away from it.  A hit on a
    // diagonal already extended past it is part of an HSP that has been found.
    int32_t *diagonalEnd = st_malloc(sizeof(int32_t) * (lX + lY));
    for (int64_t i = 0; i < lX + lY; i++) {
        diagonalEnd[i] = INT32_MIN;
    }
    int64_t hspCapacity = 64;
    SeedHsp *hsps = st_malloc(sizeof(SeedHsp) * hspCapacity);
    for (int64_t x = x0; x <= x1 - SEED_SPAN; x++) {
        uint32_t w;
        if (!seed_word(sX, x, repeatMask, &w)) {
            continue;
        }
        for (int64_t v = 0; v <= SEED_WEIGHT; v++) {
            uint32_t w2 = v == 0 ? w : w ^ (2u << (2 * (v - 1))); // A<->G and C<->T are the same bit flip
            int64_t b = w2 >> shift, lo = bucketStart[b], hi = bucketStart[b + 1];
            while (lo < hi && entries[lo].word < w2) {
                lo++;
            }
            int64_t e = lo;
            while (e < hi && entries[e].word == w2) {
                e++;
            }
            if (e - lo > SEED_MAX_WORD_COUNT) {
                continue;
            }
            for (int64_t k = lo; k < e; k++) {
                int64_t y = entries[k].pos;
                int64_t d = (x - x0) - (y - y0) + lY - 1;
                if (diagonalEnd[d] > x) {
                    continue;
                }
                SeedHsp hsp = seed_extend(sX->codes, x0, x1, sY->codes, y0, y1, x, y, xDrop);
                diagonalEnd[d] = hsp.x + hsp.length > x + 1 ? (int32_t) (hsp.x + hsp.length) : (int32_t) (x + 1);
                if (hsp.length > 0 && hsp.score >= threshold) { // a zero-length "HSP", possible only with a threshold of 0, would break the chain's sweep
                    if (*hspNumber == hspCapacity) {
                        hspCapacity *= 2;
                        hsps = st_realloc(hsps, sizeof(SeedHsp) * hspCapacity);
                    }
                    hsps[(*hspNumber)++] = hsp;
                }
            }
        }
    }
    free(entries);
    free(bucketStart);
    free(diagonalEnd);
    return hsps;
}

static int seedHsp_cmpByX(const void *a, const void *b) {
    const SeedHsp *h = a, *i = b;
    if (h->x != i->x) {
        return h->x < i->x ? -1 : 1;
    }
    return h->y < i->y ? -1 : (h->y > i->y ? 1 : 0);
}

static int seed_cmpInt64(const void *a, const void *b) {
    int64_t k = *(int64_t *) a, l = *(int64_t *) b;
    return k < l ? -1 : (k > l ? 1 : 0);
}

/*
 * Replaces hsps with the heaviest chain of them that is colinear and non-overlapping in both
 * sequences, in order.  A sweep along X, with a Fenwick tree of the best chain ending at or
 * before each Y position.
 */
static void seed_chain(SeedHsp *hsps, int64_t *hspNumber) {
    int64_t n = *hspNumber;
    if (n <= 1) {
        return;
    }
    qsort(hsps, n, sizeof(SeedHsp), seedHsp_cmpByX);
    int64_t *ends = st_malloc(sizeof(int64_t) * n), *byEnd = st_malloc(sizeof(int64_t) * n);
    int64_t *yEnds = st_malloc(sizeof(int64_t) * n);
    for (int64_t i = 0; i < n; i++) {
        ends[i] = hsps[i].x + hsps[i].length;
        byEnd[i] = i;
        yEnds[i] = hsps[i].y + hsps[i].length;
    }
    // The HSPs in order of where they end in X, by sorting (end, index) pairs
    int64_t *pairs = st_malloc(sizeof(int64_t) * 2 * n);
    for (int64_t i = 0; i < n; i++) {
        pairs[2 * i] = ends[i];
        pairs[2 * i + 1] = i;
    }
    qsort(pairs, n, 2 * sizeof(int64_t), seed_cmpInt64); // compares the end, the first of each pair
    for (int64_t i = 0; i < n; i++) {
        byEnd[i] = pairs[2 * i + 1];
    }
    free(pairs);
    qsort(yEnds, n, sizeof(int64_t), seed_cmpInt64);
    int64_t m = 0;
    for (int64_t i = 0; i < n; i++) {
        if (m == 0 || yEnds[m - 1] != yEnds[i]) {
            yEnds[m++] = yEnds[i];
        }
    }
    int64_t *fenwickScore = st_calloc(m + 1, sizeof(int64_t)), *fenwickIndex = st_malloc(sizeof(int64_t) * (m + 1));
    int64_t *best = st_malloc(sizeof(int64_t) * n), *prev = st_malloc(sizeof(int64_t) * n);
    int64_t e = 0, bestEnd = -1;
    for (int64_t i = 0; i < n; i++) {
        // Make the chains ending at or before this HSP's start in X available
        while (e < n && ends[byEnd[e]] <= hsps[i].x) {
            int64_t j = byEnd[e++];
            int64_t yEnd = hsps[j].y + hsps[j].length;
            int64_t lo = 0, hi = m; // rank of yEnd, 1 based
            while (lo < hi) {
                int64_t mid = (lo + hi) / 2;
                if (yEnds[mid] < yEnd) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            for (int64_t r = lo + 1; r <= m; r += r & -r) {
                if (best[j] > fenwickScore[r]) {
                    fenwickScore[r] = best[j];
                    fenwickIndex[r] = j;
                }
            }
        }
        // The best of them ending at or before its start in Y
        int64_t lo = 0, hi = m; // number of yEnds <= hsps[i].y
        while (lo < hi) {
            int64_t mid = (lo + hi) / 2;
            if (yEnds[mid] <= hsps[i].y) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        int64_t s = 0, p = -1;
        for (int64_t r = lo; r > 0; r -= r & -r) {
            if (fenwickScore[r] > s) {
                s = fenwickScore[r];
                p = fenwickIndex[r];
            }
        }
        best[i] = hsps[i].score + s;
        prev[i] = p;
        if (bestEnd == -1 || best[i] > best[bestEnd]) {
            bestEnd = i;
        }
    }
    int64_t chainLength = 0;
    for (int64_t i = bestEnd; i != -1; i = prev[i]) {
        chainLength++;
    }
    SeedHsp *chain = st_malloc(sizeof(SeedHsp) * chainLength);
    for (int64_t i = bestEnd, k = chainLength - 1; i != -1; i = prev[i], k--) {
        chain[k] = hsps[i];
    }
    memcpy(hsps, chain, sizeof(SeedHsp) * chainLength);
    *hspNumber = chainLength;
    free(chain);
    free(ends);
    free(byEnd);
    free(yEnds);
    free(fenwickScore);
    free(fenwickIndex);
    free(best);
    free(prev);
}

/*
 * The HSP score needed in a box of the given area.  How high a score chance alone reaches grows
 * with the log of the area searched, so a whole megabase-square problem needs a high one -- on
 * unrelated random sequence, 2000 still lets 16 HSPs into the chain of a 1 Mb x 1 Mb search, and
 * 2500 none -- while the gap between two anchors, a few hundred bases square and already known to
 * sit between homologous sequence, can take a much lower one: 1000 finds nothing by chance in a
 * 1 kb square.  Rises 250 per tenfold of area from minThreshold at 1e6, capped at maxThreshold.
 */
static int64_t seed_threshold(int64_t area, int64_t minThreshold, int64_t maxThreshold) {
    double t = minThreshold + (area > 1000000 ? 250.0 * log10((double) area / 1000000.0) : 0.0);
    return t > maxThreshold ? maxThreshold : (int64_t) t;
}

/*
 * Appends to anchors, in order, the anchors between X[x0, x1) and Y[y0, y1): the chained HSPs,
 * less constraintDiagonalTrim at each end, and, down to seedRecursionDepth levels, the anchors
 * of the gaps between them that are big enough to be worth anchoring themselves.  The whole box
 * takes seedHspThreshold; the gaps within it, flanked by anchors, the lower seed_threshold.
 */
static void seed_anchorBox(SeedSequence *sX, int64_t x0, int64_t x1, SeedSequence *sY, int64_t y0, int64_t y1,
                           PairwiseAlignmentParameters *p, bool repeatMask, int64_t depth, stList *anchors) {
    int64_t hspNumber;
    int64_t threshold = depth == 0 ? p->seedHspThreshold :
                        seed_threshold((x1 - x0) * (y1 - y0), p->seedHspThresholdMin, p->seedHspThreshold);
    SeedHsp *hsps = seed_getHsps(sX, x0, x1, sY, y0, y1, repeatMask, threshold, p->seedXDrop, &hspNumber);
    seed_chain(hsps, &hspNumber);
    int64_t pX = x0, pY = y0;
    for (int64_t i = 0; i <= hspNumber; i++) {
        int64_t nX = i < hspNumber ? hsps[i].x : x1, nY = i < hspNumber ? hsps[i].y : y1;
        int64_t gapSize = (nX - pX) * (nY - pY);
        if (depth < p->seedRecursionDepth && gapSize > p->anchorMatrixBiggerThanThis && hspNumber > 0) { // not the whole box again
            seed_anchorBox(sX, pX, nX, sY, pY, nY, p, gapSize > p->repeatMaskMatrixBiggerThanThis, depth + 1, anchors);
        }
        if (i < hspNumber) {
            for (int64_t k = p->constraintDiagonalTrim; k < hsps[i].length - p->constraintDiagonalTrim; k++) {
                stList_append(anchors, stIntTuple_construct3(hsps[i].x + k, hsps[i].y + k, p->diagonalExpansion));
            }
            pX = hsps[i].x + hsps[i].length;
            pY = hsps[i].y + hsps[i].length;
        }
    }
    free(hsps);
}

stList *getSeedAnchors(const char *sX, const char *sY, int64_t lX, int64_t lY, PairwiseAlignmentParameters *p) {
    stList *anchors = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    if (lX >= INT32_MAX || lY >= INT32_MAX) { // the index keeps positions in 32 bits; bar never asks for more than a megabase
        st_logCritical("getSeedAnchors: sequences of %" PRIi64 " and %" PRIi64 " bases are too long to anchor\n", lX, lY);
        return anchors;
    }
    SeedSequence seqX = seedSequence_construct(sX, lX), seqY = seedSequence_construct(sY, lY);
    // Masked or not by the same rule as the gaps within: repeatMaskMatrixBiggerThanThis
    seed_anchorBox(&seqX, 0, lX, &seqY, 0, lY, p, lX * lY > p->repeatMaskMatrixBiggerThanThis, 0, anchors);
    seedSequence_destruct(seqX);
    seedSequence_destruct(seqY);
    return anchors;
}
