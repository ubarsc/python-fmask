/* This file is part of 'python-fmask' - a cloud masking module
* Copyright (C) 2015  Neil Flood
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 3
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
/*
Module to implement filling of local minima in a raster surface. 

The algorithm is from 
    Soille, P., and Gratin, C. (1994). An efficient algorithm for drainage network
        extraction on DEMs. J. Visual Communication and Image Representation. 
        5(2). 181-189. 
        
The algorithm is intended for hydrological processing of a DEM, but is used by the 
Fmask cloud shadow algorithm as part of its process for finding local minima which 
represent potential shadow objects. 
*/
#include <Python.h>
#include "numpy/arrayobject.h"
#include <math.h>

/* An exception object for this module */
/* created in the init function */
struct FillMinimaState
{
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct FillMinimaState*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct FillMinimaState _state;
#endif

/* Routines for handling the hierarchical pixel queue which the
   algorithm requires.
*/
typedef struct PQstruct {
    int i, j;
    struct PQstruct *next;
} PQel;
    
typedef struct {
    PQel *first, *last;
    int n;
} PQhdr;
    
typedef struct PQ {
    int hMin;
    int numLevels;
    PQhdr *q;
} PixelQueue;
    
/* A new pixel structure */
static PQel *newPix(int i, int j) {
    PQel *p;
        
    p = (PQel *)calloc(1, sizeof(PQel));
    p->i = i;
    p->j = j;
    p->next = NULL;
    if (i>20000) {
        printf("i=%d\\n", i);
        exit(1);
    }
    return p;
}
    
/* Initialize pixel queue */
static PixelQueue *PQ_init(int hMin, int hMax) {
    PixelQueue *pixQ;
    int numLevels, i;
        
    pixQ = (PixelQueue *)calloc(1, sizeof(PixelQueue));
    numLevels = hMax - hMin + 1;
    pixQ->hMin = hMin;
    pixQ->numLevels = numLevels;
    pixQ->q = (PQhdr *)calloc(numLevels, sizeof(PQhdr));
    for (i=0; i<numLevels; i++) {
        pixQ->q[i].first = NULL;
        pixQ->q[i].last = NULL;
        pixQ->q[i].n = 0;
    }
    return pixQ;
}
    
/* Add a pixel at level h */
static void PQ_add(PixelQueue *pixQ, PQel *p, int h) {
    int ndx;
    PQel *current, *newP;
    PQhdr *thisQ;
        
    /* Take a copy of the pixel structure */
    newP = newPix(p->i, p->j);
        
    ndx = h - pixQ->hMin;
    if (ndx > pixQ->numLevels) {
        printf("Level h=%d too large. ndx=%d, numLevels=%d\\n", h, ndx, pixQ->numLevels);
        exit(1);
    }
    if (ndx < 0) {
        printf("Ndx is negative, which is not allowed. ndx=%d, h=%d, hMin=%d\n", ndx, h, pixQ->hMin);
        exit(1);
    }
    thisQ = &(pixQ->q[ndx]);
    /* Add to end of queue at this level */
    current = thisQ->last;
    if (current != NULL) {
        current->next = newP;
    }
    thisQ->last = newP;
    thisQ->n++;
    /* If head of queue is NULL, make the new one the head */
    if (thisQ->first == NULL) {
        thisQ->first = newP;
    }
}
    
/* Return TRUE if queue at level h is empty */
static int PQ_empty(PixelQueue *pixQ, int h) {
    int ndx, empty, n;
    PQel *current;
        
    ndx = h - pixQ->hMin;
    current = pixQ->q[ndx].first;
    n = pixQ->q[ndx].n;
    empty = (current == NULL);
    if (empty && (n != 0)) {
        printf("Empty, but n=%d\\n", n);
        exit(1);
    }
    if ((n == 0) && (! empty)) {
        printf("n=0, but not empty\\n");
        while (current != NULL) {
            printf("    h=%d i=%d j=%d\\n", h, current->i, current->j);
            current = current->next;
        }
        exit(1);
    }
    return empty;
}
    
/* Return the first element in the queue at level h, and remove it
   from the queue */
static PQel *PQ_first(PixelQueue *pixQ, int h) {
    int ndx;
    PQel *current;
    PQhdr *thisQ;
        
    ndx = h - pixQ->hMin;
    thisQ = &(pixQ->q[ndx]);
    current = thisQ->first;
    /* Remove from head of queue */
    if (current != NULL) {
        thisQ->first = current->next;
        if (thisQ->first == NULL) {
            thisQ->last = NULL;
        }
        thisQ->n--;
        if (thisQ->n < 0) {
            printf("n=%d in PQ_first()\\n", thisQ->n);
            exit(1);
        } else if (thisQ->n == 0) {
            if (thisQ->first != NULL) {
                printf("n=0, but 'first' != NULL. first(i,j) = %d,%d\\n", 
                    thisQ->first->i, thisQ->first->j);
            }
        }
    }
    return current;
}
    
/* Return a list of neighbouring pixels to given pixel p.  */
static PQel *neighbours(PQel *p, int nRows, int nCols) {
    int ii, jj, i, j;
    PQel *pl, *pNew;
        
    pl = NULL;
    for (ii=-1; ii<=1; ii++) {
        for (jj=-1; jj<=1; jj++) {
            if ((ii != 0) && (jj != 0)) {
                i = p->i + ii;
                j = p->j + jj;
                if ((i >= 0) && (i < nRows) && (j >= 0) && (j < nCols)) {
                    pNew = newPix(i, j);
                    pNew->next = pl;
                    pl = pNew;
                }
            }
        }
    }
    return pl;
}

#define max(a,b) ((a) > (b) ? (a) : (b))

static PyObject *fillminima_fillMinima(PyObject *self, PyObject *args)
{
    PyArrayObject *pimg, *pimg2, *pBoundaryRows, *pBoundaryCols, *pNullMask;
    int hMin, hMax;
    double dBoundaryVal;

    npy_int64 r, c;
    npy_intp i, nRows, nCols;
    npy_int16 imgval, img2val;
    PixelQueue *pixQ;
    PQel *p, *nbrs, *pNbr, *pNext;
    int hCrt;
    
    if( !PyArg_ParseTuple(args, "OOiiOdOO:fillMinima", &pimg, &pimg2, &hMin, &hMax,
                            &pNullMask, &dBoundaryVal, &pBoundaryRows, &pBoundaryCols))
        return NULL;
    
    if( !PyArray_Check(pimg) || !PyArray_Check(pimg2) || !PyArray_Check(pBoundaryRows) ||
        !PyArray_Check(pBoundaryCols) || !PyArray_Check(pNullMask) )
    {
        PyErr_SetString(GETSTATE(self)->error, "parameters 0, 1, 4, 6 and 7 must be numpy arrays");
        return NULL;
    }

    if( (PyArray_TYPE(pimg) != NPY_INT16) || (PyArray_TYPE(pimg2) != NPY_INT16))
    {
        PyErr_SetString(GETSTATE(self)->error, "parameters 0 and 1 must be int16 arrays");
        return NULL;
    }

    if( PyArray_TYPE(pNullMask) != NPY_BOOL )
    {
        PyErr_SetString(GETSTATE(self)->error, "parameter 4 must be a bool array");
        return NULL;
    }

    if( (PyArray_TYPE(pBoundaryRows) != NPY_INT64) || (PyArray_TYPE(pBoundaryCols) != NPY_INT64) )
    {
        PyErr_SetString(GETSTATE(self)->error, "parameters 6 and 7 must be int64 arrays");
        return NULL;
    }

    nRows = PyArray_DIMS(pimg)[0];
    nCols = PyArray_DIMS(pimg)[1];
    
    pixQ = PQ_init(hMin, hMax);
    
    /* Initialize the boundary */
    for (i=0; i<PyArray_DIMS(pBoundaryRows)[0]; i++) {
        r = *((npy_int*)PyArray_GETPTR1(pBoundaryRows, i));
        c = *((npy_int*)PyArray_GETPTR1(pBoundaryCols, i));
        //img2(r, c) = img(r, c);
        *((npy_int16*)PyArray_GETPTR2(pimg2, r, c)) = dBoundaryVal;
        
        p = newPix(r, c);
        //h = img(r, c);
        //PQ_add(pixQ, p, img(r, c));
        PQ_add(pixQ, p, dBoundaryVal);
    }
    
    
    /* Process until stability */
    hCrt = (int)hMin;
    do {
        while (! PQ_empty(pixQ, hCrt)) {
            p = PQ_first(pixQ, hCrt);
            nbrs = neighbours(p, nRows, nCols);
            pNbr = nbrs;
            while (pNbr != NULL) {
                r = pNbr->i;
                c = pNbr->j;
                /* Exclude null area of original image */
                if (! *((npy_bool*)PyArray_GETPTR2(pNullMask, r, c))) {
                    imgval = *((npy_int16*)PyArray_GETPTR2(pimg, r, c));
                    img2val = *((npy_int16*)PyArray_GETPTR2(pimg2, r, c));
                    if (img2val == hMax) {
                        img2val = max(hCrt, imgval);
                        *((npy_int16*)PyArray_GETPTR2(pimg2, r, c)) = img2val;
                        PQ_add(pixQ, pNbr, img2val);
                    }
                }
                pNext = pNbr->next;
                free(pNbr);
                pNbr = pNext;
            }
            free(p);
        }
        hCrt++;
    } while (hCrt < hMax);
    
    free(pixQ);

    Py_RETURN_NONE;
}


// Our list of functions in this module
static PyMethodDef FillMinimaMethods[] = {
    {"fillMinima", fillminima_fillMinima, METH_VARARGS, 
"function to implement filling of local minima in a raster surface:\n"
"call signature: fillminima(inarray, outarray, hMin, hMax, nullmask, boundaryval, boundaryRows, boundaryCols)\n"
"where:\n"
"   inarray is the input array\n"
"   outarray is the output array\n"
"   hMin is the minimum of the input image (excluding null values)\n"
"   hMax is the maximum of the input image (excluding null values)\n"
"   nullmask is the mask of where null values exist in inarray\n"
"   boundaryVal is the input boundar value\n"
"   boundaryRows and boundaryCols specify the boundary of the search\n"},
    {NULL}        /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3

static int fillminima_traverse(PyObject *m, visitproc visit, void *arg) 
{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int fillminima_clear(PyObject *m) 
{
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_fillminima",
        NULL,
        sizeof(struct FillMinimaState),
        FillMinimaMethods,
        NULL,
        fillminima_traverse,
        fillminima_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC 
PyInit__fillminima(void)

#else
#define INITERROR return

PyMODINIT_FUNC
init_fillminima(void)
#endif
{
PyObject *pModule;
struct FillMinimaState *state;

    // initialize the numpy stuff
    import_array();

#if PY_MAJOR_VERSION >= 3
    pModule = PyModule_Create(&moduledef);
#else
    pModule = Py_InitModule("_fillminima", FillMinimaMethods);
#endif
    if( pModule == NULL )
        INITERROR;

    state = GETSTATE(pModule);

    // Create and add our exception type
    state->error = PyErr_NewException("_fillminima.error", NULL, NULL);
    if( state->error == NULL )
    {
        Py_DECREF(pModule);
        INITERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return pModule;
#endif
}


