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

#include <Python.h>
#include "numpy/arrayobject.h"
#include <math.h>

/* An exception object for this module */
/* created in the init function */
struct ValueIndexesState
{
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct ValueIndexesState*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct ValueIndexesState _state;
#endif

/* 2^32 - 1 */
#define MAXUINT32 4294967295

static PyObject *valueIndexes_valndxFunc(PyObject *self, PyObject *args)
{
PyArrayObject *pInput, *pIndexes, *pValLU, *pCurrentIdx;
npy_int64 nMin, nMax, arrVal = 0;
int nDone = 0, nFound = 0, nDim, nType, i, idx;
npy_uint32 j, m;
npy_intp *pDims, *pCurrIdx;

    if( !PyArg_ParseTuple(args, "OOLLOO:valndxFunc", &pInput, &pIndexes, &nMin, &nMax, &pValLU, &pCurrentIdx))
        return NULL;

    if( !PyArray_Check(pInput) || !PyArray_Check(pIndexes) || !PyArray_Check(pValLU) || !PyArray_Check(pCurrentIdx) )
    {
        PyErr_SetString(GETSTATE(self)->error, "parameters 0, 1, 4 and 5 must be numpy arrays");
        return NULL;
    }

    if( (PyArray_TYPE(pIndexes) != NPY_UINT32) || (PyArray_TYPE(pValLU) != NPY_UINT32) || (PyArray_TYPE(pCurrentIdx) != NPY_UINT32))
    {
        PyErr_SetString(GETSTATE(self)->error, "parameters 1, 4 and 5 must be uint32 arrays");
        return NULL;
    }

    nDim = PyArray_NDIM(pInput);
    pDims = PyArray_DIMS(pInput);
    nType = PyArray_TYPE(pInput);

    if( !PyTypeNum_ISINTEGER(nType) )
    {
        PyErr_SetString(GETSTATE(self)->error, "parameter 0 must be an integer array");
        return NULL;
    }

    pCurrIdx = (npy_intp*)calloc(nDim, sizeof(npy_intp));

    while( !nDone )
    {
        /* Get the value at pCurrIdx */
        switch(nType)
        {
        case NPY_INT8: arrVal = (npy_uint64) *((npy_int8*)PyArray_GetPtr(pInput, pCurrIdx)); break;
        case NPY_UINT8: arrVal = (npy_uint64) *((npy_uint8*)PyArray_GetPtr(pInput, pCurrIdx)); break;
        case NPY_INT16: arrVal = (npy_uint64) *((npy_int16*)PyArray_GetPtr(pInput, pCurrIdx)); break;
        case NPY_UINT16: arrVal = (npy_uint64) *((npy_uint16*)PyArray_GetPtr(pInput, pCurrIdx)); break;
        case NPY_INT32: arrVal = (npy_uint64) *((npy_int32*)PyArray_GetPtr(pInput, pCurrIdx)); break;
        case NPY_UINT32: arrVal = (npy_uint64) *((npy_uint32*)PyArray_GetPtr(pInput, pCurrIdx)); break;
        case NPY_INT64: arrVal = (npy_uint64) *((npy_int64*)PyArray_GetPtr(pInput, pCurrIdx)); break;
        case NPY_UINT64: arrVal = (npy_uint64) *((npy_uint64*)PyArray_GetPtr(pInput, pCurrIdx)); break;
        }
        
        nFound = 0;
        j = 0;
        if( (arrVal >= nMin) && (arrVal <= nMax))
        {
            j = *((npy_uint32*)PyArray_GETPTR1(pValLU, arrVal - nMin));
            nFound = (j < MAXUINT32);
        }

        if( nFound )
        {
            m = *((npy_uint32*)PyArray_GETPTR1(pCurrentIdx, j));
            for( i = 0; i < nDim; i++ )
            {
                *((npy_uint32*)PyArray_GETPTR2(pIndexes, m, i)) = pCurrIdx[i];
            }
            *((npy_uint32*)PyArray_GETPTR1(pCurrentIdx, j)) = m + 1;
        }

        /* code that updates curridx - incs the next dim
         if we have done all the elements in the current 
         dim */
        idx = nDim - 1;
        while( idx >= 0 )
        {
            pCurrIdx[idx]++;
            if( pCurrIdx[idx] >= pDims[idx] )
            {
                pCurrIdx[idx] = 0;
                idx--;
            }
            else
            {
                break;
            }
        }

        /* if we are done we have run out of dims */
        nDone = (idx < 0);
    }

    free(pCurrIdx);

    Py_RETURN_NONE;
}


// Our list of functions in this module
static PyMethodDef ValueIndexesMethods[] = {
    {"valndxFunc", valueIndexes_valndxFunc, METH_VARARGS, 
"function to go through an array and create a lookup table of the array indexes for each distinct value in the data array:\n"
"call signature: valndxFunc(input, indexes, min, max, ValLU, CurrentIdx)\n"
"where:\n"
"   input is the input array\n"
"   indexes is the output array which will be filled with the indexes of each value\n"
"   min is the minimum value of input\n"
"   max is the minimum value of input\n"
"   ValLU is the lookup array to index into indexes and CurrentIdx\n"
"   CurrentIdx has the current index\n"
"   \n"},
    {NULL}        /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3

static int valueIndexes_traverse(PyObject *m, visitproc visit, void *arg) 
{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int valueIndexes_clear(PyObject *m) 
{
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_valueindexes",
        NULL,
        sizeof(struct ValueIndexesState),
        ValueIndexesMethods,
        NULL,
        valueIndexes_traverse,
        valueIndexes_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC 
PyInit__valueindexes(void)

#else
#define INITERROR return

PyMODINIT_FUNC
init_valueindexes(void)
#endif
{
PyObject *pModule;
struct ValueIndexesState *state;

    /* initialize the numpy stuff */
    import_array();

#if PY_MAJOR_VERSION >= 3
    pModule = PyModule_Create(&moduledef);
#else
    pModule = Py_InitModule("_valueindexes", ValueIndexesMethods);
#endif
    if( pModule == NULL )
        INITERROR;

    state = GETSTATE(pModule);

    /* Create and add our exception type */
    state->error = PyErr_NewException("_valueindexes.error", NULL, NULL);
    if( state->error == NULL )
    {
        Py_DECREF(pModule);
        INITERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return pModule;
#endif
}


