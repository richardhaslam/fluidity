/*  Copyright (C) 2006 Imperial College London and others.

Please see the AUTHORS file in the main source directory for a full list
of copyright holders.

Prof. C Pain
Applied Modelling and Computation Group
Department of Earth Science and Engineering
Imperial College London

amcgsoftware@imperial.ac.uk

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation,
version 2.1 of the License.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
USA
*/



#include "confdefs.h"
#include "string.h"

#ifdef HAVE_PYTHON
#include "Python.h"
#if PY_MAJOR_VERSION >= 3
#define PyInt_FromLong PyLong_FromLong
#define PyInt_AsLong PyLong_AsLong
#define PyString_Size PyUnicode_GET_SIZE
#define PyString_AsString PyUnicode_AsUTF8
#endif
#endif
#ifdef HAVE_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#endif

void set_field_from_python(char *function, int function_len, int dim, int nodes,
			   double *x, double *y, double *z, double t, double *dt,
			   int *stat, int *result_dim, void **result,
			   int (*set_result)(int, int *, void **, PyObject *))
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t)function_len);
  for (i = 0; i < function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }

  *stat=1;
  return;
#else

  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pPos, *px, *pT, *pDT;
  char *function_c;
  int i, res;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(function_len+1);
  memcpy(function_c, function, function_len);
  function_c[function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");
  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals = PyDict_New();

  // Execute the user's code and clean up allocated string.
  pCode = PyRun_String(function_c, Py_file_input, pGlobals, pLocals);
  free(function_c);

  // Extract the function from the code.
  pFunc = PyDict_GetItemString(pLocals, "val");
  if (pFunc == NULL) {
      printf("Couldn't find a 'val' function in your Python code.\n");
      *stat = 1;
      return;
  }

  // Check for errors in executing user code.
  if (PyErr_Occurred()) {
    PyErr_Print();
    *stat = 1;
    return;
  }

  // Create Python objects for function arguments
  pT = PyFloat_FromDouble(t);
  pPos = PyTuple_New(dim);
  pArgs = PyTuple_New(2);
  // note that PyTuple_SetItem steals refs, so we don't have to decrement
  // the refcounts of pT or pPos
  PyTuple_SetItem(pArgs, 0, pPos);
  PyTuple_SetItem(pArgs, 1, pT);

  // particle routines need dt too
  if (dt != NULL) {
    pDT = PyFloat_FromDouble(*dt);
    PyTuple_SetItem(pArgs, 2, pDT);
  }

  // Check for a Python error in the function call
  if (PyErr_Occurred()) {
    PyErr_Print();
    *stat = 1;
    return;
  }

  // populate position tuple
  for (i = 0; i < nodes; i++) {
    px = PyFloat_FromDouble(x[i]);
    PyTuple_SetItem(pPos, 0, px);

    if (dim > 1) {
      px = PyFloat_FromDouble(y[i]);
      PyTuple_SetItem(pPos, 1, px);

      if (dim > 2) {
        px = PyFloat_FromDouble(z[i]);
        PyTuple_SetItem(pPos, 2, px);
      }
    }

    pResult = PyObject_CallObject(pFunc, pArgs);

    // Check for a Python error in the function call
    if (PyErr_Occurred()) {
      PyErr_Print();
      *stat = 1;
      return;
    }

    res = set_result(i, result_dim, result, pResult);
    if (res) {
      *stat = res;
      return;
    }
    // Check for a Python error in result.
    if (PyErr_Occurred()) {
      PyErr_Print();
      *stat = 1;
      return;
    }

    Py_DECREF(pResult);
  }

  // clean up
  Py_DECREF(pCode);
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  PyGC_Collect();
  *stat = 0;
  #endif
}

int set_scalar_result_double(int i, int *dim, double **result, PyObject *pResult)
{
  *result[i] = PyFloat_AsDouble(pResult);
  return 0;
}

#define set_scalar_field_from_python F77_FUNC(set_scalar_field_from_python, SET_SCALAR_FIELD_FROM_PYTHON)
void set_scalar_field_from_python(char *function, int function_len, int dim,
                                  int nodes, double *x, double *y, double *z, double t,
				  double *result, int *stat)
{
  set_field_from_python(function, function_len, dim, nodes, x, y, z, t, NULL, stat,
			NULL, &result, set_scalar_result_double);
}

#define set_scalar_particles_from_python F77_FUNC(set_scalar_particles_from_python, SET_SCALAR_PARTICLES_FROM_PYTHON)
void set_scalar_particles_from_python(char *function, int function_len, int dim, int ndete,
				      double *x, double *y, double *z, double t, double dt,
				      double *result, int *stat)
{
  set_field_from_python(function, function_len, dim, ndete, x, y, z, t, &dt, stat,
			NULL, &result, set_scalar_result_double);
}

int set_scalar_result_integer(int i, int *dim, int **result, PyObject *pResult)
{
  *result[i] = PyLong_AsLong(pResult);
  return 0;
}

#define set_integer_array_from_python F77_FUNC(set_integer_array_from_python, SET_INTEGER_ARRAY_FROM_PYTHON)
void set_integer_array_from_python(char* function, int function_len, int dim,
                                   int nodes, double *x, double *y, double *z, double t,
                                   int* result, int* stat)
{
  set_field_from_python(function, function_len, dim, nodes, x, y, z, t, NULL, stat,
			NULL, &result, set_scalar_result_integer);
}

int set_vector_result_double(int i, int *dim, double **result, PyObject *pResult)
{
  if (PyObject_Length(pResult) != *dim) {
      fprintf(stderr, "Error: length of object returned from python (%d) does not match the allocated dimension of the vector field (%d).\n",
              PyObject_Length(pResult), *dim);
      return 1;
  }

  PyObject *px;
  for (int d = 0; d < *dim; d++) {
    px = PySequence_GetItem(pResult, d);
    result[d][i] = PyFloat_AsDouble(px);
    Py_DECREF(px);
  }

  return 0;
}

#define set_vector_field_from_python F77_FUNC(set_vector_field_from_python, SET_VECTOR_FIELD_FROM_PYTHON)
void set_vector_field_from_python(char *function, int function_len, int dim,
                                  int nodes, double *x, double *y, double *z, double t,
                                  int result_dim, double *result_x, double *result_y, double *result_z,
				  int *stat)
{
  double *results[] = {result_x, result_y, result_z};

  set_field_from_python(function, function_len, dim, nodes, x, y, z, t, NULL, stat,
			&result_dim, &results, set_vector_result_double);
}

int set_tensor_result_double(int i, int *dim, double **result, PyObject *pResult)
{
  PyArrayObject *pArray;

  // get a 2d array from result
  pArray = (PyArrayObject *)PyArray_ContiguousFromObject(pResult, NPY_DOUBLE, 2, 2);
  if (PyErr_Occurred()) {
    PyErr_Print();
    return 1;
  }
  if (PyArray_DIMS(pArray)[0] != dim[0] || PyArray_DIMS(pArray)[1] != dim[1]) {
    fprintf(stderr, "Error: dimensions of array returned from python ([%d, %d]) do not match allocated dimensions of the tensor_field ([%d, %d])).\n",
	    (int) PyArray_DIMS(pArray)[0], (int) PyArray_DIMS(pArray)[1], dim[0], dim[1]);
    return 1;
  }

  for (int ii = 0; ii < dim[0]; ii++) {
    for (int jj = 0; jj < dim[1]; jj++) {
      *result[i*(dim[0]+dim[1]) + jj*dim[0] + ii] = *((double*)PyArray_GETPTR2(pArray, ii, jj));
    }
  }
  Py_DECREF(pArray);

}

#define set_tensor_field_from_python F77_FUNC(set_tensor_field_from_python, SET_TENSOR_FIELD_FROM_PYTHON)
void set_tensor_field_from_python(char *function, int function_len, int dim,
                                  int nodes, double *x, double *y, double *z, double t,
				  int *result_dim, double *result, int *stat)
{
#ifndef HAVE_NUMPY
  int i;
  strncpy(function, "No Numpy support!\n", (size_t)function_len);
  for (i=0; i < function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat = 1;
  return;
#else
  import_array();

  set_field_from_python(function, function_len, dim, nodes, x, y, z, t, NULL, stat,
			result_dim, &result, set_tensor_result_double);

#endif
}


#define set_detectors_from_python F77_FUNC(set_detectors_from_python, SET_DETECTORS_FROM_PYTHON)
void set_detectors_from_python(char *function, int *function_len, int *dim,
                               int *ndete, double *t,
                               int *result_dim,
                               double result_x[], double result_y[],
                               double result_z[], int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult, *pResultItem,
    *pArgs, *px, *pT;
  char *function_c;
  int i;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");
  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
   if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }


  for (i = 0; i < *ndete; i++){
    pResultItem = PySequence_GetItem(pResult, i);

    px=PySequence_GetItem(pResultItem, 0);

    result_x[i]=PyFloat_AsDouble(px);
    // Check for a Python error in unpacking tuple.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }
    Py_DECREF(px);

    if (*result_dim>1) {
      px=PySequence_GetItem(pResultItem, 1);
      result_y[i]=PyFloat_AsDouble(px);
      // Check for a Python error in unpacking tuple.
      if (PyErr_Occurred()){
         PyErr_Print();
         return;
      }

      Py_DECREF(px);
      if (*result_dim>2) {
        px=PySequence_GetItem(pResultItem, 2);
        result_z[i]=PyFloat_AsDouble(px);
      // Check for a Python error in unpacking tuple.
       if (PyErr_Occurred()){
          PyErr_Print();
          return;
       }
        Py_DECREF(px);
      }
    }

    Py_DECREF(pResultItem);
  }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();

  *stat=0;
  return;
#endif
}

//#define set_scalar_particles_from_python_fields F77_FUNC(set_scalar_particles_from_python_fields, SET_SCALAR_PARTICLES_FROM_PYTHON_FIELDS)
void set_scalar_particles_from_python_fields(char *function, int function_len, int dim, int ndete,
				      double x[], double y[], double z[], double t, double dt,
				      int FIELD_NAME_LEN, int nfields[],
				      char field_names[nfields[0]+nfields[1]+nfields[2]][FIELD_NAME_LEN],
				      double field_vals[ndete][nfields[0]+dim*nfields[1]+dim*dim*nfields[2]],
				      int old_nfields[],
				      char old_field_names[old_nfields[0]+old_nfields[1]+old_nfields[2]][FIELD_NAME_LEN],
				      double old_field_vals[ndete][old_nfields[0]+dim*old_nfields[1]+dim*dim*old_nfields[2]],
				      int old_nattributes[],
				      char old_att_names[old_nattributes[0]+old_nattributes[1]+old_nattributes[2]][FIELD_NAME_LEN],
				      double old_attributes[ndete][old_nattributes[0]+dim*old_nattributes[1]+dim*dim*old_nattributes[2]],
				      double result[], int* stat)

{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }

  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pPos, *px, *pT, *pdT, *pField, *pNames;
  PyObject *vec, *tens;
  char *function_c;
  int i, j, k, l;
  npy_intp dims[] = {dim, dim};
  import_array();

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(function_len+1);
  memcpy(function_c, function, function_len);
  function_c[function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals = PyDict_New();

  // Execute the user's code.
  pCode = PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc = PyDict_GetItemString(pLocals, "val");
  if (pFunc == NULL) {
      printf("Couldn't find a 'val' function in your Python code.\n");
      *stat=1;
      return;
  }

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // allocate Python object for function arguments
  // Field variable dictionary space
  pNames = PyDict_New();
  // Python form of time variable.
  pT = PyFloat_FromDouble(t);
  // Python form of timestep variable.
  pdT = PyFloat_FromDouble(dt);
  // Tuple containing the current position vector.
  pPos = PyTuple_New(dim);

  // Tuple of arguments to function;
  pArgs = PyTuple_New(4);
  PyTuple_SetItem(pArgs, 3, pNames);
  PyTuple_SetItem(pArgs, 2, pdT);
  PyTuple_SetItem(pArgs, 1, pT);
  PyTuple_SetItem(pArgs, 0, pPos);

  //loop over all particles
  for (i = 0; i < ndete; i++)
    {
      // Set values for position vector.
      px = PyFloat_FromDouble(x[i]);
      PyTuple_SetItem(pPos, 0, px);

      if (dim>1) {
	px = PyFloat_FromDouble(y[i]);
	PyTuple_SetItem(pPos, 1, px);

	if (dim>2) {
	  px = PyFloat_FromDouble(z[i]);
	  PyTuple_SetItem(pPos, 2, px);
	}
      }

      //Set values for fields dictionary
      //Set field values (scalar, then vector, then tensor)
      for (j=0; j < nfields[0]; j++)//scalar fields
      	{
      	  pField = PyFloat_FromDouble(field_vals[i][j]);
      	  PyDict_SetItemString(pNames, field_names[j], pField);
	  Py_DECREF(pField);
      	}

      for (j=0; j < nfields[1]; j++)//vector fields
      	{
	  vec = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
      	  for (k=0; k < dim; k++)
      	    {
      	      ((double*)PyArray_DATA(vec))[k] = field_vals[i][nfields[0] + (j*dim) + k];
      	    }
      	  PyDict_SetItemString(pNames, field_names[j + nfields[0]], vec);
	  Py_DECREF(vec);
      	}

      for (j=0; j < nfields[2]; j++)//tensor fields
      	{
	  tens = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
      	  for (k=0; k < dim; k++)
      	    {
      	      for (l=0; l < dim; l++)
      		{
		  //l and k swapped in ptField so tensor is correctly read into python
      		  ((double*)PyArray_DATA(tens))[k + l*dim] = field_vals[i][nfields[0] + nfields[1]*dim + (j*dim*dim) + (k*dim) + l];
      		}
      	    }
	   PyDict_SetItemString(pNames, field_names[j + nfields[0] + nfields[1]], tens);
	   Py_DECREF(tens);
      	}


      //Set old_field values (scalar, then vector, then tensor)
      for (j=0; j < old_nfields[0]; j++)//scalar old_fields
      	{
      	  pField = PyFloat_FromDouble(old_field_vals[i][j]);
      	  PyDict_SetItemString(pNames, old_field_names[j], pField);
	  Py_DECREF(pField);
      	}

      for (j=0; j < old_nfields[1]; j++)//vector old_fields
      	{
	  vec = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
      	  for (k=0; k < dim; k++)
      	    {
      	      ((double*)PyArray_DATA(vec))[k] = old_field_vals[i][old_nfields[0] + (j*dim) + k];
      	    }
	  PyDict_SetItemString(pNames, old_field_names[j+old_nfields[0]], vec);
	  Py_DECREF(vec);
      	}

      for (j=0; j < old_nfields[2]; j++)//tensor old_fields
      	{
	  tens = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
      	  for (k=0; k < dim; k++)
      	    {
      	      for (l=0; l < dim; l++)
      		{
		  //l and k swapped in ptField so tensor is correctly read into python
      		  ((double*)PyArray_DATA(tens))[k + l*dim] = old_field_vals[i][old_nfields[0] + old_nfields[1]*dim + (j*dim*dim) + (k*dim) + l];
      		}
      	    }
	  PyDict_SetItemString(pNames, old_field_names[j + old_nfields[0] + old_nfields[1]], tens);
	  Py_DECREF(tens);
      	}

      //Set old_attribute values (scalar, then vector, then tensor)
      for (j=0; j < old_nattributes[0]; j++)//scalar old_attributes
      	{
      	  pField = PyFloat_FromDouble(old_attributes[i][j]);
      	  PyDict_SetItemString(pNames, old_att_names[j], pField);
	  Py_DECREF(pField);
      	}

      for (j=0; j < old_nattributes[1]; j++)//vector old_attributes
      	{
	  vec = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
      	  for (k=0; k < dim; k++)
      	    {
      	      ((double*)PyArray_DATA(vec))[k] = old_attributes[i][old_nattributes[0] + (j*dim) + k];
      	    }

	  PyDict_SetItemString(pNames, old_att_names[j+old_nattributes[0]], vec);
	  Py_DECREF(vec);
      	}

      for (j=0; j < old_nattributes[2]; j++)//tensor old_attributes
      	{
	  tens = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
      	  for (k=0; k < dim; k++)
      	    {
      	      for (l=0; l < dim; l++)
      		{
		  //l and k swapped in ptField so tensor is correctly read into python
      		  ((double*)PyArray_DATA(tens))[k + l*dim] = old_attributes[i][old_nattributes[0] + old_nattributes[1]*dim + (j*dim*dim) + (k*dim) + l];
      		}
      	    }
	  PyDict_SetItemString(pNames, old_att_names[j + old_nattributes[0] + old_nattributes[1]], tens);
	  Py_DECREF(tens);
      	}

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      pResult=PyObject_CallObject(pFunc, pArgs);

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      result[i]=PyFloat_AsDouble(pResult);

      // Check for a Python error in result.
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

    }
  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);
  Py_DECREF(pNames);

  // Force a garbage collection
  PyGC_Collect();

  *stat=0;
  return;
#endif
}

#define set_vector_particles_from_python F77_FUNC(set_vector_particles_from_python, SET_VECTOR_PARTICLES_FROM_PYTHON)
void set_vector_particles_from_python(char *function, int *function_len, int *dim, int *ndete,
			       double x[], double y[], double z[], double *t, double *dt,
				      double result[*dim][*ndete], int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }

  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pPos, *px, *pT, *pdT;

  char *function_c;
  int i;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");
  if (pFunc == NULL) {
      printf("Couldn't find a 'val' function in your Python code.\n");
      *stat=1;
      return;
  }

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Python form of timestep variable.
  pdT=PyFloat_FromDouble(*dt);

  // Tuple containing the current position vector.
  pPos=PyTuple_New(*dim);

  //Tuple containing the Arguments
  pArgs=PyTuple_New(3);
  PyTuple_SetItem(pArgs, 2, pdT);
  PyTuple_SetItem(pArgs, 1, pT);
  PyTuple_SetItem(pArgs, 0, pPos);

  for (i = 0; i < *ndete; i++)
    {

      // Set values for position vector.

      px=PyFloat_FromDouble(x[i]);
      PyTuple_SetItem(pPos, 0, px);

      if (*dim>1) {
	px=PyFloat_FromDouble(y[i]);
	PyTuple_SetItem(pPos, 1, px);

	if (*dim>2) {
	  px=PyFloat_FromDouble(z[i]);
	  PyTuple_SetItem(pPos, 2, px);
	}
      }

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      pResult=PyObject_CallObject(pFunc, pArgs);

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      if (PyObject_Length(pResult) != *dim)
      	{
      	  fprintf(stderr, "Error: length of object returned from python (%d) does not match the allocated dimension of the vector field (%d).\n",
      		  (int) PyObject_Length(pResult), *dim);
      	  *stat = 1;
      	  return;
      	}

      px=PySequence_GetItem(pResult, 0);
      result[0][i]=PyFloat_AsDouble(px);
      // Check for a Python error in unpacking tuple.
      if (PyErr_Occurred()){
      	PyErr_Print();
      	*stat=1;
      	return;
      }

      if (*dim>1) {
      	px=PySequence_GetItem(pResult, 1);
      	result[1][i]=PyFloat_AsDouble(px);
      	// Check for a Python error in unpacking tuple.
      	if (PyErr_Occurred()){
      	  PyErr_Print();
      	  return;
      	}

      	if (*dim>2) {
      	  px=PySequence_GetItem(pResult, 2);
      	  result[2][i]=PyFloat_AsDouble(px);
      	  // Check for a Python error in unpacking tuple.
      	  if (PyErr_Occurred()){
      	    PyErr_Print();
      	    return;
      	  }
      	}
      }

    }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);
  Py_DECREF(px);

  // Force a garbage collection
  PyGC_Collect();

  *stat=0;
  return;
#endif
}

//#define set_vector_particles_from_python_fields F77_FUNC(set_vector_particles_from_python_fields, SET_VECTOR_PARTICLES_FROM_PYTHON_FIELDS)
void set_vector_particles_from_python_fields(char *function, int *function_len, int *dim, int *ndete,
				      double x[], double y[], double z[], double *t, double *dt,
				      int *FIELD_NAME_LEN, int nfields[], char field_names[nfields[0]+nfields[1]+nfields[2]][*FIELD_NAME_LEN],
				      double field_vals[], int old_nfields[], char old_field_names[old_nfields[0]+old_nfields[1]+old_nfields[2]][*FIELD_NAME_LEN],
				      double old_field_vals[], int old_nattributes[],
				      char old_att_names[old_nattributes[0]+old_nattributes[1]+old_nattributes[2]][*FIELD_NAME_LEN],
				      double old_attributes[], double result[], int* stat)

{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }

  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pPos, *px, *pT, *pdT, *pField, *pNames;
  PyObject *vec, *tens;
  double *pvField = malloc(sizeof(double[*dim]));
  double *ptField = malloc(sizeof(double[*dim * *dim]));
  char *function_c;
  int i, j, k, l;
  int *dims[] = {*dim, *dim};
  import_array();

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");
  if (pFunc == NULL) {
      printf("Couldn't find a 'val' function in your Python code.\n");
      *stat=1;
      return;
  }

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Field variable dictionary space
  pNames=PyDict_New();

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Python form of timestep variable.
  pdT=PyFloat_FromDouble(*dt);

  // Tuple containing the current position vector.
  pPos=PyTuple_New(*dim);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(4);
  PyTuple_SetItem(pArgs, 3, pNames);
  PyTuple_SetItem(pArgs, 2, pdT);
  PyTuple_SetItem(pArgs, 1, pT);
  PyTuple_SetItem(pArgs, 0, pPos);

  //loop over all particles

  for (i = 0; i < *ndete; i++)
    {

      // Set values for position vector.
      px=PyFloat_FromDouble(x[i]);
      PyTuple_SetItem(pPos, 0, px);

      if (*dim>1) {
	px=PyFloat_FromDouble(y[i]);
	PyTuple_SetItem(pPos, 1, px);

	if (*dim>2) {
	  px=PyFloat_FromDouble(z[i]);
	  PyTuple_SetItem(pPos, 2, px);
	}
      }

      //Set values for fields library.

      //Set field values (scalar, then vector, then tensor)
      for (j=0; j < nfields[0]; j++)//scalar fields
      	{
      	  pField=PyFloat_FromDouble(field_vals[i * (nfields[0]+(*dim * nfields[1])+(*dim * *dim * nfields[2])) + j]);
      	  PyDict_SetItemString(pNames, field_names[j], pField);
      	}

      for (j=0; j < nfields[1]; j++)//vector fields
      	{
      	  for (k=0; k < *dim; k++)
      	    {
      	      pvField[k]=field_vals[i * (nfields[0]+(*dim * nfields[1])+(*dim * *dim * nfields[2])) + nfields[0] + (j * *dim) + k];
	      printf("%s\n", field_names[j+nfields[0]]);
	      printf("%f\n",pvField[k]);
      	    }
      	  vec=PyArray_SimpleNewFromData(1,dim,NPY_DOUBLE,pvField);
      	  PyDict_SetItemString(pNames, field_names[j+nfields[0]], vec);
      	}

      for (j=0; j < nfields[2]; j++)//tensor fields
      	{
      	  for (k=0; k < *dim; k++)
      	    {
      	      for (l=0; l < *dim; l++)
      		{
		  //l and k swapped in ptField so tensor is correctly read into python
      		  ptField[k+(l * *dim)]=field_vals[i * (nfields[0]+(*dim * nfields[1])+(*dim * *dim * nfields[2])) + nfields[0] + (nfields[1] * *dim) + (j * *dim * *dim) + (k * *dim) + l];
      		}
      	    }
	   tens = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ptField);
	   PyDict_SetItemString(pNames, field_names[j + nfields[0] + nfields[1]], tens);
      	}


      //Set old_field values (scalar, then vector, then tensor)
      for (j=0; j < old_nfields[0]; j++)//scalar old_fields
      	{
      	  pField=PyFloat_FromDouble(old_field_vals[i * (old_nfields[0]+(*dim * old_nfields[1])+(*dim * *dim * old_nfields[2])) + j]);
      	  PyDict_SetItemString(pNames, old_field_names[j], pField);
      	}

      for (j=0; j < old_nfields[1]; j++)//vector old_fields
      	{
      	  for (k=0; k < *dim; k++)
      	    {
      	      pvField[k]=old_field_vals[i * (old_nfields[0]+(*dim * old_nfields[1])+(*dim * *dim * old_nfields[2])) + old_nfields[0] + (j * *dim) + k];
      	    }
	  vec=PyArray_SimpleNewFromData(1,dim,NPY_DOUBLE,pvField);
	  PyDict_SetItemString(pNames, old_field_names[j+old_nfields[0]], vec);
      	}

      for (j=0; j < old_nfields[2]; j++)//tensor old_fields
      	{
      	  for (k=0; k < *dim; k++)
      	    {
      	      for (l=0; l < *dim; l++)
      		{
		  //l and k swapped in ptField so tensor is correctly read into python
      		  ptField[k+(l * *dim)]=old_field_vals[i * (old_nfields[0]+(*dim * old_nfields[1])+(*dim * *dim * old_nfields[2])) + old_nfields[0] + (old_nfields[1] * *dim) + (j * *dim * *dim) + (k * *dim) + l];
      		}
      	    }
	  tens = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ptField);
	  PyDict_SetItemString(pNames, old_field_names[j + old_nfields[0] + old_nfields[1]], tens);
      	}

      //Set old_attribute values (scalar, then vector, then tensor)
      for (j=0; j < old_nattributes[0]; j++)//scalar old_attributes
      	{
      	  pField=PyFloat_FromDouble(old_attributes[i * (old_nattributes[0]+(*dim * old_nattributes[1])+(*dim * *dim * old_nattributes[2])) + j]);
      	  PyDict_SetItemString(pNames, old_att_names[j], pField);
      	}

      for (j=0; j < old_nattributes[1]; j++)//vector old_attributes
      	{
      	  for (k=0; k < *dim; k++)
      	    {
      	      pvField[k]=old_attributes[i * (old_nattributes[0]+(*dim * old_nattributes[1])+(*dim * *dim * old_nattributes[2])) + old_nattributes[0] + (j * *dim) + k];
	      printf("%s\n", old_att_names[j+old_nattributes[0]]);
	      printf("%f\n",pvField[k]);
      	    }

	  vec=PyArray_SimpleNewFromData(1,dim,NPY_DOUBLE,pvField);
	  PyDict_SetItemString(pNames, old_att_names[j+old_nattributes[0]], vec);
      	}

      for (j=0; j < old_nattributes[2]; j++)//tensor old_attributes
      	{
      	  for (k=0; k < *dim; k++)
      	    {
      	      for (l=0; l < *dim; l++)
      		{
		  //l and k swapped in ptField so tensor is correctly read into python
      		  ptField[k+(l * *dim)]=old_attributes[i * (old_nattributes[0]+(*dim * old_nattributes[1])+(*dim * *dim * old_nattributes[2])) + old_nattributes[0] + (old_nattributes[1] * *dim) + (j * *dim * *dim) + (k * *dim) + l];
      		}
      	    }
	  tens = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ptField);
	  PyDict_SetItemString(pNames, old_att_names[j + old_nattributes[0] + old_nattributes[1]], tens);
      	}

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      printf("Results for dictionary\n");
      PyObject *key, *value, *py;
      Py_ssize_t pos = 0;

      while (PyDict_Next(pNames,&pos, &key, &value))
	{
	  printf("%s\n", key);
	  /* if (PyObject_Size(value) == 1) */
	  /*   { */
	  /*     printf("scalar value\n"); */
	  /*     printf("%f\n", value); */
	  /*   } */
	  /* if (PyObject_Size(value) == *dim) */
	  /*   { */
	  /*     printf("vector value\n"); */
	  /*     py=PySequence_GetItem(value, 0); */
	  /*     printf("%f\n", py); */
	  /*     py=PySequence_GetItem(value, 1); */
	  /*     printf("%f\n", py); */
	  /*   } */
	  /* if (PyObject_Size(value) == *dim * *dim) */
	  /*   { */
	  /*     printf("tensor value\n"); */
	  /*     py=PySequence_GetItem(value, 0); */
	  /*     printf("%f\n", py); */
	  /*     py=PySequence_GetItem(value, 1); */
	  /*     printf("%f\n", py); */
	  /*     py=PySequence_GetItem(value, 2); */
	  /*     printf("%f\n", py); */
	  /*     py=PySequence_GetItem(value, 3); */
	  /*     printf("%f\n", py); */
	  /*   } */

	  printf("%lf\n", value);
	  printf("%d\n", pos);
	}

      pResult=PyObject_CallObject(pFunc, pArgs);

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      if (PyObject_Length(pResult) != *dim)
	{
	  fprintf(stderr, "Error: length of object returned from python (%d) does not match the allocated dimension of the vector field (%d).\n",
		  (int) PyObject_Length(pResult), *dim);
	  *stat = 1;
	  return;
	}

      px=PySequence_GetItem(pResult, 0);
      result[(i * *dim)]=PyFloat_AsDouble(px);
      // Check for a Python error in unpacking tuple.
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }
      Py_DECREF(px);

      if (*dim>1) {
	px=PySequence_GetItem(pResult, 1);
	result[(i * *dim) + 1]=PyFloat_AsDouble(px);
	// Check for a Python error in unpacking tuple.
	if (PyErr_Occurred()){
	  PyErr_Print();
	  return;
	}
	Py_DECREF(px);

	if (*dim>2) {
	  px=PySequence_GetItem(pResult, 2);
	  result[(i * *dim) + 2]=PyFloat_AsDouble(px);
	  // Check for a Python error in unpacking tuple.
	  if (PyErr_Occurred()){
	    PyErr_Print();
	    return;
	  }
	  Py_DECREF(px);
	}
      }

    }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  //Free allocated memory
  free(ptField);
  free(pvField);
  Py_DECREF(tens);
  Py_DECREF(vec);

  // Force a garbage collection
  PyGC_Collect();

  *stat=0;
  return;
#endif
}


#define set_tensor_particles_from_python F77_FUNC(set_tensor_particles_from_python, SET_TENSOR_PARTICLES_FROM_PYTHON)
void set_tensor_particles_from_python(char *function, int *function_len, int *dim, int *ndete,
			       double x[], double y[], double z[], double *t, double *dt,
                                  double result[*dim][*dim][*ndete], int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }

  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pPos, *px, *pT, *pdT;
  PyArrayObject *pArray;
  char *function_c;
  int i, j, k;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");
  if (pFunc == NULL) {
      printf("Couldn't find a 'val' function in your Python code.\n");
      *stat=1;
      return;
  }

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Python form of timestep variable.
  pdT=PyFloat_FromDouble(*dt);

  // Tuple containing the current position vector.
  pPos=PyTuple_New(*dim);

  //Tuple containing the Arguments
  pArgs=PyTuple_New(3);
  PyTuple_SetItem(pArgs, 2, pdT);
  PyTuple_SetItem(pArgs, 1, pT);
  PyTuple_SetItem(pArgs, 0, pPos);

  for (i = 0; i < *ndete; i++)
    {

      // Set values for position vector.

      px=PyFloat_FromDouble(x[i]);
      PyTuple_SetItem(pPos, 0, px);

      if (*dim>1) {
	px=PyFloat_FromDouble(y[i]);
	PyTuple_SetItem(pPos, 1, px);

	if (*dim>2) {
	  px=PyFloat_FromDouble(z[i]);
	  PyTuple_SetItem(pPos, 2, px);
	}
      }

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      pResult=PyObject_CallObject(pFunc, pArgs);

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      pArray = (PyArrayObject *)
	PyArray_ContiguousFromObject(pResult, NPY_DOUBLE, 2, 2);

      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      if (PyArray_DIMS(pArray)[0] != dim || PyArray_DIMS(pArray)[1] != dim)
	{
	  fprintf(stderr, "Error: dimensions of array returned from python ([%d, %d]) do not match allocated dimensions of the tensor_field ([%d, %d])).\n",
		  (int) PyArray_DIMS(pArray)[0], (int) PyArray_DIMS(pArray)[1], dim, dim);
	  *stat=1;
	  return;
	}

      for (j = 0; j < dim; j++){
	for (k = 0; k < dim; k++){

	  // Note the transpose for fortran.
	  double tmp;
	  tmp = *(double*)(PyArray_DATA(pArray) + j * PyArray_STRIDES(pArray)[0] + k * PyArray_STRIDES(pArray)[1]);
	  result[j][k][i] = tmp;
	}
      }

      Py_DECREF(pArray);

      // Check for a Python error in result.
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

    }
  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();

  *stat=0;
  return;
#endif
}

//#define set_tensor_particles_from_python_fields F77_FUNC(set_tensor_particles_from_python_fields, SET_TENSOR_PARTICLES_FROM_PYTHON_FIELDS)
void set_tensor_particles_from_python_fields(char *function, int *function_len, int *dim, int *ndete,
				      double x[], double y[], double z[], double *t, double *dt,
				      int *FIELD_NAME_LEN, int nfields[], char field_names[*nfields][*FIELD_NAME_LEN],
				      double *field_vals, int old_nfields[], char old_field_names[*old_nfields][*FIELD_NAME_LEN],
				      double *old_field_vals, int old_nattributes[],
				      char old_att_names[*old_nattributes][*FIELD_NAME_LEN],
				      double *old_attributes, double result[*dim][*dim][*ndete], int* stat)

{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }

  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pPos, *px, *pT, *pdT, *pField, *pNames;
  PyArrayObject *pArray;
  double (*fields_new)[*nfields] = malloc(sizeof(double[*ndete][nfields[0]+(*dim * nfields[1])+(*dim * *dim * nfields[2])])); //now have nfields(nscalar,nvector,ntensor) set up with this
  double (*fields_old)[*old_nfields] = malloc(sizeof(double[*ndete][old_nfields[0] + (*dim * old_nfields[1]) + (*dim * *dim * old_nfields[2])]));
  double (*attributes_old)[*old_nattributes] = malloc(sizeof(double[*ndete][old_nattributes[0] + (*dim * old_nattributes[1]) + (*dim * *dim * old_nattributes[2])]));
  char *function_c;
  int i, j, k, l;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");
  if (pFunc == NULL) {
      printf("Couldn't find a 'val' function in your Python code.\n");
      *stat=1;
      return;

  }

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  for(i = 0; i < *ndete; i++) {
    for(j = 0; j < *nfields; j++) {
      fields_new[i][j] = field_vals[i * (*nfields) + j];
    }
    for(j = 0; j < *old_nfields; j++) {
      fields_old[i][j] = old_field_vals[i * (*old_nfields) + j];
    }
    for(j = 0; j < *old_nattributes; j++) {
      attributes_old[i][j] = old_attributes[i * (*old_nattributes) + j];
    }
  }

  // Field variable dictionary space
  pNames=PyDict_New();

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Python form of timestep variable.
  pdT=PyFloat_FromDouble(*dt);

  // Tuple containing the current position vector.
  pPos=PyTuple_New(*dim);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(4);
  PyTuple_SetItem(pArgs, 3, pNames);
  PyTuple_SetItem(pArgs, 2, pdT);
  PyTuple_SetItem(pArgs, 1, pT);
  PyTuple_SetItem(pArgs, 0, pPos);

  //Set values for fields vector.

  for (i = 0; i < *ndete; i++)
    {

      // Set values for position vector.
      px=PyFloat_FromDouble(x[i]);
      PyTuple_SetItem(pPos, 0, px);

      if (*dim>1) {
	px=PyFloat_FromDouble(y[i]);
	PyTuple_SetItem(pPos, 1, px);

	if (*dim>2) {
	  px=PyFloat_FromDouble(z[i]);
	  PyTuple_SetItem(pPos, 2, px);
	}
      }

      //Set values for fields library.

      for (j=0; j < *nfields; j++)
      {
        pField=PyFloat_FromDouble(fields_new[i][j]);
        PyDict_SetItemString(pNames, field_names[j], pField);
      }

      for (j=0; j < *old_nfields; j++)
      {
        pField=PyFloat_FromDouble(fields_old[i][j]);
        PyDict_SetItemString(pNames, old_field_names[j], pField);
      }

      for (j=0; j < *old_nattributes; j++)
      {
        pField=PyFloat_FromDouble(attributes_old[i][j]);
        PyDict_SetItemString(pNames, old_att_names[j], pField);
      }

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      pResult=PyObject_CallObject(pFunc, pArgs);

      // Check for a Python error in the function call
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      pArray = (PyArrayObject *)
	PyArray_ContiguousFromObject(pResult, NPY_DOUBLE, 2, 2);

      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

      if (PyArray_DIMS(pArray)[0] != dim || PyArray_DIMS(pArray)[1] != dim)
	{
	  fprintf(stderr, "Error: dimensions of array returned from python ([%d, %d]) do not match allocated dimensions of the tensor_field ([%d, %d])).\n",
		  (int) PyArray_DIMS(pArray)[0], (int) PyArray_DIMS(pArray)[1], dim, dim);
	  *stat=1;
	  return;
	}

      for (j = 0; j < dim; j++){
	for (k = 0; k < dim; k++){

	  // Note the transpose for fortran.
	  double tmp;
	  tmp = *(double*)(PyArray_DATA(pArray) + j * PyArray_STRIDES(pArray)[0] + k * PyArray_STRIDES(pArray)[1]);
	  result[j][k][i] = tmp;
	}
      }

      Py_DECREF(pArray);

      // Check for a Python error in result.
      if (PyErr_Occurred()){
	PyErr_Print();
	*stat=1;
	return;
      }

    }
  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  //Free allocated memory
  free(fields_new);
  free(fields_old);
  free(attributes_old);

  // Force a garbage collection
  PyGC_Collect();

  *stat=0;
  return;
#endif
}

#define real_from_python F77_FUNC(real_from_python, REAL_FROM_PYTHON)
void real_from_python(char* function, int* function_len,
                        double* t,
                        double* result, int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = PyFloat_AsDouble(pResult);

  // Check for a Python error in result.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}

void free_c_vector(void** vector)
{
  free(*vector);
}

void real_vector_from_python(char* function, int* function_len,
                             double* t,
                             void** result,
                             int* result_len,
                             int* stat)
{
 int i;
#ifndef HAVE_PYTHON
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT, *pResultItem;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result_len = PySequence_Length(pResult);

  // Check for a Python error in result_dim.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = malloc(*result_len * sizeof(double));

  // Unpack tuple to pointer
  for (i = 0; i < *result_len; i++){
    pResultItem = PySequence_GetItem(pResult, i);
    // Check for a Python error in unpacking tuple.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    ((double*)*result)[i]=PyFloat_AsDouble(pResultItem);

    // Check we really got a float.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    Py_DECREF(pResultItem);
  }


  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}

void integer_vector_from_python(char* function, int* function_len,
                             double* t,
                             void** result,
                             int* result_len,
                             int* stat)
{
 int i;
#ifndef HAVE_PYTHON
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT, *pResultItem;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result_len = PySequence_Length(pResult);

  // Check for a Python error in result_dim.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = malloc(*result_len * sizeof(long));

  // Unpack tuple to pointer
  for (i = 0; i < *result_len; i++){
    pResultItem = PySequence_GetItem(pResult, i);
    // Check for a Python error in unpacking tuple.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    ((long*)*result)[i]=PyInt_AsLong(pResultItem);

    // Check we really got a float.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    Py_DECREF(pResultItem);
  }


  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}

#define integer_from_python F77_FUNC(integer_from_python, INTEGER_FROM_PYTHON)
void integer_from_python(char* function, int* function_len,
                        double* t,
                        int* result, int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = PyLong_AsLong(pResult);

  // Check for a Python error in result.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}

#define string_from_python F77_FUNC(string_from_python, STRING_FROM_PYTHON)
void string_from_python(char* function, int* function_len,
                        int* result_len,
                        double* t,
                        char* result, int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT;
  int pResult_len;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult_len = PyString_Size(pResult);
  if(pResult_len > *result_len){
    fprintf(stderr, "In string_from_python\n");
    fprintf(stderr, "Warning: Truncating returned string\n");
    fflush(stderr);
    memcpy(result, PyString_AsString(pResult), *result_len * sizeof(char));
  }else{
    memcpy(result, PyString_AsString(pResult), pResult_len * sizeof(char));
    *result_len = pResult_len;
  }

  // Check for a Python error in result.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}
