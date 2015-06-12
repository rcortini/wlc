#include <Python.h>
#include <wlc/wlc.h>

// wlc_g_F
static PyObject* pywlc_g_F (PyObject* self, PyObject *args)
{
  double F, lpb;
  if (!PyArg_ParseTuple(args, "dd", &F, &lpb))
    return NULL;

  // calculate wlc_g_F
  double g = wlc_g_F (F, lpb);
  return Py_BuildValue("d", g);
}

static char pywlc_g_F_docs[] = "wlc_g_F: calculate Gibbs free energy as a function of force\n";

// wlc_f_rho
static PyObject* pywlc_f_rho (PyObject* self, PyObject *args)
{
  double rho, lpb;
  if (!PyArg_ParseTuple(args, "dd", &rho, &lpb))
    return NULL;

  // calculate wlc_g_F
  double f = wlc_f_rho (rho, lpb);
  return Py_BuildValue("d", f);
}

static char pywlc_f_rho_docs[] = "wlc_f_rho: calculate Helmholtz free energy as a function of z/L\n";

// now we initialize the module
static PyMethodDef pywlc_funcs[] = {
  {"wlc_g_F", (PyCFunction)pywlc_g_F, METH_VARARGS, pywlc_g_F_docs},
  {"wlc_f_rho", (PyCFunction)pywlc_f_rho, METH_VARARGS, pywlc_f_rho_docs},
  {NULL}
};

void initpywlc(void)
{
  Py_InitModule3 ("pywlc", pywlc_funcs, "pywlc module, wraps wlc C library");
}
