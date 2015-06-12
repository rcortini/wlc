#include <Python.h>
#include <wlc/wlc.h>

static PyObject* g_F (PyObject* self, PyObject *args)
{
  double F, lpb, g;
  if (!PyArg_ParseTuple(args, "dd", &F, &lpb))
    return NULL;

  // calculate wlc_g_F
  g = wlc_g_F (F, lpb);
  return Py_BuildValue("d", g);
}

static char g_F_docs[] = "g_F: calculate Gibbs free energy as a function of force\n";

/* now we initialize the module */
static PyMethodDef wlc_funcs[] = {
  {"g_F", (PyCFunction)g_F, METH_VARARGS, g_F_docs},
  {NULL}
};

void initwlc(void)
{
  Py_InitModule3 ("wlc", wlc_funcs, "wlc module, wraps wlc C library");
}
