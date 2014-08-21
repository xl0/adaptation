
#include <stdio.h>
#include <Python.h>


static inline unsigned long find_match(const char *seq, const char *template)
{
	unsigned long i;

	i = 0;

	while (seq[i] && template[i]) {
		if (seq[i] != 'N' && seq[i] != template[i])
			break;
		i++;
	}

	return !seq[i];
}

static PyObject* find_seq(PyObject* self, PyObject* args)
{
	const char *seq;
	const char *template;
	unsigned long n = 0;
	PyObject *res_list;


	if (!PyArg_ParseTuple(args, "ss", &seq, &template))
	        return NULL;

	res_list = PyList_New(0);
	if (!res_list)
		return NULL;

	while (template[n]) {
		if (find_match(seq, template + n))
			if (-1 == PyList_Append(res_list, Py_BuildValue("k", n)))
				return NULL;
		n++;
	}


	return res_list;
}

static PyMethodDef FindSeqMethods[] =
{
     {"find_seq", find_seq, METH_VARARGS, "Find sub-sequence in larger sequence"},
     {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC

initfind_seq_module(void)
{
     (void) Py_InitModule("find_seq_module", FindSeqMethods);
}

