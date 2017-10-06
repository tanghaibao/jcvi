#define PY_NEW(T) \
     (((PyTypeObject*)(T))->tp_new( \
             (PyTypeObject*)(T), __pyx_empty_tuple, NULL))
