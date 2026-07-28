#include <stddef.h>

#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

void attribute_visible
R_init_Rduckhts(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
