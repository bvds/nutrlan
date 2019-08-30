#include "user_ortho.h"
//
//  Make dummy version of a user orthogonalization routine.
//
#ifdef USER_ORTHO
// Note that "in" and "out" may overlap.
void USER_ORTHO(int nrow, double *in, double *out){
    int i;
    for(i=0; i<nrow; i++){
        out[i] = in[i];
    }
}
#endif
