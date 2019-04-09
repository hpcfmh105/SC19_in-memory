#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

int
main(int argc, char **argv)
{
    void *handle;
    double (*cosine)(double);
    char *error;
	printf ("start load so\n");

	printf("test\n");
	const char* s = getenv("LD_LIBRARY_PATH");
	printf("LD PATH :%s\n",(s!=NULL)? s : "getenv returned NULL");
   handle = dlopen("/lustre/atlas/proj-shared/csc143/dhuang/workflow-bench-master/build/lib/libmod_lammps_decaf_shared.so", RTLD_LAZY);
	   
	printf("end test\n");
 if (!handle) {

		 printf("undefine symbols end test exit\n");

        fprintf(stderr, "error in testdl%s\n", dlerror());
        exit(EXIT_FAILURE);
    }

   dlerror();    /* Clear any existing error */

   /* Writing: cosine = (double (*)(double)) dlsym(handle, "cos");
       would seem more natural, but the C99 standard leaves
       casting from "void *" to a function pointer undefined.
       The assignment used below is the POSIX.1-2003 (Technical
       Corrigendum 1) workaround; see the Rationale for the
       POSIX specification of dlsym(). */


   if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        exit(EXIT_FAILURE);
    }

    dlclose(handle);
    exit(EXIT_SUCCESS);
}
