/* Definition of memory allocation and error checking macros */


/* Print s together with the source file name and the current line number. */
#define LOGSTATUS(s)  printf("[%s:%d]%s\n", __FILE__, __LINE__, s)

/* Check whether s is NULL or not. Quit this program if it is NULL. */
#define MYASSERT(s)  if (!(s))   {                                      \
    printf("General Assert Error at %s:line%d\n", __FILE__, __LINE__);  \
    exit(-1); \
  }

/* Check whether s is NULL or not on a memory allocation.
// Quit this program if it is NULL. */
#define MALLOC_CHECK(s)  if ((s) == NULL)   {                     \
    printf("No enough memory at %s:line%d ", __FILE__, __LINE__); \
    perror(":");                                                  \
    exit(-1);                                                     \
  }

#define SAFE_FREE(s)  if ((s) != NULL) {   \
    free(s);                               \
    s = NULL;                              \
  }

/* Set memory space starts at pointer n of size m to zero. */
#define BZERO(n,m)  memset(n, 0, m)