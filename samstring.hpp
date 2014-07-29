#ifndef _SAMSTRING_H_
#define _SAMSTRING_H_

#include <string>
#include <sstream>
#include "sam.h"

std::string get_seq(bam1_t* read);
std::string get_qual(bam1_t* read);

#endif /* _SAMSTRING_H_ */
