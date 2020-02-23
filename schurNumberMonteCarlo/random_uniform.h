//
//  random_uniform.h
//  SchurNumber
//
//  Created by wcn-cs on 10/05/2019.
//  Ce fichier en-tête déclare une macro RANDOM_UNIFORM.
//  Si le système possède une distribution BSD, cette macro redirige vers
//  la fonction arc4random_uniform, sinon elle est remplacée par random().
//

#ifndef random_uniform_h
#define random_uniform_h

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))

#include <sys/param.h>
#ifdef BSD
#define RANDOM_UNIFORM(n) arc4random_uniform(n)
#else
#define RANDOM_UNIFORM(n) rand()%n
#endif

#else
#define RANDOM_UNIFORM(n) rand()%n
#endif

#endif /* random_uniform_h */
