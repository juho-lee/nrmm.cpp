#ifndef DEF_H_
#define DEF_H_

#ifdef _WIN32
#if _MSC_VER >= 1800
#include <cmath>
#define INF INFINITY
#else
#include <limits>
#define INF std::numeric_limits<double>::infinity()
#endif
#endif

#ifdef _WIN32
#if (_MSC_VER >= 1800)
#define foreach(it, container) for (auto &it : container)
#else
#define foreach(it, container) for each (auto &it in container)
#endif
#endif

#define ElapsedSecs(bc) ((double)(clock() - bc)/CLOCKS_PER_SEC)

#define BUFFSIZE 1000
#define PRECISION 10

#define PI 3.141592653589793238462

static const char *datapath = "../data";
static const char *resultspath = "../results";

#endif