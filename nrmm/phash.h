#ifndef PHASH_H_
#define PHASH_H_

#include <xhash>

// http://stackoverflow.com/a/16473277
template <class T>
void hash_combine(size_t &seed, const T &v)
{
#ifdef WIN32
	seed ^= stdext::hash_value<T>(v) +0x9e3779b9 + (seed << 6) + (seed >> 2);
#else
	seed ^= hash<T>()(v)+0x9e3779b9 + (seed << 6) + (seed >> 2);
#endif
}

template <class T>
size_t phash(const T &v0, const T &v1)
{
#ifdef WIN32
	size_t retval = stdext::hash_value<T>(v0);
#else
	size_t retval = hash<T>()(v0);
#endif
	hash_combine(retval, v1);
	return retval;
}

#endif