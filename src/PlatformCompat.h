#pragma once
// Small shim for the handful of MSVC CRT calls Giraffe uses, so the same
// call sites work on macOS/Linux too.
#ifdef _WIN32
#include <direct.h>
#define GIRAFFE_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#include <cmath>
#include <cstdio>
#define GIRAFFE_MKDIR(path) mkdir(path, 0755)
#define _copysign copysign
// MSVC's _itoa(value, buffer, 10); the codebase only ever uses base 10.
static inline char* _itoa(int value, char* buffer, int base)
{
	sprintf(buffer, "%d", value);
	return buffer;
}
#endif
