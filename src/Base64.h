#pragma once
#ifndef _BASE64_H_
#define _BASE64_H_

//#include <vector>
#include <string>

std::string b64encode(const void* data, const size_t len);
std::string b64decode(const void* data, const size_t len);
std::string b64decode(const std::string& str64);

#endif