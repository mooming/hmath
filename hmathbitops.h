#pragma once

#include "hmathconfig.h"
#include "hmathtypes.h"

#include <iostream>
#include <string>
#include <vector>


namespace hmath
{
namespace bitops
{

template <typename T>
typename std::enable_if<std::is_integral<T>::value, std::string>::type getBitsStrings(T value)
{
	constexpr int numBits = sizeof(value) * 8;
	constexpr T mask = static_cast<T>(static_cast<T>(1) << (numBits - 1));
	
	std::string text;
	text.reserve(numBits + 1);

	for (int i = 0; i < numBits; ++i)
	{
		char ch = (value & mask) == mask ? '1' : '0';
		text.append(1, ch);
		value = value << 1;
	}

	text.append(1, 'b');

	return text;
}

std::string getBitsStrings(float value);
std::string getBitsStrings(double value);
std::string getBitsStrings(long double value);

template <typename T>
typename std::enable_if<std::is_integral<T>::value, T>::type abs(T value)
{
	constexpr int numBits = sizeof(value) * CHAR_BIT;
	const T mask = value >> (numBits - 1);

	return (value ^ mask) - mask;
}

float abs(float value);
double abs(double value);
long double abs(long double value);

#if DO_TEST
int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages);
#endif // DO_TEST

} // bitops
} // hmath