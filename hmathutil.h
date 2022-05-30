#pragma once

#include "hmathconfig.h"
#include "hmathconstants.h"
#include "hmathtypes.h"

#include <optional>


namespace hmath
{

namespace util
{
	TFunc1 composite(const TFunc1& func1, const TFunc1& func2);
	TFunc1 operator+(const TFunc1& left, const TFunc1& right);
	TFunc1 operator-(const TFunc1& left, const TFunc1& right);
	TFunc1 operator*(const TFunc1& left, HReal right);
	TFunc1 operator*(HReal left, const TFunc1& right);
	
	HReal compare(const TFunc1& func1, const TFunc1& func2, HReal start, HReal end, HReal step = SMALL_NUMBER);
	std::optional<HReal> solveLinearEquation(HReal a, HReal b);
	std::optional<std::pair<HReal, HReal>> solveQuadraticEquation(HReal a, HReal b, HReal c);
#if DO_TEST
	int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages);
#endif // DO_TEST
} // util

} // hmath
