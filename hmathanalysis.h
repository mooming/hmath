#pragma once

#include "hmathconfig.h"
#include "hmathconstants.h"
#include "hmathtypes.h"

#include <optional>
#include <string>
#include <vector>


namespace hmath
{
namespace analysis
{
	static constexpr HReal DERIVATIVE_STEP = 1e-4;

	HReal derivativeFromBelow(const TFunc1& func, HReal x, HReal epsilon = DERIVATIVE_STEP);
	HReal derivativeFromAbove(const TFunc1& func, HReal x, HReal epsilon = DERIVATIVE_STEP);
	HReal derivative(const TFunc1& func, HReal x, HReal epsilon = DERIVATIVE_STEP);

	HReal secondOrderDerivativeFromBelow(const TFunc1& func, HReal x, HReal epsilon = DERIVATIVE_STEP);
	HReal secondOrderDerivativeFromAbove(const TFunc1& func, HReal x, HReal epsilon = DERIVATIVE_STEP);
	HReal secondOrderDerivative(const TFunc1& func, HReal x, HReal epsilon = DERIVATIVE_STEP);

	TFunc1 getDerivativeFromBelow(const TFunc1& func, HReal epsilon = DERIVATIVE_STEP);
	TFunc1 getDerivativeFromAbove(const TFunc1& func, HReal epsilon = DERIVATIVE_STEP);
	TFunc1 getDerivative(const TFunc1& func, HReal epsilon = DERIVATIVE_STEP);

	TFunc1 getSecondOrderDerivativeFromBelow(const TFunc1& func, HReal epsilon = DERIVATIVE_STEP);
	TFunc1 getSecondOrderDerivativeFromAbove(const TFunc1& func, HReal epsilon = DERIVATIVE_STEP);
	TFunc1 getSecondOrderDerivative(const TFunc1& func, HReal epsilon = DERIVATIVE_STEP);

	HReal getError(HReal approximateValue, HReal trueValue);
	HReal getRelativeError(HReal approximateValue, HReal trueValue);

	// conditions
	// The given function should be continous on range [start, end].
	// The sign of f(start) and f(end) should be different.
	// If f(start) is positive, then f(end) should be negative.
	std::optional<HRoot> bisectionMethod(int& outIterationCount,
		const TFunc1& continuousFunc, HReal start, HReal end,
		int maxCount = 30, HReal epsilon = SMALL_NUMBER);

	// conditions
	// The given function should be differentiable for every point.
	// y`(start) should not be zero.
	std::optional<HRoot> newtonRaphsonMethod(int& outIterationCount,
		const TFunc1& func, const TFunc1& derivativeFunc, HReal start,
		int maxCount = 30, HReal epsilon = SMALL_NUMBER);

	std::optional<HRoot> newtonRaphsonMethod(int& outIterationCount,
		const TFunc1& differentiableFunc, HReal start,
		int maxCount = 30, HReal epsilon = SMALL_NUMBER);

	// conditions
	// The given function should be differentiable for every point.
	// y`(start) should not be zero.
	std::optional<HRoot> secantMethod(int& outIterationCount,
		const TFunc1& differentiableFunc, HReal start, HReal start2,
		int maxCount = 30, HReal epsilon = SMALL_NUMBER);

#if DO_TEST
	int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages);
#endif // DO_TEST

} // analysis

} // hmath
