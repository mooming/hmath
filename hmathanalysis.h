#pragma once

#include "hmathconfig.h"
#include "hmathconstants.h"
#include "hmathtypes.h"

#include <cmath>
#include <optional>
#include <string>
#include <vector>


namespace hmath
{
	namespace analysis
	{
		TFunc1 composite(TFunc1 func1, TFunc1 func2);
		HReal derivativeFromBelow(TFunc1 func, HReal x, HReal epsilon = EPSILON);
		HReal derivativeFromAbove(TFunc1 func, HReal x, HReal epsilon = EPSILON);
		HReal derivative(TFunc1 func, HReal x, HReal epsilon = EPSILON);

		HReal secondOrderDerivativeFromBelow(TFunc1 func, HReal x, HReal epsilon = EPSILON);
		HReal secondOrderDerivativeFromAbove(TFunc1 func, HReal x, HReal epsilon = EPSILON);
		HReal secondOrderDerivative(TFunc1 func, HReal x, HReal epsilon = EPSILON);

		TFunc1 getDerivativeFromBelow(TFunc1 func, HReal epsilon = EPSILON);
		TFunc1 getDerivativeFromAbove(TFunc1 func, HReal epsilon = EPSILON);
		TFunc1 getDerivative(TFunc1 func, HReal epsilon = EPSILON);

		TFunc1 getSecondOrderDerivativeFromBelow(TFunc1 func, HReal epsilon = EPSILON);
		TFunc1 getSecondOrderDerivativeFromAbove(TFunc1 func, HReal epsilon = EPSILON);
		TFunc1 getSecondOrderDerivative(TFunc1 func, HReal epsilon = EPSILON);

		namespace solver
		{
			struct HRoot final
			{
				const HReal value = ZERO;
				const HReal error = ZERO;
			};

			struct HQuadraticRoots final
			{
				const HReal first = ZERO;
				const HReal second = ZERO;
			};

			std::optional<HReal> solveLinearEquation(HReal a, HReal b);
			std::optional<HQuadraticRoots> solveQuadraticEquation(HReal a, HReal b, HReal c);

			// conditions
			// The given function should be continous on range [start, end].
			// The sign of f(start) and f(end) should be different.
			// If f(start) is positive, then f(end) should be negative.
			std::optional<HRoot> bisectionMethod(int& outIterationCount,
				TFunc1 continuousFunc, HReal start, HReal end,
				int maxCount = 30, HReal epsilon = SMALL_NUMBER);

			// conditions
			// The given function should be differentiable for every point.
			// y`(start) should not be zero.
			std::optional<HRoot> newtonRaphsonMethod(int& outIterationCount,
				TFunc1 func, TFunc1 derivativeFunc, HReal start,
				int maxCount = 30, HReal epsilon = SMALL_NUMBER);

			std::optional<HRoot> newtonRaphsonMethod(int& outIterationCount,
				TFunc1 differentiableFunc, HReal start,
				int maxCount = 30, HReal epsilon = SMALL_NUMBER);

			// conditions
			// The given function should be differentiable for every point.
			// y`(start) should not be zero.
			std::optional<HRoot> secantMethod(int& outIterationCount,
				TFunc1 differentiableFunc, HReal start, HReal start2,
				int maxCount = 30, HReal epsilon = SMALL_NUMBER);
		}

		HReal getError(HReal approximateValue, HReal trueValue);
		HReal getRelativeError(HReal approximateValue, HReal trueValue);

#if DO_TEST
		int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages);
#endif // DO_TEST

	} // analysis
} // hmath
