#pragma once

#include "hmathconfig.h"
#include "hmathconstants.h"
#include "hmathtypes.h"

#include <cmath>
#include <functional>
#include <optional>
#include <string>
#include <vector>


namespace hmath
{
	namespace analysis
	{
		// Standard function with a single parameter and having the same domain and range.
		// f:x -> y, where x and y are real numbers.
		using TFunc1 = std::function<HReal(HReal)>;

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

			std::optional<HRoot> bisectionMethod(int& outPerformanceCount,
				TFunc1 continuousFunc, HReal start, HReal end,
				int maxCount = 30, HReal epsilon = SMALL_NUMBER);
		}

		HReal getError(HReal approximateValue, HReal trueValue);
		HReal getRelativeError(HReal approximateValue, HReal trueValue);

#if DO_TEST
		int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages);
#endif // DO_TEST

	} // analysis
} // hmath
