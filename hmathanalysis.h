#pragma once

#include "hmathconstants.h"
#include "hmathtypes.h"

#include <cmath>
#include <functional>


namespace hmath
{
	// Standard function with a single parameter and having the same domain and range.
	// f:x -> y, where x and y are real numbers.
	using TFunc1 = std::function<HReal(HReal)>;

	HReal derivativeFromBelow(TFunc1 func, HReal x, HReal epsilon = Epsilon);
	HReal derivativeFromAbove(TFunc1 func, HReal x, HReal epsilon = Epsilon);
	HReal derivative(TFunc1 func, HReal x, HReal epsilon = Epsilon);
	
	HReal secondOrderDerivativeFromBelow(TFunc1 func, HReal x, HReal epsilon = Epsilon);
	HReal secondOrderDerivativeFromAbove(TFunc1 func, HReal x, HReal epsilon = Epsilon);
	HReal secondOrderDerivative(TFunc1 func, HReal x, HReal epsilon = Epsilon);

	TFunc1 getDerivativeFromBelow(TFunc1 func, HReal epsilon = Epsilon);
	TFunc1 getDerivativeFromAbove(TFunc1 func, HReal epsilon = Epsilon);
	TFunc1 getDerivative(TFunc1 func, HReal epsilon = Epsilon);
	
	TFunc1 getSecondOrderDerivativeFromBelow(TFunc1 func, HReal epsilon = Epsilon);
	TFunc1 getSecondOrderDerivativeFromAbove(TFunc1 func, HReal epsilon = Epsilon);
	TFunc1 getSecondOrderDerivative(TFunc1 func, HReal epsilon = Epsilon);

	HReal getError(HReal approximateValue, HReal trueValue);
	HReal getRelativeError(HReal approximateValue, HReal trueValue);
}
