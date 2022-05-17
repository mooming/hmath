#include "hmathanalysis.h"

#include <cassert>


namespace hmath
{
	TFunc1 Composite(TFunc1 func1, TFunc1 func2)
	{
		return [func1, func2](HReal value)
		{
			auto y1 = func1(value);
			return func2(y1);
		};
	}

	HReal derivativeFromBelow(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		auto y = func(x) - func(x - epsilon);
		y = y / epsilon;
		
		return y;
	}

	HReal derivativeFromAbove(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		auto y = func(x + epsilon) - func(x);
		y = y / epsilon;

		return y;
	}

	HReal derivative(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		const auto halfEpsilon = epsilon * HALF;
		auto y = func(x + halfEpsilon) - func(x - halfEpsilon);
		y = y / epsilon;

		return y;
	}

	HReal secondOrderDerivativeFromBelow(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		auto dy = getDerivativeFromBelow(func, epsilon);
		auto ddy = derivativeFromBelow(dy, x, epsilon);

		return ddy;
	}

	HReal secondOrderDerivativeFromAbove(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		auto dy = getDerivativeFromAbove(func, epsilon);
		auto ddy = derivativeFromAbove(dy, x, epsilon);

		return ddy;
	}

	HReal secondOrderDerivative(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		auto dy = getDerivative(func, epsilon);
		auto ddy = derivative(dy, x, epsilon);

		return ddy;
	}

	TFunc1 getDerivativeFromBelow(TFunc1 func, HReal epsilon)
	{
		assert(epsilon > 0);

		auto dy = [func, epsilon](HReal value) -> HReal
		{
			return derivativeFromBelow(func, value, epsilon);
		};

		return dy;
	}

	TFunc1 getDerivativeFromAbove(TFunc1 func, HReal epsilon)
	{
		assert(epsilon > 0);

		auto dy = [func, epsilon](HReal value) -> HReal
		{
			return derivativeFromAbove(func, value, epsilon);
		};

		return dy;
	}

	TFunc1 getDerivative(TFunc1 func, HReal epsilon)
	{
		assert(epsilon > 0);

		auto derivativeFunc = [func, epsilon](HReal value) -> HReal
		{
			return derivative(func, value, epsilon);
		};

		return derivativeFunc;
	}

	TFunc1 getSecondOrderDerivativeFromBelow(TFunc1 func, HReal epsilon)
	{
		assert(epsilon > 0);

		auto derivativeFunc = [func, epsilon](HReal value) -> HReal
		{
			return secondOrderDerivativeFromBelow(func, value, epsilon);
		};

		return derivativeFunc;
	}

	TFunc1 getSecondOrderDerivativeFromAbove(TFunc1 func, HReal epsilon)
	{
		assert(epsilon > 0);

		auto derivativeFunc = [func, epsilon](HReal value) -> HReal
		{
			return secondOrderDerivativeFromAbove(func, value, epsilon);
		};

		return derivativeFunc;
	}

	TFunc1 getSecondOrderDerivative(TFunc1 func, HReal epsilon)
	{
		assert(epsilon > 0);

		auto derivativeFunc = [func, epsilon](HReal value) -> HReal
		{
			return secondOrderDerivative(func, value, epsilon);
		};

		return derivativeFunc;
	}

	HReal getError(HReal approximateValue, HReal trueValue)
	{
		return abs(approximateValue - trueValue);
	}

	HReal getRelativeError(HReal approximateValue, HReal trueValue)
	{
		auto diff = abs(approximateValue - trueValue);
		if (diff < MIN_NUMBER)
			return ZERO;

		if (abs(trueValue) < SMALL_NUMBER)
			return MAX_NUMBER;

		return diff / trueValue;
	}
}
