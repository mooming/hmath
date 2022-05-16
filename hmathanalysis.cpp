#include "hmathanalysis.h"

#include <cassert>


namespace hmath
{
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

		const auto halfEpsilon = epsilon * Half;
		auto y = func(x + halfEpsilon) - func(x - halfEpsilon);
		y = y / epsilon;

		return y;
	}

	HReal secondOrderDerivativeFromBelow(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		const auto halfEpsilon = epsilon * Half;
		const auto dy = derivative(func, x, epsilon);
		const auto dyBelow = derivative(func, (x - halfEpsilon), epsilon);
		auto ddy = dy - dyBelow;
		ddy = ddy / epsilon;

		return ddy;
	}

	HReal secondOrderDerivativeFromAbove(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		const auto halfEpsilon = epsilon * Half;
		const auto dyAbove = derivative(func, (x + halfEpsilon), epsilon);
		const auto dy = derivative(func, x, epsilon);
		auto ddy = dyAbove - dy;
		ddy = ddy / epsilon;

		return ddy;
	}

	HReal secondOrderDerivative(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		const auto halfEpsilon = epsilon * Half;
		const auto dyAbove = derivative(func, (x + halfEpsilon), epsilon);
		const auto dyBelow = derivative(func, (x - halfEpsilon), epsilon);
		auto ddy = dyAbove - dyBelow;
		ddy = ddy / epsilon;

		return ddy;
	}

	TFunc1 getDerivativeFromBelow(TFunc1 func, HReal epsilon)
	{
		assert(epsilon > 0);

		auto derivativeFunc = [func, epsilon](HReal value) -> HReal
		{
			return derivativeFromBelow(func, value, epsilon);
		};

		return derivativeFunc;
	}

	TFunc1 getDerivativeFromAbove(TFunc1 func, HReal x, HReal epsilon)
	{
		assert(epsilon > 0);

		auto derivativeFunc = [func, epsilon](HReal value) -> HReal
		{
			return derivativeFromAbove(func, value, epsilon);
		};

		return derivativeFunc;
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
		if (diff < SmallestNumber)
			return Zero;

		if (abs(trueValue) < SmallNumber)
			return BiggestNumber;

		return diff / trueValue;
	}
}
