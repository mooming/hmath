#include "hmathanalysis.h"

#include <cassert>
#include <iostream>


namespace hmath
{
namespace analysis
{
TFunc1 composite(TFunc1 func1, TFunc1 func2)
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

namespace solver
{

std::optional<HReal> solveLinearEquation(HReal a, HReal b)
{
	if (abs(a) < MIN_NUMBER)
		return std::optional<HReal>();

	return -b / a;
}

std::optional<HQuadraticRoots> solveQuadraticEquation(HReal a, HReal b, HReal c)
{
	HReal d = b * b - (4 * a * c);
	if (signbit(d))
		return std::optional<HQuadraticRoots>();

	HReal sqrtD = sqrt(d);
	HReal factor = HALF / a;
	HReal root1 = (-b + d) * factor;
	HReal root2 = (-b - d) * factor;

	return HQuadraticRoots{ root1, root2 };
}

std::optional<HRoot> bisectionMethod(int& outIterationCount,
	TFunc1 continuousFunc, HReal start, HReal end,
	int maxCount, HReal epsilon)
{
	outIterationCount = 0;

	HReal sy = continuousFunc(start);
	HReal ey = continuousFunc(end);
	
	HReal error = abs(sy);
	if (error < epsilon)
		return HRoot{ start, error };

	error = abs(ey);
	if (error < epsilon)
		return HRoot{ end, error };

	bool bSyNegative = signbit(sy);
	bool bEyNegative = signbit(ey);

	if (signbit(sy) == signbit(ey))
		return std::optional<HRoot>();

	HReal deltaRange = end - start;
	
	while (!signbit(deltaRange) && abs(deltaRange) > epsilon && outIterationCount < maxCount)
	{
		++outIterationCount;

		HReal x = start + (deltaRange * HALF);
		HReal y = continuousFunc(x);

		error = abs(y);
		if (error < epsilon)
			return HRoot{ x, error };

		bool bYNegative = signbit(y);
		if (bSyNegative != bYNegative)
		{
			end = x;
			ey = y;
			bEyNegative = bYNegative;
		}
		else
		{
			start = x;
			sy = y;
			bEyNegative = bYNegative;
		}

		deltaRange = end - start;
	} 

	return std::optional<HRoot>();
}

std::optional<HRoot> newtonRaphsonMethod(int& outIterationCount,
	TFunc1 func, TFunc1 derivativeFunc, HReal start,
	int maxCount, HReal epsilon)
{
	outIterationCount = 0;

	auto x = start;
	auto y = func(x);
	auto error = abs(y);
	if (error < epsilon)
		return HRoot{ x, error };

	while (outIterationCount < maxCount)
	{
		++outIterationCount;

		auto dy = derivativeFunc(x);
		if (abs(dy) < DIV_EPSILON)
			return std::optional<HRoot>();

		x = x - (y / dy);
		y = func(x);

		error = abs(y);
		if (error < epsilon)
			return HRoot{ x, error };
	}

	return std::optional<HRoot>();
}

std::optional<HRoot> newtonRaphsonMethod(int& outIterationCount,
	TFunc1 differentiableFunc, HReal start,
	int maxCount, HReal epsilon)
{
	return newtonRaphsonMethod(outIterationCount, differentiableFunc,
		getDerivative(differentiableFunc), start, maxCount, epsilon);
}

std::optional<HRoot> secantMethod(int& outIterationCount,
	TFunc1 func, HReal start, HReal start2,
	int maxCount, HReal epsilon)
{
	outIterationCount = 0;

	auto x1 = start;
	auto x2 = start2;
	auto y1 = func(x1);
	auto y2 = func(x2);

	auto error = abs(y1);
	if (error < epsilon)
		return HRoot{ x1, error };

	error = abs(y2);
	if (error < epsilon)
		return HRoot{ x2, error };

	while (outIterationCount < maxCount)
	{
		++outIterationCount;

		auto dx = x2 - x1;
		if (abs(dx) < DIV_EPSILON)
			return std::optional<HRoot>();
		
		auto dy = (y2 - y1) / dx;
		if (abs(dy) < DIV_EPSILON)
			return std::optional<HRoot>();

		auto oldX = x1;
		auto oldY = y1;

		x1 = x2;
		y1 = y2;

		x2 = x2 - (y2 / dy);
		y2 = func(x2);

		error = abs(y2);
		if (error < epsilon)
			return HRoot{ x2, error };
	}

	return std::optional<HRoot>();
}

} // namespace solver
		
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

#if DO_TEST
int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages)
{
	using namespace std;

	int errorCount = 0;

	{
		++inOutTestCount;

		cout << endl << "[Analysis][TC" << inOutTestCount
			<< "] Derivative tests on trigonometric functions" << endl;

		auto func = [](HReal value) -> HReal { return sin(value); };

		for (int i = 0; i <= 360; i += 30)
		{
			const HReal radians = degreesToRadians(i);

			cout << "[Analysis][TC" << inOutTestCount
				<< "] Sin(" << i << ") = " << func(radians) << endl;

			HReal value = derivative(func, radians);
			HReal trueValue = cos(radians);
			HReal errorValue = getError(value, trueValue) * 0.5;

			cout << "[Analysis][TC" << inOutTestCount
				<< "] Sin`(" << i << ") = " << derivativeFromBelow(func, radians)
				<< ", " << derivativeFromAbove(func, radians) << ", " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;

			if (abs(errorValue) > SMALL_NUMBER)
			{
				++errorCount;
				cerr << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << errorValue << " is huge than expect "
					<< SMALL_NUMBER << endl << endl;
			}

			value = secondOrderDerivative(func, radians);
			trueValue = -sin(radians);

			// Range is [-1, 1]. Hence rational error is absolute error / 2.
			errorValue = getError(value, trueValue) * 0.5;

			cout << "[Analysis][TC" << inOutTestCount
				<< "] Sin``(" << i << ") = " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;

			if (abs(errorValue) > SMALL_NUMBER)
			{
				++errorCount;
				cerr << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << errorValue << " is huge than expect "
					<< SMALL_NUMBER << endl << endl;
			}
		}
	}

	{
		constexpr HReal maxError = 1e-3;

		++inOutTestCount;

		cout << endl << "[Analysis][TC" << inOutTestCount
			<< "] Derivative tests on exponential functions" << endl;

		constexpr auto a = PI;
		auto func = [a](HReal value) -> HReal { return exp(a * value); };

		for (int i = 0; i < 10; ++i)
		{
			const HReal x = i;

			cout << "[Analysis][TC" << inOutTestCount
				<< "] Exp(" << i << ") = " << func(x) << endl;

			HReal value = derivative(func, x);
			HReal trueValue = exp(a * x) * a;
			HReal errorValue = getRelativeError(value, trueValue);

			cout << "[Analysis][TC" << inOutTestCount
				<< "] Exp`(" << i << ") = " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;
					
			if (abs(errorValue) > maxError)
			{
				++errorCount;
				cerr << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << errorValue << " is huge than expect "
					<< maxError << endl << endl;
			}

			auto derivativeFunc = getDerivative(func);
			value = secondOrderDerivative(func, x);

			trueValue = exp(a * x) * a * a;
			errorValue = getRelativeError(value, trueValue);

			cout << "[Analysis][TC" << inOutTestCount
				<< "] Exp``(" << i << ") = " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;

			if (abs(errorValue) > maxError)
			{
				++errorCount;
				cerr << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << errorValue << " is huge than expect "
					<< maxError << endl << endl;
			}
		}
	}

	{
		++inOutTestCount;

		cout << endl << "[Analysis][TC" << inOutTestCount
			<< "] Solver: Bisectional Method" << endl;

		auto func = [](HReal x) -> HReal
		{
			return 2 * (x * x * x) + 1;
		};

		int iterationCount = 0;
		auto root = solver::bisectionMethod(iterationCount, func, -10, 10);

		if (root)
		{
			cout << "[Analysis][TC" << inOutTestCount
				<< "] root = " << root->value << ", error = " << root->error
				<< ", count = " << iterationCount << endl;

			auto value = func(root->value);
			if (abs(func(root->value)) > SMALL_NUMBER)
			{
				++errorCount;
				cerr << "[Analysis][TC" << inOutTestCount
					<< "][Error] f(" << root->value << ") = " << value
					<< " has a bigger error than expected " << SMALL_NUMBER << '!' << endl;
			}
		}
		else
		{
			++errorCount;
			cerr << "[Analysis][TC" << inOutTestCount
				<< "][Error] failed to find a root in [-10, 10] with " << iterationCount << " try." << endl;
		}
		
	}

	{
		++inOutTestCount;

		cout << endl << "[Analysis][TC" << inOutTestCount
			<< "] Solver: Newton-Raphson Method" << endl;

		auto func = [](HReal x) -> HReal
		{
			return 2 * x * x * x + 1;
		};

		auto dy = [](HReal x) -> HReal
		{
			return 6 * x * x;
		};

		int iterationCount = 0;
		auto root = solver::newtonRaphsonMethod(iterationCount, func, dy, -10);

		if (root)
		{
			cout << "[Analysis][TC" << inOutTestCount
				<< "] root = " << root->value << ", error = " << root->error
				<< ", count = " << iterationCount << endl;

			auto value = func(root->value);
			if (abs(func(root->value)) > SMALL_NUMBER)
			{
				++errorCount;
				cerr << "[Analysis][TC" << inOutTestCount
					<< "][Error] f(" << root->value << ") = " << value
					<< " has a bigger error than expected " << SMALL_NUMBER << '!' << endl;
			}
		}
		else
		{
			++errorCount;
			cerr << "[Analysis][TC" << inOutTestCount
				<< "][Error] failed to find a root with " << iterationCount << " try." << endl;
		}

	}

	{
		++inOutTestCount;

		cout << endl << "[Analysis][TC" << inOutTestCount
			<< "] Solver: Newton-Raphson Method (2)" << endl;

		auto func = [](HReal x) -> HReal
		{
			return 2 * x * x * x + 1;
		};

		int iterationCount = 0;
		auto root = solver::newtonRaphsonMethod(iterationCount, func, -10);

		if (root)
		{
			cout << "[Analysis][TC" << inOutTestCount
				<< "] root = " << root->value << ", error = " << root->error
				<< ", count = " << iterationCount << endl;

			auto value = func(root->value);
			if (abs(func(root->value)) > SMALL_NUMBER)
			{
				++errorCount;
				cerr << "[Analysis][TC" << inOutTestCount
					<< "][Error] f(" << root->value << ") = " << value
					<< " has a bigger error than expected " << SMALL_NUMBER << '!' << endl;
			}
		}
		else
		{
			++errorCount;
			cerr << "[Analysis][TC" << inOutTestCount
				<< "][Error] failed to find a root with " << iterationCount << " try." << endl;
		}
	}

	{
		++inOutTestCount;

		cout << endl << "[Analysis][TC" << inOutTestCount
			<< "] Solver: Secant Method" << endl;

		auto func = [](HReal x) -> HReal
		{
			return 2 * x * x * x + 1;
		};

		int iterationCount = 0;
		auto root = solver::secantMethod(iterationCount, func, -10, -10 + EPSILON);

		if (root)
		{
			cout << "[Analysis][TC" << inOutTestCount
				<< "] root = " << root->value << ", error = " << root->error
				<< ", count = " << iterationCount << endl;

			auto value = func(root->value);
			if (abs(func(root->value)) > SMALL_NUMBER)
			{
				++errorCount;
				cerr << "[Analysis][TC" << inOutTestCount
					<< "][Error] f(" << root->value << ") = " << value
					<< " has a bigger error than expected " << SMALL_NUMBER << '!' << endl;
			}
		}
		else
		{
			++errorCount;
			cerr << "[Analysis][TC" << inOutTestCount
				<< "][Error] failed to find a root with " << iterationCount << " try." << endl;
		}
	}

	return errorCount;
}
#endif // DO_TEST
} // namespace analysis
} // namespace hmath
