#include "hmathanalysis.h"

#include "hmathbitops.h"
#include "hmathutil.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>


namespace hmath
{
namespace analysis
{

HReal derivativeFromBelow(const TFunc1& func, HReal x, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return ZERO;
	}

	assert(epsilon > 0);

	auto y = func(x) - func(x - epsilon);
	y = y / epsilon;

	return y;
}

HReal derivativeFromAbove(const TFunc1& func, HReal x, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return ZERO;
	}

	assert(epsilon > 0);

	auto y = func(x + epsilon) - func(x);
	y = y / epsilon;

	return y;
}

HReal derivative(const TFunc1& func, HReal x, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return ZERO;
	}

	assert(epsilon > 0);

	const auto halfEpsilon = epsilon * HALF;
	auto y = func(x + halfEpsilon) - func(x - halfEpsilon);
	y = y / epsilon;

	return y;
}

HReal secondOrderDerivativeFromBelow(const TFunc1& func, HReal x, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return ZERO;
	}

	assert(epsilon > 0);

	auto dy = getDerivativeFromBelow(func, epsilon);
	auto ddy = derivativeFromBelow(dy, x, epsilon);

	return ddy;
}

HReal secondOrderDerivativeFromAbove(const TFunc1& func, HReal x, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return ZERO;
	}

	assert(epsilon > 0);

	auto dy = getDerivativeFromAbove(func, epsilon);
	auto ddy = derivativeFromAbove(dy, x, epsilon);

	return ddy;
}

HReal secondOrderDerivative(const TFunc1& func, HReal x, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return ZERO;
	}

	assert(epsilon > 0);

	auto dy = getDerivative(func, epsilon);
	auto ddy = derivative(dy, x, epsilon);

	return ddy;
}

TFunc1 getDerivativeFromBelow(const TFunc1& func, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return TFunc1();
	}

	assert(epsilon > 0);

	auto dy = [func, epsilon](HReal value) -> HReal
	{
		return derivativeFromBelow(func, value, epsilon);
	};

	return dy;
}

TFunc1 getDerivativeFromAbove(const TFunc1& func, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return TFunc1();
	}

	assert(epsilon > 0);

	auto dy = [func, epsilon](HReal value) -> HReal
	{
		return derivativeFromAbove(func, value, epsilon);
	};

	return dy;
}

TFunc1 getDerivative(const TFunc1& func, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return TFunc1();
	}

	assert(epsilon > 0);

	auto derivativeFunc = [func, epsilon](HReal value) -> HReal
	{
		return derivative(func, value, epsilon);
	};

	return derivativeFunc;
}

TFunc1 getSecondOrderDerivativeFromBelow(const TFunc1& func, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return TFunc1();
	}

	assert(epsilon > 0);

	auto derivativeFunc = [func, epsilon](HReal value) -> HReal
	{
		return secondOrderDerivativeFromBelow(func, value, epsilon);
	};

	return derivativeFunc;
}

TFunc1 getSecondOrderDerivativeFromAbove(const TFunc1& func, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return TFunc1();
	}

	assert(epsilon > 0);

	auto derivativeFunc = [func, epsilon](HReal value) -> HReal
	{
		return secondOrderDerivativeFromAbove(func, value, epsilon);
	};

	return derivativeFunc;
}

TFunc1 getSecondOrderDerivative(const TFunc1& func, HReal epsilon)
{
	if (!func)
	{
		using namespace std;
		cerr << "[hmath][analysis][Error] " << __func__ << ": func is null." << endl;
		return TFunc1();
	}

	assert(epsilon > 0);

	auto derivativeFunc = [func, epsilon](HReal value) -> HReal
	{
		return secondOrderDerivative(func, value, epsilon);
	};

	return derivativeFunc;
}

std::optional<HRoot> bisectionMethod(int& outIterationCount,
	const TFunc1& continuousFunc, HReal start, HReal end,
	int maxCount, HReal epsilon)
{
	using namespace bitops;

	outIterationCount = 0;

	HReal sy = continuousFunc(start);
	HReal ey = continuousFunc(end);
	
	HReal error = abs(sy);
	if (error < epsilon)
		return HRoot{ start, error };

	error = abs(ey);
	if (error < epsilon)
		return HRoot{ end, error };

	bool bSyNegative = isNegative(sy);
	bool bEyNegative = isNegative(ey);

	if (isNegative(sy) == isNegative(ey))
		return std::optional<HRoot>();

	HReal deltaRange = end - start;
	
	while (!isNegative(deltaRange) && abs(deltaRange) > epsilon && outIterationCount < maxCount)
	{
		++outIterationCount;

		HReal x = start + (deltaRange * HALF);
		HReal y = continuousFunc(x);

		error = abs(y);
		if (error < epsilon)
			return HRoot{ x, error };

		bool bYNegative = isNegative(y);
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
	const TFunc1& func, const TFunc1& derivativeFunc, HReal start,
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
	const TFunc1& differentiableFunc, HReal start,
	int maxCount, HReal epsilon)
{
	return newtonRaphsonMethod(outIterationCount, differentiableFunc,
		getDerivative(differentiableFunc), start, maxCount, epsilon);
}

std::optional<HRoot> secantMethod(int& outIterationCount,
	const TFunc1& func, HReal start, HReal start2,
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

	constexpr HReal MAX_ERROR = EPSILON;

	{
		cout << endl << "[Analysis][TC" << ++inOutTestCount
			<< "] Derivative tests on trigonometric functions" << endl;

		auto func = [](HReal value) -> HReal
		{
			return sin(value);
		};

		auto testFunc = [func](HReal value) -> HReal
		{
			return derivative(func, value);
		};

		auto trueFunc = [](HReal value) -> HReal
		{
			return cos(value);
		};

		auto error = util::compare(testFunc, trueFunc,
			-TWO_PI, TWO_PI, TWO_PI * 0.001);

		cout << "[Analysis][TC" << inOutTestCount
			<< "] Test the derivative of Sin: error = " << error << endl;

		if (error > MAX_ERROR)
		{
			++errorCount;

			ostringstream msg;
			msg << "[Analysis][TC" << inOutTestCount
				<< "][Error] " << __LINE__ << ": "
				<< error << " is huge than expect "
				<< SMALL_NUMBER << endl << endl;

			const auto errorMsg = msg.view();
			cerr << errorMsg;

			outErrorMessages.emplace_back(errorMsg);
		}
	}

	{
		cout << endl << "[Analysis][TC" << ++inOutTestCount
			<< "] Derivative tests on trigonometric functions" << endl;

		auto func = [](HReal value) -> HReal
		{
			return sin(value);
		};

		auto testFunc = [func](HReal value) -> HReal
		{
			return secondOrderDerivative(func, value);
		};

		auto trueFunc = [](HReal value) -> HReal
		{
			return -sin(value);
		};

		auto error = util::compare(testFunc, trueFunc,
			-TWO_PI, TWO_PI, TWO_PI / 360);

		cout << "[Analysis][TC" << inOutTestCount
			<< "] Test the 2nd. order derivative of Sin: error = " << error << endl;

		if (error > MAX_ERROR)
		{
			++errorCount;

			ostringstream msg;
			msg << "[Analysis][TC" << inOutTestCount
				<< "][Error] " << __LINE__ << ": "
				<< error << " is huge than expect "
				<< SMALL_NUMBER << endl << endl;

			const auto errorMsg = msg.view();
			cerr << errorMsg;

			outErrorMessages.emplace_back(errorMsg);
		}
	}

	{
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
					
			if (abs(errorValue) > MAX_ERROR)
			{
				++errorCount;

				ostringstream msg;
				msg << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << __LINE__ << ": "
					<< errorValue << " is huge than expect "
					<< MAX_ERROR << endl << endl;

				const auto errorMsg = msg.view();
				cerr << errorMsg;

				outErrorMessages.emplace_back(errorMsg);
			}

			auto derivativeFunc = getDerivative(func);
			value = secondOrderDerivative(func, x);

			trueValue = exp(a * x) * a * a;
			errorValue = getRelativeError(value, trueValue);

			cout << "[Analysis][TC" << inOutTestCount
				<< "] Exp``(" << i << ") = " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;

			if (abs(errorValue) > MAX_ERROR)
			{
				++errorCount;

				ostringstream msg;

				msg << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << __LINE__ << ": "
					<< errorValue << " is huge than expect "
					<< MAX_ERROR << endl << endl;

				const auto errorMsg = msg.view();
				cerr << errorMsg;

				outErrorMessages.emplace_back(errorMsg);
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
		auto root = bisectionMethod(iterationCount, func, -10, 10);

		if (root)
		{
			cout << "[Analysis][TC" << inOutTestCount
				<< "] root = " << root->value << ", error = " << root->error
				<< ", count = " << iterationCount << endl;

			auto value = func(root->value);
			if (abs(func(root->value)) > SMALL_NUMBER)
			{
				++errorCount;

				ostringstream msg;
				msg << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << __LINE__
					<< ": f(" << root->value << ") = " << value
					<< " has a bigger error than expected " << SMALL_NUMBER << '!' << endl;

				const auto errorMsg = msg.view();
				cerr << errorMsg;

				outErrorMessages.emplace_back(errorMsg);
			}
		}
		else
		{
			++errorCount;

			ostringstream msg;
			msg << "[Analysis][TC" << inOutTestCount
				<< "][Error] " << __LINE__
				<< ": failed to find a root in [-10, 10] with "
				<< iterationCount << " try." << endl;

			const auto errorMsg = msg.view();
			cerr << errorMsg;

			outErrorMessages.emplace_back(errorMsg);
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
		auto root = newtonRaphsonMethod(iterationCount, func, dy, -10);

		if (root)
		{
			cout << "[Analysis][TC" << inOutTestCount
				<< "] root = " << root->value << ", error = " << root->error
				<< ", count = " << iterationCount << endl;

			auto value = func(root->value);
			if (abs(func(root->value)) > SMALL_NUMBER)
			{
				++errorCount;

				ostringstream msg;
				msg << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << __LINE__
					<< ": f(" << root->value << ") = " << value
					<< " has a bigger error than expected "
					<< SMALL_NUMBER << '!' << endl;

				const auto errorMsg = msg.view();
				cerr << errorMsg;

				outErrorMessages.emplace_back(errorMsg);
			}
		}
		else
		{
			++errorCount;

			ostringstream msg;
			msg << "[Analysis][TC" << inOutTestCount
				<< "][Error] " << __LINE__
				<< ": failed to find a root with "
				<< iterationCount << " try." << endl;

			const auto errorMsg = msg.view();
			cerr << errorMsg;

			outErrorMessages.emplace_back(errorMsg);
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
		auto root = newtonRaphsonMethod(iterationCount, func, -10);

		if (root)
		{
			cout << "[Analysis][TC" << inOutTestCount
				<< "] root = " << root->value << ", error = " << root->error
				<< ", count = " << iterationCount << endl;

			auto value = func(root->value);
			if (abs(func(root->value)) > SMALL_NUMBER)
			{
				++errorCount;

				ostringstream msg;
				msg << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << __LINE__ << ": f(" << root->value << ") = " << value
					<< " has a bigger error than expected " << SMALL_NUMBER << '!' << endl;

				const auto errorMsg = msg.view();
				cerr << errorMsg;

				outErrorMessages.emplace_back(errorMsg);
			}
		}
		else
		{
			++errorCount;

			ostringstream msg;
			msg << "[Analysis][TC" << inOutTestCount
				<< "][Error] " << __LINE__
				<< ": failed to find a root with "
				<< iterationCount << " try." << endl;

			const auto errorMsg = msg.view();
			cerr << errorMsg;

			outErrorMessages.emplace_back(errorMsg);
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
		auto root = secantMethod(iterationCount, func, -10, -10 + EPSILON);

		if (root)
		{
			cout << "[Analysis][TC" << inOutTestCount
				<< "] root = " << root->value << ", error = " << root->error
				<< ", count = " << iterationCount << endl;

			auto value = func(root->value);
			if (abs(func(root->value)) > SMALL_NUMBER)
			{
				++errorCount;

				ostringstream msg;
				msg << "[Analysis][TC" << inOutTestCount
					<< "][Error] " << __LINE__
					<< ": f(" << root->value << ") = " << value
					<< " has a bigger error than expected " << SMALL_NUMBER << '!' << endl;

				const auto errorMsg = msg.view();
				cerr << errorMsg;

				outErrorMessages.emplace_back(errorMsg);
			}
		}
		else
		{
			++errorCount;

			ostringstream msg;
			msg << "[Analysis][TC" << inOutTestCount
				<< "][Error] " << __LINE__
				<< ": failed to find a root with "
				<< iterationCount << " try." << endl;

			const auto errorMsg = msg.view();
			cerr << errorMsg;

			outErrorMessages.emplace_back(errorMsg);
		}
	}

	return errorCount;
}
#endif // DO_TEST
} // namespace analysis
} // namespace hmath
