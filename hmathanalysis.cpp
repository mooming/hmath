#include "hmathanalysis.h"

#include <cassert>
#include <iostream>


namespace hmath
{
	namespace analysis
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

#if DO_TEST
		int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages)
		{
			using namespace std;

			int errorCount = 0;

			{
				++inOutTestCount;

				cout << "[Analysis][TC" << inOutTestCount
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

				cout << "[Analysis][TC" << inOutTestCount
					<< "] Derivative tests on exponential functions" << endl;

				constexpr float a = PI;
				auto func = [a](HReal value) -> HReal { return exp(a * value); };

				for (int i = 0; i < 100; ++i)
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

			return errorCount;
		}
#endif // DO_TEST
	}
}
