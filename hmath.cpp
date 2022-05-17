
#include "hmathanalysis.h"
#include "hmathconfig.h"
#include "hmathconstants.h"
#include "hmathpolynomial.h"

#include <iostream>
#include <limits>
#include <string>
#include <vector>

#if DO_TEST

int main()
{
	using namespace std;
	using namespace hmath;

	int testCount = 0;
	int errorCount = 0;
	vector<string> errorMessages;
	
	cout << "[HMath] Test Started! ===" << endl;

	cout << "[HMath] Derivative tests on trigonometric functions" << endl;
	{	
		auto func = [](HReal value) -> HReal { return sin(value); };

		for (int i = 0; i <= 360; i+=30)
		{
			const HReal radians = degreesToRadians(i);

			cout << "[HMath] Sin(" << i << ") = " << func(radians) << endl;

			HReal value = derivative(func, radians);
			HReal trueValue = cos(radians);
			HReal errorValue = getError(value, trueValue) * 0.5;

			cout << "[HMath] Sin`(" << i << ") = " << derivativeFromBelow(func, radians)
				<< ", " << derivativeFromAbove(func, radians) << ", " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;
			
			value = secondOrderDerivative(func, radians);
			trueValue = -sin(radians);
			errorValue = getError(value, trueValue) * 0.5; // Range is [-1, 1]. Hence rational error is absolute error / 2.

			cout << "[HMath] Sin``(" << i << ") = " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;
		}
	}

	cout << "[HMath] Derivative tests on exponential functions" << endl;
	{
		constexpr float a = PI;
		auto func = [a](HReal value) -> HReal { return exp(a * value); };

		for (int i = 0; i < 100; ++i)
		{
			
			const HReal x = i;

			cout << "[HMath] Exp(" << i << ") = " << func(x) << endl;

			HReal value = derivative(func, x);
			HReal trueValue = exp(a * x) * a;
			HReal errorValue = getRelativeError(value, trueValue);

			cout << "[HMath] Exp`(" << i << ") = " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;

			auto derivativeFunc = getDerivative(func);
			value = secondOrderDerivative(func, x);
			
			trueValue = exp(a * x) * a * a;
			errorValue = getRelativeError(value, trueValue);

			cout << "[HMath] Exp``(" << i << ") = " << value
				<< ", true value = " << trueValue
				<< ", error = " << errorValue << endl;
		}
	}

	Polynomial::DoTest(testCount, errorCount, errorMessages);

	cout << endl;
	cout << "[HMath] Test Finished! ===" << endl;
	cout << "[HMath] Test Report" << endl;
	cout << "[HMath] Errors " << errorCount << endl;
	cout << endl;
	
	return 0;
}

#endif // DO_TEST
