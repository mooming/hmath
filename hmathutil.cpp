#include "hmathutil.h"

#include "hmathbitops.h"
#include <cmath>
#include <iostream>
#include <sstream>


namespace hmath
{
namespace util
{
	TFunc1 composite(const TFunc1& func1, const TFunc1& func2)
	{
		if (!func1)
		{
			using namespace std;
			cerr << "[hmath][Error] " << __func__ << ": func1 is null." << endl;

			return TFunc1();
		}

		if (!func2)
		{
			using namespace std;
			cerr << "[hmath][Error] " << __func__ << ": func2 is null." << endl;

			return TFunc1();
		}

		return [func1, func2](HReal value)
			{
				auto y1 = func1(value);

				return func2(y1);
			};
	}

	TFunc1 operator+(const TFunc1& left, const TFunc1& right)
	{
		return [left, right](HReal x) -> HReal
		{
			return left(x) + right(x);
		};
	}

	TFunc1 operator-(const TFunc1& left, const TFunc1& right)
	{
		return [left, right](HReal x) -> HReal
		{
			return left(x) - right(x);
		};
	}

	TFunc1 operator*(const TFunc1& left, HReal right)
	{
		return [left, right](HReal x) -> HReal
		{
			return left(x) * right;
		};
	}

	TFunc1 operator*(HReal left, const TFunc1& right)
	{
		return [left, right](HReal x) -> HReal
		{
			return left * right(x);
		};
	}

	HReal compare(const TFunc1& func1, const TFunc1& func2, HReal start, HReal end, HReal step)
	{
		HReal error = ZERO;

		for (HReal x = start; x < end; x += step)
		{
			auto y = func1(x);
			auto y2 = func2(x);
			auto delta = abs(y2 - y);

			error = std::max(error, delta);
		}

		return error;
	}

	std::optional<HReal> solveLinearEquation(HReal a, HReal b)
	{
		if (abs(a) < MIN_NUMBER)
			return std::optional<HReal>();

		return -b / a;
	}

	std::optional<std::pair<HReal, HReal>> solveQuadraticEquation(HReal a, HReal b, HReal c)
	{
		using TRoot = std::pair<HReal, HReal>;

		HReal d = b * b - (4 * a * c);
		if (bitops::isNegative(d))
			return std::optional<TRoot>();

		HReal sqrtD = sqrt(d);
		HReal factor = HALF / a;
		HReal root1 = (-b + d) * factor;
		HReal root2 = (-b - d) * factor;

		return TRoot{ root1, root2 };
	}

#if DO_TEST
	int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages)
	{
		using namespace std;

		int errorCount = 0;

		do
		{
			cout << "[hmathutil][TC" << ++inOutTestCount << "] Copy Function Test" << endl;

			static int copyCount = 0;
			static int moveCount = 0;

			struct CopyCounter final
			{
				int x;

				CopyCounter()
					: x(0)
				{
				}

				CopyCounter(const CopyCounter& rhs)
					: x(rhs.x)
				{
					++copyCount;
				}

				CopyCounter(CopyCounter&& rhs)
					: x(rhs.x)
				{
					++moveCount;
				}

				CopyCounter& operator= (const CopyCounter& rhs)
				{
					x = rhs.x;
					++copyCount;

					return *this;
				}

				CopyCounter& operator= (CopyCounter&& rhs) noexcept
				{
					x = rhs.x;
					++moveCount;

					return *this;
				}
			};

			CopyCounter testCounter;
			testCounter.x = 10;

			int expectedCopyCount = 0;

			cout << "[hmathutil][TC" << inOutTestCount << "] Copy Count = " << copyCount << endl;
			if (copyCount != expectedCopyCount)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount 
					<< "][Error] Incorrect copy count = " << copyCount
					<< ", expected = " << expectedCopyCount << endl;
				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);

				break;
			}

			auto func = [testCounter]()
			{
				return testCounter.x;
			};
			
			++expectedCopyCount;

			cout << "[hmathutil][TC" << inOutTestCount << "] Copy Count = " << copyCount << endl;
			if (copyCount != expectedCopyCount)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "][Error] Incorrect copy count = " << copyCount
					<< ", expected = " << expectedCopyCount << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);
				break;
			}

			func();

			cout << "[hmathutil][TC" << inOutTestCount << "] Copy Count = " << copyCount << endl;
			if (copyCount != expectedCopyCount)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "][Error] Incorrect copy count = " << copyCount
					<< ", expected = " << expectedCopyCount << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);

				break;
			}

			auto func2 = func;
			++expectedCopyCount;

			cout << "[hmathutil][TC" << inOutTestCount << "] Copy Count = " << copyCount << endl;
			if (copyCount != expectedCopyCount)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "][Error] Incorrect copy count = " << copyCount
					<< ", expected = " << expectedCopyCount << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);

				break;
			}

			func2();

			auto func3 = func2;
			++expectedCopyCount;

			cout << "[hmathutil][TC" << inOutTestCount << "] Copy Count = " << copyCount << endl;
			if (copyCount != expectedCopyCount)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "][Error] Incorrect copy count = " << copyCount
					<< ", expected = " << expectedCopyCount << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);

				break;
			}

			const auto& func4 = func3;
			func4();

			cout << "[hmathutil][TC" << inOutTestCount << "] Copy Count = " << copyCount << endl;
			if (copyCount != expectedCopyCount)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "][Error] Incorrect copy count = " << copyCount
					<< ", expected = " << expectedCopyCount << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);

				break;
			}
		} while (false);

		do
		{
			cout << "[hmathutil][TC" << ++inOutTestCount << "] Composite functiont test" << endl;

			auto func1 = [](HReal x) -> HReal { return x * x; };
			auto func2 = [](HReal x) -> HReal { return x + 3; };

			if (composite(nullptr, func2))
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "] SHOULD failed to composite two functions if func1 is null." << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);

				break;
			}

			if (composite(func1, nullptr))
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "] SHOULD failed to composite two functions if func2 is null." << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);

				break;
			}

			auto func3 = composite(func1, func2);
			if (!func3)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "] failed to composite two functions" << endl;
				
				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);

				break;
			}

			auto trueFunc = [](HReal x) -> HReal { return (x * x) + 3; };

			constexpr HReal step = 0.01;
			for (int i = 0; i < 100; ++i)
			{
				HReal x = (i - 50) * step;
				auto value = func3(x);
				auto trueValue = trueFunc(x);

				auto error = ::abs(value - trueValue);

				if (error > EPSILON)
				{
					++errorCount;

					ostringstream msg;
					msg << "[hmathutil][TC" << inOutTestCount
						<< "] composite function has the invaild outcome at " << x
						<< ", y = " << value << ", but " << trueValue << " expected" << endl;

					auto msgStr = msg.view();
					cerr << msgStr;

					outErrorMessages.emplace_back(msgStr);

					break;
				}
			}

			cout << "[hmathutil][TC" << inOutTestCount << "] Composite functiont test : Done" << endl;
		} while (false);

		{
			cout << "[hmathutil][TC" << ++inOutTestCount << "] Compare functiont pass test" << endl;

			auto func1 = [](HReal x) -> HReal { return x * x; };
			auto func2 = [](HReal x) -> HReal { return x + 3; };

			auto func3 = composite(func1, func2);
			auto trueFunc = [](HReal x) -> HReal { return (x * x) + 3; };

			auto error = compare(func3, trueFunc, -1, 1, 0.001);
			if (error > EPSILON)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "] compare function test failed with error = " << error << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);
			}

			cout << "[hmathutil][TC" << inOutTestCount << "] Compare function test: Done, error = "
				<< error << endl;
		}

		{
			cout << "[hmathutil][TC" << ++inOutTestCount << "] Compare functiont fail test" << endl;

			auto func = [](HReal x) -> HReal { return x; };
			auto trueFunc = [](HReal x) -> HReal { return x * x; };

			auto error = compare(func, trueFunc, -1, 1, 0.001);
			if (error < 2)
			{
				++errorCount;

				ostringstream msg;
				msg << "[hmathutil][TC" << inOutTestCount
					<< "] compare function test passed with error = " << error
					<< ". It SHOULD be failed with error 2" << endl;

				auto msgStr = msg.view();
				cerr << msgStr;

				outErrorMessages.emplace_back(msgStr);
			}

			cout << "[hmathutil][TC" << inOutTestCount << "] Compare function test: Done, error = "
				<< error << endl;
		}

		return errorCount;
	}
#endif // DO_TEST

} // util

} // hmath
