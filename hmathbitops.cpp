#include "hmathbitops.h"


namespace hmath
{
namespace bitops
{

std::string getBitsStrings(float value)
{
	using FloatType = decltype(value);
	using UintType = uint32_t;
	static_assert(sizeof(FloatType) == sizeof(UintType));
	
	UintType iValue = reinterpret_cast<UintType&>(value);
	constexpr int numBits = sizeof(iValue) * CHAR_BIT;
	constexpr UintType mask = 1 << (numBits - 1);

	std::string text;
	text.reserve(numBits + 1);

	for (int i = 0; i < numBits; ++i)
	{
		char ch = (iValue & mask) == mask ? '1' : '0';
		text.append(1, ch);
		iValue = iValue << 1;
	}

	text.append(1, 'b');

	return text;
}

std::string getBitsStrings(double value)
{
	using FloatType = decltype(value);
	using UintType = uint64_t;
	static_assert(sizeof(FloatType) == sizeof(UintType));

	UintType iValue = reinterpret_cast<UintType&>(value);
	constexpr int numBits = sizeof(iValue) * CHAR_BIT;
	constexpr UintType mask = static_cast<UintType>(1) << (numBits - 1);

	std::string text;
	text.reserve(numBits + 1);

	for (int i = 0; i < numBits; ++i)
	{
		char ch = (iValue & mask) == mask ? '1' : '0';
		text.append(1, ch);
		iValue = iValue << 1;
	}

	text.append(1, 'b');

	return text;
}

std::string getBitsStrings(long double value)
{
	using FloatType = decltype(value);
	using UintType = uint64_t;
	static_assert(sizeof(FloatType) == sizeof(UintType));

	UintType iValue = reinterpret_cast<UintType&>(value);
	constexpr int numBits = sizeof(iValue) * CHAR_BIT;
	constexpr UintType mask = static_cast<UintType>(1) << (numBits - 1);

	std::string text;
	text.reserve(numBits + 1);

	for (int i = 0; i < numBits; ++i)
	{
		char ch = (iValue & mask) == mask ? '1' : '0';
		text.append(1, ch);
		iValue = iValue << 1;
	}

	text.append(1, 'b');

	return text;
}

float abs(float value)
{
	using FloatType = decltype(value);
	using UintType = uint32_t;

	static_assert(sizeof(FloatType) == sizeof(UintType));

	constexpr int numBits = sizeof(FloatType) * CHAR_BIT;
	constexpr int shift = numBits - 1;
	constexpr UintType mask = ~(static_cast<UintType>(1) << shift);
	UintType bitValues = reinterpret_cast<UintType&>(value) & mask;

	return reinterpret_cast<FloatType&>(bitValues);
}

double abs(double value)
{
	using FloatType = decltype(value);
	using UintType = uint64_t;

	static_assert(sizeof(FloatType) == sizeof(UintType));

	constexpr int numBits = sizeof(FloatType) * CHAR_BIT;
	constexpr int shift = numBits - 1;
	constexpr UintType mask = ~(static_cast<UintType>(1) << shift);
	UintType bitValues = reinterpret_cast<UintType&>(value) & mask;

	return reinterpret_cast<FloatType&>(bitValues);
}

long double abs(long double value)
{
	using FloatType = decltype(value);
	using UintType = uint64_t;

	static_assert(sizeof(FloatType) == sizeof(UintType));

	constexpr int numBits = sizeof(FloatType) * CHAR_BIT;
	constexpr int shift = numBits - 1;
	constexpr UintType mask = ~(static_cast<UintType>(1) << shift);
	UintType bitValues = reinterpret_cast<UintType&>(value) & mask;

	return reinterpret_cast<FloatType&>(bitValues);
}

#if DO_TEST
int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages)
{
	int errorCount = 0;

	using namespace std;

	{
		cout << "[bitops][TC" << ++inOutTestCount << "] Print bits" << endl;

		auto printBits = [](auto value)
		{
			constexpr auto numBits = sizeof(value) * CHAR_BIT;
			cout << value << "(" << numBits << ") = " << getBitsStrings(value) <<endl;
		};

		printBits(-1);
		printBits(1);
		printBits(0 >> 1);
		printBits(1 >> 1);
		printBits(1 << 31);
		printBits((1 << 31) >> 1);
		printBits((-1 << 1) >> 1);

		printBits(abs(static_cast<int8_t>(-1)));
		printBits(abs(static_cast<int16_t>(-1)));
		printBits(abs(static_cast<int32_t>(-1)));
		printBits(abs(static_cast<int64_t>(-1L)));
		
		printBits(0.0f);
		printBits(-0.0f);
		printBits(0.0);
		printBits(-0.0);
		printBits(1.0f);
		printBits(2.0f);
		printBits(4.0f);
		printBits(8.0f);
		printBits(16.0f);
		printBits(10.0f);
		printBits(100.0f);
		printBits(1000.0f);
		printBits(-1.0f);
		
		printBits(-2.0f);
		printBits(1.0);
		printBits(-1.0);
	}

	{
		cout << "[bitops][TC" << ++inOutTestCount << "] Absolute Value Tests" << endl;

		{
			auto value = -17;
			auto trueValue = ::abs(value);
			auto absValue = bitops::abs(value);

			cout << "[bitops][TC" << ++inOutTestCount << "] absolute value of " << value << " = " << absValue << endl;
			if (absValue != trueValue)
			{
				cout << "[bitops][TC" << ++inOutTestCount << "][Error] Incorrect absolute value "
					<< absValue << ", " << trueValue << " is expected." << endl;
				++errorCount;
			}
		}
	
		{
			auto value = -19L;
			auto trueValue = ::abs(value);
			auto absValue = bitops::abs(value);

			cout << "[bitops][TC" << ++inOutTestCount << "] absolute value of " << value << " = " << absValue << endl;
			if (absValue != trueValue)
			{
				cout << "[bitops][TC" << ++inOutTestCount << "][Error] Incorrect absolute value "
					<< absValue << ", " << trueValue << " is expected." << endl;
				++errorCount;
			}
		}

		{
			auto value = -3.14f;
			auto trueValue = ::abs(value);
			auto absValue = bitops::abs(value);

			cout << "[bitops][TC" << ++inOutTestCount << "] absolute value of " << value << " = " << absValue << endl;
			if (absValue != trueValue)
			{
				cout << "[bitops][TC" << ++inOutTestCount << "][Error] Incorrect absolute value "
					<< absValue << ", " << trueValue << " is expected." << endl;
				++errorCount;
			}
		}

		{
			auto value = -3.14;
			auto trueValue = ::abs(value);
			auto absValue = bitops::abs(value);

			cout << "[bitops][TC" << ++inOutTestCount << "] absolute value of " << value << " = " << absValue << endl;
			if (absValue != trueValue)
			{
				cout << "[bitops][TC" << ++inOutTestCount << "][Error] Incorrect absolute value "
					<< absValue << ", " << trueValue << " is expected." << endl;
				++errorCount;
			}
		}
	}
	
	return errorCount;
}
#endif // DO_TEST
} // bitops
} // hmath
