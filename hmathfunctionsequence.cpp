#include "hmathfunctionsequence.h"

#include <algorithm>
#include <iostream>


namespace hmath
{

FunctionSequence::FunctionSequence(std::initializer_list<TFunc1> list)
{
	functions.reserve(list.size());

	for (auto func : list)
	{
		if (!func)
		{
			using namespace std;
			cerr << "[FunctionSequnece][Error] ctr: Null function found!" << endl;
			continue;
		}

		functions.push_back(func);
	}
}

void FunctionSequence::add(TFunc1 func)
{
	if (!func)
	{
		using namespace std;
		cerr << "[FunctionSequnece][Error] add: Null function found!" << endl;
		return;
	}

	functions.push_back(func);
}

void FunctionSequence::clear()
{
	functions.clear();
}

void FunctionSequence::empty()
{
	std::vector<TFunc1>().swap(functions);
}

TFunc1 FunctionSequence::asFunction() const
{
	return[*this](HReal value)->HReal
	{
		return get(value);
	};
}

HReal FunctionSequence::get(HReal x) const
{
	HReal result = x;

	for (auto func : functions)
	{
		result = func(result);
	}

	return result;
}

HReal FunctionSequence::getSub(HReal x, int start, int end) const
{
	start = truncateIndex(start);
	end = truncateIndex(end);

	HReal result = x;

	for (int i = start; i < end; ++i)
	{
		result = functions[i](result);
	}

	return result;
}

int FunctionSequence::truncateIndex(int index) const
{
	return std::clamp(index, 0, static_cast<int>(functions.size()));
}

#if DO_TEST
int FunctionSequence::DoTest(int& inOutTestCount,
	std::vector<std::string>& outErrorMessages)
{
	using namespace std;

	int errorCount = 0;

	{
		cout << "[hmath][FunctionSequence][TC" << ++inOutTestCount << "] Basic tests " << endl;

		auto squareFunc = [](HReal x) -> HReal { return x * x; };
		FunctionSequence fSeqs;

		for (int i = 0; i < 20; ++i)
			fSeqs.add(squareFunc);
	}

	return errorCount;
}

#endif // DO_TEST

} // hmath