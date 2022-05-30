#pragma once

#include "hmathconfig.h"
#include "hmathtypes.h"

#include <vector>


namespace hmath
{
	class FunctionSequence
	{
	private:
		std::vector<TFunc1> functions;

	public:
		FunctionSequence() = default;
		FunctionSequence(std::initializer_list<TFunc1> list);
		~FunctionSequence() = default;

		void add(TFunc1 func);
		void clear();
		void empty();

		TFunc1 asFunction() const;
		HReal get(HReal x) const;
		HReal getSub(HReal x, int start, int end) const;

#if DO_TEST
		static int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages);
#endif // DO_TEST

	private:
		int truncateIndex(int index) const;
	};
} // hmath
