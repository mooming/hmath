#pragma once

#include "hmathconfig.h"
#include "hmathconstants.h"
#include "hmathtypes.h"

#include <initializer_list>
#include <vector>
#include <ostream>


namespace hmath
{
	class Polynomial final
	{
		using TOrder = int;

	private:
		std::vector<HReal> coefficients;

	public:
		Polynomial() = default;
		Polynomial(std::initializer_list<HReal> inCoefficients);
		explicit Polynomial(const std::vector<HReal>& inCoefficients);
		explicit Polynomial(std::vector<HReal>&& inCoefficients);
		~Polynomial() = default;

		Polynomial operator+ (const Polynomial& rhs) const;
		Polynomial operator- (const Polynomial& rhs) const;
		Polynomial operator* (const Polynomial& rhs) const;
		Polynomial operator* (HReal value) const;
		void operator*= (HReal value);	

		inline bool operator== (const Polynomial& rhs) const { return coefficients == rhs.coefficients; }
		inline bool operator!= (const Polynomial& rhs) const { return coefficients != rhs.coefficients; }

		friend std::ostream& operator<< (std::ostream& stream, const Polynomial& polynomial);

	public:
		TFunc1 AsFunction() const;

		TOrder numCoefficients() const;
		TOrder getOrder() const;
		HReal getCoefficient(TOrder index) const;
		HReal evaluate(HReal value) const;

		void shiftUp(unsigned int numShift);
		void shiftDown(unsigned int numShift);
		void defferentiate();
		void integrate(HReal constant = ZERO);

		void print() const;
		void print(HReal value) const;

#if DO_TEST
		static int DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages);
#endif // DO_TEST
	};
}

