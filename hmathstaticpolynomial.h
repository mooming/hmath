
#include "hmathtypes.h"

#include <initializer_list>


namespace hmath
{
	template <int Order = 1>
	class StaticPolynomial final
	{
		static_assert(Order >= 0);

	private:
		HReal coefficients[Order];

	public:
		StaticPolynomial() = default;
		StaticPolynomial(std::initializer_list inCoefficients)
			: coefficients(inCoefficients)
		{
		}

		~StaticPolynomial() = default;

		bool operator== (const StaticPolynomial& rhs) const
		{
			for (int i = 0; i < Order; ++i)
			{
				if (coefficients[i] != rhs.coefficients[i])
					return false;
			}

			return true;
		}

		bool operator!= (const StaticPolynomial& rhs) const
		{
			return !(*this == rhs);
		}

	public:
		constexpr int getOrder() const { return Order; }

		HReal getCoefficient(int index) const
		{
			if (index < 0 || index >= getOrder())
				return ZERO;

			return coefficients.at(index);
		}

		HReal evaluate(HReal value) const
		{
			HReal y = ZERO;

			for (auto coeff : coefficients)
			{
				y = y * value + coeff;
			}

			return y;
		}
	};
}