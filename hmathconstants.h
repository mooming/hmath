#pragma once

#include "hmathtypes.h"

#include <limits>


namespace hmath
{
	static constexpr HReal Epsilon = 1.0e-5;

	static constexpr HReal Zero = 0;
	static constexpr HReal One = 1;
	static constexpr HReal Two = 2;

	static constexpr HReal Half = 0.5;
	static constexpr HReal Centi = 1.0e-2;
	static constexpr HReal Mili = 1.0e-3;
	static constexpr HReal Micro = 1.0e-6;
	static constexpr HReal Nano = 1.0e-9;

	static constexpr HReal Pi = 3.141592653589793;
	static constexpr HReal TwoPi = Pi * Two;

	static constexpr HReal SmallNumber = 1.0e-5;
	static constexpr HReal SmallestNumber = std::numeric_limits<HReal>::min();
	static constexpr HReal BiggestNumber = std::numeric_limits<HReal>::max();


	constexpr HReal degreesToRadians(HReal value)
	{
		return value / 180 * Pi;
	};

	constexpr HReal radiansToDegree(HReal value)
	{
		return value / Pi * 180;
	};
}
