#pragma once

#include "hmathtypes.h"

#include <limits>


namespace hmath
{
	// Epsilon for derivative and basic math functions
	static constexpr HReal EPSILON = 1.0e-6;
	static constexpr HReal MACHINE_EPSILON = std::numeric_limits<HReal>::epsilon();

	static constexpr HReal ZERO = 0;
	static constexpr HReal ONE = 1;
	static constexpr HReal ONE_TENTH = 0.1;
	static constexpr HReal MINUS_ONE = -1;
	static constexpr HReal TWO = 2;
	static constexpr HReal TEN = 10;

	static constexpr HReal HALF = 0.5;
	static constexpr HReal CENTI = 1.0e-2;
	static constexpr HReal MILI = 1.0e-3;
	static constexpr HReal MICRO = 1.0e-6;
	static constexpr HReal NANO = 1.0e-9;

	static constexpr HReal PI = 3.141592653589793;
	static constexpr HReal TWO_PI = PI * TWO;

	// Epsilon for physical calculation in KMS (represents 0.1 mm)
	static constexpr HReal SMALL_NUMBER = 1.0e-4;
	static constexpr HReal MIN_NUMBER = std::numeric_limits<HReal>::min();
	static constexpr HReal MAX_NUMBER = std::numeric_limits<HReal>::max();


	constexpr HReal degreesToRadians(HReal value)
	{
		return value / 180 * PI;
	};

	constexpr HReal radiansToDegree(HReal value)
	{
		return value / PI * 180;
	};
}
