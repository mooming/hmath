#pragma once

#include "hmathconfig.h"

#include <cstdint>
#include <functional>

namespace hmath
{
#if USE_HIGH_PRECISION
	using HReal = double;
#else // USE_HIGH_PRECIOSION
	using HReal = float;
#endif // USE_HIGH_PRECISION

	// Standard function with a single parameter and having the same domain and range.
	// f:x -> y, where x and y are real numbers.
	using TFunc1 = std::function<HReal(HReal)>;
}
