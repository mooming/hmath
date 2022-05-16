#pragma once

#include "hmathconfig.h"

#include <cstdint>


namespace hmath
{
#if USE_HIGH_PRECISION
		using HReal = double;
#else // USE_HIGH_PRECIOSION
		using HReal = float;
#endif // USE_HIGH_PRECISION
}
