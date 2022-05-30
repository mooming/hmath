#include "hmathtypes.h"

#include "hmathconstants.h"


namespace hmath
{
	HRoot::HRoot()
		: value(ZERO), error(ZERO)
	{
	}

	HRoot::HRoot(HReal inValue, HReal inError)
		: value(inValue), error(inError)
	{
	}
}