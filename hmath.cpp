
#include "hmath.h"

#include "hmathanalysis.h"
#include "hmathbitops.h"
#include "hmathconfig.h"
#include "hmathconstants.h"
#include "hmathfunctionsequence.h"
#include "hmathpolynomial.h"
#include "hmathutil.h"

#include <algorithm>
#include <iostream>


namespace hmath
{

#if DO_TEST
bool DoTest()
{
	using namespace std;
	
	int testCount = 0;
	int errorCount = 0;
	vector<string> errorMessages;

	cout << "[HMath] Test Started! ===" << endl;

	errorCount += bitops::DoTest(testCount, errorMessages);
	errorCount += util::DoTest(testCount, errorMessages);
	errorCount += Polynomial::DoTest(testCount, errorMessages);
	errorCount += analysis::DoTest(testCount, errorMessages);

	cout << endl;
	cout << "[HMath] Test Finished! ===" << endl;
	cout << "[HMath] Test Report" << endl;
	cout << "[HMath] Errors " << errorCount << endl;
	cout << endl;

	int errorIndex = 0;
	for (auto& msg : errorMessages)
	{
		cerr << "[Error" << errorIndex++ << "] " << msg;
	}

	return errorCount == 0;
}
#endif // DO_TEST

} // hmath