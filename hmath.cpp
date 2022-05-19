
#include "hmathanalysis.h"
#include "hmathconfig.h"
#include "hmathconstants.h"
#include "hmathpolynomial.h"

#include <iostream>
#include <limits>
#include <string>
#include <vector>

#if DO_TEST

int main()
{
	using namespace std;
	using namespace hmath;

	int testCount = 0;
	int errorCount = 0;
	vector<string> errorMessages;
	
	cout << "[HMath] Test Started! ===" << endl;

	errorCount += analysis::DoTest(testCount, errorMessages);
	errorCount += Polynomial::DoTest(testCount, errorMessages);

	cout << endl;
	cout << "[HMath] Test Finished! ===" << endl;
	cout << "[HMath] Test Report" << endl;
	cout << "[HMath] Errors " << errorCount << endl;
	cout << endl;
	
	return 0;
}

#endif // DO_TEST
