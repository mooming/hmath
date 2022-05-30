#include "hmathpolynomial.h"

#include "hmathanalysis.h"
#include "hmathutil.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>


namespace hmath
{
    Polynomial::Polynomial(std::initializer_list<HReal> inCoefficients)
        : coefficients(inCoefficients)
    {
    }

    Polynomial::Polynomial(const std::vector<HReal>& inCoefficients)
        : coefficients(inCoefficients)
    {
    }

	Polynomial::Polynomial(std::vector<HReal>&& inCoefficients)
        : coefficients(std::move(inCoefficients))
    {
    }

    Polynomial::Polynomial(TFunc1 smoothFunc, HReal point, int depth, HReal epsilon)
    {
        // Taylor Series at the given point
        auto y = smoothFunc;
        auto dy = analysis::getDerivative(y, epsilon);

        coefficients.reserve(depth);
        
        for (int i = 0; i < depth; ++i)
        {
            
        }
    }

    Polynomial Polynomial::operator+ (const Polynomial& rhs) const
    {
        TOrder sizeDiff = numCoefficients() - rhs.numCoefficients();

        const Polynomial* bigger = nullptr;
        const Polynomial* smaller = nullptr;

        if (sizeDiff < 0)
        {
            bigger = &rhs;
            smaller = this;
            sizeDiff = -sizeDiff;
        }
        else
        {
            bigger = this;
            smaller = &rhs;
        }

        assert(bigger != nullptr);
        assert(smaller != nullptr);

        const auto bigSize = bigger->numCoefficients();
        const auto smallSize = smaller->numCoefficients();
        
        Polynomial outcome;
        auto& outCoeffs = outcome.coefficients;
        outCoeffs.reserve(bigSize);

        TOrder index = 0;
        for (index = 0; index < sizeDiff; ++index)
        {
            outCoeffs.push_back(bigger->coefficients[index]);
        }

        for (TOrder i = 0; i < smallSize; ++i, ++index)
        {
            outCoeffs.push_back(bigger->coefficients[index] + smaller->coefficients[i]);
        }

        return outcome;
    }

    Polynomial Polynomial::operator- (const Polynomial& rhs) const
    {
        TOrder sizeDiff = numCoefficients() - rhs.numCoefficients();

        const Polynomial* bigger = nullptr;
        const Polynomial* smaller = nullptr;

        TOrder selfIndex = 0;
        TOrder rhsIndex = 0;

        if (sizeDiff < 0)
        {
            bigger = &rhs;
            smaller = this;
            sizeDiff = -sizeDiff;
            rhsIndex = sizeDiff;
        }
        else
        {
            bigger = this;
            smaller = &rhs;
            selfIndex = sizeDiff;
        }

        assert(bigger != nullptr);
        assert(smaller != nullptr);

        const auto bigSize = bigger->numCoefficients();
        const auto smallSize = smaller->numCoefficients();
        
        Polynomial outcome;
        auto& outCoeffs = outcome.coefficients;
        outCoeffs.reserve(bigSize);

        TOrder index = 0;
        for (index = 0; index < sizeDiff; ++index)
        {
            outCoeffs.push_back(bigger->coefficients[index]);
        }

        for (TOrder i = 0; i < smallSize; ++i, ++selfIndex, ++rhsIndex)
        {
            outCoeffs.push_back(getCoefficient(selfIndex) - rhs.getCoefficient(rhsIndex));
        }

        return outcome;
    }

    Polynomial Polynomial::operator* (const Polynomial& rhs) const
    {
        Polynomial result;

        TOrder sizeRhs = rhs.numCoefficients();
        TOrder orderRhs = sizeRhs - 1;

        for (TOrder i = 0; i < sizeRhs; ++i)
        {
            const TOrder exponent = orderRhs - i;
            assert(exponent >= 0);

            auto coeff = rhs.getCoefficient(i);
            Polynomial tmp(*this);
            tmp *= coeff;
            tmp.shiftUp(exponent);

            result = result + tmp;
        }

        return result;
    }
    
    Polynomial Polynomial::operator* (HReal value) const
    {
        Polynomial result(*this);
        result *= value;

        return result;
    }

    void Polynomial::operator*= (HReal value)
    {
        for (auto& coeff : coefficients)
        {
            coeff *= value;
        }
    }

    TFunc1 Polynomial::AsFunction() const
    {
        auto func = [*this](HReal value)
        {
            return evaluate(value);
        };

        return func;
    }

    Polynomial::TOrder Polynomial::numCoefficients() const
    {
        const auto order = static_cast<TOrder>(coefficients.size());
        return order;
    }

    Polynomial::TOrder Polynomial::getOrder() const
    {
        const auto size = static_cast<TOrder>(coefficients.size());

        return size - 1;
    }

    HReal Polynomial::getCoefficient(TOrder index) const
    {
        if (index < 0 || index >= numCoefficients())
            return ZERO;

        return coefficients.at(index);
    }

    HReal Polynomial::evaluate(HReal value) const
    {
        // Horner's method
        HReal y = ZERO;
        
        for (auto coeff : coefficients)
        {
            y = y * value + coeff;
        }

        return y;
    }

    void Polynomial::shiftUp(unsigned int numShift)
    {
        if (numShift == 0)
            return;
            
        std::vector<HReal> tmp;
        tmp.reserve(coefficients.size() + numShift);

        for (auto coeff : coefficients)
        {
            tmp.push_back(coeff);
        }

        for (unsigned int i = 0; i < numShift; ++i)
        {
            tmp.push_back(ZERO);
        }

        std::swap(tmp, coefficients);
    }
	
    void Polynomial::shiftDown(unsigned int numShift)
    {
        if (numShift == 0)
            return;

        std::vector<HReal> tmp;

        auto newSize = coefficients.size() - numShift;
        tmp.reserve(newSize);

        for (size_t i = 0; i < newSize; ++i)
        {
            tmp.push_back(coefficients[i]);
        }

        std::swap(tmp, coefficients);
    }

    void Polynomial::defferentiate()
    {
        auto order = getOrder();
        for (TOrder i = 0; i < order; ++i)
        {
            TOrder exponent = (order - i);
            coefficients[i] *= exponent;
        }

        coefficients.pop_back();
    }

    void Polynomial::integrate(HReal constant)
    {
        auto size = numCoefficients();

        for (TOrder i = 0; i < size; ++i)
        {
            TOrder newExponent = (size - i);
            coefficients[i] /= newExponent;
        }

        coefficients.push_back(constant);
    }

    void Polynomial::print() const
    {
        using namespace std;

        cout << "y = ";

        const TOrder lastIndex = getOrder();
        for (TOrder i = 0; i < lastIndex; ++i)
        {
            cout << coefficients[i] << "x^" << (lastIndex - i) << " + ";
        }

        cout << coefficients[lastIndex] << endl;
    }

    void Polynomial::print(HReal value) const
    {
        using namespace std;

        const TOrder lastIndex = getOrder();
        for (TOrder i = 0; i < lastIndex; ++i)
        {
            cout << coefficients[i] << "x^" << (lastIndex - i) << " + ";
        }

        cout << coefficients[lastIndex] << " = " << evaluate(value) << endl;
    }

    std::ostream& operator<< (std::ostream& stream, const Polynomial& polynomial)
    {
        using TOrder = Polynomial::TOrder;
        
        stream << "y = ";

        const TOrder lastIndex = polynomial.getOrder();
        if (lastIndex < 0)
            return stream;

        for (TOrder i = 0; i < lastIndex; ++i)
        {
            stream << polynomial.coefficients[i] << "x^" << (lastIndex - i) << " + ";
        }

        stream << polynomial.coefficients[lastIndex];

        return stream;
    }

#if DO_TEST
	int Polynomial::DoTest(int& inOutTestCount, std::vector<std::string>& outErrorMessages)
    {
        using namespace std;

        int errorCount = 0;

        cout << endl << "[Polynomial] TestCase " << ++inOutTestCount << ") Simple Polynomial" << endl;
        {
            auto func = [](HReal x) { return x * x * x + 2 * x * x - 3 * x + 1; };
            Polynomial p({1, 2, -3, 1});
            p.print();

            auto error = util::compare(p.AsFunction(), func, -10, 10, 0.01);
            cout << "[Polynomial][TC" << inOutTestCount
                << "] P value check: PASS, error = " << error << endl;

            if (error > EPSILON)
            {
                ++errorCount;

                ostringstream msg;
                msg << "[Polynomial][TC" << inOutTestCount << "][Error] error " << error
                    << " is bigger than expected " << SMALL_NUMBER << endl;

                auto errorMsg = msg.view();
                cerr << errorMsg;

                outErrorMessages.emplace_back(errorMsg);
            }
        }

        cout << endl << "[Polynomial] TestCase " << ++inOutTestCount << ") Add two polynomials" << endl;
        {
            auto func = [](HReal x) { return x * x * x + 2 * x * x + x + 2; };
            
            Polynomial p1({ 1, 0, 1, 0 });
            Polynomial p2({ 2, 0, 2 });
            Polynomial p3 = p1 + p2;

            p3.print();

            auto error = util::compare(p3.AsFunction(), func, -10, 10, 0.01);
            cout << "[Polynomial][TC" << inOutTestCount
                << "] P value check: PASS, error = " << error << endl;

            if (error > EPSILON)
            {
                ++errorCount;

                ostringstream msg;
                msg << "[Polynomial][TC" << inOutTestCount << "][Error] error " << error
                    << " is bigger than expected " << SMALL_NUMBER << endl;

                auto errorMsg = msg.view();
                cerr << errorMsg;

                outErrorMessages.emplace_back(errorMsg);
            }
        }

        cout << endl << "[Polynomial] TestCase " << ++inOutTestCount << ") Subtract two polynomials" << endl;
        {
            auto func = [](HReal x) { return x * x * x - 2 * x * x + x - 2; };

            Polynomial p1({ 1, 0, 1, 0 });
            Polynomial p2({ 2, 0, 2 });
            Polynomial p3 = p1 - p2;

            p3.print();

            auto error = util::compare(p3.AsFunction(), func, -10, 10, 0.01);
            cout << "[Polynomial][TC" << inOutTestCount
                << "] P value check: PASS, error = " << error << endl;

            if (error > EPSILON)
            {
                ++errorCount;

                ostringstream msg;
                msg << "[Polynomial][TC" << inOutTestCount << "][Error] error " << error
                    << " is bigger than expected " << SMALL_NUMBER << endl;

                auto errorMsg = msg.view();
                cerr << errorMsg;

                outErrorMessages.emplace_back(errorMsg);
            }
        }

        cout << endl << "[Polynomial] TestCase " << ++inOutTestCount << ") Multiply two polynomials" << endl;
        {
            Polynomial p1({ 1, 1});
            Polynomial p2({ 1, 1 });
            Polynomial p3 = p1 * p2;
            Polynomial answer({1, 2, 1});

            cout << "[Polynomial][TC" << inOutTestCount << "] (" << p1 << ") X (" << p2 << ") = (" << p3 << ')' << endl; 

            if (p3 != answer)
            {
                ++errorCount;

                ostringstream msg;
                msg << "[Polynomial][TC" << inOutTestCount << "][Error] error " << p3
                    << " doesn't coincide with " <<  answer << endl;

                auto errorMsg = msg.view();
                cerr << errorMsg;

                outErrorMessages.emplace_back(errorMsg);
            }
        }

        cout << endl << "[Polynomial] TestCase " << ++inOutTestCount << ") Multiply two polynomials (2)" << endl;
        {
            Polynomial p1({ 1, 2});
            Polynomial p2({ 2, 1 });
            Polynomial p3 = p1 * p2;
            Polynomial answer({2, 5, 2});

            cout << "[Polynomial][TC" << inOutTestCount << "] (" << p1 << ") X (" << p2 << ") = (" << p3 << ')' << endl; 

            if (p3 != answer)
            {
                ++errorCount;

                ostringstream msg;
                msg << "[Polynomial][TC" << inOutTestCount << "][Error] error " << p3
                    << " doesn't coincide with " <<  answer << endl;

                auto errorMsg = msg.view();
                cerr << errorMsg;

                outErrorMessages.emplace_back(errorMsg);
            }
        }

        cout << endl << "[Polynomial] TestCase " << ++inOutTestCount << ") Multiply three polynomials" << endl;
        {
            Polynomial p1({ 1, 1});
            Polynomial p2 = p1 * p1 * p1;
            Polynomial answer({1, 3, 3, 1});

            cout << "[Polynomial][TC" << inOutTestCount << "] (" << p1 << ") X (" << p1 << ") X (" << p1 << ") = (" << p2 << ')' << endl; 

            if (p2 != answer)
            {
                ++errorCount;

                ostringstream msg;
                msg << "[Polynomial][TC" << inOutTestCount << "][Error] error " << p2
                    << " doesn't coincide with " <<  answer << endl;

                auto errorMsg = msg.view();
                cerr << errorMsg;

                outErrorMessages.emplace_back(errorMsg);
            }
        }

        cout << endl << "[Polynomial] TestCase " << ++inOutTestCount << ") Differentiation" << endl;
        {

            Polynomial p1({ 3, 2, 1});
            
            Polynomial dp = p1;
            dp.defferentiate();

            Polynomial answer({6, 2});

            cout << "[Polynomial][TC" << inOutTestCount << "] (" << p1 << ")` = (" << dp << ')' << endl; 

            if (dp != answer)
            {
                ++errorCount;

                ostringstream msg;
                msg << "[Polynomial][TC" << inOutTestCount << "][Error] error " << dp
                    << " doesn't coincide with " <<  answer << endl;

                auto errorMsg = msg.view();
                cerr << errorMsg;

                outErrorMessages.emplace_back(errorMsg);
            }
        }

        cout << endl << "[Polynomial] TestCase " << ++inOutTestCount << ") Integration" << endl;
        {

            Polynomial p1({ 6, 2});
            
            Polynomial integral = p1;
            integral.integrate(ONE);

            Polynomial answer({3, 2, 1});

            cout << "[Polynomial][TC" << inOutTestCount << "] Integral of (" << p1 << ") = (" << integral << ')' << endl; 

            if (integral != answer)
            {
                ++errorCount;

                ostringstream msg;
                msg << "[Polynomial][TC" << inOutTestCount << "][Error] error " << integral
                    << " doesn't coincide with " <<  answer << endl;

                auto errorMsg = msg.view();
                cerr << errorMsg;

                outErrorMessages.emplace_back(errorMsg);
            }
        }

        return errorCount;
    }
#endif // DO_TEST
}
