#include <algorithm>
#include "BezierCurve.h"

namespace
{
    std::vector<size_t> g_primeNumbers;
}

const std::vector<size_t>& PrimeNumbers()
{
    return g_primeNumbers;
}

bool IsPrime(size_t number)
{
    const size_t root = static_cast<size_t>(std::sqrt(number));
    for (size_t i = 2; i <= root; ++i)
    {
        if (number %i == 0)
        {
            return false;
        }
    }

    return true;
}

bool IsPrimeCached(size_t number)
{
    //n is highest checked number
    static size_t n = 1;

    if (number > n)
    {
        for (size_t i = n + 1; i <= number; ++i)
        {
            if (IsPrime(i))
            {
                g_primeNumbers.push_back(i);
            }
        }

        n = number;
    }

    return std::binary_search(g_primeNumbers.begin(), g_primeNumbers.end(), number);
}

void GetNumberPrimeFactors(size_t number, std::vector<std::pair<size_t, size_t>>* result)
{
    assert(result != nullptr);

    result->clear();

    if (IsPrimeCached(number))
    {
        result->emplace_back(number, 1);
        return;
    }

    const size_t half = static_cast<size_t>(number / 2);

    for (size_t prime : PrimeNumbers())
    {
        std::pair<size_t, size_t> numberAndPower(prime, 0);

        while (number % prime == 0)
        {
            number /= prime;
            ++numberAndPower.second;
        }

        if (numberAndPower.second > 0)
        {
            result->push_back(numberAndPower);
        }

        if (prime > half)
        {
            assert(number == 1);
            return;
        }
    }
}

/* The number of combinations from n by i */
size_t Combinations(size_t n, size_t i)
{
    std::vector<std::pair<size_t, size_t>> primeFactors;

    std::unordered_map<size_t, size_t> numerator;
    std::unordered_map<size_t, size_t> denominator;

    auto factorial = [&primeFactors](size_t n, std::unordered_map<size_t, size_t>& c)
    {
        for (size_t i = 2; i <= n; ++i)
        {
            GetNumberPrimeFactors(i, &primeFactors);

            for (auto& factorAndPow : primeFactors)
            {
                const size_t factor = factorAndPow.first;
                const size_t pow = factorAndPow.second;

                auto iVal = c.find(factor);
                if (iVal == c.end())
                {
                    c.insert(factorAndPow);
                }
                else
                {
                    iVal->second += pow;
                }
            }
        }
    };

    factorial(n, numerator);
    factorial(i, denominator);
    factorial(n - i, denominator);

    struct
    {
        size_t numerator;
        size_t denominator;
    } result{ 1, 1 };

    for (auto& numeratorFactor : numerator)
    {
        auto iDenominatorFactor = denominator.find(numeratorFactor.first);

        if (iDenominatorFactor != denominator.end())
        {
            /* Result is always integral number */
            assert(numeratorFactor.second >= iDenominatorFactor->second);

            /* Reduce */
            numeratorFactor.second -= iDenominatorFactor->second;
            iDenominatorFactor->second = 0;
        }

        result.numerator *= std::pow(numeratorFactor.first, numeratorFactor.second);
    }

    return result.numerator;

    // n! / (i!(n-i)!)
    //return Factorial(n) / (Factorial(i) * Factorial(n - i));
}