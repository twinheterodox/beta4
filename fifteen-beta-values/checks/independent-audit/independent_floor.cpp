#include <algorithm>
#include <iostream>
#include <vector>

using Integer = long long;

Integer floor_divide(Integer numerator, Integer denominator) {
    Integer quotient = numerator / denominator;
    if (numerator % denominator < 0) --quotient;
    return quotient;
}

Integer residue(Integer numerator, Integer denominator) {
    return numerator - denominator * floor_divide(numerator, denominator);
}

int main() {
    Integer a, b;
    while (std::cin >> a >> b) {
        const Integer denominator = 96 * b;
        const Integer cap = a < 2 * b ? 96 * (a-b) : denominator;
        std::vector<Integer> points{0, cap};
        for (int j=0; j<=30; ++j) {
            points.push_back(residue(j*a, denominator));
            points.push_back(residue((96-j)*a, denominator));
        }
        points.push_back(residue(-72*a-48*b, denominator));
        points.push_back(residue(168*a-48*b, denominator));
        std::sort(points.begin(), points.end());
        points.erase(std::unique(points.begin(), points.end()), points.end());
        points.erase(std::upper_bound(points.begin(), points.end(), cap), points.end());
        auto value = [&](Integer y) {
            int count = 2;
            for (int j=0; j<=30; ++j) {
                count += floor_divide(y-2*j*a, 2*denominator)
                       + floor_divide(2*(96-j)*a-y, 2*denominator) + 1;
            }
            count -= 3*(floor_divide(y+144*a+96*b, 2*denominator)
                      - floor_divide(y-336*a+96*b, 2*denominator));
            return count;
        };
        int maximum = 0;
        for (std::size_t j=0; j<points.size(); ++j) {
            maximum = std::max(maximum, value(2*points[j]));
            if (j+1<points.size())
                maximum = std::max(maximum, value(points[j]+points[j+1]));
        }
        std::cout << 16*maximum-221*(a/b) << '\n';
    }
}
