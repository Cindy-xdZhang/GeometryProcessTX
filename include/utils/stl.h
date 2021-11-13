#ifndef UTILS_STL_H
#define UTILS_STL_H

#include <vector>
#include <cassert>
#include <iterator>
#include <algorithm>


////////////////////////////////////////////////////////////////////////////////
/// STL extensions
////////////////////////////////////////////////////////////////////////////////

namespace std
{
    template<typename T1, typename T2, typename T3>
    struct triplet
    {
        T1 first;
        T2 second;
        T3 third;

        triplet()
        {
            first = T1();
            second = T2();
            third = T3();
        }

        triplet(const T1& m1, const T2& m2, const T3& m3)
        {
            first = m1;
            second = m2;
            third = m3;
        }
    };

    template<typename T1, typename T2, typename T3>
    triplet<T1, T2, T3> make_triplet(const T1& m1, const T2& m2, const T3& m3)
    {
        triplet<T1, T2, T3> ans;
        ans.first = m1;
        ans.second = m2;
        ans.third = m3;
        return ans;
    }

    template<typename T1, typename T2, typename T3, typename T4>
    struct quadruplet
    {
        T1 first;
        T2 second;
        T3 third;
        T4 fourth;

        quadruplet()
        {
            first = T1();
            second = T2();
            third = T3();
            fourth = T4();
        }

        quadruplet(const T1& m1, const T2& m2, const T3& m3, const T4& m4)
        {
            first = m1;
            second = m2;
            third = m3;
            fourth = m4;
        }
    };

    template<typename T1, typename T2, typename T3, typename T4>
    quadruplet<T1, T2, T3, T4> make_quadruplet(const T1& m1, const T2& m2, const T3& m3, const T4& m4)
    {
        quadruplet<T1, T2, T3, T4> ans;
        ans.first = m1;
        ans.second = m2;
        ans.third = m3;
        ans.fourth = m4;
        return ans;
    }
}


////////////////////////////////////////////////////////////////////////////////
/// C++ stl utils
////////////////////////////////////////////////////////////////////////////////

namespace utils
{
    inline std::vector<int> intersection(std::vector<int>& v1, std::vector<int>& v2)
    {
        std::vector<int> v3;

        std::sort(v1.begin(), v1.end());
        std::sort(v2.begin(), v2.end());

        std::set_intersection(v1.begin(), v1.end(),
                              v2.begin(), v2.end(),
                              back_inserter(v3));
        return v3;
    }

    inline std::vector<std::pair<int, int>> zip(const std::vector<int>& v1, const std::vector<int>& v2)
    {
        assert(v1.size() == v2.size());

        std::vector<std::pair<int, int>> v3;
        v3.reserve(v1.size());

        std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(v3),
                       [](int a, int b)
                       {
                           return std::make_pair(a, b);
                       });

        return v3;
    }

    inline std::vector<std::pair<int, int>> zip(int v, const std::vector<int>& v2)
    {
        std::vector<int> v1(v2.size(), v);
        return utils::zip(v1, v2);
    }
}

#endif // UTILS_STL_H
