#ifndef LI_HPP
#define LI_HPP

#include<map>
#include<algorithm>
#include<functional>

// mapに基づいた線形補間クラス
template<typename X = double, typename Y = double>
class LI {
public:
    typedef std::map<X, Y> Points;

    LI(){}

    template<class Iter>
    LI(Iter first, Iter last) : points_(first, last) {}

    void append(X x, Y y){
        points_[x] = y;
    }

    Y operator[] (X x) const{
        assert(!points_.empty());
        //xがinsert済みのものと一致していればその値を返し，そうでなければ補間計算の結果を返す
        return points_.find(x) == points_.end() ? interpolate_(x) : points_.at(x);
    }

private:
    static bool greater_(typename Points::value_type p, X x){
        return p.first > x ? true : false;
    }

    Y interpolate_(X x) const{
        //xを超える最小のpointを指すイテレータ
        auto it_h = std::find_if(points_.begin(), points_.end(), 
                                 std::bind2nd(std::ptr_fun(greater_), x));

        //補外が必要になるようなケースはassertで弾く
        assert(it_h != points_.begin() && it_h != points_.end());

        //xを下回る最大のpointを指すイテレータ
        auto it_l = it_h;
        it_l--;

        return linear_((*it_l).first, (*it_l).second, (*it_h).first, (*it_h).second, x);
    }

    Y linear_(X x1, Y y1, X x2, Y y2, X x) const{
        return (x - x1) * (y2 - y1) / (x2 - x1) + y1;
    }

    Points points_;
};

#endif