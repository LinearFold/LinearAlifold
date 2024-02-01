#include <cassert>
#define VALUE_MIN std::numeric_limits<double>::lowest()

typedef double value_type;

// log space: borrowed from CONTRAfold
inline value_type Fast_LogExpPlusOne(value_type x) {

    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    // Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    // Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    // Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    // Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    // Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    // 8 polynomials needed.

    assert(value_type(0.0000000000) <= x && x <= value_type(11.8624794162) && "Argument out-of-range.");
    if (x < value_type(3.3792499610)) {
        if (x < value_type(1.6320158198)) {
            if (x < value_type(0.6615367791))
                return ((value_type(-0.0065591595) * x + value_type(0.1276442762)) * x + value_type(0.4996554598)) * x + value_type(0.6931542306);
            return ((value_type(-0.0155157557) * x + value_type(0.1446775699)) * x + value_type(0.4882939746)) * x + value_type(0.6958092989);
        }
        if (x < value_type(2.4912588184))
            return ((value_type(-0.0128909247) * x + value_type(0.1301028251)) * x + value_type(0.5150398748)) * x + value_type(0.6795585882);
        return ((value_type(-0.0072142647) * x + value_type(0.0877540853)) * x + value_type(0.6208708362)) * x + value_type(0.5909675829);
    }
    if (x < value_type(5.7890710412)) {
        if (x < value_type(4.4261691294))
            return ((value_type(-0.0031455354) * x + value_type(0.0467229449)) * x + value_type(0.7592532310)) * x + value_type(0.4348794399);
        return ((value_type(-0.0010110698) * x + value_type(0.0185943421)) * x + value_type(0.8831730747)) * x + value_type(0.2523695427);
    }
    if (x < value_type(7.8162726752))
        return ((value_type(-0.0001962780) * x + value_type(0.0046084408)) * x + value_type(0.9634431978)) * x + value_type(0.0983148903);
    return ((value_type(-0.0000113994) * x + value_type(0.0003734731)) * x + value_type(0.9959107193)) * x + value_type(0.0149855051);
}

inline void Fast_LogPlusEquals(value_type &x, value_type y)

{
    // x = log(exp(x) + exp(y));
    // std::cout<<"x y in Fast_LogPlusEquals "<<x<<' '<<y<<std::endl;
    if (x < y)
        std::swap(x, y);
    if (y > value_type(VALUE_MIN / 2) && x - y < value_type(11.8624794162))
        x = Fast_LogExpPlusOne(x - y) + y;

    // cout << "x " << x << endl;
    // cout << "y " << y << endl;
    // cout<<"x "<<x<<' '<<log(exp(x)+exp(y))<<endl;
}

inline value_type Fast_Exp(value_type x) {
    // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    // Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    // Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    // Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    // Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    // Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    // Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    // 6 polynomials needed.

    if (x < value_type(-2.4915033807)) {
        if (x < value_type(-5.8622823336)) {
            if (x < value_type(-9.91152))
                return value_type(0);
            return ((value_type(0.0000803850) * x + value_type(0.0021627428)) * x + value_type(0.0194708555)) * x + value_type(0.0588080014);
        }
        if (x < value_type(-3.8396630909))
            return ((value_type(0.0013889414) * x + value_type(0.0244676474)) * x + value_type(0.1471290604)) * x + value_type(0.3042757740);
        return ((value_type(0.0072335607) * x + value_type(0.0906002677)) * x + value_type(0.3983111356)) * x + value_type(0.6245959221);
    }
    if (x < value_type(-0.6725053211)) {
        if (x < value_type(-1.4805375919))
            return ((value_type(0.0232410351) * x + value_type(0.2085645908)) * x + value_type(0.6906367911)) * x + value_type(0.8682322329);
        return ((value_type(0.0573782771) * x + value_type(0.3580258429)) * x + value_type(0.9121133217)) * x + value_type(0.9793091728);
    }
    if (x < value_type(0))
        return ((value_type(0.1199175927) * x + value_type(0.4815668234)) * x + value_type(0.9975991939)) * x + value_type(0.9999505077);
    return (x > value_type(46.052) ? value_type(1e20) : expf(x));
}