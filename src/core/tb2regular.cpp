/*----------------------------------------------------------------------
 *
 * Incremental regular constraint
 *
 */

#include "tb2regular.hpp"
#include "utils/tb2automaton.hpp"
#include "search/tb2clusters.hpp"

//--------------------------------------------------
// output operator for vector, pair
namespace std {

template <typename T> ostream& operator<<(ostream& os, vector<T> const& v)
{
    os << "v(sz=" << v.size() << ")[";
    bool first = true;
    for (auto&& t : v) {
        if (first)
            first = false;
        else
            os << ",";
        os << t;
    }
    os << "]";
    return os;
}

template <typename U, typename T>
ostream& operator<<(ostream& os, pair<U, T> const& p)
{
    return os << "p{" << p.first << "," << p.second << "}";
}

}

WRegular::WRegular(WCSP* wcsp, EnumeratedVariable** scope_in,
                                   int arity_in)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , lb(MIN_COST)
{
    for (int i = 0; i != arity_; ++i) {
        conflictWeights.push_back(0);
    }
}

WRegular::~WRegular()
{
}

std::ostream& WRegular::printstate(std::ostream& os)
{
    os << "lb = " << lb << " connected = " << connected()
       << " depth = " << Store::getDepth() << "\n";
    for (int i = 0; i != arity_; ++i) {
        auto *x = scope[i];
        if (connected(i))
            os << " * ";
        else
            os << "   ";
        os << "var " << i << " ";
        x->print(os);
        if (x->assigned() && connected(i))
            os << "*****";
        os << "\n";
    }
    return  os;
}

void WRegular::propagate()
{
 
}

void WRegular::assign(int idx)
{
    static const bool debug{false};

    if (debug)
        cout << "After assign of " << idx << " state = \n"
             << state{this} << "\n";
}


void WRegular::remove(int idx)
{
    static const bool debug{false};

    if (debug)
        cout << "In remove " << idx << " run = " << run << "\n";

}

void WRegular::increase(int idx)
{
}

void WRegular::decrease(int idx)
{
}

void WRegular::projectFromZero(int idx)
{
    static const bool debug{false};

    if (!connected(idx))
        return;

    if (debug)
        cout << "In projectFromZero state = " << state{this} << "\n";

}


double WRegular::computeTightness()
{
    return 1.0;
}

void WRegular::read(istream& is)
{

}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

