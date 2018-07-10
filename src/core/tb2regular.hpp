#ifndef TB2REG_HPP_
#define TB2REG_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

/* Enforces the weighted regular constraint. Closely related to weighted MDD 
   except for the absence of long edges
*/
class WRegular : public AbstractNaryConstraint
{
public:
    WRegular(WCSP *wcsp, EnumeratedVariable **scope_in, int arity_in,
                     vector<vector<int>> clq_in, int rhs_in);
    WRegular(WCSP *wcsp, EnumeratedVariable **scope_in, int arity_in);
    ~WRegular();

    void read(istream& file);

    bool extension() const FINAL { return false; }

    void assign(int idx) override;
    void remove(int idx) override;
    void increase(int idx) override;
    void decrease(int idx) override;
    void projectFromZero(int idx) override;

    void propagate() override;

    Cost eval( const String& s ) override {
         return MIN_COST;
    }

    vector<Long> conflictWeights;   // used by weighted degree heuristics
    Long getConflictWeight(int varIndex) const override {assert(varIndex>=0);assert(varIndex<arity_);return conflictWeights[varIndex]+Constraint::getConflictWeight();}
    void incConflictWeight(Constraint *from) override {
        //assert(fromElim1==NULL);
        //assert(fromElim2==NULL);
        if (from==this) {
            Constraint::incConflictWeight(1);
        } else if (deconnected()) {
            for (int i=0; i<from->arity(); i++) {
                int index = getIndex(from->getVar(i));
                if (index>=0) { // the last conflict constraint may be derived from two binary constraints (boosting search), each one derived from an n-ary constraint with a scope which does not include parameter constraint from
                    assert(index < arity_);
                    conflictWeights[index]++;
                }
            }
        }
    }
    double computeTightness() override;
    void dump(ostream&, bool) override {}
private:
    // ----------------------------------------------------------------------
    // definition

    //----------------------------------------------------------------------
    // operation

    // amount we have already projected to c_zero
    StoreCost lb;

public:
 
};

inline std::ostream& operator<<(std::ostream& os, WRegular::state s)
{
    return s.print(os);
}

#endif /* TB2CLQ_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

