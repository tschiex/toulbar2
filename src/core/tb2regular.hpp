#ifndef TB2REG_HPP_
#define TB2REG_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "utils/tb2automaton.hpp"

/* 
  Includes a structure to describe an unrolled automaton anfd enforce local consistencies incrementally
  This is a layered weighted acyclic DAG with one layer per variable in the scope.
  Targeted structure: 
     - layers numbered from 0, DAC order
     - each layer starts with a sequence of states (indexed from 0): not stored (only max)
     - arcs in a layer connect states from the layer to the next layer. 
        - the source and target are state numbers in their respective layer
        - the arc has a value and a weight (Cost)
        - the arc can become inactive if it belongs to no feasible path, will be removed from the list of active arcs
    - All arcs are accessible per layer/nodeidx/value (not backtrackable)
    - All arcs are accessible by layer and nodeidx (not backtrackable)
    - Backtrackabkle degrees can be stored in side vectors (in/out) if needed
    - We have a vector of backtrackable lists of arcs per layer/value that plays the role of the Qij (reverse?)
*/

class WRegular : public AbstractNaryConstraint
{
public:
    WRegular(WCSP *wcsp, EnumeratedVariable **scope_in, int arity_in, WFA &automata);
    //WRegular(WCSP *wcsp, EnumeratedVariable **scope_in, int arity_in), sequence, distance, matrix);
    ~WRegular();

    struct Arc {
        Arc(int layer, int source, int value, Cost weight, int target);
        int layer; // index to the variable associated in arc in the DAC order
        int source;
        int target;  // backtrackable?
        Value value; // value associated with arc
        Cost weight;
        bool active;// in Qlayer/value

        int get_source();
        int get_target();
        int get_weight();
    };

    int get_width();
    int get_layer_width(int layer);
    int get_layer_num();
    int get_arc_num();

    Arc* get_arc(int layer, int source, Value value) {
        // asserts!
        return allArcs[layer][source][value];
    }

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
    std::ostream& printstate(std::ostream& os);

    private:
    vector<EnumeratedVariable*> DACScope;
    vector<vector<vector<Arc*>>> allArcs; // For each layer/source-node/value

    vector<vector<BTList<Arc>>> arcsAtLayerValue; /* arcs in feasible paths in a given layer/value: better than Qij ? */
    StoreCost lb;         // amount we have already projected to c_zero
    vector<vector<StoreCost>> delta; // one per variable and value
    vector<vector<StoreCost>> alpha; // shortest path length from source to layer/node for AC
    vector<vector<StoreCost>> alphap; // shortest path length from source to layer/node with unaries for DAC
    vector<vector<StoreCost>> beta; // shortest path length from sink to layer/node for AC
    vector<vector<StoreInt>> Pred; // nodeidx of the previous state that gives alpha (or is the arc more useful?)
    vector<vector<StoreInt>> Predp; // nodeidx of the previous state that gives alphap (or is the arc more useful?)
    vector<vector<StoreInt>> Succ; // nodeidx of the next state that gives beta (or is the arc more useful?)
};

#endif /* TB2REG_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */