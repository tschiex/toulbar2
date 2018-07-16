/*----------------------------------------------------------------------
 *
 * Incremental regular constraint
 *
 */

#include "tb2regular.hpp"
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

bool DACCompare(const EnumeratedVariable* v1, const EnumeratedVariable* v2)
{
    return (v1->getDACOrder() < v2->getDACOrder());
}


// ---------------------------------------------------
// ------------ The WRegular class -------------------
// ---------------------------------------------------
WRegular::WRegular(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, WFA &automaton)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , lb(MIN_COST)
{
    // Copy the scope and DAC order it
    DACScope.assign(scope_in,scope_in+arity_in);
    sort(DACScope.begin(),DACScope.end(),DACCompare);

    // Unroll the automata
    // Only DFA for now: one initial state, one transition per value
    if (automaton.getAcceptingStates().size() != 1) {
        cerr << "Only single accepting states automata can be used with regular.\n";
        exit(EXIT_FAILURE);
    }
    map<unsigned int, unsigned int> stateTranslate;
    stateTranslate[automaton.getAcceptingStates()[0].first] = 0;
    vector<unsigned int> stateRevTranslate;
    stateRevTranslate.push_back(automaton.getAcceptingStates()[0].first);

    Cost initialWeight = automaton.getAcceptingStates()[0].second;
    vector<map<int,pair<int,Cost>>> transitions; //per state and value: target and weight
    transitions.resize(automaton.getNbStates());
    for (const auto &transition : automaton.getTransitions()) {
        if (transitions[transition->start].find(transition->symbol) ==  transitions[transition->start].end())
            transitions[transition->start][transition->symbol]= make_pair(transition->end,transition->weight);
        else {
            cerr << "Only deterministic automata can be used with regular.";
            exit(EXIT_FAILURE);
        }
    }

    // map of final states to prune / weight last layer 
    map <int,Cost> finalStates;
    for (auto accepting : automaton.getAcceptingStates())
        finalStates[accepting.first] = accepting.second;
    
    //Start layering
    map<unsigned int, unsigned int> nextTranslate; // maps automaton states to nodeIdx in layer
    vector<unsigned int> nextRevTranslate; // the other way around
    // Refs to swap at each layer
    map<unsigned int, unsigned int> &currentTranslateRef = stateTranslate;
    map<unsigned int, unsigned int> &nextTranslateRef = nextTranslate;
    vector<unsigned int> &currentRevTranslateRef = stateRevTranslate;
    vector<unsigned int> &nextRevTranslateRef = nextRevTranslate;
    // variables starting with "a/u" are resp. automaton/unrolled automaton states.
    // second letter may be c or n for current/next layer
    // Translate got from a to u, RevTranslate from u to a

    // allocate for all layers
    allArcs.resize(DACScope.size());
    outDegrees.resize(DACScope.size());

    for (unsigned int layer = 0; layer < DACScope.size(); layer ++){
        unsigned int unStateId = 0;
        // allocate for all layer states
        allArcs[layer].resize(stateRevTranslate.size());
        outDegrees[layer].resize(stateRevTranslate.size());
        for (unsigned int ucState = 0; ucState < stateRevTranslate.size();  ucState++) {
            unsigned int acState = currentRevTranslateRef[ucState];
            // allocate for all symbols/values: use domain size
            allArcs[layer][ucState].resize(DACScope[layer]->getDomainInitSize());

            for (const auto& transition : transitions[acState]) {
                unsigned int anTarget = transition.second.first;
                if (layer != DACScope.size()-1 || finalStates.count(anTarget)) { // last layer: prune with accepting states
                    if (!nextTranslate.count(anTarget)) // the target is unknown and must be mapped
                        nextTranslateRef[anTarget] = unStateId++; 
                    // let's create arcs
                    unsigned int unTarget =  nextTranslateRef[anTarget]; // get the target's id in the new layer
                    Cost arcWeight = transition.second.second;
                    if (layer == 0) arcWeight += initialWeight;
                    if (layer == DACScope.size()-1 && finalStates.count(anTarget))
                        arcWeight += finalStates[anTarget];
                    Arc* myArc = new Arc(layer,ucState,transition.first,transition.second.second,unTarget);
                    allArcs[layer][ucState][transition.first] = myArc;
                }
            }
        }
        // all states done, lets update translators: clear current and swap
        currentTranslateRef.clear();
        currentRevTranslateRef.clear();
        swap(currentTranslateRef, nextTranslateRef);
        swap(currentRevTranslateRef, nextRevTranslateRef);
    }

    // Now we should clean states with no successors (last layer done)
    // TODO

    // Init conflict weights
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
        cout << "After assign of " << idx << "\n";
}


void WRegular::remove(int idx)
{
    static const bool debug{false};

    if (debug)
        cout << "In remove " << idx << endl;

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
        cout << "In projectFromZero " << "\n";

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

