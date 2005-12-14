/** \file tb2queue.hpp
 *  \brief Propagation queue with time stamping.
 * 
 */

#ifndef TB2QUEUE_HPP_
#define TB2QUEUE_HPP_

#include "tb2btlist.hpp"

typedef enum {NOTHING_EVENT=0, INCREASE_EVENT=1, DECREASE_EVENT=2} EventType;

struct CostVariableWithTimeStamp 
{
    CostVariable *var;
    long long timeStamp;
    int incdec;
};

class Queue : private BTList<CostVariableWithTimeStamp>
{  
    // make it private because we don't want copy nor assignment
    Queue(const Queue &s);
    Queue& operator=(const Queue &s);
    
public:
    Queue() : BTList<CostVariableWithTimeStamp>(NULL) {}
    
    int getSize() const {return BTList<CostVariableWithTimeStamp>::getSize();}
    bool empty() const {return BTList<CostVariableWithTimeStamp>::empty();}

    void clear() {BTList<CostVariableWithTimeStamp>::clear();}
    
    void push(DLink<CostVariableWithTimeStamp> *elt, long long curTimeStamp) {
        if (elt->content.timeStamp < curTimeStamp) {
            elt->content.timeStamp = curTimeStamp;
            push_back(elt, false);
        }
    }
    
    void push(DLink<CostVariableWithTimeStamp> *elt, EventType incdec, long long curTimeStamp) {
        elt->content.incdec |= incdec;
        push(elt, curTimeStamp);
    }
    
    CostVariable *pop() {
        assert(!empty());
        DLink<CostVariableWithTimeStamp> *elt = pop_back(false);
        elt->content.timeStamp = -1;
        elt->content.incdec = NOTHING_EVENT;
        return elt->content.var;
    }
    
    CostVariable *pop(int *incdec) {
        assert(!empty());
        *incdec = (*rbegin()).incdec;
        return pop();
    }
    
    CostVariable *pop_min();
    CostVariable *pop_min(int *incdec);
    CostVariable *pop_max();
    CostVariable *pop_max(int *incdec);
};

#endif /*TB2QUEUE_HPP_*/
