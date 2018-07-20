/** \file tb2btlist.hpp
 *  \brief Backtrackable double-linked list.
 * 
 * Convention: 
 * 
 * elements can be inserted at the end of the list only
 * these insertions can be undone in the reverse order of their insertion
 * 
 * elements can be removed in any order
 * these removals can be undone in the reverse order of their removal.
 * 
 */

#ifndef TB2BTLIST_HPP_
#define TB2BTLIST_HPP_

#include "tb2store.hpp"

// One cell of a doubly linked list of <T>
template <class T>
class DLink {
public:
    bool removed; // true if the corresponding element has been removed
    DLink* next;
    DLink* prev;
    T content;

public:
    DLink<T>()
        : removed(true)
        , next(NULL)
        , prev(NULL)
    {
    }
};

// A backtrackable doubly linked list with a store stack
// that remembers which DLink cell must be restored where
template <class T>
class BTList {
    StoreStack<BTList, DLink<T>*>* storeUndo;
    int size;
    DLink<T>* head;
    DLink<T>* last;

public:
    BTList(StoreStack<BTList, DLink<T>*>* s)
        : storeUndo(s)
        , size(0)
        , head(NULL)
        , last(NULL)
    {
    }

    int getSize() const { return size; }
    bool empty() const { return size == 0; }

    // Warning! clear() is not a backtrackable operation
    void clear()
    {
        size = 0;
        head = NULL;
        last = NULL;
    }

    bool inBTList(DLink<T>* elt)
    {
        for (iterator iter = begin(); iter != end(); ++iter) {
            if (elt == iter.getElt())
                return !elt->removed;
        }
        return false;
    }

    void push_back(DLink<T>* elt, bool backtrack)
    {
        assert(!inBTList(elt));
        size++;
        elt->removed = false;
        if (last != NULL) {
            last->next = elt;
            elt->prev = last;
        } else {
            head = elt;
            elt->prev = NULL;
        }
        last = elt;
        last->next = NULL;
        if (backtrack)
            storeUndo->store(this, NULL);
    }

    void undoPushBack()
    {
        assert(last != NULL);
        size--;
        last->removed = true;
        if (last->prev != NULL) {
            last = last->prev;
            last->next->prev = NULL;
            last->next = NULL;
        } else {
            head = NULL;
            last = NULL;
        }
    }

    void erase(DLink<T>* elt, bool backtrack)
    {
        assert(!elt->removed);
        size--;
        elt->removed = true;
        if (elt->prev != NULL) {
            assert(!elt->prev->removed);
            assert(elt->prev->next == elt);
            elt->prev->next = elt->next;
        } else
            head = elt->next;
        if (elt->next != NULL) {
            assert(!elt->next->removed);
            assert(elt->next->prev == elt);
            elt->next->prev = elt->prev;
        } else
            last = elt->prev;
        if (backtrack) {
            storeUndo->store(this, elt->prev);
            storeUndo->store(this, elt);
        }
    }

    void undoErase(DLink<T>* elt, DLink<T>* prev)
    {
        assert(elt->removed);
        size++;
        elt->removed = false;
        if (prev != NULL) {
            assert(!prev->removed);
            elt->prev = prev;
            elt->next = prev->next;
            if (prev->next != NULL)
                prev->next->prev = elt;
            else
                last = elt;
            prev->next = elt;
        } else {
            if (head != NULL)
                head->prev = elt;
            else
                last = elt;
            elt->prev = NULL;
            elt->next = head;
            head = elt;
        }
    }

    // deprecated method to be used with erase(..) storing just one element
    //    void undoErase(DLink<T> *elt) {
    //        assert(elt->removed);
    //        size++;
    //        elt->removed = false;
    //        if (elt->prev != NULL) {
    //            assert(!elt->prev->removed);
    //            assert(elt->prev->next == elt->next);
    //            elt->prev->next = elt;
    //        } else head = elt;
    //        if (elt->next != NULL) {
    //            assert(!elt->next->removed);
    //            assert(elt->next->prev == elt->prev);
    //            elt->next->prev = elt;
    //        } else last = elt;
    //    }

    DLink<T>* pop_back(bool backtrack)
    {
        assert(last != NULL);
        DLink<T>* oldlast = last;
        erase(last, backtrack);
        return oldlast;
    }

    class iterator {
        DLink<T>* elt;

    public:
        iterator() { elt = NULL; }
        iterator(DLink<T>* e)
            : elt(e)
        {
        }

        T operator*() const
        {
            assert(elt != NULL);
            return elt->content;
        }

        DLink<T>* getElt() const { return elt; }

        iterator& operator++()
        { // Prefix form
            if (elt != NULL) {
                while (elt->next != NULL && elt->next->removed) {
                    elt = elt->next;
                }
                elt = elt->next;
            }
            assert(elt == NULL || !elt->removed);
            return *this;
        }

        iterator& operator--()
        { // Prefix form
            if (elt != NULL) {
                while (elt->prev != NULL && elt->prev->removed) {
                    elt = elt->prev;
                }
                elt = elt->prev;
            }
            assert(elt == NULL || !elt->removed);
            return *this;
        }

        // To see if you're at the end:
        bool operator==(const iterator& iter) const { return elt == iter.elt; }
        bool operator!=(const iterator& iter) const { return elt != iter.elt; }
    };

    iterator begin() { return iterator(head); }
    iterator end() { return iterator(NULL); }
    iterator rbegin() { return iterator(last); }
    iterator rend() { return end(); }
};

typedef BTList<ConstraintLink> ConstraintList;
typedef BTList<Variable*> VariableList;
typedef BTList<Separator*> SeparatorList;

/*
 * For internal use only! Interaction between tb2store and tb2btlist
 * 
 */

template <class T, class V>
template <class Q>
void StoreStack<T, V>::restore(BTList<Q>** l, DLink<Q>** elt, ptrdiff_t& x)
{
    if (elt[x] == NULL) {
        l[x]->undoPushBack();
    } else {
        assert(l[x] == l[x - 1]);
        l[x]->undoErase(elt[x], elt[x - 1]);
        x--;
    }
}

// BTList Wrapper for easier usage
template <typename T>
class DLinkStore {
private:
    int blkSize;
    StoreInt curEmpty;
    StoreInt curUsingBlkIndex;
    vector<DLink<T>*> blockStore;

public:
    typedef typename vector<DLink<T>*>::iterator iterator;

    DLinkStore(int blkSize_)
        : blkSize(blkSize_)
        , curEmpty(0)
        , curUsingBlkIndex(0)
    {
        blockStore.push_back(new DLink<T>[blkSize]);
    }
    ~DLinkStore()
    {
        for (vector<DLink<int>*>::iterator it = blockStore.begin();
             it != blockStore.end(); it++) {
            delete[] * it;
        }
        blockStore.clear();
    }

    DLink<T>* allocate(const T& ele)
    {
        if (curEmpty >= blkSize) {
            curEmpty = 0;
            curUsingBlkIndex = curUsingBlkIndex + 1;
            if (curUsingBlkIndex >= (int)blockStore.size()) {
                blockStore.push_back(new DLink<T>[blkSize]);
            }
        }
        DLink<int>* container = &(blockStore[curUsingBlkIndex][curEmpty]);
        container->content = ele;
        curEmpty = curEmpty + 1;
        return container;
    }
    iterator begin() { return blockStore.begin(); }
    iterator end() { return blockStore.end(); }
};

template <typename T>
class BTListWrapper {
private:
    BTList<T> list;
    DLinkStore<T>* dlinkStore;

public:
    typedef typename BTList<T>::iterator iterator;
    typedef typename vector<DLink<T>*>::iterator allIterator;

    BTListWrapper(DLinkStore<int>* dlinkStore)
        : list(&Store::storeIndexList)
        , dlinkStore(dlinkStore)
    {
    }

    ~BTListWrapper() {}

    size_t size()
    {
        return list.getSize();
    }
    bool empty()
    {
        return list.getSize() == 0;
    }
    void push_back(const T& ele)
    {
        DLink<T>* container = dlinkStore->allocate(ele);
        list.push_back(container, true);
    }
    void remove(const T& ele)
    {
        DLink<T>* target = NULL;
        for (typename BTList<T>::iterator it = list.begin(); it != list.end()
             && (target == NULL);
             ++it) {
            if (*it == ele) {
                target = it.getElt();
            }
        }
        if (target != NULL) {
            list.erase(target, true);
        }
    }
    void erase(iterator& it)
    {
        list.erase(it.getElt(), true);
    }
    iterator begin() { return list.begin(); }
    iterator end() { return list.end(); }
    allIterator allBegin() { return dlinkStore->begin(); }
    allIterator allEnd() { return dlinkStore->end(); }
};

#endif /*TB2BTLIST_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
