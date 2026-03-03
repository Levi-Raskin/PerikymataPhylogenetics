#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
class RandomVariable;


class Node {

    public:
                            Node(void);
        void                addNeighbor(Node* p);
        Node*               chooseNeighborAtRandom(RandomVariable* rng, Node* excludingNode);
        Node*               chooseNeighborAtRandom(RandomVariable* rng, std::set<Node*> excludedNodes);
        Node*               getAncestor(void) { return ancestor; }
        std::vector<Node*>& getDescendants(void);
        bool                getFlag(void) { return flag; }
        bool                getHasData(void){ return hasData;}
        int                 getIndex(void) { return index; }
        bool                getIsFossil(void) { return isFossil; }
        bool                getIsTip(void) { return isTip; }
        double              getLength(void){return branchLength; }
        std::string         getName(void) { return name; }
        std::set<Node*>&    getNeighbors(void) { return neighbors; }
        int                 getOffset(void) { return offset; }
        int                 getScratchInt(void) { return scratchInt; }
        double              getTime(void) { return time;}
        void                removeNeighbor(Node* p);
        void                removeAllNeighbors(void) { neighbors.clear(); }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setFlag(bool tf) { flag = tf; }
        void                setHasData(bool tf) { hasData = tf;}
        void                setIndex(int x) { index = x; }
        void                setIsFossil(bool tf) { isFossil = tf;}
        void                setIsTip(bool tf) { isTip = tf; }
        void                setLength(double x) {branchLength = x;}
        void                setName(std::string s) { name = s; }
        void                setOffset(int x) { offset = x; }
        void                setScratchInt(int x) { scratchInt = x; }
        void                setTime(double x) {time = x;}
    
    private:
        Node*               ancestor;
        double              branchLength;
        std::vector<Node*>  descendantsVector;
        bool                flag;
        bool                hasData;
        int                 index;
        bool                isFossil;
        bool                isTip;
        std::string         name;
        std::set<Node*>     neighbors;
        int                 offset;
        int                 scratchInt;
        double              time;
};

#endif
