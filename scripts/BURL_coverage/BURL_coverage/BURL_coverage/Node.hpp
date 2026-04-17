#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
class RandomVariable;


class Node {

    public:
                            Node(void);
        void                addNeighbor(Node* p);
        Node*               getAncestor(void) { return ancestor; }
        std::vector<Node*>& getDescendants(void);
        bool                getFlag(void) { return flag; }
        int                 getIndex(void) { return index; }
        bool                getIsTip(void) { return isTip; }
        std::string         getName(void) { return name; }
        std::set<Node*>&    getNeighbors(void) { return neighbors; }
        int                 getOffset(void) { return offset; }
        void                removeNeighbor(Node* p);
        void                removeAllNeighbors(void) { neighbors.clear(); }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setFlag(bool tf) { flag = tf; }
        void                setIndex(int x) { index = x; }
        void                setIsTip(bool tf) { isTip = tf; }
        void                setName(std::string s) { name = s; }
        void                setOffset(int x) { offset = x; }
    
    private:
        std::vector<Node*>  descendantsVector;
        std::set<Node*>     neighbors;
        Node*               ancestor;
        std::string         name;
        bool                isTip;
        bool                flag;
        int                 index;
        int                 offset;
};

#endif
