#ifndef Node_hpp
#define Node_hpp

#include <string>
#include <vector>

class Node{
    public:
                                    Node(void);
        void                        addNeighbor(Node* p);
        Node*                       getAncestor(void) { return ancestor; }
        std::vector<Node*>          getDescendants(void);
        int                         getIndex(void) { return index; }
        bool                        getIsTip(void) { return isTip; }
        double                      getBranchLength(void){return branchLength; }
        std::string                 getName(void) { return name; }
        const std::vector<Node*>&   getNeighbors(void) { return neighbors; }
        void                        removeNeighbor(Node* p);
        void                        removeAllNeighbors(void) { neighbors.clear(); }
        void                        setAncestor(Node* p) { ancestor = p; }
        void                        setIndex(int x) { index = x; }
        void                        setIsTip(bool tf) { isTip = tf; }
        void                        setBranchLength(double x) {branchLength = x;}
        void                        setName(std::string s) { name = s; }
    
    private:
        //Functions
        
        //Objects organized from greatest amount of memory required to least
        std::vector<Node*>          neighbors;
        Node*                       ancestor;
        std::string                 name;
        double                      branchLength; //subtending branch
        int                         index;
        bool                        isTip;
};

#endif
