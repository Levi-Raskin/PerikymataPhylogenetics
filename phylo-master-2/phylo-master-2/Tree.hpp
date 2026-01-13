#ifndef Tree_hpp
#define Tree_hpp

#include <vector>

class Node;

class Tree{
    public:
                                    Tree(void) = delete;
                                    Tree(int nt, double lambda);
                                    ~Tree(void);
        void                        initializeDownPassSequence(void);
        const std::vector<Node*>&   getDownPassSequence(void) { return downPassSequence; }
        int                         getNumTaxa(void) { return numTaxa; }
        Node*                       getRoot(void) {return root;}
        void                        print(void);
    private:
        //Functions
        Node*                       addNode(void);
        void                        setNeighbors(Node* p, Node* q);
        void                        showNode(Node* p, int indent);
        void                        passDown(Node* p, Node* from);
        void                        removeAsNeighbors(Node* p, Node* q);
        //Objects, ordered by memory footprint
        std::vector<Node*>          downPassSequence;
        std::vector<Node*>          nodes;
        Node*                       root;
        int                         numTaxa;
};

#endif
