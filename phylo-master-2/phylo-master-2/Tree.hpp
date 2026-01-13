#ifndef Tree_hpp
#define Tree_hpp

#include <vector>

class Node;

class Tree{
    public:
                                    Tree(void) = delete;
                                    Tree(int nt, double lambda);
                                    Tree(std::vector<std::string> taxonNames, double lambda);
                                    Tree(std::string newick);
                                    Tree(const Tree& t);
        Tree&                       operator=(const Tree& t);
                                    ~Tree(void);
        void                        initializeDownPassSequence(void);
        const std::vector<Node*>&   getDownPassSequence(void) { return downPassSequence; }
        int                         getNumTaxa(void) { return numTaxa; }
        Node*                       getRoot(void) {return root;}
        void                        print(void);
    private:
        //Functions
        Node*                       addNode(void);
        void                        clone(const Tree& t);
        void                        deleteNodes(void);
        void                        setNeighbors(Node* p, Node* q);
        void                        setOffsets(void);
        void                        showNode(Node* p, int indent);
        void                        passDown(Node* p, Node* from);
        std::vector<std::string>    parseNewickString(std::string);
        void                        removeAsNeighbors(Node* p, Node* q);
        //Objects, ordered by memory footprint
        std::vector<Node*>          downPassSequence;
        std::vector<Node*>          nodes;
        Node*                       root;
        int                         numTaxa;
};

#endif
