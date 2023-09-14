#pragma once

#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <cmath>
#include <limits>

//#define TREE_DEBUG

#define TREE_PRINTF(dst, ...)                                       \
    {                                                               \
        fprintf(dst, "%s(%d): %s\n", __FILE__, __LINE__, __func__); \
        fprintf(dst, __VA_ARGS__);                                  \
        fflush(NULL);                                               \
    }

#ifdef TREE_DEBUG
#define TREE_PRINT_INFO(...)              \
    {                                     \
        fprintf(stdout, "[INFO] ");       \
        TREE_PRINTF(stdout, __VA_ARGS__); \
    }
#else
#define TREE_PRINT_INFO(...)
#endif

#define TREE_PRINT_ERROR(dst, ...)     \
    {                                  \
        fprintf(stdout, "[ERROR] ");   \
        TREE_PRINTF(dst, __VA_ARGS__); \
        exit(EXIT_FAILURE);            \
    }

namespace Tree {
    enum class Type : unsigned char {
        Body = 0,
        Cell,
    };

    struct Node {
        Type type;
        bool update;
        double position[3];
        Node* next;
    };

    struct Body:public Node {
        unsigned int id;
    };

    struct Cell:public Node {
        bool flag = false;
        Node* more = nullptr;
        Node** subP;

        Cell() {}
        
        Cell(unsigned char DIM):flag(true) {
            subP = new Node*[1<<DIM];
        }

        void Construct(unsigned int DIM) {
            subP = new Node*[1<<DIM];
            flag = true;
        }

        ~Cell() {
            if(flag) {   
                delete[] subP;
            }
        }

        void ClearPosition(unsigned char DIM) {
            for(int i = 0;i<DIM;++i)
                position[i] = 0;
        }
    };

    class NeighborParticleSearchTree {
    public:
        NeighborParticleSearchTree(unsigned int _DIM, unsigned int reserve_number):m_reserve_num(reserve_number), DIM(_DIM), NSUB(1<<_DIM){
            TREE_PRINT_INFO("start\n");
            bodies = new Body[reserve_number];
            boundary_length_ = new double[_DIM];

            for(int i = 0;i<reserve_number;++i) {
                bodies[i].type   = Type::Body;
                bodies[i].update = true;
                bodies[i].id     = i;
            }
            TREE_PRINT_INFO("finish\n");
        }

        ~NeighborParticleSearchTree() {
            TREE_PRINT_INFO("start\n");
            NewTree();
            Cell* c;
            while(freeCell != nullptr) {
                c = static_cast<Cell*>(freeCell);
                freeCell = c->next;
                delete c;
            }
            if(freeCell != nullptr)
                delete freeCell;
        
            delete[] bodies;
            delete[] boundary_length_;
            TREE_PRINT_INFO("finish\n");
        }

        void Resize(int size) {
            if(size <= m_reserve_num )
                m_size = size;
            else {
                TREE_PRINT_ERROR(stdout, "Size exceeds reserved size\n");
            }
        }

        void CopyPos(double pos_x, unsigned int id, unsigned int dim) {
            if(id < m_size && dim < DIM)
                bodies[id].position[dim] = pos_x;
            else
                TREE_PRINT_ERROR(stdout, "Acces to outside of memory\n");
        }

        double GetPos(unsigned int id, unsigned int dim) {
            if(id < m_size && dim < DIM)
                return bodies[id].position[dim];
            else
                TREE_PRINT_ERROR(stdout, "Acces to outside of memory\n");
        }
        
        void UpdateTree() {
            TREE_PRINT_INFO("start\n");
            Body* p;
            NewTree();
            root = MakeCell();
            root->ClearPosition(DIM);
            ExpandBox();
            for(p = bodies;p<bodies+m_size;++p)
                LoadBody(p);
            tdepth = 0;
            PropagateInfo(root, rsize, 0);
            ThreadTree(root,nullptr);
            TREE_PRINT_INFO("finish\n");
        }

        void FindNeighborParticle(const double* pos, const double radius, std::vector<unsigned int>& interaction_list, bool clear = true) {
            TREE_PRINT_INFO("start\n");
            if(clear)
                interaction_list.clear();
            WalkTree(pos,radius,interaction_list,root,rsize);
            TREE_PRINT_INFO("finish\n");
        }

        void FindNeighborParticleWithPeriodicBoundary(const double* pos, const double radius, const double* boundary_length, std::vector<unsigned int>& interaction_list, bool clear = true) {
            TREE_PRINT_INFO("start\n");
            if(clear)
                interaction_list.clear();
            for(int dim = 0;dim < DIM; ++dim)
                boundary_length_[dim] = boundary_length[dim];
            WalkTreeWithPeriodicBoundary(pos,radius,interaction_list,root,rsize);
            TREE_PRINT_INFO("finish\n");
        }

    private:
        int m_reserve_num;
        int m_size;

        Cell* root;//root pointer
        Body* bodies;
        Node* freeCell = nullptr;
        Node** active;

        int DIM;
        int NSUB;
        int tdepth;//木の高さ
        double length;
        double rsize = 1;//root size
        double* boundary_length_;

        bool first_call = true;

        void NewTree() {
            Node* p;
            if(!first_call) {
                p = root;
                while(p != nullptr) {
                    if(p->type == Type::Cell) {
                        p->next = freeCell;
                        freeCell = p;
                        p = static_cast<Cell*>(p)->more;
                    }else {
                        p = p->next;
                    }
                }
            }else {
                first_call = false;
            }
            root = nullptr;
        }

        Cell* MakeCell() {
            Cell* c;
            int i;
            if(freeCell == nullptr) {
                c = new Cell(DIM);
            }
            else {
                c = static_cast<Cell*>(freeCell);
                freeCell = c->next;
            }
            c->type = Type::Cell;
            c->update = false;
            for(int i = 0;i<NSUB;++i)
                c->subP[i] = nullptr;
            
            return c;
        }

        void ExpandBox() {
            rsize = 1;
            double dmax = 0,d;
            Body* p;
            for(p = bodies;p<bodies + m_size;++p) {
                for(int dim = 0;dim<DIM;++dim) {
                    d = std::abs(p->position[dim]-root->position[dim]);
                    if(d > dmax)
                        dmax = d;
                }
            }

            while(rsize < 2*dmax)
                rsize = 2*rsize;
        }

        void LoadBody(Body* p) {
            Cell* q; Cell* c;
            int qind,k;
            double qsize, dist2;

            q = root;
            qind = SubIndex(p,q);
            qsize = rsize;
            
            while(q->subP[qind] != nullptr) {
                if(q->subP[qind]->type == Type::Body) {
                    double dist2 = 0;
                    for(int d = 0;d<DIM;++d)
                        dist2 += (p->position[d]-q->subP[qind]->position[d])*(p->position[d]-q->subP[qind]->position[d]);

                    if (dist2 == 0.0)
                        TREE_PRINT_ERROR(stdout, "Particles in the same position exist.\n");

                    c = MakeCell();
                    for(int k = 0;k<DIM;++k)
                        c->position[k] = q->position[k] + (p->position[k] < q->position[k] ? -qsize:qsize)/4;
                    c->subP[SubIndex(static_cast<Body*>(q->subP[qind]),c)] = q->subP[qind];
                    q->subP[qind] = c;
                }
                q = static_cast<Cell*>(q->subP[qind]);
                qind = SubIndex(p,q);
                qsize = qsize/2;
                if(qsize == 0)
                    TREE_PRINT_ERROR(stdout, "Tree is so deep that a cell size reaches zero\n");
            }
            q->subP[qind] = p;
        }

        int SubIndex(Body* p, Cell* q) {
            int ind, k;
            ind = 0;
            for(k = 0;k<DIM;++k) {
                if(q->position[k] <= p->position[k])
                    ind += NSUB >> (k+1);
            }
            return ind;
        }

        void PropagateInfo(Cell* p, double psize, int lev) {
            Node* q;
            tdepth = std::max(tdepth,lev);
            for(int i = 0;i<NSUB;++i) {
                if((q=p->subP[i]) != nullptr) {
                    if(q->type == Type::Cell)
                        PropagateInfo(static_cast<Cell*>(q), psize/2, lev+1);
                    p->update |= q->update;
                }
            }
        }

        void ThreadTree(Node* p, Node* n) {
            int ndesc, i;
            Node* desc[NSUB+1];

            p->next = n;
            if(p->type == Type::Cell)
            {
                ndesc = 0;
                for(i = 0;i<NSUB;++i)
                    if(static_cast<Cell*>(p)->subP[i] != nullptr)
                        desc[ndesc++] = static_cast<Cell*>(p)->subP[i];
                
                static_cast<Cell*>(p)->more = desc[0];
                desc[ndesc] = n;
                for(i = 0;i<ndesc;++i)
                    ThreadTree(desc[i], desc[i+1]);
            }
        }

        void WalkTree(const double* pos,const double radius,std::vector<unsigned int>& interaction_list,Cell* p,double psize) {
            Node* q;
            //search all of p's direct descendants
            for(q = p->more;q != p->next;q = q->next) {
                if(q->type == Type::Cell && isNearTarget(pos,radius,q->position,psize/2))
                    WalkTree(pos,radius,interaction_list,static_cast<Cell*>(q),psize/2);
                else if(q->type == Type::Body)
                {
                    double length = 0;

                    for(int dim = 0;dim<DIM;++dim)
                        length += (pos[dim] - q->position[dim])*(pos[dim] - q->position[dim]);

                    if(length <= radius*radius)
                        interaction_list.emplace_back(static_cast<Body*>(q)->id);
                }
            }
        }

        bool isNearTarget(const double* pos, double radius,double* posCell, double cellSize) {
            double dx, farLen;

            farLen = cellSize + radius;

            for(int dim = 0;dim < DIM;++dim) {
                dx = pos[dim] - posCell[dim];
                if(std::abs(dx) > farLen)
                    return false;
            }

            return true;
        }

        void WalkTreeWithPeriodicBoundary(const double* pos,const double radius,std::vector<unsigned int>& interaction_list,Cell* p,double psize) {
            Node* q;
            //search all of p's direct descendants
            for(q = p->more;q != p->next;q = q->next) {
                if(q->type == Type::Cell && isNearTargetWithPeriodicBoundary(pos,radius,q->position,psize/2))
                    WalkTreeWithPeriodicBoundary(pos,radius,interaction_list,static_cast<Cell*>(q),psize/2);
                else if(q->type == Type::Body) {
                    double length = 0;

                    for(int dim = 0;dim<DIM;++dim) {
                        double dx = PeriodicDistance(pos[dim], q->position[dim], dim);
                        length += dx * dx;
                    }

                    if(length <= radius*radius)
                        interaction_list.emplace_back(static_cast<Body*>(q)->id);
                }
            }
        }

        double PeriodicDistance(double x1, double x2, int dim) {
            double X = x1-x2;
            if(X>0.5*boundary_length_[dim]){
                X -= boundary_length_[dim];
            } else if(X<-0.5*boundary_length_[dim]){
                X += boundary_length_[dim];
            }
            return X;
        }

        bool isNearTargetWithPeriodicBoundary(const double* pos, double radius,double* posCell, double cellSize) {
            double dx, farLen;

            farLen = cellSize + radius;

            for(int dim = 0;dim < DIM;++dim) {
                dx = PeriodicDistance(pos[dim],posCell[dim],dim);
                if(std::abs(dx) > farLen)
                    return false;
            }

            return true;
        }
    };
}