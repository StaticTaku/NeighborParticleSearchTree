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

    enum class SearchMode : unsigned char {
        GATHER,
        SYMMETRY
    };

    struct Node {
        Type type;
        double position[3];
        Node* next;
    };

    struct Body:public Node {
        unsigned int id;
        double search_radius = 0;
    };

    struct Cell:public Node {
        double max_search_radius = 0;
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
        NeighborParticleSearchTree(unsigned int _DIM, unsigned int reserve_number):m_reserve_num(reserve_number), DIM(_DIM), m_NSUB(1<<_DIM){
            TREE_PRINT_INFO("start\n");
            m_bodies = new Body[reserve_number];
            m_boundary_length = new double[_DIM];

            for(int i = 0;i<reserve_number;++i) {
                m_bodies[i].type   = Type::Body;
                m_bodies[i].id     = i;
            }
            TREE_PRINT_INFO("finish\n");
        }

        NeighborParticleSearchTree(const NeighborParticleSearchTree&) = delete;
        NeighborParticleSearchTree& operator=(const NeighborParticleSearchTree&) = delete;
        
        NeighborParticleSearchTree(NeighborParticleSearchTree&& o):
        m_reserve_num(o.m_reserve_num), m_size(o.m_size), DIM(o.DIM), m_NSUB(o.m_NSUB),
        m_rsize(o.m_rsize), m_first_call((o.m_first_call)), m_root(o.m_root), 
        m_bodies(o.m_bodies), m_free_cell(o.m_free_cell), m_boundary_length(o.m_boundary_length) {
            o.m_root            = nullptr;
            o.m_bodies          = nullptr;
            o.m_free_cell       = nullptr;
            o.m_boundary_length = nullptr;
        }

        NeighborParticleSearchTree& operator=(NeighborParticleSearchTree&& o) {
            if (&o == this)
                return *this;
            delete [] m_root, delete [] m_bodies, delete [] m_free_cell, delete [] m_boundary_length;

            m_root            = o.m_root;
            m_bodies          = o.m_bodies;
            m_free_cell       = o.m_free_cell;
            m_boundary_length = o.m_boundary_length;
            
            m_reserve_num = o.m_reserve_num, m_size = o.m_size, DIM = o.DIM, m_NSUB= o.m_NSUB,
            m_rsize= o.m_rsize, m_first_call = o.m_first_call, m_root= o.m_root, 
            m_bodies= o.m_bodies, m_free_cell= o.m_free_cell;

            o.m_root             = nullptr;
            o.m_bodies           = nullptr;
            o.m_free_cell        = nullptr;
            o.m_boundary_length  = nullptr;

            return *this;
        }

        ~NeighborParticleSearchTree() {
            TREE_PRINT_INFO("start\n");
            NewTree();
            Cell* c;
            while(m_free_cell != nullptr) {
                c = static_cast<Cell*>(m_free_cell);
                m_free_cell = c->next;
                delete c;
            }

            delete m_free_cell;
            delete[] m_bodies;
            delete[] m_boundary_length;
            TREE_PRINT_INFO("finish\n");
        }

        void Resize(int size) {
            if(size <= m_reserve_num)
                m_size = size;
            else
                TREE_PRINT_ERROR(stdout, "Size exceeds reserved size\n");
        }

        void CopyPos(double pos_x, unsigned int id, unsigned int dim) {
            if(id < m_size && dim < DIM)
                m_bodies[id].position[dim] = pos_x;
            else
                TREE_PRINT_ERROR(stdout, "Acces to outside of memory\n");
        }

        void CopySearchRadius(double search_radius, unsigned int id) {
            if(id < m_size)
                m_bodies[id].search_radius = search_radius;
            else
                TREE_PRINT_ERROR(stdout, "Acces to outside of memory\n");
        }

        double GetPos(unsigned int id, unsigned int dim) {
            if(id < m_size && dim < DIM)
                return m_bodies[id].position[dim];
            else
                TREE_PRINT_ERROR(stdout, "Acces to outside of memory\n");
        }
        
        void UpdateTree() {
            TREE_PRINT_INFO("start\n");
            Body* p;
            NewTree();
            m_root = MakeCell();
            m_root->ClearPosition(DIM);
            ExpandBox();
            for(p = m_bodies;p<m_bodies+m_size;++p)
                LoadBody(p);
            //PropagateInfo(m_root, m_rsize, 0);
            ThreadTree(m_root,nullptr);
            TREE_PRINT_INFO("finish\n");
        }
        
        template <SearchMode SEARCH_MODE = SearchMode::GATHER>
        void FindNeighborParticle(const double* pos, const double radius, std::vector<unsigned int>& interaction_list, bool clear = true) {
            TREE_PRINT_INFO("start\n");
            if(clear)
                interaction_list.clear();
            WalkTree<SEARCH_MODE>(pos,radius,interaction_list,m_root,m_rsize);
            TREE_PRINT_INFO("finish\n");
        }

        template <SearchMode SEARCH_MODE = SearchMode::GATHER>
        void FindNeighborParticleWithPeriodicBoundary(const double* pos, const double radius, const double* boundary_length, std::vector<unsigned int>& interaction_list, bool clear = true) {
            TREE_PRINT_INFO("start\n");
            if(clear)
                interaction_list.clear();
            for(int dim = 0;dim < DIM; ++dim)
                m_boundary_length[dim] = boundary_length[dim];
            WalkTreeWithPeriodicBoundary<SEARCH_MODE>(pos,radius,interaction_list,m_root,m_rsize);
            TREE_PRINT_INFO("finish\n");
        }

    private:
        int m_reserve_num;
        int m_size;
        Cell* m_root;//root pointer
        Body* m_bodies;
        Node* m_free_cell = nullptr;
        double* m_boundary_length;
        int DIM;
        int m_NSUB;
        double m_rsize = 1;//m_root size
        bool m_first_call = true;

        void NewTree() {
            Node* p;
            if(!m_first_call) {
                p = m_root;
                while(p != nullptr) {
                    if(p->type == Type::Cell) {
                        p->next = m_free_cell;
                        m_free_cell = p;
                        p = static_cast<Cell*>(p)->more;
                    }else {
                        p = p->next;
                    }
                }
            }else {
                m_first_call = false;
            }
            m_root = nullptr;
        }

        Cell* MakeCell() {
            Cell* c;
            int i;
            if(m_free_cell == nullptr) {
                c = new Cell(DIM);
            }
            else {
                c = static_cast<Cell*>(m_free_cell);
                m_free_cell = c->next;
            }
            c->type = Type::Cell;
            for(int i = 0;i<m_NSUB;++i)
                c->subP[i] = nullptr;
            c->max_search_radius = 0;
            return c;
        }

        void ExpandBox() {
            m_rsize = 1;
            double dmax = 0,d;
            Body* p;
            for(p = m_bodies;p<m_bodies + m_size;++p) {
                for(int dim = 0;dim<DIM;++dim) {
                    d = std::abs(p->position[dim]-m_root->position[dim]);
                    if(d > dmax)
                        dmax = d;
                }
            }

            while(m_rsize < 2*dmax)
                m_rsize = 2*m_rsize;
        }

        void LoadBody(Body* p) {
            Cell* q; Cell* c;
            int qind,k;
            double qsize, dist2;

            q = m_root;
            qind = SubIndex(p,q);
            qsize = m_rsize;
            q->max_search_radius = std::max(q->max_search_radius, p->search_radius);
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
                    c->max_search_radius = p->search_radius; //initialization of max_search_radius
                    q->subP[qind] = c;
                }
                q = static_cast<Cell*>(q->subP[qind]);
                q->max_search_radius = std::max(q->max_search_radius, p->search_radius);
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
                    ind += m_NSUB >> (k+1);
            }
            return ind;
        }

        void PropagateInfo(Cell* p, double psize, int lev) {
            Node* q;
            for(int i = 0;i<m_NSUB;++i) {
                if((q=p->subP[i]) != nullptr) {
                    if(q->type == Type::Cell)
                        PropagateInfo(static_cast<Cell*>(q), psize/2, lev+1);
                }
            }
        }

        void ThreadTree(Node* p, Node* n) {
            int ndesc, i;
            Node* desc[m_NSUB+1];

            p->next = n;
            if(p->type == Type::Cell)
            {
                ndesc = 0;
                for(i = 0;i<m_NSUB;++i)
                    if(static_cast<Cell*>(p)->subP[i] != nullptr)
                        desc[ndesc++] = static_cast<Cell*>(p)->subP[i];
                
                static_cast<Cell*>(p)->more = desc[0];
                desc[ndesc] = n;
                for(i = 0;i<ndesc;++i)
                    ThreadTree(desc[i], desc[i+1]);
            }
        }

        template <SearchMode SEARCH_MODE>
        void WalkTree(const double* pos,const double radius,std::vector<unsigned int>& interaction_list,Cell* p,double psize) {
            Node* q;
            //search all of p's direct descendants
            if constexpr (SEARCH_MODE == SearchMode::GATHER) {
                for(q = p->more;q != p->next;q = q->next) {
                    if(q->type == Type::Cell && isNearTarget(pos,radius,q->position,psize/2))
                        WalkTree<SEARCH_MODE>(pos,radius,interaction_list,static_cast<Cell*>(q),psize/2);
                    else if(q->type == Type::Body)
                    {
                        double length = 0;

                        for(int dim = 0;dim<DIM;++dim)
                            length += (pos[dim] - q->position[dim])*(pos[dim] - q->position[dim]);
                        
                        if(length <= radius*radius)
                            interaction_list.emplace_back(static_cast<Body*>(q)->id);
                    }
                }
            }else if constexpr (SEARCH_MODE == SearchMode::SYMMETRY) {
                for(q = p->more;q != p->next;q = q->next) {
                    if(q->type == Type::Cell && (isNearTargetWithPeriodicBoundary(pos,radius,q->position,psize/2) || isNearTargetWithPeriodicBoundary(pos, static_cast<Cell*>(q)->max_search_radius, q->position, psize/2)))
                        WalkTree<SEARCH_MODE>(pos,radius,interaction_list,static_cast<Cell*>(q),psize/2);
                    else if(q->type == Type::Body)
                    {
                        double length = 0;

                        for(int dim = 0;dim<DIM;++dim)
                            length += (pos[dim] - q->position[dim])*(pos[dim] - q->position[dim]);
                        
                        double search_radius = static_cast<Body*>(q)->search_radius;
                        if((length <= radius*radius) || (length <= search_radius * search_radius))
                            interaction_list.emplace_back(static_cast<Body*>(q)->id);
                    }
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

        template <SearchMode SEARCH_MODE>
        void WalkTreeWithPeriodicBoundary(const double* pos,const double radius,std::vector<unsigned int>& interaction_list,Cell* p,double psize) {
            Node* q;
            //search all of p's direct descendants
            if constexpr (SEARCH_MODE == SearchMode::GATHER) {
                for(q = p->more;q != p->next;q = q->next) {
                    if(q->type == Type::Cell && isNearTargetWithPeriodicBoundary(pos,radius,q->position,psize/2))
                        WalkTreeWithPeriodicBoundary<SearchMode::GATHER>(pos,radius,interaction_list,static_cast<Cell*>(q),psize/2);
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
            }else if constexpr (SEARCH_MODE == SearchMode::SYMMETRY) {
                for(q = p->more;q != p->next;q = q->next) {
                    if(q->type == Type::Cell && (isNearTargetWithPeriodicBoundary(pos,radius,q->position,psize/2) || isNearTargetWithPeriodicBoundary(pos, static_cast<Cell*>(q)->max_search_radius, q->position, psize/2)))
                        WalkTreeWithPeriodicBoundary<SearchMode::SYMMETRY>(pos,radius,interaction_list,static_cast<Cell*>(q),psize/2);
                    else if(q->type == Type::Body) {
                        double length = 0;

                        for(int dim = 0;dim<DIM;++dim) {
                            double dx = PeriodicDistance(pos[dim], q->position[dim], dim);
                            length += dx * dx;
                        }

                        double search_radius = static_cast<Body*>(q)->search_radius;
                        if((length <= radius*radius) || (length <= search_radius * search_radius))
                            interaction_list.emplace_back(static_cast<Body*>(q)->id);
                    }
                }
            }
        }

        double PeriodicDistance(double x1, double x2, int dim) {
            double X = x1-x2;
            if(X>0.5*m_boundary_length[dim]){
                X -= m_boundary_length[dim];
            } else if(X<-0.5*m_boundary_length[dim]){
                X += m_boundary_length[dim];
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