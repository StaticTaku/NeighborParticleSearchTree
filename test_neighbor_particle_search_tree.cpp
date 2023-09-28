#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>

#include "neighbor_particle_search_tree.hpp"

int main() {
{
    constexpr int DIM              = 3;    //spatial dimension
    constexpr int reserved_size    = 1000; //Number of particles inside tree allocated in memory in advance
    constexpr int num_of_particles = 100;  //Actual number of particles to be searched for. Note that num_of_particles <= reserved_size has to be true.
    double search_radius = 2;
    double point_of_search[DIM] = {5, 5, 5};
    std::vector<unsigned int> interaction_list;

    //////////////////////////////////////////////TEST1///////////////////////////////////////////////////
    std::cout << "TEST1: (Check for neighbor particle search)\n";
    std::vector<int> ans_list = {4,5,6};
    Tree::NeighborParticleSearchTree tree(DIM, reserved_size);
    tree.Resize(num_of_particles); //change actual number of particles to be searched for.
    for(int i = 0;i<num_of_particles;++i) //Arrange particles in a cube diagonal. The size of cube: size_x, size_y, size_z = [0,100]
        for(int dim = 0;dim<3;++dim)
            tree.CopyPos(i, i, dim); //copy pos[dim] to ith_particle[dim] inside the tree 
        
    tree.UpdateTree(); //construct tree. Needs to be called after update particle positions
    tree.FindNeighborParticle(point_of_search, search_radius, interaction_list); //Store the index of the particle in the region of search_radius from point_of_search into interaction_list

    if(interaction_list.size() != ans_list.size()) {
        std::cout << "TEST1 FAILED. there are more particles than answer\n";
        std::exit(EXIT_FAILURE);
    }
    std::sort(interaction_list.begin(), interaction_list.end());

    printf("particles found from (%f, %f, %f): \n", point_of_search[0], point_of_search[1], point_of_search[2]);
    //Answer will be id=4,5,6 if the proggram is not bugged.
    for(int i=0;i<interaction_list.size();++i) {
        int id = interaction_list[i];
        printf("id=%d, (%f,%f,%f)\n",id, tree.GetPos(id,0), tree.GetPos(id,1), tree.GetPos(id,2));
        if(id != ans_list[i]) {
            std::cout << "TEST1 FAILED. particle id: " << id << " is not neighobor particle.\n";
            std::exit(EXIT_FAILURE);
        }
    }
    std::cout << "TEST1 PASSED\n";
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n";

    //////////////////////////////////////////////TEST2///////////////////////////////////////////////////
    std::cout << "TEST2: (Check for neighbor particle search with periodic boundary):\n";
    ans_list = {0,1,2,3,4,94,95,96,97,98,99};
    point_of_search[0] = 1, point_of_search[1] = 1, point_of_search[2] = 1;
    search_radius = 10;
    constexpr double periodic_boundary_length[3] = {100,100,100}; //Box size of periodic boundary box. The size of box: size_x, size_y, size_z = [0,100]

    for(int i = 0;i<num_of_particles;++i) //Arrange particles in a cube diagonal
        for(int dim = 0;dim<DIM;++dim)
            tree.CopyPos(-i, i, dim); //copy pos[dim] to ith_particle[dim] inside the tree 

    tree.UpdateTree(); //construct tree. Needs to be called after update particle positions
    tree.FindNeighborParticleWithPeriodicBoundary(point_of_search, search_radius, periodic_boundary_length, interaction_list); //Store the index of the particle in the region of search_radius from point_of_search into interaction_list

    if(interaction_list.size() != ans_list.size()) {
        std::cout << "TEST2 FAILED. there are more particles than answer\n";
        std::exit(EXIT_FAILURE);
    }
    std::sort(interaction_list.begin(), interaction_list.end());
    printf("particles found from (%f, %f, %f): \n", point_of_search[0], point_of_search[1], point_of_search[2]);
    //Answer will be id=99, 98, 97, 96, 95, 94, 4, 3, 2, 1, 0 if the proggram is not bugged
    for(int i=0;i<interaction_list.size();++i) {
        int id = interaction_list[i];
        printf("id=%d, (%f,%f,%f)\n",id, tree.GetPos(id,0), tree.GetPos(id,1), tree.GetPos(id,2));
        if(id != ans_list[i]) {
            std::cout << "TEST2 FAILED. particle id: " << id << " is not neighobor particle.\n";
            std::exit(EXIT_FAILURE);
        }
    }
    std::cout << "TEST2 PASSED\n";
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "\n";

    //////////////////////////////////////////////TEST3///////////////////////////////////////////////////
    std::cout << "TEST3: (Check for Tree::SearchMode::SYMMETRY):\n";
    ans_list = {40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60};
    point_of_search[0] = -50, point_of_search[1] = -50, point_of_search[2] = -50;
    search_radius = 1;

    for(int i = 0;i<num_of_particles;++i) {//Arrange particles in a cube diagonal
        tree.CopySearchRadius(10.1*std::sqrt(3), i);
        for(int dim = 0;dim<DIM;++dim)
            tree.CopyPos(-i, i, dim); //copy pos[dim] to ith_particle[dim] inside the tree 
    }

    tree.UpdateTree(); //construct tree. Needs to be called after update particle positions
    tree.FindNeighborParticleWithPeriodicBoundary<Tree::SearchMode::SYMMETRY>(point_of_search, 0.01, periodic_boundary_length, interaction_list); //Store the index of the particle in the region of search_radius from point_of_search into interaction_list

    if(interaction_list.size() != ans_list.size()) {
        std::cout << "TEST3 FAILED. there are more particles than answer\n";
        std::exit(EXIT_FAILURE);
    }
    std::sort(interaction_list.begin(), interaction_list.end());
    printf("particles found from (%f, %f, %f): \n", point_of_search[0], point_of_search[1], point_of_search[2]);
    //Answer will be id=40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60 if the proggram is not bugged
    for(int i=0;i<interaction_list.size();++i) {
        int id = interaction_list[i];
        printf("id=%d, (%f,%f,%f)\n",id, tree.GetPos(id,0), tree.GetPos(id,1), tree.GetPos(id,2));
        if(id != ans_list[i]) {
            std::cout << "TEST3 FAILED. particle id: " << id << " is not neighobor particle.\n";
            std::exit(EXIT_FAILURE);
        }
    }
    std::cout << "TEST3 PASSED\n";
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "\n";

    //////////////////////////////////////////////TEST3///////////////////////////////////////////////////
    std::cout << "TEST4 (Check for move semantics): \n";
    Tree::NeighborParticleSearchTree tree_move_test = std::move(tree);
    tree_move_test.FindNeighborParticleWithPeriodicBoundary<Tree::SearchMode::SYMMETRY>(point_of_search, search_radius, periodic_boundary_length, interaction_list); //Store the index of the particle in the region of search_radius from point_of_search into interaction_list

    if(interaction_list.size() != ans_list.size()) {
        std::cout << "TEST4 FAILED. there are more particles than answer\n";
        std::cout << "TEST4 FAILED. Move semantics is not correctly implemented\n";
        std::exit(EXIT_FAILURE);
    }
    std::sort(interaction_list.begin(), interaction_list.end());
    printf("particles found from (%f, %f, %f): \n", point_of_search[0], point_of_search[1], point_of_search[2]);
    //Answer will be id=99, 98, 97, 96, 95, 94, 4, 3, 2, 1, 0 if the proggram is not bugged
    for(int i=0;i<interaction_list.size();++i) {
        int id = interaction_list[i];
        printf("id=%d, (%f,%f,%f)\n",id, tree_move_test.GetPos(id,0), tree_move_test.GetPos(id,1), tree_move_test.GetPos(id,2));
        if(id != ans_list[i]) {
            std::cout << "TEST4 FAILED. particle id: " << id << " is not neighobor particle.\n";
            std::cout << "Move semantics is not correctly implemented\n";
            std::exit(EXIT_FAILURE);
        }
    }

    Tree::NeighborParticleSearchTree tree_move_test2(2, 1000);
    tree_move_test = std::move(tree_move_test);
    tree = std::move(tree_move_test2);
}

    std::cout << "TEST4 PASSED\n";
    //////////////////////////////////////////////////////////////////////////////////////////////////////
}