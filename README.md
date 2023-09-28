# Extremely Simple Neighbor Particle Search Tree
header-only C++ library. \
Periodic boundary is supported. \
Test (Demo) is shown in test_neighbor_particle_search_tree.cpp

## Instalation
```sh
git clone https://github.com/StaticTaku/neighbor_particle_search_tree.git
```

## Usage
Create an object, resize, copy particle positions, update tree, then you can find neighbor particles near any points! \
If positions of particles change, then you need to resizse (if number of particles changes), copy particle positions that has changed its positions, update tree, then you can find neighbor particles near any points!

### Example
```c++
    Tree::NeighborParticleSearchTree tree(DIM, reserved_size); //reserved_size is number of particles allocated inside tree in advance
    tree.Resize(num_of_particles); //change actual number of particles to be searched for. Note that num_of_particles <=  reserved_size.
    for(int i = 0;i<num_of_particles;++i) //Arrange particles in a cube diagonal. The size of cube: size_x, size_y, size_z = [0,100]
        for(int dim = 0;dim<DIM;++dim)
            tree.CopyPos(pos[i][dim], i, dim); //copy pos[i][dim] to ith_particle[dim] inside the tree 
        
    tree.UpdateTree(); //construct tree. Needs to be called after update particle positions
    tree.FindNeighborParticle(point_of_search, search_radius, interaction_list); //Store the index of the particle in the region of search_radius from point_of_search into interaction_list
```

## Search Option
### GATHER (Default)
Search radius is finite. \
Particles within a radius `search_radius` centered on `point_of_search` are registered as neighboring particles in interaction_list.
```c++
tree.FindNeighborParticle(point_of_search, search_radius, interaction_list);

tree.FindNeighborParticleWithPeriodicBoundary(point_of_search, search_radius, periodic_boundary_length, interaction_list);
```
or
```c++
tree.FindNeighborParticle<Tree::SearchMode::GATHER>(point_of_search, search_radius, interaction_list);

tree.FindNeighborParticle<Tree::SearchMode::GATHER>(point_of_search, search_radius, periodic_boundary_length, interaction_list);
```
### SYMMETRY
Search radius is finite. \
Need to set search radius of all particles `search_radius_list` before `tree.UpdateTree()`.
```c++
for(int i = 0;i<num_of_particles;++i)
    tree.CopySearchRadius(search_radius_list[i], i);
```
1. Particles within a radius `search_radius` centered on `point_of_search` and 
2. Any particles `j` that covers the position `point_of_search` within a radius `search_radius_list[j]` centered on `position[j]` 

are registered as neighboring particles in interaction_list.
```c++
tree.FindNeighborParticle<Tree::SearchMode::SYMMETRY>(point_of_search, search_radius, interaction_list);

tree.FindNeighborParticleWithPeriodicBoundary<Tree::SearchMode::SYMMETRY>(point_of_search, search_radius, periodic_boundary_length, interaction_list);
```

## how to run test (or demo)
```sh
g++ -std=c++17 -O3 test_neighbor_particle_search_tree.cpp && ./a.out
```
