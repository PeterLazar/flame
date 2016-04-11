#ifndef FLAME_H
#define FLAME_H

#define EPSILON 1E-9

#include <utility>

enum FlameObjectTypes
{
    OBT_NORMAL ,
    OBT_SUPPORT ,
    OBT_OUTLIER
};

class Flame
{
public:
    void partialQuickSort( int first, int last, int part );
    void defineSupports( int knn, float thd );
    void localApproximation(int p_steps, float p_epsilon );
    void makeClusters( float p_thd );
    int getCsoCount() { return m_cso_count; }
    int getN() { return m_N; }
    int* getClusters() { return m_clusters; }

    Flame() {}// a bunch of nullptr's
    Flame(float** p_data, int p_n , int p_m);
    ~Flame();

private:
    int m_N;
    int m_M;
    int m_K;
    int m_KMAX;

    int   **m_graph;
    float **m_dists;

    int    *m_nncounts;
    float **m_weights;

    int m_cso_count;
    char *m_obtypes;

    float **m_fuzzyships;

    int m_count;
    int *m_clusters;

    std::pair<int, float>* m_indexFloatArray;
};



#endif // FLAME_H
