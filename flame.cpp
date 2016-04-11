
#include"stdio.h"
#include"stdlib.h"
#include"string.h"
#include"math.h"
#include <utility>

#include"flame.h"

float Flame_Euclidean2( float *x, float *y, int m )
{
    float d = 0;
    int i;
    for(i=0; i<m; i++ ) d += ( x[i] - y[i] ) * ( x[i] - y[i] );
    return sqrt( d );
}

Flame::Flame(float **p_data, int p_n, int p_m )
{
    m_N = p_n;
    m_M = p_m;
//    printf( "m_N % i \n", m_N);

    m_KMAX = sqrt( m_N ) + 10;
    if( m_KMAX >= m_N ) m_KMAX = m_N - 1;
    m_indexFloatArray = new std::pair<int, float> [ m_N ];

    m_graph = new int*[ m_N ];
    m_dists = new float*[ m_N ];

    for( int i = 0; i < m_N; i++){
        m_graph[ i ] = new int[ m_KMAX ];
        m_dists[ i ] = new float[ m_KMAX ];
//            weights[i] = new float*[MAX];

        for( int j = 0; j < m_N; j++){
            m_indexFloatArray[ j ].first = j;
            m_indexFloatArray[ j ].second = Flame_Euclidean2( p_data[i], p_data[j], m_M );
        }

        partialQuickSort( 0, m_N-1, m_KMAX+1 );
        /* Store MAX number of nearest neighbors. */
        for( int j = 0; j < m_KMAX; j++){
            m_graph[ i ][ j ] = m_indexFloatArray[ j+1 ].first;
            m_dists[ i ][ j ] = m_indexFloatArray[ j+1 ].second;
        }
    }
    delete[] m_indexFloatArray;
}

Flame::~Flame()
{
    for ( int i = 0; i < m_N; i++ ) {
        delete[] m_graph[i];
        delete[] m_dists[i];
    }
    delete[] m_graph;
    delete[] m_dists;


}

void Flame::partialQuickSort( int first, int last, int part )
{
    int lower=first+1, upper=last;
    float pivot;
    std::pair<int, float> val;
    if( first >= last ) return;
    val = m_indexFloatArray[first];
    m_indexFloatArray[first] = m_indexFloatArray[ (first+last)/2 ];
    m_indexFloatArray[ (first+last)/2 ] = val;
    pivot = m_indexFloatArray[ first ].second;

    while( lower <= upper ){
        while( lower <= last && m_indexFloatArray[lower].second < pivot ) lower ++;
        while( pivot < m_indexFloatArray[upper].second ) upper --;
        if( lower < upper ){
            val = m_indexFloatArray[lower];
            m_indexFloatArray[lower] = m_indexFloatArray[upper];
            m_indexFloatArray[upper] = val;
            upper --;
        }
        lower ++;
    }
    val = m_indexFloatArray[first];
    m_indexFloatArray[first] = m_indexFloatArray[upper];
    m_indexFloatArray[upper] = val;
    if( first < upper-1 ) partialQuickSort( first, upper-1, part );
    if( upper >= part ) return;
    if( upper+1 < last ) partialQuickSort( upper+1, last, part );
}

void Flame::defineSupports( int p_knn, float p_thd )
{
    float *density = new float[ m_N ];
    float d, sum, sum2, fmin, fmax = 0.0;

    delete[] m_nncounts;
    m_nncounts = new int[ m_N ];
    m_weights = new float*[ m_N ];

    if( p_knn > m_KMAX ) p_knn = m_KMAX;
    m_K = p_knn;
    for( int i = 0; i < m_N; i++){
        /* To include all the neighbors that have distances equal to the
         * distance of the most distant one of the K-Nearest Neighbors */
        int k = p_knn;
        d = m_dists[i][p_knn-1];
        for( int j = p_knn; j < m_KMAX; j++)
            if( m_dists[i][j] == d ) k++; else break;
        m_nncounts[i] = k;
        m_weights[i] = new float[k];

        /* The definition of weights in this implementation is
         * different from the previous implementations where distances
         * or similarities often have to be transformed in some way.
         *
         * But in this definition, the weights are only dependent on
         * the ranking of distances of the neighbors, so it is more
         * robust against distance transformations. */
        sum = 0.5*k*(k+1.0);
        for( int j = 0; j < k; j++)
            m_weights[i][j] = (k-j) / sum;

        sum = 0.0;
        for( int j = 0; j < k; j++) sum += m_dists[i][j];
        density[i] = 1.0 / (sum + EPSILON);
    }
    sum = 0.0;
    sum2 = 0.0;
    for(int i = 0; i < m_N; i++){
        sum += density[i];
        sum2 += density[i] * density[i];
    }
    sum = sum / m_N;
    /* Density threshold for possible outliers. */
    p_thd = sum + p_thd * sqrt( sum2 / m_N - sum * sum );

    if ( m_obtypes != nullptr )
        delete m_obtypes;
    m_obtypes = new char[m_N];
    for( int i = 0; i < m_N; i++)
        m_obtypes[i] = OBT_NORMAL;


    m_cso_count = 0;
    for( int i = 0; i < m_N; i++){
        int k = m_nncounts[i];
        fmax = 0.0;
        fmin = density[i] / density[ m_graph[i][0] ];
        for( int j = 1; j < k; j++){
            d = density[i] / density[ m_graph[i][j] ];
            if( d > fmax ) fmax = d;
            if( d < fmin ) fmin = d;
            /* To avoid defining neighboring objects or objects close
             * to an outlier as CSOs.  */
            if( m_obtypes[ m_graph[i][j] ] ) fmin = 0.0;
        }
//        printf(" for el %i we got fmin %f and fmax %f \n", i, fmin, fmax);
        if( fmin >= 1.0 ){
            m_cso_count ++;
            m_obtypes[i] = OBT_SUPPORT;
        }else if( fmax <= 1.0 && density[i] < p_thd ){
            m_obtypes[i] = OBT_OUTLIER;
        }
    }
    delete[] density;
}



void Flame::localApproximation( int p_steps, float p_epsilon )
{
    m_fuzzyships = new float*[ m_N ];
    float **fuzzyships2 = new float*[ m_N ];
    char even = 0;
    double dev;

    int k = 0;
    for( int i = 0; i < m_N; i++){
        m_fuzzyships[i] = new float[ m_cso_count+1 ];
        fuzzyships2[i] = new float[ m_cso_count+1 ];
        for( int j = 0; j <= m_cso_count; j++)
            m_fuzzyships[i][j] = fuzzyships2[i][j] = 0;

        if( m_obtypes[i] == OBT_SUPPORT ){
            /* Full membership to the cluster represented by itself. */
            m_fuzzyships[i][k] = 1.0;
            fuzzyships2[i][k] = 1.0;
            k ++;
        }else if( m_obtypes[i] == OBT_OUTLIER ){
            /* Full membership to the outlier group. */
            m_fuzzyships[i][m_cso_count] = 1.0;
            fuzzyships2[i][m_cso_count] = 1.0;
        }else{
            /* Equal memberships to all clusters and the outlier group.
             * Random initialization does not change the results. */
            for( int j = 0; j <= m_cso_count; j++)
                m_fuzzyships[i][j] = fuzzyships2[i][j] = 1.0/(m_cso_count+1);
        }
    }

    for( int t = 0; t < p_steps; t++){
        printf("iter %i", t);
        dev = 0;
        for( int i = 0; i < m_N; i++){
            int knn = m_nncounts[i];
            int *ids = m_graph[i];
            float *wt = m_weights[i];
            float *fuzzy = m_fuzzyships[i];
            float **fuzzy2 = fuzzyships2;
            double sum = 0.0;
            if( m_obtypes[i] != OBT_NORMAL ) continue;
            if( even ){
                fuzzy = fuzzyships2[i];
                fuzzy2 = m_fuzzyships;
            }
            /* Update membership of an object by a linear combination of
             * the memberships of its nearest neighbors. */
            for( int j = 0; j <= m_cso_count; j++){
                fuzzy[j] = 0.0;
                for(k=0; k<knn; k++) fuzzy[j] += wt[k] * fuzzy2[ ids[k] ][j];
                dev += (fuzzy[j] - fuzzy2[i][j]) * (fuzzy[j] - fuzzy2[i][j]);
                sum += fuzzy[j];
            }
            for( int  j = 0; j <= m_cso_count; j++) fuzzy[j] = fuzzy[j] / sum;
        }
        even = ! even;
        if( dev < p_epsilon ) break;
    }

    /* update the membership of all objects to remove clusters
     * that contains only the CSO. */
    for( int i = 0; i < m_N; i++){
        int knn = m_nncounts[i];
        int *ids = m_graph[i];
        float *wt = m_weights[i];
        float *fuzzy = m_fuzzyships[i];
        float **fuzzy2 = fuzzyships2;
        for( int j = 0; j <= m_cso_count; j++){
            fuzzy[j] = 0.0;
            for(k=0; k<knn; k++) fuzzy[j] += wt[k] * fuzzy2[ ids[k] ][j];
            dev += (fuzzy[j] - fuzzy2[i][j]) * (fuzzy[j] - fuzzy2[i][j]);
        }
    }
    for( int i = 0; i < m_N; i++) delete[] fuzzyships2[i];
    delete[] fuzzyships2;
}


void Flame::makeClusters( float p_thd )
{
    int imax;
    int C = m_cso_count+1;
    float fmax;

//    IntArray *clust;
    if ( m_clusters != nullptr ) delete[] m_clusters;
    m_clusters = new int[ m_N ];
    m_indexFloatArray = new std::pair<int, float>[ m_N ];

    /* Sort objects based on the "entropy" of fuzzy memberships. */
    for( int i = 0; i < m_N; i++){
        m_indexFloatArray[i].first = i;
        m_indexFloatArray[i].second = 0.0;
        for( int j = 0; j < C; j++){
            float fs = m_fuzzyships[i][j];
            if( fs > EPSILON ) m_indexFloatArray[i].second -= fs * log( fs );
        }
    }
    partialQuickSort( 0, m_N-1, m_N );


//    self->clusters = (IntArray*) calloc( C, sizeof(IntArray) );
    if( p_thd <0 || p_thd > 1.0 ){
        /* Assign each object to the cluster
         * in which it has the highest membership. */
        for( int i = 0; i < m_N; i++ ){
            int id = m_indexFloatArray[i].first;
            fmax = 0;
            imax = -1;
            for( int j = 0; j < C; j++ ){
                if( m_fuzzyships[id][j] > fmax ){
                    imax = j;
                    fmax = m_fuzzyships[id][j];
                }
            }
            m_clusters[i] = imax;
        }
//    }else{
//        /* Assign each object to all the clusters
//         * in which it has membership higher than thd,
//         * otherwise, assign it to the outlier group.*/
//        for(i=0; i<m_N; i++){
//            int id = m_indexFloatArray[i].first;
//            imax = -1;
//            for(j=0; j<C; j++){
//                if( m_fuzzyships[id][j] > p_thd || ( j == C-1 && imax <0 ) ){
//                    imax = j;
//                    clust = self->clusters + j;
//                    IntArray_Push( self->clusters + j, id );
//                }
//            }
//        }
    }
    /* removing empty clusters */

//    C = 0;
//    for( int i = 0; i < m_cso_count; i++){
//        if( self->clusters[i].size >0 ){
//            self->clusters[C] = self->clusters[i];
//            C ++;
//        }
//    }
    /* keep the outlier group, even if its empty */
//    self->clusters[C] = self->clusters[m_cso_count];
//    C ++;
//    for( int i = C; i < m_cso_count+1; i++)
//        memset( self->clusters+i, 0, sizeof(IntArray) );
//    m_count = C;
    delete[] m_indexFloatArray;
}

